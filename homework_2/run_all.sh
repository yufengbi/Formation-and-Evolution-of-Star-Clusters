#!/bin/bash
# =============================================================================
#  作业 2：对比"有无黑洞"两类星团系统的演化
#  ---------------------------------------------------------------------------
#  课堂练习要求（见 slide）：
#    * 用 mcluster 生成两组系统，Plummer 模型，R_h = 1 pc
#        - 单质量系统（质量全为 1 Msun）              -> "无黑洞"参照组
#        - 双质量系统（质量分两组 1 Msun / 10 Msun）   -> "含黑洞"组
#          调节轻/重占比：1:9 与 5:5
#    * PeTar 关闭恒星演化与外引力势
#    * 模拟 2-3 t_rh
#    * 绘制 r_h 与 r_c 演化图，对比两类系统：
#      是否发生质量迁移（质量分层）与核心塌缩
#
#  本脚本职责：生成初始条件 -> N 体积分 -> 能量自检 -> 后处理。
#  分析绘图由 analyze.py / make_report.py 完成。
#
#  【为什么需要"看门狗"】
#  含 10 Msun 重星的系统会形成极紧的多体系统，PeTar 的 SDAR 时间变换可能退化：
#  日志里 |Int_err_cum/E| 会从 1e-9 一路涨到 0.5 以上，Step_count 触到
#  --ar-max-nstep 上限（1e6），此后模拟几乎不再前进——程序不报错、不退出，
#  只是一直烧 CPU。仅靠"跑完再检查能量"发现不了这种情况，所以每次积分都有
#  看门狗盯着：①累计积分误差超限；②超过 STALL_SEC 秒没有产生新快照。
#  触发任一条件就杀掉该次积分并重抽种子重跑。
#
#  运行：
#     bash run_all.sh                                    # 三种配置 x 3 个实现
#     CONFIGS="B_heavy10 C_heavy50" bash run_all.sh      # 只跑指定配置
#     NREAL=1 T_END=5 bash run_all.sh                    # 快速冒烟测试
#     SEEDS="1 2 3" bash run_all.sh                      # 指定种子（可复现）
# =============================================================================

set -u

# ---------------------------------------------------------------------------
# 0. 参数
# ---------------------------------------------------------------------------
N_PART="${N_PART:-1000}"          # 粒子（恒星）数
R_HALF="${R_HALF:-1}"             # 半质量半径 [pc]
# 积分时长：默认按每个系统各自的弛豫时间自动确定，t_end = T_RH_FACTOR * t_rh。
#   题目要求模拟 2-3 t_rh；本作业取 2 t_rh（下限），原因见 README §2.6：
#   含 10 Msun 重星的系统在约 2 t_rh 之后 SDAR 时间变换开始频繁退化，
#   更长的积分窗口在数值上已不可信（实测失败率从 t_end=20 Myr 时的 ~44% 
#   升到 t_end=50 Myr 时的 ~90%）。
T_RH_FACTOR="${T_RH_FACTOR:-2.0}" # t_end = T_RH_FACTOR * t_rh
T_END="${T_END:-}"                # 若显式设置则覆盖自动值（冒烟测试用）
T_OUT="${T_OUT:-0.5}"             # 快照输出间隔 [Myr]
NTHREAD="${NTHREAD:-4}"           # 每个任务的 OpenMP 线程数
NREAL="${NREAL:-3}"               # 每个配置的随机实现数
MAX_TRY="${MAX_TRY:-40}"          # 单个实现的尝试上限
NJOBS="${NJOBS:-4}"               # 同时并行的任务数
DT_SOFT="${DT_SOFT:-0}"           # 0 = PeTar 确定性自适应公式（见 README）
G_AU=0.00449830997959438          # 天文单位制引力常数 (Msun, pc, Myr)

MFRAC="${MFRAC:-0.01,0.05,0.1,0.2,0.3,0.5,0.7,0.9,0.99}"

# 能量守恒自检判据
TOL_DE="${TOL_DE:-0.05}"           # |dE/E| 上限
TOL_MODIFY="${TOL_MODIFY:-0.01}"   # max|Modify| / |E0| 上限
MAX_DUMP="${MAX_DUMP:-10}"         # 诊断文件数兜底上限

# 看门狗
STALL_SEC="${STALL_SEC:-360}"      # 多少秒没有新快照就判定为卡死
ERR_CUM_MAX="${ERR_CUM_MAX:-0.1}"  # |Int_err_cum/E| 上限

SEEDS="${SEEDS:-}"                 # 可选：空格分隔种子列表，按 run 序号分配
CONFIGS="${CONFIGS:-A_single B_heavy10 C_heavy50}"

WORKDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTROOT="${WORKDIR}/output"
EVIDIR="${WORKDIR}/evidence"
LOGDIR="${WORKDIR}/log"
ATTEMPTS="${WORKDIR}/attempts.tsv"

# ---------------------------------------------------------------------------
# 1. 配置表
# ---------------------------------------------------------------------------
#  A_single  : 1000 x 1 Msun                M = 1000 Msun  （"无黑洞"）
#  B_heavy10 : 900 x 1 + 100 x 10 Msun      M = 1900 Msun  （重星占数量 10%）
#  C_heavy50 : 500 x 1 + 500 x 10 Msun      M = 5500 Msun  （重星占数量 50%）
#
#  mcluster 关键选项：
#    -P 0   Plummer 剖面        -R 1   半质量半径 1 pc
#    -f 0   无 IMF，全为 1 Msun
#    -f 6   多组分：-m 给各组质量，-a 给各组的**数量占比**
#           （源码 main.c: noffset[j+1] = noffset[j] + N*alpha[j]）
#    -C 5   输出 NBODY6++ 格式（PeTar 可读）
#    -u 1   天文单位 (Msun, pc, km/s)
#    -t 0   无潮汐外场（与 PeTar 端"关闭外引力势"一致）
mcluster_args() {
    case "$1" in
        A_single)
            echo "-N $N_PART -R $R_HALF -P 0 -f 0 -C 5 -u 1 -t 0" ;;
        B_heavy10)
            echo "-N $N_PART -R $R_HALF -P 0 -f 6 -m 1.0 -a 0.9 -m 10.0 -a 0.1 -C 5 -u 1 -t 0" ;;
        C_heavy50)
            echo "-N $N_PART -R $R_HALF -P 0 -f 6 -m 1.0 -a 0.5 -m 10.0 -a 0.5 -C 5 -u 1 -t 0" ;;
        *)  echo "" ;;
    esac
}

config_title() {
    case "$1" in
        A_single)  echo "单质量 1 Msun（无黑洞参照）" ;;
        B_heavy10) echo "双质量 1/10 Msun，重星占数量 10%" ;;
        C_heavy50) echo "双质量 1/10 Msun，重星占数量 50%" ;;
    esac
}

# ---------------------------------------------------------------------------
# 2. 环境检查
# ---------------------------------------------------------------------------
for cmd in mcluster petar petar.init petar.select petar.data.gether \
           petar.data.process python3 awk stat; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "错误: 找不到命令 '$cmd'。请先安装 PeTar 并把 \$prefix/bin 加入 PATH。" >&2
        exit 1
    fi
done

if [ "$(ulimit -s)" != "unlimited" ]; then
    ulimit -s unlimited 2>/dev/null || echo "警告: 无法设置 ulimit -s unlimited" >&2
fi
export OMP_STACKSIZE="${OMP_STACKSIZE:-256M}"

mkdir -p "$OUTROOT" "$EVIDIR" "$LOGDIR"
[ -f "$ATTEMPTS" ] || printf "config\trun\tseed\tverdict\tt_end\tdEoverE\tmaxModify\tnote\n" > "$ATTEMPTS"

# ---------------------------------------------------------------------------
# 3. 辅助函数
# ---------------------------------------------------------------------------
record_attempt() {   # config run seed verdict t_end de maxModify note
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" >> "$ATTEMPTS"
}

save_evidence() {    # dir tag note
    local dir="$1" tag="$2" note="$3"
    local ev="${EVIDIR}/${tag}"
    mkdir -p "$ev"
    for f in output mc.info REJECT.txt; do
        [ -f "${dir}/${f}" ] && cp -f "${dir}/${f}" "$ev/"
    done
    cp -f "${dir}"/data.hard_large_energy_* "${ev}/" 2>/dev/null
    cp -f "${dir}"/data.dump_large_step_*    "${ev}/" 2>/dev/null
    [ -f "${ev}/output" ] && mv -f "${ev}/output" "${ev}/petar.log" 2>/dev/null
    printf "%s\n" "$note" > "${ev}/NOTE.txt"
    echo "      证据 -> ${ev}"
}

newest_snapshot_mtime() {
    local m=0 f t
    for f in data.[0-9]*; do
        [ -f "$f" ] || continue
        t=$(stat -c %Y "$f" 2>/dev/null) || continue
        [ "$t" -gt "$m" ] && m=$t
    done
    echo "$m"
}

# ---------------------------------------------------------------------------
# 4. 单次尝试：生成 -> 积分（带看门狗）-> 能量自检
#    返回 0 = 通过；1 = 自检不合格；2 = petar 崩溃；3 = 看门狗中止
# ---------------------------------------------------------------------------
run_attempt() {
    local cfg="$1" ridx="$2" seed="$3"

    # 安全护栏：确保 rm -rf 只作用于 output/<已知配置>/run<数字>
    case " $CONFIGS " in *" ${cfg} "*) ;; *) echo "错误: 未知配置 '${cfg}'"; return 1 ;; esac
    case "$ridx" in ''|*[!0-9]*) echo "错误: 非法 run 序号 '${ridx}'"; return 1 ;; esac

    local dir="${OUTROOT}/${cfg}/run${ridx}"
    local tag="${cfg}_run${ridx}_seed${seed}"

    rm -rf "$dir"; mkdir -p "$dir"; cd "$dir" || return 1

    # --- 4.1 初始条件 -------------------------------------------------------
    if ! mcluster $(mcluster_args "$cfg") -s "$seed" -o mc >mc.log 2>&1; then
        echo "      mcluster 失败"
        record_attempt "$cfg" "$ridx" "$seed" "crash" - - - "mcluster 失败"; return 1
    fi
    local infile
    infile=$(ls mc.dat.* 2>/dev/null | head -1)
    [ -z "$infile" ] && { echo "      mcluster 未生成 mc.dat.*"; return 1; }

    # --- 4.2 单位换算 (km/s -> pc/Myr) -------------------------------------
    if ! petar.init -v kms2pcmyr -f input "$infile" >init.log 2>&1; then
        echo "      petar.init 失败"; return 1
    fi

    # --- 4.3 由初始条件确定积分时长 t_end = T_RH_FACTOR * t_rh ------------
    local FEAT TRH TEND
    FEAT=$(python3 "${WORKDIR}/calc_trh.py" input "$T_RH_FACTOR") || {
        echo "      calc_trh.py 失败"; return 1; }
    TRH=$(echo "$FEAT" | awk '{print $5}')
    TEND=$(echo "$FEAT" | awk '{print $6}')
    [ -n "$T_END" ] && TEND="$T_END"
    echo "$FEAT" > features.txt        # N mbar Rh sigma3D trh tend
    echo "      t_rh = ${TRH} Myr  ->  t_end = ${TEND} Myr"

    # --- 4.4 选择最优二进制 ------------------------------------------------
    petar.select --optional mpi,omp,avx512,avx2 >/dev/null 2>&1

    # --- 4.5 N 体积分（后台 + 看门狗） -------------------------------------
    #  -i 1 : 输入 ASCII（petar.init 产物）、快照输出 ASCII，便于自定义分析
    #  -s 0 : dt_soft 由 PeTar 按 2.6e-4*G*M/sigma_3D^3 自动确定并规整到 2^-n
    echo "      积分中（dt_soft=${DT_SOFT}, OMP=${NTHREAD}）... $(date +%H:%M:%S)"
    OMP_NUM_THREADS="$NTHREAD" OMP_STACKSIZE=256M \
        petar -u 1 -t "$TEND" -o "$T_OUT" -s "$DT_SOFT" -i 1 input >output 2>&1 &
    local ppid=$!

    local last_m last_change now ec nm gone="" note=""
    last_m=$(newest_snapshot_mtime)
    last_change=$(date +%s)

    while kill -0 "$ppid" 2>/dev/null; do
        sleep 5

        # ① 累计积分误差超限 —— SDAR 时间变换已退化，继续跑没有意义
        ec=$(grep -oE '\|Int_err_cum/E\|:[ ]*[0-9.eE+-]+' output 2>/dev/null \
             | tail -1 | awk '{print $2}')
        if [ -n "$ec" ] && awk -v v="$ec" -v m="$ERR_CUM_MAX" 'BEGIN{exit !(v>m)}'; then
            gone="watchdog"; note="积分器发散 |Int_err_cum/E| = ${ec} > ${ERR_CUM_MAX}"
            break
        fi

        # ② 长时间没有新快照 —— 卡死
        nm=$(newest_snapshot_mtime)
        if [ "$nm" -gt "$last_m" ]; then last_m="$nm"; last_change=$(date +%s); fi
        now=$(date +%s)
        if [ $((now - last_change)) -gt "$STALL_SEC" ]; then
            gone="watchdog"; note="卡死：${STALL_SEC} s 内无新快照"
            break
        fi
    done

    if [ -n "$gone" ]; then
        kill -TERM "$ppid" 2>/dev/null; sleep 2; kill -KILL "$ppid" 2>/dev/null
        wait "$ppid" 2>/dev/null
        echo "      ⏱  ${note}"
        record_attempt "$cfg" "$ridx" "$seed" "watchdog" - - - "$note"
        save_evidence "$dir" "$tag" "$note"
        return 3
    fi

    local pst=0
    wait "$ppid" || pst=$?
    if [ "$pst" -ne 0 ]; then
        echo "      ❌ petar 非正常退出（退出码 ${pst}）"
        record_attempt "$cfg" "$ridx" "$seed" "crash" - - - "petar 退出码 ${pst}"
        save_evidence "$dir" "$tag" "petar 异常退出，退出码 ${pst}"
        return 2
    fi

    # --- 4.6 能量守恒自检 ---------------------------------------------------
    #  日志格式: Time: <t> ...
    #            Physic: <Error/Total> <Error> <Error_cum> <Total E> <Kinetic>
    #                    <Potential> <Modify> ...
    #  -> awk 字段 $5 = <Total E>，$8 = <Modify>（与作业 1 一致）
    local CHK
    CHK=$(awk '
        /^Time:/{t=$2}
        /^Physic:/{
            if (!init) { e0=$5; init=1 }
            e=$5; m=$8; if (m<0) m=-m; if (m>mmax) mmax=m; et=t
        }
        END{
            if (init) printf "%.6f %.6f %.3f %.6g\n", e0, e, et, mmax
            else print "0 0 -1 0"
        }' output)

    local E0 E1 TT MM ndump
    E0=$(echo "$CHK" | cut -d' ' -f1)
    E1=$(echo "$CHK" | cut -d' ' -f2)
    TT=$(echo "$CHK" | cut -d' ' -f3)
    MM=$(echo "$CHK" | cut -d' ' -f4)
    ndump=$(ls data.hard_large_energy_* data.dump_large_step_* 2>/dev/null | wc -l)

    local VERDICT
    VERDICT=$(awk -v e0="$E0" -v e="$E1" -v m="$MM" -v t="$TT" -v nd="$ndump" \
        -v tend="$TEND" -v tde="$TOL_DE" -v tmod="$TOL_MODIFY" -v ndmax="$MAX_DUMP" 'BEGIN{
        if (e0==0) { print "REJECT@@初始能量为零"; exit }
        de = (e-e0)/e0; if (de<0) de=-de
        mr = m/(e0<0?-e0:e0)
        if (t < tend-0.5) { printf "REJECT@@未跑完 t=%.1f < %.1f Myr", t, tend; exit }
        if (nd > ndmax)   { printf "REJECT@@诊断文件过多 %d", nd; exit }
        if (mr > tmod)    { printf "REJECT@@maxModify/|E0| = %.3g > %.3g", mr, tmod; exit }
        if (de > tde)     { printf "REJECT@@|dE/E| = %.3g > %.3g", de, tde; exit }
        printf "PASS@@t=%.1f Myr  |dE/E|=%.3e  maxModify/|E0|=%.3e  dump=%d", t, de, mr, nd
    }')

    printf "      t=%.1f  E0=%.2f  E=%.2f  maxModify=%.3g  dump=%d\n" \
           "$TT" "$E0" "$E1" "$MM" "$ndump"

    local de
    de=$(awk -v a="$E0" -v b="$E1" 'BEGIN{d=(b-a)/a; if(d<0)d=-d; printf "%.3e", d}')

    case "$VERDICT" in
        PASS*)
            echo "      ✅ ${VERDICT#*@@}"
            record_attempt "$cfg" "$ridx" "$seed" "pass" "$TT" "$de" "$MM" "dump=$ndump"
            return 0 ;;
        *)
            echo "      ❌ ${VERDICT#*@@}"
            record_attempt "$cfg" "$ridx" "$seed" "reject" "$TT" "$de" "$MM" "${VERDICT#*@@}"
            save_evidence "$dir" "$tag" "${VERDICT#*@@}"
            return 1 ;;
    esac
}

# ---------------------------------------------------------------------------
# 5. 后处理
# ---------------------------------------------------------------------------
postprocess() {
    local cfg="$1" ridx="$2"
    local dir="${OUTROOT}/${cfg}/run${ridx}"
    cd "$dir" || return 1
    rm -f data.lagr data.core data.esc data.esc_single data.esc_binary \
          data.snap.lst data.parallel.*
    petar.data.gether data >gether.log 2>&1
    petar.data.process -G "$G_AU" -s ascii -e -m "$MFRAC" data.snap.lst >process.log 2>&1
    rm -f data.parallel.*
    [ -f data.lagr ] && [ -f data.core ]
}

# ---------------------------------------------------------------------------
# 6. 驱动：对 (配置, 实现序号) 反复尝试直到通过自检
# ---------------------------------------------------------------------------
drive() {
    local cfg="$1" ridx="$2"
    local dir="${OUTROOT}/${cfg}/run${ridx}"
    local attempt=0

    echo "===== [$(date +%H:%M:%S)] ${cfg}_run${ridx} : $(config_title "$cfg") ====="

    while [ "$attempt" -lt "$MAX_TRY" ]; do
        attempt=$((attempt+1))

        local seed
        if [ -n "$SEEDS" ]; then
            seed=$(echo "$SEEDS" | awk -v k="$ridx" '{print $k}')
            [ -z "$seed" ] && seed=$(echo "$SEEDS" | awk '{print $1}')
            case "$cfg" in B_*) seed=$((seed+1000));; C_*) seed=$((seed+2000));; esac
        else
            seed=$(od -An -N2 -tu2 /dev/urandom | awk '{print $1 % 30000 + 1}')
        fi

        echo "   [${attempt}/${MAX_TRY}] seed = ${seed}   $(date +%H:%M:%S)"
        run_attempt "$cfg" "$ridx" "$seed"
        local st=$?

        if [ "$st" -eq 0 ]; then
            local chk dt feat TRH TEND
            chk=$(awk '/^Time:/{t=$2} /^Physic:/{if(!i){e0=$5;i=1} e=$5; if($8<0)m=-$8; else m=$8; if(m>mm){mm=m}} END{printf "%.6f %.6f %.3f %.6g", e0, e, t, mm}' "${dir}/output")
            dt=$(grep -m1 '^F s ' "${dir}/data.par" | awk '{print $3}')
            feat=$(cat "${dir}/features.txt" 2>/dev/null || echo "- - - - - -")
            TRH=$(echo "$feat" | awk '{print $5}')
            TEND=$(echo "$feat" | awk '{print $6}')
            cat > "${dir}/RUN_INFO.txt" <<EOF
config          = ${cfg}
title           = $(config_title "$cfg")
seed            = ${seed}
N               = ${N_PART}
R_half          = ${R_HALF} pc (Plummer)
t_rh            = ${TRH} Myr   (由本实初始条件计算)
t_end           = ${TEND} Myr  (= ${T_RH_FACTOR} t_rh，题目要求 2-3 t_rh)
output_interval = ${T_OUT} Myr
dt_soft(hex)    = ${dt}
OMP_NUM_THREADS = ${NTHREAD}
attempts        = ${attempt}
E0 / E(t_end)   = $(echo "$chk" | awk '{printf "%.4f / %.4f", $1, $2}')
max|Modify|     = $(echo "$chk" | awk '{printf "%.6g", $4}')
EOF
            echo "      后处理 ..."
            postprocess "$cfg" "$ridx"
            echo "      ✅ 完成 (seed=${seed}, 共尝试 ${attempt} 次)   $(date +%H:%M:%S)"
            return 0
        fi
    done

    echo "      ❌ ${cfg}_run${ridx}: ${MAX_TRY} 次尝试均未通过" >&2
    return 1
}

# ---------------------------------------------------------------------------
# 7. 主流程：所有 (配置, 实现) 任务进队列，按 NJOBS 并发
# ---------------------------------------------------------------------------
echo "============================================================"
echo " 作业2 : 有无黑洞两类星团系统的演化对比"
echo " 配置          = ${CONFIGS}"
echo " 每组实现数    = ${NREAL}"
echo " N             = ${N_PART}     R_half = ${R_HALF} pc (Plummer)"
echo " t_end         = ${T_RH_FACTOR} * t_rh（每个系统各自计算）${T_END:+, 覆盖为 ${T_END} Myr}"
echo " 输出间隔      = ${T_OUT} Myr"
echo " dt_soft       = ${DT_SOFT} (0 = PeTar 自适应)"
echo " 线程/任务     = ${NTHREAD}    并行任务 = ${NJOBS}"
echo " 看门狗        = 卡死 ${STALL_SEC} s / |Int_err_cum/E| > ${ERR_CUM_MAX}"
echo "============================================================"

running=0
pids=""
for cfg in $CONFIGS; do
    for r in $(seq 1 "$NREAL"); do
        drive "$cfg" "$r" > "${LOGDIR}/${cfg}_run${r}.log" 2>&1 &
        pids="$pids $!"
        running=$((running+1))
        if [ "$running" -ge "$NJOBS" ]; then
            wait -n 2>/dev/null || true
            running=$((running-1))
        fi
    done
done

for p in $pids; do wait "$p" || true; done

echo ""
echo "===== 尝试统计（attempts.tsv） ====="
awk -F'\t' 'NR>1{n[$1"_"$4]++; if($5!="-"){if(!(($1) in mx) || $5>mx[$1]) mx[$1]=$5}} END{}' "$ATTEMPTS" 2>/dev/null
printf "%-14s %6s %8s %8s %8s\n" "配置" "通过" "自检拒绝" "看门狗" "崩溃"
for cfg in $CONFIGS; do
    p=$(awk -F'\t' -v c="$cfg" '$1==c && $4=="pass"' "$ATTEMPTS" | wc -l)
    rj=$(awk -F'\t' -v c="$cfg" '$1==c && $4=="reject"' "$ATTEMPTS" | wc -l)
    wd=$(awk -F'\t' -v c="$cfg" '$1==c && $4=="watchdog"' "$ATTEMPTS" | wc -l)
    cr=$(awk -F'\t' -v c="$cfg" '$1==c && $4=="crash"' "$ATTEMPTS" | wc -l)
    printf "%-14s %6d %8d %8d %8d\n" "$cfg" "$p" "$rj" "$wd" "$cr"
done

echo ""
ok=1
for cfg in $CONFIGS; do
    p=$(awk -F'\t' -v c="$cfg" '$1==c && $4=="pass"' "$ATTEMPTS" | wc -l)
    [ "$p" -ge "$NREAL" ] || ok=0
done
if [ "$ok" -eq 1 ]; then
    echo "✅ 全部完成。产物在 ${OUTROOT}/"
else
    echo "⚠️  有实现未完成，见 ${ATTEMPTS} 与 ${LOGDIR}/"
fi
