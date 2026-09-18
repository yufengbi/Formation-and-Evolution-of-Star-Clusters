#!/bin/bash
# =============================================================================
#  作业 1：疏散星团（star cluster）数值模拟
#  ---------------------------------------------------------------------------
#  要求：
#    * 借助 PeTar / sample 文件夹中的案例 script 实现星团数值模拟
#    * N = 1000
#    * Plummer profile（Plummer 密度分布）
#    * 演化到 200 Myr
#    * 用 petar.movie 绘制星团动画与拉格朗日半径演化
#
#  本脚本在 PeTar 官方示例 sample/star_cluster_plummer_N1k.sh 的基础上改写。
#
#  【关于随机种子与自动重试】
#  官方示例用 mcluster 默认的 -s 0，即"按时间随机"，每次得到不同的星团实现，
#  结果不可复现（连作者也无法重现上一次的结果）。
#
#  但这里有一个必须处理的数值问题：本 setup（N=1000、Kroupa IMF 至 150 Msun
#  的重星、无恒星演化、无潮汐场）会让星团核心深度坍缩并形成稳定的少体（多重）
#  系统，而 PeTar 的 SDAR 方法对这类系统**不保证经典能量守恒**（PeTar README
#  "Significant Hard Energy" 一节明确说明）。失稳时系统会凭空获得能量、由束缚
#  变成不束缚，但程序不会报错、动画照样能出——只有检查能量才能发现。
#
#  实测 7 个随机实现中只有约 2 个可用（成功率 ~30%）。
#
#  因此本脚本采用：**随机抽取种子 -> 模拟 -> 能量守恒自检 -> 不合格自动重抽
#  重跑**，直到拿到一个数值上可信的实现。实际使用的种子会写入 output/SEED.txt，
#  因此结果仍可复现（设 SEED_OVERRIDE=<种子> 即可重现同一次运行）。
#
#  运行方式：  bash run_cluster.sh
#  固定种子：  SEED_OVERRIDE=3 bash run_cluster.sh
# =============================================================================

set -u

# ---------------------------------------------------------------------------
# 0. 参数设置
# ---------------------------------------------------------------------------
N_PARTICLE=1000          # 粒子（恒星）数
R_HALF=1                 # 半质量半径 [pc]
T_END=200.0              # 积分终止时间 [Myr]  <-- 作业要求 200 Myr
T_OUT=1.0                # 快照输出间隔 [Myr] -> 共 201 个快照
NTHREAD=4                # OpenMP 线程数
PREFIX=data              # 输出文件名前缀
MAX_TRY=10               # 自动重试上限

# 树时间步 dt_soft [Myr]
# 官方示例用 petar.find.dt 自动调优，但该工具是按**墙钟时间**挑选最优值的：
# 机器上有其它负载时它会给出偏小很多的结果（实测从 0.0039 退化到 0.00024，
# 即慢 16 倍），而且不可复现。因此这里固定为官方示例在空载机器上得到的最优
# 值的一半（与 sample 脚本的最终取值一致）。可用 DT_SOFT 环境变量覆盖。
DT_SOFT="${DT_SOFT:-0.001953125}"

# 能量守恒自检判据
#  物理有效性的判据是**全局能量守恒**（见文件开头说明）。
#  PeTar 的 hard_large_energy 诊断文件表示硬积分中出现了较大的局部能量误差，
#  但官方文档指出：这类误差若能偶发出现，"通常不会显著影响系统的整体动力学演化"。
#  因此把 dump 文件数作为**警告**报告，只用一个较宽松的上限做兜底。
TOL_DE=0.05              # |dE/E| 上限
TOL_MODIFY=0.01          # max|Modify| / |E0| 上限
MAX_DUMP=10              # hard_large_energy 诊断文件数的兜底上限

SEED_OVERRIDE="${SEED_OVERRIDE:-}"   # 可选：固定种子

# 天文单位制下的引力常数 (Msun, pc, km/s) -> 后处理保持单位一致
G_AU=0.00449830997959438

WORKDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${WORKDIR}/output"
EVIDIR="${WORKDIR}/evidence"

# ---------------------------------------------------------------------------
# 1. 环境检查
# ---------------------------------------------------------------------------
for cmd in mcluster petar petar.init petar.select petar.find.dt \
           petar.data.gether petar.data.process petar.movie python3; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "错误: 找不到命令 '$cmd'。请先安装 PeTar 并把 \$prefix/bin 加入 PATH。" >&2
        exit 1
    fi
done

# 栈上限：PeTar 深递归需要足够大的栈，否则会 segmentation fault
if [ "$(ulimit -s)" != "unlimited" ]; then
    ulimit -s unlimited 2>/dev/null || echo "警告: 无法设置 ulimit -s unlimited" >&2
fi
export OMP_STACKSIZE="${OMP_STACKSIZE:-128M}"

# ---------------------------------------------------------------------------
# 2. 辅助函数
# ---------------------------------------------------------------------------
# 保存失稳实现的诊断证据，便于在报告中说明
save_evidence() {
    local tag="$1"
    local dst="${EVIDIR}/${tag}"
    mkdir -p "$dst"
    cp -f ${PREFIX}.hard_large_energy_* "$dst"/ 2>/dev/null
    cp -f mc.info output "$dst"/ 2>/dev/null
    [ -f output ] && mv -f "$dst/output" "$dst/petar.log" 2>/dev/null
    echo "      证据已保存到 ${dst}"
}

# 抽取随机种子：1..30000。
# 相比 mcluster 的 -s 0（按时间随机），这样抽出的种子会被记录，运行仍可复现。
draw_seed() {
    if [ -n "$SEED_OVERRIDE" ]; then
        echo "$SEED_OVERRIDE"
    else
        od -An -N2 -tu2 /dev/urandom | awk '{print $1 % 30000 + 1}'
    fi
}

echo "============================================================"
echo " 作业1: 星团数值模拟"
echo " N            = $N_PARTICLE"
echo " 半质量半径   = $R_HALF pc (Plummer)"
echo " 积分时间     = $T_END Myr"
echo " 输出间隔     = $T_OUT Myr"
echo " 随机种子     = ${SEED_OVERRIDE:-随机抽取（自动重试直到通过能量自检）}"
echo " 输出目录     = $OUTDIR"
echo "============================================================"

# ---------------------------------------------------------------------------
# 3. 主循环：生成 -> 积分 -> 能量自检 -> 通过则继续 / 失败则重抽
# ---------------------------------------------------------------------------
SEED=""
OK=0

for attempt in $(seq 1 "$MAX_TRY"); do
    SEED=$(draw_seed)
    echo ""
    echo "------------------------------------------------------------"
    echo "[尝试 $attempt/$MAX_TRY]  随机种子 seed = $SEED"
    echo "------------------------------------------------------------"

    rm -rf "$OUTDIR"
    mkdir -p "$OUTDIR"
    cd "$OUTDIR"

    # --- 3.1 生成初始条件：mcluster ---------------------------------------
    #  -N 1000   粒子数 = 1000
    #  -R 1      半质量半径 = 1 pc
    #  -P 0      密度剖面 = Plummer (0=Plummer, 1=King, 2=Subr+2007, ...)
    #  -f 1      质量函数 = Kroupa (2001)
    #  -C 5      输出 NBODY6++GPU 格式（PeTar 可直接读取）
    #  -u 1      使用天文单位 (Msun, pc, km/s)
    #  -s <seed> 随机种子
    #  -o mc     输出文件名前缀
    echo "[1/6] mcluster 生成 Plummer 球初始条件 ..."
    if ! mcluster -N "$N_PARTICLE" -R "$R_HALF" -P 0 -f 1 -C 5 -u 1 \
                  -s "$SEED" -o mc > mc.log 2>&1; then
        echo "      mcluster 失败，重试"; continue
    fi
    INFILE=$(ls mc.dat.* 2>/dev/null | head -1)
    if [ -z "$INFILE" ]; then
        echo "      mcluster 未生成 mc.dat.*，重试"; continue
    fi

    # --- 3.2 转换成 PeTar 输入文件 ----------------------------------------
    #  mcluster 输出的速度单位是 km/s，PeTar 内部要求 pc/Myr，
    #  '-v kms2pcmyr' 完成单位换算。
    echo "[2/6] 转换为 PeTar 输入文件 input ..."
    if ! petar.init -v kms2pcmyr -f input "$INFILE" > init.log 2>&1; then
        echo "      petar.init 失败，重试"; continue
    fi

    # --- 3.3 选择最优二进制 ------------------------------------------------
    echo "[3/6] petar.select 选择最优二进制版本 ..."
    petar.select --optional mpi,omp,avx512,avx2 > /dev/null 2>&1

    # --- 3.4 树时间步 ------------------------------------------------------
    #  见文件开头 DT_SOFT 的说明：官方用 petar.find.dt 自动调优，但那依赖
    #  墙钟时间、不可复现，故这里固定取值。
    echo "[4/6] 树时间步 dt_soft = $DT_SOFT Myr"

    # --- 3.5 运行 N 体积分 -------------------------------------------------
    #  用 if ! 包裹：PeTar 遇到 SDAR 内部错误会 abort，此时应当重试而不是退出
    echo "[5/6] 运行 PeTar N 体模拟 (0 -> $T_END Myr) ..."
    if ! OMP_NUM_THREADS="$NTHREAD" OMP_STACKSIZE=128M \
         petar -u 1 -t "$T_END" -o "$T_OUT" -s "$DT_SOFT" input > output 2>&1; then
        echo "      ❌ petar 非正常退出（SDAR 内部报错，见上）"
        save_evidence "seed${SEED}_crash"
        continue
    fi

    # --- 3.6 能量守恒自检 --------------------------------------------------
    #  从日志提取每个输出时刻的 Total energy ($5) 与 Modify ($8，
    #  硬积分中粒子"修改"引起的累计能量变化)。SDAR 失稳时 Modify 会增长到
    #  与 |E| 同量级。
    #  判据（以全局能量守恒为准）：
    #    * 跑到 t_end
    #    * max|Modify| / |E0| < TOL_MODIFY
    #    * |dE/E| < TOL_DE
    #  hard_large_energy 诊断文件只作警告（官方文档：偶发的局部硬积分误差
    #  不影响整体动力学演化），但用一个宽松上限 MAX_DUMP 兜底。
    echo "      --- 能量守恒自检 ---"
    CHK_E0=0; CHK_E=0; CHK_T=-1; CHK_MOD=0
    eval "$(awk '
        /^Time:/{t=$2}
        /^Physic:/{
            if (!init) { e0=$5; init=1 }
            e=$5; m=$8; if (m<0) m=-m; if (m>mmax) mmax=m; et=t
        }
        END{
            if (init) printf "CHK_E0=%.6f CHK_E=%.6f CHK_T=%.3f CHK_MOD=%.6g\n",
                             e0, e, et, mmax
            else print "CHK_T=-1"
        }' output)"

    ndump=$(ls ${PREFIX}.hard_large_energy_* 2>/dev/null | wc -l)
    printf "      跑到 t=%.1f Myr；E: %.2f -> %.2f；max|Modify| = %.3g；dump 文件 %d 个\n" \
           "$CHK_T" "$CHK_E0" "$CHK_E" "$CHK_MOD" "$ndump"

    VERDICT=$(awk -v e0="$CHK_E0" -v e="$CHK_E" -v m="$CHK_MOD" -v t="$CHK_T" \
                  -v tend="$T_END" -v nd="$ndump" \
                  -v tde="$TOL_DE" -v tmod="$TOL_MODIFY" -v ndmax="$MAX_DUMP" 'BEGIN{
        if (e0==0) { print "REJECT|初始能量为零"; exit }
        de = (e-e0)/e0; if (de<0) de=-de
        mr = m/(e0<0?-e0:e0)
        if (t < tend-0.5) { printf "REJECT|未跑完 (t=%.1f < %.1f Myr)\n", t, tend; exit }
        if (nd > ndmax)   { printf "REJECT|hard_large_energy 诊断文件过多 (%d > %d)\n", nd, ndmax; exit }
        if (mr > tmod)    { printf "REJECT|max|Modify|/|E0| = %.3g > %.3g\n", mr, tmod; exit }
        if (de > tde)     { printf "REJECT||dE/E| = %.3g > %.3g\n", de, tde; exit }
        if (nd > 0) { printf "PASS|dE/E = %.3e，max|Modify|/|E0| = %.3e（%d 个偶发 hard_large_energy 诊断，属可接受范围）\n", de, mr, nd; exit }
        printf "PASS|dE/E = %.3e，max|Modify|/|E0| = %.3e，无诊断文件\n", de, mr
    }')

    case "$VERDICT" in
        PASS*)
            echo "      ✅ $VERDICT"
            OK=1
            break
            ;;
        *)
            echo "      ❌ $VERDICT  -> 该实现数值失稳，重抽种子"
            save_evidence "seed${SEED}_fail"
            ;;
    esac
done

if [ "$OK" -ne 1 ]; then
    echo ""
    echo "错误: $MAX_TRY 次尝试均未通过能量守恒自检。" >&2
    exit 1
fi

echo ""
echo "============================================================"
echo " 选定实现: seed = $SEED （已通过能量守恒自检）"
echo "============================================================"

# 记录种子，保证可复现
cat > SEED.txt <<EOF
seed            = $SEED
N               = $N_PARTICLE
R_half          = $R_HALF pc (Plummer)
t_end           = $T_END Myr
tree_dt         = $DT_SOFT Myr
OMP_NUM_THREADS = $NTHREAD
E(0)            = $CHK_E0
E(t_end)        = $CHK_E
|dE/E|          = $(awk -v a="$CHK_E0" -v b="$CHK_E" 'BEGIN{d=(b-a)/a; if(d<0)d=-d; printf "%.3e", d}')
max|Modify|     = $CHK_MOD
hard_large_energy 诊断文件数 = $ndump

复现命令:  SEED_OVERRIDE=$SEED bash run_cluster.sh
EOF
cat SEED.txt

# ---------------------------------------------------------------------------
# 4. 后处理
#    petar.data.gether : 合并各线程写出的分块数据
#    petar.data.process: 检测双星、计算拉格朗日半径与核心半径
#    '-G' 保持天文单位一致
# ---------------------------------------------------------------------------
echo ""
echo "[6/6] 后处理 (gether + process) ..."
petar.data.gether "$PREFIX" > gether.log 2>&1
petar.data.process -G "$G_AU" "${PREFIX}.snap.lst" > process.log 2>&1

# ---------------------------------------------------------------------------
# 5. 动画 + 拉格朗日半径演化
#    -i none          : petar 未开启 interrupt mode
#    -t none          : 未使用外部势场
#    -m x-y -R 20     : 左面板画 x-y 投影，视野 ±20 pc
#    -L data.lagr     : 右面板画拉格朗日半径演化
#    --rlagr-*        : 纵轴 0.01-100 pc，对数坐标（膨胀跨度大）
# ---------------------------------------------------------------------------
echo ""
echo "生成动画与拉格朗日半径演化图 ..."
petar.movie -i none -t none \
    -m x-y -R 20 \
    -L "${PREFIX}.lagr" --rlagr-min 0.01 --rlagr-max 1000 --rlagr-scale log \
    "${PREFIX}.snap.lst" > movie.log 2>&1

echo ""
echo "============================================================"
echo " 全部完成！主要产物："
echo "   ${OUTDIR}/movie.mp4        <- 星团动画 + 拉格朗日半径演化"
echo "   ${OUTDIR}/${PREFIX}.lagr   <- 拉格朗日半径数据"
echo "   ${OUTDIR}/${PREFIX}.*.png  <- 各时刻快照图"
echo "   ${OUTDIR}/SEED.txt         <- 使用的随机种子（可复现）"
echo "   ${OUTDIR}/output           <- PeTar 运行日志"
echo "============================================================"
