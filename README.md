# Formation and Evolution of Star Clusters — 作业与练习记录

本仓库用于存放 FESC（Formation and Evolution of Star Clusters）相关的作业、
练习与模拟结果，方便在**本地机器**和**集群**之间同步。

工作方式：

```bash
git pull            # 开始工作前先同步
# ... 做作业、跑模拟 ...
git add -A
git commit -m "描述这次做了什么"
git push            # 提交到 GitHub
```

> 注意：`FDPS` / `PeTar` / `SDAR` / `mcluster` 这 4 个软件包**不放在本仓库里**。
> 它们各自是独立的 git 仓库（有自己的 origin），用 `git clone` 单独管理。

---

## 目录

| 目录 | 内容 |
|------|------|
| `smoke/` | PeTar 冒烟测试（第一轮，极端初始条件，遇到已知保护性中止） |
| `smoke/simple/` | PeTar 冒烟测试（第二轮，简单初始条件，**跑通**） |
| `homework_1/` | 作业 1：疏散星团数值模拟（N=1000, Plummer, 200 Myr） |
| `homework_2/` | 作业 2：对比有无黑洞两类星团系统的演化（单质量 vs 含 10 M☉ 重星） |

---

## 一、PeTar 冒烟测试

### 目的

验证 PeTar 从安装到后处理的**整条链路**是否可用。

### 运行环境

| 项目 | 值 |
|------|-----|
| 集群 | Loong（登录节点 `mgt`，计算节点 `cn1-cn5` / `gn1`） |
| PeTar 版本 | `1759e`（分支 `experiment`） |
| SDAR 版本 | `447e`（分支 `experiment`） |
| FDPS 版本 | 8.0 |
| 编译配置 | `--with-mpi=no`（纯 OpenMP），`interrupt=off`，`external=off` |
| 二进制 | `petar.omp.avx512`（avx512 自动探测） |
| 单位制 | `-u 1`，即 Msun / pc / pc/Myr |

### 复现步骤

```bash
# 加载环境（PATH / PYTHONPATH / OMP_STACKSIZE / ulimit）
source <PeTar安装前缀>/petar_env.sh

cd smoke/simple

# 1) 生成初始条件（PeTar 仓库自带的确定性 IC 生成器）
python3 <PeTar>/test/validation/make_ic.py \
        --case functional_smoke --output initial.dat
#    → 1 对双星 (m=1.0, 0.8 Msun; a=0.01 pc; e=0.2) + 14 颗背景星，共 16 粒子
#    → 速度单位已是 pc/Myr，故下面 petar.init 不传 -v

# 2) 转成 PeTar 输入格式
petar.init -f data.init initial.dat

# 3) 运行求解器：t_end = 10 Myr，快照间隔 1
petar -u 1 -t 10 -o 1 -b 1 data.init >output 2>&1

# 4) 后处理
petar.data.process -G 0.00449830997959438 data.snap.lst
```

### 结果

| 项目 | 结果 |
|------|------|
| 退出码 | `0`（正常结束：`FDPS has successfully finished.`） |
| 耗时 | 约 20 秒 |
| 快照 | `data.0` ~ `data.10`（共 11 个） |
| 后处理产物 | `data.lagr`（拉格朗日半径）、`data.core`、`data.esc_single`、`data.esc_binary` |
| 硬积分器调用 | `AR_step_sum = 186.52`、`H4_step_sum = 67.92`（说明 SDAR 的 AR / Hermite 路径确实被用到） |
| 双星识别 | `Cluster = 2` |

### 已知问题

1. **第一轮（`smoke/`）失败**：使用 PeTar 自带的 `functional_mcluster_dual_merge_smoke`
   初始条件时，求解器在最初几步就中止（退出码 134）：

   ```
   Error! negative integrated time step persists after ds reductions in integrateToTime
         (time-transformation gauge issue, not ds size)
     n_negative_dt_accepted = 2
   ```

   该检查位于 `SDAR/src/AR/symplectic_integrator.h`，是**故意加入的保护性中止**
   （防止负时间步无限缩减）。SDAR 的 `lessons-learned.md` 记载该保护在其
   生产用初始条件上「0 次触发」，本轮的极端初始条件（双星周期 ~2e-8 Myr，
   半长轴 ~1e-8 pc ≈ 2 km）超出了它被验证过的范围。

2. 第二轮（`smoke/simple/`）虽然跑通，但留下了两个诊断转储文件：
   ```
   data.dump_large_step_h4_n11_g1_t4.281250_...
   data.hard_large_energy_h4_n11_g1_t4.250000_...
   ```
   说明 t ≈ 4.25–4.28 时 Hermite 积分器出现过大步长/能量误差，但运行自行恢复、未中断。

### 命令记录

每个运行目录下的 `commands.log` 记录了当时执行的确切命令行与时间戳。

---

## 二、作业 1：疏散星团数值模拟

详见 [`homework_1/README.md`](homework_1/README.md)。

### 要求与实现

| 项目 | 要求 | 实现 |
|------|------|------|
| 粒子数 | N = 1000 | 1000 |
| 密度剖面 | Plummer profile | `mcluster -P 0` |
| 演化时间 | 200 Myr | 200 Myr（每 1 Myr 一个快照，共 201 个） |
| 输出 | `petar.movie` 动画 + 拉格朗日半径演化 | `output/movie.mp4` |

### 运行环境

| 项目 | 值 |
|------|-----|
| 机器 | WSL Ubuntu 26.04，16 核 |
| 编译器 | gcc / gfortran 15.2 |
| PeTar 版本 | `1762e`（FDPS 8.0 + SDAR 447e） |
| 编译配置 | `--with-mpi=no`（纯 OpenMP），`interrupt=off`，`external=off` |
| 二进制 | `petar.omp.avx2` |
| 单位制 | `-u 1`，即 Msun / pc / pc/Myr |

### 结果

| 量 | 初始 | 200 Myr |
|------|------|------|
| 总质量 | 575.9 Msun | 575.9 Msun |
| 粒子数 | 1000 | 1000 |
| 总能量 E | −255.33 | −258.21（仍束缚） |
| 半质量半径 R50 | 1.05 pc | 7.60 pc（×7.3） |
| 核心半径 Rc | 0.477 pc | 0.280 pc（×0.6） |
| R90 | 3.26 pc | 304 pc |

即典型的孤立星团弛豫演化：**核心收缩 + 外层膨胀**同时发生，
200 Myr 内总能量始终为负，星团尚未瓦解，但已在明显朝瓦解（疏散）方向演化。

### 重要发现：随机实现高度敏感

这个 setup（N=1000 + Kroupa IMF 至 150 Msun 的重星 + 无恒星演化 + 无潮汐场）
会让星团核心深度坍缩并形成稳定的少体（多重）系统，触发 PeTar **SDAR 方法
对多重系统不保证经典能量守恒**的已知限制（官方 README *Significant Hard Energy*）。

失稳时系统会**凭空获得能量、由束缚变为不束缚，而程序不报错**——动画照样生成，
只有检查能量才能发现。实测 11 个随机实现中只有约 30% 可用：

| seed | 末态 E / 初态 E | max&#124;Modify&#124; | 判定 |
|------|------|------|------|
| 3 | −231.0 / −234.1 | 4e−3 | ✅ |
| 8 | −336.1 / −336.1 | 7e−6 | ✅ |
| 19046 | −258.2 / −255.3 | 0.068 | ✅ **采用** |
| 1 | +67.0 / −351.7 | 295 | ❌ |
| 5 | +837.2 / −274.6 | 650 | ❌ |
| 10 | +961.7 / −327.0 | 822 | ❌ |
| 2 | 崩溃（SDAR 报 gauge error 退出） | — | ❌ |

已用对照实验排除可控因素：OpenMP 线程数（1/4 都失稳）、树时间步（减小 2×/8× 无效）、
SDAR 步长（`--ar-ds-scale 0.1` 无效）、初始条件维里化（实测 E_kin/|E_pot| = 0.499）、
mcluster 参数（比对 `.info` 文件逐项一致）。

因此 `run_cluster.sh` 采用 **随机种子 + 能量守恒自检 + 自动重试**：
随机抽种子 → 跑 → 检查能量 → 不合格自动重抽重跑 → 合格后才做后处理与动画，
并把实际用到的种子写入 `output/SEED.txt`（因此仍可复现）。

### 目录内容

| 文件/目录 | 说明 |
|------|------|
| `run_cluster.sh` | 主脚本（含能量自检与自动重试） |
| `analyze.py` | 分析绘图脚本 |
| `README.md` | 完整实验报告 |
| `summary.txt` | 数值汇总 |
| `fig_lagr.png` | 拉格朗日半径演化图 |
| `fig_energy.png` | 能量守恒检验 |
| `fig_energy_compare.png` | 通过自检 vs 失稳实现的对比 |
| `output/` | 正式运行结果（快照、后处理、动画、日志、SEED.txt） |
| `evidence/`、`evidence_seed10_failed/` | 被拒绝的失稳实现的诊断证据 |
| `commands.log` | 命令行记录 |

> 说明：为控制仓库体积，`output/data.*.png`（201 张动画帧图，18 MB，
> 即 `movie.mp4` 的帧）未纳入版本控制。

### 命令记录

见 [`homework_1/commands.log`](homework_1/commands.log)。

---

## 三、作业 2：对比有无黑洞两类星团系统的演化

详见 [`homework_2/README.md`](homework_2/README.md)，
报告 PDF 在 [`homework_2/report/homework_2_report.pdf`](homework_2/report/homework_2_report.pdf)。

### 要求与实现

| 项目 | 要求（slide） | 实现 |
|------|------|------|
| 初始条件 | mcluster，Plummer，$R_h$ = 1 pc | `mcluster -N 1000 -P 0 -R 1 -t 0` |
| "无黑洞"组 | 单质量系统，全为 1 M☉ | `-f 0` → 1000 × 1 M☉，M = 1000 M☉ |
| "有黑洞"组 | 双质量系统 1 M☉/10 M☉，占比 1:9 与 5:5 | `-f 6 -m 1.0 -a 0.9 -m 10.0 -a 0.1` 与 `-a 0.5 / 0.5` |
| PeTar 设置 | 关闭恒星演化与外引力势 | 均为编译期关闭（本机构建不含外场） |
| 演化时间 | 2–3 $t_{rh}$ | $t_{end}=2\,t_{rh}$：A/B/C = 27.5 / 20 / 12 Myr |
| 输出 | $r_h$、$r_c$ 演化对比 + 质量迁移/核心塌缩判断 | `homework_2/figs/fig_rh.png`、`fig_rc.png` |
| 提交 | 报告 PDF（思路 + 参数 + 结果对比解说） | `homework_2/report/homework_2_report.pdf`（7 页，只含题目要求的三节） |

> 报告严格按 slide 的提交要求组织：**研究思路 / 模拟具体配置参数 / 结果绘图对比与解说**，
> 结果部分只画题目要求的两张图（$r_h$、$r_c$ 演化对比）。
> 其他诊断图（质量分层、平均质量剖面、能量守恒等）不进报告，
> 放在 `homework_2/figs/extra/` 备查。

### 结果（3 组配置 × 3 个随机实现的中位数）

| 量 | A 无黑洞 | B 重星 10% | C 重星 50% |
|---|---|---|---|
| $r_h$ 变化 | ×1.00 | ×1.10 | ×0.96 |
| $r_c$ 变化 | ×0.78 | ×**0.34** | ×0.84 |
| $r_c/r_h$ | 0.32 | **0.14** | 0.37 |
| $r_{h,\rm heavy}/r_{h,\rm light}$ | 1.000 | **0.27** | **0.53** |
| $\langle m\rangle(<r_h)/\bar m$ | 1.000 | **1.72** | **1.28** |

- **整体尺度不变**（$r_h$ 都在 1 pc 附近），但**内部发生了强烈的质量重分布**。
- **发生质量迁移**：重星子系统比轻星子系统紧凑 2–4 倍；无黑洞对照组的两个分层指标恒为 1.000。
- **发生核心收缩**：B 组 $r_c$ 掉到初值的 34%（三个实现分别是 0.156 / 0.031 / 0.172 pc，
  其中一个已进入深度塌缩）。
- **反直觉但有物理意义**：重星只占 10% 的 B 组比占 50% 的 C 组效应**更剧烈**——
  少量重星形成一个"小而密"的独立子系统时，内部弛豫极快，反而造成更极端的核心集中。

### 重要工程经验：SDAR 的两种失效模式

含 10 M☉ 重星的系统会形成成员数 20+ 的极紧多体系统，触发 PeTar SDAR 的已知限制
（官方 README *Significant Hard Energy*：慢化哈密顿量与经典哈密顿量不同，
**不保证物理能量守恒**）。实测到两种失效：

1. **能量凭空增加**——系统从束缚变成不束缚，但**程序不报错、正常退出、图照样出**。
   只能靠检查能量诊断（`Physic` 行的 `Modify` 列）发现。
2. **SDAR 时间变换退化（"卡死"）**——`|Int_err_cum/E|` 从 $10^{-9}$ 涨到 0.5 以上，
   `Step_count` 撞到 `--ar-max-nstep` 上限（$10^6$），模拟几乎不再前进，
   但**程序不报错、不退出，只是一直占着 CPU**（实测可持续几十分钟）。

因此 `run_all.sh` 采用 **随机种子 + 能量自检 + 看门狗 + 自动重试**：
- 能量自检：$|\Delta E/E| < 5\%$ 且 $\max|{\rm Modify}|/|E_0| < 1\%$
- 看门狗：`|Int_err_cum/E| > 0.1` 或 360 s 无新快照 → 杀掉重跑

### 一个关键教训：积分窗口要落在数值可靠区内

第一轮用 50 Myr（= 4–9 $t_{rh}$）的窗口跑，B 组 16 次尝试**只有 2 次通过（12.5%）**；
失败点几乎都在 $t \approx 17$–25 Myr（$\approx$ 1.7–2.5 $t_{rh}$）。
把窗口收到题目要求的 2 $t_{rh}$ 后成功率升到 33%，配合自动重试即可拿到 3 个干净实现；
A、C 组则一次通过。**结论：题目给的 2–3 $t_{rh}$ 窗口本身就是安全的，
跑更久反而会掉进数值失效区。**
第一阶段的尝试记录与失败日志保留在
`homework_2/attempts_phase1_50Myr.tsv` 与 `homework_2/evidence/_phase1_50Myr/`。

> 仓库体积：`homework_2/output/` 里的 ASCII 快照约 210 MB，已由 `.gitignore` 排除
> （可由 `bash run_all.sh` 重新生成）；仓库内保留 `data.lagr`、`data.core`、
> 运行日志与 `RUN_INFO.txt`（种子）等用于核对与复现的轻量产物。
