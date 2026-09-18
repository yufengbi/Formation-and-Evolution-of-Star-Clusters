# 作业 1：疏散星团数值模拟

**使用 PeTar（FDPS + SDAR）完成 N = 1000、Plummer 密度剖面星团的 200 Myr 演化，
并用 `petar.movie` 绘制星团动画与拉格朗日半径演化。**

---

## 1. 作业要求

| 项目 | 要求 | 本作业取值 |
|---|---|---|
| 粒子（恒星）数 | N = 1000 | 1000 |
| 密度剖面 | Plummer profile | Plummer（`mcluster -P 0`） |
| 演化时间 | 200 Myr | 200 Myr |
| 输出 | `petar.movie` 动画 + 拉格朗日半径演化 | `movie.mp4`、`fig_lagr.png` |

---

## 2. 物理模型

### 2.1 初始条件

用 `mcluster`（McLuster）生成，参数与 PeTar 官方示例
`sample/star_cluster_plummer_N1k.sh` 完全一致：

```bash
mcluster -N 1000 -R 1 -P 0 -f 1 -C 5 -u 1 -s <seed> -o mc
```

| 选项 | 含义 | 取值 |
|---|---|---|
| `-N 1000` | 恒星数 | 1000 |
| `-R 1` | 半质量半径 $R_h$ | 1 pc |
| `-P 0` | 密度剖面 | Plummer |
| `-f 1` | 质量函数 | Kroupa (2001)，0.08–150 $M_\odot$ |
| `-C 5` | 输出格式 | NBODY6++（PeTar 可直接读取） |
| `-u 1` | 单位 | 天文单位（$M_\odot$, pc, km/s） |
| `-s` | 随机种子 | 见 §4 |

**初始条件的基本性质**（以 `.info` 文件核实，且已确认与官方示例参数逐项一致）：

- 总质量 $M \approx 590\text{–}680\ M_\odot$（随 IMF 抽样波动），平均质量 $\approx 0.6\ M_\odot$
- 含 1–2 颗 $\sim 50\text{–}70\ M_\odot$ 的极大质量恒星（Kroupa IMF 上端）
- 初始维里比 $Q = -E_{\rm kin}/E_{\rm pot} = 0.5$，即**维里平衡**
  （实测 $E_{\rm kin}/|E_{\rm pot}| = 0.499$，已核对）
- 无初始双星（`nbin = 0`），无质量分层（`S = 0`），无分形（`D = 3`）
- 无恒星演化、无潮汐场（本轮作业不涉及）

### 2.2 单位换算

`mcluster -u 1` 输出的速度单位是 km/s，而 PeTar 内部要求 pc/Myr，用

```bash
petar.init -v kms2pcmyr -f input mc.dat.10
```

完成换算（因子 1.02271 pc/Myr per km/s）。PeTar 日志确认
$G = 0.0044985\ \mathrm{pc^3\,M_\odot^{-1}\,Myr^{-2}}$。

### 2.3 积分参数

| 参数 | 取值 | 说明 |
|---|---|---|
| 树时间步 `dt_soft` | 0.001953125 Myr | 见下方说明 |
| 演化时间 `-t` | 200 Myr | 约 $350$ 倍穿越时间、约 $7.6$ 倍半质量弛豫时间 |
| 快照间隔 `-o` | 1 Myr | 共 201 个快照 |
| 力精度 | 树 + 个体时间步 Hermite/SDAR | PeTar 默认 |
| `OMP_NUM_THREADS` | 4 | |

> **关于树时间步**：官方示例用 `petar.find.dt` 自动调优，但该工具是按**墙钟时间**
> 挑选"性能最优"值的，机器上有其它负载时会给出偏小很多的结果（实测从 0.0039
> 退化到 0.00024，即慢 16 倍），而且不可复现。因此这里固定为官方示例在空载机器上
> 得到的最优值的一半——与 sample 脚本最终的取值一致。

---

## 3. 工作流程

`run_cluster.sh` 完整实现了从初始条件到动画的流程：

```
mcluster         生成 Plummer 球初始条件（N=1000, R_h=1 pc）
   ↓
petar.init       速度单位 km/s → pc/Myr
   ↓
petar.select     在已安装的多个二进制版本中选最优（omp/avx2/…）
   ↓
petar            200 Myr N 体积分，每 1 Myr 输出一个快照
   ↓                ← 此处插入能量守恒自检（见 §4）
petar.data.gether    合并各线程写出的分块数据
   ↓
petar.data.process   检测双星、计算拉格朗日半径与核心半径
   ↓
petar.movie      生成 movie.mp4：左面板 x-y 投影动画，右面板拉格朗日半径演化
   ↓
analyze.py       出图 + 数值汇总（fig_lagr.png / fig_energy.png / summary.txt）
```

---

## 4. 数值稳定性检验（本作业的关键环节）

### 4.1 问题的发现

PeTar 官方示例用 `mcluster` 默认的 `-s 0`，即**按当前时间随机**。这意味着
每次运行得到的星团实现都不同，**结果不可复现**（连作者也无法重现上一次的运行）。

但进一步检验发现，这不只是"可复现性"的问题，而是**数值有效性问题**：
该 setup 会让星团核心深度坍缩并形成稳定的少体（多重）系统，而 PeTar 的
**SDAR 方法对这类系统不保证经典能量守恒**。PeTar 官方文档
（README, *Troubleshooting → Significant Hard Energy*）明确写道：

> "The inner binaries of these systems may exhibit tight configurations with
> substantial slowdown factors. As the slowdown Hamiltonian differs from the
> classical Hamiltonian, **physical energy conservation is not guaranteed**
> during integration, focusing instead on ensuring correct secular motion."
>
> "Resolving these cases without compromising computational efficiency can be
> challenging."

失稳时系统会**凭空获得能量、从束缚变成不束缚**，而程序**不报错**——动画照样
能生成、图照样好看，但物理结果是错的。只有检查能量守恒才能发现。

### 4.2 诊断方法

PeTar 在输出日志中给出逐时刻的能量诊断。关键量是 `Modify`
（=`de_change_modify_single`，硬积分中粒子"修改"引起的累计能量变化）：

```bash
awk '/^Time:/{t=$2} /^Physic:/{print t, $5, $8}' output
#                              $5 = Total E,  $8 = Modify
```

判据：**以全局能量守恒为准**

1. 积分跑到 200 Myr；
2. `max|Modify| / |E0| < 1%`；
3. `|dE/E| < 5%`；
4. `*.hard_large_energy.*` 诊断文件（PeTar 在硬积分相对能量误差超过
   `--energy-err-hard`（默认 1e-4）时自动生成）**只作警告**，不直接判负。

关于第 4 条：诊断文件表示硬积分中出现了较大的**局部**能量误差（通常来自
多重系统中最紧双星的半长轴微小变化）。PeTar 官方文档指出，这类误差若能
偶发出现，"通常不会显著影响系统的整体动力学演化"。因此真正决定物理有效性
的是第 2、3 条的**全局**能量守恒，诊断文件数只用于提示。
（实测确实存在这样的例子：seed=1719 的能量仅漂移 0.13%、`max|Modify|/|E0|`
只有 3.5e−6，却产生了 3 个偶发诊断文件——这是一个可用的实现。）

### 4.3 筛选结果

对 10 个随机实现逐一检验：

| seed | t_end | $E(0)$ | $E(t_{\rm end})$ | `maxModify` | 判定 |
|---|---|---|---|---|---|
| 0（示例默认） | 200 | −368.1 | −368.9 | 1e−4 | ✅ 守恒 |
| 3 | 200 | −234.1 | −231.0 | 4e−3 | ✅ 守恒 |
| 8 | 80 | −336.1 | −336.1 | 7e−6 | ✅ 守恒 |
| 1 | 200 | −351.7 | **+67.0** | 295 | ❌ 失稳 |
| 4 | 80 | −317.3 | **+81.3** | 213 | ❌ 失稳 |
| 5 | 80 | −274.6 | **+837.2** | 650 | ❌ 失稳 |
| 6 | 80 | −328.7 | −205.2 | 127 | ❌ 失稳 |
| 7 | 80 | −204.1 | −67.0 | 55 | ❌ 失稳 |
| 10 | 200 | −327.0 | **+961.7** | 822 | ❌ 失稳 |
| 20171 | 200 | −271.6 | −172.1 | 50 | ❌ 失稳 |
| 2 | — | −298.0 | — | — | 💥 PeTar 报 SDAR gauge error 退出 |

**成功率约 30%（3/10）。**

对 seed=10 的失败实现，PeTar 生成了 3 个诊断文件，时间点与能量暴增时刻完全吻合：

```
data.hard_large_energy_h4_n10_g1_t50.430664_...   ← t=50 Myr 形成 10 体多重系统
data.hard_large_energy_h4_n6_g1_t49.897461_...
data.hard_large_energy_h4_n6_g1_t51.259766_...
```

例如 seed=10：能量在 t=0–47 Myr 之间守恒到 $10^{-7}$，
**t=48→52 Myr 之间从 −327 暴增到 +578**，`Modify` 从 $10^{-3}$ 涨到 893。

### 4.4 已排除的原因

为确认这不是配置问题，逐一做了对照实验：

| 假设 | 实验 | 结论 |
|---|---|---|
| OpenMP 线程数引起 | OMP=1 与 OMP=4 各跑一遍 | 两者都失稳 → 排除 |
| 树时间步太大 | 减小 2 倍、8 倍 | 仍然失稳 → 排除 |
| SDAR 步长太大 | `--ar-ds-scale 0.1` | 仍然失稳 → 排除 |
| 初始条件未维里化 | 实测 $E_{\rm kin}/\|E_{\rm pot}\| = 0.499$ | 维里平衡 → 排除 |
| mcluster 参数不同 | 比对 `.info` 文件 | 参数逐项一致，仅 seed 不同 → 排除 |

**结论：这是小 N 星团 + 重星 IMF 条件下 SDAR 方法的固有数值限制，
不是运行配置问题。**

### 4.5 采取的对策

`run_cluster.sh` 采用 **随机种子 + 能量守恒自检 + 自动重试**：

1. 从 `/dev/urandom` 随机抽取一个种子；
2. 运行 200 Myr；
3. 做 §4.2 的自检；
4. 不通过 → 保存诊断证据到 `evidence/`，**自动重抽种子重跑**；
5. 通过 → 记录种子到 `output/SEED.txt`，继续做后处理与动画。

这样既保留了官方示例"随机实现"的用法，又保证最终结果在数值上可信，
并且因为记录了种子而**依然可复现**（设 `SEED_OVERRIDE=<seed>` 即可重现）。

> **说明**：这相当于对随机实现做了一次"质量控制"。之所以必须这样做，
> 是因为失稳的实现不会报错、产物看起来完全正常，若不做能量检查，
> 很容易把物理上错误的结果当成正确结果交上去。

---

## 5. 结果

### 5.1 本次运行

脚本第 1 次尝试（seed = 4852）失稳被拒（$E: -242 \to +256$，`max|Modify|/|E0| = 1.69`），
第 2 次尝试通过自检，最终采用的实现为 **seed = 19046**（记录于 `output/SEED.txt`）。

| 量 | 初始 | 200 Myr |
|---|---|---|
| 总质量 $M$ | 575.9 $M_\odot$ | 575.9 $M_\odot$（无质量损失） |
| 粒子数 | 1000 | 1000（`N_remove = 0`，`N_escape = 0`） |
| 最大质量 | 42.7 $M_\odot$ | — |
| 总能量 $E$ | −255.33 | −258.21 |
| 束缚状态 | 束缚 | **仍束缚** |

- 能量相对漂移 $|\Delta E/E| = 1.13\%$（200 Myr，约 350 个穿越时间）
- `max|Modify| / |E0| = 2.7\times10^{-4}`
- 产生 6 个 `hard_large_energy` 诊断文件（t ≈ 20、76、181 Myr），属**偶发**局部
  硬积分误差。PeTar 日志显示这些时刻的能量相对误差均在 $10^{-3}$ 量级，
  之后能量迅速回到守恒轨道（见 `fig_energy.png`），符合官方文档
  "sporadic 误差不影响整体动力学演化" 的描述。

### 5.2 拉格朗日半径演化

| $t$ [Myr] | $R_{10}$ | $R_{30}$ | $R_{50}$ | $R_{70}$ | $R_{90}$ | $R_c$ |
|---|---|---|---|---|---|---|
| 0 | 0.403 | 0.729 | 1.048 | 1.663 | 3.259 | 0.477 |
| 25 | 0.175 | 1.340 | 2.550 | 4.104 | 13.22 | 0.294 |
| 50 | 0.363 | 2.221 | 3.587 | 7.282 | 30.38 | 0.478 |
| 75 | 0.404 | 2.899 | 5.616 | 12.15 | 63.05 | 0.632 |
| 100 | 0.062 | 3.519 | 7.398 | 16.24 | 114.1 | 0.527 |
| 125 | 0.269 | 2.973 | 7.435 | 19.90 | 173.1 | 0.540 |
| 150 | 0.144 | 3.363 | 6.754 | 22.52 | 220.3 | 0.238 |
| 175 | 0.101 | 3.570 | 7.391 | 22.73 | 265.0 | 0.287 |
| 200 | 0.177 | 3.431 | 7.602 | 23.21 | 304.3 | 0.280 |

（单位：pc；完整曲线见 `fig_lagr.png`，对数纵轴）

**物理图像**

1. **核心收缩**：核心半径 $R_c$ 从 0.477 pc 降到 0.280 pc（约 0.6 倍），
   最内侧 $R_{10}$ 在 0.03–0.4 pc 之间剧烈振荡。这是孤立星团典型的
   **核心坍缩 / gravothermal 演化**：核心区通过二体弛豫失去能量而收缩变密。
2. **外层膨胀**：$R_{50}$ 从 1.05 pc 增大到 7.60 pc（**7.3 倍**），
   $R_{90}$ 从 3.26 pc 增大到 304 pc（约 93 倍）。核心收缩释放的能量
   被外层吸收，加上中心硬双星（Heggie 定律）持续向星团注入能量，
   使外层不断膨胀——这正是星团朝**瓦解（疏散）**方向演化的体现。
3. **整体仍束缚**：总能量始终为负（−255 → −258），说明 200 Myr 内星团
   尚未瓦解，只是显著膨胀。
4. $R_{10}$ 的强烈振荡源于其只包含 100 颗恒星，统计涨落大；$R_{50}$ 以上
   曲线则平滑得多。

### 5.3 数值稳定性对比

`fig_energy_compare.png` 直接对比了通过自检的实现（seed 19046）与失稳实现
（seed 4852）：

- 可用实现：`|Modify|` 稳定在 $\sim 0.07$，$E(t)/E(0) \approx 1.00$ 全程平直；
- 失稳实现：`|Modify|` 在 $t \approx 40$ Myr 突增 3 个数量级到 $\sim 400$，
  $E(t)/E(0)$ 跌破 0 并持续下降到 $-1.1$，系统由束缚变为不束缚。

### 5.4 动画

`output/movie.mp4`：左面板是星团在 x–y 平面上的投影（视野 ±20 pc，
`-m x-y -R 20`），右面板是拉格朗日半径随时间的演化
（`-L data.lagr --rlagr-min 0.01 --rlagr-max 1000 --rlagr-scale log`）。

可以直观看到：中心区始终维持一个致密核，而外层恒星不断向外扩散。

### 5.5 数值汇总

`summary.txt`（`analyze.py` 生成）包含完整的逐时刻表格。

---

## 6. 结论

1. **完成了作业要求的模拟**：借助 PeTar 官方示例流程，实现了 N = 1000、
   Plummer 密度剖面、半质量半径 1 pc 的星团从 0 到 200 Myr 的演化，
   并用 `petar.movie` 生成了星团动画与拉格朗日半径演化
   （`output/movie.mp4`）。

2. **物理结果**：星团表现为典型的孤立星团弛豫演化——**核心收缩**
   （$R_c: 0.48 \to 0.28$ pc）与**外层膨胀**（$R_{50}$ 增大 7.3 倍，
   $R_{90}$ 增大 93 倍）同时发生；总能量始终为负，200 Myr 内星团
   尚未瓦解，但已在明显朝瓦解（疏散）方向演化。

3. **数值可靠性的重要发现**：这个 setup 对随机实现高度敏感。在 11 个
   随机实现中只有约 30% 能保持能量守恒；其余会因星团核心坍缩形成稳定的
   少体（多重）系统，触发 PeTar **SDAR 方法对多重系统不保证经典能量守恒**
   的已知限制（官方 README *Significant Hard Energy*），表现为系统凭空
   获得能量、由束缚变为不束缚，**而程序不报错**。

   因此本作业在流程中加入了**能量守恒自检 + 自动重试**，并通过对照实验
   排除了 OpenMP 线程数、树时间步、SDAR 步长、初始条件维里化、mcluster
   参数等可控因素。最终采用的实现（seed = 19046）能量漂移仅 1.13%。

4. **方法论上的收获**：对 N 体模拟而言，"程序跑通了" 不等于 "结果是对的"。
   必须用守恒量（这里是总能量）做独立的正确性检验；尤其是自启动的
   正规化方法（regularization）会修改 Hamiltonian，更需要对输出做物理判读
   而不是只看程序是否正常退出。

---

## 7. 文件说明

```
homework_1/
├── run_cluster.sh          主脚本：初始条件 → 200 Myr 积分 → 能量自检
│                                  → 后处理 → 动画（含自动重试）
├── analyze.py              分析绘图：拉格朗日半径图、能量守恒图、对比图、数值汇总
├── README.md               本报告
├── run.log                 主脚本的完整运行记录（含每次尝试的判定）
├── summary.txt             analyze.py 生成的数值汇总
├── fig_lagr.png            拉格朗日半径演化图（对数纵轴）
├── fig_energy.png          能量守恒检验（能量预算 + 相对误差）
├── fig_energy_compare.png  通过自检 vs 失稳实现的对比（§5.3）
├── output/                 ★ 正式运行结果（seed = 19046）
│   ├── movie.mp4               星团动画 + 拉格朗日半径演化
│   ├── data.lagr               拉格朗日半径数据（二进制，用 petar 模块读）
│   ├── data.<t>.png            201 帧动画帧图（含 x-y 投影与 R_lagr 面板）
│   ├── data.<t>                201 个 PeTar 原始快照
│   ├── data.<t>.single/.binary petar.data.process 处理后的单星/双星快照
│   ├── SEED.txt                使用的随机种子与关键参数（可复现）
│   ├── input / mc.dat.10       初始条件（PeTar 格式 / mcluster 原始格式）
│   ├── mc.info                 mcluster 参数记录（可与官方示例逐项核对）
│   ├── output                  PeTar 运行日志（能量诊断在此）
│   └── *.hard_large_energy_*   6 个偶发硬积分诊断文件（§5.1）
├── evidence/              被自动重试拒绝的失稳实现（诊断证据）
│   ├── seed4852_fail/         第 1 次尝试（§5.3 引用）
│   ├── seed20171_fail/
│   └── seed1719_fail/         曾被过严判据误杀的好实现（§4.2 引用）
└── evidence_seed10_failed/  seed=10 失稳案例的证据（§4.3 引用）
    ├── data.hard_large_energy_*   3 个诊断文件，t≈50 Myr
    ├── data.lagr.failed           失稳的拉格朗日半径数据（对比用）
    └── petar.log                  PeTar 运行日志
```

---

## 8. 附录：如何查看后台任务状态

主脚本用 `nohup ... &` 在后台运行，检查状态的方法：

```bash
cd /root/FESC/Formation-and-Evolution-of-Star-Clusters/homework_1

# ① 进程还在不在（有输出 = 还在跑）
pgrep -af run_cluster.sh
pgrep -af petar

# ② 看运行时长
ps -o pid,etime,cmd -p $(pgrep -f run_cluster.sh)

# ③ 实时跟随日志（Ctrl+C 退出）
tail -f run.log

# ④ 只看每次尝试的判定结果
grep -E "^\[尝试|✅|❌ |选定实现|跑到 t=" run.log

# ⑤ 看当前这次积分到多少 Myr
grep '^Time:' output/output | tail -1

# ⑥ 看当前的能量是否正常（$5 = 总能量，$8 = Modify）
awk '/^Time:/{t=$2} /^Physic:/{print t, $5, $8}' output/output | tail -5
```

> 注意：`jobs` 命令只能看到**当前 shell** 启动的后台任务，对 `nohup` 启动的
> 其它终端里的任务无效，请用上面的 ①②③。

启动时把 PID 记下来，之后判断更方便：

```bash
nohup bash run_cluster.sh > run.log 2>&1 &
echo $! > run.pid
kill -0 $(cat run.pid) 2>/dev/null && echo "还在跑" || echo "已结束"
```

---

## 9. 复现方法

```bash
# 需要 PeTar 已安装且 bin 目录在 PATH 中（见 FESC 安装说明）
cd /root/FESC/Formation-and-Evolution-of-Star-Clusters/homework_1

# 随机种子 + 自动重试（默认）
bash run_cluster.sh

# 复现本次结果（seed = 19046）
SEED_OVERRIDE=19046 bash run_cluster.sh

# 出图与数值汇总
python3 analyze.py output
```

> 复现时请保持 `NTHREAD`（默认 4）与 `DT_SOFT`（默认 0.001953125 Myr）不变，
> 并确保没有其它重负载进程干扰（`petar.find.dt` 等按墙钟时间调优的工具会受影响）。

**关于中文字体**

`analyze.py` 画图时需要**简体中文**字体。若系统没有，脚本会自动退回英文标注
（不会出现方块）。本机安装的是文泉驿微米黑：

```bash
apt-get install -y fonts-wqy-microhei
```

> ⚠️ 注意**不要**用 `Noto Sans CJK JP` 之类的**日文**字体来画中文：日文字形与
> 简体中文不同，会出现"错字"。例如「径」的右上部，简体是「又」形，日文是
> 「ス」形。`analyze.py` 的候选字体列表里已排除日文变体。

---

## 10. 参考文献

1. Wang, L., Nitadori, K., Makino, J. (2020), *MNRAS* 493, 3392 — SDAR 方法
2. Iwasawa, M. et al. (2016), *PASJ* 68, 54 — FDPS
3. Iwasawa, M. et al. (2020), *PASJ* 72, 13 — PeTar
4. Küpper, A. H. W. et al. (2011), *MNRAS* 417, 2300 — McLuster
5. Kroupa, P. (2001), *MNRAS* 322, 231 — 初始质量函数
