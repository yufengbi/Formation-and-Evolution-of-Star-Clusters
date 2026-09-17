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
