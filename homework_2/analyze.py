#!/usr/bin/env python3
"""
作业2 分析脚本：对比"有无黑洞"两类星团系统的演化

输入：output/<config>/run<N>/ 下的
      data.snap.lst    快照文件清单
      data.<t>         ASCII 快照（-i 1 输出）
      data.lagr        Lagrangian 半径 / 各半径内的平均质量、天体数（peTar 后处理）
      data.core        密度中心与核心半径 Rc
      output           PeTar 运行日志（能量诊断）

输出：figs/*.png 与 summary.txt

用法:  python3 analyze.py [输出根目录, 默认 ./output]
"""

import os
import sys
import glob
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import petar
except ImportError:  # pragma: no cover
    petar = None

HERE = os.path.dirname(os.path.abspath(__file__))
OUTROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "output")
FIGDIR = os.path.join(HERE, "figs")
CACHE = os.path.join(HERE, "analysis_metrics.npz")

G_AU = 0.00449830997959438     # pc^3 / (Msun * Myr^2)
MASS_FRACTION = np.array([0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99])
M_HEAVY_THRESHOLD = 5.0        # 质量 > 5 Msun 归为"重星（黑洞）"

# ---------------------------------------------------------------------------
# 配置表
# ---------------------------------------------------------------------------
CONFIGS = [
    dict(key="A_single",  label="无黑洞 (单质量 1 M$_\\odot$)",       color="C0",
         title="单质量系统  N=1000 × 1 M$_\\odot$，M = 1000 M$_\\odot$"),
    dict(key="B_heavy10", label="黑洞组 10% (1/10 M$_\\odot$)",      color="C1",
         title="双质量系统  900 × 1 + 100 × 10 M$_\\odot$，M = 1900 M$_\\odot$"),
    dict(key="C_heavy50", label="黑洞组 50% (1/10 M$_\\odot$)",      color="C3",
         title="双质量系统  500 × 1 + 500 × 10 M$_\\odot$，M = 5500 M$_\\odot$"),
]

# ---------------------------------------------------------------------------
# 中文字体（只接受简体中文字体；找不到就退回英文，避免出现日文字形）
# ---------------------------------------------------------------------------
_CJK_CANDIDATES = (
    "WenQuanYi Micro Hei", "WenQuanYi Zen Hei",
    "Noto Sans CJK SC", "Source Han Sans SC", "Source Han Sans CN",
    "Microsoft YaHei", "SimHei", "PingFang SC", "Heiti SC", "AR PL UMing CN",
)
_CJK_FONT_FILES = (
    "/usr/share/fonts/truetype/wqy/wqy-microhei.ttc",
    "/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc",
)


def _setup_cjk_font():
    from matplotlib import font_manager
    for path in _CJK_FONT_FILES:
        if os.path.exists(path):
            try:
                font_manager.fontManager.addfont(path)
            except Exception:
                pass
    for name in _CJK_CANDIDATES:
        try:
            font_manager.findfont(name, fallback_to_default=False)
        except Exception:
            continue
        plt.rcParams["font.sans-serif"] = [name] + list(plt.rcParams["font.sans-serif"])
        plt.rcParams["axes.unicode_minus"] = False
        return name
    return None


CJK_FONT = _setup_cjk_font()
CJK = CJK_FONT is not None
print(f"[字体] {CJK_FONT if CJK else '未找到简体中文字体 -> 图注用英文'}")


def L(zh, en):
    return zh if CJK else en


# ===========================================================================
# 1. 读取原始数据
# ===========================================================================
def read_energy_log(path):
    """从 PeTar 日志提取逐时刻能量诊断。

    格式:  Physic: <Error/Total> <Error> <Error_cum> <Total E> <Kinetic>
                   <Potential> <Modify> ...
    """
    t, etot, ekin, epot, eerr, emod = [], [], [], [], [], []
    cur_t = None
    if not os.path.isfile(path):
        return tuple(np.array([]) for _ in range(6))
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith("Time:"):
                cur_t = float(line.split()[1])
            elif line.startswith("Physic:") and cur_t is not None:
                v = line.split()[1:]
                if len(v) >= 7:
                    eerr.append(float(v[0]))
                    etot.append(float(v[3]))
                    ekin.append(float(v[4]))
                    epot.append(float(v[5]))
                    emod.append(float(v[6]))
                    t.append(cur_t)
    return (np.array(t), np.array(etot), np.array(ekin),
            np.array(epot), np.array(eerr), np.array(emod))


def read_snapshot(path):
    """读取 ASCII 快照。

    列: mass, x, y, z, vx, vy, vz, bin_stat, r_search, id, mass_bk, status,
        r_in, r_out, acc_x, acc_y, acc_z, pot_tot, pot_soft, n_b
    """
    d = np.loadtxt(path, skiprows=1)
    if d.ndim == 1:
        d = d[None, :]
    return dict(mass=d[:, 0], pos=d[:, 1:4], vel=d[:, 4:7], pot=d[:, 17])


def snapshot_list(run_dir):
    """返回按时间排序的 (t, 文件路径) 列表"""
    lst = os.path.join(run_dir, "data.snap.lst")
    files = []
    if os.path.isfile(lst):
        with open(lst) as fh:
            for line in fh:
                fn = line.strip()
                if not fn:
                    continue
                p = fn if os.path.isabs(fn) else os.path.join(run_dir, fn)
                if os.path.isfile(p):
                    files.append(p)
    else:
        cand = [p for p in glob.glob(os.path.join(run_dir, "data.?"))
                + glob.glob(os.path.join(run_dir, "data.?.?"))
                if os.path.isfile(p)]
        files = cand
    out = []
    for p in files:
        try:
            with open(p) as fh:
                head = fh.readline().split()
            out.append((float(head[2]), p))
        except Exception:
            continue
    out.sort(key=lambda x: x[0])
    return out


def read_lagr(run_dir):
    """读取 peTar 后处理给出的 Lagrangian 数据（按时间排序）"""
    if petar is None:
        return None
    p = os.path.join(run_dir, "data.lagr")
    if not os.path.isfile(p):
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lag = petar.LagrangianMultiple(
            mass_fraction=MASS_FRACTION, calc_energy=True)
        lag.fromfile(p)
    o = np.argsort(lag.time)
    return dict(time=lag.time[o], r=lag.all.r[o], m=lag.all.m[o],
                n=lag.all.n[o], vr=lag.all.vr[o],
                mass_fraction=MASS_FRACTION)


def read_core(run_dir):
    """读取密度中心位置与核心半径 Rc（按时间排序）"""
    if petar is None:
        return None
    p = os.path.join(run_dir, "data.core")
    if not os.path.isfile(p):
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        core = petar.Core()
        core.fromfile(p)
    o = np.argsort(core.time)
    return dict(time=core.time[o], pos=core.pos[o], rc=core.rc[o])


def enclosed_radius(r, m, frac, center_zero=True):
    """给定相对中心的距离 r 与质量 m，返回累计质量达到 frac*M 的半径"""
    o = np.argsort(r)
    cs = np.cumsum(m[o])
    target = frac * m.sum()
    if cs[-1] < target:
        return np.nan
    return float(np.interp(target, cs, r[o]))


def snapshot_metrics(snap, center):
    """从单张快照计算半径/分层指标"""
    m = snap["mass"]
    d = snap["pos"] - center
    r = np.sqrt((d * d).sum(axis=1))
    heavy = m > M_HEAVY_THRESHOLD
    res = dict(
        M=m.sum(), N=len(m),
        rh=enclosed_radius(r, m, 0.5),
        r10=enclosed_radius(r, m, 0.10),
        r90=enclosed_radius(r, m, 0.90),
        mbar=m.mean(),
    )
    # 重/轻子系统的半质量半径
    #  若某一组为空（例如单质量系统里没有 >5 Msun 的重星），则两组的
    #  半质量半径都取整体值 —— 此时"重星/轻星"本就没有区别，比值恒为 1，
    #  这正是等质量对照组的正确极限，也便于在同一张图里画出参照线。
    for nm, sel in (("heavy", heavy), ("light", ~heavy)):
        if sel.sum() > 0:
            res[f"rh_{nm}"] = enclosed_radius(r[sel], m[sel], 0.5)
            res[f"n_{nm}"] = int(sel.sum())
            res[f"mbar_{nm}"] = float(m[sel].mean())
        else:
            res[f"rh_{nm}"] = res["rh"]
            res[f"n_{nm}"] = 0
            res[f"mbar_{nm}"] = np.nan
    # 半质量半径以内的平均质量（直接的质量分层度量）
    if np.isfinite(res["rh"]):
        sel = r <= res["rh"]
        res["mbar_rh"] = float(m[sel].mean())
        res["f_heavy_rh"] = float(heavy[sel].mean())
    else:
        res["mbar_rh"] = np.nan
        res["f_heavy_rh"] = np.nan
    return res


# ===========================================================================
# 2. 汇总一个 run 的所有时刻
# ===========================================================================
METRIC_KEYS = ("rh", "r10", "r90", "mbar", "rh_heavy", "rh_light",
               "mbar_rh", "f_heavy_rh", "n_heavy", "n_light", "M", "N")


def analyse_run(cfg_key, run_idx):
    run_dir = os.path.join(OUTROOT, cfg_key, f"run{run_idx}")
    if not os.path.isdir(run_dir):
        return None

    core = read_core(run_dir)
    snaps = snapshot_list(run_dir)
    if not snaps:
        return None

    # 密度中心随时间的插值函数（快照时刻与 data.core 时刻一致）
    if core is not None:
        ct, cpos = core["time"], core["pos"]
    else:
        ct, cpos = np.array([snaps[0][0]]), np.zeros((1, 3))

    times, metrics = [], {k: [] for k in METRIC_KEYS}
    for t, path in snaps:
        try:
            snap = read_snapshot(path)
        except Exception:
            continue
        center = np.array([np.interp(t, ct, cpos[:, i]) for i in range(3)])
        mt = snapshot_metrics(snap, center)
        times.append(t)
        for k in METRIC_KEYS:
            metrics[k].append(mt.get(k, np.nan))

    out = dict(cfg=cfg_key, run=run_idx, dir=run_dir,
               time=np.array(times),
               **{k: np.array(v, dtype=float) for k, v in metrics.items()})
    if core is not None:
        out["core_time"] = core["time"]
        out["rc"] = core["rc"]
    if petar is not None:
        out["lagr"] = read_lagr(run_dir)
    e = read_energy_log(os.path.join(run_dir, "output"))
    out["energy"] = dict(zip(("t", "etot", "ekin", "epot", "eerr", "emod"), e))
    return out


def collect_runs(cfg):
    runs = []
    for d in sorted(glob.glob(os.path.join(OUTROOT, cfg["key"], "run*"))):
        m = os.path.basename(d)
        if not m.startswith("run") or not m[3:].isdigit():
            continue
        r = analyse_run(cfg["key"], int(m[3:]))
        if r is not None:
            runs.append(r)
    return runs


# ===========================================================================
# 3. 出图
# ===========================================================================
def _save(fig, name, extra=False):
    """保存图。extra=True 的图放入 figs/extra/（题目只要 r_h 与 r_c 两张主图）"""
    d = os.path.join(FIGDIR, "extra") if extra else FIGDIR
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, name)
    fig.tight_layout()
    fig.savefig(p, dpi=160)
    plt.close(fig)
    print(f"[写入] {p}")


def _ref_grid(runs, key="time"):
    """公共时间网格：取快照最多的那个实现的时刻序列"""
    if not runs:
        return np.array([])
    return max((r[key] for r in runs if key in r), key=len)


def align(runs, key, tref, tkey="time"):
    """把各实现的曲线插值到公共时间网格上（超出该实现时间范围的填 NaN）"""
    out = []
    for r in runs:
        t = r.get(tkey)
        v = r.get(key)
        if t is None or v is None or len(t) == 0:
            continue
        vi = np.interp(tref, t, v)
        vi[tref > t[-1] + 1e-9] = np.nan
        vi[tref < t[0] - 1e-9] = np.nan
        out.append(vi)
    return np.vstack(out) if out else None


def value_at(run, key, t, tkey="time"):
    """取 run 在时刻 t 的值（线性插值；超出范围返回 nan）"""
    tt, vv = run.get(tkey), run.get(key)
    if tt is None or vv is None or len(tt) == 0 or t > tt[-1] + 1e-9:
        return np.nan
    return float(np.interp(t, tt, vv))


def _band(ax, runs, key, color, label, alpha=0.18, lw=2.0, ls="-", xscale=None):
    """把所有实现的曲线画成细线 + 中位数粗线 + 1σ 带

    xscale: 可选，把时间坐标做变换（例如换成归一化时间 t/t_rh）
    """
    tref = _ref_grid(runs)
    v = align(runs, key, tref)
    if v is None:
        return
    x = xscale(tref) if xscale else tref
    med = np.nanmedian(v, axis=0)
    lo = np.nanpercentile(v, 16, axis=0)
    hi = np.nanpercentile(v, 84, axis=0)
    for row in v:
        ax.plot(x, row, color=color, lw=0.7, alpha=0.45, zorder=2)
    ax.fill_between(x, lo, hi, color=color, alpha=alpha, lw=0, zorder=1)
    ax.plot(x, med, color=color, lw=lw, ls=ls, label=label, zorder=3)


def _xhdr():
    return L(r"$t\,/\,t_{rh}$", "$t/t_{rh}$")


def _norm(t, cfg, t_rh):
    """把时刻换算成归一化弛豫时间 t/t_rh"""
    trh = t_rh.get(cfg["key"]) if t_rh else None
    return t / trh if trh else t


def _mark_trh(ax):
    for xv, ls in ((1.0, ":"), (2.0, "--")):
        ax.axvline(xv, color="0.45", ls=ls, lw=1.0)


def plot_rh(all_runs, t_rh, extra=False):
    """半质量半径 r_h 演化对比

    主图（extra=False）横轴用绝对时间 t [Myr]；
    归一化版本（extra=True）横轴用 t/t_rh，存入 figs/extra/。
    """
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    for cfg in CONFIGS:
        runs = all_runs.get(cfg["key"], [])
        if extra:
            _band(ax, runs, "rh", cfg["color"], cfg["label"],
                  xscale=lambda t, c=cfg: _norm(t, c, t_rh))
        else:
            _band(ax, runs, "rh", cfg["color"], cfg["label"])
    if extra:
        ax.set_xlabel(_xhdr(), fontsize=12)
        _mark_trh(ax)
    else:
        ax.set_xlabel("$t$  [Myr]", fontsize=12)
        for cfg in CONFIGS:
            trh = t_rh.get(cfg["key"])
            if trh:
                ax.axvline(2 * trh, color=cfg["color"], ls="--", lw=1.0, alpha=0.6)
    ax.set_ylabel("$r_h$  [pc]", fontsize=12)
    ax.set_title(L("半质量半径 $r_h$ 的演化", "Half-mass radius $r_h$"), fontsize=12)
    ax.grid(alpha=0.3, ls=":")
    ax.legend(fontsize=10)
    _save(fig, "fig_rh_normalised.png" if extra else "fig_rh.png", extra)


def plot_rc(all_runs, t_rh, extra=False):
    """核心半径 r_c 演化对比；横轴同 plot_rh"""
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    for cfg in CONFIGS:
        runs = all_runs.get(cfg["key"], [])
        tref = _ref_grid(runs, "core_time")
        v = align(runs, "rc", tref, tkey="core_time")
        if v is None:
            continue
        x = _norm(tref, cfg, t_rh) if extra else tref
        for row in v:
            ax.plot(x, row, color=cfg["color"], lw=0.7, alpha=0.45)
        ax.fill_between(x, np.nanpercentile(v, 16, axis=0),
                        np.nanpercentile(v, 84, axis=0),
                        color=cfg["color"], alpha=0.18, lw=0)
        ax.plot(x, np.nanmedian(v, axis=0), color=cfg["color"], lw=2.0,
                label=cfg["label"])
        if not extra:
            trh = t_rh.get(cfg["key"])
            if trh:
                ax.axvline(2 * trh, color=cfg["color"], ls="--", lw=1.0, alpha=0.6)
    if extra:
        ax.set_xlabel(_xhdr(), fontsize=12)
        _mark_trh(ax)
    else:
        ax.set_xlabel("$t$  [Myr]", fontsize=12)
    ax.set_ylabel("$r_c$  [pc]", fontsize=12)
    ax.set_title(L("核心半径 $r_c$ 的演化", "Core radius $r_c$"), fontsize=12)
    ax.grid(alpha=0.3, ls=":")
    ax.legend(fontsize=10)
    _save(fig, "fig_rc_normalised.png" if extra else "fig_rc.png", extra)


def plot_norm(all_runs, extra=False):
    """归一化演化 + 核心集中度 r_c/r_h"""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

    for cfg in CONFIGS:
        runs = all_runs.get(cfg["key"], [])
        if not runs:
            continue
        for i, r in enumerate(runs):
            t = r["time"]
            axes[0].plot(t, r["rh"] / r["rh"][0], color=cfg["color"],
                         lw=0.8, alpha=0.5,
                         label=cfg["label"] if i == 0 else None)
            if "rc" in r:
                axes[1].plot(r["core_time"], r["rc"] / r["rh"][0],
                             color=cfg["color"], lw=0.8, alpha=0.5)

    axes[0].set_xlabel("$t$ [Myr]"); axes[0].set_ylabel("$r_h(t)/r_h(0)$")
    axes[0].set_title(L("半质量半径（归一化）", "Normalised $r_h$"))
    axes[1].set_xlabel("$t$ [Myr]"); axes[1].set_ylabel("$r_c(t)/r_h(0)$")
    axes[1].set_title(L("核心半径（以初始 $r_h$ 为单位）", "Core radius"))
    for ax in axes:
        ax.grid(alpha=0.3, ls=":")
    axes[0].legend(fontsize=9)
    _save(fig, "fig_rh_rc_normalised.png", extra)


def plot_segregation(all_runs, t_rh, extra=False):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

    for cfg in CONFIGS:
        runs = all_runs.get(cfg["key"], [])
        if not runs:
            continue
        for i, r in enumerate(runs):
            t = _norm(r["time"], cfg, t_rh)
            lbl = cfg["label"] if i == 0 else None
            # 重/轻子系统的半质量半径之比（<1 表示重星更向中心集中）
            axes[0].plot(t, r["rh_heavy"] / r["rh_light"], color=cfg["color"],
                         lw=1.4, alpha=0.8, label=lbl)
            # 半质量半径以内的平均质量 / 整体平均质量
            axes[1].plot(t, r["mbar_rh"] / r["mbar"], color=cfg["color"],
                         lw=1.4, alpha=0.8, label=lbl)

    axes[0].axhline(1.0, color="k", lw=0.9, ls="--")
    axes[0].set_ylabel("$r_{h,\\rm heavy} / r_{h,\\rm light}$")
    axes[0].set_title(L("重星与轻星的半质量半径之比  (<1 = 重星向心集中)",
                        "Mass segregation ratio"))
    axes[1].axhline(1.0, color="k", lw=0.9, ls="--")
    axes[1].set_ylabel("$\\langle m\\rangle(<r_h) / \\bar m$")
    axes[1].set_title(L("$r_h$ 以内平均质量 / 整体平均质量", "Mean mass inside $r_h$"))
    for ax in axes:
        ax.set_xlabel(_xhdr())
        ax.grid(alpha=0.3, ls=":")
        ax.legend(fontsize=9)
    _save(fig, "fig_segregation.png", extra)


def plot_massprofile(all_runs, extra=False):
    """不同时刻：累计平均质量 <m>(<r) 沿半径的分布（质量分层剖面）"""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), sharey=True)
    for ax, cfg in zip(axes, [CONFIGS[1], CONFIGS[2]]):
        runs = all_runs.get(cfg["key"], [])
        if not runs:
            continue
        run = runs[0]
        core = read_core(run["dir"])
        snaps = snapshot_list(run["dir"])
        times = [0.0, run["time"][-1] / 2.0, run["time"][-1]]
        cols = plt.cm.viridis(np.linspace(0.15, 0.85, len(times)))
        for tt, col in zip(times, cols):
            i = int(np.argmin(np.abs(run["time"] - tt)))
            path = min(snaps, key=lambda s: abs(s[0] - run["time"][i]))[1]
            snap = read_snapshot(path)
            c = np.array([np.interp(snaps[i][0], core["time"], core["pos"][:, k])
                          for k in range(3)])
            d = snap["pos"] - c
            r = np.sqrt((d * d).sum(axis=1))
            o = np.argsort(r)
            rr, mm = r[o], snap["mass"][o]
            cum_m = np.cumsum(mm)
            cum_n = np.arange(1, len(mm) + 1)
            ax.plot(rr, cum_m / cum_n, color=col, lw=1.6,
                    label=f"$t$ = {run['time'][i]:.0f} Myr")
        ax.axhline(1.0, color="k", lw=0.8, ls="--")
        ax.set_xscale("log")
        ax.set_xlim(0.01, 30)
        ax.set_xlabel(L("距中心的半径 $r$ [pc]", "radius [pc]"))
        ax.set_title(L("累计平均质量 $\\langle m\\rangle(<r)$", "Mean mass profile")
                     + f"  —  {cfg['key']}")
        ax.grid(alpha=0.3, which="both", ls=":")
        ax.legend(fontsize=9)
    axes[0].set_ylabel("$\\langle m\\rangle(<r)$  [M$_\\odot$]")
    _save(fig, "fig_massprofile.png", extra)


def plot_energy(all_runs, extra=False):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for cfg in CONFIGS:
        for r in all_runs.get(cfg["key"], []):
            e = r["energy"]
            if e["t"].size == 0:
                continue
            de = np.abs((e["etot"] - e["etot"][0]) / e["etot"][0])
            axes[0].plot(e["t"], de, color=cfg["color"], lw=1.2, alpha=0.85)
            axes[1].plot(e["t"], np.abs(e["emod"]) / abs(e["etot"][0]),
                         color=cfg["color"], lw=1.2, alpha=0.85)
    axes[0].axhline(0.05, color="k", ls="--", lw=0.9)
    axes[0].text(0.02, 0.055, "5% tol.", transform=axes[0].get_yaxis_transform(),
                 fontsize=8)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("$t$ [Myr]"); axes[0].set_ylabel("$|\\Delta E/E|$")
    axes[0].set_title(L("全局能量守恒", "Global energy conservation"))
    axes[1].axhline(0.01, color="k", ls="--", lw=0.9)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("$t$ [Myr]"); axes[1].set_ylabel("$\\max|{\\rm Modify}|/|E_0|$")
    axes[1].set_title(L("SDAR 累计能量修正", "Cumulative SDAR correction"))
    for ax in axes:
        ax.grid(alpha=0.3, which="both", ls=":")
    _save(fig, "fig_energy.png", extra)


def plot_panel(all_runs, t_rh, extra=False):
    """报告用 2x2 汇总面板（横轴统一为归一化弛豫时间 t/t_rh）"""
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))

    ax = axes[0, 0]
    for cfg in CONFIGS:
        _band(ax, all_runs.get(cfg["key"], []), "rh", cfg["color"],
              cfg["label"], lw=1.8, xscale=lambda t, c=cfg: _norm(t, c, t_rh))
    ax.set_ylabel("$r_h$ [pc]")
    ax.set_title(L("(a) 半质量半径 $r_h$", "(a) half-mass radius"))
    ax.grid(alpha=0.3, ls=":"); ax.legend(fontsize=9)

    ax = axes[0, 1]
    for cfg in CONFIGS:
        runs = all_runs.get(cfg["key"], [])
        for i, r in enumerate(runs):
            if "rc" in r:
                ax.plot(_norm(r["core_time"], cfg, t_rh), r["rc"],
                        color=cfg["color"], lw=1.3, alpha=0.8,
                        label=cfg["label"] if i == 0 else None)
    ax.set_ylabel("$r_c$ [pc]")
    ax.set_title(L("(b) 核心半径 $r_c$", "(b) core radius"))
    ax.grid(alpha=0.3, ls=":"); ax.legend(fontsize=9)

    ax = axes[1, 0]
    for cfg in CONFIGS:
        for i, r in enumerate(all_runs.get(cfg["key"], [])):
            ax.plot(_norm(r["time"], cfg, t_rh),
                    r["rh_heavy"] / r["rh_light"], color=cfg["color"],
                    lw=1.3, alpha=0.8, label=cfg["label"] if i == 0 else None)
    ax.axhline(1.0, color="k", ls="--", lw=0.9)
    ax.set_ylabel("$r_{h,\\rm heavy}/r_{h,\\rm light}$")
    ax.set_title(L("(c) 质量分层：重/轻半质量半径之比", "(c) mass segregation"))
    ax.grid(alpha=0.3, ls=":"); ax.legend(fontsize=9)

    ax = axes[1, 1]
    for cfg in CONFIGS:
        for i, r in enumerate(all_runs.get(cfg["key"], [])):
            ax.plot(_norm(r["time"], cfg, t_rh), r["mbar_rh"] / r["mbar"],
                    color=cfg["color"], lw=1.3, alpha=0.8,
                    label=cfg["label"] if i == 0 else None)
    ax.axhline(1.0, color="k", ls="--", lw=0.9)
    ax.set_ylabel("$\\langle m\\rangle(<r_h)/\\bar m$")
    ax.set_title(L("(d) $r_h$ 以内平均质量 / 整体平均质量", "(d) mean mass inside $r_h$"))
    ax.grid(alpha=0.3, ls=":"); ax.legend(fontsize=9)

    _save(fig, "fig_panel.png", extra)


# ===========================================================================
# 4. 数值汇总
# ===========================================================================
def relaxation_times(all_runs):
    """各系统的 t_rh（Spitzer 半质量弛豫时间）

    t_rh = 0.138 sqrt(N R_h^3 /(G mbar)) / ln(0.11 N)

    优先用 run_all.sh 在**初始条件**上算出的值（存在 features.txt 里），
    保证与自动确定的 t_end 一致；缺失时退回用 t=0 快照的实测值估算。
    """
    out = {}
    for cfg in CONFIGS:
        runs = all_runs.get(cfg["key"], [])
        if not runs:
            continue
        val = np.nan
        feat = os.path.join(runs[0]["dir"], "features.txt")
        if os.path.isfile(feat):
            try:
                val = float(open(feat).read().split()[4])
            except Exception:
                val = np.nan
        if not np.isfinite(val):
            r = runs[0]
            N, mbar, rh = r["N"][0], r["mbar"][0], r["rh"][0]
            val = 0.138 * np.sqrt(N * rh ** 3 / (G_AU * mbar)) / np.log(0.11 * N)
        out[cfg["key"]] = val
    return out


def write_summary(all_runs, t_rh):
    lines = []
    A = "=" * 74
    lines += [A, "作业2  有无黑洞两类星团系统演化对比  ——  数值汇总", A, ""]
    lines += ["参数：N = 1000，Plummer 剖面，R_h = 1 pc，维里平衡 Q = 0.5，",
              "      无潮汐外场，无恒星演化，PeTar 树 + SDAR 个体时间步积分", ""]

    def med_at(cfg_key, key, t, tkey="time"):
        """各实现插值到时刻 t 后取中位数"""
        vals = [value_at(r, key, t, tkey) for r in all_runs.get(cfg_key, [])]
        vals = [v for v in vals if np.isfinite(v)]
        return float(np.nanmedian(vals)) if vals else np.nan

    lines += ["-" * 74, "1. 各系统的弛豫时间与积分时长", "-" * 74]
    lines += [f"  {'配置':<12}{'M [Msun]':>10}{'mbar':>8}{'t_rh [Myr]':>12}"
              f"{'t_end [Myr]':>13}{'t_end/t_rh':>12}"]
    for cfg in CONFIGS:
        runs = all_runs.get(cfg["key"], [])
        if not runs:
            continue
        r = runs[0]
        trh = t_rh.get(cfg["key"], np.nan)
        tend = max(x["time"][-1] for x in runs)
        lines.append(f"  {cfg['key']:<12}{r['M'][0]:>10.0f}{r['mbar'][0]:>8.2f}"
                     f"{trh:>12.2f}{tend:>13.1f}{tend / trh:>12.1f}")
    lines.append("")

    lines += ["-" * 74, "2. 半质量半径 r_h 与核心半径 r_c 的演化（中位数）", "-" * 74]
    lines += [f"  {'配置':<12}{'r_h(0)':>9}{'r_h(end)':>10}{'倍数':>8}"
              f"{'r_c(0)':>9}{'r_c(end)':>10}{'倍数':>8}{'r_c/r_h(end)':>13}"]
    for cfg in CONFIGS:
        runs = all_runs.get(cfg["key"], [])
        if not runs:
            continue
        key = cfg["key"]
        tend = max(x["time"][-1] for x in runs)
        rh0 = med_at(key, "rh", 0.0)
        rhe = med_at(key, "rh", tend)
        rc0 = med_at(key, "rc", 0.0, "core_time")
        rce = med_at(key, "rc", tend, "core_time")
        lines.append(f"  {key:<12}{rh0:>9.3f}{rhe:>10.3f}{rhe / rh0:>8.2f}"
                     f"{rc0:>9.3f}{rce:>10.3f}{rce / rc0:>8.2f}{rce / rhe:>13.3f}")
    lines.append("")

    lines += ["-" * 74, "3. 质量分层（质量迁移）诊断", "-" * 74]
    lines += ["  r_h,heavy/r_h,light < 1 表示重星比轻星更向中心集中（发生质量分层）",
              f"  {'配置':<12}{'比值(0)':>10}{'比值(end)':>12}"
              f"{'<m>(<r_h)/mbar (0)':>21}{'(end)':>9}"]
    for cfg in CONFIGS:
        runs = all_runs.get(cfg["key"], [])
        if not runs:
            continue
        key = cfg["key"]
        tend = max(x["time"][-1] for x in runs)
        ratio0 = med_at(key, "rh_heavy", 0.0) / med_at(key, "rh_light", 0.0)
        ratioe = med_at(key, "rh_heavy", tend) / med_at(key, "rh_light", tend)
        mr0 = med_at(key, "mbar_rh", 0.0) / med_at(key, "mbar", 0.0)
        mre = med_at(key, "mbar_rh", tend) / med_at(key, "mbar", tend)
        lines.append(f"  {key:<12}{ratio0:>10.3f}{ratioe:>12.3f}"
                     f"{mr0:>21.3f}{mre:>9.3f}")
    lines.append("")

    lines += ["-" * 74, "4. 能量守恒检验（每个随机实现）", "-" * 74]
    lines += [f"  {'配置':<12}{'run':>5}{'seed':>8}{'|dE/E|':>12}"
              f"{'maxModify/|E0|':>18}{'末态束缚':>12}"]
    for cfg in CONFIGS:
        for r in all_runs.get(cfg["key"], []):
            e = r["energy"]
            if e["t"].size == 0:
                continue
            de = abs((e["etot"][-1] - e["etot"][0]) / e["etot"][0])
            mm = np.abs(e["emod"]).max() / abs(e["etot"][0])
            seed = "-"
            info = os.path.join(r["dir"], "RUN_INFO.txt")
            if os.path.isfile(info):
                with open(info) as fh:
                    for ln in fh:
                        if ln.startswith("seed"):
                            seed = ln.split("=")[1].strip()
                            break
            bound = "束缚" if e["etot"][-1] < 0 else "不束缚"
            lines.append(f"  {cfg['key']:<12}{r['run']:>5}{seed:>8}{de:>12.3e}"
                         f"{mm:>18.3e}{bound:>12}")
    lines.append("")
    lines.append(A)

    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(HERE, "summary.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(f"\n[写入] {os.path.join(HERE, 'summary.txt')}")


# ===========================================================================
# 5. 缓存指标（便于绘图与报告脚本复用）
# ===========================================================================
def save_cache(all_runs, t_rh):
    arrays = {}
    for cfg in CONFIGS:
        runs = all_runs.get(cfg["key"], [])
        if not runs:
            continue
        arrays[f"{cfg['key']}_time"] = _ref_grid(runs)
        tref = arrays[f"{cfg['key']}_time"]
        for k in METRIC_KEYS:
            v = align(runs, k, tref)
            if v is not None:
                arrays[f"{cfg['key']}_{k}"] = v
        if "rc" in runs[0]:
            ctref = _ref_grid(runs, "core_time")
            arrays[f"{cfg['key']}_ctime"] = ctref
            arrays[f"{cfg['key']}_rc"] = align(runs, "rc", ctref, "core_time")
    for k, v in t_rh.items():
        arrays[f"trh_{k}"] = np.array([v])
    np.savez(CACHE, **arrays)
    print(f"[写入] {CACHE}")


def main():
    all_runs = {}
    for cfg in CONFIGS:
        runs = collect_runs(cfg)
        all_runs[cfg["key"]] = runs
        print(f"[数据] {cfg['key']}: {len(runs)} 个实现, "
              f"{len(runs[0]['time']) if runs else 0} 个快照")

    if not any(all_runs.values()):
        print("错误: 未找到任何数据，请先运行 run_all.sh", file=sys.stderr)
        return 1

    t_rh = relaxation_times(all_runs)
    print("[弛豫时间] " + ", ".join(f"{k} = {v:.2f} Myr" for k, v in t_rh.items()))

    os.makedirs(FIGDIR, exist_ok=True)
    # 主图（题目要求的两个：r_h 与 r_c 演化对比，横轴为时间 Myr）
    plot_rh(all_runs, t_rh)
    plot_rc(all_runs, t_rh)
    # 额外诊断图 -> figs/extra/
    plot_rh(all_runs, t_rh, extra=True)
    plot_rc(all_runs, t_rh, extra=True)
    plot_panel(all_runs, t_rh, extra=True)
    plot_norm(all_runs, extra=True)
    plot_segregation(all_runs, t_rh, extra=True)
    plot_energy(all_runs, extra=True)
    plot_massprofile(all_runs, extra=True)
    write_summary(all_runs, t_rh)
    save_cache(all_runs, t_rh)
    return 0


if __name__ == "__main__":
    sys.exit(main())
