#!/usr/bin/env python3
"""
作业1 结果分析脚本

读取 PeTar 后处理生成的数据，产出：
  1. fig_lagr.png     拉格朗日半径演化图 (R10/R30/R50/R70/R90 + 核心半径)
  2. fig_energy.png   能量守恒检验图
  3. summary.txt      关键数值汇总

用法:  python3 analyze.py [数据目录, 默认 ./output]
"""

import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import petar

# ---------------------------------------------------------------------------
# 中文字体设置
#  注意：不能用 Noto Sans CJK JP 之类的**日文**字体——日文字形与简体中文有差异
#  （例如「径」的右半在日文字体里是「ス+土」，中文是「又+土」）。
#  这里只接受简体中文字体；找不到就退回英文，避免出现错字或方块。
# ---------------------------------------------------------------------------
_CJK_CANDIDATES = (
    "WenQuanYi Micro Hei",      # 文泉驿微米黑（简体）
    "WenQuanYi Zen Hei",        # 文泉驿正黑（简体）
    "Noto Sans CJK SC",         # 思源黑体简体
    "Source Han Sans SC",
    "Source Han Sans CN",
    "Microsoft YaHei",          # 微软雅黑（Windows）
    "SimHei",
    "PingFang SC",              # macOS
    "Heiti SC",
    "AR PL UMing CN",
)

# matplotlib 对 .ttc 字体集合只索引第一个字面，常只暴露出日文变体，
# 因此优先直接注册简体中文字体文件
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
if CJK:
    print(f"[字体] 图中文字使用简体中文字体: {CJK_FONT}")
else:
    print("[字体] 未找到简体中文字体，图注改用英文")


def L(zh, en):
    """有简体中文字体时返回中文，否则返回英文"""
    return zh if CJK else en

DATADIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "output")
FIGDIR = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# 1. 拉格朗日半径演化
# ---------------------------------------------------------------------------
def plot_lagrangian():
    lagr = petar.LagrangianMultiple()
    lagr.fromfile(os.path.join(DATADIR, "data.lagr"))

    t = lagr.time                       # (n_time,)
    mfrac = lagr.initargs["mass_fraction"]     # [0.1, 0.3, 0.5, 0.7, 0.9]
    R = lagr.all.r                      # (n_time, n_frac+1), 最后一列是核心半径 Rc

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    cmap = plt.cm.viridis(np.linspace(0.05, 0.85, len(mfrac)))
    for i, f in enumerate(mfrac):
        lbl = f"$R_{{{int(round(f * 100))}}}$"
        ax.plot(t, R[:, i], color=cmap[i], lw=1.6, label=lbl)
    ax.plot(t, R[:, -1], color="crimson", lw=1.4, ls="--", label="$R_c$ (core)")

    ax.set_yscale("log")
    ax.set_xlabel("$t$ [Myr]", fontsize=12)
    ax.set_ylabel(L("半径 [pc]", "Radius [pc]"), fontsize=12)
    ax.set_title(L("拉格朗日半径演化  (N=1000, Plummer, $R_h$ = 1 pc)",
                   "Lagrangian radii evolution  (N=1000, Plummer, $R_h$ = 1 pc)"),
                 fontsize=12)
    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(0.008, 600)
    ax.grid(alpha=0.3, which="both", ls=":")
    ax.legend(ncol=2, fontsize=10, frameon=True, loc="lower right",
              framealpha=0.95)

    # 标注半质量半径 R50 的初值与终值
    i50 = int(np.argmin(np.abs(mfrac - 0.5)))
    ann_kw = dict(textcoords="axes fraction", fontsize=9, color=cmap[i50],
                  bbox=dict(boxstyle="round,pad=0.25", fc="white",
                            ec=cmap[i50], lw=0.7, alpha=0.9),
                  arrowprops=dict(arrowstyle="->", color=cmap[i50], lw=0.9))
    ax.annotate(f"$R_{{50}}$ = {R[0, i50]:.2f} pc", xy=(t[0], R[0, i50]),
                xytext=(0.30, 0.13), **ann_kw)
    ax.annotate(f"$R_{{50}}$ = {R[-1, i50]:.2f} pc", xy=(t[-1], R[-1, i50]),
                xytext=(0.60, 0.60), **ann_kw)

    fig.tight_layout()
    out = os.path.join(FIGDIR, "fig_lagr.png")
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[写入] {out}")

    return t, mfrac, R


# ---------------------------------------------------------------------------
# 2. 能量守恒检验
# ---------------------------------------------------------------------------
def read_energy(logfile):
    """从 PeTar 日志中提取每个输出时刻的能量诊断。

    日志格式::

        Time: <t>  N_real(loc): ... N_escape(glb): ...
        Energy:  Error/Total  Error  Error_cum  Total  Kinetic  Potential  ...
        Physic:  <误差/总能量> <误差> <累计误差> <总能量> <动能> <势能> <Modify> ...
    """
    t, etot, ekin, epot, eerr, emod = [], [], [], [], [], []
    cur_t = None
    with open(logfile, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith("Time:"):
                cur_t = float(line.split()[1])
            elif line.startswith("Physic:") and cur_t is not None:
                v = line.split()[1:]
                if len(v) >= 7:
                    eerr.append(float(v[0]))    # Error/Total
                    etot.append(float(v[3]))    # Total
                    ekin.append(float(v[4]))    # Kinetic
                    epot.append(float(v[5]))    # Potential
                    emod.append(float(v[6]))    # Modify
                    t.append(cur_t)
    return (np.array(t), np.array(etot), np.array(ekin),
            np.array(epot), np.array(eerr), np.array(emod))


def plot_energy():
    logfile = os.path.join(DATADIR, "output")
    t, etot, ekin, epot, eerr, emod = read_energy(logfile)
    if t.size == 0:
        print("[跳过] 未在日志中找到能量诊断")
        return None

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    ax = axes[0]
    ax.plot(t, etot, "k-", lw=1.6, label="Total $E$")
    ax.plot(t, ekin, "C3-", lw=1.0, alpha=0.8, label="Kinetic $T$")
    ax.plot(t, epot, "C0-", lw=1.0, alpha=0.8, label="Potential $W$")
    ax.set_xlabel("$t$ [Myr]", fontsize=11)
    ax.set_ylabel("Energy  [$M_\\odot\\,$pc$^2$\\,Myr$^{-2}$]", fontsize=11)
    ax.set_title("Energy budget", fontsize=11)
    ax.grid(alpha=0.3, ls=":")
    ax.legend(fontsize=9)

    ax = axes[1]
    ax.plot(t, np.abs(eerr), "C4-", lw=1.2)
    ax.axhline(1e-3, color="r", ls="--", lw=0.9, label="1e-3 tolerance")
    ax.set_yscale("log")
    ax.set_xlabel("$t$ [Myr]", fontsize=11)
    ax.set_ylabel("$|\\Delta E / E|$", fontsize=11)
    ax.set_title("Relative energy error", fontsize=11)
    ax.grid(alpha=0.3, which="both", ls=":")
    ax.legend(fontsize=9)

    fig.tight_layout()
    out = os.path.join(FIGDIR, "fig_energy.png")
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[写入] {out}")

    return t, etot, ekin, epot, eerr, emod


# ---------------------------------------------------------------------------
# 3. 数值稳定性对比图：通过的实现 vs 被拒绝的实现
# ---------------------------------------------------------------------------
def _latest_evidence_log():
    """在 evidence/ 下找一个".hard_large_energy"或失稳的运行日志"""
    root = os.path.join(FIGDIR, "evidence")
    if not os.path.isdir(root):
        return None, None
    best, label = None, None
    for d in sorted(os.listdir(root)):
        p = os.path.join(root, d, "petar.log")
        if os.path.isfile(p) and os.path.getsize(p) > 20000:
            best, label = p, d
    return best, label


def plot_energy_compare():
    good_log = os.path.join(DATADIR, "output")
    bad_log, bad_label = _latest_evidence_log()
    if bad_log is None:
        return

    print(f"[对比] 失稳实现: {bad_label}")

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for log, lab, col in ((good_log, L("通过自检", "accepted"), "C0"),
                          (bad_log, L("失稳", "unstable") + f" ({bad_label})", "C3")):
        t, etot, ekin, epot, eerr, emod = read_energy(log)
        if t.size == 0:
            continue
        axes[0].plot(t, np.abs(emod), col, lw=1.4, label=lab)
        axes[1].plot(t, etot / etot[0], col, lw=1.6, label=lab)

    axes[0].set_yscale("log")
    axes[0].set_xlabel("$t$ [Myr]", fontsize=11)
    axes[0].set_ylabel("$|{\\rm Modify}|$", fontsize=11)
    axes[0].set_title(L("SDAR 累计能量修正", "Cumulative SDAR energy correction"),
                      fontsize=11)
    axes[0].grid(alpha=0.3, which="both", ls=":")
    axes[0].legend(fontsize=9)

    axes[1].axhline(0.0, color="k", lw=0.8)
    axes[1].set_xlabel("$t$ [Myr]", fontsize=11)
    axes[1].set_ylabel("$E(t)/E(0)$", fontsize=11)
    axes[1].set_title(L("总能量演化（0 以下 = 束缚）", "Total energy  ($E<0$ = bound)"),
                      fontsize=11)
    axes[1].grid(alpha=0.3, ls=":")
    axes[1].legend(fontsize=9)

    fig.tight_layout()
    out = os.path.join(FIGDIR, "fig_energy_compare.png")
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[写入] {out}")


# ---------------------------------------------------------------------------
# 4. 数值汇总
# ---------------------------------------------------------------------------
def write_summary(t, mfrac, R, ener):
    lines = []
    lines.append("=" * 68)
    lines.append("作业1  星团数值模拟结果汇总")
    lines.append("=" * 68)
    lines.append("")
    lines.append("参数:")
    lines.append("  N (粒子数)        = 1000")
    lines.append("  密度剖面          = Plummer (mcluster -P 0)")
    lines.append("  质量函数          = Kroupa (2001)  (mcluster -f 1)")
    lines.append("  半质量半径 R_h    = 1 pc")
    lines.append("  积分时间          = 200 Myr")
    lines.append("  快照间隔          = 1 Myr")
    lines.append("")
    lines.append("-" * 68)
    lines.append("1. 拉格朗日半径演化  [pc]")
    lines.append("-" * 68)
    header = "  t [Myr] " + "".join(f"  R{int(round(f*100)):<8}" for f in mfrac) \
             + "     Rc"
    lines.append(header)
    for tt in [0, 25, 50, 75, 100, 125, 150, 175, 200]:
        i = int(np.argmin(np.abs(t - tt)))
        row = f"  {t[i]:7.1f} " + "".join(f"  {R[i, j]:<9.4g}" for j in range(len(mfrac))) \
              + f"  {R[i, -1]:.4g}"
        lines.append(row)

    lines.append("")
    lines.append("  半质量半径 R50 :  %.3f pc  ->  %.3f pc   (增大 %.1f 倍)"
                 % (R[0, 2], R[-1, 2], R[-1, 2] / R[0, 2]))
    lines.append("  核心半径   Rc  :  %.3f pc  ->  %.3f pc   (增大 %.1f 倍)"
                 % (R[0, -1], R[-1, -1], R[-1, -1] / R[0, -1]))

    if ener is not None:
        t_e, etot, ekin, epot, eerr, emod = ener
        lines.append("")
        lines.append("-" * 68)
        lines.append("2. 能量守恒检验")
        lines.append("-" * 68)
        lines.append(f"  初始总能量 E0            = {etot[0]:.4f}")
        lines.append(f"  末态总能量 E(t_end)      = {etot[-1]:.4f}")
        lines.append(f"  总能量相对漂移 |dE/E|    = {abs((etot[-1]-etot[0])/etot[0]):.3e}")
        lines.append(f"  单步最大相对误差         = {np.abs(eerr).max():.3e}")
        lines.append(f"  系统末态束缚状态         = "
                     f"{'束缚 (E<0)' if etot[-1] < 0 else '不束缚 (E>0)'}")

    lines.append("")
    lines.append("=" * 68)
    txt = "\n".join(lines)
    print(txt)
    out = os.path.join(FIGDIR, "summary.txt")
    with open(out, "w") as fh:
        fh.write(txt + "\n")
    print(f"\n[写入] {out}")


def main():
    t, mfrac, R = plot_lagrangian()
    ener = plot_energy()
    plot_energy_compare()
    write_summary(t, mfrac, R, ener)


if __name__ == "__main__":
    main()
