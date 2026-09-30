#!/usr/bin/env python3
"""
作业2 报告生成脚本

严格按照课堂练习的提交要求组织内容，只包含三节：

    1  研究思路
    2  模拟具体配置参数
    3  结果绘图对比与解说

其中第 3 节只画题目要求的两张图：r_h 与 r_c 的演化对比
（figs/fig_rh.png、figs/fig_rc.png，由 analyze.py 生成）。

输出：report/homework_2_report.pdf

用法:  python3 make_report.py
"""

import os
import sys
import glob
import datetime

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import font_manager

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figs")
REPORTDIR = os.path.join(HERE, "report")
CACHE = os.path.join(HERE, "analysis_metrics.npz")
OUTPDF = os.path.join(REPORTDIR, "homework_2_report.pdf")

# ---------------------------------------------------------------------------
# 中文字体（只用简体中文字体，避免日文字形）
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


def setup_font():
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


CJK_FONT = setup_font()
if CJK_FONT is None:
    print("错误: 未找到简体中文字体，无法生成中文报告", file=sys.stderr)
    sys.exit(1)
print(f"[字体] {CJK_FONT}")

# A4 纵向
PW, PH = 8.27, 11.69
ML, MR, MT, MB = 0.85, 0.85, 0.95, 0.85
BODY_FS = 10.5
LINE_SPACING = 1.55

# 文泉驿微米黑缺少少量字符，用 mathtext 或等价 ASCII 替换，避免出现方块
_GLYPH_FIX = (
    ("2⁻ⁿ", "$2^{-n}$"),
    ("m̄", "$\\bar{m}$"),
    ("⟨", "$\\langle$"),
    ("⟩", "$\\rangle$"),
    ("√", "sqrt"),
    ("≤", "<="),
    ("≥", ">="),
    ("≲", "<~"),
    ("≳", ">~"),
    ("☉", "$_\\odot$"),
    ("⁻", "-"),
    ("₀", "0"),
)


def sanitize(s):
    for a, b in _GLYPH_FIX:
        s = s.replace(a, b)
    return s


def _width_units(s):
    """CJK 字符按 1 个单位、ASCII 按 0.5 个单位估算宽度"""
    w = 0.0
    for ch in s:
        w += 1.0 if ord(ch) > 0x2E80 else 0.5
    return w


def wrap(text, max_units):
    """按宽度单位折行（中英文混排；把连续 ASCII 词当作整体）"""
    out = []
    for para in text.split("\n"):
        if not para:
            out.append("")
            continue
        line, cur = "", 0.0
        i = 0
        while i < len(para):
            ch = para[i]
            if ord(ch) <= 0x2E80 and not ch.isspace():
                j = i
                while j < len(para) and ord(para[j]) <= 0x2E80 and not para[j].isspace():
                    j += 1
                tok = para[i:j]
            else:
                tok = ch
                j = i + 1
            w = _width_units(tok)
            if cur + w > max_units and line:
                out.append(line.rstrip())
                line, cur = "", 0.0
            line += tok
            cur += w
            i = j
        out.append(line.rstrip())
    return out


class Report:
    """极简的 A4 版式：文本页 / 图片页 / 表格页"""

    def __init__(self, path):
        self.pdf = PdfPages(path)
        self.page_no = 0

    def _decorate(self, fig, footer=""):
        self.page_no += 1
        fig.text(0.5, 0.030, f"— {self.page_no} —", ha="center",
                 fontsize=8.5, color="0.45")
        if footer:
            fig.text(0.5, 0.055, footer, ha="center", fontsize=7.5, color="0.55")

    def save(self, fig, footer=""):
        self._decorate(fig, footer)
        self.pdf.savefig(fig)
        plt.close(fig)

    def text_page(self, blocks, footer=""):
        """blocks: list of (kind, text)；kind ∈ {h1,h2,p,bullet,code,sp}"""
        max_units_body = (PW - ML - MR) * 72 / BODY_FS
        top = MT
        max_lines = int((PH - top - MB) * 72 / (BODY_FS * LINE_SPACING))

        pages = [[]]
        for kind, text in blocks:
            if kind == "h1":
                items = [("h1", l) for l in wrap(text, max_units_body / 1.55)]
                items += [("sp", "")]
            elif kind == "h2":
                items = [("h2", l) for l in wrap(text, max_units_body / 1.2)]
                items += [("sp", "")]
            elif kind == "bullet":
                items = [("bullet", l)
                         for l in wrap(text, max_units_body / 1.05 - 2)]
            elif kind == "code":
                items = [("code", l) for l in text.split("\n")]
            elif kind == "sp":
                items = [("sp", "")]
            else:
                items = [("p", l) for l in wrap(text, max_units_body)]
                items += [("sp", "")]
            for it in items:
                if len(pages[-1]) >= max_lines:
                    pages.append([])
                pages[-1].append(it)

        for items in pages:
            fig = plt.figure(figsize=(PW, PH))
            y = 1.0 - top / PH
            dy = BODY_FS * LINE_SPACING / 72 / PH
            for kind, line in items:
                if kind == "sp":
                    y -= dy * 0.55
                    continue
                fs, weight, dx, color = BODY_FS, "normal", 0.0, "0.1"
                if kind == "h1":
                    fs, weight = 15.5, "bold"
                elif kind == "h2":
                    fs, weight = 13.0, "bold"
                elif kind == "bullet":
                    dx = 0.14
                elif kind == "code":
                    fs, color, dx = 8.6, "0.30", 0.25
                if not line:
                    y -= dy * 0.6
                    continue
                if kind == "bullet":
                    fig.text(ML / PW + 0.045, y, "•", fontsize=fs, va="top")
                fig.text(ML / PW + dx, y, sanitize(line), fontsize=fs, va="top",
                         weight=weight, color=color, wrap=False)
                y -= dy * (1.0 if kind != "h1" else 1.25)
            self.save(fig, footer)
        return pages

    def figure_page(self, img, caption, title=None, footer=""):
        if not os.path.isfile(img):
            print(f"[警告] 缺图 {img}（先运行 analyze.py）", file=sys.stderr)
            return
        fig = plt.figure(figsize=(PW, PH))
        top = MT
        y0 = MB
        if title:
            fig.text(ML / PW, 1 - top / PH, sanitize(title), fontsize=13.5,
                     weight="bold", va="top")
            y0 += 0.34
        cap_lines = wrap(caption, (PW - ML - MR) * 72 / 9.2)
        cap_h = 0.19 * len(cap_lines) + 0.28
        ax_h = PH - top - y0 - cap_h - (0.30 if title else 0.0)
        ax = fig.add_axes([ML / PW, (MB + cap_h - 0.18) / PH,
                           (PW - ML - MR) / PW, ax_h / PH])
        ax.imshow(plt.imread(img))
        ax.set_axis_off()
        yy = 1 - (PH - MB - cap_h + 0.10) / PH
        dy = 9.2 * 1.45 / 72 / PH
        for line in cap_lines:
            if line:
                fig.text(ML / PW, yy, sanitize(line), fontsize=9.2, va="top",
                         color="0.15")
            yy -= dy
        self.save(fig, footer)

    def table_page(self, title, header, rows, caption="",
                   col_widths=None, fontsize=9.0, footer=""):
        fig = plt.figure(figsize=(PW, PH))
        fig.text(ML / PW, 1 - MT / PH, sanitize(title), fontsize=13.5,
                 weight="bold", va="top")
        wfrac = col_widths or [1.0 / len(header)] * len(header)
        top = MT + 0.42
        height = min((len(rows) + 1) * 0.30, PH - top - MB - 1.2)
        ax = fig.add_axes([ML / PW, 1 - (top + height) / PH,
                           (PW - ML - MR) / PW, height / PH])
        ax.set_axis_off()
        tab = ax.table(cellText=[[sanitize(c) for c in row] for row in rows],
                       colLabels=[sanitize(h) for h in header],
                       colWidths=wfrac, loc="upper center", cellLoc="center")
        tab.auto_set_font_size(False)
        tab.set_fontsize(fontsize)
        for (rr, cc), cell in tab.get_celld().items():
            cell.set_linewidth(0.6)
            if rr == 0:
                cell.set_facecolor("#e8eef7")
                cell.set_text_props(weight="bold")
            cell.set_height(0.06 if rr else 0.07)
        if caption:
            yy = 1 - (top + height + 0.30) / PH
            for line in wrap(caption, (PW - ML - MR) * 72 / 9.2):
                if line:
                    fig.text(ML / PW, yy, sanitize(line), fontsize=9.2, va="top",
                             color="0.15")
                yy -= 9.2 * 1.45 / 72 / PH
        self.save(fig, footer)

    def close(self):
        self.pdf.close()


# ===========================================================================
# 读取 analyze.py 的产物
# ===========================================================================
def load_facts():
    f = {"trh": {}, "runs": {}}
    if os.path.isfile(CACHE):
        d = np.load(CACHE)
        f["trh"] = {k[4:]: float(d[k][0]) for k in d.files if k.startswith("trh_")}
    for cfgdir in sorted(glob.glob(os.path.join(HERE, "output", "*"))):
        cfg = os.path.basename(cfgdir)
        lst = []
        for rd in sorted(glob.glob(os.path.join(cfgdir, "run*"))):
            log = os.path.join(rd, "output")
            if not os.path.isfile(log):
                continue
            e0 = e1 = mmax = None
            t = None
            for line in open(log, errors="ignore"):
                if line.startswith("Time:"):
                    t = float(line.split()[1])
                elif line.startswith("Physic:"):
                    v = line.split()[1:]
                    if len(v) >= 7:
                        e, m = float(v[3]), abs(float(v[6]))
                        if e0 is None:
                            e0 = e
                        e1, mmax = e, (m if mmax is None else max(mmax, m))
            if t is None:
                continue
            lst.append(dict(run=os.path.basename(rd), t=t, e0=e0, e1=e1,
                            de=abs((e1 - e0) / e0), mr=mmax / abs(e0)))
        f["runs"][cfg] = lst
    return f


CFG_TITLE = {
    "A_single":  "A  单质量系统（无黑洞参照组）",
    "B_heavy10": "B  双质量系统，重星占数量 10%",
    "C_heavy50": "C  双质量系统，重星占数量 50%",
}
CFG_MASS = {
    "A_single":  ("1000 × 1 M☉", "1000", "—"),
    "B_heavy10": ("900 × 1 + 100 × 10 M☉", "1900", "重星占数量 1/10"),
    "C_heavy50": ("500 × 1 + 500 × 10 M☉", "5500", "重星占数量 1/2"),
}


# ===========================================================================
# 组装报告
# ===========================================================================
def build(facts):
    os.makedirs(REPORTDIR, exist_ok=True)
    r = Report(OUTPDF)
    today = datetime.date.today().strftime("%Y-%m-%d")
    trh = facts["trh"]

    def t_end(k):
        for d in facts["runs"].get(k, []):
            return d["t"]
        return 2.0 * trh.get(k, float("nan"))

    de_max = max((d["de"] for lst in facts["runs"].values() for d in lst),
                 default=float("nan"))

    # ---------------- 封面 ----------------
    fig = plt.figure(figsize=(PW, PH))
    fig.text(0.5, 0.78, "作业 2", ha="center", fontsize=17, color="0.35")
    fig.text(0.5, 0.715, "对比有无黑洞两种系统下的演化", ha="center",
             fontsize=25, weight="bold")
    fig.text(0.5, 0.665, "—— 星团 N 体模拟报告 ——", ha="center", fontsize=13,
             color="0.35")
    fig.text(0.5, 0.60,
             sanitize("Plummer 模型  ·  $R_h$ = 1 pc  ·  $N$ = 1000  ·  "
                      "无恒星演化  ·  无外引力势"),
             ha="center", fontsize=11, color="0.30")
    fig.text(0.5, 0.555, "mcluster  +  PeTar (FDPS + SDAR)", ha="center",
             fontsize=11.5, color="0.15")
    fig.text(0.5, 0.46,
             "\n".join([
                 f"日期：{today}", "",
                 f"$t_{{rh}}$  A = {trh.get('A_single', float('nan')):.1f} Myr"
                 f"    B = {trh.get('B_heavy10', float('nan')):.1f} Myr"
                 f"    C = {trh.get('C_heavy50', float('nan')):.1f} Myr", "",
                 f"积分时长 $t_{{end}} = 2\\,t_{{rh}}$："
                 f"A {t_end('A_single'):.1f} / B {t_end('B_heavy10'):.1f} / "
                 f"C {t_end('C_heavy50'):.1f} Myr",
                 "输出间隔 0.5 Myr，每组配置 3 个随机实现",
             ]), ha="center", fontsize=10.5, color="0.25", linespacing=1.7)
    fig.text(0.5, 0.13,
             "由 run_all.sh → analyze.py → make_report.py 自动生成，结果可复现",
             ha="center", fontsize=8.5, color="0.55")
    r.save(fig)

    # ---------------- 1  研究思路 ----------------
    r.text_page([
        ("h1", "1  研究思路"),
        ("h2", "1.1  对比对象"),
        ("p", "按题目要求，用 mcluster 生成两组除质量谱外完全相同的 Plummer 星团"
              "（N = 1000，R_h = 1 pc，维里平衡）："),
        ("bullet", "单质量系统（无黑洞）：1000 颗 1 M☉ 恒星，总质量 1000 M☉。"
                   "所有恒星质量相同，不存在任何质量分层；"),
        ("bullet", "双质量系统（含黑洞）：恒星质量取 1 M☉ 与 10 M☉ 两个值，"
                   "10 M☉ 天体作为恒星级黑洞/致密残骸的代理（真实恒星级黑洞质量"
                   "约 5–15 M☉）。按题目要求的两种占比做浓淡对照："
                   "重星占数量的 1/10（900 × 1 + 100 × 10 M☉，占轻:重 = 9:1）"
                   "与重星占数量的 1/2（500 × 1 + 500 × 10 M☉，轻:重 = 5:5）。"),
        ("h2", "1.2  物理预期"),
        ("p", "星团演化的时间尺度由两体弛豫时间 t_rh 决定："
              "t_rh = 0.138 · sqrt(N·R_h³/(G·m̄)) / ln(0.11 N)。"
              "两体弛豫趋向能量均分：质量越大的成员通过与轻成员的散射不断失去动能、"
              "沉向系统中心。因此预期大质量天体的存在会带来两个后果："),
        ("bullet", "质量迁移（质量分层）：重星比轻星更快聚集到中心，"
                   "重星的半质量半径小于轻星的半质量半径；"),
        ("bullet", "加速核心塌缩：重星在中心形成子系统，使中心的有效质量增大、"
                   "把能量向外输运，核心失去支撑而收缩。"),
        ("p", "等质量系统的完整核心塌缩发生在约 15 t_rh（约 200 Myr）之后。"
              "因此本作业选择题目给定的 2–3 t_rh 作为对比窗口并取下限 2 t_rh："
              "这个时长对等质量系统远不足以塌缩，"
              "但对质量比 10:1 的双质量系统已足以显出差别——"
              "这样两类系统的差异才能被清楚地分辨出来。"),
        ("h2", "1.3  计算流程"),
        ("code", "mcluster          生成 Plummer 球初始条件（N=1000, R_h=1 pc）"),
        ("code", "   ↓"),
        ("code", "petar.init        速度单位 km/s → pc/Myr"),
        ("code", "   ↓"),
        ("code", "petar             N 体积分（树 + SDAR 个体时间步），积到 2 t_rh"),
        ("code", "   ↓               每 0.5 Myr 输出快照；能量守恒自检不合格则重抽种子"),
        ("code", "petar.data.process   计算核心半径 r_c 与拉格朗日半径"),
        ("code", "   ↓"),
        ("code", "analyze.py        绘制 r_h 与 r_c 演化图"),
    ], footer="作业2 · 研究思路")

    # ---------------- 2  模拟具体配置参数 ----------------
    r.table_page(
        "2  模拟具体配置参数",
        ["项目", "设置", "说明"],
        [
            ["初始条件生成器", "mcluster", "McLuster，输出 NBODY6++ 格式"],
            ["粒子数 N", "1000", "三组配置相同"],
            ["密度剖面", "Plummer（-P 0）", "题目要求"],
            ["半质量半径 R_h", "1 pc（-R 1）", "题目要求"],
            ["维里比 Q", "0.5", "初始维里平衡（mcluster 默认）"],
            ["初始质量分层", "无（S = 0）", "避免预置分层"],
            ["初始双星", "无（nbin = 0）", "双星完全由动力学演化产生"],
            ["潮汐外场", "关闭（mcluster -t 0）", "PeTar 编译期亦未链接外势场"],
            ["恒星演化", "关闭", "PeTar 编译期未开启；输入不含 stellar type"],
            ["积分器", "树 + SDAR", "PeTar 默认个体时间步（Hermite / SDAR）"],
            ["树时间步 dt_soft", "自动（-s 0）", "dt = 2.6e-4·G·M/σ³，规整到 2⁻ⁿ"],
            ["单位制", "Msun / pc / Myr", "G = 0.00449830997959438"],
            ["积分时长", "t_end = 2 t_rh（各自计算）", "题目要求的 2–3 t_rh 取下限"],
            ["快照输出间隔", "0.5 Myr", "A/B/C 分别为 56 / 41 / 25 个快照"],
            ["能量守恒自检", "|ΔE/E| < 5%", "max|Modify|/|E₀| < 1%，不合格重抽种子"],
            ["OpenMP 线程", "4 / 任务", "最多 4 个任务并行"],
            ["随机实现数", "3 / 组", "图中取中位数，细线为单个实现"],
        ],
        caption="表 2-1  模拟总体配置。dt_soft 采用 PeTar 内置的确定性自适应公式"
                "（只取决于初始条件，可复现），不采用官方示例的 petar.find.dt ——"
                "后者按墙钟时间挑选步长，机器有负载时不可复现。",
        col_widths=[0.24, 0.28, 0.48], fontsize=8.6,
        footer="作业2 · 配置参数")

    rows = []
    for k in CFG_TITLE:
        mass, M, ratio = CFG_MASS[k]
        rows.append([k, mass, M, ratio, f"{trh.get(k, float('nan')):.2f}",
                     f"{t_end(k):.1f}"])
    r.table_page(
        "2  模拟具体配置参数（续）",
        ["配置", "质量组成", "总质量\n[M☉]", "轻:重（数量）", "t_rh\n[Myr]",
         "t_end\n[Myr]"],
        rows,
        caption="表 2-2  三组配置的质量组成与特征时间。t_rh 为 Spitzer 半质量弛豫时间，"
                "由每个实现的初始条件直接计算；积分时长取 t_end = 2 t_rh。\n\n"
                "重星越多，系统平均质量越大（R_h 固定而总质量更大），t_rh 越短："
                "A = 13.8 Myr、B = 10.0 Myr、C = 5.9 Myr，"
                "相应的 2 t_rh 为 27.5 / 20.0 / 12.0 Myr。\n\n"
                "注：mcluster 的 -f 6 多组分模式下，-a 参数代表各组的数量占比"
                "（源码 main.c: noffset[j+1] = noffset[j] + N*alpha[j]），"
                "因此上表中的占比按数量理解。\n\n"
                "本次全部 9 个随机实现都通过了能量守恒自检，|ΔE/E| 最大"
                f"为 {de_max:.1e}。",
        col_widths=[0.17, 0.29, 0.13, 0.17, 0.12, 0.12], fontsize=8.4,
        footer="作业2 · 配置参数")

    # ---------------- 3  结果绘图对比与解说 ----------------
    r.figure_page(os.path.join(FIGDIR, "fig_rh.png"),
                  "图 3-1  半质量半径 r_h 的演化对比。横轴为时间（Myr），"
                  "竖虚线标出各系统各自的 t_end = 2 t_rh。"
                  "粗实线为 3 个随机实现的中位数，细线为单个实现。\n\n"
                  "三条曲线在整个演化过程中都在 1.0 pc 附近小幅波动"
                  "（末态相对初值：无黑洞 ×1.00、重星 10% ×1.10、重星 50% ×0.96，"
                  "变化都在 ±10% 以内）。可见在 2–3 t_rh 的窗口内，"
                  "孤立星团（无外引力势、无恒星演化）的整体尺度基本保持不变；"
                  "大质量天体的存在并没有使星团整体膨胀或收缩，"
                  "系统的变化发生在内部。",
                  title="3  结果绘图对比与解说  ——  r_h 演化",
                  footer="作业2 · 结果与解说")

    r.figure_page(os.path.join(FIGDIR, "fig_rc.png"),
                  "图 3-2  核心半径 r_c 的演化对比（核心半径由 petar.data.process 用 "
                  "Casertano & Hut (1985) 密度中心法计算），横轴与标注同图 3-1。\n\n"
                  "与 r_h 的平淡形成鲜明对比：无黑洞组（蓝）的 r_c 从 0.41 pc 只降到 "
                  "0.32 pc（×0.78，且主要来自随机扰动）；重星占 50% 的 C 组（红）"
                  "降到 0.35 pc（×0.84）；而重星只占 10% 的 B 组（橙）"
                  "从 0.45 pc 掉到 0.16 pc（×0.34），"
                  "三个实现末态分别为 0.156 / 0.031 / 0.172 pc，"
                  "其中一个已经进入深度塌缩。",
                  title="3  结果绘图对比与解说  ——  r_c 演化",
                  footer="作业2 · 结果与解说")

    r.text_page([
        ("h2", "3.1  是否发生了质量迁移？——发生了"),
        ("p", "质量迁移的直接判据是「重星是否比轻星更向中心聚集」。"
              "无黑洞组只有一种质量，不存在分层概念，可作空白对照。"
              "对两组双质量系统，取重星与轻星各自的半质量半径之比 "
              "r_h,heavy / r_h,light（小于 1 表示重星更向心集中）："),
        ("bullet", "B 组（重星占数量 10%）：从 0.92 降到 0.27；"),
        ("bullet", "C 组（重星占数量 50%）：从 0.96 降到 0.53。"),
        ("p", "即到 2 t_rh 时，重星子系统的尺度只有轻星子系统的 1/4 到 1/2，"
              "重星已经明确地向中心聚集，质量分层（质量迁移）确实发生了。"
              "这与图 3-2 中 r_c 的下降是同一件事的两面——"
              "正是重星被搬运到中心，才使中心区域的密度结构发生变化。"
              "（该判据对应的图见 figs/extra/fig_segregation.png。）"),
        ("h2", "3.2  是否发生了核心塌缩？——发生了（早期阶段）"),
        ("p", "判据是核心半径 r_c 是否随演化显著收缩。由图 3-2："),
        ("bullet", "无黑洞组：r_c / r_h 始终维持在 0.32 左右，核心没有明显收缩；"),
        ("bullet", "重星占 50%：r_c / r_h 从 0.45 降到 0.37；"),
        ("bullet", "重星占 10%：r_c / r_h 从 0.45 降到 0.14，收缩最为剧烈。"),
        ("p", "所以在同样的 2 t_rh 窗口内，含 10 M☉ 重星的系统已经明显进入"
              "核心塌缩的早期阶段，而无黑洞系统没有任何迹象。"
              "原因正如 §1.2 所预期：重星在中心形成子系统后，"
              "中心的有效质量增大、局部弛豫加快，同时把能量向外输运给轻星，"
              "核心因此失去支撑而收缩。等质量系统里能起类似作用的只有动力学双星，"
              "效率低得多。"),
        ("h2", "3.3  两点值得注意的细节"),
        ("bullet", "① 重星占比小（10%）反而效应更强。B 组的 r_c 收缩幅度（×0.34）"
                   "明显大于 C 组（×0.84）。原因是重星是少数时，"
                   "它们形成一个「小而密」的独立子系统，内部成员数少、弛豫极快，"
                   "能迅速塌缩成极紧凑的核心；而当重星占一半时，"
                   "重星本身就是星团主体，「重星 vs 轻星」的对比反而没那么悬殊。"
                   "也就是说，少量黑洞造成的核心集中更极端。"),
        ("bullet", "② r_c 的离散度很大。B 组三个实现的末态 r_c 分别是 "
                   "0.156 / 0.031 / 0.172 pc，相差近一个数量级，"
                   "而 r_h 的离散度很小。这说明核心塌缩发生的时刻本身具有随机性，"
                   "因此需要用多个随机实现（本作业每组 3 个）来观察其典型行为。"),
    ], footer="作业2 · 结果与解说")

    r.close()
    print(f"[写入] {OUTPDF}")
    return OUTPDF


def main():
    facts = load_facts()
    if not facts["runs"]:
        print("错误: 未找到任何运行结果，请先运行 run_all.sh 与 analyze.py",
              file=sys.stderr)
        return 1
    build(facts)
    return 0


if __name__ == "__main__":
    sys.exit(main())
