#!/usr/bin/env python3
"""
由 PeTar 输入文件计算星团特征时间，用于自动确定积分时长。

用法:  python3 calc_trh.py <petar_input_file> [T_RH_FACTOR]

输出（一行，空格分隔）:
    N mbar Rh sigma3D trh tend

其中
    t_rh = 0.138 * sqrt(N * R_h^3 / (G * mbar)) / ln(0.11 N)   (Spitzer 半质量弛豫时间)
    t_end = T_RH_FACTOR * t_rh，四舍五入到 0.5 Myr（T_RH_FACTOR 默认 2.0）

R_h 用质量加权的三维半质量半径（相对质量中心）。
"""

import sys

import numpy as np

G_AU = 0.00449830997959438          # pc^3 / (Msun * Myr^2)


def main():
    if len(sys.argv) < 2:
        print("用法: python3 calc_trh.py <petar_input_file> [T_RH_FACTOR]",
              file=sys.stderr)
        return 1
    path = sys.argv[1]
    factor = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0

    d = np.loadtxt(path, skiprows=1)
    if d.ndim == 1:
        d = d[None, :]
    m, pos, vel = d[:, 0], d[:, 1:4], d[:, 4:7]

    N = m.size
    M = m.sum()
    mbar = M / N

    # 质量加权的三维半质量半径（相对质量中心）
    com = (m[:, None] * pos).sum(axis=0) / M
    r = np.linalg.norm(pos - com, axis=1)
    o = np.argsort(r)
    rh = float(np.interp(0.5 * M, np.cumsum(m[o]), r[o]))

    sigma3 = float(np.sqrt((m * (vel ** 2).sum(axis=1)).sum() / M))
    trh = 0.138 * np.sqrt(N * rh ** 3 / (G_AU * mbar)) / np.log(0.11 * N)
    tend = round(factor * trh * 2.0) / 2.0     # 取到 0.5 Myr

    print(f"{N:d} {mbar:.4f} {rh:.4f} {sigma3:.4f} {trh:.3f} {tend:.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
