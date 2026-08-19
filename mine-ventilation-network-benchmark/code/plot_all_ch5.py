# -*- coding: utf-8 -*-
"""Chapter 5 data figures (Fig. 11/12/13) - unified style, Times New Roman, ~9 pt at print size.

Fig. 11  solving time vs number of edges n (line + marker, full width)
Fig. 12  convergence histories (one panel per network, three method curves each)
Fig. 13  pairwise differences              (grouped bars, full width)
"""
import os, sys, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pubstyle as ps

OUT = HERE
ENG = os.path.normpath(os.path.join(OUT, "..", "实例数据", "engineering_networks"))

# canonical order A B C D F
NODES = {"hongshagang_1": 38, "hongshagang_difficult": 113, "data_mine": 363,
         "daxing": 588, "jinchuan": 795}
EDGES = {"hongshagang_1": 46, "hongshagang_difficult": 155, "data_mine": 495,
         "daxing": 783, "jinchuan": 1001}
NET_ORDER = ["hongshagang_1", "hongshagang_difficult", "data_mine", "daxing", "jinchuan"]
NET_TAG = {"hongshagang_1": "A", "hongshagang_difficult": "B", "data_mine": "C",
           "daxing": "D", "jinchuan": "E"}
NET_EDGES = {k: EDGES[k] for k in NET_ORDER}

# =====================================================================
# Fig. 11 : solving time vs number of edges n (display width 5.49 in)
# =====================================================================
def fig11():
    fs = ps.apply(5.0, 5.49)
    _csv = os.path.join(ENG, "summary_split.csv")
    _solves = {}
    with open(_csv, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            _solves.setdefault(r["network"], {})[r["method"]] = float(r["solve_s"])
    rows = [(NET_TAG[k], EDGES[k], _solves[k]["loop"], _solves[k]["node"],
             _solves[k]["dcp"]) for k in NET_ORDER]
    rows = sorted(rows, key=lambda r: r[1])
    xs = [r[1] for r in rows]
    fig, ax = plt.subplots(figsize=(5.0, 2.8), dpi=300)
    ax.semilogy(xs, [r[2] for r in rows], "o-", color=ps.C_LOOP,
                label="Loop method", ms=4.5, lw=1.4)
    ax.semilogy(xs, [r[3] for r in rows], "s--", color=ps.C_NODE,
                label="Node pressure method", ms=4.5, lw=1.4)
    ax.semilogy(xs, [r[4] for r in rows], "^-.", color=ps.C_DCP,
                label="Edge resistance-law method", ms=5.0, lw=1.4)
    ax.set_xlim(0, 1100)
    ax.set_xticks([0, 200, 400, 600, 800, 1000, 1100])
    ax.set_ylim(5e-3, 0.5)
    ax.set_yticks([0.01, 0.1])
    ax.set_xlabel("Number of edges $n$")
    ax.set_ylabel("Solving time (s)")
    for x, y, tag in zip(xs, [max(r[2], r[3], r[4]) for r in rows],
                         [r[0] for r in rows]):
        ax.annotate(tag, (x, y), textcoords="offset points",
                    xytext=(0, 7), ha="center", va="bottom",
                    fontsize=fs, fontname="Times New Roman")
    fig.tight_layout(rect=[0, 0, 1, 0.80])
    fig.legend(loc="upper center", bbox_to_anchor=(0.5, 0.965), ncol=3,
               frameon=False, fontsize=fs, columnspacing=1.0, handletextpad=0.4)
    fig.savefig(os.path.join(OUT, "fig11_scaling.png"), dpi=400)
    plt.close(fig)
    print("saved fig11_scaling.png (fs=%.2f)" % fs)

# =====================================================================
# Fig. 12 : 3D point-line convergence histories (one panel per method)
# =====================================================================
def fig12():
    """Fig. 12: convergence histories, one panel per network (A-E), three method
    curves per panel (tolerance-normalized residual vs iteration, log-y)."""
    fs = ps.apply(7.6, 5.77)
    eps_h, eps_q = 0.1, 1e-4
    conv = {tag: np.load(os.path.join(OUT, "conv_%s.npz" % tag))
            for tag in ("A", "B", "C", "D", "E")}

    def series(tag, method):
        d = conv[tag]
        if method == "loop":
            arr = d["nt"]; y = np.maximum(arr[:, 1], 1e-12) / eps_h if len(arr) else None
        elif method == "node":
            arr = d["node"]; y = np.maximum(arr[:, 1], 1e-12) / eps_q if len(arr) else None
        else:
            arr = d["dcp"]; y = np.maximum(arr[:, 1], 1e-12) if len(arr) else None
        if y is None or len(y) == 0:
            return None, None
        return arr[:, 0], y

    nets = ("A", "B", "C", "D", "E")
    mets = [("loop", "Loop method", ps.C_LOOP, "o"),
            ("node", "Node pressure method", ps.C_NODE, "s"),
            ("dcp", "Edge resistance-law method", ps.C_DCP, "^")]

    fig = plt.figure(figsize=(7.6, 4.2), dpi=300)
    # manual layout: top row A B C, bottom row D E centered
    panel_w, panel_h, gap = 0.285, 0.305, 0.026
    left = 0.075
    top_y, bot_y = 0.565, 0.100
    xA, xB, xC = left, left + panel_w + gap, left + 2 * (panel_w + gap)
    tot_w = 2 * panel_w + gap
    bx0 = 0.5 - tot_w / 2.0
    pos = {"A": (xA, top_y), "B": (xB, top_y), "C": (xC, top_y),
           "D": (bx0, bot_y), "E": (bx0 + panel_w + gap, bot_y)}
    axes = {}
    for tag in nets:
        ax = fig.add_axes([pos[tag][0], pos[tag][1], panel_w, panel_h])
        axes[tag] = ax
    # legend row on top
    leg_y = 0.955
    for tag in nets:
        ax = axes[tag]
        for m, _, col, mk in mets:
            xs, ys = series(tag, m)
            if xs is None:
                continue
            ax.semilogy(xs, ys, color=col, marker=mk, ms=3.2, lw=1.15)
        ax.axhline(1.0, color=ps.C_TOL, ls="--", lw=0.8)
        ax.set_title("Network %s" % tag, fontsize=fs, pad=2)
        ax.set_xlim(0, 40)
        ax.set_xticks([0, 10, 20, 30, 40])
        ax.set_ylim(1e-6, 1e10)
        ax.set_yticks([1e-6, 1e-2, 1e2, 1e6, 1e10])
        ax.set_yticklabels([r"$10^{-6}$", r"$10^{-2}$", r"$10^{2}$",
                            r"$10^{6}$", r"$10^{10}$"])
        ax.tick_params(labelsize=fs - 0.4)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, p: "%d" % v))
        if tag not in ("A", "D"):
            ax.set_yticklabels([])
        if tag in ("D", "E"):
            ax.set_xlabel("Iteration $k$", fontsize=fs, labelpad=1)
        if tag in ("A", "D"):
            ax.set_ylabel("Tolerance-normalized residual", fontsize=fs, labelpad=1)
    handles = [plt.Line2D([], [], color=col, marker=mk, ms=4.5, lw=1.2, label=lab)
               for m, lab, col, mk in mets]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.990),
               ncol=3, frameon=False, fontsize=fs, columnspacing=1.1,
               handletextpad=0.4)
    fig.savefig(os.path.join(OUT, "fig12_conv.png"), dpi=300)
    plt.close(fig)
    print("saved fig12_conv.png (fs=%.2f)" % fs)
CASE_OF = {"hongshagang_1": "红沙岗一号",
           "hongshagang_difficult": "红沙岗困难",
           "data_mine": "网络C",
           "daxing": "大兴",
           "jinchuan": "金川"}

def pairwise_all(folder):
    """Max pairwise |dq| / |dp| of one network from tab5_pairwise_canonical.npz,
    as {(a, b): (dq, dp)} for the three method pairs."""
    d = np.load(os.path.join(ENG, "tab5_pairwise_canonical.npz"), allow_pickle=True)
    cases = list(d["case"])
    name = CASE_OF[folder]
    if name not in cases:
        raise KeyError("case %s not in tab5_pairwise_canonical.npz" % name)
    i = cases.index(name)
    return {("loop", "node"): (float(d["dq_ln"][i]), float(d["dp_ln"][i])),
            ("loop", "dcp"): (float(d["dq_ld"][i]), float(d["dp_ld"][i])),
            ("node", "dcp"): (float(d["dq_nd"][i]), float(d["dp_nd"][i]))}


def fig13():
    fs = ps.apply(5.2, 5.0)
    nets = [(k, NET_TAG[k]) for k in NET_ORDER]
    pairs = [("loop", "node"), ("loop", "dcp"), ("node", "dcp")]
    pair_colors = [ps.C_PAIR_LN, ps.C_PAIR_LD, ps.C_PAIR_ND]
    pair_labels = ["loop\u2013node", "loop\u2013edge", "node\u2013edge"]
    data = {tag: pairwise_all(folder) for folder, tag in nets}
    fig, axes = plt.subplots(1, 2, figsize=(5.2, 2.2), dpi=300)
    xs = np.arange(len(nets))
    bw = 0.26
    for pi, (a, b) in enumerate(pairs):
        dq = [data[t][(a, b)][0] for _, t in nets]
        dp = [data[t][(a, b)][1] for _, t in nets]
        off = (pi - 1) * bw
        axes[0].bar(xs + off, dq, width=bw, color=pair_colors[pi], label=pair_labels[pi],
                    edgecolor="none", log=True)
        axes[1].bar(xs + off, dp, width=bw, color=pair_colors[pi], label=pair_labels[pi],
                    edgecolor="none", log=True)
    for ax, eps, ylim, yticks, ylab in [
        (axes[0], 1e-4, (1e-7, 1e-2), [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2],
         "max |$\\Delta q$| (m$^3$/s)"),
        (axes[1], 0.1, (1e-7, 1e0), [1e-7, 1e-5, 1e-3, 1e-1], "max |$\\Delta p$| (Pa)"),
    ]:
        ax.axhline(eps, color=ps.C_TOL, ls="--", lw=0.9)
        ax.text(4.32, eps, "$\\varepsilon_q$" if eps < 1 else "$\\varepsilon_h$",
                va="center", ha="left", fontsize=fs)
        ax.set_xticks(xs)
        ax.set_xticklabels([t for _, t in nets])
        ax.set_ylim(*ylim)
        ax.set_yticks(yticks)
        ax.set_xlabel("Network")
        ax.set_ylabel(ylab)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in pair_colors]
    fig.legend(handles, pair_labels, loc="upper center", bbox_to_anchor=(0.5, 0.99),
               ncol=3, frameon=False, fontsize=fs, columnspacing=1.2, handletextpad=0.4)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(os.path.join(OUT, "fig13_pairwise.png"), dpi=300)
    plt.close(fig)
    print("saved fig13_pairwise.png (fs=%.2f)" % fs)


if __name__ == "__main__":
    fig11()
    fig12()
    fig13()
