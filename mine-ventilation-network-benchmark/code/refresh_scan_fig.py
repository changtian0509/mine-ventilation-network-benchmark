# -*- coding: utf-8 -*-
"""Fig. 10(b) / fig10_scan.png: clean 8x8 success/failure block maps.

Unified with the chapter-5 figure style (pubstyle): Times New Roman, one
uniform font size, and xiaowu (9 pt) at the in-paper display width (the
7.8 x 3.0 in export is placed at 5.77 in wide in the paper).
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pubstyle as ps

OUT = HERE
NPZ = os.path.join(OUT, "fig10_scan.npz")

def cell_edges(vals, log=False):
    """Cell boundaries so that each of the 8 grid values is the center of its cell."""
    v = np.asarray(vals, dtype=float)
    if log:
        lv = np.log(v)
        d = (lv[-1] - lv[0]) / (len(v) - 1)
        edges = np.exp(np.linspace(lv[0] - d/2, lv[-1] + d/2, len(v) + 1))
    else:
        d = (v[-1] - v[0]) / (len(v) - 1)
        edges = np.linspace(v[0] - d/2, v[-1] + d/2, len(v) + 1)
    return edges

def draw_panel(ax, M, x_edges, y_edges, xlog):
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            color = ps.C_OK if M[i, j] else ps.C_FAIL
            if xlog:
                lx0, lx1 = np.log(x_edges[j]), np.log(x_edges[j+1])
                ix = (lx1 - lx0) * 0.05
                xx0, xx1 = np.exp(lx0 + ix), np.exp(lx1 - ix)
            else:
                ix = (x_edges[j+1] - x_edges[j]) * 0.05
                xx0, xx1 = x_edges[j] + ix, x_edges[j+1] - ix
            iy = (y_edges[i+1] - y_edges[i]) * 0.05
            yy0, yy1 = y_edges[i] + iy, y_edges[i+1] - iy
            ax.add_patch(mpatches.Rectangle((xx0, yy0), xx1 - xx0, yy1 - yy0,
                                            facecolor=color, edgecolor="none"))
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1])

def main():
    # export 7.8 in wide -> placed at 4.84 in wide in the paper => 9 pt (xiaowu)
    fs = ps.apply(7.8, 5.77)
    d = np.load(NPZ)
    A2, R = d["A2"], d["R"]
    methods = [("loop", "Loop method"), ("node", "Node pressure method"),
               ("dcp", "Edge resistance-law method")]
    r_edges = cell_edges(R, log=True)
    a2_edges = cell_edges(A2, log=False)

    fig, axes = plt.subplots(1, 3, figsize=(7.8, 3.0), dpi=300)
    for ax, (key, name) in zip(axes, methods):
        M = np.asarray(d[key])
        draw_panel(ax, M, r_edges, a2_edges, xlog=True)
        ax.set_xscale("log")
        # 0.1 / 0.2 / 0.5 = grid values R[0], R[3], R[7]; short labels avoid
        # overlaps at the uniform 9 pt (on-paper) font size
        ax.set_xticks(R[[0, 3, 7]])
        ax.set_xticklabels(["0.1", "0.2", "0.5"])
        ax.xaxis.set_minor_locator(plt.NullLocator())  # drop log minor ticks entirely
        ax.xaxis.set_minor_formatter(plt.NullFormatter())  # (no 2x10^-1 etc.)
        ax.set_yticks(A2)
        ax.set_yticklabels(["%.2f" % v for v in A2])
        ax.set_xlabel(r"$r$  (N$\cdot$s$^2$/m$^8$)")
        ax.set_ylabel(r"$a_2$  (Pa$\cdot$s$^2$/m$^6$)")
        ax.set_title(name)
        ax.set_aspect("auto")

    legend = [mpatches.Patch(facecolor=ps.C_OK, edgecolor="none", label="Success"),
              mpatches.Patch(facecolor=ps.C_FAIL, edgecolor="none", label="Failure")]
    # Legend centered in the top band; explicit layout so that the full panel
    # titles (e.g. "Edge resistance-law method") are never clipped at the
    # right edge of the figure.
    fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, 1.0),
               ncol=2, frameon=False, columnspacing=1.0, handletextpad=0.5,
               borderaxespad=0.2)
    fig.subplots_adjust(left=0.135, right=0.97, top=0.77, bottom=0.20, wspace=0.42)
    if __name__ == "__main__":
        fig.savefig(os.path.join(OUT, "fig10_scan.png"), dpi=300)
        plt.close(fig)
        print("saved fig10_scan.png (fs=%.2f)" % fs)
    return fig

if __name__ == "__main__":
    main()