# -*- coding: utf-8 -*-
"""Unified publication style for all paper figures (Times New Roman, ~9 pt at print size)."""
import matplotlib.pyplot as plt

# ---- method identity (Okabe-Ito, color-blind safe, used in EVERY data figure) ----
C_LOOP = "#E69F00"    # amber     : loop method (global Newton)
C_LOOP_HC = "#D55E00" # vermillion: loop method (Hardy-Cross variant)
C_NODE = "#CC79A7"    # orchid    : node pressure method
C_DCP = "#009E73"     # teal      : edge resistance-law method (hero)

# ---- physical elements (Figs 5-8: fan edge / friction edge / solution) ----
C_FAN = "#B5504B"     # fan edge (muted red)
C_FRIC = "#4A72A8"    # friction edge (muted blue)
C_SOL = "#2C3E50"     # solution marker (dark slate)

# ---- semantic ----
C_OK = "#5E9B61"      # success (scan maps)
C_FAIL = "#C4554D"    # failure (scan maps)
C_TOL = "#4D4D4D"     # tolerance reference line
C_GRAY = "#7F7F7F"    # neutral annotations / reference lines
C_LIGHT = "#BFBFBF"

# ---- pairwise differences (Fig. 13) ----
C_PAIR_LN = "#7C6CCF" # loop-node pair
C_PAIR_LD = "#33B5A5" # loop-dcp pair
C_PAIR_ND = "#E28E2C" # node-dcp pair


def apply(export_w_in, display_w_in, target_pt=9.0, lw=1.2, title_bold=False):
    """Apply the unified rcParams.

    export_w_in : figure width at export (inches)
    display_w_in: width at which the figure is placed in the paper (inches)
    target_pt   : desired on-paper font size (xiaowu = 9 pt)
    Returns the base font size used in the export.
    """
    fs = target_pt * export_w_in / display_w_in
    plt.rcParams.update({
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
        "font.size": fs,
        "axes.labelsize": fs,
        "axes.titlesize": fs,
        "xtick.labelsize": fs,
        "ytick.labelsize": fs,
        "legend.fontsize": fs,
        "lines.linewidth": lw,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 3.2,
        "ytick.major.size": 3.2,
        "axes.grid": False,
        "legend.frameon": False,
    })
    return fs