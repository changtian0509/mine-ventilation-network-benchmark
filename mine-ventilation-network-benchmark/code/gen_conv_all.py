# -*- coding: utf-8 -*-
"""Generate convergence-history npz files (conv_A..conv_E) for the five Chapter-5
networks (A=hongshagang I, B=hongshagang diff, C=data.json, D=daxing, E=jinchuan).

Run with:  $env:PYTHONHASHSEED="0"; python -X utf8 gen_conv_all.py
"""
import pathlib, sys, io, os, inspect, json
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
ROOT = r"D:\software_and_algo\software_development\network_solution"
sys.path.insert(0, ROOT); os.chdir(ROOT)
from utils.r_w_xls import xls_file_to_json
from algo.Initial_vent import VentNetwork, build_network
from algo import Vent_node, Vent_iterate, Vent_dcp

OUT = r"D:\博士\论文\论文9\figs_ch5"
NS, NT, HC, DM = [], [], [], []

def inject_all():
    src = inspect.getsource(Vent_node.solve_node_pressures)
    m = "        F, J = build_system(p_root)\n"
    inj = (m + "        NS.append(float(np.max(np.abs(F))) if nvar else 0.0)\n")
    ns = dict(Vent_node.solve_node_pressures.__globals__); ns["NS"] = NS
    exec(compile(src.replace(m, inj, 1), "<inj_node>", "exec"), ns)
    Vent_node.solve_node_pressures = ns["solve_node_pressures"]

    src = inspect.getsource(Vent_iterate._iterate_newton)
    m = "    while k <= count:\n"
    inj = m + "        NT.append(float(np.max(np.abs(F))) if L else 0.0)\n"
    ns = dict(Vent_iterate._iterate_newton.__globals__); ns["NT"] = NT
    exec(compile(src.replace(m, inj, 1), "<inj_newton>", "exec"), ns)
    Vent_iterate._iterate_newton = ns["_iterate_newton"]

    src = inspect.getsource(Vent_iterate.iterate_hardy_cross)
    m = "            # 判断精度是否满足要求(风量精度和阻力精度)\n"
    inj = ("            HC.append(float(max(DH)) if DH else float('nan'))\n" + m)
    ns = dict(Vent_iterate.iterate_hardy_cross.__globals__); ns["HC"] = HC
    exec(compile(src.replace(m, inj, 1), "<inj_hc>", "exec"), ns)
    Vent_iterate.iterate_hardy_cross = ns["iterate_hardy_cross"]

    src = inspect.getsource(Vent_dcp.solve_dcp_pressures)
    m = "    for it in range(maxIter):\n"
    inj = m + "        DM.append(float(merit(r, e)))\n"
    ns = dict(Vent_dcp.solve_dcp_pressures.__globals__); ns["DM"] = DM
    exec(compile(src.replace(m, inj, 1), "<inj_dcp>", "exec"), ns)
    Vent_dcp.solve_dcp_pressures = ns["solve_dcp_pressures"]

def clear():
    for L in (NS, NT, HC, DM): del L[:]

def load_gd(fn, kind):
    if kind == "json":
        with io.open(os.path.join(ROOT, "data", fn), encoding="utf-8-sig") as f:
            return json.load(f)
    return xls_file_to_json(os.path.join(ROOT, "data", fn))

def build(fn, kind):
    gd = load_gd(fn, kind)
    vnet = VentNetwork(); build_network(gd, vnet.graph()); return vnet

def run_case(fn, kind, tag, outname):
    clear()
    vnet = build(fn, kind); vnet.addVirtualST()
    ok_n = Vent_node.vno_node_NoadjustNegativeEdges(vnet)
    node = np.array([[i, v] for i, v in enumerate(NS)], dtype=float)
    clear()
    vnet = build(fn, kind); vnet.addVirtualST()
    ok_l = Vent_iterate.vno_NoadjustNegativeEdges(vnet)
    nt = np.array([[i, v] for i, v in enumerate(NT)], dtype=float) if NT else np.empty((0,2))
    hc = np.array([[i, v] for i, v in enumerate(HC)], dtype=float) if HC else np.empty((0,2))
    clear()
    vnet = build(fn, kind); vnet.addVirtualST()
    ok_d = Vent_dcp.vno_dcp_NoadjustNegativeEdges(vnet)
    dcp = np.array([[i, v] for i, v in enumerate(DM)], dtype=float) if DM else np.empty((0,2))
    np.savez(os.path.join(OUT, outname), node=node, nt=nt, hc=hc, dcp=dcp)
    print("%s ok(n,l,d)=%s node=%d(末%.3g) nt=%d(末%.3g) hc=%d dcp=%d(末%.3g)" % (
        tag, (ok_n, ok_l, ok_d), len(node), node[-1,1] if len(node) else float('nan'),
        len(nt), nt[-1,1] if len(nt) else float('nan'), len(hc), len(dcp),
        dcp[-1,1] if len(dcp) else float('nan')))

inject_all()
cases = [
    ("红沙岗一号煤矿-论文.xls", "xls", "A", "conv_A.npz"),
    ("红沙岗一号煤矿困难通风-论文5.xls", "xls", "B", "conv_B.npz"),
    ("data_paper5.json", "json", "C", "conv_C.npz"),
    ("daxing_paper5.json", "json", "D", "conv_D.npz"),
    ("jinchuan.xls", "xls", "E", "conv_E.npz"),
]
for fn, kind, tag, out in cases:
    run_case(fn, kind, tag, out)
print("DONE")
