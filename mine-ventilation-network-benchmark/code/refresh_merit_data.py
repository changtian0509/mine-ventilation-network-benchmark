# -*- coding: utf-8 -*-
"""Refresh fig9/12/14 merit data: node trace at OUTER iteration level + corrected tolerances."""
import sys, io, os, inspect
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
ROOT = r"D:\software_and_algo\software_development\network_solution"
sys.path.insert(0, ROOT); os.chdir(ROOT)
from utils.r_w_xls import xls_file_to_json
from algo.Initial_vent import VentNetwork, build_network
from algo import Vent_node, Vent_iterate, Vent_dcp

OUT = r"D:\博士\论文\论文9\figs_ch5"
os.makedirs(OUT, exist_ok=True)

NS, NT, HC, DM = [], [], [], []

def inject_all():
    src = inspect.getsource(Vent_node.solve_node_pressures)
    m = "        F, J = build_system(p_root)\n"
    inj = (m +
           "        NS.append(float(np.max(np.abs(F))) if nvar else 0.0)\n")
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

def make_mini_dg():
    rf, a0, a1, a2, R, q0 = 0.05, 300.0, 20.0, -0.5, 1.5, 0.001
    import igraph as ig
    dg = ig.Graph(directed=True)
    dg.add_vertices(2)
    dg.vs[0]["id"]="n1"; dg.vs[1]["id"]="n2"
    dg.add_edge(0,1); dg.add_edge(1,0)
    dg["maxCount"]=800; dg["q_precise"]=1e-2; dg["h_precise"]=5.0
    dg["st_hard"]=False; dg["Q"]=0.0
    fan,fric=dg.es[0],dg.es[1]
    for e in (fan,fric):
        e["id"]="e"; e["r"]=0.0; e["adjust_r"]=0.0; e["adjust_h"]=0.0
        e["a0"]=0.0; e["a1"]=0.0; e["a2"]=0.0; e["fq"]=0.0; e["q"]=q0; e["cs_area"]=1.0
    fan["id"]="e1"; fan["r"]=rf; fan["a0"]=a0; fan["a1"]=a1; fan["a2"]=a2
    fric["id"]="e2"; fric["r"]=R
    return dg

class MiniVnet:
    def __init__(self,dg): self.dg=dg
    def graph(self): return self.dg
    def vSource(self): return 0
    def vAir(self): return -1
    def totalFlow(self): return 0.001

def clear():
    for L in (NS, NT, HC, DM): del L[:]

def build(fn):
    gd = xls_file_to_json(os.path.join(ROOT, "data", fn))
    vnet = VentNetwork(); build_network(gd, vnet.graph()); return vnet

def run_case(fn, tag, outname):
    clear()
    vnet = build(fn); vnet.addVirtualST()
    ok_n = Vent_node.vno_node_NoadjustNegativeEdges(vnet)
    node = np.array([[i, v] for i, v in enumerate(NS)], dtype=float)
    clear()
    vnet = build(fn); vnet.addVirtualST()
    ok_l = Vent_iterate.vno_NoadjustNegativeEdges(vnet)
    nt = np.array([[i, v] for i, v in enumerate(NT)], dtype=float) if NT else np.empty((0,2))
    hc = np.array([[i, v] for i, v in enumerate(HC)], dtype=float) if HC else np.empty((0,2))
    clear()
    vnet = build(fn); vnet.addVirtualST()
    ok_d = Vent_dcp.vno_dcp_NoadjustNegativeEdges(vnet)
    dcp = np.array([[i, v] for i, v in enumerate(DM)], dtype=float) if DM else np.empty((0,2))
    np.savez(os.path.join(OUT, outname), node=node, nt=nt, hc=hc, dcp=dcp)
    print("%s ok(n,l,d)=%s node=%d(末%.3g) nt=%d hc=%d dcp=%d(末%.3g)" % (
        tag, (ok_n, ok_l, ok_d), len(node), node[-1,1] if len(node) else float('nan'),
        len(nt), len(hc), len(dcp), dcp[-1,1] if len(dcp) else float('nan')))

def run_mini():
    clear()
    dg = make_mini_dg()
    ok_n = Vent_node.solve_node_pressures(MiniVnet(dg))
    node = np.array([[i, v] for i, v in enumerate(NS)], dtype=float)
    clear()
    dg = make_mini_dg(); vnet = MiniVnet(dg)
    from algo.Vent_iterate import find_cycle_fast
    cl = find_cycle_fast(vnet, False)
    ok_l = Vent_iterate.iterate(vnet, cl, 40, False, True)
    nt = np.array([[i, v] for i, v in enumerate(NT)], dtype=float) if NT else np.empty((0,2))
    hc = np.array([[i, v] for i, v in enumerate(HC)], dtype=float) if HC else np.empty((0,2))
    clear()
    dg = make_mini_dg()
    ok_d = Vent_dcp.solve_dcp_pressures(MiniVnet(dg))
    dcp = np.array([[i, v] for i, v in enumerate(DM)], dtype=float) if DM else np.empty((0,2))
    np.savez(os.path.join(OUT, "fig9_mini.npz"), node=node, nt=nt, hc=hc, dcp=dcp)
    print("mini ok(n,l,d)=%s node=%d(末%.3g) nt=%d hc=%d dcp=%d(末%.3g)" % (
        (ok_n, ok_l, ok_d), len(node), node[-1,1] if len(node) else float('nan'),
        len(nt), len(hc), len(dcp), dcp[-1,1] if len(dcp) else float('nan')))

inject_all()
run_mini()
run_case("syn_C_n620.xls", "合成C", "fig12_synC.npz")
run_case("jinchuan.xls", "金川", "fig14_jinchuan.npz")
print("DONE")
