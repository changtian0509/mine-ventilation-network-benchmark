# -*- coding: utf-8 -*-
"""Canonical pass for Chapter 5: best-of-3 timings + iteration counts + merit traces (single process)."""
import sys, io, os, inspect, copy, time
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
    inj = m + "        NS.append(float(np.max(np.abs(F))) if nvar else 0.0)\n"
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
    inj = "            HC.append(float(max(DH)) if DH else float('nan'))\n" + m
    ns = dict(Vent_iterate.iterate_hardy_cross.__globals__); ns["HC"] = HC
    exec(compile(src.replace(m, inj, 1), "<inj_hc>", "exec"), ns)
    Vent_iterate.iterate_hardy_cross = ns["iterate_hardy_cross"]

    src = inspect.getsource(Vent_dcp.solve_dcp_pressures)
    m = "    for it in range(maxIter):\n"
    inj = m + "        DM.append(float(merit(r, e)))\n"
    ns = dict(Vent_dcp.solve_dcp_pressures.__globals__); ns["DM"] = DM
    exec(compile(src.replace(m, inj, 1), "<inj_dcp>", "exec"), ns)
    Vent_dcp.solve_dcp_pressures = ns["solve_dcp_pressures"]

inject_all()

def clear():
    for L in (NS, NT, HC, DM): del L[:]

def build(fn):
    gd = xls_file_to_json(os.path.join(ROOT, "data", fn))
    vnet = VentNetwork(); build_network(gd, vnet.graph()); return vnet

SOLVERS = [("loop", Vent_iterate.vno_NoadjustNegativeEdges),
           ("node", Vent_node.vno_node_NoadjustNegativeEdges),
           ("dcp",  Vent_dcp.vno_dcp_NoadjustNegativeEdges)]

def counts():
    return (len(NS), len(NT), len(HC), len(DM))

def run_best3(fn, solver):
    best = None
    for _ in range(3):
        clear()
        vnet = build(fn); vnet.addVirtualST()
        t0 = time.perf_counter()
        ret, _ = solver(vnet)
        dt = time.perf_counter() - t0
        c = counts()
        if best is None or dt < best[0]:
            best = (dt, ret, vnet, c)
    return best

def diff_stats(dgA, dgB):
    maxq = max((abs(dgA.es[e.index]["q"] - dgB.es[e.index]["q"])
                for e in dgA.es if e["id"] != "ST"), default=0.0)
    maxp = max((abs(dgA.vs[v.index]["p"] - dgB.vs[v.index]["p"])
                for v in dgA.vs if v["id"] not in ("S", "T")), default=0.0)
    return maxq, maxp

CASES = [
    ("jinchuan.xls", "Jinchuan", "金川"),
    ("红沙岗一号煤矿-论文.xls", "Hongshagang I", "红沙岗一号"),
    ("红沙岗一号煤矿困难通风-论文5.xls", "Hongshagang diff", "红沙岗困难"),
    ("syn_A_n240.xls", "Synthetic A", "合成A"),
    ("syn_B_n430.xls", "Synthetic B", "合成B"),
    ("syn_C_n620.xls", "Synthetic C", "合成C"),
]

tab5 = []
for fn, en, cn in CASES:
    res = {}
    for name, solver in SOLVERS:
        res[name] = run_best3(fn, solver)
    dref = res["dcp"][2].graph()
    neg = len([e for e in dref.es if e["id"] != "ST" and e["q"] < 0])
    row = {"case": cn, "n": dref.vcount() - 2, "m": dref.ecount() - 1}
    for name in ("loop", "node", "dcp"):
        dt, ret, vnet, c = res[name]
        dg = vnet.graph()
        maxq = maxp = float("nan")
        if name != "dcp" and ret:
            maxq, maxp = diff_stats(dg, dref)
        row[name] = dict(t=dt, ok=bool(ret), it=(c[1] + c[2] if name == "loop" else (c[0] if name == "node" else c[3])),
                         maxq=maxq, maxp=maxp)
    row["neg"] = neg
    tab5.append(row)
    print("%s: n=%d m=%d neg=%d" % (cn, row["n"], row["m"], neg))
    for name in ("loop", "node", "dcp"):
        r = row[name]
        print("   %-4s ok=%s t=%.3fs it=%d maxdq=%s maxdp=%s" %
              (name, r["ok"], r["t"], r["it"],
               ("%.2e" % r["maxq"]) if np.isfinite(r["maxq"]) else "-",
               ("%.2e" % r["maxp"]) if np.isfinite(r["maxp"]) else "-"))

np.savez(os.path.join(OUT, "ch5_tab5.npz"), **{
    "case": np.array([r["case"] for r in tab5]),
    "n": np.array([r["n"] for r in tab5]),
    "m": np.array([r["m"] for r in tab5]),
    "neg": np.array([r["neg"] for r in tab5]),
    "t_loop": np.array([r["loop"]["t"] for r in tab5]),
    "t_node": np.array([r["node"]["t"] for r in tab5]),
    "t_dcp":  np.array([r["dcp"]["t"] for r in tab5]),
    "it_loop": np.array([r["loop"]["it"] for r in tab5]),
    "it_node": np.array([r["node"]["it"] for r in tab5]),
    "it_dcp":  np.array([r["dcp"]["it"] for r in tab5]),
    "dq_loop": np.array([r["loop"]["maxq"] if np.isfinite(r["loop"]["maxq"]) else 0.0 for r in tab5]),
    "dp_loop": np.array([r["loop"]["maxp"] if np.isfinite(r["loop"]["maxp"]) else 0.0 for r in tab5]),
    "dq_node": np.array([r["node"]["maxq"] if np.isfinite(r["node"]["maxq"]) else 0.0 for r in tab5]),
    "dp_node": np.array([r["node"]["maxp"] if np.isfinite(r["node"]["maxp"]) else 0.0 for r in tab5]),
})
print("saved ch5_tab5.npz")
print("DONE-PART1")
