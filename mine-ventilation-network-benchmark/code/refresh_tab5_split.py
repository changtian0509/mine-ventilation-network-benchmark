# -*- coding: utf-8 -*-
"""Split-phase timing for Chapter 5 Table 5: build / max-flow allocation / solve."""
import sys, io, os, time, inspect, csv
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
ROOT = r"D:\software_and_algo\software_development\network_solution"
sys.path.insert(0, ROOT); os.chdir(ROOT)
from utils.r_w_xls import xls_file_to_json
from algo.Initial_vent import VentNetwork, build_network
from algo import Vent_node, Vent_iterate, Vent_dcp

OUT = r"D:\博士\论文\论文9\实例数据\engineering_networks"
os.makedirs(OUT, exist_ok=True)

CASES = [
    ("jinchuan.xls", "Jinchuan", "金川", "jinchuan"),
    ("红沙岗一号煤矿-论文.xls", "Hongshagang I", "红沙岗一号", "hongshagang_1"),
    ("红沙岗一号煤矿困难通风-论文5.xls", "Hongshagang diff", "红沙岗困难", "hongshagang_difficult"),
    ("syn_A_n240.xls", "Synthetic A", "合成A", "synthetic_A"),
    ("syn_B_n430.xls", "Synthetic B", "合成B", "synthetic_B"),
    ("syn_C_n620.xls", "Synthetic C", "合成C", "synthetic_C"),
]
SOLVERS = [("loop", Vent_iterate.vno_NoadjustNegativeEdges),
           ("node", Vent_node.vno_node_NoadjustNegativeEdges),
           ("dcp",  Vent_dcp.vno_dcp_NoadjustNegativeEdges)]

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
def counts():
    return (len(NS), len(NT), len(HC), len(DM))

ORIG_IFQ = Vent_iterate.ifq
FLOW = []
def timed_ifq(vnet):
    t0 = time.perf_counter()
    ret = ORIG_IFQ(vnet)
    FLOW.append(time.perf_counter() - t0)
    return ret
def patch(mode):
    if mode == "timed":
        Vent_iterate.ifq = Vent_node.ifq = Vent_dcp.ifq = timed_ifq
    else:
        Vent_iterate.ifq = Vent_node.ifq = Vent_dcp.ifq = ORIG_IFQ

def build_net(fn):
    t0 = time.perf_counter()
    gd = xls_file_to_json(os.path.join(ROOT, "data", fn))
    vnet = VentNetwork(); build_network(gd, vnet.graph())
    vnet.addVirtualST()
    tb = time.perf_counter() - t0
    return vnet, tb

def run3(fn, solver):
    runs = []
    for _ in range(3):
        clear()
        vnet, tb = build_net(fn)
        patch("timed"); del FLOW[:]
        t0 = time.perf_counter()
        ret, _ = solver(vnet)
        ttot = time.perf_counter() - t0
        patch("restore")
        tf = sum(FLOW)
        ts = ttot - tf
        c = counts()
        runs.append(dict(tb=tb, tf=tf, ts=ts, ttot=ttot, ret=bool(ret), c=c))
    return runs

rows, meta = [], {}
for fn, en, cn, key in CASES:
    per_net, all_build, all_flow = {}, [], []
    for name, solver in SOLVERS:
        runs = run3(fn, solver)
        best_total = min(runs, key=lambda r: r["ttot"])
        c = best_total["c"]
        it = (c[1] + c[2] if name == "loop" else (c[0] if name == "node" else c[3]))
        per_net[name] = dict(tb=min(r["tb"] for r in runs), tf=min(r["tf"] for r in runs),
                             ts=min(r["ts"] for r in runs), ttot=min(r["ttot"] for r in runs),
                             ok=all(r["ret"] for r in runs), it=it)
        all_build += [r["tb"] for r in runs]
        all_flow  += [r["tf"] for r in runs]
    bnet, fnet = min(all_build), min(all_flow)
    meta[key] = dict(cn=cn, bnet=bnet, fnet=fnet)
    for name in ("loop", "node", "dcp"):
        r = per_net[name]
        rows.append(dict(network=key, network_cn=cn, method=name,
                         build_s=bnet, flow_s=fnet, solve_s=r["ts"],
                         total_s=round(bnet + fnet + r["ts"], 6),
                         iterations=r["it"], converged=int(r["ok"])))
        print("%s %-4s build=%.4f flow=%.4f solve=%.4f total=%.4f it=%d" %
              (cn, name, bnet, fnet, r["ts"], bnet + fnet + r["ts"], r["it"]))

with open(os.path.join(OUT, "summary_split.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["network","network_cn","method","build_s","flow_s","solve_s","total_s","iterations","converged"])
    w.writeheader()
    for r in rows:
        w.writerow(r)

# also refresh summary.csv with the same phase times + total
nmap = {"jinchuan":("金川",795,1001,207,1006.8),
        "hongshagang_1":("红沙岗一号",38,46,9,151.295),
        "hongshagang_difficult":("红沙岗困难",113,155,43,154.09),
        "synthetic_A":("合成A",240,312,73,302.4),
        "synthetic_B":("合成B",430,556,127,497.8),
        "synthetic_C":("合成C",620,800,181,703.5)}
with open(os.path.join(OUT, "summary.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["network","network_cn","n","m","L","Q_m3_s","method","converged","iterations","build_s","flow_s","solve_s","time_s_best"])
    for r in rows:
        cn,n,m,L,Q = nmap[r["network"]]
        w.writerow([r["network"], r["network_cn"], n, m, L, Q, r["method"],
                    r["converged"], r["iterations"], r["build_s"], r["flow_s"], r["solve_s"], r["total_s"]])
print("saved summary_split.csv and summary.csv")
print("DONE")
