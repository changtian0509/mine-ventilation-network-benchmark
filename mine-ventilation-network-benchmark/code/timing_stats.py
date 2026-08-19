# -*- coding: utf-8 -*-
import sys, io, os, time, inspect, gc, json, csv
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
ROOT = r"D:\software_and_algo\software_development\network_solution"
sys.path.insert(0, ROOT); os.chdir(ROOT)
from utils.r_w_xls import xls_file_to_json
from algo.Initial_vent import VentNetwork, build_network
from algo import Vent_node, Vent_iterate, Vent_dcp
NS, NT, HC, DM = [], [], [], []
def inject_all():
    src = inspect.getsource(Vent_node.solve_node_pressures)
    m = "        F, J = build_system(p_root)\n"
    ns = dict(Vent_node.solve_node_pressures.__globals__); ns["NS"] = NS
    exec(compile(src.replace(m, m + "        NS.append(float(np.max(np.abs(F))) if nvar else 0.0)\n", 1), "<inj>", "exec"), ns)
    Vent_node.solve_node_pressures = ns["solve_node_pressures"]
    src = inspect.getsource(Vent_iterate._iterate_newton)
    m = "    while k <= count:\n"
    ns = dict(Vent_iterate._iterate_newton.__globals__); ns["NT"] = NT
    exec(compile(src.replace(m, m + "        NT.append(float(np.max(np.abs(F))) if L else 0.0)\n", 1), "<inj>", "exec"), ns)
    Vent_iterate._iterate_newton = ns["_iterate_newton"]
    src = inspect.getsource(Vent_iterate.iterate_hardy_cross)
    m = "            # 判断精度是否满足要求(风量精度和阻力精度)\n"
    ns = dict(Vent_iterate.iterate_hardy_cross.__globals__); ns["HC"] = HC
    exec(compile(src.replace(m, "            HC.append(float(max(DH)) if DH else float('nan'))\n" + m, 1), "<inj>", "exec"), ns)
    Vent_iterate.iterate_hardy_cross = ns["iterate_hardy_cross"]
    src = inspect.getsource(Vent_dcp.solve_dcp_pressures)
    m = "    for it in range(maxIter):\n"
    ns = dict(Vent_dcp.solve_dcp_pressures.__globals__); ns["DM"] = DM
    exec(compile(src.replace(m, m + "        DM.append(float(merit(r, e)))\n", 1), "<inj>", "exec"), ns)
    Vent_dcp.solve_dcp_pressures = ns["solve_dcp_pressures"]
inject_all()
def clear(): [x.clear() for x in (NS, NT, HC, DM)]
def counts(): return (len(NS), len(NT), len(HC), len(DM))
ORIG_IFQ = Vent_iterate.ifq
FLOW = []
def timed_ifq(vnet):
    t0 = time.perf_counter(); ret = ORIG_IFQ(vnet); FLOW.append(time.perf_counter() - t0); return ret
def patch(on):
    Vent_iterate.ifq = Vent_node.ifq = Vent_dcp.ifq = (timed_ifq if on else ORIG_IFQ)
CASES = [
    ("红沙岗一号煤矿-论文.xls", "xls", "A", "hongshagang_1"),
    ("红沙岗一号煤矿困难通风-论文5.xls", "xls", "B", "hongshagang_difficult"),
    ("data_paper5.json", "json", "C", "data_mine"),
    ("daxing_paper5.json", "json", "D", "daxing"),
    ("jinchuan.xls", "xls", "E", "jinchuan"),
]
SOLVERS = [("loop", Vent_iterate.vno_NoadjustNegativeEdges, True),
           ("node", Vent_node.vno_node_NoadjustNegativeEdges, False),
           ("dcp",  Vent_dcp.vno_dcp_NoadjustNegativeEdges, False)]
def load_gd(fn, kind):
    if kind == "json":
        with io.open(os.path.join(ROOT, "data", fn), encoding="utf-8-sig") as f:
            return json.load(f)
    return xls_file_to_json(os.path.join(ROOT, "data", fn))
N_ROUNDS = 10
out_rows = []
for fn, kind, lab, key in CASES:
    gd0 = load_gd(fn, kind)
    acc = {name: dict(tb=[], tf=[], ts=[], it=None, ok=True) for name, _, _ in SOLVERS}
    # warm-up all solvers on this network
    for name, solver, is_loop in SOLVERS:
        for _ in range(2):
            clear()
            vnet = VentNetwork(); build_network(gd0, vnet.graph()); vnet.addVirtualST()
            patch(True); del FLOW[:]
            solver(vnet); patch(False)
    gc.collect()
    order = [s[0] for s in SOLVERS]
    for rnd in range(N_ROUNDS):
        rotated = order[rnd % 3:] + order[:rnd % 3]
        for name in rotated:
            solver = next(s[1] for s in SOLVERS if s[0] == name)
            is_loop = next(s[2] for s in SOLVERS if s[0] == name)
            clear()
            t0 = time.perf_counter()
            vnet = VentNetwork(); build_network(gd0, vnet.graph()); vnet.addVirtualST()
            tb = time.perf_counter() - t0
            patch(True); del FLOW[:]
            t1 = time.perf_counter()
            ret, _ = solver(vnet)
            ttot = time.perf_counter() - t1
            patch(False)
            c = counts()
            it = (c[1] + c[2]) if is_loop else (c[0] if name == "node" else c[3])
            a = acc[name]
            a["tb"].append(tb); a["tf"].append(sum(FLOW)); a["ts"].append(ttot - sum(FLOW))
            if a["it"] is None: a["it"] = it
            a["ok"] = a["ok"] and bool(ret)
    for name, _, _ in SOLVERS:
        a = acc[name]
        r = dict(network=key, method=name,
                 tb_mean=float(np.mean(a["tb"])), tb_std=float(np.std(a["tb"])),
                 tf_mean=float(np.mean(a["tf"])), tf_std=float(np.std(a["tf"])),
                 ts_mean=float(np.mean(a["ts"])), ts_std=float(np.std(a["ts"])),
                 it=a["it"], ok=a["ok"])
        out_rows.append(r)
        print("%s %-4s build=%.4f±%.4f flow=%.4f±%.4f solve=%.4f±%.4f it=%d ok=%s" % (
            lab, name, r["tb_mean"], r["tb_std"], r["tf_mean"], r["tf_std"],
            r["ts_mean"], r["ts_std"], r["it"], r["ok"]))
OUT = r"D:\博士\论文\论文9\实例数据\engineering_networks"
with open(os.path.join(OUT, "timing_stats.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
    w.writeheader()
    for r in out_rows: w.writerow(r)
print("saved timing_stats.csv (N_ROUNDS=%d, interleaved)" % N_ROUNDS)
