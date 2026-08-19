# -*- coding: utf-8 -*-
"""Regenerate engineering_networks data (solutions/pressures/topology/summary/pairwise)
for the five real-mine networks of Chapter 5 (networks A/B/C/D/F).

Run with:  $env:PYTHONHASHSEED="0"; python -X utf8 refresh_eng_data.py
"""
import sys, io, os, time, inspect, csv, json
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
ROOT = r"D:\software_and_algo\software_development\network_solution"
sys.path.insert(0, ROOT); os.chdir(ROOT)
from utils.r_w_xls import xls_file_to_json
from algo.Initial_vent import VentNetwork, build_network
from algo import Vent_node, Vent_iterate, Vent_dcp

OUT = r"D:\博士\论文\论文9\实例数据\engineering_networks"
os.makedirs(OUT, exist_ok=True)

# (file, kind, english label, chinese label, key, m_nodes, n_edges, L, Q_m3_s)
CASES = [
    ("红沙岗一号煤矿-论文.xls", "xls", "Hongshagang I", "红沙岗一号",
     "hongshagang_1", 38, 46, 9, 151.295),
    ("红沙岗一号煤矿困难通风-论文5.xls", "xls", "Hongshagang diff", "红沙岗困难",
     "hongshagang_difficult", 113, 155, 43, 154.09),
    ("data_paper5.json", "json", "Network C", "网络C（data.json）",
     "data_mine", 363, 495, 133, 100.0),
    ("daxing_paper5.json", "json", "Daxing", "大兴",
     "daxing", 588, 783, 196, 1000.0),
    ("jinchuan.xls", "xls", "Jinchuan", "金川",
     "jinchuan", 795, 1001, 207, 1006.8),
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

def load_gd(fn, kind):
    if kind == "json":
        with io.open(os.path.join(ROOT, "data", fn), encoding="utf-8-sig") as f:
            return json.load(f)
    return xls_file_to_json(os.path.join(ROOT, "data", fn))

def build_net(fn, kind):
    t0 = time.perf_counter()
    gd = load_gd(fn, kind)
    vnet = VentNetwork(); build_network(gd, vnet.graph())
    vnet.addVirtualST()
    tb = time.perf_counter() - t0
    return vnet, gd, tb

def collect(vnet):
    dg = vnet.graph()
    q, h, p = {}, {}, {}
    for e in dg.es:
        eid = str(e["id"])
        if eid.upper() == "ST" or eid.startswith("S0") or eid.startswith("T0"):
            continue
        q[eid] = float(e["q"])
        h[eid] = float(dg.vs[e.source]["p"] - dg.vs[e.target]["p"])
    for v in dg.vs:
        nid = str(v["id"])
        if nid.upper() in ("S", "T"):
            continue
        p[nid] = float(v["p"])
    return q, h, p

def kcl_real(vnet):
    dg = vnet.graph()
    bal = {}
    for e in dg.es:
        u, v = e.tuple
        bal.setdefault(u, 0.0); bal.setdefault(v, 0.0)
        bal[u] -= e["q"]; bal[v] += e["q"]
    real = [v.index for v in dg.vs if str(v["id"]).upper() not in ("S", "T")]
    return max((abs(bal[i]) for i in real), default=0.0)

def reversed_count(vnet):
    dg = vnet.graph()
    return sum(1 for e in dg.es
               if str(e["id"]).upper() != "ST"
               and not str(e["id"]).startswith(("S0", "T0"))
               and float(e["q"]) < 0)

def run3(fn, kind, solver):
    runs = []
    for _ in range(3):
        clear()
        vnet, gd, tb = build_net(fn, kind)
        patch("timed"); del FLOW[:]
        t0 = time.perf_counter()
        ret, _ = solver(vnet)
        ttot = time.perf_counter() - t0
        patch("restore")
        tf = sum(FLOW)
        ts = ttot - tf
        c = counts()
        runs.append(dict(tb=tb, tf=tf, ts=ts, ttot=ttot, ret=bool(ret), c=c,
                         vnet=vnet, kcl=kcl_real(vnet)))
    return runs

summary_rows, pair = [], {}
for fn, kind, en, cn, key, m_nodes, n_edges, L, Q in CASES:
    per_net = {}
    sols = {}
    rev_dcp = None
    for name, solver in SOLVERS:
        runs = run3(fn, kind, solver)
        best = min(runs, key=lambda r: r["ttot"])
        c = best["c"]
        it = (c[1] + c[2] if name == "loop" else (c[0] if name == "node" else c[3]))
        per_net[name] = dict(tb=min(r["tb"] for r in runs), tf=min(r["tf"] for r in runs),
                             ts=min(r["ts"] for r in runs), ttot=min(r["ttot"] for r in runs),
                             ok=all(r["ret"] for r in runs), it=it, kcl=best["kcl"])
        q, h, p = collect(best["vnet"])
        sols[name] = (q, h, p)
        if name == "dcp":
            rev_dcp = reversed_count(best["vnet"])
        print("%s %-4s build=%.4f flow=%.4f solve=%.4f total=%.4f it=%d kcl=%.2e" %
              (cn, name, per_net[name]["tb"], per_net[name]["tf"], per_net[name]["ts"],
               per_net[name]["tb"] + per_net[name]["tf"] + per_net[name]["ts"], it, per_net[name]["kcl"]))
    bnet = min(per_net[n]["tb"] for n in per_net)
    fnet = min(per_net[n]["tf"] for n in per_net)
    rev = rev_dcp
    print("%s reversed-flow edges at steady state: %d" % (cn, rev))
    for name in ("loop", "node", "dcp"):
        r = per_net[name]
        summary_rows.append(dict(network=key, network_cn=cn, m_nodes=m_nodes, n_edges=n_edges,
                                 L=L, Q_m3_s=Q, reversed_edges=rev,
                                 method=name, converged=int(r["ok"]), iterations=r["it"],
                                 build_s=bnet, flow_s=fnet, solve_s=r["ts"],
                                 time_s_best=round(bnet + fnet + r["ts"], 6)))
    # pairwise
    ids = sorted(set(sols["loop"][0]) & set(sols["node"][0]) & set(sols["dcp"][0]),
                 key=lambda x: (len(x), x))
    nids = sorted(set(sols["loop"][2]) & set(sols["node"][2]) & set(sols["dcp"][2]),
                  key=lambda x: (len(x), x))
    pair[key] = {}
    for a, b in (("loop", "node"), ("loop", "dcp"), ("node", "dcp")):
        qa, ha, pa = sols[a]; qb, hb, pb = sols[b]
        pair[key][(a, b)] = (max(abs(qa[e] - qb[e]) for e in ids),
                             max(abs(pa[nd] - pb[nd]) for nd in nids))
    folder = os.path.join(OUT, key)
    os.makedirs(folder, exist_ok=True)
    # topology
    gd = load_gd(fn, kind)
    fan_ids = {str(f["id"]) for f in gd.get("fans", [])}
    with open(os.path.join(folder, "topology.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["edge_id", "name", "from", "to", "r_N_s2_m8", "adjust_r",
                    "is_fan", "a0_Pa", "a1_Pa_s_m3", "a2_Pa_s2_m6"])
        for e in sorted(gd["edges"], key=lambda e: (len(str(e["id"])), str(e["id"]))):
            eid = str(e["id"])
            fan = next((x for x in gd.get("fans", []) if str(x["id"]) == eid), None)
            w.writerow([eid, e.get("name", "巷道"), str(e["s"]), str(e["t"]),
                        repr(float(e["r"])), repr(float(e.get("adjust_r", 0.0))),
                        1 if fan else 0,
                        fan["a0"] if fan else 0.0, fan["a1"] if fan else 0.0,
                        fan["a2"] if fan else 0.0])
    # solutions + pressures
    with open(os.path.join(folder, "solutions.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["method", "edge_id", "q_m3_s", "h_Pa"])
        for name in ("loop", "node", "dcp"):
            q, h, _ = sols[name]
            for eid in ids:
                w.writerow([name, eid, repr(q[eid]), repr(h[eid])])
    with open(os.path.join(folder, "pressures.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["method", "node_id", "p_Pa"])
        for name in ("loop", "node", "dcp"):
            _, _, p = sols[name]
            for nid in nids:
                w.writerow([name, nid, repr(p[nid])])

with open(os.path.join(OUT, "summary.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["network", "network_cn", "m_nodes", "n_edges", "L",
                                      "Q_m3_s", "reversed_edges",
                                      "method", "converged", "iterations", "build_s",
                                      "flow_s", "solve_s", "time_s_best"])
    w.writeheader()
    for r in summary_rows:
        w.writerow(r)
with open(os.path.join(OUT, "summary_split.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["network", "network_cn", "method", "build_s",
                                      "flow_s", "solve_s", "total_s", "iterations", "converged"])
    w.writeheader()
    for r in summary_rows:
        w.writerow(dict(network=r["network"], network_cn=r["network_cn"], method=r["method"],
                        build_s=r["build_s"], flow_s=r["flow_s"], solve_s=r["solve_s"],
                        total_s=r["time_s_best"], iterations=r["iterations"],
                        converged=r["converged"]))
# pairwise npz (canonical order matches Table 3: A B C D F)
order = ["hongshagang_1", "hongshagang_difficult", "data_mine", "daxing", "jinchuan"]
cn = {"hongshagang_1": "红沙岗一号", "hongshagang_difficult": "红沙岗困难",
      "data_mine": "网络C", "daxing": "大兴", "jinchuan": "金川"}
arr = {k: np.array([pair[key][(a, b)][i] for key in order])
       for k, (a, b), i in (("dq_ln", ("loop", "node"), 0), ("dp_ln", ("loop", "node"), 1),
                            ("dq_ld", ("loop", "dcp"), 0), ("dp_ld", ("loop", "dcp"), 1),
                            ("dq_nd", ("node", "dcp"), 0), ("dp_nd", ("node", "dcp"), 1))}
np.savez(os.path.join(OUT, "tab5_pairwise_canonical.npz"),
         case=np.array([cn[k] for k in order]), **arr)
print("saved summary.csv / summary_split.csv / solutions / pressures / topology / npz")
for key in order:
    d = pair[key]
    print("%-18s dq ln=%.3e ld=%.3e nd=%.3e | dp ln=%.3e ld=%.3e nd=%.3e" % (
        key, d[("loop","node")][0], d[("loop","dcp")][0], d[("node","dcp")][0],
        d[("loop","node")][1], d[("loop","dcp")][1], d[("node","dcp")][1]))
print("DONE")
