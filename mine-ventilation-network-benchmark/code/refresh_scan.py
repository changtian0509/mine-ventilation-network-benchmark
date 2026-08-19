# -*- coding: utf-8 -*-
"""Re-run Fig.10 scan: 8x8 (a2, r) grid; record success + iteration counts + Jacobian PD checks."""
import sys, io, os, inspect, contextlib
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
ROOT = r"D:\software_and_algo\software_development\network_solution"
sys.path.insert(0, ROOT); os.chdir(ROOT)
from utils.r_w_xls import xls_file_to_json
from algo.Initial_vent import VentNetwork, build_network
from algo import Vent_iterate, Vent_node, Vent_dcp
from algo.Vent_iterate import find_cycle_fast

OUT = r"D:\博士\论文\论文9\figs_ch5"
TEMPLATE = xls_file_to_json(os.path.join(ROOT, "data", "dcp优势-7边角联小网络.xls"))
EDGES = [(e["s"], e["t"], e["r"]) for e in TEMPLATE["edges"]]

# ---------------- count recorders ----------------
NS, NT, HC, DM = [], [], [], []
def inject_all():
    src = inspect.getsource(Vent_node.solve_node_pressures)
    m = "        F, J = build_system(p_root)\n"
    inj = m + "        NS.append(0.0)\n"
    ns = dict(Vent_node.solve_node_pressures.__globals__); ns["NS"] = NS
    exec(compile(src.replace(m, inj, 1), "<inj_node>", "exec"), ns)
    Vent_node.solve_node_pressures = ns["solve_node_pressures"]

    src = inspect.getsource(Vent_iterate._iterate_newton)
    m = "    while k <= count:\n"
    inj = m + "        NT.append(0.0)\n"
    ns = dict(Vent_iterate._iterate_newton.__globals__); ns["NT"] = NT
    exec(compile(src.replace(m, inj, 1), "<inj_newton>", "exec"), ns)
    Vent_iterate._iterate_newton = ns["_iterate_newton"]

    src = inspect.getsource(Vent_iterate.iterate_hardy_cross)
    m = "            # 判断精度是否满足要求(风量精度和阻力精度)\n"
    inj = "            HC.append(0.0)\n" + m
    ns = dict(Vent_iterate.iterate_hardy_cross.__globals__); ns["HC"] = HC
    exec(compile(src.replace(m, inj, 1), "<inj_hc>", "exec"), ns)
    Vent_iterate.iterate_hardy_cross = ns["iterate_hardy_cross"]

    src = inspect.getsource(Vent_dcp.solve_dcp_pressures)
    m = "    for it in range(maxIter):\n"
    inj = m + "        DM.append(0.0)\n"
    ns = dict(Vent_dcp.solve_dcp_pressures.__globals__); ns["DM"] = DM
    exec(compile(src.replace(m, inj, 1), "<inj_dcp>", "exec"), ns)
    Vent_dcp.solve_dcp_pressures = ns["solve_dcp_pressures"]

inject_all()
def clear():
    for L in (NS, NT, HC, DM): del L[:]

def make_gd(a2, r, maxCount=800, a0=200.0, a1=30.0, Q=12.0):
    gd = {"config": [{"vnet_method":"out","q_precise":1e-4,"maxCount":maxCount,"h_precise":0.1,"Q":Q,"solver_method":"dcp"}],
          "nodes":[{"id":i,"x":(i-1)*100.0,"y":0.0,"q_sensor":1,"p_sensor":1} for i in range(1,7)],
          "edges":[],"fans":[{"id":1,"a0":a0,"a1":a1,"a2":a2}]}
    for k,(s,t,r0) in enumerate(EDGES,1):
        rr = r0 if k==1 else r
        gd["edges"].append({"name":"e%d"%k,"id":k,"s":s,"t":t,"r":rr,"adjust_r":0,"q":0,"v":0,"p":0,"adjust_h":0,
            "influence":0.99,"influenced":0.99,"stable":0.99,"rel":0.99,"rel_sen":0,"low_v":0.25,"up_v":25,
            "cs_area":10,"win_area":0,"incl":0,"padra":0.5,"pacra":0.5,"prtara":0.5,"picara":0.5,"pilfba":1,
            "q_sensor":1,"p_sensor":1})
    return gd

def build(gd):
    vnet = VentNetwork(); build_network(gd, vnet.graph()); return vnet

def run(gd, name):
    vnet = build(gd); vnet.addVirtualST()
    with contextlib.redirect_stdout(io.StringIO()):
        if name=="loop":
            ok,_ = Vent_iterate.vno_NoadjustNegativeEdges(vnet)
            it = len(NT) + len(HC)
        elif name=="node":
            ok,_ = Vent_node.vno_node_NoadjustNegativeEdges(vnet)
            it = len(NS)
        else:
            ok,_ = Vent_dcp.vno_dcp_NoadjustNegativeEdges(vnet)
            it = len(DM)
    dg = vnet.graph()
    qs = [abs(e["q"]) for e in dg.es if e["id"] != "ST" and not e["id"].startswith(("S0","T0"))]
    good = bool(ok) and (max(qs) < 1e4)
    return good, it, dg

def steady_state_metrics(dg):
    """At the DCP steady state: fan slope, loop/node Jacobian min eigenvalues."""
    import scipy.sparse as sp
    from algo.some_tool import R as Rfn, is_fan_edge
    edges = [e for e in dg.es if e["id"] != "ST" and not e["id"].startswith(("S0","T0"))]
    n = dg.vcount()
    # reduced incidence A (drop reference node = source 0)
    rows, cols, vals = [], [], []
    for e in edges:
        u, v = e.tuple
        if u != 0: rows.append(u-1); cols.append(e.index); vals.append(1.0)
        if v != 0: rows.append(v-1); cols.append(e.index); vals.append(-1.0)
    A = sp.coo_matrix((vals,(rows,cols)), shape=(n-1, len(edges))).tocsc()
    d = np.zeros(len(edges))
    for k, e in enumerate(edges):
        r = Rfn(e)
        q = e["q"]
        if is_fan_edge(e):
            d[k] = 2.0*(r - e["a2"])*abs(q) - e["a1"]
        else:
            d[k] = 2.0*r*abs(q)
    min_slope = float(np.min(d))
    D = sp.diags(d)
    # node Jacobian A D^{-1} A^T
    try:
        Din = sp.diags(1.0/np.where(np.abs(d) > 1e-300, d, 1e300))
        Jn = (A @ Din @ A.T).toarray()
        lam_n = float(np.linalg.eigvalsh((Jn + Jn.T)/2).min())
    except Exception:
        lam_n = float('nan')
    # loop Jacobian via fundamental loops
    try:
        vnet = VentNetwork()
        dg2 = dg
        # build from current graph (already solved) - construct loop matrix manually via null space of A
        # A (n-1 x m), loop space dimension L = m - (n-1)
        Z = sp.linalg.null_space(A.toarray())
        if Z.shape[1] > 0:
            Jl = (Z.T @ D @ Z).toarray() if hasattr(Z.T @ D @ Z, 'toarray') else (Z.T @ D @ Z)
            Jl = (Z.T @ D @ Z)
            lam_l = float(np.linalg.eigvalsh((Jl + Jl.T)/2).min())
        else:
            lam_l = float('nan')
    except Exception:
        lam_l = float('nan')
    return min_slope, lam_n, lam_l

A2 = np.linspace(-0.8, -0.45, 8)
R  = np.geomspace(0.1, 0.5, 8)
res = {}
counts = {}
min_slope = np.zeros((8,8)); lam_n = np.zeros((8,8)); lam_l = np.zeros((8,8))
for name in ("loop","node","dcp"):
    M = np.zeros((8,8), dtype=int); C = np.zeros((8,8), dtype=int)
    for i,a2 in enumerate(A2):
        for j,r in enumerate(R):
            clear()
            good, it, dg = run(make_gd(a2,r), name)
            M[i,j] = int(good); C[i,j] = it if good else 0
            if name == "dcp" and good:
                ms, ln, ll = steady_state_metrics(dg)
                min_slope[i,j] = ms; lam_n[i,j] = ln; lam_l[i,j] = ll
    res[name]=M; counts[name]=C
    print(name, "rate=%.3f" % M.mean())
    succ = C[M.astype(bool)]
    print("   mean iters (conv)=%.1f  max=%d" % (succ.mean() if len(succ) else float('nan'), succ.max() if len(succ) else 0))

np.savez(os.path.join(OUT,"fig10_scan.npz"), A2=A2, R=R,
         loop=res["loop"], node=res["node"], dcp=res["dcp"],
         loop_iter=counts["loop"], node_iter=counts["node"], dcp_iter=counts["dcp"],
         min_slope=min_slope, lam_n=lam_n, lam_l=lam_l)

both = ((res["loop"]==0) | (res["node"]==0))
neg = min_slope < 0
print("instances where at least one classical method fails: %d/64" % both.sum())
print("instances with negative fan-edge slope at steady state: %d/64" % neg.sum())
print("instances where node Jacobian non-PD: %d/64" % (lam_n <= 0).sum())
print("instances where loop Jacobian non-PD: %d/64" % (lam_l <= 0).sum())
print("saved fig10_scan.npz")
