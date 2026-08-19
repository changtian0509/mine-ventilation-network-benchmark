# -*- coding: utf-8 -*-
"""Recompute steady-state metrics for the 64 scan instances (DCP only)."""
import sys, io, os, contextlib
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
ROOT = r"D:\software_and_algo\software_development\network_solution"
sys.path.insert(0, ROOT); os.chdir(ROOT)
from utils.r_w_xls import xls_file_to_json
from algo.Initial_vent import VentNetwork, build_network
from algo import Vent_dcp
from algo.some_tool import R as Rfn, is_fan_edge
OUT = r"D:\博士\论文\论文9\figs_ch5"

TEMPLATE = xls_file_to_json(os.path.join(ROOT, "data", "dcp优势-7边角联小网络.xls"))
EDGES = [(e["s"], e["t"], e["r"]) for e in TEMPLATE["edges"]]

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

def metrics(dg):
    edges = [e for e in dg.es if e["id"] != "ST" and not e["id"].startswith(("S0","T0"))]
    n = 6  # real nodes
    rows, cols, vals = [], [], []
    for e in edges:
        u, v = e.tuple
        if u != 0 and u < n: rows.append(u-1); cols.append(e.index); vals.append(1.0)
        if v != 0 and v < n: rows.append(v-1); cols.append(e.index); vals.append(-1.0)
    A = sp.coo_matrix((vals,(rows,cols)), shape=(n-1, len(edges))).tocsc()
    d = np.zeros(len(edges)); qvec = np.zeros(len(edges))
    for k, e in enumerate(edges):
        q = e["q"]; qvec[k] = q
        r = Rfn(e)
        if is_fan_edge(e):
            d[k] = 2.0*(r - e["a2"])*abs(q) - e["a1"]
        else:
            d[k] = 2.0*r*abs(q)
    D = sp.diags(d)
    Din = sp.diags(np.where(np.abs(d) > 1e-300, 1.0/d, 0.0))
    Jn = (A @ Din @ A.T).toarray()
    lam_n = float(np.linalg.eigvalsh((Jn+Jn.T)/2).min())
    Z = sla.null_space(A.toarray())  # (m x L) basis of ker A
    Jl = (Z.T @ D @ Z)
    lam_l = float(np.linalg.eigvalsh((Jl+Jl.T)/2).min())
    return qvec, d, lam_n, lam_l

A2 = np.linspace(-0.8, -0.45, 8)
R  = np.geomspace(0.1, 0.5, 8)
lam_n = np.zeros((8,8)); lam_l = np.zeros((8,8)); min_slope = np.zeros((8,8))
qfan = np.zeros((8,8)); slopefan = np.zeros((8,8))
for i,a2 in enumerate(A2):
    for j,r in enumerate(R):
        vnet = VentNetwork(); build_network(make_gd(a2,r), vnet.graph()); vnet.addVirtualST()
        with contextlib.redirect_stdout(io.StringIO()):
            ok,_ = Vent_dcp.vno_dcp_NoadjustNegativeEdges(vnet)
        if not ok:
            lam_n[i,j] = lam_l[i,j] = min_slope[i,j] = float('nan')
            continue
        qvec, d, ln, ll = metrics(vnet.graph())
        lam_n[i,j] = ln; lam_l[i,j] = ll
        min_slope[i,j] = d.min(); qfan[i,j] = qvec[0]; slopefan[i,j] = d[0]

old = np.load(os.path.join(OUT,"fig10_scan.npz"))
np.savez(os.path.join(OUT,"fig10_scan.npz"),
         A2=A2, R=R,
         loop=old["loop"], node=old["node"], dcp=old["dcp"],
         loop_iter=old["loop_iter"], node_iter=old["node_iter"], dcp_iter=old["dcp_iter"],
         lam_n=lam_n, lam_l=lam_l, min_slope=min_slope, qfan=qfan, slopefan=slopefan)
print("min fan slope at steady state: min=%.3f max=%.3f (all positive?)" % (np.nanmin(min_slope), np.nanmax(min_slope)))
print("node Jacobian lam_min: min=%.3f max=%.3f  (PD if >0)" % (np.nanmin(lam_n), np.nanmax(lam_n)))
print("loop Jacobian lam_min: min=%.3f max=%.3f  (PD if >0)" % (np.nanmin(lam_l), np.nanmax(lam_l)))
print("fan q range: %.2f .. %.2f" % (qfan.min(), qfan.max()))
print("fan slope range: %.2f .. %.2f" % (slopefan.min(), slopefan.max()))
print("saved")
