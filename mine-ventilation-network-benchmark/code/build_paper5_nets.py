# -*- coding: utf-8 -*-
"""Build the two real-mine paper networks with embedded local fans.

Replaces the three synthetic networks (C/D/E) of the previous paper draft with
two real-mine networks:
  - daxing.xls  -> daxing_paper5.json   (4 added local fans; the strong original
                             fan edge is left passive to keep a unique steady state)
  - data.json   -> data_paper5.json     (3 added local fans)

Each added fan is designed at the baseline airflow q0 of its branch so that
a1 > 0, a0 > 0, the operating point lies past the fan-characteristic peak
(q0 > q_peak = a1/(2|a2|)), and the composite (net) slope at q0 is positive.
A negative-resistance edge in daxing (unphysical data artifact) is taken with
abs(r).  All three solvers must converge on the final networks.

Run with:  $env:PYTHONHASHSEED="0"; python -X utf8 build_paper5_nets.py
"""
import sys, io, os, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
ROOT = r"D:\software_and_algo\software_development\network_solution"
sys.path.insert(0, ROOT); os.chdir(ROOT)
from utils.r_w_xls import xls_file_to_json
from algo.Initial_vent import VentNetwork, build_network
from algo import Vent_iterate, Vent_node, Vent_dcp

DATA = os.path.join(ROOT, "data")
SOLVERS = [("loop", Vent_iterate.vno_NoadjustNegativeEdges),
           ("node", Vent_node.vno_node_NoadjustNegativeEdges),
           ("dcp",  Vent_dcp.vno_dcp_NoadjustNegativeEdges)]

def build(gd):
    vnet = VentNetwork(); build_network(gd, vnet.graph()); vnet.addVirtualST()
    return vnet

def solve(vnet, name):
    ret, _ = dict(SOLVERS)[name](vnet)
    return bool(ret)

def edge_flows(vnet):
    dg = vnet.graph()
    out = {}
    for e in dg.es:
        eid = str(e['id'])
        if eid.upper() == 'ST' or eid.startswith('S0') or eid.startswith('T0'):
            continue
        out[eid] = dict(q=float(e['q']), r=float(e['r']), s=int(e.source), t=int(e.target))
    return out

def fan_design(q0, r, net_slope=30.0, boost=1.4, a2=-0.20):
    # B: net-law slope at q0 when a1 = 0 (upper bound on s to keep a1 > 0)
    B = (2.0*r + 2.0*abs(a2))*q0
    # s must exceed these lower bounds so that q0 lies past the fan peak
    # (a1 < 2|a2| q0) and a0 > 0
    S_peak = 2.0*r*q0*1.1
    S_a0 = (2.0 - boost)*r*q0 + abs(a2)*q0
    lo = max(S_peak, S_a0, 0.5)
    hi = max(0.95*B, lo + 0.1)
    s = 0.5*(lo + hi)
    s = min(s, max(net_slope, lo))
    s = max(s, lo)
    a1 = B - s
    a0 = boost*r*q0*q0 - a1*q0 + abs(a2)*q0*q0
    assert a1 > 0 and a0 > 0 and s > 2.0*r*q0, (q0, r, a1, a0, s)
    return dict(a0=round(a0, 4), a1=round(a1, 4), a2=a2)

def pick_edges(flows, fan_ids, n_add, lo_frac, hi_frac=0.95):
    cand = []
    for eid, d in flows.items():
        if eid in fan_ids: continue
        if d['q'] > 0 and d['r'] > 0:
            cand.append((d['q'], d['r'], eid, d['s'], d['t']))
    cand.sort(key=lambda x: -x[0])
    if len(cand) < 2: return []
    # deterministic contiguous window over the flow-ranked candidates
    pool = cand[int(len(cand)*lo_frac): int(len(cand)*hi_frac)]
    picked, used = [], set()
    for c in pool:
        if len(picked) >= n_add: break
        if c[3] in used or c[4] in used: continue
        picked.append(c); used.add(c[3]); used.add(c[4])
    return picked

def diagnose(vnet, tag):
    dg = vnet.graph()
    rev = 0; neg = []
    for e in dg.es:
        eid = str(e['id'])
        if eid.upper() == 'ST' or eid.startswith('S0') or eid.startswith('T0'): continue
        q = float(e['q']); r = float(e['r'])
        a0 = float(e.attributes().get('a0', 0.0) or 0.0)
        a1 = float(e.attributes().get('a1', 0.0) or 0.0)
        a2 = float(e.attributes().get('a2', 0.0) or 0.0)
        is_fan = (a0 != 0 or a1 != 0 or a2 != 0)
        if q < 0: rev += 1
        slope = 2*r*abs(q) - a1 - 2*a2*abs(q) if is_fan else 2*r*abs(q)
        if slope < 1e-9:
            neg.append((eid, q, r, is_fan, slope))
    print('[%s] rev=%d neg_slope_edges=%d' % (tag, rev, len(neg)))
    for x in neg[:8]:
        print('   neg edge %s q=%.3f r=%s is_fan=%s slope=%.3e' % x)
    return rev

def run_network(tag, gd, n_fans_total, keep_orig_fans=True):
    gd = json.loads(json.dumps(gd))
    cfg = gd['config'][0]
    cfg['q_precise'] = 1e-4; cfg['h_precise'] = 0.1
    # unphysical negative-resistance edges -> abs(r)
    for e in gd['edges']:
        if float(e['r']) < 0:
            print('[%s] fixing negative r on edge %s: %s -> %s' % (tag, e['id'], e['r'], abs(float(e['r']))))
            e['r'] = abs(float(e['r']))
    vnet = build(gd)
    ok0 = solve(vnet, 'dcp')
    print('[%s] baseline dcp converged=%s' % (tag, ok0))
    if not ok0: return None
    diagnose(vnet, tag + ' baseline')
    flows = edge_flows(vnet)
    fans0 = gd.get('fans', [])
    n_add = n_fans_total - (len(fans0) if keep_orig_fans else 0)
    fan_ids = set(str(f['id']) for f in (fans0 if keep_orig_fans else []))
    base_fans = list(fans0) if keep_orig_fans else []

    def candidate(lo_frac):
        picked = pick_edges(flows, fan_ids, n_add, lo_frac)
        fans = list(base_fans)
        for q0, r, eid, s, t in picked:
            d = fan_design(q0, r)
            d['id'] = int(eid); fans.append(d)
        gd2 = json.loads(json.dumps(gd)); gd2['fans'] = fans
        return gd2, picked

    def accepted(gd2):
        nt, hc = [], []
        vn_loop = build(gd2)
        ok_l = _solve_loop_instrumented(vn_loop, nt, hc)
        vn_node = build(gd2); ok_n, _ = Vent_node.vno_node_NoadjustNegativeEdges(vn_node)
        vn_dcp = build(gd2); ok_d, _ = Vent_dcp.vno_dcp_NoadjustNegativeEdges(vn_dcp)
        if not (ok_l and ok_n and ok_d):
            return None
        if len(hc) > 0:              # loop fell back to Hardy-Cross
            return None
        qL = _flow_map(vn_loop); qN = _flow_map(vn_node); qD = _flow_map(vn_dcp)
        ids = set(qL) & set(qN) & set(qD)
        dqLD = max(abs(qL[e]-qD[e]) for e in ids)
        dqND = max(abs(qN[e]-qD[e]) for e in ids)
        if dqLD > 5e-5 or dqND > 5e-5:
            return None
        return dict(newton=len(nt), dqLD=dqLD, dqND=dqND, rev=sum(1 for v in qD.values() if v < 0))

    for lo_frac in [0.30, 0.25, 0.35, 0.20, 0.40, 0.15, 0.45, 0.10, 0.05, 0.0, 0.50, 0.55]:
        gd2, picked = candidate(lo_frac)
        print('[%s] try window lo=%.2f: picked %d fan edges' % (tag, lo_frac, len(picked)))
        if len(picked) < n_add:
            print('   WARNING only %d fan edges available' % len(picked))
            continue
        for q0, r, eid, s, t in picked:
            d = fan_design(q0, r)
            print('   +fan edge %s q0=%.3f r=%.4f -> a0=%.2f a1=%.3f a2=%.3f' % (eid, q0, r, d['a0'], d['a1'], d['a2']))
        info = accepted(gd2)
        if info is not None:
            print('[%s] ACCEPTED window lo=%.2f: %s' % (tag, lo_frac, info))
            vnet2 = build(gd2)
            res = {name: solve(vnet2, name) for name, _ in SOLVERS}
            rev = diagnose(vnet2, tag + ' with-fans')
            return gd2, res, rev
        else:
            print('[%s]   rejected window lo=%.2f' % (tag, lo_frac))
    print('[%s] NO acceptable fan configuration found' % tag)
    return None


# ---- loop-method instrumentation helpers ----
LOOP_NT, LOOP_HC = [], []
def _instrument_loop():
    import inspect as _inspect
    src = _inspect.getsource(Vent_iterate._iterate_newton)
    m = '    while k <= count:\n'
    ns = dict(Vent_iterate._iterate_newton.__globals__); ns['LOOP_NT'] = LOOP_NT
    exec(compile(src.replace(m, m + '        LOOP_NT.append(k)\n', 1), '<n>', 'exec'), ns)
    Vent_iterate._iterate_newton = ns['_iterate_newton']
    src = _inspect.getsource(Vent_iterate.iterate_hardy_cross)
    m = '            # 判断精度是否满足要求(风量精度和阻力精度)\n'
    ns = dict(Vent_iterate.iterate_hardy_cross.__globals__); ns['LOOP_HC'] = LOOP_HC
    exec(compile(src.replace(m, '            LOOP_HC.append(k)\n' + m, 1), '<h>', 'exec'), ns)
    Vent_iterate.iterate_hardy_cross = ns['iterate_hardy_cross']
_instrument_loop()
def _solve_loop_instrumented(vnet, nt, hc):
    del LOOP_NT[:]; del LOOP_HC[:]
    ok, _ = Vent_iterate.vno_NoadjustNegativeEdges(vnet)
    nt.extend(LOOP_NT); hc.extend(LOOP_HC)
    return bool(ok)
def _flow_map(vnet):
    return {str(e['id']): float(e['q']) for e in vnet.graph().es}
gd_dax = xls_file_to_json(os.path.join(DATA, 'daxing.xls'))
out = run_network('daxing', gd_dax, n_fans_total=4, keep_orig_fans=False)
if out:
    gd, res, rev = out
    with io.open(os.path.join(DATA, 'daxing_paper5.json'), 'w', encoding='utf-8') as f:
        json.dump(gd, f, ensure_ascii=False, indent=1)
    print('saved daxing_paper5.json')

with io.open(os.path.join(DATA, 'data.json'), encoding='utf-8-sig') as f:
    gd_data = json.load(f)
out = run_network('data.json', gd_data, n_fans_total=3, keep_orig_fans=False)
if out:
    gd, res, rev = out
    with io.open(os.path.join(DATA, 'data_paper5.json'), 'w', encoding='utf-8') as f:
        json.dump(gd, f, ensure_ascii=False, indent=1)
    print('saved data_paper5.json')
