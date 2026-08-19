import io, sys
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
d = np.load('fig10_scan.npz')
loop, node, dcp = d['loop'], d['node'], d['dcp']
li, ni, di = d['loop_iter'], d['node_iter'], d['dcp_iter']
print('loop ok=%d/64 rate=%.1f%% mean=%.1f max=%d fail=%d' % (loop.sum(), 100*loop.mean(), li[loop.astype(bool)].mean(), li.max(), (1-loop).sum()))
print('node ok=%d/64 rate=%.1f%% mean=%.1f max=%d fail=%d' % (node.sum(), 100*node.mean(), ni[node.astype(bool)].mean(), ni.max(), (1-node).sum()))
print('dcp  ok=%d/64 rate=%.1f%% mean=%.1f max=%d' % (dcp.sum(), 100*dcp.mean(), di[dcp.astype(bool)].mean(), di.max()))
lf, nf = (1-loop).astype(bool), (1-node).astype(bool)
print('at least one classical fails:', (lf | nf).sum(), '/64')
print('both fail:', (lf & nf).sum())
print('loop fail only:', (lf & ~nf).sum(), '| node fail only:', (~lf & nf).sum())
