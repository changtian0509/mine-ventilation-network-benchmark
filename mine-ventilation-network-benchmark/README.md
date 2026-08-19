# Mine Ventilation Network Benchmark Data

Data and results accompanying the manuscript on the *edge resistance-law method*
(resistance-law iteration under dual conservation) for steady-state solving of
mine ventilation networks, including the locally pressure-facility-dominated
negative-resistance recirculation (LPFNR) network class.

This repository provides the network data, the computed solutions, and the
data-generation scripts used in the numerical experiments of the manuscript:
five real-mine ventilation networks (A–E) and a 64-instance parameter sweep of
small fan-dominated recirculation (LPFNR) networks.

Data repository: https://github.com/changtian0509/mine-ventilation-network-benchmark

The original mine data are provided from the research project of the corresponding author, Jinzhang Jia (Liaoning Technical University).

## Repository structure

| Path | Description |
|---|---|
| `engineering_networks/` | Five test networks A–E: topologies and steady-state solutions of three solvers |
| `fan_family_64/`       | 64-instance sweep over the 6-node / 7-edge LPFNR benchmark family |
| `topologies/`          | Network topology figures (SVG/PNG) |
| `raw/`                 | Original mine data files and the processed inputs used to build the paper networks |
| `code/`                | Python scripts used to build the networks and regenerate the published data |

## Solvers

Three ventilation-network solvers are compared throughout:

- `loop` — the loop (mesh) method, global Newton with Hardy–Cross fallback;
- `node` — the node-pressure method, global Newton;
- `dcp` — the edge resistance-law method proposed in the manuscript.

All results were produced with a single solver implementation and common
settings: stopping tolerances (ε_q, ε_h) = (10⁻⁴ m³/s, 0.1 Pa), the max-flow
initial airflow allocation, at most 1000 iterations, single-threaded execution
on Windows 10 / Python 3.7, and timings taken as the best of three runs.

## Engineering networks (`engineering_networks/`)

| Network | Key | Source | n | m | L | Fan edges | Adjust edges | Reversed-flow edges | Q (m³/s) |
|---|---|---|---|---|---|---|---|---|---|
| A | `hongshagang_1` | Hongshagang No. 1 mine | 38 | 46 | 9 | 0 | 6 | 0 | 151.295 |
| B | `hongshagang_difficult` | Hongshagang mine, difficult-ventilation conditions | 113 | 155 | 43 | 0 | 39 | 0 | 154.09 |
| C | `data_mine` | Actual mine network (from `data.json`), 3 local fans added | 363 | 495 | 133 | 3 | 0 | 10 | 100.0 |
| D | `daxing` | Daxing mine, 4 local fans added | 588 | 783 | 196 | 4 | 0 | 184 | 1000.0 |
| E | `jinchuan` | Jinchuan mine | 795 | 1001 | 207 | 7 | 0 | 229 | 1006.8 |

Notation: n = number of nodes, m = number of edges, L = number of independent
loops. `Fan edges` counts fan branches with a quadratic characteristic
h_fan(q) = a0 + a1·q + a2·q²; `Adjust edges` counts branches with an added
regulating resistance; `Reversed-flow edges` counts edges whose solved airflow
is opposite to the assigned orientation, indicating the scale of reverse flow
and recirculation in each network.

Each network folder contains:

- `topology.csv` — columns `edge_id, name, from, to, r_N_s2_m8, adjust_r,
  is_fan, a0_Pa, a1_Pa_s_m3, a2_Pa_s2_m6`; non-fan edges have zero fan
  coefficients;
- `solutions.csv` — columns `method, edge_id, q_m3_s, h_Pa`: steady-state
  airflow and pressure drop of every edge for every solver;
- `pressures.csv` — columns `method, node_id, p_Pa`: steady-state node
  pressure for every solver.

Summary files at the `engineering_networks/` level:

- `summary.csv` — per-network/per-method size parameters, convergence flag,
  iteration count, and build / max-flow / solve timings (best of three);
- `summary_split.csv` — per-network/per-method timing components;
- `timing_stats.csv` — mean and standard deviation of repeated timings;
- `tab5_pairwise_canonical.npz` — the maximum pairwise differences between the
  three solver solutions (max |Δq| in m³/s and max |Δp| in Pa), i.e. the numbers
  reported in the results table of the manuscript (numpy `.npz`).

## LPFNR benchmark family (`fan_family_64/`)

A 6-node / 7-edge base topology with one fan branch and L = 2 independent
loops (Fig. 9 of the manuscript). Total airflow Q = 12 m³/s and the same
stopping tolerances as above. The fan coefficients a0 = 200 Pa and
a1 = 30 Pa·s/m³ are fixed, and the grid a2 ∈ [−0.8, −0.45] Pa·s²/m⁶
(8 levels) × r ∈ [0.1, 0.5] N·s²/m⁸ (8 log-spaced levels) is swept, giving
64 instances whose common max-flow starting point lies in the
negative-stiffness region of the fan characteristic.

- `instances.csv` — the 64 (a2, r) parameter pairs;
- `summary.csv` — per instance and method: convergence flag (1/0) and
  iteration count;
- `edges.csv` / `nodes.csv` — per instance and method: edge airflow and
  pressure drop / node pressures;
- `topology.csv` — the fixed base topology (a2 is scanned per instance).

Convergence summary (PYTHONHASHSEED = 0): loop 30/64, node 44/64, and the edge
resistance-law method 64/64; see `summary.csv` for details.

## Raw inputs (`raw/`)

| Repository file | Original file | Used for |
|---|---|---|
| `hongshagang_1.xls` | 红沙岗一号煤矿-论文.xls | Network A |
| `hongshagang_difficult.xls` | 红沙岗一号煤矿困难通风-论文5.xls | Network B |
| `data_mine_original.json` | data.json | Network C (raw) |
| `data_mine_paper5.json` | data_paper5.json | Network C (3 local fans added) |
| `daxing_original.xls` | daxing.xls | Network D (raw) |
| `daxing_paper5.json` | daxing_paper5.json | Network D (4 local fans added) |
| `jinchuan.xls` | jinchuan.xls | Network E |

## Code (`code/`)

The Python scripts that build the paper networks from the raw inputs and
regenerate the published data (e.g. `build_paper5_nets.py`,
`refresh_eng_data.py`, `refresh_scan.py`). They import the solver
implementation; adjust the `ROOT` / `OUT` paths in each script before running.

Reproducibility note: the node-pressure iterations on the near-critical
fan-family instances depend on Python hash randomization. All published numbers
were generated with

```powershell
$env:PYTHONHASHSEED = "0"; python -X utf8 refresh_scan.py
```

## License

The data and code are made available under the
[Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)
license. See `LICENSE`.

## Citation

If you use this repository, please cite the manuscript:

> [Authors], *[Title]*, [Journal] (in preparation).
> Data: mine-ventilation-network-benchmark, https://github.com/changtian0509/mine-ventilation-network-benchmark (Zenodo DOI to be assigned).
