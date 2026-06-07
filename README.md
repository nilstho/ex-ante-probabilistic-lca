# Ex-ante Probabilistic LCA — one Python pipeline (Monte Carlo + GSA on Brightway2)

A single, reproducible **Python** pipeline for **ex-ante probabilistic life cycle
assessment**: parameter presampling, Monte Carlo uncertainty propagation on a parameterised
[Brightway2](https://brightway.dev/) model, and global sensitivity analysis (GSA) — all in
one place.

It consolidates a workflow that originally spanned **three tools and three languages**:

> **R** (presampling) → **Activity Browser** (LCA calculation) → **MATLAB** (GSA)

into a **single Jupyter notebook**, removing the manual file hand-offs between the steps and
making the whole analysis seedable, version-controlled, and reproducible end-to-end.

Built for emerging-technology assessment, where parameters are uncertain *a priori* and you
want to (a) propagate that uncertainty to the impact result and (b) identify which
parameters drive it. It implements the Monte Carlo + GSA steps of the
Safe-and-Sustainable-by-Design framework of **Blanco et al. (2024)** (see
[Acknowledgments](#acknowledgments)), applied to a III-V/Si tandem PV case study, and is
written so it can be adapted to **any** parameterised Brightway2 project — see
[`docs/WORKFLOW.md`](docs/WORKFLOW.md).

---

## From three tools to one

| Step | Originally | Here (one Python pipeline) |
|---|---|---|
| **Presampling** of uncertain parameters | R | Python — `scipy.stats` (PERT, triangular, normal, lognormal, uniform, Bernoulli, beta) |
| **LCA calculation** per scenario | Activity Browser (manual scenario runs + Excel export) | Python — Brightway2 (`bw2calc`) driven in a loop |
| **GSA** (Borgonovo δ) | MATLAB | Python — [SALib](https://salib.readthedocs.io/) |

```
Presampling  ─►  Monte Carlo simulation  ─►  Global Sensitivity Analysis
(sample X)       (evaluate formulas,            (Borgonovo δ via SALib)
                  run LCA per scenario → Y)
```

Why consolidate: one language and environment instead of three; no CSV/Excel hand-offs
between tools; a single fixed random seed makes the entire chain reproducible; and the whole
analysis lives in version control.

### What surfaced while porting the calculation into Python

Reproducing the Activity Browser numbers in Python initially gave different answers. Running
every variant against a clean restore of the project pinned down why — and these findings
are baked into the notebook as guards:

| Problem | Symptom | Handling |
|---|---|---|
| **NumPy 2.0 kernel** | `np.NaN was removed in the NumPy 2.0 release` crash in `lci()` | Use a **numpy < 2** kernel (`bw2data 3.6.x` is incompatible with numpy ≥ 2); Layer 1 warns you |
| **Only one database scanned** | The 4 Si‑supply‑chain activity parameters had no effect; their GSA sensitivity was a spurious zero | Evaluate formula exchanges in **all** relevant databases; a validation gate flags unused parameters |
| **Database contamination** | Re‑running gave drifting results | Overwrite **every** parameterised exchange every iteration → the result becomes state‑independent |
| **AB scenario export not reproducible** | AB produced **different numbers on identical inputs** between runs, matching no consistent parameter application | A further reason to script the calculation: this Python pipeline is deterministic given the sample matrix |

On a clean database the Python pipeline reproduces the original parameter approach
**exactly** (`0.136787, 0.188661, 0.105696, …`) and gives the same answer on every run.

---

## Example results — III‑V/Si tandem PV

Two cases, simple → full (the notebook is organised the same way): first the 4‑parameter module
case, then the complete 26‑factor model. Impact = global warming of 1 kWh of electricity
(ILCD 2.0 climate change total).

### Case A — simple: module performance (4 parameters)

The reduced exercise varying only the generation parameters `Eff_PV`, `PR_PV`, `LT`, `Irrad`.
The score is purely multiplicative in the four, so the δ **ranking is reproducible regardless of
background state** (only the absolute scale is baseline‑dependent): **`LT ≫ Irrad > Eff_PV > PR_PV`**.
A clean starting point — lifetime dominates because it sets the per‑kWh denominator.

![Simplified 4-parameter case](docs/gsa_elec_4param.png)

### Case B — full model (26 factors)

All 21 direct parameters **plus** the five `pₓ` success‑probabilities (the binary design choices'
chances), which enter indirectly via the switches they parameterise (the thesis's full t=0 model).

**Output distribution** (kg CO₂‑eq / kWh): mean **0.200**, median **0.175**, P5 **0.106**,
P95 **0.405** — right‑skewed, with a long tail toward higher impacts.

![Monte Carlo distribution](docs/mc_distribution.png)

**Global sensitivity (Borgonovo δ)** — which factors drive that spread (`*` = indirect `pₓ`):

![GSA delta heatmap](docs/gsa_delta_heatmap.png)

| Rank | Factor | δ | S1 | Meaning |
|---|---|---|---|---|
| 1 | `LT` | 0.220 | 0.174 | Panel lifetime (sets the per‑kWh denominator) |
| 2 | `bin_CuSint` | 0.194 | 0.380 | Cu laser‑vs‑chemical sintering choice (largest first‑order variance share) |
| 3 | `P_movpe_tool` | 0.137 | 0.069 | MOVPE tool power load |
| 4 | `RT_movpe` | 0.126 | 0.045 | MOVPE runtime |
| 5 | `bin_FM` | 0.121 | 0.053 | Front‑metal nanoink choice |
| 6 | `Irrad` | 0.107 | 0.030 | Irradiation |
| 7 | `* pi_FM` | 0.082 | 0.010 | P(Cu vs Ag nanoink) |
| 9 | `* pi_CuSint` | 0.080 | 0.006 | P(laser sintering Cu) |

Lifetime and the discrete copper‑sintering choice dominate; the rest contribute modestly. Full
table: [`examples/gsa_delta_results.csv`](examples/gsa_delta_results.csv).

**Output vs. each input** (scatter screening — the visual companion to δ; a clear trend = influential,
a flat cloud = negligible, two vertical stripes = a binary choice):

![GSA scatter grid](docs/gsa_scatter.png)

### How this compares to the thesis (Blanco 2022, Ch. 6)

**Consistent, not identical.** The impact magnitude (~0.20 kg CO₂‑eq/kWh) matches the thesis's
t=0 snapshot, and the same *set* of factors dominates (lifetime, MOVPE power & runtime, the
front‑metal choice). The **ordering differs** — the thesis t=0 ranks `P_movpe_tool` first and finds
the `pₓ` more influential than the binaries, whereas we get `LT` first and the binaries above their
`pₓ`. This is structural, not a parameter‑selection issue:

- **Functional unit.** Per‑kWh impact ∝ `1/(LT·Irrad·Eff·PR)`, so lifetime/irradiation scale the
  *whole* footprint while `P_movpe_tool` scales only the MOVPE branch — elevating `LT` here.
- **Discrete‑jump size.** `bin_CuSint` injects a large discrete impact (laser‑sintering reagents),
  making the binary dominate and diluting its `pₓ`. The thesis itself notes this binary‑vs‑`pₓ`
  ordering flips between its case studies, i.e. it is model‑specific.
- **Different model.** This runs the reconstructed `SiTaSol_F2v1a`/`IEA_PVPS_2020` databases with
  ILCD 2.0, not the thesis's exact integrated LCA — so exact δ values aren't expected to coincide.

---

## Quick start

```bash
# 1. Create an environment with numpy < 2 (critical!)
conda env create -f environment.yml
conda activate bw2-mc-gsa

# 2. Launch Jupyter and open the notebook
jupyter lab MC_workflow.ipynb
```

### Two notebooks

| Notebook | Use it when |
|---|---|
| **`MC_workflow.ipynb`** (recommended) | Structured/layered. All case-specific settings live in two layers at the top (environment binding + a declarative **parameter registry**); a validation gate auto-checks parameters against the model formulas. Best for adapting to new cases. |
| `MC_workflow_simple.ipynb` | Minimal, everything inline. Good for a quick read of the bare mechanism. |

Run it top to bottom — it goes from **simple to full**:

| Section | Purpose | Edit? |
|---|---|---|
| 1 | Imports, **numpy‑version check** | — |
| 2 | Environment binding — project / databases / FU / method | per project |
| 3 | **Shared engine** — samplers + the Monte‑Carlo loop + formula‑exchange index | never |
| **Part A** | **Simple case: module performance (4 parameters)** — presample, run, MC + scatter, GSA | start here |
| **Part B** | **Full model (26 factors):** registry → presample → load → **validation gate** + introspection → **verification** → full MC run → MC + scatter → GSA heatmap | **registry edited most** |

> **In Part B, always run the verification first.** If the 5 scores match the reference, the full
> run will too. If they don't, check the kernel/numpy version and that every parameterised database
> is listed in `PARAMETERISED_DBS`.

The **validation gate** (Part B) cross-checks your parameter registry against the actual model
formulas and **fails loudly** if a parameter is unresolved or silently unused — preventing the
"I forgot a database, so some parameters did nothing" class of bug. New users should start with the
two self‑contained notebooks in [`examples/`](examples/).

---

## The one principle that makes it reproducible

> **Every parameterised exchange is overwritten from the parameters in every iteration.**

Because nothing is left at a "default", the LCA score depends only on the sampled
parameter vector — **not** on whatever state the database was left in by a previous run.
This is why the workflow is immune to the database contamination that made earlier
attempts (and the AB scenario engine) non‑reproducible.

---

## Repository layout

```
.
├── MC_workflow.ipynb         # structured/layered workflow (recommended)
├── MC_workflow_simple.ipynb  # minimal all-inline variant
├── environment.yml           # conda env (numpy < 2, brightway2, SALib, …)
├── requirements.txt          # pip alternative
├── examples/
│   ├── tea_demo.ipynb         # self-contained teaching example (cup of tea)
│   ├── ev_demo.ipynb          # self-contained teaching example (electric vehicle)
│   ├── README.md              # what the examples teach
│   ├── figs/                  # rendered analysis dashboards
│   ├── PROB_X_t0_short.csv    # 10-scenario sample matrix for a fast verification run
│   └── gsa_delta_results.csv  # GSA δ table from the full 10k run
├── docs/
│   ├── WORKFLOW.md            # how to adapt this to your own Brightway2 project
│   ├── gsa_delta_heatmap.png  # GSA results (figures used in this README)
│   ├── gsa_scatter.png
│   ├── gsa_elec_4param.png
│   └── mc_distribution.png
└── LICENSE
```

The Brightway2 project backup, ecoinvent databases and course materials are **not**
included (licensing + size). Point the notebook at your own project.

### Educational examples

New to the workflow? Start with [`examples/`](examples/) — two **self‑contained** notebooks
(cup of tea, electric vehicle) that run the identical presample → Monte Carlo → GSA machinery
on toy models with **no ecoinvent or project backup required**, each with the full set of
visuals. See [`examples/README.md`](examples/README.md).

---

## Requirements

- Python ≥ 3.9
- **numpy < 2** (hard requirement for `bw2data` 3.6.x)
- `brightway2`, `bw2data`, `bw2calc`, `bw2io`
- `SALib`, `scipy`, `pandas`, `matplotlib`, `seaborn`

See [`environment.yml`](environment.yml).

---

## Acknowledgments

The III‑V/Si tandem photovoltaics case study, the parameterised model, and the underlying
**Safe‑and‑Sustainable‑by‑Design (SSbD) global‑sensitivity framework** are the work of
**Carlos Felipe Blanco** and colleagues at the Institute of Environmental Sciences (CML),
Leiden University. Carlos developed the original analysis across **R** (presampling),
the **Activity Browser** (LCA calculation), and **MATLAB** (GSA).

This repository ports that same method into a **single Python pipeline** on Brightway2 —
the contribution here is the consolidation and reproducibility, not the method itself. Full
credit for the framework, the model, and the case study goes to Carlos and co‑authors.

## References

- Blanco, C. F., Behrens, P., Vijver, M. G., Peijnenburg, W. J. G. M., Quik, J. T. K., &
  Cucurachi, S. (2024). **A framework for guiding safe and sustainable‑by‑design innovation.**
  *Journal of Industrial Ecology, 29*(1).
  https://doi.org/10.1111/jiec.13609
- Blanco Rocha, C. F. (2022). *Guiding safe and sustainable technological innovation under
  uncertainty: a case study of III‑V/silicon photovoltaics* (PhD thesis, Leiden University).
  https://hdl.handle.net/1887/3455392

If you use this workflow in academic work, please cite the framework paper above.

## License

[MIT](LICENSE)
