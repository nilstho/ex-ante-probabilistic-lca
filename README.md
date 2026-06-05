# Ex-ante Probabilistic LCA — Monte Carlo & Global Sensitivity Analysis in Brightway2

A reproducible, scriptable workflow for **ex-ante probabilistic life cycle assessment**:
parameter-based uncertainty propagation and global sensitivity analysis (GSA) on a
parameterised [Brightway2](https://brightway.dev/) LCA model — **entirely in Python, with
no dependence on the Activity Browser (AB) scenario engine**.

Built for emerging-technology assessment, where parameters are uncertain *a priori* and
you want to (a) propagate that uncertainty to the impact result and (b) identify which
parameters drive it.

It was built for a III‑V/Si tandem photovoltaics case study but is written so it can be
adapted to **any** parameterised Brightway2 project. See
[`docs/WORKFLOW.md`](docs/WORKFLOW.md) for the generalisation guide.

---

## Why this exists

The original course workflow ran the uncertainty analysis through the Activity Browser's
*scenario analysis* feature and exported the results to Excel. Reproducing those numbers
in Python kept giving different answers. After running every variant against a clean
restore of the project, the root causes turned out to be:

| Problem | Symptom | Fix |
|---|---|---|
| **NumPy 2.0 kernel** | `np.NaN was removed in the NumPy 2.0 release` crash in `lci()` | Run with a **numpy < 2** kernel (`bw2data 3.6.x` is incompatible with numpy ≥ 2) |
| **Only one database scanned** | The 4 Si‑supply‑chain activity parameters had no effect; their GSA sensitivity was a spurious zero | Evaluate formula exchanges in **all** relevant databases, not just the foreground |
| **Database contamination** | Re‑running gave drifting results | Overwrite **every** parameterised exchange every iteration → the result becomes state‑independent |
| **AB scenario export is itself unreliable** | AB produced **different numbers on identical inputs** between runs, matching no consistent parameter application | Use this Python workflow as the reproducible reference instead |

The key empirical finding: on a clean database this workflow reproduces the original
parameter approach **exactly** (`0.136787, 0.188661, 0.105696, …`), whereas the AB
scenario export drifted run‑to‑run. **The Python workflow is the trustworthy reference.**

---

## What it does

```
Presampling  ─►  Monte Carlo simulation  ─►  Global Sensitivity Analysis
(sample X)       (evaluate formulas,            (Borgonovo δ via SALib)
                  run LCA per scenario → Y)
```

1. **Presampling** — draw `N` samples for each uncertain parameter from its distribution
   (PERT, triangular, normal, lognormal, uniform, Bernoulli, beta) and save to
   `PROB_X_t0.csv`.
2. **Monte Carlo** — for each scenario, write the sampled values into the parameterised
   exchange amounts (via their stored `formula`), run the LCA, and collect the score `Y`.
3. **GSA** — compute the Borgonovo δ (delta) sensitivity measure with
   [SALib](https://salib.readthedocs.io/) to rank which parameters drive the output
   uncertainty.

---

## Quick start

```bash
# 1. Create an environment with numpy < 2 (critical!)
conda env create -f environment.yml
conda activate bw2-mc-gsa

# 2. Launch Jupyter and open the notebook
jupyter lab MC_workflow.ipynb
```

Then run the cells top to bottom:

| Section | Purpose |
|---|---|
| 0 | Imports, **numpy‑version check**, project/database/method setup |
| 1 | Presampling → `PROB_X_t0.csv` (skip if you already have it) |
| 2 | Load the sample matrix `X` |
| 3 | Pre‑index every formula exchange in **all** databases |
| 4 | **Verification** — 5 scenarios vs the reproducible reference (must show `Diff ≈ 0`) |
| 5 | Full Monte Carlo run (10 000 iterations) → `Model_results_python.csv` |
| 6 | Histogram + boxplot of the output distribution |
| 7 | GSA (Borgonovo δ) → `GSA_results.csv` + `GSA_heatmap.png` |

> **Always run Section 4 first.** If the 5 verification scores match the reference, the
> full run will too. If they don't, stop and check the kernel/numpy version and that all
> databases are scanned.

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
├── MC_workflow.ipynb        # the complete, verified workflow
├── environment.yml          # conda env (numpy < 2, brightway2, SALib, …)
├── requirements.txt         # pip alternative
├── examples/
│   └── PROB_X_t0_short.csv   # 10‑scenario sample matrix for a fast verification run
├── docs/
│   └── WORKFLOW.md          # how to adapt this to your own Brightway2 project
└── LICENSE
```

The Brightway2 project backup, ecoinvent databases and course materials are **not**
included (licensing + size). Point the notebook at your own project.

---

## Requirements

- Python ≥ 3.9
- **numpy < 2** (hard requirement for `bw2data` 3.6.x)
- `brightway2`, `bw2data`, `bw2calc`, `bw2io`
- `SALib`, `scipy`, `pandas`, `matplotlib`, `seaborn`

See [`environment.yml`](environment.yml).

---

## Citing / context

Developed as part of an LCA uncertainty course (III‑V/Si tandem PV case study). If this
workflow is useful in academic work, a mention is appreciated.

## License

[MIT](LICENSE)
