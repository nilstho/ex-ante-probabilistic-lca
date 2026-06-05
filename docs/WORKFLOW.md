# Adapting the workflow to your own Brightway2 case

This guide explains how the workflow is structured and exactly what to change to run it
on a **different** parameterised Brightway2 project. No part of it is specific to the PV
case study except the concrete names and numbers, which are all isolated in a few cells.

---

## 0. Mental model

The workflow rests on one idea:

> A parameterised Brightway2 exchange stores a **`formula`** string (e.g. `RT_movpe*P_movpe_tool`).
> To run a Monte Carlo scenario you simply **evaluate that formula** with a set of scalar
> parameter values and write the result into the exchange's `amount`, then run the LCA.

Doing this for **every** parameterised exchange, **every** iteration, makes the result a
pure function of the parameter vector `X[i]`. That is what guarantees reproducibility and
makes the result independent of any prior database state.

```
for each scenario i:
    params = { name: value for name, value in zip(param_names, X[i]) }
    for each parameterised exchange:
        exchange['amount'] = eval(exchange['formula'], params)
        exchange.save()
    score_i = LCA(functional_unit, method).score
```

Everything else (sampling, plotting, SALib) is standard.

---

## 1. Prerequisites in your Brightway2 project

Your model must already be **parameterised inside Brightway2**, i.e. the relevant
exchanges carry a `formula` field and the parameter names used in those formulas are
registered as project/activity/database parameters. You can confirm this:

```python
from bw2data.parameters import ParameterizedExchange
for pe in ParameterizedExchange.select():
    print(pe.formula)
```

If that prints your formulas, you are ready. (This is exactly what the Activity Browser
sets up when you define parameters there — you can keep using AB to *build* the model and
this workflow to *run* the analysis.)

---

## 2. The cells you must edit

All case‑specific values live in **Section 0** and **Section 1** of the notebook.

### Section 0 — project, databases, functional unit, method

```python
bd.projects.set_current("FINAL")                       # ← your project name

pv   = bw.Database("SiTaSol_F2v1a")                    # ← your foreground database
iea  = bw.Database("IEA_PVPS_2020")                    # ← any other DB with formula exchanges
elec = pv.get("ba8e2884bda6251e448284b8d4e9cc15_copy1") # ← your functional-unit activity
ipcc = bw.Method(('ILCD 2.0 2018 midpoint',
                  'climate change', 'climate change total'))  # ← your LCIA method
```

How to find these:

```python
bd.databases                              # list database names
list(pv)[:5]                              # browse activities; copy the code you want
[m for m in bw.methods if 'climate' in str(m).lower()]   # find the method tuple
```

### Section 1 — the parameter samples

Replace the sampling block with **your** parameters and distributions. Each line draws
`replications` samples; the variable name **must match the name used in the formulas**.

```python
np.random.seed(42)            # keep a fixed seed for reproducibility
replications = 10000

MyParam = norm.rvs(mean, sd, size=replications)        # normal
Other   = uniform.rvs(low, width, size=replications)   # uniform on [low, low+width]
Frac    = rpert(replications, min=.2, mode=.3, max=.7) # PERT (helper defined in the cell)
Choice  = bernoulli.rvs(0.5, size=replications)        # 0/1 switch
# ... lognormal, triangular, beta are also imported
```

Then list the parameters that enter the model (exclude purely intermediate ones such as
the `pi_*` success‑probabilities used only to draw Bernoulli switches):

```python
param_names = ['MyParam', 'Other', 'Frac', 'Choice', ...]
```

`var_level` (project vs activity) is only used to label the exported CSV in the
AB‑compatible format; it does not affect the calculation. Set every entry to `'project'`
if you don't care about that distinction.

### Section 3 — which databases to scan

```python
for db in (pv, iea):          # ← include EVERY database that contains formula exchanges
    for act in db:
        for exc in act.exchanges():
            if 'formula' in exc:
                formula_exchanges.append(exc)
```

**This is the cell people get wrong.** If a parameter lives in an exchange in a database
you don't scan, that parameter silently does nothing and its GSA sensitivity will be a
meaningless zero. When in doubt, scan all databases:

```python
for db in (bw.Database(n) for n in bd.databases):
    ...
```

(Skipping `biosphere3` is fine — it has no formula exchanges — but scanning it is harmless.)

### Section 4 — verification reference

Put a few known‑good scores here (e.g. from a trusted AB run *on a freshly restored
project*, or from a previous validated Python run) so every future run is checked
automatically. If you have no reference yet, just run 5 iterations once and record the
numbers.

---

## 3. Distribution cheat‑sheet (`scipy.stats`)

| Distribution | Call | Notes |
|---|---|---|
| Normal | `norm.rvs(loc, scale, size=N)` | `scale` = standard deviation |
| Uniform | `uniform.rvs(low, width, size=N)` | support is `[low, low+width]` |
| Triangular | `triang.rvs(c, loc, scale, size=N)` | `c` = relative mode in `[0,1]` |
| Lognormal | `lognorm.rvs(s, scale=exp(μ), size=N)` | `s` = σ of underlying normal; median = `scale` |
| Beta | `beta.rvs(a, b, size=N)` | on `[0,1]` |
| Bernoulli | `bernoulli.rvs(p, size=N)` | 0/1 switch |
| PERT | `rpert(N, min, mode, max, shape=4)` | helper defined in the notebook (beta‑based) |

The PERT helper:

```python
def rpert(size, min, mode, max, shape=4):
    alpha  = 1 + shape * (mode - min) / (max - min)
    beta_p = 1 + shape * (max - mode) / (max - min)
    return beta.rvs(alpha, beta_p, size=size) * (max - min) + min
```

---

## 4. Reading the GSA output

`delta.analyze` returns two measures per parameter:

- **`delta`** — Borgonovo's δ: a moment‑independent importance measure (uses the whole
  output distribution, not just the variance). Higher = more influential. Robust to
  non‑linear / non‑monotonic relationships, which is why it suits LCA models with switches.
- **`S1`** — first‑order Sobol index (variance‑based), provided for comparison.

`*_conf` columns are bootstrap confidence intervals. Rank parameters by `delta`; the
heatmap in Section 7 visualises this.

> GSA needs enough samples to be stable. For ~20 parameters, use **≥ 1000** iterations
> (the case study uses 10 000). Check that the `delta_conf` intervals are small relative
> to `delta`.

---

## 5. Performance notes

- Each iteration does several SQLite writes (`exc.save()`) plus a full `lci()/lcia()`.
  10 000 iterations on ~30 parameterised exchanges typically takes tens of minutes.
- The dominant cost is `lci()` rebuilding the technosphere matrix each time. If you need
  it faster, look into building the matrix once and modifying it in place with
  `bw2calc`'s lower‑level API, or using `presamples` — but the simple approach here is the
  most transparent and is what guarantees an exact match with the formula semantics.
- Do **not** parallelise by sharing one Brightway2 project across processes — each writes
  to the same SQLite file. Use separate project copies (or separate `BRIGHTWAY2_DIR`s) if
  you parallelise.

---

## 6. Reproducibility checklist

- [ ] Kernel uses **numpy < 2** (Section 0 prints a warning otherwise).
- [ ] `param_names` order matches the columns of `PROB_X_t0.csv` and the order you want in
      the GSA.
- [ ] **All** databases containing formula exchanges are scanned in Section 3.
- [ ] Section 4 verification shows `Diff ≈ 0` before launching the full run.
- [ ] A fixed random seed is set in Section 1.

---

## 7. Why not just trust the Activity Browser scenario export?

During development, AB's scenario analysis produced **different results for the same
inputs on different runs**, and the numbers matched no consistent application of the
parameters (the best‑fitting subset of the affected activity parameters differed per
scenario). The Python workflow, by contrast, is deterministic given `X` and reproduces the
documented parameter approach exactly. Keep using AB to **build and parameterise** your
model; use this workflow to **run** the uncertainty/GSA analysis.
