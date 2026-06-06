# Adapting the workflow to your own Brightway2 case

This guide explains how the workflow is structured and exactly what to change to run it
on a **different** parameterised Brightway2 project. No part of it is specific to the PV
case study except the concrete names and numbers, which are all isolated in a few cells.

> Context: this pipeline unifies three previously separate steps — presampling (originally
> **R**), LCA calculation (originally **Activity Browser**), and GSA (originally **MATLAB**) —
> into one Python notebook. The three steps map to Layers 2–3 (presampling), 5–7 (calculation),
> and 9 (GSA) below.

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

In the structured notebook (`MC_workflow.ipynb`) all case-specific values live in
**Layer 1** and **Layer 2**. The simple notebook has the same content inline in Sections 0–1.

### Layer 1 — environment binding

```python
PROJECT           = "FINAL"                                   # ← your project name
PARAMETERISED_DBS = ["SiTaSol_F2v1a", "IEA_PVPS_2020"]        # ← EVERY DB with formula exchanges
FU_DB, FU_CODE    = "SiTaSol_F2v1a", "ba8e2884...copy1"       # ← your functional-unit activity
METHOD            = ('ILCD 2.0 2018 midpoint',
                     'climate change', 'climate change total') # ← your LCIA method
```

How to find these:

```python
bd.databases                                             # list database names
list(bw.Database("SiTaSol_F2v1a"))[:5]                   # browse activities; copy the code you want
[m for m in bw.methods if 'climate' in str(m).lower()]   # find the method tuple
```

> `PARAMETERISED_DBS` is the cell people get wrong. If a parameter lives in an exchange in
> a database you don't list, that parameter silently does nothing and its GSA sensitivity
> becomes a meaningless zero. **Layer 5's validation gate now catches this for you** (see §2.1),
> but list every relevant database anyway. When in doubt, use `list(bd.databases)` (skipping
> `biosphere3` is fine — it has no formula exchanges — but scanning it is harmless).

### Layer 2 — the parameter registry (the main thing you edit)

One declarative table is the single source of truth. Add/remove/reorder a parameter by
editing **one row**; the sampling, CSV columns, GSA labels and `var_level` are all derived.

```python
PARAMETERS = [
  dict(name="MyParam", dist="normal",  args=dict(loc=30, scale=5),               group="project",  desc="..."),
  dict(name="Frac",    dist="pert",    args=dict(min=.2, mode=.3, max=.7),       group="project",  desc="..."),
  dict(name="Energy",  dist="lognormal", args=dict(s=np.log(1.22), scale=110),   group="activity", desc="..."),
  dict(name="pi_x",    dist="pert",    args=dict(min=.5, mode=.7, max=.8),       group=None,       desc="helper"),
  dict(name="Switch",  dist="bernoulli", args=dict(p="pi_x"),                    group="project",  desc="dependent draw"),
  # ...
]
```

Rules:

* **`name` must match the variable used in the Brightway2 exchange `formula`.**
* `dist` is a key of the `SAMPLERS` dict; `args` are that sampler's keyword arguments.
* A **string** `args` value (`p="pi_x"`) means *use the already-sampled array of that
  parameter* — this is how dependent draws (a Bernoulli whose probability is itself
  uncertain) are expressed without leaving the table.
* `group="project"`/`"activity"` only labels the exported CSV. Set `group=None` for
  intermediate helpers (the `pi_*`) so they're excluded from the model output and the GSA.
* **List order = draw order.** Put a dependent entry after the helper it references so the
  random-number stream is reproducible.

Need another distribution? Add one line to `SAMPLERS` (a `lambda size, **args: ...`), no
other change required.

### 2.1 Layer 5 — validation gate (you don't edit it, but rely on it)

After indexing the formula exchanges, the notebook parses every formula and cross-checks it
against your registry:

```python
unresolved = used - declared    # formula needs a parameter you didn't supply  -> raises
unused     = declared - used    # parameter affects nothing (wrong/missing DB)  -> warns
```

`unresolved` raises immediately; `unused` prints a warning. Either message points straight
at a mis-typed name or a forgotten database — the bugs that are otherwise invisible until
the numbers come out subtly wrong.

### Layer 6 — verification reference

Put a few known-good scores in `ref` (from a trusted run on a freshly restored project, or a
previously validated Python run) so every future run is checked automatically. If you have
no reference yet, run 5 iterations once and record the numbers.

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
