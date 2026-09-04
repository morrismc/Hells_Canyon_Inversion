# Multi-basin inversion framework

Runs the Hells Canyon inversion over several basins with **globally shared
settings**, so that cross-basin agreement in capture timing is a result
rather than an artifact of per-basin tuning.

## Files

| File | Role |
|---|---|
| `hc_basin_defaults.m` | Everything shared: burn-in, dt, bounds, proposal tuning, cave priors |
| `hc_basin_list.m` | One entry per basin — only what differs |
| `run_hc_inversion.m` | The engine: the MCMC as a function, safe inside `parfor` |
| `run_all_basins.m` | Driver: parallel loop, per-basin save, error isolation |
| `hc_multibasin_summary.m` | QA table + cross-basin comparison figure |

## Quick start

```matlab
addpath(genpath('C:\path\to\topotoolbox'))
summary = run_all_basins();                       % all configured basins
summary = run_all_basins('n_postburn', 5e4);      % shorter test
summary = run_all_basins('parallel', false);      % serial, for debugging
```

Any field of `hc_basin_defaults` can be overridden globally as a
name/value pair.

## Design decisions worth knowing

**Basins are inverted independently.** Each gets its own `t_capture`. Five
tributaries independently converging on ~2.2 Ma is far stronger evidence
than fitting one shared parameter by construction. A hierarchical model
with shared timing is the natural follow-up *if* the independent runs
justify it.

**Cave priors apply to every basin, including cave-less tributaries.**
These streams are graded to the Snake, so their base-level history *is*
the Snake's incision history, which the caves record. Without them,
`U_pre` and `U_post` are unidentifiable in a cave-less basin: the profile
alone cannot separate `U` from `K`, because relict steepness is
`(U_pre/K)^(1/n)` and, with `ksn_ref` sampled, `U_pre` and `K` trade off
exactly.

**Likelihood weighting is held identical with or without caves.**
`hc_loglikelihood` applies the `N_eff` down-weighting only when cave data
is present, so a cave-less basin would otherwise carry ~50× more stream
weight and would not be comparable. The engine folds the same factor into
`stream_sigma` instead — which is how Gallen balances datasets anyway.

**Parallelism is free.** Basins are independent, so `parfor` needs no
coordination. Without Parallel Computing Toolbox, MATLAB runs `parfor`
serially and the same code still works. Plotting happens after the loop,
never inside it.

## Before trusting a multi-basin run

1. **Regression-test the engine.** Run basin 1 (the already-validated trunk
   basin) with a fixed `rng_seed` through `run_hc_inversion` and confirm it
   reproduces what `main_hc_bayes_inversion.m` gives. The refactor must not
   change the science.
2. **Read the QA table, not the figures.** `hc_multibasin_summary` flags
   `RAILED` / `ACCEPT` / `NONFINITE` / `FWDFAIL`. A basin flagged `RAILED`
   is reporting its prior, not its data.
3. **Expect to drop basins.** Small, steep tributaries may be debris-flow
   dominated, where stream power does not apply, and some will have either
   fully adjusted (knickpoint already exited) or not responded at all —
   either way carrying no timing information. Decide the exclusion
   criterion *before* seeing the answers.

## ⚠️ Sync obligation

`run_hc_inversion.m` currently **duplicates** the MCMC loop from
`main_hc_bayes_inversion.m`. This was deliberate: it means this branch adds
files only and touches nothing on the single-basin branch, so the merge is
clean while that recipe is still being iterated on.

The cost is that **any change to the model or sampler in
`main_hc_bayes_inversion.m` must be mirrored here**, or multi-basin runs
will silently use stale science.

Once the single-basin recipe is frozen, do the one-time reconciliation:
make `main_hc_bayes_inversion.m` a thin wrapper that builds a `cfg` and
calls `run_hc_inversion`, leaving one implementation. That is a small,
deliberate commit best made when nothing else is in flight.
