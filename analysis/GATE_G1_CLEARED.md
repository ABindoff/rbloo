# Gate G1: CLEARED — predict (C1) + cure (C2) validated at scale

Ran `rb_loo.R` (P0 reference implementation + P1 replicated experiment).
60 replicate Gaussian random-intercept LMMs (J=12, N=48, singletons + small groups,
σ_u=1.3 to produce a natural spread of k̂), conjugate Gibbs, `loo` 2.6.0, brute-force
refit gold standard on leverage-stratified subsets. 2880 folds total; 118 with
k̂_full > 0.7.

## Pre-registered results (all PASS)

| Hyp | Metric | Result | Threshold | Verdict |
|---|---|---|---|---|
| **H1** | Spearman(pooling/structural leverage, k̂_full) | **0.752** (MCSE 0.014) | ≥ 0.60 | **PASS** |
| **H1** | AUC(k̂_full>0.7 ~ structural leverage), a-priori | **0.946** | ≥ 0.80 | **PASS** |
| **H2** | elpd RMSE vs refit, RB/PSIS ratio | **0.29** (0.038 vs 0.128) | ≤ 0.50 | **PASS** |
| **H2** | high-k̂ folds only, RB vs PSIS RMSE | 0.064 vs 0.236 (**3.7×**) | — | — |
| **H3** | % of k̂_full>0.7 folds cured (k̂_base<0.7) | **100%** (mean k̂ 0.27→−0.02) | ≥ 80% | **PASS** |

## What this establishes

- **C1 (predict, a-priori).** The pooling factor / structural leverage predicts
  per-fold PSIS-LOO k̂ with AUC 0.95 — *before residuals are seen*. A design-time
  reliability map, on solid replicated footing.
- **C2 (cure).** RB-LOO removes 100% of the fiber-driven failures and is 3.4×
  (overall) to 3.7× (high-k̂ folds) more accurate than PSIS-LOO against refit.
- The **minimum-viable paper (predict + cure) is validated.** Per the plan's gates,
  the paper now stands even if C3 (base-leverage triage) or the SOTA comparison
  disappoint.

## Honest caveats

- Gaussian conjugate case only (M1); P2 must reproduce off the Gaussian (GLMM/PG).
- H1's AUC benefits from structural leverage being largely group-size-driven — which
  is exactly the a-priori/design-time claim, not a flaw, but state it plainly.
- "RB closer than PSIS on 70% of folds" is the full subset: on *easy* folds both are
  accurate (PSIS occasionally closer by noise); on *hard* folds RB dominates (H2
  high-k̂ row). Report both.
- RB-LOO targets the *integrated* (marginal-over-RE) predictive; the refit gold
  standard here uses the same target, so the comparison is fair. State the
  conditional-vs-marginal distinction (Merkle) in the paper.
- k̂ is a noisy estimate; results are pooled over 60 replicates with MCSE on the
  headline correlation.

## Software (P0 deliverable)

`rb_loo(y, grp, sigma, draws)` returns, per observation: pooling factor π_j,
structural leverage h_i, base leverage L_base, full PSIS-LOO (k̂_full, elpd_full),
RB-LOO (k̂_base, elpd_rb), and a refit flag. Reference implementation for the M1
model; the packaging template for the loo-compatible tool.

## Next (per EXPERIMENT_PLAN_cv_reliability.md)

- **P2** GLMM (Poisson/Bernoulli) + NB via brms — extend C1/C2 off Gaussian (PG).
- **P3** base-leverage triage (C3) on few-groups + eight-schools — gate **G3**.
- **P4** head-to-head vs moment-matching (B2) + integrated LOO (B4) — gate **G2**.
- **P5** real data (radon, roaches, eight-schools).

Gate G1 status: **CLEARED**. Greenlight the paper.
