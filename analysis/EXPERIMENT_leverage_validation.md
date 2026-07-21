# Validation experiment: does the leverage decomposition predict PSIS-LOO reliability?

**Purpose.** Turn the base-leverage theory (`THEORY_base_leverage.md`) into a
falsifiable test. The σ_u-shift check already done validates the *influence*
algebra against refit; this validates the **operational claim** — that the two
leverage terms predict PSIS-LOO k̂ and that RB-LOO cures exactly the folds the
fiber term flags. Pre-register the thresholds below before running.

## The two hypotheses (falsifiable)

- **H1 (fiber / pooling).** A fold's full-θ PSIS-LOO k̂ rises with its total
  influence L_i, and RB-LOO (marginalising the fiber analytically) *reduces* k̂ on
  exactly the folds with a large fiber term (large local 1−π_j × residual²).
- **H2 (base / horizontal).** After RB, the residual base-only importance weights'
  k̂ rises with the base leverage L_i^base = g̃_i′M⁻¹g̃_i; the folds where RB-LOO
  still fails (base k̂ > 0.7) are the high-base-leverage ones, i.e. RB does *not*
  rescue outlier-in-small-group folds — and the base leverage says so a priori.

## Pre-registered success / falsification thresholds

| Test | Success | Falsified if |
|---|---|---|
| H1a | Spearman(L_i total, k̂_full) ≥ 0.7 | < 0.4 (leverage doesn't predict k̂) |
| H1b | RB-LOO drops k̂ below 0.7 on ≥80% of folds with high fiber term & small base term | RB fails to reduce k̂ there |
| H2a | Spearman(L_base, k̂_base) ≥ 0.6 | < 0.3 (base leverage doesn't predict residual failure) |
| H2b | AUC ≥ 0.8 for "k̂_base > 0.7" predicted by top-decile L_base | AUC ≤ 0.65 |
| Acc | RB-LOO elpd matches refit-LOO within 2·MCSE on low-L_base folds | systematic bias beyond MCSE |

A clean falsification is itself a result: it would mean the linearised influence
is too crude a k̂ proxy (push to a second-order or exact-refit leverage), or that
k̂ is driven by something the metric split misses.

## Design

For each model: fit the full-data posterior once; then per observation compute
(a) analytic leverages from the Fisher blocks / fibr, (b) full-θ PSIS-LOO k̂ and
elpd (`loo`), (c) RB-LOO — analytic fiber-marginalised pointwise predictive with
base-only importance weights, giving k̂_base and elpd_rb, (d) brute-force refit
LOO on a stratified subset (gold standard: low, mid, high leverage folds).

Quantities to correlate/plot:
- k̂_full vs L_i (total)                          — H1a
- (k̂_full − k̂_base) vs fiber term               — H1b (does RB help where fiber is big?)
- k̂_base vs L_base                               — H2a
- ROC of {k̂_base > 0.7} from L_base              — H2b
- elpd_rb vs elpd_refit                          — accuracy
- structural (residual-free) leverage vs realized k̂ — bonus: how much predictive
  power is a-priori (design only) vs needs residuals

## Model ladder (cheapest first)

1. **Gaussian random-intercept LMM** — analytic, no sampler needed (conjugate;
   direct posterior draws). Full control: everything (leverages, RB-LOO, refit) is
   exact, so this isolates whether k̂ itself tracks the leverage. Include
   singleton and outlier groups (reuse the `verify_base_leverage.R` design).
   *This one can run in the current cloud R immediately.*
2. **GLMM (Poisson & Bernoulli random intercept)** — via brms/rstan. Tests the
   Pólya-Gamma conditional-Gaussian extension (the info term becomes the PG
   working weight). Confirms the leverage formulas and RB-LOO under a non-Gaussian
   outcome where the conditional is still exactly Gaussian.
3. **bjlm smooth change-point** — the boundary case. Here the fiber conditional is
   only linearised (Phase 2), so expect the leverage/k̂ link to degrade at extreme
   sharpness. Deliberately map *where* it breaks; this is the honest limit and
   ties back to the change-point RB finding. Needs the bjlm sampler (later).

## Instrumentation notes

- k̂ from `loo::psis`; use its own MCSE for the accuracy test.
- For H1b/H2, define "high fiber term" and "high base leverage" by the analytic
  quantities *before* looking at k̂ (no post-hoc thresholds).
- Report the *structural* vs *influence* split explicitly: how much of k̂ is
  predictable from design alone (structural h_i, computable pre-data) vs needs the
  realised residuals (influence L_i). This is the practically important line — a
  purely structural predictor enables design-time triage.
- Stratify the refit gold-standard sample across the leverage range; do not sample
  uniformly (the high-leverage tail is rare and is the whole point).

## What a pass would establish

That fibr's Fisher blocks yield, in closed form, a per-observation reliability map
for PSIS-LOO on hierarchical models — flag folds PSIS will fail (fiber term), cure
the curable ones (RB-LOO), and identify the genuinely refit-only ones (base term)
— validated against k̂ and refit. That is the empirical spine of both the LOO
paper and the fibr-paper repositioning.

## First action

Build Model 1 (Gaussian LMM) in the cloud R now: it needs no Stan, reuses the
verified leverage code, and settles H1a/H2a/accuracy under full analytic control
before spending MCMC on the GLMM. If Model 1 falsifies H1a or H2a, stop and
rethink before touching brms or bjlm.
