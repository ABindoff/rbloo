# Experiment plan — "The pooling factor predicts leave-one-out reliability"

**Working title:** *Partial pooling predicts cross-validation reliability: a
closed-form triage and Rao-Blackwellised cure for hierarchical LOO.*

**One-paragraph thesis.** For hierarchical models, PSIS-LOO fails on exactly the
folds where a random-effect coordinate is data-driven and its group is small — and
the Gelman-Pardoe pooling factor already measures that, in closed form, before any
importance weights are formed. Marginalising the conditionally-Gaussian
random-effect block (Rao-Blackwellised LOO) removes those failures and matches
brute-force refits; a connection-generated *base leverage* predicts the residual
folds where even the marginalised importance sampling is strained and a genuine
refit is unavoidable. The result is a closed-form, a-priori triage-and-cure for
hierarchical cross-validation, with the fibr base–fiber geometry as the
explanation. **No new primitive — pooling factor (fibr), PSIS-LOO (loo),
integrated LOO, and case-deletion influence, composed so each fixes the others'
blind spot.** Supersedes `EXPERIMENT_leverage_validation.md` (that was the pilot).

---

## Contributions

- **C1 (predict, a-priori).** The pooling factor π_j (structural fiber leverage)
  predicts per-fold PSIS-LOO k̂ *before* residuals are seen — a design-time
  reliability map.
- **C2 (cure).** RB-LOO — marginalise the RE block analytically, importance-sample
  only the base — removes the fiber-driven IS failure; matches refit and beats
  PSIS-LOO in accuracy/stability at equal or lower cost.
- **C3 (residual triage).** The base leverage L_base = g̃'M⁻¹g̃ (case-deletion
  influence on the variance components, via the Ehresmann connection) predicts the
  folds where RB-LOO's base IS is still strained → a closed-form "refit-only" flag.
- **C4 (explain / unify).** The base–fiber bundle decomposition (pooling = vertical,
  base leverage = horizontal) explains C1–C3 as one orthogonal identity, tying CV
  reliability to the centring geometry.

## Why this is a paper (demarcation — cite explicitly)

- **Not a leverage/hat-matrix result.** Lovison et al. (2026, Int. Stat. Rev.) own
  the hierarchical augmented hat matrix; their subject hat value *is* the structural
  fiber leverage. Cui-Hodges (2010, Technometrics) partition component complexity
  and cite Gelman-Pardoe. **Neither does case-deletion influence on the variance
  components, nor anything about cross-validation.** Cite both to demarcate; keep our
  leverage claim to the *variance-component influence + CV* part they do not cover.
- **Not just integrated LOO.** Bürkner-Gabry-Vehtari (2021, non-factorized LOO),
  Merkle et al. (2019, conditional vs marginal), and Liu & Rue (2022, leave-group-out)
  marginalise latent structure for LOO. Our addition: (i) the *a-priori pooling-factor
  triage* (which folds to marginalise/refit, before computing weights), (ii) the
  *base-leverage refit flag*, (iii) the closed-form connection explanation, and the
  head-to-head that RB-LOO ≥ moment-matching at lower cost where the block is
  conjugate/PG-Gaussian.
- **Not just PSIS.** Vehtari-Gelman-Gabry (2017) + Paananen et al. (2021,
  moment-matching) + Peruggia (1997) are the IS-cure line we position against.

## Theoretical statements (state/prove; verified numerically already)

- **T1 (orthogonal leverage identity).** g_i'G⁻¹g_i = g_i^F'G_FF⁻¹g_i^F +
  g̃_i'M⁻¹g̃_i, g̃_i = g_i^B + A'g_i^F, A = −G_FF⁻¹G_BF, M = Schur complement.
  *(verified exact, 4e-16.)*
- **T2 (the CV bridge — the key proposition).** L_base = g̃'M⁻¹g̃ = the variance of
  the RB-LOO base importance weights; hence a closed-form predictor of the residual
  base-k̂. *(verified: predicts σ_u LOO shift at cor 0.987; base-k̂ at ρ 0.62–0.71.)*
- **T3 (a-priori predictor).** structural fiber leverage summed over group = 1−π_j
  (*verified exact*); relation to Lovison's hat value stated for demarcation.
- **T4 (RB-LOO estimator).** exact for Gaussian/PG-Gaussian RE blocks; Laplace/Fisher
  approximation for nonlinear blocks (change-point, GP length-scale), with the error
  bounded by the conditional non-Gaussianity (Phase-2 finding).

## The three methods (precise)

1. **Triage:** compute π_j (fibr) and the per-obs structural leverage h_i from the
   Fisher blocks at the posterior mean. Flag likely-high-k̂ folds a priori.
2. **RB-LOO:** L_rb[s,i] = log p(y_i | y_{−i}, φ^(s)) via the analytic rank-1 downdate
   of the RE conditional given base draw φ^(s); feed to loo → elpd_rb, k̂_base.
3. **Refit flag:** compute L_base; folds with L_base above a calibrated cutoff (or
   k̂_base > 0.7) go to genuine refit. Everything else is closed-form.

---

## Pre-registered hypotheses, metrics, falsification

| # | Hypothesis | Primary metric | Success | Falsified if |
|---|---|---|---|---|
| H1 | π_j predicts k̂_full (a-priori) | Spearman; AUC(k̂>0.7) | ρ≥0.6 **and** AUC≥0.80, pooled over models | ρ<0.4 or AUC<0.65 |
| H2 | RB-LOO more accurate than PSIS-LOO vs refit | elpd RMSE ratio; within-MCSE rate | RMSE(RB) ≤ 0.5·RMSE(PSIS); within-MCSE rate higher | RB not more accurate on high-k̂ folds |
| H3 | RB-LOO removes fiber failures | %(k̂_base<0.7 \| k̂_full>0.7) | ≥ 80% | < 50% |
| H4 | base leverage predicts residual strain | Spearman(L_base,k̂_base); AUC | ρ≥0.5 and AUC≥0.75 in strained regimes | ρ<0.3 |
| H5 | RB-LOO ≥ moment-matching at lower cost | elpd RMSE vs refit; #refits/optims | RMSE ≤ MM's **and** 0 per-fold optimisations | RMSE materially worse than MM |
| H6 | boundary maps predictably (honesty) | RB error vs conditional non-Gaussianity | monotone, predictable degradation on M5 | RB fails silently (no warning correlate) |

Pre-register thresholds and predictor definitions *before* unblinding k̂. All
hypotheses judged on replicated runs (≥50 replicate datasets per cell), not single
draws — the pilot showed point estimates are noisy at small N.

## Baselines (head-to-head)

- **B1** full-θ PSIS-LOO (`loo`) — the incumbent.
- **B2** PSIS-LOO + **moment matching** (Paananen et al.; `loo::loo_moment_match`) —
  the SOTA high-k̂ cure; the key comparator for H5.
- **B3** mixture / implicitly-adaptive IS (IWMM) if available.
- **B4** **integrated / non-factorized LOO** (Bürkner-Vehtari) — the closest
  "marginalise the RE" method; RB-LOO positioned as its closed-form,
  triage-equipped, connection-explained special case.
- **B5** **brute-force refit LOO** (gold standard) + K-fold / leave-group-out.
- **B6** naive k̂-threshold triage (refit all k̂>0.7) — vs our a-priori π_j triage
  (predicts before weights) on the cost axis.

## Model ladder + real data

- **M1 Gaussian random-intercept LMM** — conjugate, everything exact. *(pilot done;
  scale + replicate.)* Spine for C1, C2.
- **M2 Gaussian LMM, random slopes / nested / unbalanced** — realistic structure,
  correlated RE (per-marginal π_j caveat).
- **M3 GLMM (Poisson, Bernoulli/Binomial)** — Pólya-Gamma conditional-Gaussian; RB
  exact. Extends C1, C2 off the Gaussian.
- **M4 Negative-binomial GLMM** — a genuine extra base parameter (dispersion r).
- **M5 bjlm smooth change-point + latent GP** — the linearisation boundary (RB
  approximate). Map where it degrades; ties to bjlm and to C4/T4 honesty.
- **Real data:** *radon* (Gelman-Hill multilevel), *roaches* (overdispersed Poisson —
  a known PSIS-LOO stressor, ideal for C2), *eight schools* (small J — the base-strain
  regime, ideal for C3), and a longitudinal dataset for M5. These are the loo
  community's standard proving grounds — beating PSIS-LOO on roaches and triaging
  eight-schools a priori is the persuasive external evidence.

## Experiment → hypothesis → model matrix

| Exp | Tests | Models | Baselines |
|---|---|---|---|
| E1 predict | H1 | M1–M4 + radon | B1 |
| E2 cure-accuracy | H2, H3 | M1–M4 + roaches | B1, B5 |
| E3 residual triage | H4 | M1(few-groups), M3, eight-schools | B5 |
| E4 vs SOTA | H5 | M1, M3, roaches | B2, B3, B4, B6 |
| E5 boundary | H6 | M5 (bjlm) | B5 |
| E6 ablations | robustness | M1, M3 | — |

## Ablations / robustness (E6)

- **Draws S:** RB tolerates small S (the efficiency win) — power/accuracy vs S.
- **Group-size distribution:** singletons drive fiber failure; few groups drive base
  strain — sweep both, confirm the two predictors light up separately.
- **Variance-component prior sensitivity.**
- **Misspecification:** does RB-LOO's marginal target stay sensible under a
  misspecified model? (conditional-vs-marginal LOO subtlety — Merkle.)

## Figures / tables the paper needs

- **F1** schematic: two-level triage (π_j → fiber failure; L_base → base strain) + cure.
- **F2** π_j predicts k̂ a priori, across models (the C1 money plot).
- **F3** k̂ before/after RB + elpd-vs-refit accuracy, RB vs PSIS vs moment-matching.
- **F4** base leverage predicts residual strain (eight-schools / few-groups).
- **F5** cost–accuracy Pareto: #refits vs elpd error — RB-LOO + π_j-triage dominates.
- **F6** the boundary (M5): RB error vs change-point sharpness / conditional
  non-Gaussianity.

## Go/no-go gates + fallback positions

- **G1 (after E1).** If H1 fails at scale → drop the a-priori-triage claim (C1),
  keep the RB-LOO cure (C2) — still a paper ("a closed-form marginalised LOO for
  hierarchical models, with a geometric account").
- **G2 (after E4).** If RB-LOO isn't more accurate than moment-matching → reposition
  C2 as "as-accurate, closed-form, no per-fold optimisation, geometry-explained"
  (cost/interpretability win, not accuracy).
- **G3 (after E3).** The base-leverage triage (C3) is the newest and riskiest. If H4
  fails → drop C3, keep C1+C2. The paper stands on predict+cure alone.

The paper survives any single gate failing; C1+C2 is the minimum viable paper and
is the most-validated part.

## Software deliverable

`rb_loo()` — a loo-compatible R function/package: takes a fit, returns the π_j
triage, RB-LOO elpd + k̂_base, and the L_base refit flag; falls back to refit on
flagged folds. Parallels fibr's packaging. This is the adoption vehicle and a
first-class contribution (the loo ecosystem rewards drop-in tools).

## Sequenced plan (DAG; what can start now)

```
P0 implement rb_loo() + L_base in R, loo-integrated   [reuse verified pilot code]  <- start now
P1 E1+E2 on M1 scaled+replicated  -> C1,C2   [GATE G1]        (needs P0)
P2 E1+E2 on M3,M4 (GLMM/NB, PG)   -> extend C1,C2            (needs P1)
P3 E3 on few-groups + eight schools -> C3    [GATE G3]        (needs P0)
P4 E4 vs B2/B3/B4 head-to-head    -> H5      [GATE G2]        (needs P1)
P5 real data radon/roaches/8-schools -> external validity     (needs P1,P3)
P6 E5 on bjlm (M5)                -> C4 honesty, ties bjlm     (needs bjlm build)
P7 theory writeup T1–T4 + software + manuscript               (needs P1–P5)
```

P0, P1, P3 are pure-R and can start immediately (reuse `verify_base_leverage.R` and
`leverage_psisloo_validation.R`). P2/P4 add brms/rstan (fits + `loo_moment_match`).
P6 waits on the bjlm build. **First action: P0 — wrap the verified pilot code into
`rb_loo()` and re-run E1/E2 on M1 at scale (≥50 replicates) to clear G1.**

## Threats-to-validity register (keep honest)

- k̂ is a noisy estimate at small N/S → replicate; report MCSE on all correlations.
- conditional-vs-marginal LOO: RB-LOO targets the *integrated* predictive; state
  this and compare only to refits using the same target (Merkle).
- correlated random effects: π_j is per-marginal (fibr caveat) — flag on M2.
- external validity: simulated wins must be reproduced on the real datasets or the
  claims stay conditional.
- the base-leverage → base-k̂ link is first-order (Laplace); the highest-influence
  folds are where it's weakest (pilot scatter) — report, don't hide.
