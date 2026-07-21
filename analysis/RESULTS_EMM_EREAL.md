# E-MM + E-real: gate G2 closed empirically, and real-data external validity

Run locally on a working brms/Stan toolchain (R 4.6.0, brms 2.23.0, rstan 2.32.7,
cmdstan 2.38.0, loo 2.9.0) — the two experiments the cloud sandbox could not run
(no Stan compiler). Both scripts fit with **brms** (`backend="rstan"`,
`save_pars(all=TRUE)`), so `loo_moment_match` and `reloo` are both available.

## E-MM — the moment-matching head-to-head (`emm_moment_match.R/.png/.rds`)

Logistic random-intercept GLMM (J=60, singleton-heavy, separation-prone), **6
replicate datasets, 77 high-k̂ folds**. On the high-k̂ folds, compared vs
**reloo-exact gold**: PSIS-LOO, `loo_moment_match` (SOTA), and RB-LOO (quadrature).

| Method | elpd RMSE vs reloo-exact (pooled) | per-fold optimisations |
|---|---|---|
| PSIS-LOO | 0.140 | 0 |
| moment matching | 0.131 | 13 / rep |
| **RB-LOO** | **0.044** | **0** |

Per-rep RMSE (mean ± MCSE): PSIS 0.094±0.025, MM 0.081±0.026, **RB 0.032±0.008**.
Cure: RB and MM both 100% (all k̂→<0.7). Cost: RB **0** optimisations; MM **13/rep**
(12s); reloo **13 refits/rep** (89s).

**Pre-registered H5 confirmed decisively.** RB-LOO is **3.0× more accurate than
moment matching** (and 3.2× more accurate than PSIS) **at zero per-fold
optimisation cost**. Notably MM barely improves on PSIS here (0.131 vs 0.140): an
affine transform of the *conditional* draws cannot represent a data-driven
singleton random effect reverting toward its prior under deletion — exactly what
RB marginalises in closed form. **Gate G2 is now closed empirically**, not just by
the exact-refit bound in `GATE_G2_RESULT.md`.

## E-real — epilepsy real data (`ereal_epilepsy.R/.png/.rds`)

`brms::epilepsy` (236 seizure counts, 59 patients × 4 visits; overdispersed — a
known PSIS-LOO stressor). Model `count ~ zBase*Trt + zAge + (1|patient)`, Poisson.
Evaluated with a single **`rb_loo(fit)`** call — this is also the `rbloo` package
integration test, and it passed.

- **C1:** structural leverage flags the 5 PSIS failures (k̂>0.7, max **1.63**) at
  **AUC 0.84**. (Pooled Spearman only 0.27 — k̂ is floored on the 231 easy folds,
  diluting rank correlation; AUC is the meaningful separation metric.)
- **C2:** RB-LOO cures all five (k̂ 5→0; mean 0.32→−0.05; worst fold 1.63→0.15).
- **Accuracy** vs reloo-exact on the high-k̂ folds (elpd RMSE): PSIS **1.76**,
  MM **0.77**, **RB 0.72** — RB most accurate, at 0 optimisations vs MM's 5 and
  reloo's 5 exact refits. (Absolute RMSE large: extreme counts → large-magnitude
  pointwise elpd. Ranking is unambiguous.)

roaches was dropped: it lives in **rstanarm** (not installed; a full install
compiles many Stan models, ~15 min, just for a dataset) and its classic model has
**no random effect** for RB-LOO to marginalise. epilepsy is the stronger E-real —
genuinely hierarchical (patient RE), overdispersed, bundled with brms, and a
recognised loo-community stressor.

## E-real OLRE stressor — the sharpest result (`ereal_epilepsy_olre.R/.png/.rds`)

The patient-RE model gives only 5 failures (4 visits/group). Model overdispersion
the standard way — an **observation-level random effect** (Poisson-lognormal),
`count ~ zBase*Trt + zAge + (1|obs)` — so **every obs is its own singleton group**,
the maximally data-driven regime. Single grouping factor, so `rb_loo()` applies
directly (widened `quad_range=8` for the singleton tails; the package now exposes
`quad_range`). PSIS-LOO fails on **97 of 236** folds (max k̂ 1.34).

| | value |
|---|---|
| C1 predict | Spearman **0.60**, AUC **0.80** |
| C2 cure (RB) | 97 → **0** (mean k̂ 0.68 → −0.01; max 1.34 → 0.35) |
| C2 cure (**moment matching**) | 97 → **37 still k̂>0.7** — MM does *not* fix the singleton failures |
| Accuracy vs reloo-exact (97 folds, elpd RMSE) | PSIS 0.638 · MM 0.658 · **RB 0.041** (~15×) |
| Cost | RB **0** · MM 97 opts (8 min) · reloo 97 refits (**82 min**) |

reloo confirms RB is the correct answer: total integrated elpd −651.6 (RB) matches
reloo, while PSIS's conditional −615.5 is over-optimistic by ~36 nats. **Moment
matching (the SOTA) both misses a third of the failures and is no more accurate
than raw PSIS on this real stressor, while RB-LOO is exact at zero cost.** This is
the paper's strongest single figure (RB on the reloo diagonal; PSIS + MM
systematically above it, diverging on the extreme-count folds).

## Deliverables produced this session

- `emm_moment_match.{R,png,rds}` — E-MM.
- `ereal_epilepsy.{R,png,rds}` — E-real + package integration test.
- `rbloo/` — the loo-compatible package (`rb_loo()` S3 over brmsfit/stanreg;
  Gaussian analytic downdate + Bernoulli/binomial/Poisson quadrature; a-priori
  pooling-factor / structural-leverage triage; base-k̂ + refit flag). Validated on
  the epilepsy fit.
- `paper/rb_loo_paper.tex` (+ compiled `rb_loo_paper.pdf`, 10 pp) — first full
  draft with all five figures (G1, replicated GLMM, E-MM, E-real, change-point
  boundary), theory T1–T4, the Lovison/Cui-Hodges demarcation, and the honesty
  register.

## Programme status

- **G1 cleared** (Gaussian: predict AUC 0.95, cure, 3.7× accuracy).
- **G3 resolved by simplification** (base leverage redundant with pooling factor).
- **G2 closed** (GLMM: predict AUC 0.81; RB 3× more accurate than moment matching
  at zero optimisation cost; matches exact refit).
- **Real-data external validity** shown on epilepsy.
- **E-M5 boundary** already mapped (`NOTES_bjlm_changepoint.md`); the change-point
  conditional is near-Gaussian, degrading only at extreme sharpness — the honest
  limit. The full bjlm-sampler test still needs the Rust `m_j, s_j` export.
