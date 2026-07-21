# rbloo

**Rao-Blackwellised leave-one-out for hierarchical models** — a closed-form,
a-priori triage-and-cure for hierarchical cross-validation.

For hierarchical models, PSIS-LOO fails on exactly the folds where a
random-effect coordinate is data-driven and its group is small. The
Gelman–Pardoe **pooling factor** (the [fibr](https://github.com/) package)
already measures that, in closed form, *before* any importance weights are
formed. `rbloo` uses it to triage folds, then **Rao-Blackwellises** the LOO —
marginalising the random-effect block analytically (Gaussian) or by cheap 1-D
quadrature (Bernoulli / binomial / Poisson) — so the fiber-driven failures
disappear at **zero refit cost**, matching exact refits.

## Install

```r
# install.packages("remotes")
remotes::install_github("ABindoff/rbloo")
# or from a local clone of this repo:
# devtools::install(".")
```

The paper and the scripts that reproduce every figure live in [`analysis/`](analysis/)
(see `analysis/REVIEW_synthesis.md` and the manuscript `analysis/paper/`).

## Use

```r
library(brms); library(rbloo)
fit <- brm(count ~ zBase * Trt + zAge + (1 | patient),
           data = epilepsy, family = poisson(),
           save_pars = save_pars(all = TRUE))

rb <- rb_loo(fit)
rb
#> rb_loo  (poisson GLMM; N=236 obs, G=59 groups)
#>   elpd_rb = ...    elpd_full(PSIS) = ...
#>   PSIS-LOO : #(k>0.7) = ...
#>   RB-LOO   : #(k>0.7) = 0   <- fiber failures removed
#>   refit_flag: 0 fold(s) still k_base>0.70 -> send to reloo
```

The returned object carries, per observation:

| field | meaning |
|---|---|
| `pooling_factor` | Gelman–Pardoe pooling factor π_j of the obs's group (a-priori) |
| `structural_leverage` | structural fiber leverage h_i; sums to 1−π_j over a group (a-priori) |
| `pointwise$elpd_rb` | Rao-Blackwellised pointwise elpd (the cure) |
| `diagnostics$pareto_k` | base Pareto-k̂ after RB (should be < 0.7) |
| `diagnostics$pareto_k_full` | full conditional PSIS-LOO k̂ (the incumbent) |
| `refit_flag` | TRUE where k_base still > `base_cut` → send to `reloo` |

## What it targets

RB-LOO targets the **integrated** (marginal-over-RE) predictive; compare only to
refits / `reloo` using the same target (Merkle, Furr & Rabe-Hesketh 2019). RB is
**exact** for Gaussian / Pólya-Gamma-Gaussian RE blocks and **approximate**
(1-D quadrature) for the other GLMM families.

## Supported models

RB-LOO applies to a fit with **one grouping factor**, on the default link, in the
families `gaussian`, `bernoulli`, `binomial`, `poisson`. The random-effect
structure it handles depends on the family:

- **`gaussian`: any random-effect design** — random intercept, random slopes,
  correlated or not (`(1|g)`, `(1+x|g)`, `(0+x|g)`, `(x1+x2|g)`). The conditional
  is exactly Gaussian, so RB-LOO is the closed-form `p × p` matrix downdate and
  is **exact** (validated fold-by-fold against brute-force refits).
- **`bernoulli` / `binomial` / `poisson`: a single random intercept**, via 1-D
  quadrature (numerically exact to grid resolution).

Anything outside this scope is **not** Rao-Blackwellised — `rb_loo()` never
returns a wrong RB number. Instead it emits a **loud warning** naming the reason
and **falls back to plain PSIS-LOO** (`brms::loo()`), returning an `rb_loo` object
unmistakably marked as PSIS-only: the RB-specific fields (`pooling_factor`,
`structural_leverage`, `pointwise$elpd_rb`, `diagnostics$pareto_k`) are all `NA`,
only the PSIS/full fields are populated, and `print()` leads with a
`PSIS-LOO fallback (RB-LOO not applied)` banner plus the Pareto-k count and advice.
The out-of-scope cases that trigger this are:

- multiple / crossed / nested grouping factors (`(1|g1)+(1|g2)`, `(1|g1/g2)`);
- **random slopes on a non-Gaussian family** (multivariate GLMM quadrature is not
  yet implemented — the Gaussian case above is exact);
- unsupported families (e.g. `negbinomial`, `Gamma`, `student`, ordinal,
  zero-inflated / hurdle);
- distributional / heteroscedastic models (where `sigma` or another parameter is
  itself regressed);
- non-default links, and fits that use observation weights.

PSIS-LOO is a valid computation for these models, but it can be unreliable on
influential observations: inspect the Pareto-k̂ diagnostics, and for folds with
k̂ > 0.7 use `reloo()` or `loo_moment_match()`. Malformed arguments and non-fit
objects remain hard errors — those are usage errors, not model-scope limitations.

`brmsfit` is fully supported and tested. The `stanreg` (rstanarm) path shares the
same estimator and the same warn-and-fall-back behaviour, but is not yet validated
end-to-end; it emits a warning to that effect.

## References

- Gelman & Pardoe (2006), *Technometrics* — the pooling factor.
- Bindoff (2026), *fibr* — structural leverage, the base–fiber bundle.
- Vehtari, Gelman & Gabry (2017), *Stat. Comput.* — PSIS-LOO and k̂.
- Paananen, Piironen, Bürkner & Vehtari (2021) — moment matching (the comparator).
- Merkle, Furr & Rabe-Hesketh (2019), *Psychometrika* — conditional vs marginal LOO.
