# Phase 2 — RB-PIT in the bjlm setting: the non-linear change-point

Moving the RB-PIT into bjlm (smooth logistic change-point outcome + latent GP).
The H0 pipeline map settles the structure:
- Intercept RE `u_b0` is EXACTLY Gaussian given hypers (Gaussian / Polya-Gamma) ->
  replace-RB is exact there; `m_j, s_j` are already computed in Rust's
  `sample_random_effects_weighted`, and `z_j=(u_j-m_j)/s_j ~ N(0,1)` is already
  the planned ASIS non-centring handle. The sampler machinery and the RB-PIT
  accelerator are the SAME object.
- The change-point enters `mu` NON-LINEARLY (HMC, not conjugate). This is the one
  coordinate where "RB breaks" was the open worry, so it is what I tested.

## What was tested

Controlled isolation (all other params known = best case for RB): draw
`omega~ ~ prior`, simulate `y_i = mu(tau_i;omega~)+N(0,sigma^2)` with bjlm's exact
mean, compute the EXACT conditional `p(omega|y)` on a dense grid, and compare
three PITs to uniform:
- exact           -- positive control (uniform by construction; harness check)
- moment-matched  -- best possible Gaussian (exact cond. mean+sd); miscalibration
                     = pure NON-GAUSSIANITY of the conditional
- Fisher / pf     -- bjlm's approximation: mode + Fisher precision Glik =
                     sum(dmu^2)/sigma^2, using pooling_factor.R's exact gradient
                     dmu/domega = -(b1 + delta*s*(1+d*rho*(1-s))).
Swept identification strength (n_obs) and sharpness (rho); S=800, 0.95 KS crit 0.048.

## Findings

1. **The change-point conditional is close to Gaussian.** The moment-matched PIT
   is calibrated across every regime (KS 0.018-0.043, all in band), including
   sparse (n_obs=3) and sharp (rho=20). Non-Gaussianity is NOT the dominant risk.

2. **bjlm's Fisher/pf approximation is adequate, with mild degradation at extreme
   sharpness.** In band for rho<=8 and all n_obs; edges OUT only at rho=20
   (KS 0.060 vs 0.048) where the Gauss-Newton curvature (dropping mu'') and the
   mode!=mean gap start to bite in the data-dominated sharp kink. The effect is
   modest and near the S=800 noise floor.

3. **Sparsity and sharpness alone do not break it.** The prior-dominated regimes
   are trivially Gaussian (conditional ~ prior); the interesting failures are
   data-dominated + very sharp, and even there the error is small.

## Reading

RB on the change-point is VIABLE in this isolated setting. Practical guidance:
prefer a **moment-matched** conditional (exact conditional mean + sd) over
mode+Fisher for the change-point -- it is uniformly the more robust of the two and
costs nothing extra once you have the conditional. Reserve mode+Fisher (the raw
pooling_factor precision) for a quick screen.

## The caveat that matters

This test fixes ALL other parameters at truth. It therefore excludes the real
bjlm risks, which are interactions:
- uncertainty in b1, delta, rho, sigma feeding the change-point conditional;
- mode!=mean amplified once the base (hypers, sigma_u, GP length-scale funnel)
  is itself uncertain and moving m_j, s_j;
- the latent GP field, which the current SBC harness does not rank at all
  (sampled by elliptical slice in Rust, not exposed in R draws).
The decisive test needs the actual sampler with the base varying -- i.e. the
Rust per-draw m_j, s_j surfaced to R.

## Integration seam (concrete)

Current rank: `.sbc_score_functional` (R/sbc.R:294) does `mean(vt < truth)` over
ESS-thinned draws. RB-PIT drop-in:
    u_rb = mean_over_draws( pnorm(truth, m_draw, s_draw) )
where (m_draw, s_draw) are the per-draw conditional moments.
- Intercept RE: m_j, s_j exact, already in `sample_random_effects_weighted`.
- Change-point: use the moment-matched Gaussian conditional (mean+sd), not
  mode+Fisher.
Both require surfacing the moments from Rust -- the SAME `m_j, s_j` the ASIS
non-centring (z_j) already needs. One export unlocks both the sampler fix and the
calibration accelerator.

## Recommended next step

Export `m_j, s_j` (or `z_j`) from `sample_random_effects_weighted` to the R draws,
then: (a) replace-RB the intercept RE in `.sbc_score_functional` and confirm the
E2-style L-slashing win on the real sampler; (b) add the `z_j` within-fit residual
check to catch fiber-update bugs replace-RB is blind to (Phase 1b); (c) extend RB
to the GP field via its conditional-Gaussian-given-hypers posterior, which the
harness currently does not test at all.

Run: `Rscript bjlm_changepoint_rb.R` (base R only; ~1.7 min).
