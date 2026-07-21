# Phase 1b — Control-variate PIT: a clarifying negative result

Motivation (from Phase 1 E4): replace-RB is blind to fiber-update bugs because it
discards the sampler's fiber draws. Could a control-variate PIT keep the fiber
draws (stay fiber-sensitive) yet cut variance toward replace-RB by subtracting the
fibr conditional CDF as a correlated control?

    u_cv = (1/M) sum c(base_m)  +  (1/L) sum_{sub} ( g_l - c_l )

with M >> L emulating the target model's cost asymmetry (c cheap on the low-dim
base; fiber = latent GP field, expensive).

## Verdict: no free lunch at small L

**E2 — the CV does NOT reduce variance.** Within-simulation SD: CV = counting to
within 1-4% (variance ratio x1.0), while replace-RB is x500-x15000 smaller. The
hoped-for "best of both" does not materialise.

**E5 — why, decomposed on a fixed dataset.** Per-draw PIT variance splits into
  - fiber-sampling noise  E[c(1-c)] = 0.12 - 0.25   (dominant)   <- CV KEEPS this
  - base-marginal noise   Var(c)    = 0.0004 - 0.009 (tiny)      <- CV removes this
Given the base, the indicator is Bernoulli(c), so g|base carries variance c(1-c).
CV keeps g (that is exactly what makes it fiber-sensitive), so it keeps the
dominant noise. Replace-RB substitutes c for g -> removes it entirely -> which is
precisely why it goes blind. The two goals (low variance, fiber sensitivity) are
in fundamental tension at small L for any indicator-based statistic.

## What CV *does* buy (honestly)

- **Fiber sensitivity retained (E3):** under a fiber bug (alpha 20% too narrow) CV
  detects it (D=0.079) just like counting (D=0.093); replace-RB is blind (D=0.031).
- **Base sensitivity (E4):** CV also flags the sigma^2 bug via its base term.
- **Continuity:** unlike the discrete counting rank (L+1 atoms) CV is continuous, so
  it works with the continuous ECDF band and has no rank-discreteness floor.
- **Caveat:** CV is mildly over-dispersed vs uniform under the correct sampler
  (within-sim noise; E1 D=0.031 just over the 0.027 band), so its null needs the
  extra noise folded in, or more fiber draws, before its power is trustworthy.

## Consequence for the design

There is no shortcut to slashing L for FIBER coordinates: detecting a fiber-
conditional error requires enough fiber draws to resolve the g-distribution, full
stop. The L-slashing win is a BASE-side phenomenon (replace-RB, Phase 1 E2). So:

  base / hyperparameter coords  -> replace-RB, tiny L        (big compute win)
  fiber coords (latent field)   -> genuine fiber draws;
                                   use CV for continuity + base-noise removal,
                                   but budget L for the fiber-noise floor.

fibr's pi_j is the dial that says which coordinates fall on which side.

## Not pursued (dead ends this ruled out)

- CV with M=L collapses algebraically to the counting rank.
- Averaging many fiber draws per base draw to denoise g just reconstructs c
  (Gaussian conditional) -> back to replace-RB -> blind. The circularity is
  intrinsic: in the conditionally-Gaussian block the only fiber calibration signal
  is whether the sampler's conditional matches the analytic one, which needs fiber
  draws to see.

Run: `Rscript cv_pit_prototype.R` (base R only; ~25 s).
