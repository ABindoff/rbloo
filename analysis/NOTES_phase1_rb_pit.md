# Phase 1 prototype — Rao-Blackwellised PIT for SBC

Conjugate normal-normal hierarchy (sigma_y known), exact Gibbs sampler, base-R.
Two SBC statistics from the *same* chain: the standard counting rank, and a
Rao-Blackwellised PIT that integrates each fiber coordinate out analytically
given the base draws (mu, sigma^2), using exactly the moments fibr returns —
prior precision 1/sigma^2, lik_information n_j/sigma_y^2, giving pi_j and the
conditional sd s_j = 1/sqrt(G_FF,j). In production the moment block is a call to
`fibr::prior_fraction()`.

## What the prototype established

**E1 — validity.** Under a correct sampler both statistics are uniform (inside
the simultaneous band): counting D=0.018, RB D=0.012, crit 0.028. The RB PIT is a
legitimate PIT.

**E2 — variance reduction (the headline).** Fixing one dataset and rerunning the
sampler, the within-simulation SD of the PIT estimate falls by 1–2 orders of
magnitude, and the reduction *grows as the coordinate becomes data-dominated*
(pi_j -> 0): effective-L multiplier x24 (n=2) up to x592 (n=80). This is the real
compute lever — RB decouples calibration resolution from L, so you can run a
handful of posterior draws per fit where counting ranks need hundreds. It helps
*most* exactly in the informative-data regime where the prior-IS path fails.

**E3 — RB is NOT a free power boost (the corrective finding).** Once the counting
rank is tested against its correct discrete-uniform null (not the continuous one —
that mistake inflates small-L counting to spurious power=1.00), RB's per-simulation
detection power for a variance-component (sigma^2) bug is comparable to, even
slightly below, counting. The speed-up is L (per-fit cost), not S (number of fits).

**E4 — RB relocates what is under test (the honest caveat).** For a bug in the
fiber (alpha) update — credible intervals 20% too narrow — the counting rank
detects it clearly (D=0.065, out of band, classic S-shaped signature) while the
replace-RB PIT is *blind* (D=0.014, in band): it substitutes the correct analytic
conditional for the sampler's alpha draws, so it cannot see an error in those
draws. RB tests the base coordinates it conditions on (weighted by pi_j), not the
fiber it marginalises.

## Refined thesis for the paper

RB-PIT's value is (1) drastic within-sim variance reduction -> slash L, biggest
where pi_j is small; and (2) it directs sensitivity onto the base/hyperparameters,
blind to pure fiber-update errors. It is not a uniform power increase. Whether it
helps or hurts a given bug depends on where the bug sits relative to the base/fiber
split that pi_j quantifies. Use RB for base calibration; keep a fiber-sensitive
statistic for the fiber.

## Most promising next estimator (revealed by E4)

A **control-variate PIT** that keeps the sampler's fiber draws (so it still tests
them) but subtracts the analytic conditional CDF as a correlated control variate to
cut Monte Carlo noise. Naively it collapses to the counting rank, because the only
estimate of the control's mean is the same L draws. It pays off *only* under an
asymmetry: the analytic conditional CDF is cheap to evaluate, so its mean can be
pinned down on many more base draws than the expensive fiber draws — then the
residual (indicator − conditional CDF) carries the fiber-bug signal at reduced
variance. Worth prototyping next; it would recover E4 sensitivity without losing
the E2 gain.

## Caveats / to firm up

- The pooled KS test ignores within-simulation dependence across the J groups, so
  counting size runs a little high (0.03–0.12). Replace with the Säilynoja
  discrete-uniform simultaneous band per coordinate.
- Gaussian outcome makes the fiber conditional exact. For the non-linear outcome
  in the joint longitudinal GP model the conditional is only Laplace/EP-Gaussian,
  so RB there is approximate — which is also where the real calibration risk lives.

Run: `Rscript rb_pit_sbc_prototype.R` (base R only; ~5 min).
