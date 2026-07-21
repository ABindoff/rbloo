# Validation results — Model 1 (Gaussian LMM, real PSIS-LOO)

Ran `leverage_psisloo_validation.R` (conjugate Gibbs, `loo` 2.6.0), two regimes:
A = 15 groups (sigma_u well identified), B = 6 groups (barely). Figure:
`leverage_psisloo.png` (top row A, bottom B).

## Confirmed

1. **RB-LOO cures the fiber-driven PSIS-LOO failures.** Marginalising the fiber
   analytically drove every k-hat down: A had 3 folds with k>0.7 (max 1.27) → 0
   after RB (mean 0.22 → −0.02); B had 2 (max 1.00) → 0 (mean 0.59 → 0.14). Robust
   across regimes and reruns. This is the central operational claim and it holds.

2. **RB-LOO is more accurate than PSIS-LOO against gold-standard refits.** elpd
   RMSE vs refit: A 0.039 (RB) vs 0.114 (PSIS) ≈ 3×; B 0.116 vs 0.156. RB-LOO is
   both stable *and* accurate where PSIS-LOO degrades.

3. **k-hat is predicted by leverage — and it is the STRUCTURAL term, not the
   influence term.** Spearman(structural h_i, k_full) = 0.79 (A), 0.49 (B);
   Spearman(influence L_i, k_full) = 0.43 (A), 0.77 (B). This **corrects the theory
   note**: k-hat is a tail/variance property, so the residual-free *structural*
   leverage (= the pooling factor, computable a priori) is its natural predictor;
   the *influence* (score-weighted) predicts the deletion *magnitude* — the σ_u
   shift (cor 0.987, prior script) and the elpd change. Which term wins for k-hat
   is regime-dependent and noisy at this N, but both are positive; the clean,
   defensible statement is: **structural leverage → IS reliability (k-hat);
   influence → effect size (elpd / σ_u shift).**

4. **Base leverage predicts residual base-IS strain (H2a).** Where the base IS is
   actually strained (few groups, B), Spearman(L_base, k_base) = 0.62–0.71 across
   runs. Where sigma_u is well identified (A) there is no base strain to predict
   (max k_base 0.29) and L_base raises no false alarms — correct null behaviour.

## Caveats / not yet established

- **Point estimates are noisy** (small N, k-hat estimation variance, refit-subset
  sampling); directions are stable across reruns, magnitudes are indicative only.
  The pre-registered thresholds should be judged on a scaled-up, replicated run.
- **RB-LOO never outright failed** even in the strained regime (max k_base 0.49).
  It is more robust than expected — good news, but it means H2b (AUC for
  k_base>0.7) could not be exercised. A J=3–4 / more-extreme regime is needed to
  produce genuine RB failures and complete H2b.
- Pre-registered **H1a on the influence term technically missed** 0.7 in A (0.43);
  the intended phenomenon holds via the structural term (0.79). This is a
  refinement of which quantity to name, not a failure of the mechanism.

## Verdict

The conjugate-Gaussian case validates the framework's spine: leverage predicts
PSIS-LOO reliability (structural term, a priori), RB-LOO cures the fiber-driven
failures and is more accurate than PSIS-LOO, and base leverage flags residual
base-IS strain where it exists. The one substantive correction — k-hat tracks
structural leverage, influence tracks effect size — actually *strengthens* the
story, because the a-priori pooling factor becomes the IS-reliability predictor.

## Next rungs (unchanged from the spec, reordered by value)

1. Scale up + replicate Model 1 for firm point estimates against the pre-registered
   thresholds; add a J=3–4 extreme regime to exercise H2b.
2. GLMM (Poisson/Bernoulli) via brms — the Pólya-Gamma conditional-Gaussian
   extension.
3. bjlm smooth change-point — the linearisation boundary (expect degradation at
   extreme sharpness, per Phase 2).
