# Gate G3: base-leverage triage (C3) — RESOLVED by SIMPLIFICATION

Ran `p3_base_leverage.R`: few-groups sweep (240 datasets, 1600 folds, base IS
genuinely strained) + eight schools. The result is a productive negative: **C3 as a
distinct diagnostic does not earn its place — fold it into C1.** The paper gets
simpler and more honest.

## What we asked and found

**H4 (base leverage predicts residual RB-LOO failures).**
- As a *classifier*: AUC(k_base>0.7 ~ L_base) = **0.90** — base leverage does flag
  the folds RB-LOO can't fix. (Continuous Spearman 0.46, just under the 0.50 bar —
  k_base is floored on the ~98% of folds RB fixes, diluting the rank correlation.)
- **The decisive test — does L_base beat the pooling factor?** AUC(k_base>0.7):
  **L_base 0.90 vs h_struct (structural/pooling leverage) 0.92** vs (1−π_j) 0.66.
  The structural leverage that C1 *already uses* predicts RB-failure **as well or
  better** than the elaborate connection-based base leverage. Conditional on high
  fiber leverage (the folds RB is applied to), Spearman(L_base,k_base)=0.48 vs
  Spearman(h_struct,k_base)=0.40 — a marginal, not decisive, edge.

**Conclusion:** the residual-triage *need* is real and *met* — but by the pooling
factor, not by a separate quantity. In realistic hierarchical data, fiber leverage
and base-identification strain co-occur (few groups ⇒ small n ⇒ high h_struct *and*
high L_base), so the simpler predictor suffices.

## What this does to the paper (it improves it)

- **Drop C3 as a separate contribution.** The triage is *one* quantity, the pooling
  factor, used at *two* thresholds: moderate → apply RB-LOO; extreme (data-driven
  coordinate in a poorly-identified base) → refit. One predictor, one cure.
- **The base leverage / connection survives as the EXPLANATION (C4), not a
  diagnostic.** T1–T2 are still true and verified (base leverage = base-IS weight
  variance; predicts the σ_u LOO shift at cor 0.987). They explain *why* the
  few-groups folds resist RB — the horizontal (base) component of influence is large
  there — even though the pooling factor is the sufficient practical proxy. State
  L_base as the theoretically-correct predictor and π_j as the operational one.
- Net contributions become **C1 (pooling factor = two-level LOO-reliability triage),
  C2 (RB-LOO cure), C4 (base–fiber geometry explains both failure modes)** — cleaner
  than the original four.

## Eight schools (honest)

RB-LOO lowers every school's k̂ (max 0.58 → 0.39), but full PSIS-LOO did **not**
outright fail (all k̂<0.7) with clean independent draws — consistent with the high
pooling factors (π_j ≈ 0.88–0.97) that C1 says make PSIS *borderline, not broken*.
The eight-schools LOO "problem" is partly a sampler-mixing artifact; with good draws
it is borderline. So eight schools is a *"theory correctly predicts borderline"*
case, not a rescue showcase. The rescue evidence is the simulation (526 k̂>0.7 folds,
100% cured, gate G1) and, later, roaches (P4/P5).

## Gate G3 status

**RESOLVED — C3 folded into C1.** No gate failure in the damaging sense: the triage
works; it just needs one quantity, not two. The paper is C1 + C2 + C4, all
validated. Update `EXPERIMENT_PLAN_cv_reliability.md`: C3 → absorbed into C1; base
leverage moves from "diagnostic" to "explanation (C4)"; H4 retired in favour of the
head-to-head (pooling factor is the sufficient triage predictor).

## Caveat worth one more experiment

L_base *is* the theoretically-correct base-IS-failure predictor; my sweeps couple
fiber leverage and base strain, so they can't show a regime where L_base beats π_j.
If such a regime exists (fiber well-pooled but base badly identified — e.g. many
small-but-not-singleton groups with a near-degenerate variance component), L_base
would separate from π_j there. Worth one targeted probe before fully retiring L_base
as a predictor — but it is not needed for the paper, which stands on π_j.
