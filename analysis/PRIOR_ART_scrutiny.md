# Prior-art scrutiny (post-read): the base-leverage / connection result

Adversarial search + full read of the two closest papers (Lovison 2026;
Cui-Hodges 2010). **Verdict: the structural-leverage half is now firmly spoken
for — Lovison (2026) owns the hierarchical hat matrix — but the case-deletion
influence on the variance components via the connection, and the whole
cross-validation-reliability application, are untouched by this literature. Frame
the contribution around those, NOT around a new leverage/hat-matrix.**

## Lovison, Sciandra, Albano & Di Maria (2026), Int. Stat. Rev., doi:10.1111/insr.70030

Read in full. The **augmented hat matrix** for HGLMs, via the Lee–Nelder
h-likelihood "augmented GLM": augmented response y_a=(y,ψ), design T=[[X,Z],[0,I]],
GLS projector H*_a = Σ_a^{-1/2}T(TᵀΣ_a^{-1}T)^{-1}TᵀΣ_a^{-1/2}, block-diagonal into
subject-level (n×n) and cluster-level (r×r) hat values. New thresholds t_s=2h̄',
t_c=2h̄ for high-leverage subjects/clusters.

What this means for us:
- **The subject-level hat value IS my *structural* fiber leverage** (the
  shrinkage-weighted observation leverage; their σ_u²→0 limit → h̄=1 is exactly the
  pooling behaviour). So "structural fiber leverage ≈ 1−π_j" is **fully
  anticipated** — do not claim it. Predecessors: Demidenko & Stukel (2005), Nobre &
  Singer (2011), Loy & Hofmann (2013), Wei-Hu-Fung (1998).
- BUT Lovison is **leverage only, explicitly not influence** — they *adjust Y to
  keep modified points non-influential*. No case-deletion.
- Variance components (φ, λ) are **held known / plugged in** (Cases 1–2). The hat
  matrix is over (β, b) *given* the variance components. So the influence **on** the
  variance components — my base leverage — is **outside their object entirely.**
- No pooling-factor naming, no Gelman-Pardoe, no cross-validation, no PSIS/k̂, no
  Ehresmann connection, no Schur complement (their block-diagonality comes from
  Σ_a being block-diagonal, not a Schur reduction of a hyperparameter block).

## Cui, Hodges, Kong & Carlin (2010), Technometrics 52(1):124–136

Read in full. **Partitions total DF (model complexity / effective #params)** across
model *effects*: DF(X_2j) = tr[X_2jΓ_2jX_2jᵀ(X_1Γ_1X_1ᵀ+X_2Γ_2X_2ᵀ+Γ_0)^{-1}], the
ratio of an effect's modelled variance to total variance; DF(δ)+DF(ξ)+DF(ε)=total.

What this means for us:
- It is an **effect-level complexity** measure (a variance-ratio-weighted trace),
  **not** a per-observation case-deletion influence and **not** a leverage for
  outlier detection. It answers "how many effective parameters does clustering
  use", not "which observation destabilises σ_u".
- It **explicitly cites Gelman-Pardoe (2006)** pooling ("reminiscent… although
  their intent and results were rather different"). So the pooling ↔
  component-complexity bridge is partly made — the pooling-factor-is-interesting
  claim is not virgin territory, though the CV use is.
- The DF is a function of the variance ratios (hence shrinkage/pooling), like π_j.
- No case-deletion, no CV, no PSIS/k̂, no connection, no per-observation influence.

## Net demarcation — what is anticipated vs what survives

**Anticipated (do not claim):**
- Structural leverage of (β, b) given variance components, incl. the shrinkage
  reading and cluster-level leverage — Lovison (2026) and predecessors.
- Pooling factor ↔ effective-complexity/DF of a component — Cui-Hodges (2010),
  Gelman-Pardoe (2006).
- Integrated/marginal LOO cure — Bürkner-Gabry-Vehtari; Merkle; Liu & Rue.
- Case-deletion IS heavy tails / influential = high k̂ — Peruggia (1997); Vehtari.

**Untouched by this literature (the real opening):**
1. **Case-deletion influence ON the variance components (base leverage), via the
   Ehresmann connection.** Lovison conditions on the variance components; Cui-Hodges
   does complexity, not influence. The object g̃_iᵀM⁻¹g̃_i with g̃_i = g_i^B + A′g_i^F —
   how deleting an observation propagates through the connection into the
   hyperparameter posterior — is in neither. This is where the singleton/funnel
   pathology lives, and it is the piece our σ_u-shift verification (cor 0.987)
   targets.
2. **The cross-validation-reliability application.** Nobody links the pooling factor
   to PSIS-LOO k̂ reliability, nor the connection-generated base leverage to
   RB-LOO's residual base-IS stability, nor uses either to triage/cure LOO. This is
   the freshest and most defensible angle, and it is completely open.
3. **The unifying Ehresmann-connection framing** tying pooling (vertical), the
   variance-component funnel (horizontal), case influence, and CV reliability into
   one orthogonal identity on the base–fiber bundle.

## Recommendation (changed by the read)

Do **not** pitch this as "a new hierarchical leverage / hat-matrix decomposition" —
that now collides head-on with Lovison (2026), published in the same window. Pitch
it as **a cross-validation-reliability result**: *the pooling factor predicts where
PSIS-LOO fails, the connection-generated base leverage predicts where the RB cure
still fails, and RB-LOO fixes the rest* — validated (RB-LOO 3–5× more accurate than
PSIS-LOO vs refit; k̂ tracks structural leverage at ρ≈0.79; base leverage tracks
base-k̂ at ρ≈0.62–0.71). The geometry (base–fiber bundle, connection) is the
*explanation* of why, and the tie to the fibr paper; it is not the headline.

For the fibr paper itself: a *short* addition — "the connection also generates a
case-deletion influence on the hyperparameters, giving the pooling factor an
operational role in cross-validation reliability" — cite Lovison (2026) and
Cui-Hodges (2010) to demarcate, and keep the leverage claim to the *influence /
variance-component / CV* part that they do not cover. A separate focused methods
paper on the CV-reliability result is the higher-upside path.
