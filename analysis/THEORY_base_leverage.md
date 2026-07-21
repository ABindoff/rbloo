# The base leverage: a connection-orthogonal decomposition of case-deletion influence

**Status:** derivation + numerical verification. Feeds the PSIS-LOO acceleration
question and, potentially, a repositioning of the fibr paper. Notation follows
`fibr_paper_jcgs.tex`. Verification script: `verify_base_leverage.R`.

---

## 1. Claim

The metric-orthogonal Ehresmann connection A = −G_FF⁻¹G_BF does more than set the
centring weight π_j. Applied to a single observation's score, it produces an
**exact orthogonal decomposition of that observation's generalized (case-deletion)
leverage** into a vertical (fiber) term and a horizontal (base) term:

    g_i' G⁻¹ g_i  =  g_i^F' G_FF⁻¹ g_i^F   +   g̃_i' M⁻¹ g̃_i
                     \___ vertical: pooling __/   \___ horizontal: base leverage __/

with g̃_i = g_i^B − G_BF G_FF⁻¹ g_i^F = g_i^B + A' g_i^F  and  M = G_BB − G_BF G_FF⁻¹ G_FB
(the Schur complement = marginal base precision). The horizontal term is the
"base-space pathology that π_j does not govern" (the paper's own phrase) — it is
π_j's missing partner, and it is built from the same connection.

---

## 2. Setup and the two objects (a distinction that matters)

Base θ (hyperparameters: σ_u, population β, GP hypers), fiber α (group
coordinates), metric G = [[G_BB, G_BF],[G_FB, G_FF]] with G_FF block-diagonal.
One observation i touches only its own group, so its score splits as
g_i = (g_i^B, g_i^F) with g_i^F nonzero only in the α_{j(i)} block.

There are **two** leverage-type quantities and they must not be conflated (this
was the correction numerical scrutiny forced):

- **Structural leverage** h_i = (Fisher info of obs i)/G_FF,j — residual-free,
  computable a priori from the design. Summed over a group it equals **1 − π_j**
  exactly (verified to 1e-16). This is fibr's pooling factor: the *capacity* for
  influence.
- **Influence** L_i = g_i' G⁻¹ g_i — the score-weighted, data-dependent version,
  = structural leverage × realized standardized-residual².

Both decompose by the identity above (structural via E[g g'] = the info of obs i;
influence via the realized g_i). π_j is the vertical *structural* term; the
horizontal partner exists at both levels.

**Empirical correction (see VALIDATION_RESULTS.md).** An earlier draft claimed the
*influence* predicts a fold's PSIS-LOO k̂. Model-1 validation overturns that: k̂ is
a tail/variance property, so the residual-free **structural** leverage (= the
pooling factor) is its natural predictor (Spearman up to 0.79 vs k̂_full), while
the **influence** predicts the deletion *magnitude* — the σ_u shift (cor 0.987) and
the elpd change. Net: structural leverage → IS reliability; influence → effect
size. The base *influence* term does predict the residual base-IS k̂ after RB
(Spearman 0.62–0.71 where the base IS is strained).

---

## 3. Why it is "already in the connection"

The identity is the Pythagorean theorem for the tangent-space split the connection
defines. g_i^F is the vertical component of the score; g̃_i = g_i^B + A′g_i^F is
its horizontal component — the base score with the fiber part **transported up
through A**. Their squared metric-norms (G_FF⁻¹ and M⁻¹ respectively) are the two
leverages, and they add to the total because horizontal ⟂ vertical in G. So the
base leverage is nothing but ‖g_i‖²_horizontal; the connection that gives the
centring weight is the *same* connection that lifts a case into the base.

The decomposition is purely algebraic (Schur complement) and holds regardless of
flatness. What flatness buys — your Prop. on the flat true connection — is that A
is path-independent, so the horizontal lift, hence the base leverage, is a
coherent global quantity rather than a holonomy-ambiguous one. And it must use the
**true** Fisher metric: computed with a sampler's working metric the split is the
curved one and the leverage is wrong — the same true-vs-working distinction the
paper already draws, reused.

---

## 4. Worked example: influence on σ_u flows entirely through A

Gaussian random intercept, base = (μ, τ=log σ_u), fiber = α, σ known. For
observation i in group j with residual r_i = y_i − μ − α_j:

    g_i^B = ( r_i/σ² , 0 ),   g_i^F = e_j · r_i/σ² .

The τ-component of the *direct* base score is **zero** — an observation carries no
first-order information about σ_u on its own. Yet its leave-one-out influence on
σ_u is real, and the decomposition says exactly where it comes from:

    g̃_i,μ  = (r_i/σ²) · π_j                         (grand mean sees it, shrunk by pooling)
    g̃_i,τ  = (r_i/σ²) · 2 α_j e^{−2τ} / G_FF,j   ∝  r_i · α_j / (σ_u² G_FF,j)

The entire σ_u-influence is the connection term A′g_i^F — proportional to the
observation's residual **times its group effect α_j**. Large exactly when a
high-residual point sits in an extreme group; amplified by M⁻¹_{ττ} ~ 1/J (few
groups). This is the singleton/outlier-group pathology, derived, not assumed.

---

## 5. The PSIS-LOO / k̂ bridge

In RB-LOO you marginalise the fiber analytically and importance-sample only the
base. The log-weight for dropping obs i is that observation's *fiber-integrated*
score contribution; its variance under the base posterior (covariance M⁻¹) is
exactly g̃_i' M⁻¹ g̃_i = the base influence. So:

    base influence  =  variance of the RB-LOO base importance weights  =  k̂ predictor.

That gives a two-number, closed-form triage per observation:

| fiber term (1−π_j locally) | base term (horizontal) | verdict |
|---|---|---|
| small | small | plain PSIS-LOO fine |
| large | small | **RB-LOO's sweet spot** — the fiber marginalisation cures it |
| any | large | refit; RB won't save it (outlier-in-small-group folds) |

The pooling factor says how much RB helps; the base leverage says whether what's
left is stable.

---

## 6. Verification (verify_base_leverage.R)

- **Part A (identity):** random SPD metric, |LHS − RHS| = 4.4e-16. The
  decomposition is exact.
- **Part B (Gaussian random-intercept LMM, 12 groups incl. singletons + 2 outlier
  groups):**
  - one-step base influence predicts the **actual** σ_u leave-one-out shift (full
    refit) at **cor 0.987**; base influence vs |Δτ| at 0.966.
  - structural fiber leverage summed per group equals **1 − π_j to 1e-16**.
  - base leverage concentrates on near-singleton/outlier groups and its share of
    total leverage falls monotonically with group size — the derived pathology.
- **Known caveat, visible in the scatter:** the one-step (linearised) influence
  *under*-shoots the two largest σ_u shifts. First-order case-deletion is a floor,
  not exact, for the highest-influence points — the same Laplace boundary as the
  change-point in Phase 2. Refit or a second-order correction is needed there.

---

## 7. Prior-art scrutiny (required before any novelty claim)

**The adversarial search has now been run — see `PRIOR_ART_scrutiny.md`. It is
sobering.** The generalized leverage g' G⁻¹ g and its Schur decomposition are
classical case-influence machinery (Cook 1986; Wei, Hu & Fung 1998); the
*hierarchical* leverage/DOF partition is anticipated by Hodges & Sargent (2001),
Cui-Hodges-Kong-Carlin (2010, "Partitioning Degrees of Freedom…"), and — most
dangerously — Lovison (2025, "Augmented Hat-Matrix of Hierarchical GLMs"), which I
could not read and which may already give the fixed/random split; the RB-LOO cure
is anticipated by the integrated / non-factorized / leave-group-out LOO line
(Bürkner-Gabry-Vehtari; Merkle et al.; Liu & Rue); and "influential = high k̂" is
loo folklore. What plausibly survives, pending a full read of Lovison and
Cui-Hodges:

1. **Framing the case-influence split as the horizontal/vertical decomposition of
   the metric-orthogonal Ehresmann connection** of the hierarchical bundle — i.e.
   identifying the base leverage as A′g_i^F transported through *your* connection,
   the same object that gives the centring weight. The bundle-geometric reading is
   the paper's, not the influence literature's.
2. **π_j as the vertical/structural half of a single influence identity whose
   other half is the variance-component (funnel) leverage** — unifying the two
   pathologies the paper currently presents as separate.
3. **The operational payoff: base leverage = RB-LOO base-IS weight variance**, a
   closed-form a-priori predictor of PSIS-LOO reliability and a triage rule. This
   application is not in the influence-geometry literature.

Do NOT claim the bare decomposition; claim the connection framing, the
two-pathologies unification, and the LOO application, and cite the case-influence
precedents explicitly.

---

## 8. What this could do for the fibr paper

The current paper's honest position (per the priority reviews) is under some
pressure: the pooling-factor-as-centring-weight is anticipated (Tan & Nott 2013),
the flat connection has precedents (Amari 2001; Cox & Reid 1987), and the paper's
own conclusion is slightly deflationary — "π_j is a diagnostic, not a faster
parameterisation," and the funnel is "a separate pathology π_j does not govern."

The base leverage changes that arc:

- It gives the connection a **second, less-anticipated job** beyond the centring
  weight: generating the base-space influence — precisely the pathology the paper
  flags as outside π_j's reach. The connection *closes its own gap*.
- It **unifies the two pathologies** into one orthogonal identity: partial pooling
  (vertical) and the variance-component funnel (horizontal) are the two halves of
  a single case-influence decomposition in the bundle metric. That is a structural
  result, not a reframing.
- It gives π_j an **operational role** — PSIS-LOO reliability and the RB-LOO cure —
  answering the paper's own "what is the pooling factor for?" with something
  sharper than a diagnostic flag.

A plausible repositioning: *from* "a flat connection whose leftover is the
(known) pooling factor" *to* "the hierarchical bundle's metric-orthogonal split
decomposes case-deletion influence into a vertical pooling term and a horizontal
base-leverage term; this unifies partial pooling with the variance-component
funnel and yields a closed-form predictor of cross-validation reliability."

**Post-read verdict (`PRIOR_ART_scrutiny.md`, both papers now read in full).** The
*structural* leverage half is spoken for: Lovison (2026) is the hierarchical
augmented hat matrix — the subject-level hat value *is* the structural fiber
leverage — and Cui-Hodges (2010) partitions component complexity and already cites
Gelman-Pardoe pooling. So do **not** pitch this as a new leverage/hat-matrix
decomposition; that collides head-on with Lovison (2026). But both papers hold the
variance components fixed (Lovison) or measure complexity not influence
(Cui-Hodges), and **neither touches case-deletion influence on the variance
components, nor cross-validation reliability.** Those are the real opening. Reframe
the contribution as a **CV-reliability result** — the pooling factor predicts where
PSIS-LOO fails, the connection-generated base leverage predicts where the RB cure
still fails, RB-LOO fixes the rest (all validated) — with the base–fiber geometry as
the *explanation*, not the headline. For the fibr paper: at most a short addition
tying the connection to CV reliability, citing Lovison (2026) and Cui-Hodges (2010)
to demarcate. The higher-upside path is a separate focused methods paper on the
CV-reliability result.
