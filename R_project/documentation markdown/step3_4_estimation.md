# Steps 3–4: Demand Estimation, Markups, and Marginal Costs
## Miller & Weinberg (2017) — Replication in R

---

## 1. Model Specification

### 1.1 Estimating equations

**Simple logit (Berry 1994)**

$$\ln(s_{jt}) - \ln(s_{0t}) = \alpha \, p_{jt} + \mathbf{x}_j' \boldsymbol{\beta} + \xi_{jt}$$

**Nested logit (Berry 1994; MW2017 eq. 7)**

$$\ln(s_{jt}) - \ln(s_{0t}) = \alpha \, p_{jt} + \sigma \ln(s_{j|g,t}) + \mathbf{x}_j' \boldsymbol{\beta} + \xi_{jt}$$

where $s_{j|g,t} = s_{jt} / \sum_{k \in g} s_{kt}$ is the within-group conditional share and $\sigma \in [0,1)$ is the nesting parameter. All inside goods form a single nest; the outside good is the only alternative group.

**Paper's full model: Random Coefficients Nested Logit (RCNL)**

$$U_{ijt} = \delta_{jt} + \mu_{ijt} + \epsilon_{ijt}$$

where $\delta_{jt} = \alpha \, p_{jt} + \mathbf{x}_j' \boldsymbol{\beta} + \xi_{jt}$ is the mean utility and $\mu_{ijt} = \mathbf{x}_{jt}' \Pi \, d_i$ captures consumer heterogeneity through income interactions with price, a constant, and calories. Estimated by GMM with 500 simulated consumers per market.

### 1.2 Dependent variable

$$y_{jt} = \ln(s_{jt}) - \ln(s_{0t})$$

- $s_{jt} = q_{jt} / \text{msize}_t$ — raw share (not normalized).
- $s_{0t} = 1 - \sum_j s_{jt}$ — outside good share.
- Range in data: $[-5.38, -1.73]$. All finite, no zeros.

### 1.3 Endogenous variables

| Variable | Column | Endogeneity source |
|----------|--------|--------------------|
| $p_{jt}$ | `price` | Simultaneous: firms observe $\xi_{jt}$ when setting prices |
| $\ln(s_{j|g,t})$ | `logcondshr` | Mechanical function of equilibrium shares (nested logit only) |

### 1.4 Instrument sets

**Baseline (Spec D): 5 excluded instruments**

| ID | Variable | Type | Identifies |
|----|----------|------|-----------|
| Z1 | `dist` (miles × diesel) | Cost-side | $\alpha$ (price) |
| Z2 | `coalpost` (ABI+MC post-merger) | Cost-side / structure | $\alpha$ (price) |
| Z4 | `sum_dist` (market-total distance) | BLP-style | $\sigma$ (nesting) |
| Z7 | `sumdist_abi` (sum_dist × ABI) | BLP-style interaction | $\sigma$ (nesting) |
| Z8 | `sumdist_mc` (sum_dist × MC) | BLP-style interaction | $\sigma$ (nesting) |

**Extended (rejected): adds Z9–Z12**

| ID | Variable | Type | Status in small sample |
|----|----------|------|----------------------|
| Z9 | mean income × 1 | Demographics | Produces wrong-sign alpha |
| Z10 | mean income × calories | Demographics | Sargan rejects |
| Z11 | mean income × size | Demographics | Identifies off perturbation |
| Z12 | mean income × import | Demographics | Excluded from baseline |

**Degenerate (dropped): Z3, Z5, Z6** — zero variation because $J=8$ in all markets.

### 1.5 Fixed effects

| Set | Dummies | Role |
|-----|---------|------|
| Product FE (`prodfe`) | 8 | Absorb brand × size mean utility |
| Date FE (`datefe`) | 15 (16 periods, 1 reference) | Absorb common time shocks |

City FEs are **not** included in the demand specification — cross-city price variation identifies $\alpha$. City FEs enter the supply-side cost regression only.

---

## 2. Estimation Results

### 2.1 Comparison table

| | OLS (NL) | Spec B (logit) | **Spec D (NL)** | Paper (RCNL) |
|---|---|---|---|---|
| Model | Nested logit | Simple logit | **Nested logit** | RCNL |
| Method | OLS | 2SLS | **2SLS** | GMM |
| Excluded IVs | 0 | 5 | **5** | 12 |
| $\hat{\alpha}$ | +0.025 | **−0.059** | **−0.027** | **−0.109** |
| (clustered se) | (0.062) | (0.070) | **(0.047)** | — |
| $\hat{\sigma}$ / $\hat{\rho}$ | 0.975 | — | **0.940** | **0.778** |
| (clustered se) | (0.048) | — | **(0.107)** | — |
| Mean own-price elasticity | — | −0.72 | **−5.01** | **−4.65** |
| 1st-stage F (price) | — | 15.7 | **15.7** | — |
| 1st-stage F (logcondshr) | — | — | **21.3** | — |
| Wu-Hausman p | — | — | 0.085 | — |
| Sargan p | — | 0.000 | 0.003 | — |
| N | 640 | 640 | **640** | 94,656 |

### 2.2 Key coefficients

- **$\hat{\alpha} = -0.027$** (Spec D): correct sign, statistically insignificant at 5% but economically meaningful. Attenuated relative to the paper's $-0.109$, consistent with the simpler specification and smaller perturbed sample.
- **$\hat{\sigma} = 0.940$** (Spec D): within $[0,1)$, indicating strong within-nest correlation. Higher than the paper's $\hat{\rho} = 0.778$, implying less substitution to the outside good in our estimate.
- **OLS bias**: OLS produces $\hat{\alpha} > 0$ (wrong sign) in the nested logit — classic upward simultaneity bias.

### 2.3 First-stage diagnostics

- Price first-stage F = 15.7: above the Staiger-Stock threshold of 10. Strongest contributors: `coalpost` (t = 9.6), `sum_dist` (t = −3.1).
- Logcondshr first-stage F = 21.3: comfortably strong. BLP-style instruments (sum_dist and interactions) are the primary drivers.
- Sargan rejects at p = 0.003: expected given data perturbation (see Section 3 below).

---

## 3. IV Sensitivity Analysis

### 3.1 Specification comparison

| Spec | Model | IVs | $\hat{\alpha}$ | Sign OK? | $\hat{\sigma}$ | $\sigma \in [0,1)$? | F(price) | Sargan p |
|------|-------|-----|----------------|----------|-----------------|---------------------|----------|----------|
| A | Logit | 2 (cost) | +0.097 | No | — | — | 30.6 | 0.273 |
| B | Logit | 5 (cost+BLP) | −0.059 | Yes | — | — | 15.7 | 0.000 |
| C | NL | 2 (just-id) | −0.715 | Yes | 8.128 | No | 30.6 | — |
| **D** | **NL** | **5 (cost+BLP)** | **−0.027** | **Yes** | **0.940** | **Yes** | **15.7** | **0.003** |
| E | Logit | 9 (+income) | +0.148 | No | — | — | 41.7 | 0.000 |

### 3.2 Why income instruments fail

- The replication data was intentionally perturbed: random noise on prices/quantities and 1.5× inflation of ABI/Miller volumes (`smalldata.m` lines 31–38).
- With only 5 cities, market-level income variation is minimal. The instruments identify off perturbation artifacts rather than true demand heterogeneity.
- The first-stage F rises from 15.7 to 41.7 — spurious relevance that does not imply validity.

### 3.3 Why Spec D is preferred

Spec D is the unique specification satisfying all four criteria simultaneously:

1. $\hat{\alpha} < 0$ (correct sign on price)
2. $\hat{\sigma} \in [0, 1)$ (valid nesting parameter)
3. Both first-stage F-statistics exceed 10
4. Instruments have clear economic justification (cost shifters + market-structure variation)

The Sargan rejection (p = 0.003) is noted but attributed to the data perturbation rather than a fundamental identification failure.

---

## 4. Elasticities

### 4.1 Own-price elasticities

| Model | Mean | Range |
|-------|------|-------|
| Spec B (simple logit) | −0.72 | [−1.03, −0.38] |
| **Spec D (nested logit)** | **−5.01** | **[−7.52, −2.42]** |
| **Paper RCNL** | **−4.65** | **[−7.37, −2.46]** |

The simple logit produces unrealistically low elasticities (the IIA problem). The nested logit and RCNL are in close agreement: mean own-price elasticities around −4.7 to −5.0, with a range spanning −2.4 to −7.5.

### 4.2 Own-price elasticities by product (RCNL, 80 markets)

| Product | Mean | Range |
|---------|------|-------|
| Bud Light 12pk | −5.05 | [−6.65, −3.26] |
| Bud Light 24pk | −4.24 | [−6.02, −2.82] |
| Coors Light 12pk | −4.33 | [−5.39, −3.54] |
| Coors Light 24pk | −3.60 | [−4.87, −2.46] |
| Corona 12pk | −5.08 | [−7.13, −3.01] |
| Heineken 12pk | −4.94 | [−6.72, −2.60] |
| Miller Lite 12pk | −5.37 | [−7.35, −3.27] |
| Miller Lite 24pk | −4.60 | [−6.05, −3.23] |

Pattern: 12-packs are more elastic than 24-packs within the same brand. Import brands (Corona, Heineken) have elasticities comparable to domestic 12-packs.

### 4.3 Cross-elasticity patterns (RCNL, market 1)

| From \ To | Bud 24 (same brand) | Coors 12 (same size) | Corona 12 (import) |
|-----------|-------|--------|--------|
| Bud Light 12pk | 0.78 | 0.31 | 0.28 |

Cross-elasticities are highest within the same brand (different size), moderate across domestic brands at the same size, and lowest between domestic and import products. This is the differentiated substitution pattern that BLP/RCNL is designed to capture and that the simple logit misses entirely (IIA).

---

## 5. Markups and Marginal Costs

### 5.1 Bertrand markups (RCNL parameters, FOC inversion)

Marginal costs are recovered by inverting the Bertrand first-order condition:

$$mc_j = p_j - \underbrace{\left[ -(\mathbf{\Omega} \circ \frac{\partial \mathbf{s}}{\partial \mathbf{p}}')^{-1} \mathbf{s} \right]_j}_{\text{markup}_j}$$

where $\Omega_{jk} = 1$ if products $j$ and $k$ are owned by the same firm.

| Firm | Price | Markup | MC | Lerner |
|------|-------|--------|-----|--------|
| ABI | $13.63 | $4.33 | $9.30 | 32.0% |
| Corona | $13.90 | $3.08 | $10.83 | 22.2% |
| Heineken | $13.96 | $3.20 | $10.76 | 23.0% |
| Coors | $8.96 | $2.73 | $6.23 | 31.0% |
| MillerCoors | $12.08 | $3.47 | $8.61 | 29.7% |

- ABI has the highest Lerner index (32%), consistent with its dominant domestic position and multi-product portfolio.
- Imports have lower Lerner indices (22–23%) despite higher prices, reflecting higher marginal costs (transportation, import duties).
- All implied marginal costs are positive — no negative MC pathology.

### 5.2 Validation against Matlab (`sres_bertrand.mat`)

| Metric | Value |
|--------|-------|
| Correlation (R vs Matlab MC) | 0.982 |
| Mean |difference| | $0.32 |
| Max |difference| | $2.74 |
| Negative MC | 0 (both) |

The residual difference arises because the Matlab code runs a 3-step coalition/fringe FOC procedure (`f_impute_mc.m`) with derivative recomputation via 500-draw simulation at each step, while our R implementation evaluates the full ownership FOC directly. The 0.98 correlation confirms economic equivalence.

### 5.3 Ownership matrices

| Period | Firms | Co-owned pairs (of 28) | Key feature |
|--------|-------|----------------------|-------------|
| Pre-merger (FY 2006–07) | ABI, Coors, Corona, Heineken, MillerCoors | 3 | Coors and Miller are separate firms |
| Post-merger (FY 2010–11) | ABI, Corona, Heineken, MillerCoors | 7 | Coors Light joins Miller Lite under MillerCoors |

The merger increases co-owned pairs from 3 to 7, creating a 4-product firm (MillerCoors) that internalizes pricing externalities across Coors Light and Miller Lite in both sizes. This structural change is the basis for the merger simulation in Step 5.

---

## 6. Interpretation and Limitations

### 6.1 Small sample vs full data

- The paper estimates RCNL on 94,656 observations (37 regions, 39 products). We have 640 observations (5 cities, 8 products).
- In the full data, $J$ varies across markets, providing identification for the nesting parameter via BLP-style instruments. In our sample, $J = 8$ everywhere, shutting off this channel.
- Random coefficient parameters ($\Pi$) cannot be credibly estimated with 5 cities of income variation. The paper's RCNL leverages 37 regions.

### 6.2 Data perturbation

- The replication sample was constructed by adding random noise to prices and quantities and inflating ABI/Miller volumes by 1.5× (`smalldata.m`).
- This perturbation breaks the exclusion restrictions that income-based instruments rely on, causing Sargan rejection and sign reversals in expanded IV specifications.
- The pre-computed derivative and elasticity matrices (`derMat_2`, `elasMat_2`) reflect the perturbed data, so our supply-side analysis is internally consistent.

### 6.3 Differences from paper estimates

| Parameter | Our 2SLS (Spec D) | Paper RCNL | Explanation |
|-----------|-------------------|------------|-------------|
| $\alpha$ | −0.027 | −0.109 | Attenuation: simpler model, smaller sample, perturbation |
| $\sigma / \rho$ | 0.940 | 0.778 | Higher nesting: fewer products reduce outside-good substitution |
| Mean own-elas | −5.01 | −4.65 | Broadly consistent despite parameter differences |

The supply-side analysis (markups, merger simulation) uses the paper's pre-estimated RCNL parameters from `dresgmm2.mat`, not our 2SLS estimates. This is the intended workflow: the replication package itself takes demand as given.

---

## 7. Takeaways

- The nested logit 2SLS with cost-side and BLP-style instruments (Spec D) is the only specification producing both $\hat{\alpha} < 0$ and $\hat{\sigma} \in [0,1)$ in this sample. It serves as the baseline demand estimate for the problem set.
- Income-based instruments, despite high first-stage F-statistics, produce wrong-sign coefficients and fail overidentification tests. This is a direct consequence of the data perturbation and limited cross-sectional variation (5 cities), not a flaw in the research design.
- Own-price elasticities from the nested logit (−5.0) and RCNL (−4.7) are in close agreement, both far more plausible than the simple logit (−0.7). The nested structure is essential for generating realistic substitution patterns in differentiated product markets.
- Bertrand markups average $3.51 (Lerner index 28.7%), with ABI earning the highest margin (32%). All implied marginal costs are positive and economically sensible.
- The RCNL cross-elasticity matrix reveals that within-brand substitution (e.g., Bud 12pk → Bud 24pk: 0.78) exceeds cross-brand substitution at the same size (Bud 12pk → Coors 12pk: 0.31), which in turn exceeds domestic-to-import substitution (Bud 12pk → Corona 12pk: 0.28). This differentiated pattern is the economic foundation for the merger analysis.
- Results should be interpreted as qualitative replication. Exact numerical agreement with the paper is not achievable — nor expected — given the perturbed small sample.
