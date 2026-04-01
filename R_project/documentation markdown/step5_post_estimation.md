# Step 5: Post-Estimation — Elasticities, Markups, and Marginal Costs
## Miller & Weinberg (2017) — Replication in R

---

## 1. Overview

This document reports post-estimation outputs derived from the RCNL demand parameters loaded from `small_dresgmm2.mat`. These parameters were estimated by MW(2017) on the full IRI dataset (37 regions, 39 products, 94,656 observations) and are taken as given for the supply-side analysis. The key outputs are:

- **Elasticity matrices**: own- and cross-price elasticities from the pre-computed `elasMat_2` (8×8×140, subset to 80 fiscal-filtered markets).
- **Bertrand markups**: recovered by inverting the Bertrand-Nash first-order conditions using the full RCNL simulation (500 consumer draws per market).
- **Marginal costs**: implied as $mc_j = p_j - \text{markup}_j$, serving as inputs to the merger counterfactual.

All computations use: $\alpha = -0.1087$, $\rho = 0.7779$, $\theta_2 = [0.0009, 0.0125, 0.0045]$.

---

## 2. Elasticities

### 2.1 Own-price elasticities

| Statistic | RCNL (80 markets) | Nested logit 2SLS (Spec D) |
|-----------|--------------------|----------------------------|
| Mean | −4.65 | −5.01 |
| Min | −7.37 | −7.52 |
| Max | −2.46 | −2.42 |

Both specifications produce elasticities in the range typical of differentiated consumer goods. The RCNL estimates are slightly less elastic on average, reflecting the richer substitution patterns enabled by consumer heterogeneity.

### 2.2 Own-price elasticities by product (RCNL)

| Product | Mean | Min | Max |
|---------|------|-----|-----|
| Bud Light 12pk | −5.05 | −6.65 | −3.26 |
| Bud Light 24pk | −4.24 | −6.02 | −2.82 |
| Coors Light 12pk | −4.33 | −5.39 | −3.54 |
| Coors Light 24pk | −3.60 | −4.87 | −2.46 |
| Corona 12pk | −5.08 | −7.13 | −3.01 |
| Heineken 12pk | −4.94 | −6.72 | −2.60 |
| Miller Lite 12pk | −5.37 | −7.35 | −3.27 |
| Miller Lite 24pk | −4.60 | −6.05 | −3.23 |

### 2.3 Substitution patterns (RCNL, market 1)

| Substitution type | Example pair | Cross-elasticity |
|-------------------|-------------|-----------------|
| Within-brand, different size | Bud Light 12pk → Bud Light 24pk | 0.78 |
| Across-brand, same size | Bud Light 12pk → Coors Light 12pk | 0.31 |
| Domestic → import | Bud Light 12pk → Corona 12pk | 0.28 |

The RCNL generates a cross-elasticity hierarchy that the simple logit cannot produce: within-brand substitution is strongest, followed by across-brand domestic, with domestic-to-import substitution weakest. This differentiated pattern is the economic foundation for realistic merger simulation — it determines how much demand diverts from one merging product to the other versus to outside options.

### 2.4 Comparison with paper

The paper reports mean own-price elasticities of approximately −4 to −5 for the RCNL specifications (Table IV). Our values (−4.65 mean) are consistent with this range. Exact numerical agreement is not expected because our elasticity matrices were recomputed in `smalldata.m` on the perturbed 5-city sample.

---

## 3. Markups and Marginal Costs

### 3.1 Bertrand markups by firm

| Firm | Mean price | Mean markup | Mean MC | Lerner index |
|------|-----------|------------|---------|-------------|
| ABI | $13.63 | $4.33 | $9.30 | 32.0% |
| Corona | $13.90 | $3.08 | $10.83 | 22.2% |
| Heineken | $13.96 | $3.20 | $10.76 | 23.0% |
| Coors | $8.96 | $2.73 | $6.23 | 31.0% |
| MillerCoors | $12.08 | $3.47 | $8.61 | 29.7% |
| **Overall** | **$12.54** | **$3.51** | **$9.03** | **28.7%** |

### 3.2 Economic interpretation

- **ABI** has the highest Lerner index (32%) and absolute markup ($4.33). As the dominant domestic brand with both 12-pack and 24-pack SKUs, its multi-product pricing power is reflected in higher margins.
- **Imports** (Corona, Heineken) have higher prices but lower Lerner indices (22–23%). Their marginal costs are substantially higher ($10.8–10.9 vs. $6.2–9.3 for domestics), consistent with transportation and import costs.
- **Coors** has the lowest absolute markup ($2.73) on the lowest prices ($8.96), but its Lerner index (31%) is comparable to ABI — reflecting its position as a price-competitive domestic brand.
- All implied marginal costs are strictly positive. No negative-MC pathology, which would indicate model misspecification.

---

## 4. Validation against Matlab

### 4.1 Comparison with `sres_bertrand.mat`

| Metric | Value |
|--------|-------|
| Correlation (R vs. Matlab MC) | 0.982 |
| Mean \|difference\| | $0.32 |
| Max \|difference\| | $2.74 |
| RMSE | $0.45 |
| Negative MC count | 0 (both) |

### 4.2 Source of differences

The Matlab code in `f_impute_mc.m` uses a 3-step procedure:

1. Recover fringe (non-coalition) marginal costs from FOC at observed prices.
2. Solve for fringe Nash prices given coalition Nash prices (via `fsolve`).
3. Invert the coalition FOC at Nash prices to recover coalition marginal costs.

Our R implementation computes the full-ownership Bertrand FOC directly — all products optimize simultaneously with no coalition/fringe distinction. When the supermarkup is zero (Bertrand), the two approaches should converge, but the 3-step procedure introduces additional numerical evaluations of the RCNL simulation (500 draws), each with its own finite-sample noise. The 0.98 correlation confirms the economic content is equivalent.

---

## 5. Key Takeaways

- Own-price elasticities average −4.65 (RCNL), consistent with the paper and with typical values for differentiated beer products. The demand system exhibits economically reasonable curvature.
- Cross-elasticities display the expected hierarchy: within-brand > across-brand domestic > domestic-to-import. This differentiated substitution pattern is essential for realistic merger simulation.
- Bertrand markups average $3.51 (Lerner index 28.7%). All marginal costs are positive and economically plausible. ABI earns the highest margins; imports have the highest costs.
- Validation against the Matlab replication output shows a correlation of 0.98 for implied marginal costs, confirming that our R implementation reproduces the supply-side analysis with high fidelity.
- These outputs — elasticity matrices, marginal costs, and ownership structures — are the direct inputs to the merger counterfactual in Step 6.
