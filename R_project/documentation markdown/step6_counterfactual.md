# Step 6: Merger Counterfactual — Bertrand Simulation
## Miller & Weinberg (2017) — Replication in R

---

## 1. Simulation Setup

### 1.1 Equilibrium concept

The counterfactual simulates a **Bertrand-Nash** equilibrium under post-merger ownership. Each firm simultaneously chooses prices for its products to maximize profit, taking rival prices as given. The first-order condition for product $j$ owned by firm $f$:

$$p_j = mc_j + \underbrace{\left[ -(\mathbf{\Omega}^{\text{post}} \circ \frac{\partial \mathbf{s}}{\partial \mathbf{p}}')^{-1} \mathbf{s} \right]_j}_{\text{markup}_j(\mathbf{p})}$$

where $\mathbf{\Omega}^{\text{post}}_{jk} = 1$ if products $j$ and $k$ are co-owned under the post-merger ownership structure.

### 1.2 Ownership change

| Period | Firm 4 (Coors) | Firm 5 (MillerCoors) | Co-owned pairs (of 28) |
|--------|---------------|---------------------|----------------------|
| Pre-merger | Coors Light 12pk, 24pk | Miller Lite 12pk, 24pk | 3 |
| **Post-merger** | **Absorbed into firm 5** | **Coors Light + Miller Lite (all 4 SKUs)** | **7** |

The merger internalizes pricing externalities between Coors Light and Miller Lite. Pre-merger, a price increase on Coors Light diverts demand partly to Miller Lite (a rival product). Post-merger, this diversion is captured by the same firm, reducing the competitive cost of a price increase.

### 1.3 Demand parameters

RCNL parameters from `small_dresgmm2.mat`, taken as given:
- $\alpha = -0.1087$, $\rho = 0.7779$, $\theta_2 = [0.0009, 0.0125, 0.0045]$
- Shares and derivatives recomputed at each candidate price vector via full RCNL simulation (500 consumer draws).

### 1.4 Marginal costs

Held fixed at pre-merger Bertrand values recovered in Step 5. The assumption is that the merger does not change production costs (the pure unilateral-effects benchmark). The paper considers cost efficiencies separately in counterfactual scenarios 2 and 8.

### 1.5 Solver

`nleqslv` (Broyden method), tolerances: $\text{xtol} = 10^{-12}$, $\text{ftol} = 10^{-12}$, $\text{maxit} = 500$. Starting values: pre-merger observed prices.

---

## 2. Results

### 2.1 Price changes by firm

| Firm | $\bar{p}^{\text{pre}}$ | $\bar{p}^{\text{post}}$ | $\Delta p$ ($) | $\Delta p$ (%) |
|------|------------------------|-------------------------|---------------|---------------|
| **Coors** | $8.96 | $9.32 | **+0.365** | **+4.26%** |
| **MillerCoors** | $12.08 | $12.22 | **+0.142** | **+1.13%** |
| ABI | $13.63 | $13.70 | +0.065 | +0.49% |
| Corona | $13.90 | $13.89 | −0.013 | −0.09% |
| Heineken | $13.96 | $13.95 | −0.017 | −0.12% |
| **Overall** | **$12.54** | **$12.65** | **+0.111** | **+1.05%** |

### 2.2 Price changes by product

| Product | $\Delta p$ ($) | $\Delta p$ (%) |
|---------|---------------|---------------|
| Coors Light 24pk | +0.202 | +2.59% |
| Miller Lite 24pk | +0.270 | +2.30% |
| Coors Light 12pk | +0.162 | +1.67% |
| Miller Lite 12pk | +0.157 | +1.08% |
| Bud Light 24pk | +0.070 | +0.56% |
| Bud Light 12pk | +0.061 | +0.42% |
| Corona 12pk | −0.013 | −0.09% |
| Heineken 12pk | −0.017 | −0.12% |

Larger packs show bigger dollar-value increases, consistent with higher cross-elasticities within the merged entity's portfolio.

### 2.3 Merging firms vs. rivals

| Group | $\overline{\Delta p}$ ($) | $\overline{\Delta p}$ (%) |
|-------|--------------------------|--------------------------|
| Merging parties (Coors + MC) | +0.198 | +1.91% |
| Rival: ABI | +0.065 | +0.49% |
| Rival: imports | −0.015 | −0.10% |

### 2.4 Consumer surplus

| Metric | Value |
|--------|-------|
| Mean CS (pre-merger) | 3.203 |
| Mean CS (post-merger) | 3.180 |
| $\Delta$CS | **−0.024** |
| $\Delta$CS (%) | **−0.85%** |

Consumer surplus declines in the majority of markets, reflecting the higher post-merger prices. The magnitude is modest because the merging parties account for approximately 20% of inside-good market share, and the price increases are moderate (1–4%).

---

## 3. Validation Checks

| Check | Result |
|-------|--------|
| Markets converged | **80/80** (100%) |
| Max FOC residual norm | $9.97 \times 10^{-13}$ |
| Negative post-merger prices | **0** |
| Negative post-merger shares | **0** |
| Post-merger markups positive | **Yes** (all products) |

The solver achieves machine-precision convergence in all markets. The equilibrium is well-behaved: no boundary solutions, no corner cases.

---

## 4. Comparison with Paper

### 4.1 Paper's counterfactual structure

The paper's Tables 6–8 report merger counterfactuals under the **Price Leadership Equilibrium** (PLE) model, not standard Bertrand. The PLE model allows ABI to set a supermarkup above Bertrand-Nash prices, which all coalition members (ABI + MillerCoors) follow. The Bertrand simulation (our exercise) corresponds to the **unilateral effects component** only.

The paper's `cfmerger.txt` contains results for 6 scenarios across PLE and Bertrand. The first row of each scenario block represents the Bertrand price change; the second row adds the PLE coordination component.

### 4.2 Magnitude comparison

| Metric | Our simulation | Paper (Bertrand row, Scenario 1) |
|--------|---------------|----------------------------------|
| Merging firms $\Delta p$ | +$0.20 (1.9%) | +$0.64 (Scenario 1, col 1) |
| ABI $\Delta p$ | +$0.065 (0.5%) | +$0.50 (Scenario 1, col 1) |
| ABI/merging ratio | 0.33 | ~0.78 |

Our price effects are smaller in absolute terms. This is expected because:
- The small sample (5 cities, 8 products) has different demand curvature than the full data (37 regions, 39 products).
- The perturbed prices and quantities alter the equilibrium derivatives.
- The paper's Bertrand results condition on the PLE-estimated marginal costs (with supermarkup recovery at $\delta_f = 0.26$), while ours condition on standard Bertrand MC.

The qualitative pattern is consistent: merging firms raise prices more than rivals, and ABI responds modestly.

---

## 5. Economic Interpretation

### 5.1 Unilateral effects

The merger creates a 4-product firm (Coors Light 12pk, 24pk + Miller Lite 12pk, 24pk) that internalizes the pricing externalities between these brands. Pre-merger, a price increase on Coors Light 12pk diverts some consumers to Miller Lite 12pk — a gain captured by a rival. Post-merger, this diversion stays within the firm, reducing the competitive cost of raising prices. The result: merging firms raise prices by +1.9% on average.

### 5.2 Rival responses

- **ABI** (+0.5%): prices are strategic complements in differentiated-products Bertrand. When the merged entity raises prices, ABI faces less competition and optimally raises its own prices — but by much less than the merging parties.
- **Imports** (−0.1%): Corona and Heineken are weak substitutes for domestic beer (cross-elasticity ~0.28 vs. 0.31 for domestic-to-domestic). The merger barely affects their pricing incentives. The slight price *decrease* reflects a small share reallocation toward imports as domestic prices rise.

### 5.3 The coordination gap

The central finding of MW(2017) is that Bertrand unilateral effects **underpredict** the observed price changes. In the data, ABI's prices rose by approximately the same amount as Miller/Coors after the joint venture — a pattern inconsistent with the 0.33 ABI/merging ratio predicted by Bertrand. The paper's PLE model, in which ABI leads a coalition that coordinates on a common supermarkup, closes this gap. Our Bertrand simulation quantifies the unilateral-effects baseline against which the coordination hypothesis is measured.

---

## 6. Model Limitations

- **This is a Bertrand counterfactual, not the paper's full structural model.** The paper estimates a Price Leadership Equilibrium (PLE) with incentive compatibility constraints (ICC), in which the MillerCoors joint venture facilitates coordination among ABI, Miller, and Coors. Our simulation captures only the unilateral effects of the merger — the component due to ownership-matrix consolidation.
- **Marginal costs are held fixed.** The paper considers efficiency scenarios (Scenarios 2 and 8 in `cfmerger.txt`) in which the merger reduces MillerCoors production costs. Our baseline assumes no cost efficiencies.
- **The small perturbed sample limits quantitative precision.** Price effects are directionally correct but smaller than the paper's, reflecting the reduced cross-sectional variation and the data perturbation. Results should be interpreted qualitatively.
- **The outside good is passive.** The model assumes the outside option does not respond to inside-good price changes. In reality, consumers may substitute to craft beer, wine, or spirits — categories not in the model.

---

## 7. Key Takeaways

- The MillerCoors merger raises merging parties' prices by approximately **+2%** and rival ABI's prices by **+0.5%** under standard Bertrand-Nash competition. Consumer surplus declines by about **0.9%**.
- Import brands (Corona, Heineken) are largely unaffected, consistent with weak cross-brand substitution between domestic and imported beer.
- The ratio of ABI's price increase to the merging parties' increase is approximately **0.33** under Bertrand — substantially below the near-1.0 ratio observed in the actual data. This gap is the empirical basis for the paper's coordination hypothesis.
- The Bertrand counterfactual serves as a **lower bound** on the price effects of the merger. The paper's PLE model, which accounts for coordination facilitated by the joint venture, predicts larger and more symmetric price increases across all coalition members.
- These results are qualitatively consistent with the paper's findings and demonstrate the standard empirical IO toolkit for merger simulation: demand estimation → marginal cost recovery → ownership change → equilibrium re-solve → welfare computation.
