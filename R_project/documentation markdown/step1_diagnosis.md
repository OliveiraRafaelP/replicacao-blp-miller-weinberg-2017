# Step 1: Replication Package Diagnosis
## Miller & Weinberg (2017) — Matlab-to-R Translation Plan

---

## A. Directory Structure Map

```
replication/
├── code/
│   ├── main_starter.m              ← MASTER SCRIPT (orchestrates full pipeline)
│   ├── MTIMESX/                    ← 3D matrix multiply library (C MEX)
│   ├── functions/                  ← Core estimation code
│   │   ├── main_data.m             ← Data construction from raw .mat
│   │   ├── main_spec.m             ← Model specification & options
│   │   ├── smalldata.m             ← Creates small replication sample
│   │   ├── f_daugment.m            ← Augments demand results (delta, mu, ai)
│   │   ├── impute_bertrand.m       ← Bertrand baseline (sm=0)
│   │   ├── main_supply_bind.m      ← Supply: binding ICC (with rebalancing)
│   │   ├── main_supply_nonbind.m   ← Supply: non-binding ICC
│   │   ├── main_supply_bind_nopool.m ← Binding ICC, no pooling
│   │   ├── main_supply_bind_x.m    ← Binding ICC, two supermarkups
│   │   ├── combine_imputed.m       ← Combines year-specific results
│   │   ├── f_impute_mc.m           ← Core MC imputation (3-step FOC inversion)
│   │   ├── f_mu.m                  ← Consumer heterogeneity (mu_ij = x2 * d_i * theta2w')
│   │   ├── f_ownMat.m              ← Ownership matrix construction
│   │   ├── f_pi_m.m                ← PLE equilibrium prices/shares/profits
│   │   ├── f_loss_bind.m           ← ICC constraint evaluation (binding)
│   │   ├── f_loss_nonbind.m        ← Loss function (non-binding)
│   │   ├── f_rebalance.m           ← Iterative regional deviation adjustment
│   │   ├── f_inclusive.m           ← Inclusive value for CS computation
│   │   ├── rcnl_meanval.m          ← BLP contraction mapping (calls C MEX)
│   │   ├── rcnl_indsh.m            ← Individual choice probabilities (RCNL)
│   │   ├── rcnl_der1.m             ← Price derivatives (ds/dp)
│   │   ├── cf_foc_partial.m        ← FOC evaluation (partial equilibrium)
│   │   ├── results_baseanalysis.m  ← Tables 2, 4, 5; Figures 3, 4, 5
│   │   ├── results_costregs.m      ← Table 3 (cost regressions)
│   │   ├── results_cfanalysis.m    ← Tables 6, 7, 8; Figures 6, G3
│   │   └── results_cfmergers.m     ← Merger counterfactual engine
│   └── idcheck/                    ← Monte Carlo identification check
├── data/
│   ├── raw/
│   │   ├── small_scanner.mat       ← 1120x9: prices, quantities, IDs
│   │   ├── small_demosE.mat        ← 35x1500: income draws
│   │   └── small_dresgmm2.mat      ← Demand estimates from MW(2017)
│   └── analysis/                   ← Intermediate results
│       ├── daugfile.mat            ← Augmented demand (delta, mu, pcoefi)
│       ├── sres_bertrand.mat       ← Bertrand MC baseline
│       ├── gammapar.mat            ← Cost regression coefficients
│       ├── df_250/                 ← Results for delta_f = 0.25
│       └── ...                     ← Other discount factor folders
└── results/                        ← Tables (.txt) and Figures (.pdf)
```

---

## B. Data Architecture

### B.1 Raw Scanner Data: `small_scanner.mat` (1120 x 9)

The composite ID in column 1 encodes five identifiers via positional digits:

```
ID = firmid * 1e10 + brndid * 1e8 + sizeid * 1e6 + cityid * 1e4 + yearid * 1e2 + montid
```

| Column | Variable   | Description                         | Values                          |
|--------|-----------|-------------------------------------|---------------------------------|
| 1      | id2       | Composite ID (firm+brand+size+city+year+month) | ~1e10–5e10       |
| 2      | p_jt      | Price ($ per unit)                  | 6.69–17.63                      |
| 3      | sizeid    | Pack size class                     | 2 (12-pack), 3 (24/30-pack)    |
| 4      | q_jt      | Unit sales (144oz equiv.)           | 472–75,606                      |
| 5      | miles     | Miles to brewery (raw)              | 0–2.35                          |
| 6      | calor_raw | Calories (raw, /100 then demeaned)  | 102–166                         |
| 7      | msize     | Market size (potential consumers)   | 75,347–613,831                  |
| 8      | dist      | Distance to ABI brewery             | 0–10.31                         |
| 9      | distbutfor| Distance (but-for counterfactual)   | 0–10.31                         |

### B.2 Decoded IDs (from `main_data.m`)

| firmid | Firm          | brndid | Brand               |
|--------|---------------|--------|---------------------|
| 1      | ABI           | 1      | Bud Light (12-pk)   |
| 2      | Import A      | 4      | Coors Light (12-pk) |
| 3      | Import B      | 5      | Corona (12-pk only) |
| 4      | Miller        | 7      | Heineken (12-pk only)|
| 5      | MillerCoors   | 13     | Miller Lite (12-pk) |

- **5 cities** (cityid 1–5, from 37 original IRI regions)
- **8 products** (5 brands x 2 sizes, minus Corona 24pk and Heineken 24pk)
- **~140 market-time observations** (city x quarter, excluding 1 year post-merger)
- **8 products per market** (consistent J=8 across all markets)

### B.3 Coalition Structure

```
coalid = (firmid == 1 | firmid == 5 | firmid == 4)  →  ABI + Miller + MillerCoors
leadid = (firmid == 1)                               →  ABI is price leader
```

**Critical insight**: The paper models ABI as the price leader of a coalition that includes Miller and Coors/MillerCoors. This is NOT standard Bertrand — it's Price Leadership Equilibrium (PLE).

### B.4 Demand Estimates: `small_dresgmm2.mat`

| Variable    | Shape      | Content                                       |
|-------------|-----------|-----------------------------------------------|
| theta1_2    | (63, 1)   | Linear params: alpha (price coef) + 62 FEs    |
| theta2_2    | (3, 1)    | RC params: [0.0009, 0.0125, 0.0045]           |
| rho_2       | (1, 1)    | Nesting param = 0.778                          |
| elasMat_2   | (8,8,140) | Elasticity matrices per market                 |
| derMat_2    | (8,8,140) | Derivative matrices (ds/dp) per market         |

**theta2w construction** (RCNL2 specification):
```
theti = [1; 2; 3]    →  rows: price, constant, calories
thetj = [2; 2; 2]    →  all in column 2 (income interaction)
theta2w = sparse(theti, thetj, theta2_2)  →  3x2 matrix
```
This means: price sensitivity, base utility, and calorie preference all vary with income.

### B.5 Demographics: `small_demosE.mat` (35 x 1500)

- 35 rows = 7 years x 5 cities (year-city combinations)
- 1500 columns, but only first 500 used (ns=500)
- Demeaned household income draws from PUMS/ACS
- Expanded to observation level via `yearcityid` mapping

---

## C. Full Estimation Pipeline

### Phase 0: Data Augmentation (`f_daugment.m`) — RUN ONCE

```
Raw data → main_data() → vars, ids
                ↓
Demand params (theta1, theta2, rho) from small_dresgmm2.mat
                ↓
Contraction mapping: rcnl_meanval() → delta (mean utility)
                ↓
mu = x2 * (d_i * theta2w')       → consumer heterogeneity
ai = d_i * theta2w(1,:)          → individual price coefficient
pcoefi = alpha + ai              → full price coefficient
deltanp = delta - alpha * p_jt   → non-price mean utility
                ↓
Save → daugfile.mat (used by ALL subsequent steps)
```

**Economic logic**: Takes demand-side estimates as GIVEN from MW(2017). Does not re-estimate demand. The supply-side analysis conditions on these demand parameters.

### Phase 1: Bertrand Baseline (`impute_bertrand.m`)

```
For each market (cdid):
    f_impute_mc(sm=0, ...) → mc_bertrand
```
With sm=0 (no supermarkup), observed prices ARE Nash prices. This gives baseline MC under standard Bertrand competition.

### Phase 2: Supply-Side with Binding ICCs (`main_supply_bind.m`)

This is the CORE of the paper. For each discount factor δ ∈ {0.20, 0.25, 0.26, 0.30, 0.35, 0.40}:

```
For each fiscal year (2006, 2007, 2010, 2011):
    Iteration 0: fminsearch → base supermarkup (sm_base)
    Iterations 1-8: f_rebalance → regional deviations (sm_devs)
        Each iteration tightens ICC tolerance: 0.2 → 0.001
        
    Output: mc, pnash, snash, smbase, smdevs, smfinal
```

**ICC constraint** (Incentive Compatibility):
```
V^PLE ≥ V^Nash + δ/(1-δ) * (V^Deviation - V^Nash)
```
Where V^PLE = coalition profit under price leadership, V^Nash = Nash profit, V^Deviation = one-shot deviation profit.

### Phase 3: MC Imputation (`f_impute_mc.m`) — 3-Step Algorithm

```
Step 1: Fringe MC
    - Fringe firms play Nash-Bertrand given observed prices
    - mc_fringe = p_observed - markup (from FOC inversion)

Step 2: Fringe Nash Prices
    - Coalition Nash prices: p_nash = p_observed - sm
    - Solve fsolve() for fringe FOC given fixed coalition Nash prices
    
Step 3: Coalition MC
    - Evaluate derivatives at Nash prices
    - Invert coalition FOC: mc_coal = p_nash - (Owner .* der')^{-1} * s_nash
```

### Phase 4: Results Generation

| Function               | Outputs                    | Paper Tables/Figures |
|------------------------|----------------------------|---------------------|
| results_costregs       | Cost regressions, df selection | Table 3, G1, G2    |
| results_baseanalysis   | Summary stats, markups, profits | Tables 2, 4, 5; Figs 3-5 |
| results_cfmergers      | Merger simulations         | (intermediate)      |
| results_cfanalysis     | CF price/welfare effects   | Tables 6, 7, 8; Fig 6 |

---

## D. Validation Targets (from pre-computed results)

### Table 3, Column 1 (spec_df_coef_mod1.txt): Cost Regression Coefficients
Row 1 (mpost — Miller post-merger): ranges from 1.007 (Bertrand) to -1.046 (df=0.40)
Row 2 (cpost — Coors post-merger): ranges from 0.811 to -1.321
Row 3 (dist — distance to ABI): ranges from 0.327 to -1.849
Row 4 (apost — ABI post-merger): ranges from 0.057 to -0.018

**Key finding**: df=0.26 (col 4) gives "sensible" cost coefficients: mpost=-0.015, cpost=-0.252

### Table 2 (sumstats_pq.txt): Summary Statistics
ABI: share=15.1%, price=$14.91 (12pk); share=31.2%, price=$12.42 (24pk)

### Table 5 (mean_markups.txt): Mean Supermarkups
ABI: 5.33 (2007), 6.20 (2010) — substantial price leadership premium

### Table 4 (eqeffects.txt): Equilibrium Effects
Profit gains from PLE: 22.34–31.37% (row 1)
Consumer surplus loss: -1.34 to -1.66% (row 2)

---

## E. Key Economic Model (NOT Standard BLP)

This paper is **not** a standard BLP exercise. Critical distinctions:

### E.1 Price Leadership Equilibrium (PLE)
- ABI is the price **leader** of a tacit coalition (ABI + Miller + Coors)
- Leader sets a **supermarkup** sm above Nash-Bertrand prices
- All coalition members' prices shift by sm (or sm varies by size class)
- Fringe firms (imports) best-respond to coalition prices

### E.2 Incentive Compatibility Constraints (ICC)
- Coalition is sustained by repeated game logic
- ICC ensures no firm wants to deviate from PLE to grab one-period Nash profit
- Discount factor δ governs how much firms value future cooperation
- The paper finds δ ≈ 0.26 produces economically sensible cost estimates

### E.3 Supermarkup Structure
- `bysize=0`: Single supermarkup for all products
- `bysize=2`: Two supermarkups (12-pack vs. 24/30-pack)
- Regional deviations allow ICC to bind differentially across cities

### E.4 What the Problem Set Asks vs. What the Paper Does
| Problem Set Question | Paper's Actual Model |
|---------------------|---------------------|
| "Linear 2SLS" (Q3.1) | Can approximate with logit IV regression |
| "Full BLP random coefs" (Q3.2) | Paper uses RCNL (nested + random coefs) |
| "Merger simulation Bertrand" (Q4) | Paper does PLE + ICC, not pure Bertrand |

The problem set simplifies by asking for standard BLP/Bertrand merger simulation. The replication package implements the full PLE model. We should do BOTH.

---

## F. Proposed R Project Architecture

```
R_project/
├── 00_setup.R              ← Packages, paths, helper functions
├── 01_load_data.R          ← Load .mat files, decode IDs, build tibble
├── 02_clean_shares.R       ← Shares, outside option, trimming, validation
├── 03_instruments.R        ← BLP instruments, cost shifters, market FEs
├── 04_demand_linear.R      ← 2SLS logit demand (ivreg)
├── 05_demand_blp.R         ← Full BLP via pyblp (reticulate) or manual
├── 06_elasticities.R       ← Own/cross elasticities, comparison w/ paper
├── 07_markups_mc.R         ← FOC inversion, MC imputation
├── 08_merger_sim.R         ← Bertrand merger counterfactual
├── 09_welfare.R            ← Consumer surplus, total surplus
├── 10_bonus_efficiency.R   ← 5% MC reduction scenario
├── helpers/
│   ├── rcnl_shares.R       ← RCNL individual share computation
│   ├── rcnl_derivatives.R  ← Price derivatives (ds/dp)
│   ├── ownership_matrix.R  ← f_ownMat equivalent
│   ├── contraction_map.R   ← BLP contraction mapping
│   └── foc_inversion.R     ← FOC partial equilibrium solver
├── output/
│   ├── tables/             ← Replicated tables (CSV/LaTeX)
│   ├── figures/            ← Plots (PDF/PNG)
│   └── validation/         ← Cross-checks against Matlab results
├── data/                   ← Processed CSV files
├── problemset.Rmd          ← Final RMarkdown deliverable
└── step1_diagnosis.md      ← This document
```

---

## G. Technical Risks: Matlab → R Translation

### RISK 1: C MEX Contraction Mapping [HIGH]
**Problem**: `rcnl_meanval.m` calls `contrMap_rcnl.mexw64` — compiled C code for the BLP contraction mapping with RCNL nesting.
**Mitigation**: 
- Write pure R contraction mapping (slow but correct)
- Alternatively use `pyblp` via reticulate (has optimized contraction mapping)
- If speed matters, write Rcpp version
**Impact**: Correctness is fine; speed may be 10-100x slower in pure R

### RISK 2: Demand Parameters Are GIVEN, Not Re-estimated [MEDIUM]
**Problem**: The replication package does NOT re-estimate demand. It loads `small_dresgmm2.mat` with pre-estimated theta1, theta2, rho. The problem set asks students to estimate BLP.
**Mitigation**:
- For the problem set: use pyblp to estimate demand from the constructed dataset
- For replication fidelity: also load the given parameters and verify supply-side results
**Impact**: Two parallel tracks needed

### RISK 3: 3D Matrix Operations (mtimesx) [MEDIUM]
**Problem**: `rcnl_der1.m` uses mtimesx for 3D matrix multiplications. R has no native equivalent.
**Mitigation**: 
- Use `array()` with `apply()` or explicit loops over the 3rd dimension
- For 8x8x140, the overhead is negligible
**Impact**: Low — small data makes this manageable

### RISK 4: fsolve → R Equivalent [MEDIUM]
**Problem**: Matlab's `fsolve` (Levenberg-Marquardt) used in `cf_foc_partial` for equilibrium solving.
**Mitigation**: 
- R package `nleqslv` provides equivalent nonlinear equation solvers
- `rootSolve::multiroot` as backup
- Both support Levenberg-Marquardt and Newton methods
**Impact**: Should translate cleanly; edge cases in convergence may differ

### RISK 5: Perturbed/Masked Data [LOW but IMPORTANT]
**Problem**: `smalldata.m` reveals the replication data is intentionally perturbed:
- Prices: `price + 0.2*(rand-0.5)`; ABI/Miller prices multiplied by 1.5
- Quantities: `quant + (rand-0.5)`; ABI/Miller quantities multiplied by 1.5
- Only 5 of 37 regions, 8 of 39 products
**Implication**: Results will NOT exactly match the published paper. They should be directionally consistent. Validation targets are the pre-computed `.txt` files in `results/`, not the paper's tables.

### RISK 6: PLE vs. Standard BLP for Problem Set [LOW]
**Problem**: The problem set pseudocode assumes standard BLP + Bertrand merger. The actual replication package implements PLE + ICC — a substantially more complex model.
**Mitigation**:
- Scripts 04-08: implement the standard BLP pipeline the problem set asks for
- Supplementary analysis: replicate the full PLE model for bonus/discussion
**Impact**: The problem set is intentionally simplified; we deliver both

### RISK 7: Fiscal Year / Time Period Alignment [LOW]
**Problem**: Matlab code uses fiscal years (Oct-Sep) with RCNL2 specification shifting by `montid >= 4`. Misalignment could cause wrong sample selection.
**Mitigation**: Replicate exact `obsin` logic from `main_data.m` lines 44-63.
**Impact**: Affects which observations enter estimation; must match exactly

---

## H. Confirmed vs. Uncertain Elements

### CONFIRMED (verified in code + data)
- [x] 1120 obs, 9 columns in scanner data
- [x] 8 products per market (J=8), 5 cities, ~140 market-time obs
- [x] ID decoding formula: firm/brand/size/city/year/month from composite ID
- [x] Coalition = {ABI(1), Miller(4), MillerCoors(5)}, Leader = ABI(1)
- [x] Demand params: alpha = theta1_2(1), rho = 0.778, theta2 has 3 RC params
- [x] x2 = [price, constant, calories] for RCNL2
- [x] ns = 500 simulated consumers
- [x] Outside share = 1 - sum(inside shares per market)
- [x] Cost shifters: mpost, cpost, dist, product FE, city FE, date FE
- [x] Optimal discount factor: df ≈ 0.26

### UNCERTAIN (need verification in Step 2)
- [ ] Exact firm-to-brand mapping (need to decode all unique id2 values)
- [ ] Whether `calor` column 6 values (102-166) represent actual calories or index
- [ ] Whether shares sum correctly within markets after filtering
- [ ] Exact number of observations after obsin filtering for RCNL2
- [ ] Whether pyblp can handle the RCNL nesting parameter directly

---

## I. Immediate Next Steps (Step 2)

1. Load `small_scanner.mat` in R via `R.matlab::readMat()`
2. Decode composite IDs; build product-market tibble
3. Validate: J=8 per market, unique firms, share sums
4. Compute inside/outside shares
5. Plot: shares vs prices by firm (faceted by market sample)
6. Cross-check against `sumstats_pq.txt` targets
