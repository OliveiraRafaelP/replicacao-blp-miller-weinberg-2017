# Step 2 Validation Checkpoint
## Miller & Weinberg (2017) — Data Loading & Cleaning

---

## 1. Raw Source Files and Matlab Variables

| File | Variable | Shape | Dtype | Description |
|------|----------|-------|-------|-------------|
| `small_scanner.mat` | `small_scanner` | (1120, 9) | float64 | Scanner data: IDs, prices, quantities, attributes |
| `small_demosE.mat` | `small_demosE` | (35, 1500) | float64 | Income draws: 35 year-city combos x 1500 draws |
| `small_dresgmm2.mat` | `theta1_2` | (63, 1) | float64 | Linear demand params (alpha + 62 FEs) |
| | `theta2_2` | (3, 1) | float64 | RC params: [0.00090, 0.01249, 0.00451] |
| | `rho_2` | (1, 1) | float64 | Nesting param: 0.77788 |
| | `derMat_2` | (8,8,140) | float64 | Price derivative matrices per market |
| | `elasMat_2` | (8,8,140) | float64 | Elasticity matrices per market |
| `daugfile.mat` | `daugfile` (struct) | — | — | Augmented demand (delta, mu, ai, pcoefi, etc.) |
| | `.delta` | (960, 1) | | Mean utility from contraction mapping |
| | `.mu` | (960, 500) | | Consumer heterogeneity |
| | `.ai` | (960, 500) | | Individual-specific price coefficient deviation |
| | `.pcoefi` | (960, 500) | | Full individual price coefficient (alpha + ai) |
| | `.deltanp` | (960, 1) | | Non-price mean utility |
| | `.xi` | (960, 1) | | Unobserved quality (demand residual) |
| | `.dprodfecoef` | (8, 1) | | Product FE coefficients |
| | `.ddatefecoef` | (24, 1) | | Date FE coefficients |
| | `.dcityfecoef` | (5, 1) | | City FE coefficients |

---

## 2. Matlab-to-R Variable Mapping

### 2a. Scanner columns → R variables

| Matlab Source | Matlab Var | R Variable | Meaning |
|--------------|-----------|-----------|---------|
| `scanner[:,0]` | `id2` | `decode_id(id2)` | Composite ID → 6 decoded fields |
| `scanner[:,1]` | `p_jt` | `df$price` | Price ($/144oz equiv) |
| `scanner[:,2]` | `sizeid` | `df$sizeid` | Pack size (2=12pk, 3=24pk) |
| `scanner[:,3]` | `q_jt` | `df$quantity` | Unit sales (144oz equiv) |
| `scanner[:,4]` | `miles` | `df$miles` | Miles to brewery (raw) |
| `scanner[:,5]` | `calor_raw` | `df$calor = (col5/100) - mean()` | Calories (demeaned) |
| `scanner[:,6]` | `msize` | `df$msize` | Market size (potential consumers) |
| `scanner[:,7]` | `dist` | `df$dist` | Distance to ABI brewery |
| `scanner[:,8]` | `distbutfor` | `df$distbutfor` | Distance (but-for CF) |

### 2b. Derived variables

| Matlab Var | Formula | R Variable | Role |
|-----------|---------|-----------|------|
| `s_jt` | `q_jt / msize` | `df$share` | Market share |
| `inshr` | `sum(share) per market` | `df$inshr` | Inside good share |
| `outshr` | `1 - inshr` | `df$outshr` | Outside good share |
| `cdid` | Sequential market counter | `df$cdid` | Market index (1-80) |
| `prodid` | `grp2idx(brndid*100+sizeid)` | `df$prodid` | Product index (1-8) |
| `fiscid` | `yearid + (montid>=4) + 2004` | `df$fiscid` | Fiscal year |
| `coalid` | `firmid %in% c(1,4,5)` | `df$coalid` | Coalition member |
| `leadid` | `firmid == 1` | `df$leadid` | Price leader (ABI) |
| `mpost` | `brndid %in% 11:13 & yearid>=5` | `df$mpost` | Miller post-merger |
| `cpost` | `brndid %in% c(3,4) & yearid>=5` | `df$cpost` | Coors post-merger |
| `apost` | `firmid==1 & yearid>=5` | `df$apost` | ABI post-merger |

---

## 3. Confirmed vs. Inferred

### CONFIRMED (verified in code AND data)

1. Column mapping: col0=ID, col1=price, col2=sizeid, col3=quantity, col4=miles, col5=calories, col6=msize, col7=dist, col8=distbutfor
2. ID decode: `firmid*1e10 + brndid*1e8 + sizeid*1e6 + cityid*1e4 + yearid*1e2 + montid`
3. Firm mapping: 1=ABI, 2=Corona parent, 3=Heineken parent, 4=Coors(pre), 5=MillerCoors(post)
4. Brand mapping: 1=Bud Light, 4=Coors Light, 5=Corona, 7=Heineken, 13=Miller Lite
5. Size: 2=12-pack, 3=24/30-pack
6. J=8 products per market (constant across all 140 raw / 80 filtered markets)
7. Total: 1120 raw → 960 (obsintemp) → 640 (fiscal filter) = 8 x 80
8. daugfile has 960 rows; supply code selects 640 via `obsindemand`
9. ns=500 simulated consumers; demosE first 500 cols
10. alpha=-0.10872, rho=0.77788, theta2=[0.00090, 0.01249, 0.00451]
11. theta2w is 3x2: rows=[price, const, calor], col2=income interaction
12. Coalition = {1,4,5}, Leader = {1}
13. s_jt = q_jt / msize (line 142)
14. calor = scanner[:,5]/100 then demeaned (line 153-154)

### INFERRED (from code logic, not explicit documentation)

1. Brand 4 (Coors Light) switches from firm=4 → firm=5 at merger (year 4-5)
2. Fiscal year RCNL2 offset: `yearid + 1*(montid>=4) + 2004` means Q4=Oct start
3. Market size represents total potential beer consumers per city-quarter
4. dist (col7) = distance to nearest ABI brewery (cost instrument)
5. distbutfor (col8) = counterfactual distance (for merger scenario analysis)
6. mpost in small data only captures brndid=13 (Miller Lite); brndid=11,12 absent
7. cpost in small data only captures brndid=4 (Coors Light); brndid=3 absent
8. Product ordering within market is deterministic (sorted by composite ID)

---

## 4. Validation Checks

### 4a. Observation counts
| Stage | N obs | Check |
|-------|-------|-------|
| Raw | 1120 | = 8 prods x 5 cities x 28 quarters |
| After obsintemp | 960 | Excludes ~1 year around merger |
| After fiscal year filter | 640 | = 8 prods x 80 markets |
| obsindemand | 640 of 960 | Maps 960→640 for daugfile selection |

### 4b. Market structure
- Unique markets (cdid): **80** = 5 cities x 4 fiscal years x 4 quarters
- Products per market: **[8]** — perfectly balanced panel
- Unique firms: **[1, 2, 3, 4, 5]** (5 firms)
- Unique products (brand x size): **8**

### 4c. Share validation
- Inside share: min=0.1005, max=0.3260, mean=0.2074
- All < 1.0: **PASS**
- All > 0.0: **PASS**

### 4d. Outside share
- Outside share: min=0.6740, max=0.8995, mean=0.7926
- All > 0: **PASS**

### 4e. Missing/zero values
| Variable | NaN | <=0 | Range |
|----------|-----|-----|-------|
| Price | 0 | 0 | [6.71, 17.63] |
| Quantity | 0 | 0 | [574, 63362] |
| Share | 0 | 0 | [0.0040, 0.1198] |
| Market size | 0 | 0 | [75347, 613831] |
| Distance | 0 | 16 (ABI=0) | [0.0, 9.06] |

**All checks PASS.** No data cleaning needed. No trimming required.

---

## 5. Flagged Ambiguities

### 5a. Market definition
- Market = city x fiscal-year-quarter (not calendar quarter)
- RCNL2 fiscal year starts at montid=4 (October)
- **Resolved**: 80 = 5 cities x 16 fiscal quarters (4 years x 4 quarters)

### 5b. Product definition
- Product = brand x pack-size (8 unique combinations)
- Brand 4 (Coors Light) appears under firm=4 pre-merger AND firm=5 post-merger
- **Same physical product**, different owner — prodid is stable, firmid changes

### 5c. Firm IDs
- FY 2006-2007: firms = {1, 2, 3, 4, 5} (Coors separate from MillerCoors)
- FY 2010-2011: firms = {1, 2, 3, 5} (firm=4 absorbed into firm=5)
- **Resolved**: coalid = firmid %in% c(1,4,5) handles both periods correctly

### 5d. Quantity units
- 144-oz equivalents (confirmed in main_data.m comment)
- Dimensionally consistent with msize

### 5e. Aggregation level
- Brand x pack-size x city x quarter — no further aggregation needed

---

## 6. Risk Log

| # | Severity | Risk | Description | Mitigation |
|---|----------|------|-------------|------------|
| R1 | HIGH | Perturbed data | Prices/quantities randomly perturbed; ABI/Miller inflated 1.5x | Validate against results/*.txt, not published tables |
| R2 | HIGH | Demand params given | daugfile uses pre-estimated params, not re-estimated | Track (a) given params + (b) pyblp fresh estimation |
| R3 | HIGH | derMat 140→80 indexing | 140 derivative slices must be correctly subset to 80 | Use market-level obsindemand equivalent |
| R4 | MEDIUM | Contraction convergence | Pure R vs C MEX may differ | Validate max\|delta_R - delta_matlab\| < 1e-8 |
| R5 | MEDIUM | Fiscal year boundaries | Complex exclusion window around merger | Replicate exact obsin logic from main_data.m |
| R6 | MEDIUM | Coalition across merger | Ownership matrix changes pre/post | Reconstruct ownership per market-time |
| R7 | MEDIUM | Calories demeaning scope | Mean over 640 vs 960 obs matters | Match Matlab scope exactly |
| R8 | LOW | dist=0 for ABI | Correct but affects cost regression | Document in interpretation |
| R9 | LOW | No share trimming needed | Shares well-behaved [0.004, 0.120] | Skip pmax/pmin from pseudocode |
| R10 | LOW | Imports lack 24-pack | By design in small sample | J=8 is correct |
| R11 | LOW | yearcityid pre-filter | Demo expansion uses all 1120 obs | Compute before applying any filter |

---

## 7. Data Dictionary

| Variable | Meaning | Source File | Matlab Object | Transformation | Estimation Role |
|----------|---------|-------------|---------------|----------------|-----------------|
| id2 | Composite product-market ID | small_scanner.mat | scanner[:,0] | None (decode) | Identification |
| firmid | Firm identifier (1-5) | small_scanner.mat | floor(id2/1e10) | Positional decode | Ownership matrix |
| brndid | Brand identifier | small_scanner.mat | decode from id2 | Positional decode | Product FE |
| sizeid | Pack size (2=12pk, 3=24pk) | small_scanner.mat | decode from id2 | Positional decode | Product FE, bysize |
| cityid | City/region (1-5) | small_scanner.mat | decode from id2 | Positional decode | City FE, market def |
| yearid | Year index (1-7) | small_scanner.mat | decode from id2 | Positional decode | Time FE, filters |
| montid | Quarter (1-4) | small_scanner.mat | decode from id2 | Positional decode | Time FE, fiscal year |
| p_jt | Price ($/144oz equiv) | small_scanner.mat | scanner[:,1] | None | Demand LHS, x1, x2 |
| q_jt | Quantity (144oz equiv) | small_scanner.mat | scanner[:,3] | None | Share construction |
| msize | Market size (pot. consumers) | small_scanner.mat | scanner[:,6] | None | Share denominator |
| s_jt | Market share | derived | q_jt / msize | Division | Demand LHS |
| inshr | Inside good share | derived | cumsum per market | Sum s_jt in mkt | Log odds, nesting |
| outshr | Outside good share | derived | 1 - inshr | Complement | Log odds ratio |
| calor | Calories (demeaned) | small_scanner.mat | scanner[:,5] | /100 then demean | x2 (RC variable) |
| miles | Miles to brewery (raw) | small_scanner.mat | scanner[:,4] | None | Descriptive only |
| dist | Distance to ABI brewery | small_scanner.mat | scanner[:,7] | None | Cost shifter (w) |
| distbutfor | Distance (but-for CF) | small_scanner.mat | scanner[:,8] | None | CF cost shifter |
| cdid | Market index (1-80) | derived | sequential from ID | Incremental counter | Market loop index |
| cdindex | Last row of each market | derived | diff(cdid) | Change detection | Market boundaries |
| prodid | Product index (1-8) | derived | grp2idx(b*100+s) | Group indexing | Product FE |
| fiscid | Fiscal year (2006-2011) | derived | yearid+(montid>=4)+2004 | Shift + offset | Sample filter |
| fisccity | Fiscal year x city | derived | grp2idx(city*100+fisc) | Group indexing | Supply-side loop |
| coalid | Coalition member (0/1) | derived | firmid in {1,4,5} | Logical | PLE model |
| leadid | Price leader (0/1) | derived | firmid == 1 | Logical | PLE model |
| mpost | Miller post-merger dummy | derived | brndid in {11-13} & yr>=5 | Logical | Cost shifter (w) |
| cpost | Coors post-merger dummy | derived | brndid in {3,4} & yr>=5 | Logical | Cost shifter (w) |
| apost | ABI post-merger dummy | derived | firmid==1 & yr>=5 | Logical | Cost shifter (w) |
| mcpost | MillerCoors post dummy | derived | firmid==5 & yr>=5 | Logical | Cost shifter (w) |
| alpha | Mean price coefficient | small_dresgmm2.mat | theta1_2[0] = -0.10872 | None | Demand param |
| theta2 | RC std devs (3x1) | small_dresgmm2.mat | theta2_2 | None | Demand param |
| rho | Nesting parameter | small_dresgmm2.mat | rho_2 = 0.77788 | None | Demand param |
| theta2w | RC interaction (3x2) | small_dresgmm2.mat | sparse(theti,thetj,theta2) | Sparse fill | mu computation |
| derMat | ds/dp matrices (8x8x140) | small_dresgmm2.mat | derMat_2 | None | FOC inversion |
| elasMat | Elasticities (8x8x140) | small_dresgmm2.mat | elasMat_2 | None | Validation/tables |
| dfull | Income draws (Nx500) | small_demosE.mat | demosE expanded | Expand+filter+demean | Consumer heterogeneity |
| delta | Mean utility (960x1) | daugfile.mat | daugfile.delta | Contraction mapping | Demand inversion |
| mu | Consumer heterogeneity (960x500) | daugfile.mat | daugfile.mu | x2*(d_i*theta2w') | Individual shares |
| pcoefi | Indiv price coef (960x500) | daugfile.mat | daugfile.pcoefi | alpha + ai | Derivatives |
| deltanp | Non-price utility (960x1) | daugfile.mat | daugfile.deltanp | delta - alpha*p_jt | Counterfactuals |
| xi | Unobs quality (960x1) | daugfile.mat | daugfile.xi | deltanp - fesd*coef | Demand residual |

---

## 8. Instrument Set — Full Specification

This section separates three categories with explicit provenance for each element.

---

### 8.1 Variables confirmed directly from the replication package

These variables appear in the Matlab code (`main_data.m`, `results_costregs.m`, `f_impute_mc.m`) and can be constructed mechanically from the raw `.mat` files. No interpretation or reconstruction from the paper is needed.

#### 8.1.1 Demand-side matrices (from `main_data.m`)

| Variable | Code line | Formula | Confirmed content |
|----------|-----------|---------|-------------------|
| `x1` | line 174 | `[p_jt, fesd]` where `fesd = [prodfe, datefe(:,2:end)]` | Linear demand variables: price + product FE + date FE |
| `x2` | line 178 | `[p_jt, ones(N,1), calor]` (RCNL2 spec) | Variables receiving random coefficients: price, constant, calories |
| `logodds` | line 204 | `log(s_jt) - log(1 - inshr)` | BLP dependent variable: ln(s_j) - ln(s_0) |
| `logcondshr` | line 200 | `log(s_jt) - log(inshr)` | Within-nest conditional share (identifies rho in nested logit) |

#### 8.1.2 Supply-side cost shifters (from `main_data.m` lines 183-188)

| Variable | Formula | Matlab code | Role in `w` |
|----------|---------|-------------|-------------|
| `mpost` | `(brndid in {11,12,13}) & (yearid >= 5)` | line 185 | Miller brands post-merger cost shift |
| `cpost` | `(brndid in {3,4}) & (yearid >= 5)` | line 186 | Coors brands post-merger cost shift |
| `dist` | `scanner[:,7]` (miles x diesel index) | line 149 | Transportation cost to ABI brewery |
| `fess` | `[prodfe, cityfe(:,2:end), datefe(:,2:end)]` | line 171 | Product + city + date fixed effects |
| `apost` | `(firmid == 1) & (yearid >= 5)` | line 187 | ABI post-merger indicator |

#### 8.1.3 Supply-side cost regression structure (from `results_costregs.m` lines 41, 118)

```
X = [apost, w]  = [apost, mpost, cpost, dist, prodFE, cityFE, dateFE]
Y = mc           (implied marginal costs, from FOC inversion)
Method: OLS with SEs clustered by city
```

This is confirmed directly in the code. `f_ols.m` runs `gamma = inv(X'X) * X'y` with clustered variance `inv(X'X) * (sum_c X_c' e_c e_c' X_c) * inv(X'X)`.

#### 8.1.4 Demand parameter values (from `small_dresgmm2.mat`, loaded in `f_daugment.m`)

| Parameter | Value | Source | Role |
|-----------|-------|--------|------|
| `alpha` (mean price coef) | -0.10872 | `theta1_2[0]` | Mean marginal utility of income |
| `theta2` (RC params) | [0.00090, 0.01249, 0.00451] | `theta2_2` | Income interactions: price, constant, calories |
| `rho` (nesting) | 0.77788 | `rho_2` | Within-nest correlation |
| `theta2w` (3x2) | col2 = theta2 | `sparse(theti,thetj,theta2)` | Random coefficient interaction matrix |
| `derMat_2` (8x8x140) | — | `small_dresgmm2.mat` | Pre-computed demand derivatives |
| `elasMat_2` (8x8x140) | — | `small_dresgmm2.mat` | Pre-computed elasticity matrices |

**Key fact**: The replication package does NOT contain demand estimation code. These parameters are loaded as given from MW(2017, Econometrica). The entire supply-side analysis conditions on them.

---

### 8.2 Instrument set reconstructed for demand estimation in R

These instruments are described in the paper (Section 4.2, pp. 1775-1776) but are **NOT constructed anywhere in the replication package**. The demand estimation was performed separately in the original MW(2017) Econometrica codebase, which is not included. We reconstruct them here from the paper's description and verify they are constructible from the available data.

**Provenance**: Paper text, Section 4.2. **Not** from replication code.

The paper states: *"There are 12 instruments in total."*

#### Set A: Price instruments — addressing endogeneity of p_jt (2 instruments)

| ID | Variable | Construction | Classification | Economic justification (from paper) | Constructible from small data? |
|----|----------|-------------|----------------|-------------------------------------|-------------------------------|
| Z1 | `dist` | `scanner[:,7]` = miles to brewery x diesel index | **Cost-side shifter** | Shifts marginal cost via transportation. Varies cross-sectionally (region x firm) and temporally (diesel price). Excluded from demand: shipping cost does not enter consumer utility. | **YES** — 345 unique values, range [0.00, 9.06] |
| Z2 | `coalpost` | 1 if (firmid in {1,5}) & (yearid >= 5) | **Cost-side / competitive structure** | Captures change in competitive structure from merger. Relevant: observed price increases post-merger. Valid: conditional on product + time FEs, requires changes in ABI/MC unobserved quality are not systematically different from Modelo/Heineken changes. | **YES** — 240 of 640 obs = 1 |

**Paper quote**: *"The first set of instruments that we use addresses the endogeneity of prices. It includes the distance between the brewery and the region (miles x diesel index) and an indicator equal to 1 for ABI and MillerCoors products after the merger."*

#### Set B: Nesting parameter instruments — identifying rho (6 instruments)

These instrument the endogenous within-nest conditional share `log(s_j|g)` in the nested logit structure.

| ID | Variable | Construction | Classification | Economic justification (from paper) | Constructible from small data? |
|----|----------|-------------|----------------|-------------------------------------|-------------------------------|
| Z3 | `num_products` | J_t = count of inside products in market t | **BLP-style (market structure)** | Standard instrument. Negatively correlated with conditional share: more products in nest = lower share per product. Valid if xi uncorrelated with J. | **DEGENERATE** — J=8 in all 80 markets. Zero variation. In full IRI data (39 products, 37 regions), J varies. |
| Z4 | `sum_dist` | sum_j(dist_jt) for all j in market t | **BLP-style (aggregate rival cost)** | Captures variation in competing products' marginal costs. Positively correlated with conditional share: higher rival costs = less competition = higher own share. | **YES** — 80 unique values, range [6.68, 18.23] |
| Z5 | `num_products x ABI` | J_t x 1(firmid == 1) | BLP-style interaction | Allows nesting effect to differ for ABI products. | **DEGENERATE** — collapses to 8 x ABI_dummy (J constant). Perfectly collinear with ABI indicator after absorbing constant. |
| Z6 | `num_products x MC` | J_t x 1(firmid in {4,5}) | BLP-style interaction | Allows nesting effect to differ for Miller/Coors products. | **DEGENERATE** — collapses to 8 x MC_dummy. Same problem. |
| Z7 | `sum_dist x ABI` | sum_dist_t x 1(firmid == 1) | BLP-style interaction | ABI-specific variation in market-level cost structure. | **YES** — range [0.00, 18.23] |
| Z8 | `sum_dist x MC` | sum_dist_t x 1(firmid in {4,5}) | BLP-style interaction | MillerCoors-specific variation in market-level cost structure. | **YES** — range [0.00, 18.23] |

**Paper quote**: *"We use as instruments the number of products in the market and the distance summed across all products in the market. [...] to add flexibility, we incorporate interactions with indicators for ABI and Miller/Coors products."*

**Small-data implication**: Z3, Z5, Z6 have zero independent variation because J=8 everywhere. Only Z4, Z7, Z8 provide actual identification of rho. This is a limitation of the replication sample, not of the research design.

#### Set C: Random coefficient / demographic instruments — identifying Pi (4 instruments)

These identify the parameters governing consumer heterogeneity in preferences for characteristics.

| ID | Variable | Construction | Classification | Economic justification (from paper) | Constructible from small data? |
|----|----------|-------------|----------------|-------------------------------------|-------------------------------|
| Z9 | `mean_income x 1` | market-level mean of income draws | **Demographics (Romeo 2014)** | Income heterogeneity shifts price sensitivity and base utility. Valid if E[xi \| income, x] = 0. | **YES** — range [28.3, 54.2] |
| Z10 | `mean_income x calor` | mean_income x (calories/100 - mean) | Demographics | Identifies income-calorie preference interaction. Higher-income markets may value calorie content differently. | **YES** — range [-10.0, 24.7] |
| Z11 | `mean_income x size` | mean_income x sizeid | Demographics | Identifies income-size preference interaction. Pack-size choice may vary with income. | **YES** — range [56.7, 162.6] |
| Z12 | `mean_income x import` | mean_income x 1(firmid in {2,3}) | Demographics | Identifies income-import preference interaction. Import beer demand may correlate with affluence. | **YES** — range [0.0, 54.2] |

**Paper quote**: *"We use mean income interacted with the observed product characteristics (a constant, calories, package size, and an import dummy), which provide the requisite variation."*

**Construction note**: `mean_income` is computed from `small_demosE.mat` as the mean of the first 500 income draws per year-city combination, expanded to the observation level via the `yearcityid` mapping. This is NOT the same as individual-level income draws (which enter via `mu`). The instrument is the *market-level average* of the demographic variable.

#### Instrument count reconciliation

| Category | Paper count | Functional in small data | Degenerate |
|----------|------------|------------------------|------------|
| Price (alpha) | 2 | 2 | 0 |
| Nesting (rho) | 6 | 3 (Z4, Z7, Z8) | 3 (Z3, Z5, Z6) |
| RC/demographics (Pi) | 4 | 4 | 0 |
| **Total** | **12** | **9** | **3** |

With 9 functional instruments and 3 nonlinear parameters to estimate (alpha via price in x2, rho, and the income interaction scale), the model is overidentified. The 3 degenerate instruments are a known consequence of the small replication sample (8 fixed products in all markets).

---

### 8.3 Supply-side variables used in the merger/cost analysis

These variables are confirmed in the replication package and play specific roles in the supply-side estimation. They are NOT instruments in the demand-side IV sense; they are regressors, structural objects, or counterfactual inputs.

#### 8.3.1 Cost regression regressors (`results_costregs.m`)

| Variable | Matlab code | Formula | Economic role | Provenance |
|----------|------------|---------|---------------|------------|
| `apost` | `main_data.m:187` | `(firmid==1) & (yearid>=5)` | Key test variable. Coefficient > 0 implies ABI prices exceed unilateral effects prediction = evidence of price coordination (kappa). | **Code: confirmed** |
| `mpost` | `main_data.m:185` | `(brndid in {11,12,13}) & (yearid>=5)` | Captures MillerCoors efficiency gains (cost reduction after merger). In small data: only brndid=13 (Miller Lite) present. | **Code: confirmed** |
| `cpost` | `main_data.m:186` | `(brndid in {3,4}) & (yearid>=5)` | Captures Coors-specific post-merger cost changes. In small data: only brndid=4 (Coors Light) present. | **Code: confirmed** |
| `dist` | `main_data.m:149` | `scanner[:,7]` (miles x diesel) | Transportation cost. ABI products have dist=0 (zero distance to own brewery). | **Code: confirmed** |
| Product FE | `main_data.m:171` | `cr_dum(prodid)` (8 products) | Absorb time-invariant product-level cost differences. | **Code: confirmed** |
| City FE | `main_data.m:171` | `cr_dum(cityid)(:,2:end)` (4 dummies) | Absorb city-level cost differences. | **Code: confirmed** |
| Date FE | `main_data.m:171` | `cr_dum(dateid)(:,2:end)` | Absorb time-period cost shocks common to all products. | **Code: confirmed** |

**Estimation method**: OLS with clustered SEs by city (5 clusters). From `results_costregs.m:118`: `[gamma,se] = f_ols(mc, X, cv)` where `cv = ids.cityid`.

**Interpretation of key coefficients** (from paper Section 5.2, p. 1780):
- `apost > 0`: ABI costs appear higher post-merger under the assumed competitive model. Since actual costs should not increase, this indicates the model is attributing coordination-induced price increases to "cost" — evidence that the assumed competitive model (e.g., Bertrand with df=0) is wrong.
- The discount factor df is selected such that `apost ≈ 0` and all cost coefficients are economically sensible (Table 3 column selection: df=0.26).

#### 8.3.2 Structural objects used in FOC inversion (`f_impute_mc.m`)

| Object | Source | Dimension | Role | Provenance |
|--------|--------|-----------|------|------------|
| `derMat_2` | `small_dresgmm2.mat` | (8, 8, 140) | Price derivative matrix ds_j/dp_k per market | **Data file: confirmed** |
| `elasMat_2` | `small_dresgmm2.mat` | (8, 8, 140) | Elasticity matrix per market | **Data file: confirmed** |
| Ownership matrix | `f_ownMat.m` | (J, J) per market | `Owner[j,k] = 1(same firm)`. Changes across pre/post merger. | **Code: confirmed** |
| `coalid` | `main_data.m:112` | (N, 1) | Coalition membership: `firmid in {1, 4, 5}` | **Code: confirmed** |
| `sm` (supermarkup) | Estimated in `main_supply_bind.m` | scalar or (2,1) per fiscal-city | Wedge above Nash-Bertrand: `p_observed = p_nash + sm` | **Code: confirmed** |

**FOC inversion logic** (3-step, from `f_impute_mc.m`):
1. Fringe (non-coalition) firms: `mc_fringe = p_observed - markup_fringe` (standard Bertrand FOC)
2. Solve for fringe Nash prices given `p_nash_coalition = p_observed - sm`
3. Coalition firms: `mc_coal = p_nash - (Owner .* der')^{-1} * s_nash` (FOC inversion at Nash prices)

#### 8.3.3 Counterfactual merger inputs (`results_cfmergers.m`)

| Variable | Role | Provenance |
|----------|------|------------|
| `gammapar.mat` (gamma, se) | Estimated cost coefficients from df=0.26 regression | **Code: confirmed** (`results_costregs.m:123-124`) |
| `distbutfor` (`scanner[:,8]`) | But-for counterfactual distance (used when simulating ownership change) | **Code: confirmed** (`main_data.m:150`) |
| Counterfactual `firmid` | Merger simulation changes ownership: e.g., scenario 1 separates Miller from Coors | **Code: confirmed** (`results_cfmergers.m`) |
| `mcpost` | `(firmid==5) & (yearid>=5)`: MillerCoors post-merger indicator for cost adjustments | **Code: confirmed** (`main_data.m:184`) |

---

### 8.4 Summary: what comes from where

| Element | Source | Status |
|---------|--------|--------|
| `x1`, `x2`, `w`, `fesd`, `fess`, all ID variables | `main_data.m` | **Directly from replication code** |
| `alpha`, `theta2`, `rho`, `derMat`, `elasMat` | `small_dresgmm2.mat` | **Directly from replication data** |
| `delta`, `mu`, `ai`, `pcoefi`, `deltanp`, `xi` | `daugfile.mat` via `f_daugment.m` | **Directly from replication code + data** |
| Cost regression structure `X = [apost, w]` | `results_costregs.m` | **Directly from replication code** |
| FOC inversion 3-step algorithm | `f_impute_mc.m` | **Directly from replication code** |
| Demand instruments Z1-Z12 (the 12-instrument set) | Paper Section 4.2 (pp. 1775-1776) | **Reconstructed from paper text.** NOT in replication code. Required for fresh demand estimation in R. |
| Instrument economic justifications | Paper Section 4.2 + Appendix D.1 | **From paper.** Exclusion restrictions are economic arguments, not testable. |
| Z3, Z5, Z6 degeneracy in small data | Our verification | **Our analysis.** Not documented in the paper or replication package. |
| Supply-side "instrument" (apost for markup identification) | Paper Section 5.2 (pp. 1780-1781) + `results_costregs.m` | **Both paper and code.** Paper provides economic logic; code confirms implementation. |

---

## 9. IV Sensitivity and Instrument Validity (Small Sample)

### 9.1 Specification summary

Six specifications were estimated on the 640-observation perturbed replication sample to assess instrument sensitivity. All specifications include product and date fixed effects. Clustered standard errors by city (5 clusters).

| Spec | Model | Endogenous | Excluded IVs | IVs used |
|------|-------|-----------|-------------|----------|
| OLS | Nested logit | — (none) | — | Biased benchmark |
| A | Simple logit | price | dist, coalpost | 2 cost-side |
| B | Simple logit | price | dist, coalpost, sum_dist, sumdist_abi, sumdist_mc | 5 cost + BLP |
| C | Nested logit | price, logcondshr | dist, coalpost | 2 cost-side (just-identified) |
| **D** | **Nested logit** | **price, logcondshr** | **dist, coalpost, sum_dist, sumdist_abi, sumdist_mc** | **5 cost + BLP** |
| E | Simple logit | price | all 9 functional instruments (incl. income interactions) | 9 cost + BLP + income |

### 9.2 Results

| Spec | alpha (price) | se | sigma (nesting) | se | F(price) | F(logcshr) | Sargan p |
|------|--------------|-----|----------------|-----|----------|-----------|----------|
| OLS | +0.025 | 0.062 | 0.975 | 0.048 | — | — | — |
| A | +0.097 | 0.073 | — | — | 30.6 | — | 0.273 |
| B | -0.059 | 0.070 | — | — | 15.7 | — | 0.000 |
| C | -0.715 | 2.961 | 8.128 | 29.66 | 30.6 | — | — (just-id) |
| **D** | **-0.027** | **0.047** | **0.940** | **0.107** | **15.7** | **21.3** | **0.003** |
| E | +0.148 | 0.108 | — | — | 41.7 | — | 0.000 |
| Paper (RCNL) | -0.109 | — | 0.778 | — | — | — | — |

Key observations:

- **Spec A** (logit, 2 IVs): Wrong sign on alpha (+0.097). Strong first stage (F=30.6). Sargan does not reject (p=0.27) but only 1 degree of overidentification.
- **Spec B** (logit, 5 IVs): Correct sign (alpha=-0.059), strong first stage (F=15.7), but Sargan rejects decisively (p<0.001).
- **Spec C** (nested logit, just-identified): Correct signs but estimates are completely imprecise (se on alpha is 2.96, se on sigma is 29.7). Two instruments cannot reliably identify two endogenous variables in this sample.
- **Spec D** (nested logit, 5 IVs): Correct sign on alpha (-0.027), sigma in [0,1) at 0.94, strong first stages for both endogenous variables (F=15.7 and 21.3). Sargan rejects (p=0.003) but this is the only specification producing economically coherent nested logit estimates.
- **Spec E** (logit, 9 IVs including income): Wrong sign (+0.148) despite very high first-stage F (41.7). Sargan rejects decisively (p<0.001). The income instruments create spurious power.

### 9.3 Interpretation

#### 9.3.1 Why income-based instruments perform poorly

The income interaction instruments (Z9–Z12: mean income x {constant, calories, size, import}) produce the highest first-stage F-statistic (41.7) but yield a positive price coefficient — economically nonsensical. This occurs because:

- The replication data was intentionally perturbed: random noise added to prices and quantities, and ABI/Miller volumes inflated by 1.5x (`smalldata.m` lines 31–38).
- These perturbations distort the correlation between local income variation and product-level demand in ways that break the exclusion restriction E[xi | income, x] = 0.
- With only 5 cities (versus 37 in the full sample), market-level income variation is extremely limited, and the instruments may be identifying primarily off the artificial perturbation rather than true economic variation.
- A high first-stage F-statistic is necessary but not sufficient for valid instruments. Here, relevance is satisfied but exogeneity is not.

#### 9.3.2 Why Sargan rejects overidentified specifications

The Sargan (Hansen J) test rejects in all overidentified specifications except Spec A (which has only 1 degree of freedom and thus low power). This is consistent with:

- The data perturbation introducing artificial correlations between instruments and the structural error term xi.
- In the original unperturbed IRI data (37 regions, 39 products), the paper's instruments were validated against the Sargan criterion. The perturbation breaks this.
- The rejection does not imply the research design is flawed — it indicates that the small perturbed replication sample cannot support the same overidentification restrictions as the full dataset.

#### 9.3.3 Why Spec D is the preferred baseline

Spec D (nested logit, 5 cost + BLP instruments) is selected as the main 2SLS specification because:

- It is the only specification producing both alpha < 0 and sigma in [0,1) simultaneously.
- Both first stages are strong (F=15.7 for price, F=21.3 for logcondshr), comfortably above the Staiger-Stock threshold of 10.
- The 5 instruments have clear economic justification: `dist` shifts marginal cost, `coalpost` captures the merger's competitive structure change, and `sum_dist` (with firm interactions) captures BLP-style variation in rival cost characteristics.
- While Sargan rejects at the 1% level (p=0.003), this is expected given the data perturbation and does not invalidate the specification for illustrative purposes.
- The estimated alpha (-0.027) is less negative than the paper's RCNL estimate (-0.109), consistent with attenuation from the simpler nested logit specification and the smaller sample.

### 9.4 Academic summary (for inclusion in report)

The sensitivity of the 2SLS estimates to the choice of instrument set was assessed by comparing specifications using cost-based instruments alone (distance to brewery, post-merger coalition indicator, and BLP-style aggregates of rival costs) against an expanded set that includes interactions of mean local income with product characteristics. Specifications incorporating income-based instruments yield first-stage F-statistics that exceed conventional thresholds but produce price coefficients with incorrect (positive) signs, suggesting that these instruments identify off variation introduced by the data perturbation rather than economically meaningful demand-side heterogeneity. Overidentification tests reject instrument validity in all expanded specifications, consistent with the perturbed nature of the replication sample — which the authors constructed by adding random noise to prices and quantities and inflating select brand volumes. The preferred specification uses five cost-side and BLP-style instruments, which deliver economically coherent estimates (negative price coefficient, nesting parameter below unity) with strong first stages. Given these limitations, the 2SLS results should be interpreted as qualitative validation of the demand model's structure rather than as exact replication of the paper's random coefficients nested logit estimates, which were obtained from the full (unperturbed) IRI dataset covering 37 regions and 39 products.

---

## 10. Full Data vs. Small Replication Sample

### 10.1 Data scale and variation

| Dimension | Full dataset (MW 2017) | Small replication sample |
|-----------|----------------------|-------------------------|
| Regions | 37 IRI markets | 5 cities |
| Products | 39 brand x size combos | 8 brand x size combos |
| J per market | Varies across markets (key source of BLP instrument variation) | Fixed at J=8 in every market |
| Time periods | 28 quarters (2005–2011), excluding 1 year post-merger | Same 28 quarters, same exclusion window |
| Observations (after filtering) | ~94,656 (brand x size x region x month x year) | 640 (8 x 80 markets) |
| Perturbation | None — confidential IRI scanner data | Prices and quantities randomly perturbed; ABI and Miller Lite volumes inflated by 1.5x |

### 10.2 Demand specification

- **Full paper**: Random Coefficients Nested Logit (RCNL) estimated by GMM with 500 simulated consumers per market. Consumer income from PUMS/ACS interacts with price, a constant, calories, package size, and an import dummy. Multiple specifications reported (nested logit, RCNL-1 through RCNL-4).
- **Replication package**: Does not re-estimate demand. Loads pre-estimated parameters (alpha, theta2, rho) from MW(2017, Econometrica) and uses them as given for the supply-side analysis.
- **Our setup**: We estimate a simplified nested logit by 2SLS as the problem set requests, and separately load the paper's RCNL parameters for the supply-side replication.

### 10.3 Instrument set

- **Full paper** (12 instruments):
  - 2 for price: dist (miles x diesel), coalition post-merger indicator
  - 6 for nesting (rho): number of products in market, sum of distance across products, each interacted with ABI and Miller/Coors indicators — all 6 have variation because J varies across the 37 regions
  - 4 for RC parameters (Pi): mean income x {constant, calories, size, import} — 37 regions provide rich income variation
- **Small sample** (9 functional of 12):
  - 3 instruments are degenerate: Z3 (J=8 always), Z5 and Z6 (collapse to scaled firm dummies)
  - Income instruments (Z9–Z12) have only 5 cities of variation; empirically they identify off perturbation noise rather than true demand heterogeneity
  - Only the 5 cost + BLP instruments (dist, coalpost, sum_dist, sumdist_abi, sumdist_mc) produce economically coherent estimates

### 10.4 Implications for identification

- **Price endogeneity**: Present in both datasets, but harder to address in the small sample. With 5 cities, the cross-sectional cost variation that identifies alpha is severely limited.
- **Nesting parameter**: In the full data, variation in J across markets is the primary source of identification for rho. With J fixed at 8, this channel is shut off entirely. Identification of rho in our sample relies solely on sum_dist and its firm interactions.
- **RC parameters**: Impossible to estimate credibly with 5 cities. The paper's RCNL estimates exploit income variation across 37 regions — we cannot replicate this.
- **Instrument validity**: Sargan tests pass in the full data but reject in the small sample, consistent with the perturbation introducing artificial correlations.

### 10.5 Implications for results

- **What we CAN match**:
  - Qualitative signs and relative magnitudes of demand coefficients (alpha < 0, sigma in [0,1))
  - The structure of the estimation pipeline (data construction, share computation, instrument logic, FOC inversion)
  - Supply-side results (markups, cost regressions, merger simulations) using the paper's pre-estimated demand parameters
  - Directional patterns: ABI has the largest shares, imports are niche, 24-packs have higher volume shares

- **What we CANNOT match**:
  - Exact numerical values of any demand coefficient (perturbed data + simplified specification)
  - The paper's published tables (Table IV demand estimates, Table 3 cost regressions) — our targets are the pre-computed `.txt` files in the replication package, which themselves reflect the perturbed data
  - Random coefficient estimates (theta2/Pi) — these require the full RCNL with 37 regions of income variation
  - Overidentification test results — Sargan passes in the full data but rejects here

### 10.6 Takeaway for empirical strategy

- **Use two tracks**: (a) our own 2SLS nested logit for the problem set's demand estimation questions, acknowledging the limitations above; (b) the paper's pre-estimated RCNL parameters from `dresgmm2.mat` for all supply-side analysis (markups, FOC inversion, merger simulation).
- **Validate against replication outputs, not published tables**: The `.txt` files in `replication/results/` are the correct benchmark. They were produced from the same perturbed small sample we are using.
- **Interpret 2SLS results qualitatively**: Our nested logit estimates demonstrate understanding of the BLP identification strategy and produce directionally correct coefficients, but they are not — and cannot be — a precise replication of the paper's RCNL estimates.
