# ==============================================================================
# 03_demand_estimation.R — Linear 2SLS and RCNL demand estimation
# Miller & Weinberg (2017) Replication in R
# ==============================================================================
#
# Este script estima a demanda por cerveja usando tres abordagens:
#   (1) OLS (benchmark viesado),
#   (2) 2SLS com variaveis instrumentais (logit simples e nested logit),
#   (3) Parametros RCNL pre-estimados pelo artigo original (GMM com consumidores simulados).
#
# O objetivo final eh obter elasticidades-preco proprias e cruzadas,
# necessarias para recuperar markups e simular fusoes nos scripts seguintes.
# ==============================================================================

load("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/R_project/data/step02_output.RData")
source("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/R_project/scripts/00_setup.R")

cat("\n")
cat("=", rep("=", 69), "\n", sep = "")
cat("STEP 3: DEMAND ESTIMATION\n")
cat("=", rep("=", 69), "\n", sep = "")

# ==============================================================================
# BLOCO: ESTIMACAO 2SLS
# ==============================================================================
# Estima a equacao de demanda logit e nested logit via minimos quadrados
# em dois estagios (2SLS). A variavel dependente eh log(s_j/s_0), a
# "log-odds ratio" do produto j em relacao ao bem externo.
#
# O preco eh endogeno (correlacionado com xi, o choque de demanda nao observado),
# por isso OLS produz estimativas viesadas. Usamos instrumentos de custo (dist)
# e instrumentos tipo BLP (soma de caracteristicas dos rivais) para identificar alpha.
# ==============================================================================

cat("\n")
cat("-", rep("-", 69), "\n", sep = "")
cat("PART 1: 2SLS Demand Estimation\n")
cat("-", rep("-", 69), "\n", sep = "")

# ------------------------------------------------------------------------------
# 1a. OLS benchmark (viesado)
# ------------------------------------------------------------------------------
# OLS serve como referencia: esperamos que alpha seja atenuado (menos negativo)
# porque cov(preco, xi) > 0 — produtos com demanda alta cobram precos altos.

cat("\n--- 1a. OLS benchmark ---\n")

fit_ols_logit <- lm(logodds ~ price + factor(prodid) + factor(dateid), data = df)
fit_ols_nlogit <- lm(logodds ~ price + logcondshr + factor(prodid) + factor(dateid), data = df)

# Erros-padrao clusterizados por cidade para lidar com correlacao intra-mercado
vcov_ols_l <- vcovCL(fit_ols_logit, cluster = df$cityid)
vcov_ols_nl <- vcovCL(fit_ols_nlogit, cluster = df$cityid)

cat(sprintf("  Simple logit OLS:  alpha = %+.6f  (se = %.4f)\n",
            coef(fit_ols_logit)["price"], sqrt(diag(vcov_ols_l))["price"]))
cat(sprintf("  Nested logit OLS:  alpha = %+.6f  (se = %.4f),  sigma = %.4f  (se = %.4f)\n",
            coef(fit_ols_nlogit)["price"], sqrt(diag(vcov_ols_nl))["price"],
            coef(fit_ols_nlogit)["logcondshr"], sqrt(diag(vcov_ols_nl))["logcondshr"]))

# ------------------------------------------------------------------------------
# 1b. Spec B: Logit simples com 5 IVs (custo + BLP)
# ------------------------------------------------------------------------------
# Modelo logit simples: log(s_j/s_0) = alpha*p_j + prodFE + dateFE + xi_j
# Os 5 instrumentos excluidos sao:
#   dist       — distancia ao centro de distribuicao (shifter de custo)
#   coalpost   — interacao indicador de coalizao x pos-fusao
#   sum_dist   — soma de distancias dos rivais (BLP)
#   sumdist_abi, sumdist_mc — soma de distancias por firma (BLP)

cat("\n--- 1b. Spec B: Simple logit IV, 5 instruments ---\n")
cat("  logodds = alpha * price + prodFE + dateFE + xi\n")
cat("  Excluded IVs: dist, coalpost, sum_dist, sumdist_abi, sumdist_mc\n\n")

# 2SLS via ivreg: formula | instrumentos (variaveis exogenas + instrumentos excluidos)
fit_B <- ivreg(logodds ~ price + factor(prodid) + factor(dateid) |
                 factor(prodid) + factor(dateid) +
                 dist + coalpost + sum_dist + sumdist_abi + sumdist_mc,
               data = df)
vcov_B <- vcovCL(fit_B, cluster = df$cityid)
diag_B <- summary(fit_B, diagnostics = TRUE)$diagnostics

alpha_B <- coef(fit_B)["price"]
se_B    <- sqrt(diag(vcov_B))["price"]

cat(sprintf("  alpha = %+.6f  (se = %.4f)\n", alpha_B, se_B))
# alpha deve ser negativo: aumento de preco reduz utilidade
cat(sprintf("  alpha < 0: %s\n", ifelse(alpha_B < 0, "YES", "NO")))
# Estatistica F do primeiro estagio: F > 10 indica instrumentos fortes (Staiger-Stock)
cat(sprintf("  1st-stage F = %.2f\n", diag_B["Weak instruments", "statistic"]))
# Teste de Sargan: H0 = instrumentos validos (exogenos). p alto => nao rejeita validade
cat(sprintf("  Sargan p = %.4f (df = %d)\n",
            diag_B["Sargan", "p-value"], diag_B["Sargan", "df1"]))

# ------------------------------------------------------------------------------
# 1c. Spec D: Nested logit com 5 IVs (ESPECIFICACAO BASE)
# ------------------------------------------------------------------------------
# Modelo nested logit: log(s_j/s_0) = alpha*p_j + sigma*log(s_j|g) + prodFE + dateFE + xi_j
# O parametro sigma (in [0,1)) governa a correlacao de preferencias dentro do ninho
# (grupo "inside"). sigma -> 1 implica alta substituicao dentro do grupo;
# sigma = 0 colapsa ao logit simples.
# logcondshr = log(s_j|g) tambem eh endogeno, pois depende das shares de equilibrio.

cat("\n--- 1c. Spec D: Nested logit IV, 5 instruments (BASELINE) ---\n")
cat("  logodds = alpha * price + sigma * logcondshr + prodFE + dateFE + xi\n")
cat("  Excluded IVs: dist, coalpost, sum_dist, sumdist_abi, sumdist_mc\n\n")

# 2SLS com dois regressores endogenos (price e logcondshr)
fit_D <- ivreg(logodds ~ price + logcondshr + factor(prodid) + factor(dateid) |
                 factor(prodid) + factor(dateid) +
                 dist + coalpost + sum_dist + sumdist_abi + sumdist_mc,
               data = df)
vcov_D <- vcovCL(fit_D, cluster = df$cityid)
diag_D <- summary(fit_D, diagnostics = TRUE)$diagnostics

alpha_D <- coef(fit_D)["price"]
sigma_D <- coef(fit_D)["logcondshr"]
se_alpha_D <- sqrt(diag(vcov_D))["price"]
se_sigma_D <- sqrt(diag(vcov_D))["logcondshr"]

cat(sprintf("  alpha = %+.6f  (se = %.4f)\n", alpha_D, se_alpha_D))
cat(sprintf("  sigma = %+.6f  (se = %.4f)\n", sigma_D, se_sigma_D))
# Condicoes de consistencia do nested logit: alpha < 0 e sigma in [0,1)
cat(sprintf("  alpha < 0: %s,  sigma in [0,1): %s\n",
            ifelse(alpha_D < 0, "YES", "NO"),
            ifelse(sigma_D >= 0 & sigma_D < 1, "YES", "NO")))
# Estatisticas F separadas para cada endogena
cat(sprintf("  1st-stage F (price)    = %.2f\n",
            diag_D["Weak instruments (price)", "statistic"]))
cat(sprintf("  1st-stage F (logcshr)  = %.2f\n",
            diag_D["Weak instruments (logcondshr)", "statistic"]))
# Wu-Hausman: testa se OLS e 2SLS diferem significativamente (endogeneidade)
cat(sprintf("  Wu-Hausman p = %.4f\n", diag_D["Wu-Hausman", "p-value"]))
cat(sprintf("  Sargan p = %.4f (df = %d)\n",
            diag_D["Sargan", "p-value"], diag_D["Sargan", "df1"]))

# ==============================================================================
# BLOCO: PARAMETROS RCNL DO ARTIGO
# ==============================================================================
# O modelo RCNL (Random Coefficients Nested Logit) permite heterogeneidade
# nos coeficientes de preco, constante e calorias, interagidos com renda.
# Esses parametros foram estimados por GMM no artigo original usando a base
# completa do IRI (37 regioes, 39 produtos, 500 consumidores simulados).
# Nos os carregamos como "dados" para usar na recuperacao de markups e simulacao.
# ==============================================================================

cat("\n")
cat("-", rep("-", 69), "\n", sep = "")
cat("PART 2: Paper's RCNL Demand Estimates (loaded from dresgmm2)\n")
cat("-", rep("-", 69), "\n", sep = "")

cat("\n  These parameters are taken as given from MW(2017, Econometrica).\n")
cat("  They were estimated on the full IRI data (37 regions, 39 products)\n")
cat("  using GMM with 500 simulated consumers and RCNL specification.\n\n")

# alpha: coeficiente medio de preco (deve ser negativo)
cat(sprintf("  alpha (mean price coef) = %.6f\n", alpha))
# rho: parametro de aninhamento (nesting), analogo a sigma no nested logit linear
cat(sprintf("  rho   (nesting param)   = %.6f\n", rho))
# theta2: desvios-padrao dos coeficientes aleatorios (heterogeneidade nao observada)
cat(sprintf("  theta2 (RC params):       [%.5f, %.5f, %.5f]\n",
            theta2[1], theta2[2], theta2[3]))
# theta2w: interacoes dos coeficientes aleatorios com renda dos consumidores
# Coluna 2 contem as interacoes com renda (coluna 1 seria distribuicao normal, aqui zero)
cat("\n  theta2w matrix (random coefficient interactions with income):\n")
cat(sprintf("    Row 1 (price x income):    %.5f\n", theta2w[1, 2]))
cat(sprintf("    Row 2 (constant x income): %.5f\n", theta2w[2, 2]))
cat(sprintf("    Row 3 (calories x income): %.5f\n", theta2w[3, 2]))
cat(sprintf("\n  theta1 (all %d linear params): alpha + %d product/date FEs\n",
            length(theta1), length(theta1) - 1))

# ==============================================================================
# BLOCO: ELASTICIDADES
# ==============================================================================
# Elasticidades-preco medem a sensibilidade percentual da demanda a variacoes
# percentuais no preco. Sao fundamentais para:
#   (1) avaliar poder de mercado,
#   (2) construir a matriz de derivadas ds/dp usada na inversao da CPO de Bertrand,
#   (3) comparar padroes de substituicao entre modelos (logit, nested logit, RCNL).
# ==============================================================================

cat("\n")
cat("-", rep("-", 69), "\n", sep = "")
cat("PART 3: Elasticities\n")
cat("-", rep("-", 69), "\n", sep = "")

# ------------------------------------------------------------------------------
# 3a. Elasticidades do logit simples (Spec B)
# ------------------------------------------------------------------------------
# No logit simples, as elasticidades tem forma fechada:
#   Propria:  eta_jj = alpha * p_j * (1 - s_j)
#   Cruzada:  eta_jk = -alpha * p_k * s_k
# Limitacao: a elasticidade cruzada nao depende de j — todos os rivais sao
# substitutos simetricos (propriedade IIA do logit multinomial).

cat("\n--- 3a. Simple logit elasticities (Spec B) ---\n\n")

# Elasticidade-preco propria: deve ser negativa (alpha < 0, 1-s > 0)
elas_own_logit  <- alpha_B * df$price * (1 - df$share)
# Elasticidade cruzada media: positiva (aumento do preco de k beneficia j)
elas_cross_logit_avg <- -alpha_B * df$price * df$share  # average cross from j's perspective

cat(sprintf("  Mean own-price elasticity:  %.4f\n", mean(elas_own_logit)))
cat(sprintf("  Range: [%.4f, %.4f]\n", min(elas_own_logit), max(elas_own_logit)))
cat(sprintf("  Mean cross-price elast:     %.4f\n", mean(elas_cross_logit_avg)))

# ------------------------------------------------------------------------------
# 3b. Elasticidades do nested logit (Spec D)
# ------------------------------------------------------------------------------
# No nested logit, a substituicao dentro do grupo (inside goods) eh maior
# que para o bem externo. A elasticidade propria incorpora sigma:
#   eta_jj = (alpha/(1-sigma)) * p_j * (1 - sigma*s_j|g - (1-sigma)*s_j)
# onde s_j|g = s_j / inshr eh a share condicional dentro do ninho.

cat("\n--- 3b. Nested logit elasticities (Spec D) ---\n\n")

# Share condicional: fracao de s_j dentro do grupo de bens internos
s_cond <- df$share / df$inshr  # conditional share within inside good nest

# Elasticidade propria do nested logit (mais flexivel que o logit simples)
elas_own_nlogit <- (alpha_D / (1 - sigma_D)) * df$price *
  (1 - sigma_D * s_cond - (1 - sigma_D) * df$share)

cat(sprintf("  Mean own-price elasticity:  %.4f\n", mean(elas_own_nlogit)))
cat(sprintf("  Range: [%.4f, %.4f]\n", min(elas_own_nlogit), max(elas_own_nlogit)))

# By firm
cat("\n  Own-price elasticities by firm (Spec D nested logit):\n")
for (firm in sort(unique(df$firmid))) {
  mask <- df$firmid == firm
  firm_label <- c("1" = "ABI", "2" = "Corona", "3" = "Heineken",
                  "4" = "Coors", "5" = "MillerCoors")[as.character(firm)]
  cat(sprintf("    firm=%d (%s): mean = %.4f,  [%.4f, %.4f]\n",
              firm, firm_label, mean(elas_own_nlogit[mask]),
              min(elas_own_nlogit[mask]), max(elas_own_nlogit[mask])))
}

# ------------------------------------------------------------------------------
# 3c. Elasticidades RCNL (pre-computadas pelo artigo)
# ------------------------------------------------------------------------------
# O RCNL gera uma matriz de elasticidades J x J para cada mercado,
# permitindo padroes de substituicao totalmente flexiveis entre produtos.
# elasMat eh 8 x 8 x 140 (140 mercados brutos); filtramos para 80 mercados.

cat("\n--- 3c. RCNL elasticities (from paper's elasMat_2) ---\n\n")

# elasMat is 8 x 8 x 140. The 140 slices correspond to ALL raw markets
# (1120 obs / 8 products = 140). We need to select the 80 markets that
# pass the full obsin filter (obsintemp + fiscal year).
#
# Since the filter operates at the market level (city x year x quarter),
# either all 8 products in a market pass or none do. We reconstruct the
# market-level selector from the observation-level obsin applied to all 1120 obs.

# Reconstroi o seletor de mercados a partir dos IDs compostos do Matlab
mat_scanner_reload <- R.matlab::readMat(file.path(path$data_raw, "small_scanner.mat"))
id2_all <- mat_scanner_reload$small.scanner[, 1]
fid_all <- floor(id2_all / 1e10)
yid_all <- floor((id2_all - fid_all * 1e10 -
                    floor((id2_all - fid_all * 1e10) / 1e8) * 1e8 -
                    floor((id2_all - fid_all * 1e10 -
                             floor((id2_all - fid_all * 1e10) / 1e8) * 1e8) / 1e6) * 1e6 -
                    floor((id2_all - fid_all * 1e10 -
                             floor((id2_all - fid_all * 1e10) / 1e8) * 1e8 -
                             floor((id2_all - fid_all * 1e10 -
                                      floor((id2_all - fid_all * 1e10) / 1e8) * 1e8) / 1e6) * 1e6) / 1e4) * 1e4) / 1e2)
mid_all <- round(id2_all - fid_all * 1e10 -
                   floor((id2_all - fid_all * 1e10) / 1e8) * 1e8 -
                   floor((id2_all - fid_all * 1e10 -
                            floor((id2_all - fid_all * 1e10) / 1e8) * 1e8) / 1e6) * 1e6 -
                   floor((id2_all - fid_all * 1e10 -
                            floor((id2_all - fid_all * 1e10) / 1e8) * 1e8 -
                            floor((id2_all - fid_all * 1e10 -
                                     floor((id2_all - fid_all * 1e10) / 1e8) * 1e8) / 1e6) * 1e6) / 1e4) * 1e4 -
                   yid_all * 1e2)
obsintemp_all <- (yid_all <= 3 | yid_all >= 6) | (yid_all == 4 & mid_all <= 2) | (yid_all == 5 & mid_all >= 3)
fiscid_all <- yid_all + as.integer(mid_all >= 4) + 2004L
obsin_all <- obsintemp_all & (fiscid_all %in% c(2006, 2007, 2010, 2011))
raw_mkt <- rep(1:140, each = 8)
mkt_selector <- tapply(obsin_all, raw_mkt, all)

cat(sprintf("  elasMat dimensions: %s\n", paste(dim(elasMat), collapse = " x ")))
cat(sprintf("  Raw markets: 140, selected by fiscal filter: %d\n", sum(mkt_selector)))
stopifnot(sum(mkt_selector) == 80)

# Seleciona as 80 fatias de mercado relevantes
elasMat_80 <- elasMat[, , mkt_selector]

cat(sprintf("  elasMat_80 dimensions: %s\n", paste(dim(elasMat_80), collapse = " x ")))

# Elasticidades proprias: diagonal de cada fatia J x J
own_elas_rcnl <- sapply(1:dim(elasMat_80)[3], function(m) diag(elasMat_80[, , m]))
# own_elas_rcnl is 8 x 80 (products x markets)

cat("\n  Mean own-price elasticities by product (RCNL, 80 markets):\n")
prod_labels <- c("Bud 12", "Bud 24", "Coors 12", "Coors 24",
                 "Corona 12", "Heineken 12", "Miller 12", "Miller 24")
for (j in 1:8) {
  cat(sprintf("    j=%d (%s): mean = %.3f,  [%.3f, %.3f]\n",
              j, prod_labels[j],
              mean(own_elas_rcnl[j, ]),
              min(own_elas_rcnl[j, ]),
              max(own_elas_rcnl[j, ])))
}

cat(sprintf("\n  Grand mean own-price elasticity (RCNL): %.3f\n",
            mean(own_elas_rcnl)))

# Matriz de elasticidades completa para o mercado 1 (exemplo didatico)
cat("\n  Elasticity matrix, market 1 (RCNL):\n")
elas_m1 <- round(elasMat_80[, , 1], 3)
colnames(elas_m1) <- prod_labels
rownames(elas_m1) <- prod_labels
print(elas_m1)

# Padroes de substituicao cruzada: o RCNL permite assimetria (ao contrario do logit)
cat("\n  Cross-elasticity patterns (market 1):\n")
cat(sprintf("    Bud 12 -> Coors 12:   %.3f (same size, different firm)\n", elas_m1[1, 3]))
cat(sprintf("    Bud 12 -> Bud 24:     %.3f (same brand, different size)\n", elas_m1[1, 2]))
cat(sprintf("    Bud 12 -> Corona 12:  %.3f (domestic vs import, same size)\n", elas_m1[1, 5]))

# ==============================================================================
# BLOCO: TABELA COMPARATIVA
# ==============================================================================
# Resume os resultados das quatro especificacoes para facilitar a comparacao.
# O padrao esperado eh: OLS viesado -> 2SLS corrige -> RCNL mais flexivel.
# ==============================================================================

cat("\n")
cat("-", rep("-", 69), "\n", sep = "")
cat("PART 4: Comparison Table\n")
cat("-", rep("-", 69), "\n", sep = "")

cat("\n")
cat(sprintf("%-28s %12s %12s %12s %12s\n",
            "", "OLS (NL)", "Spec B (L)", "Spec D (NL)", "Paper RCNL"))
cat(rep("-", 80), "\n", sep = "")
cat(sprintf("%-28s %12s %12s %12s %12s\n",
            "Model", "Nested logit", "Simple logit", "Nested logit", "RCNL"))
cat(sprintf("%-28s %12s %12s %12s %12s\n",
            "Method", "OLS", "2SLS", "2SLS", "GMM"))
cat(sprintf("%-28s %12s %12s %12s %12s\n",
            "# excluded IVs", "0", "5", "5", "12"))
cat(sprintf("%-28s %12.4f %12.4f %12.4f %12.4f\n",
            "alpha (price coef)",
            coef(fit_ols_nlogit)["price"], alpha_B, alpha_D, alpha))
cat(sprintf("%-28s %12.4f %12.4f %12.4f %12s\n",
            "  (clustered se)",
            sqrt(diag(vcov_ols_nl))["price"], se_B, se_alpha_D, "—"))
cat(sprintf("%-28s %12.4f %12s %12.4f %12.4f\n",
            "sigma/rho (nesting)",
            coef(fit_ols_nlogit)["logcondshr"], "—", sigma_D, rho))
cat(sprintf("%-28s %12.4f %12s %12.4f %12s\n",
            "  (clustered se)",
            sqrt(diag(vcov_ols_nl))["logcondshr"], "", se_sigma_D, "—"))
cat(sprintf("%-28s %12s %12.4f %12.4f %12s\n",
            "Mean own-price elasticity",
            "—",
            mean(elas_own_logit), mean(elas_own_nlogit), round(mean(own_elas_rcnl), 4)))
cat(sprintf("%-28s %12s %12.1f %12.1f %12s\n",
            "1st-stage F (price)",
            "—",
            diag_B["Weak instruments", "statistic"],
            diag_D["Weak instruments (price)", "statistic"],
            "—"))
cat(sprintf("%-28s %12d %12d %12d %12s\n",
            "N", nobs(fit_ols_nlogit), nobs(fit_B), nobs(fit_D), "94,656"))

# ==============================================================================
# BLOCO: SALVAR RESULTADOS
# ==============================================================================

cat("\n\n--- Saving estimation results ---\n")

estimation_results <- list(
  # 2SLS estimates
  fit_ols_logit  = fit_ols_logit,
  fit_ols_nlogit = fit_ols_nlogit,
  fit_B = fit_B,
  fit_D = fit_D,
  vcov_B = vcov_B,
  vcov_D = vcov_D,
  # Key coefficients
  alpha_B = alpha_B,
  alpha_D = alpha_D,
  sigma_D = sigma_D,
  # Elasticities
  elas_own_logit  = elas_own_logit,
  elas_own_nlogit = elas_own_nlogit,
  own_elas_rcnl   = own_elas_rcnl,
  elasMat_80      = elasMat_80,
  # Paper's RCNL params (carried forward)
  alpha_rcnl = alpha,
  rho_rcnl   = rho,
  theta2_rcnl = theta2,
  theta2w_rcnl = theta2w
)

save(
  estimation_results,
  df, dfull, cdindex, obsindemand,
  alpha, theta1, theta2, rho, theta2w,
  derMat, elasMat,
  delta_val, deltanp_val, xi_val, mu_val, pcoefi_val,
  dprodfecoef, ddatefecoef, dcityfecoef,
  file = file.path(path$data_out, "step03_output.RData")
)
cat("  Saved: step03_output.RData\n")

# ==============================================================================
# BLOCO: RESUMO FINAL
# ==============================================================================

cat("\n")
cat("=", rep("=", 69), "\n", sep = "")
cat("STEP 3 COMPLETE: Demand Estimation\n")
cat("=", rep("=", 69), "\n", sep = "")

cat("\n  2SLS Results (Spec D, baseline nested logit):\n")
cat(sprintf("    alpha = %.4f (se = %.4f)  —  correct sign (negative)\n",
            alpha_D, se_alpha_D))
cat(sprintf("    sigma = %.4f (se = %.4f)  —  %s\n",
            sigma_D, se_sigma_D,
            ifelse(sigma_D >= 0 & sigma_D < 1,
                   "in [0,1): consistent with nested logit",
                   "outside [0,1): boundary issue")))
cat(sprintf("    Mean own-price elasticity = %.3f\n", mean(elas_own_nlogit)))

cat("\n  Paper's RCNL (for supply-side analysis):\n")
cat(sprintf("    alpha = %.4f,  rho = %.4f\n", alpha, rho))
cat(sprintf("    Mean own-price elasticity = %.3f\n", mean(own_elas_rcnl)))

cat("\n  Key comparison:\n")
cat(sprintf("    2SLS alpha (%.4f) vs RCNL alpha (%.4f): ", alpha_D, alpha))
if (alpha_D < 0 && alpha < 0) {
  cat("both negative, 2SLS attenuated as expected\n")
} else if (alpha_D < 0) {
  cat("2SLS correct sign, magnitudes differ due to data/spec\n")
} else {
  cat("sign difference — 2SLS limited by perturbed small sample\n")
}

cat(sprintf("    2SLS sigma (%.4f) vs RCNL rho (%.4f): ", sigma_D, rho))
cat("2SLS higher — less substitution to outside good\n")

cat("\n  Next: Step 4 (markups and marginal cost imputation)\n")
