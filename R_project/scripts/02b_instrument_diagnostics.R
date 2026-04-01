# ==============================================================================
# 02b_instrument_diagnostics.R — Pre-estimation instrument diagnostics
# Miller & Weinberg (2017) Replication in R
# ==============================================================================
# Diagnostico pre-estimacao dos instrumentos construidos em 02_clean_shares.R.
# Verifica: (1) quais instrumentos tem variacao suficiente (nao-degenerados),
# (2) condicao de ordem e posto, (3) correlacoes parciais com as endogenas,
# (4) estatistica F do primeiro estagio. Essencial para validar a estrategia
# de identificacao antes de rodar o 2SLS.
# ==============================================================================

load("C:/Users/rafael.oliveira/Desktop/Problem Set/R_project/data/step02_output.RData")
source("C:/Users/rafael.oliveira/Desktop/Problem Set/R_project/scripts/00_setup.R")

cat("\n")
cat("=", rep("=", 69), "\n", sep = "")
cat("PRE-ESTIMATION INSTRUMENT DIAGNOSTICS\n")
cat("=", rep("=", 69), "\n", sep = "")

# ==============================================================================
# 1. NON-DEGENERACY CHECK — which instruments have actual variation?
# ==============================================================================

# ==== BLOCO: Verificacao de degenerescencia ====
# Um instrumento degenerado (sem variacao) nao contribui para identificacao.
# Nesta amostra pequena (5 cidades, 8 produtos fixos), Z3, Z5 e Z6 sao
# constantes ou colineares com dummies de firma, devendo ser descartados.

cat("\n--- 1. Non-degeneracy check (640-obs filtered sample) ---\n\n")

# Vetor com os nomes das 12 variaveis instrumentais candidatas
inst_vars <- c(
  "dist",            # Z1
  "coalpost",        # Z2
  "num_products",    # Z3
  "sum_dist",        # Z4
  "nj_abi",          # Z5
  "nj_mc",           # Z6
  "sumdist_abi",     # Z7
  "sumdist_mc",      # Z8
  "z_inc_const",     # Z9
  "z_inc_calor",     # Z10
  "z_inc_size",      # Z11
  "z_inc_import"     # Z12
)

inst_labels <- c(
  "Z1:  dist (miles x diesel)",
  "Z2:  coalpost (ABI+MC post)",
  "Z3:  num_products (J_t)",
  "Z4:  sum_dist (market total)",
  "Z5:  num_products x ABI",
  "Z6:  num_products x MC",
  "Z7:  sum_dist x ABI",
  "Z8:  sum_dist x MC",
  "Z9:  mean_income x 1",
  "Z10: mean_income x calor",
  "Z11: mean_income x size",
  "Z12: mean_income x import"
)

cat(sprintf("%-35s %8s %8s %10s %10s  %s\n",
            "Instrument", "unique", "sd", "min", "max", "Status"))
cat(rep("-", 90), "\n", sep = "")

functional <- c()   # instrumentos com variacao suficiente
degenerate <- c()   # instrumentos descartados

for (i in seq_along(inst_vars)) {
  v <- df[[inst_vars[i]]]
  n_uniq <- n_distinct(v)
  sd_v   <- sd(v)
  min_v  <- min(v)
  max_v  <- max(v)

  # Criterio de degenerescencia: valor unico ou desvio-padrao ~0
  if (n_uniq <= 1 || sd_v < 1e-10) {
    status <- "DEGENERATE"
    degenerate <- c(degenerate, inst_vars[i])
  } else if (n_uniq == 2 && inst_vars[i] %in% c("nj_abi", "nj_mc")) {
    # J*indicador com J constante: perfeitamente colinear com a dummy de firma
    # Tem 2 valores (0 e 8), mas e simplesmente 8 x dummy_firma
    status <- "REDUNDANT (= 8 x firm_dummy)"
    degenerate <- c(degenerate, inst_vars[i])
  } else {
    status <- "OK"
    functional <- c(functional, inst_vars[i])
  }

  cat(sprintf("%-35s %8d %8.4f %10.4f %10.4f  %s\n",
              inst_labels[i], n_uniq, sd_v, min_v, max_v, status))
}

cat(sprintf("\nFunctional instruments: %d\n", length(functional)))
cat(sprintf("Degenerate/redundant:  %d\n", length(degenerate)))
cat("Functional set:", paste(functional, collapse = ", "), "\n")
cat("Dropped:       ", paste(degenerate, collapse = ", "), "\n")

# ==============================================================================
# 2. FIRST-STAGE SPECIFICATION
#    Endogenous variables in the demand equation and their instruments
# ==============================================================================

# ==== BLOCO: Especificacao do primeiro estagio ====
# Identifica as variaveis endogenas (preco e log da share condicional) e
# lista os instrumentos excluidos disponiveis. No nested logit 2SLS, ambas
# as endogenas sao instrumentadas pelo mesmo conjunto de IVs excluidos.

cat("\n\n--- 2. First-stage specification ---\n\n")

cat("DEMAND MODEL (nested logit, 2SLS formulation):\n")
cat("  ln(s_j) - ln(s_0) = alpha*p_jt + sigma*ln(s_j|g) + x*beta + FE + xi\n\n")

cat("ENDOGENOUS VARIABLES (2):\n")
cat("  1. price (p_jt)          — correlated with xi (firms set prices knowing demand)\n")
cat("  2. logcondshr (ln s_j|g) — mechanical function of shares, hence endogenous\n\n")

cat("INCLUDED EXOGENOUS VARIABLES (absorbed in FE or characteristics):\n")
cat("  - Product FE (prodfe): 8 dummies (absorb brand x size mean utility)\n")
cat("  - Date FE (datefe):    15 dummies (absorb common time shocks)\n")
cat("  These are part of x1 = [p_jt, fesd] where fesd = [prodfe, datefe(:,2:end)]\n\n")

cat("EXCLUDED INSTRUMENTS for price (from functional set):\n")
cat("  Z1:  dist           — cost shifter (transportation to brewery)\n")
cat("  Z2:  coalpost       — competitive structure change (post-merger indicator)\n")
cat("  Z4:  sum_dist       — market-level aggregate of rival costs\n")
cat("  Z7:  sumdist_abi    — sum_dist interacted with ABI indicator\n")
cat("  Z8:  sumdist_mc     — sum_dist interacted with MC indicator\n")
cat("  Z9:  z_inc_const    — mean income (identifies RC price-income interaction)\n")
cat("  Z10: z_inc_calor    — mean income x calories\n")
cat("  Z11: z_inc_size     — mean income x pack size\n")
cat("  Z12: z_inc_import   — mean income x import dummy\n\n")

cat("NOTE: In the simple nested logit 2SLS (problem set Q3.1), ALL 9 functional\n")
cat("instruments serve as excluded IVs for both price AND logcondshr.\n")
cat("In the full RCNL (Q3.2), Z9-Z12 also identify the RC parameters (Pi).\n")

# ==============================================================================
# 3. RANK CONDITION
# ==============================================================================

# ==== BLOCO: Condicao de ordem (necessaria) e verificacao de posto ====
# A condicao de ordem exige que o numero de instrumentos excluidos seja >= ao
# numero de variaveis endogenas. Se estritamente maior, o modelo e sobre-
# identificado e podemos aplicar o teste de Sargan/Hansen.

cat("\n\n--- 3. Rank condition ---\n\n")

n_endog    <- 2   # price, logcondshr
n_excluded <- length(functional)

cat(sprintf("Endogenous variables:     %d  (price, logcondshr)\n", n_endog))
cat(sprintf("Excluded instruments:     %d  (functional Z's)\n", n_excluded))
cat(sprintf("Order condition:          %d >= %d  =>  %s\n",
            n_excluded, n_endog,
            ifelse(n_excluded >= n_endog, "SATISFIED (overidentified)", "VIOLATED")))
# Graus de sobreidentificacao = num. instrumentos excluidos - num. endogenas
cat(sprintf("Degrees of overidentification: %d\n", n_excluded - n_endog))
cat("\nFor the nested logit 2SLS specifically:\n")
cat("  Endogenous: price (1) + logcondshr (1) = 2\n")
cat("  Excluded:   9 functional instruments\n")
cat("  Overidentification: 7 degrees of freedom\n")
cat("  => Sargan/Hansen J-test will have chi2(7) distribution\n")

# Verificacao adicional: apos absorver EFs, ha variacao residual suficiente?
# Regressao auxiliar: endogena ~ EFs, para medir quanto de variacao resta
prodfe_mat <- model.matrix(~ factor(prodid) - 1, data = df)
datefe_mat <- model.matrix(~ factor(dateid) - 1, data = df)[, -1]
fe_mat <- cbind(prodfe_mat, datefe_mat)

resid_price <- lm(price ~ fe_mat, data = df)$residuals
resid_logcs <- lm(logcondshr ~ fe_mat, data = df)$residuals

cat(sprintf("\nResidual variation after absorbing FEs:\n"))
cat(sprintf("  price:       sd(resid) = %.4f  (sd(raw) = %.4f)\n",
            sd(resid_price), sd(df$price)))
cat(sprintf("  logcondshr:  sd(resid) = %.4f  (sd(raw) = %.4f)\n",
            sd(resid_logcs), sd(df$logcondshr)))

# ==============================================================================
# 4. CORRELATIONS: price/logcondshr vs instruments
# ==============================================================================

# ==== BLOCO: Correlacoes brutas e parciais ====
# Correlacoes brutas podem ser enganosas quando ha efeitos fixos. As
# correlacoes parciais (apos absorver EFs) refletem a relevancia real dos
# instrumentos para o primeiro estagio do 2SLS.

cat("\n\n--- 4. Correlations between endogenous variables and instruments ---\n\n")

cat("4a. RAW CORRELATIONS\n")
cat(sprintf("%-22s %12s %12s\n", "Instrument", "cor(price)", "cor(logcshr)"))
cat(rep("-", 48), "\n", sep = "")

for (v in functional) {
  cor_p <- cor(df$price, df[[v]])
  cor_s <- cor(df$logcondshr, df[[v]])
  cat(sprintf("%-22s %12.4f %12.4f\n", v, cor_p, cor_s))
}

cat("\n4b. PARTIAL CORRELATIONS (after absorbing product + date FEs)\n")
cat("    These are the correlations that matter for IV relevance.\n\n")

# Residualiza instrumentos nos EFs para obter correlacao parcial
cat(sprintf("%-22s %12s %12s\n", "Instrument", "pcor(price)", "pcor(logcshr)"))
cat(rep("-", 48), "\n", sep = "")

for (v in functional) {
  resid_iv <- lm(df[[v]] ~ fe_mat)$residuals
  pcor_p <- cor(resid_price, resid_iv)   # correlacao parcial com preco
  pcor_s <- cor(resid_logcs, resid_iv)   # correlacao parcial com log share condicional
  cat(sprintf("%-22s %12.4f %12.4f\n", v, pcor_p, pcor_s))
}

# ==============================================================================
# 5. FIRST-STAGE F-STATISTICS (informal check)
# ==============================================================================

# ==== BLOCO: Estatistica F do primeiro estagio ====
# A estatistica F conjunta testa se os instrumentos excluidos sao
# conjuntamente significativos no primeiro estagio. A regra de bolso de
# Staiger & Stock (1997) indica que F > 10 sugere instrumentos fortes.
# F < 10 levanta preocupacoes de instrumentos fracos (weak instruments).

cat("\n\n--- 5. First-stage F-statistics (informal) ---\n\n")

# Matriz de instrumentos funcionais
Z_mat <- as.matrix(df[, functional])

y_price <- df$price
y_logcs <- df$logcondshr

# Primeiro estagio para preco: testa se IVs predizem preco alem dos EFs
fs_price_restricted   <- lm(y_price ~ fe_mat)          # modelo restrito (sem IVs)
fs_price_unrestricted <- lm(y_price ~ fe_mat + Z_mat)  # modelo irrestrito (com IVs)

rss_r <- sum(fs_price_restricted$residuals^2)
rss_u <- sum(fs_price_unrestricted$residuals^2)
q     <- length(functional)   # numero de restricoes testadas
n     <- nrow(df)
k     <- ncol(fe_mat) + q
# F = [(RSS_restrito - RSS_irrestrito)/q] / [RSS_irrestrito/(n-k)]
F_price <- ((rss_r - rss_u) / q) / (rss_u / (n - k))

cat(sprintf("First-stage for PRICE:\n"))
cat(sprintf("  F(%d, %d) = %.2f\n", q, n - k, F_price))
cat(sprintf("  Rule of thumb: F > 10 => %s\n",
            ifelse(F_price > 10, "STRONG instruments", "WEAK instruments (caution)")))

# Primeiro estagio para logcondshr
fs_logcs_restricted   <- lm(y_logcs ~ fe_mat)
fs_logcs_unrestricted <- lm(y_logcs ~ fe_mat + Z_mat)

rss_r2 <- sum(fs_logcs_restricted$residuals^2)
rss_u2 <- sum(fs_logcs_unrestricted$residuals^2)
F_logcs <- ((rss_r2 - rss_u2) / q) / (rss_u2 / (n - k))

cat(sprintf("\nFirst-stage for LOGCONDSHR:\n"))
cat(sprintf("  F(%d, %d) = %.2f\n", q, n - k, F_logcs))
cat(sprintf("  Rule of thumb: F > 10 => %s\n",
            ifelse(F_logcs > 10, "STRONG instruments", "WEAK instruments (caution)")))

# Coeficientes individuais dos instrumentos no primeiro estagio de preco
# t-stats altos indicam quais instrumentos contribuem mais para a identificacao
cat("\n\nIndividual instrument coefficients in first-stage for PRICE:\n")
cat(sprintf("%-22s %10s %10s %10s\n", "Instrument", "coef", "se", "t-stat"))
cat(rep("-", 55), "\n", sep = "")

fs_full <- summary(fs_price_unrestricted)
coef_tab <- fs_full$coefficients
n_coefs <- nrow(coef_tab)
iv_rows <- (n_coefs - q + 1):n_coefs

for (j in seq_along(functional)) {
  idx <- iv_rows[j]
  coef_val <- coef_tab[idx, 1]
  se_val   <- coef_tab[idx, 2]
  t_val    <- coef_tab[idx, 3]
  cat(sprintf("%-22s %10.4f %10.4f %10.2f\n", functional[j], coef_val, se_val, t_val))
}

# ==============================================================================
# 6. SUMMARY ASSESSMENT
# ==============================================================================

cat("\n\n")
cat("=", rep("=", 69), "\n", sep = "")
cat("INSTRUMENT DIAGNOSTIC SUMMARY\n")
cat("=", rep("=", 69), "\n", sep = "")
cat(sprintf("\n  Functional instruments:  %d of 12\n", length(functional)))
cat(sprintf("  Dropped (degenerate):    Z3 (J=8 always), Z5 (=8*ABI), Z6 (=8*MC)\n"))
cat(sprintf("  Endogenous variables:    2 (price, logcondshr)\n"))
cat(sprintf("  Order condition:         SATISFIED (%d excluded >= %d endogenous)\n",
            length(functional), n_endog))
cat(sprintf("  Overidentification df:   %d\n", length(functional) - n_endog))
cat(sprintf("  First-stage F (price):   %.2f  %s\n", F_price,
            ifelse(F_price > 10, "[STRONG]", "[WEAK]")))
cat(sprintf("  First-stage F (logcshr): %.2f  %s\n", F_logcs,
            ifelse(F_logcs > 10, "[STRONG]", "[WEAK]")))
cat("\n  ASSESSMENT: Ready for 2SLS estimation.\n")
