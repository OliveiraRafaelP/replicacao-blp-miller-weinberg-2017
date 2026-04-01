# ==============================================================================
# 02c_iv_sensitivity.R — IV sensitivity: baseline vs full instrument set
# Miller & Weinberg (2017) Replication in R
# ==============================================================================
# Analise de sensibilidade da estimacao 2SLS a escolha do conjunto de
# instrumentos. Compara duas especificacoes: (1) baseline com 5 IVs de custo
# e BLP, (2) full com 9 IVs incluindo interacoes com renda. Tambem inclui
# OLS como benchmark viesado e o teste de sobreidentificacao de Sargan.
# Objetivo: verificar se os coeficientes estruturais sao estaveis.
# ==============================================================================

load("C:/Users/rafael.oliveira/Desktop/Problem Set/R_project/data/step02_output.RData")
source("C:/Users/rafael.oliveira/Desktop/Problem Set/R_project/scripts/00_setup.R")

cat("\n")
cat("=", rep("=", 69), "\n", sep = "")
cat("IV SENSITIVITY ANALYSIS: Baseline vs Full Instrument Set\n")
cat("=", rep("=", 69), "\n", sep = "")

# ==============================================================================
# MODEL: Nested logit 2SLS
#   logodds = alpha * price + sigma * logcondshr + prodFE + dateFE + xi
#   Endogenous: price, logcondshr
# ==============================================================================

# ==== BLOCO: Definicao dos conjuntos de instrumentos ====
# Baseline: apenas deslocadores de custo e instrumentos BLP tradicionais.
# Full: adiciona interacoes renda x caracteristicas, que no modelo RCNL
# identificam os coeficientes aleatorios (heterogeneidade de preferencias).

# Especificacao (1): Baseline — custo + BLP (sem renda)
iv_baseline <- c("dist", "coalpost", "sum_dist", "sumdist_abi", "sumdist_mc")

# Especificacao (2): Full — baseline + interacoes com renda
iv_full <- c("dist", "coalpost", "sum_dist", "sumdist_abi", "sumdist_mc",
             "z_inc_const", "z_inc_calor", "z_inc_size", "z_inc_import")

# ==== BLOCO: Construcao das formulas para ivreg ====
# Sintaxe do ivreg: y ~ endogenas + exogenas | exogenas + instrumentos
# Os efeitos fixos de produto e data entram como exogenas em ambos os lados.

make_formula <- function(iv_names) {
  iv_str <- paste(iv_names, collapse = " + ")
  as.formula(paste0(
    "logodds ~ price + logcondshr + factor(prodid) + factor(dateid) | ",
    "factor(prodid) + factor(dateid) + ", iv_str
  ))
}

fml_baseline <- make_formula(iv_baseline)
fml_full     <- make_formula(iv_full)

# ==============================================================================
# A. RUN BOTH SPECIFICATIONS
# ==============================================================================

# ==== BLOCO: Estimacao 2SLS com erros-padrao clusterizados ====
# Erros-padrao clusterizados por cidade (vcovCL) corrigem correlacao
# intra-cluster nos residuos, essencial quando mercados da mesma cidade
# compartilham choques nao observados.

cat("\n--- Running 2SLS regressions ---\n")

fit_baseline <- ivreg(fml_baseline, data = df)
fit_full     <- ivreg(fml_full, data = df)

# Erros-padrao robustos clusterizados por cidade
vcov_baseline <- vcovCL(fit_baseline, cluster = df$cityid)
vcov_full     <- vcovCL(fit_full, cluster = df$cityid)

se_baseline <- sqrt(diag(vcov_baseline))
se_full     <- sqrt(diag(vcov_full))

coef_b <- coef(fit_baseline)
coef_f <- coef(fit_full)

cat("  Baseline (5 IVs): done\n")
cat("  Full     (9 IVs): done\n")

# ==============================================================================
# B. FIRST-STAGE F-STATISTICS
# ==============================================================================

# ==== BLOCO: F do primeiro estagio por especificacao ====
# Compara a forca dos instrumentos entre as duas especificacoes.
# R2 parcial mede a fracao da variacao da endogena explicada pelos IVs
# apos absorver os efeitos fixos.

cat("\n--- First-stage diagnostics ---\n\n")

fe_mat <- model.matrix(~ factor(prodid) + factor(dateid), data = df)[, -1]
y_price <- df$price
y_logcs <- df$logcondshr

compute_fs_F <- function(iv_names, label) {
  Z_mat <- as.matrix(df[, iv_names])
  q <- ncol(Z_mat)
  n <- nrow(df)

  # Primeiro estagio para preco
  fs_r <- lm(y_price ~ fe_mat)
  fs_u <- lm(y_price ~ fe_mat + Z_mat)
  rss_r <- sum(fs_r$residuals^2)
  rss_u <- sum(fs_u$residuals^2)
  k <- ncol(fe_mat) + q
  F_p <- ((rss_r - rss_u) / q) / (rss_u / (n - k))

  # Primeiro estagio para logcondshr
  fs_r2 <- lm(y_logcs ~ fe_mat)
  fs_u2 <- lm(y_logcs ~ fe_mat + Z_mat)
  rss_r2 <- sum(fs_r2$residuals^2)
  rss_u2 <- sum(fs_u2$residuals^2)
  F_s <- ((rss_r2 - rss_u2) / q) / (rss_u2 / (n - k))

  # R2 parcial: quanto da variacao residual os IVs explicam
  pR2_p <- 1 - rss_u / rss_r
  pR2_s <- 1 - rss_u2 / rss_r2

  cat(sprintf("  %s (%d IVs):\n", label, q))
  cat(sprintf("    Price:      F(%d,%d) = %7.2f   partial R2 = %.4f   %s\n",
              q, n - k, F_p, pR2_p,
              ifelse(F_p > 10, "[STRONG]", "[WEAK]")))
  cat(sprintf("    Logcondshr: F(%d,%d) = %7.2f   partial R2 = %.4f   %s\n",
              q, n - k, F_s, pR2_s,
              ifelse(F_s > 10, "[STRONG]", "[WEAK]")))

  list(F_price = F_p, F_logcs = F_s, pR2_price = pR2_p, pR2_logcs = pR2_s)
}

fs_b <- compute_fs_F(iv_baseline, "Baseline")
cat("\n")
fs_f <- compute_fs_F(iv_full, "Full")

# ==============================================================================
# C. COEFFICIENT COMPARISON
# ==============================================================================

# ==== BLOCO: Comparacao de coeficientes entre especificacoes ====
# Se alpha (preco) e sigma (nesting) mudam muito entre baseline e full,
# isso indica sensibilidade a escolha dos instrumentos — sinal de alerta.
# sigma deve estar em [0,1) para modelo bem definido; alpha deve ser < 0.

cat("\n\n--- Coefficient comparison ---\n\n")

params <- c("price", "logcondshr")

cat(sprintf("%-14s | %12s %12s | %12s %12s | %8s\n",
            "Parameter", "Baseline", "(se)", "Full", "(se)", "Delta %"))
cat(rep("-", 80), "\n", sep = "")

for (p in params) {
  b_coef <- coef_b[p]
  b_se   <- se_baseline[p]
  f_coef <- coef_f[p]
  f_se   <- se_full[p]
  delta_pct <- (f_coef - b_coef) / abs(b_coef) * 100

  cat(sprintf("%-14s | %12.6f %12.6f | %12.6f %12.6f | %+7.1f%%\n",
              p, b_coef, b_se, f_coef, f_se, delta_pct))
}

# sigma = coeficiente de logcondshr = parametro de nesting do nested logit
# sigma em [0,1): substitucao dentro do ninho; sigma -> 0: logit simples
cat("\n")
sigma_b <- coef_b["logcondshr"]
sigma_f <- coef_f["logcondshr"]
cat(sprintf("Implied sigma (nesting):  Baseline = %.4f,  Full = %.4f\n", sigma_b, sigma_f))
cat(sprintf("  sigma in [0,1)?  Baseline: %s,  Full: %s\n",
            ifelse(sigma_b >= 0 & sigma_b < 1, "YES", "NO"),
            ifelse(sigma_f >= 0 & sigma_f < 1, "YES", "NO")))

# alpha = coeficiente de preco; deve ser negativo (lei da demanda)
alpha_b <- coef_b["price"]
alpha_f <- coef_f["price"]
cat(sprintf("\nImplied alpha (price):    Baseline = %.6f,  Full = %.6f\n", alpha_b, alpha_f))
cat(sprintf("  alpha < 0 (correct sign)?  Baseline: %s,  Full: %s\n",
            ifelse(alpha_b < 0, "YES", "NO"),
            ifelse(alpha_f < 0, "YES", "NO")))

# Compara com estimativas do artigo original (RCNL completo)
cat(sprintf("\n  Paper's alpha (from dresgmm2): %.6f\n", alpha))
cat(sprintf("  Paper's rho (nesting):         %.6f\n", rho))

# ==============================================================================
# D. HAUSMAN-STYLE COMPARISON
# ==============================================================================

# ==== BLOCO: Teste tipo Hausman ====
# Compara os coeficientes de preco entre as duas especificacoes IV.
# Se os instrumentos adicionais sao validos, ambas as estimativas devem
# convergir. Uma diferenca significativa sugere que algum instrumento
# viola a restricao de exclusao na especificacao mais ampla.

cat("\n\n--- Stability assessment ---\n\n")

d_alpha <- alpha_f - alpha_b
# Variancia da diferenca (sob hipotese de independencia, para ilustracao)
v_d <- vcov_full["price", "price"] - vcov_baseline["price", "price"]
if (v_d > 0) {
  hausman_t <- d_alpha / sqrt(v_d)
  cat(sprintf("Hausman-style test (price coef):\n"))
  cat(sprintf("  Difference:  %.6f\n", d_alpha))
  cat(sprintf("  t-stat:      %.3f\n", hausman_t))
  cat(sprintf("  |t| < 1.96:  %s (no significant difference at 5%%)\n",
              ifelse(abs(hausman_t) < 1.96, "YES", "NO")))
} else {
  cat("Hausman variance negative (finite-sample issue) — using informal comparison\n")
  cat(sprintf("  Absolute difference in alpha: %.6f\n", abs(d_alpha)))
  cat(sprintf("  As %% of baseline alpha:       %.1f%%\n", abs(d_alpha / alpha_b) * 100))
}

# ==============================================================================
# E. OLS BENCHMARK (for comparison — known to be biased)
# ==============================================================================

# ==== BLOCO: Benchmark OLS (viesado) ====
# OLS ignora a endogeneidade do preco. O vies de simultaneidade faz com que
# alpha_OLS seja menos negativo que alpha_IV (vies para cima), pois produtos
# com alta demanda nao observada (xi alto) tem precos mais altos.

cat("\n\n--- OLS benchmark (biased, for reference) ---\n\n")

fit_ols <- lm(logodds ~ price + logcondshr + factor(prodid) + factor(dateid), data = df)
vcov_ols <- vcovCL(fit_ols, cluster = df$cityid)

cat(sprintf("%-14s | %12s %12s\n", "Parameter", "OLS", "(se)"))
cat(rep("-", 42), "\n", sep = "")
for (p in params) {
  cat(sprintf("%-14s | %12.6f %12.6f\n",
              p, coef(fit_ols)[p], sqrt(diag(vcov_ols))[p]))
}

cat("\nExpected bias: OLS alpha should be LESS negative than IV alpha\n")
cat("  (upward bias from simultaneity: high-demand products have high prices)\n")
cat(sprintf("  OLS alpha = %.6f vs IV baseline alpha = %.6f vs IV full alpha = %.6f\n",
            coef(fit_ols)["price"], alpha_b, alpha_f))

# ==============================================================================
# F. OVERIDENTIFICATION TEST (Sargan-Hansen J)
# ==============================================================================

# ==== BLOCO: Teste de sobreidentificacao de Sargan ====
# O teste J de Sargan/Hansen verifica se os instrumentos excluidos satisfazem
# a restricao de exclusao (E[Z'xi] = 0). Sob H0, J ~ chi2(q - k_endog).
# Rejeicao (p < 0.05) indica que ao menos um instrumento e invalido,
# ou que o modelo esta mal especificado.

cat("\n\n--- Overidentification tests ---\n\n")

sargan_b <- summary(fit_baseline, diagnostics = TRUE)$diagnostics
sargan_f <- summary(fit_full, diagnostics = TRUE)$diagnostics

cat("Baseline (5 IVs, 2 endogenous, df = 3):\n")
print(sargan_b)
cat("\nFull (9 IVs, 2 endogenous, df = 7):\n")
print(sargan_f)

# ==============================================================================
# G. SUMMARY TABLE
# ==============================================================================

cat("\n\n")
cat("=", rep("=", 69), "\n", sep = "")
cat("SUMMARY: IV SENSITIVITY ANALYSIS\n")
cat("=", rep("=", 69), "\n", sep = "")

cat("\n")
cat(sprintf("%-30s %15s %15s %15s\n", "", "OLS", "IV Baseline", "IV Full"))
cat(rep("-", 75), "\n", sep = "")
cat(sprintf("%-30s %15s %15s %15s\n", "Instruments", "—",
            "5 (cost+BLP)", "9 (+income)"))
cat(sprintf("%-30s %15.6f %15.6f %15.6f\n", "alpha (price)",
            coef(fit_ols)["price"], alpha_b, alpha_f))
cat(sprintf("%-30s %15.6f %15.6f %15.6f\n", "  (clustered se)",
            sqrt(diag(vcov_ols))["price"], se_baseline["price"], se_full["price"]))
cat(sprintf("%-30s %15.4f %15.4f %15.4f\n", "sigma (logcondshr)",
            coef(fit_ols)["logcondshr"], sigma_b, sigma_f))
cat(sprintf("%-30s %15.4f %15.4f %15.4f\n", "  (clustered se)",
            sqrt(diag(vcov_ols))["logcondshr"],
            se_baseline["logcondshr"], se_full["logcondshr"]))
cat(sprintf("%-30s %15s %15.2f %15.2f\n", "1st-stage F (price)",
            "—", fs_b$F_price, fs_f$F_price))
cat(sprintf("%-30s %15s %15.2f %15.2f\n", "1st-stage F (logcondshr)",
            "—", fs_b$F_logcs, fs_f$F_logcs))
cat(sprintf("%-30s %15s %15.4f %15.4f\n", "partial R2 (price)",
            "—", fs_b$pR2_price, fs_f$pR2_price))
cat(sprintf("%-30s %15d %15d %15d\n", "N",
            nobs(fit_ols), nobs(fit_baseline), nobs(fit_full)))

cat("\n")
cat("Paper reference (RCNL, not nested logit 2SLS):\n")
cat(sprintf("  alpha = %.6f,  rho = %.6f\n", alpha, rho))

cat("\nINTERPRETATION:\n")
diff_pct <- abs((alpha_f - alpha_b) / alpha_b) * 100
if (diff_pct < 10) {
  cat(sprintf("  Price coefficient moves %.1f%% between specifications — STABLE.\n", diff_pct))
} else if (diff_pct < 25) {
  cat(sprintf("  Price coefficient moves %.1f%% — MODERATE sensitivity to income IVs.\n", diff_pct))
} else {
  cat(sprintf("  Price coefficient moves %.1f%% — HIGH sensitivity to income IVs.\n", diff_pct))
}
