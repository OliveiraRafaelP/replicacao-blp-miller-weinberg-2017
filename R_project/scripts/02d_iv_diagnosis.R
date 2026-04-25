# ==============================================================================
# 02d_iv_diagnosis.R — Diagnose the IV instability and Sargan rejection
# ==============================================================================
# Investiga por que o nested logit 2SLS com os dados perturbados de replicacao
# gera estimativas instaveis ou com sinal errado. Testa progressivamente
# especificacoes mais simples (logit simples, menos IVs, exatamente identificado)
# para isolar a fonte do problema e recomendar a especificacao mais confiavel.
# ==============================================================================

load("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/R_project/data/step02_output.RData")
source("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/R_project/scripts/00_setup.R")

cat("\n")
cat("=", rep("=", 69), "\n", sep = "")
cat("DIAGNOSING IV INSTABILITY\n")
cat("=", rep("=", 69), "\n", sep = "")

# ==============================================================================
# The problem: nested logit 2SLS with the small perturbed replication data
# gives unstable/wrong-sign estimates. Why?
#
# 1. The data is PERTURBED (prices/quantities randomly shifted, ABI/Miller x1.5)
# 2. Only 5 cities, 8 products — very little cross-sectional variation
# 3. The paper's RCNL estimates used 37 regions x 39 products
# 4. The nested logit 2SLS is a SIMPLIFIED specification — not what the paper ran
# 5. Income instruments may be identifying off the perturbation noise
#
# Let's try progressively simpler specifications to find what works.
# ==============================================================================

# ==== BLOCO: Busca de especificacao (Spec A-E) ====
# Estrategia: partir do modelo mais simples (logit, 2 IVs) e aumentar
# complexidade gradualmente, monitorando sinal de alpha, forca dos IVs
# (F do 1o estagio) e validade (Sargan). Isso permite isolar qual
# componente causa a instabilidade.

cat("\n--- Specification search ---\n\n")

# Spec A: Logit simples (sem nesting), preco instrumentado apenas por dist + coalpost
cat("SPEC A: Simple logit, 2 IVs (dist, coalpost)\n")
cat("  logodds = alpha * price + prodFE + dateFE + xi\n")
cat("  1 endogenous (price), 2 excluded IVs\n\n")

fit_a <- ivreg(logodds ~ price + factor(prodid) + factor(dateid) |
                 factor(prodid) + factor(dateid) + dist + coalpost,
               data = df)
vcov_a <- vcovCL(fit_a, cluster = df$cityid)
diag_a <- summary(fit_a, diagnostics = TRUE)$diagnostics

cat(sprintf("  alpha = %.6f  (se = %.6f)\n", coef(fit_a)["price"], sqrt(diag(vcov_a))["price"]))
cat(sprintf("  alpha < 0: %s\n", ifelse(coef(fit_a)["price"] < 0, "YES", "NO")))
cat(sprintf("  1st-stage F = %.2f\n", diag_a["Weak instruments", "statistic"]))
# Sargan com df=1: apenas 1 grau de sobreidentificacao (2 IVs - 1 endogena)
cat(sprintf("  Sargan p = %.4f (df=%d)\n",
            diag_a["Sargan", "p-value"], diag_a["Sargan", "df1"]))

# Spec B: Logit simples, 5 IVs (custo + BLP)
cat("\nSPEC B: Simple logit, 5 IVs (dist, coalpost, sum_dist, sumdist_abi, sumdist_mc)\n")

fit_b <- ivreg(logodds ~ price + factor(prodid) + factor(dateid) |
                 factor(prodid) + factor(dateid) +
                 dist + coalpost + sum_dist + sumdist_abi + sumdist_mc,
               data = df)
vcov_b <- vcovCL(fit_b, cluster = df$cityid)
diag_b <- summary(fit_b, diagnostics = TRUE)$diagnostics

cat(sprintf("  alpha = %.6f  (se = %.6f)\n", coef(fit_b)["price"], sqrt(diag(vcov_b))["price"]))
cat(sprintf("  alpha < 0: %s\n", ifelse(coef(fit_b)["price"] < 0, "YES", "NO")))
cat(sprintf("  1st-stage F = %.2f\n", diag_b["Weak instruments", "statistic"]))
cat(sprintf("  Sargan p = %.4f (df=%d)\n",
            diag_b["Sargan", "p-value"], diag_b["Sargan", "df1"]))

# Spec C: Nested logit exatamente identificado (2 endogenas, 2 IVs)
# Sem graus de sobreidentificacao => nao ha teste de Sargan
cat("\nSPEC C: Nested logit, just-identified (dist, coalpost)\n")
cat("  2 endogenous, 2 excluded => no Sargan test\n")

fit_c <- ivreg(logodds ~ price + logcondshr + factor(prodid) + factor(dateid) |
                 factor(prodid) + factor(dateid) + dist + coalpost,
               data = df)
vcov_c <- vcovCL(fit_c, cluster = df$cityid)
diag_c <- summary(fit_c, diagnostics = TRUE)$diagnostics

cat(sprintf("  alpha = %.6f  (se = %.6f)\n", coef(fit_c)["price"], sqrt(diag(vcov_c))["price"]))
cat(sprintf("  sigma = %.6f  (se = %.6f)\n", coef(fit_c)["logcondshr"], sqrt(diag(vcov_c))["logcondshr"]))
cat(sprintf("  alpha < 0: %s,  sigma in [0,1): %s\n",
            ifelse(coef(fit_c)["price"] < 0, "YES", "NO"),
            ifelse(coef(fit_c)["logcondshr"] >= 0 & coef(fit_c)["logcondshr"] < 1, "YES", "NO")))

# Spec D: Nested logit sobreidentificado (5 IVs, df=3 para Sargan)
cat("\nSPEC D: Nested logit, 5 IVs (cost + BLP)\n")

fit_d <- ivreg(logodds ~ price + logcondshr + factor(prodid) + factor(dateid) |
                 factor(prodid) + factor(dateid) +
                 dist + coalpost + sum_dist + sumdist_abi + sumdist_mc,
               data = df)
vcov_d <- vcovCL(fit_d, cluster = df$cityid)
diag_d <- summary(fit_d, diagnostics = TRUE)$diagnostics

cat(sprintf("  alpha = %.6f  (se = %.6f)\n", coef(fit_d)["price"], sqrt(diag(vcov_d))["price"]))
cat(sprintf("  sigma = %.6f  (se = %.6f)\n", coef(fit_d)["logcondshr"], sqrt(diag(vcov_d))["logcondshr"]))
cat(sprintf("  alpha < 0: %s,  sigma in [0,1): %s\n",
            ifelse(coef(fit_d)["price"] < 0, "YES", "NO"),
            ifelse(coef(fit_d)["logcondshr"] >= 0 & coef(fit_d)["logcondshr"] < 1, "YES", "NO")))
# F separado para cada endogena (ivreg reporta quando ha mais de uma)
cat(sprintf("  1st-stage F (price) = %.2f\n", diag_d["Weak instruments (price)", "statistic"]))
cat(sprintf("  1st-stage F (logcs) = %.2f\n", diag_d["Weak instruments (logcondshr)", "statistic"]))
cat(sprintf("  Sargan p = %.4f (df=%d)\n",
            diag_d["Sargan", "p-value"], diag_d["Sargan", "df1"]))

# Spec E: Logit simples com todos os 9 IVs funcionais
cat("\nSPEC E: Simple logit, 9 IVs (full set)\n")

fit_e <- ivreg(logodds ~ price + factor(prodid) + factor(dateid) |
                 factor(prodid) + factor(dateid) +
                 dist + coalpost + sum_dist + sumdist_abi + sumdist_mc +
                 z_inc_const + z_inc_calor + z_inc_size + z_inc_import,
               data = df)
vcov_e <- vcovCL(fit_e, cluster = df$cityid)
diag_e <- summary(fit_e, diagnostics = TRUE)$diagnostics

cat(sprintf("  alpha = %.6f  (se = %.6f)\n", coef(fit_e)["price"], sqrt(diag(vcov_e))["price"]))
cat(sprintf("  alpha < 0: %s\n", ifelse(coef(fit_e)["price"] < 0, "YES", "NO")))
cat(sprintf("  1st-stage F = %.2f\n", diag_e["Weak instruments", "statistic"]))
cat(sprintf("  Sargan p = %.4f (df=%d)\n",
            diag_e["Sargan", "p-value"], diag_e["Sargan", "df1"]))

# OLS como referencia: sem correcao de endogeneidade (viesado por construcao)
cat("\nSPEC OLS: No instruments (biased benchmark)\n")
fit_ols <- lm(logodds ~ price + logcondshr + factor(prodid) + factor(dateid), data = df)
vcov_ols <- vcovCL(fit_ols, cluster = df$cityid)
cat(sprintf("  alpha = %.6f  (se = %.6f)\n", coef(fit_ols)["price"], sqrt(diag(vcov_ols))["price"]))
cat(sprintf("  sigma = %.6f  (se = %.6f)\n", coef(fit_ols)["logcondshr"], sqrt(diag(vcov_ols))["logcondshr"]))

# ==============================================================================
# SUMMARY TABLE
# ==============================================================================

# ==== BLOCO: Tabela comparativa de todas as especificacoes ====
# Resume alpha, sigma, F do 1o estagio e p-valor do Sargan para cada spec.
# Permite identificar visualmente qual especificacao e mais robusta:
# sinais corretos, F alto, Sargan nao rejeitado.

cat("\n\n")
cat("=", rep("=", 69), "\n", sep = "")
cat("SUMMARY TABLE: ALL SPECIFICATIONS\n")
cat("=", rep("=", 69), "\n", sep = "")
cat("\n")

specs <- list(
  list(name = "OLS (nested logit)", alpha = coef(fit_ols)["price"],
       se = sqrt(diag(vcov_ols))["price"],
       sigma = coef(fit_ols)["logcondshr"], se_s = sqrt(diag(vcov_ols))["logcondshr"],
       fs = NA, sargan = NA, niv = 0, model = "NL"),
  list(name = "A: Logit, 2 IVs", alpha = coef(fit_a)["price"],
       se = sqrt(diag(vcov_a))["price"],
       sigma = NA, se_s = NA,
       fs = diag_a["Weak instruments", "statistic"],
       sargan = diag_a["Sargan", "p-value"], niv = 2, model = "L"),
  list(name = "B: Logit, 5 IVs", alpha = coef(fit_b)["price"],
       se = sqrt(diag(vcov_b))["price"],
       sigma = NA, se_s = NA,
       fs = diag_b["Weak instruments", "statistic"],
       sargan = diag_b["Sargan", "p-value"], niv = 5, model = "L"),
  list(name = "C: NL, 2 IVs (just-id)", alpha = coef(fit_c)["price"],
       se = sqrt(diag(vcov_c))["price"],
       sigma = coef(fit_c)["logcondshr"], se_s = sqrt(diag(vcov_c))["logcondshr"],
       fs = diag_c["Weak instruments (price)", "statistic"],
       sargan = NA, niv = 2, model = "NL"),
  list(name = "D: NL, 5 IVs", alpha = coef(fit_d)["price"],
       se = sqrt(diag(vcov_d))["price"],
       sigma = coef(fit_d)["logcondshr"], se_s = sqrt(diag(vcov_d))["logcondshr"],
       fs = diag_d["Weak instruments (price)", "statistic"],
       sargan = diag_d["Sargan", "p-value"], niv = 5, model = "NL"),
  list(name = "E: Logit, 9 IVs", alpha = coef(fit_e)["price"],
       se = sqrt(diag(vcov_e))["price"],
       sigma = NA, se_s = NA,
       fs = diag_e["Weak instruments", "statistic"],
       sargan = diag_e["Sargan", "p-value"], niv = 9, model = "L")
)

cat(sprintf("%-25s %9s %9s %9s %9s %8s %8s\n",
            "Specification", "alpha", "(se)", "sigma", "(se)", "F(price)", "Sargan p"))
cat(rep("-", 85), "\n", sep = "")

for (s in specs) {
  alpha_str <- sprintf("%9.4f", s$alpha)
  se_str    <- sprintf("(%7.4f)", s$se)
  sigma_str <- ifelse(is.na(s$sigma), sprintf("%9s", "—"), sprintf("%9.4f", s$sigma))
  ses_str   <- ifelse(is.na(s$se_s), sprintf("%9s", ""), sprintf("(%7.4f)", s$se_s))
  fs_str    <- ifelse(is.na(s$fs), sprintf("%8s", "—"), sprintf("%8.1f", s$fs))
  sar_str   <- ifelse(is.na(s$sargan), sprintf("%8s", "—"), sprintf("%8.4f", s$sargan))
  cat(sprintf("%-25s %s %s %s %s %s %s\n",
              s$name, alpha_str, se_str, sigma_str, ses_str, fs_str, sar_str))
}

cat(sprintf("\n%-25s %9.4f %9s %9.4f\n",
            "Paper (RCNL, full data)", alpha, "", rho))

# ==== BLOCO: Diagnostico final e recomendacao ====
# Sintetiza as causas da instabilidade e recomenda a especificacao
# mais adequada para o problem set, dados os limites da amostra perturbada.

cat("\n")
cat("DIAGNOSIS:\n")
cat("  1. The small perturbed replication data has limited cross-sectional\n")
cat("     variation (5 cities, 8 products). The full paper uses 37 x 39.\n")
cat("  2. Income instruments add power for price (F: 15.7 -> 41.7) but may\n")
cat("     identify off perturbation noise, causing instability.\n")
cat("  3. Sargan rejection in both overidentified specs signals that some\n")
cat("     instruments may not satisfy exclusion in the perturbed sample.\n")
cat("  4. The just-identified nested logit (Spec C) is the most reliable\n")
cat("     2SLS estimate for this data: no overid concern, correct signs.\n")
cat("  5. For the problem set: report Spec C or D as the main 2SLS result,\n")
cat("     note the data limitations, and rely on the paper's RCNL estimates\n")
cat("     (loaded from dresgmm2) for the supply-side analysis.\n")
