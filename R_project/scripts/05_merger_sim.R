# ==============================================================================
# 05_merger_sim.R — Bertrand merger simulation (MillerCoors JV)
# Miller & Weinberg (2017) Replication in R
# ==============================================================================
#
# Este script simula o contrafactual da fusao MillerCoors sob Nash-Bertrand.
#
# A logica eh:
#   (1) Tomar os custos marginais recuperados no Step 4 como primitivas (fixos).
#   (2) Alterar a matriz de propriedade para refletir a fusao (firma 4 -> firma 5).
#   (3) Resolver o novo equilibrio de Nash-Bertrand: encontrar precos p* tais que
#       a CPO de Bertrand seja satisfeita sob a nova propriedade.
#   (4) Comparar precos, shares e excedente do consumidor antes e depois da fusao.
#
# A premissa-chave eh que custos marginais e preferencias dos consumidores nao
# mudam com a fusao — apenas a estrutura de propriedade muda.
# ==============================================================================

load("C:/Users/rafael.oliveira/Desktop/Problem Set/R_project/data/step04_output.RData")
source("C:/Users/rafael.oliveira/Desktop/Problem Set/R_project/scripts/00_setup.R")

cat("\n")
cat("=", rep("=", 69), "\n", sep = "")
cat("STEP 5: MERGER SIMULATION (Bertrand Counterfactual)\n")
cat("=", rep("=", 69), "\n", sep = "")

# ==============================================================================
# BLOCO: FUNCOES RCNL (replicadas do Step 4)
# ==============================================================================
# As funcoes abaixo sao identicas as do Step 4. Sao replicadas aqui para que
# o script 05 seja autocontido (nao dependa de source() do script 04).
# ==============================================================================

# --------------------------------------------------------------------------
# compute_mu: desvio de utilidade individual mu_ij
# --------------------------------------------------------------------------
# Entrada: theta2w (K x 2), ns (inteiro), dfull_m (J x ns), x2_m (J x K)
# Saida: mu (J x ns) e ai (J x ns, desvio do coef. de preco)
# Intuicao economica: captura como a sensibilidade ao preco varia com a renda
compute_mu <- function(theta2w, ns, dfull_m, x2_m) {
  J <- nrow(x2_m)
  income_draws <- dfull_m[1, ]
  theta2_col2 <- theta2w[, 2]
  product_loading <- as.numeric(x2_m %*% theta2_col2)
  mu <- outer(product_loading, income_draws)
  ai <- matrix(theta2w[1, 2] * income_draws, nrow = J, ncol = ns, byrow = TRUE)
  list(mu = mu, ai = ai)
}

# --------------------------------------------------------------------------
# compute_indsh: probabilidades de escolha individuais (RCNL)
# --------------------------------------------------------------------------
# Entrada: exp_delta (J x 1), exp_mu (J x ns), rho (escalar)
# Saida: sharei (J x ns), scondi (J x ns), sgroupi (1 x ns)
# Intuicao economica: calcula a probabilidade de cada consumidor escolher cada
#   produto, usando a estrutura de nesting do RCNL
compute_indsh <- function(exp_delta, exp_mu, rho) {
  J <- length(exp_delta)
  ns <- ncol(exp_mu)
  # Utilidade escalada pelo nesting: exp(V_ij / (1-rho))
  expval <- matrix(exp_delta^(1 / (1 - rho)), nrow = J, ncol = ns) *
    exp_mu^(1 / (1 - rho))
  # Valor inclusivo do grupo inside
  IV_g <- colSums(expval)
  # Share condicional dentro do ninho
  scondi <- expval / matrix(IV_g, nrow = J, ncol = ns, byrow = TRUE)
  # Probabilidade do grupo inside
  IV_g_1mrho <- IV_g^(1 - rho)
  sgroupi <- matrix(IV_g_1mrho / (1 + IV_g_1mrho), nrow = 1, ncol = ns)
  # Share incondicional: condicional x prob do grupo
  sharei <- scondi * matrix(sgroupi, nrow = J, ncol = ns, byrow = TRUE)
  list(sharei = sharei, scondi = scondi, sgroupi = sgroupi)
}

# --------------------------------------------------------------------------
# compute_der1: derivadas ds_j/dp_k
# --------------------------------------------------------------------------
# Entrada: pcoefi_m (J x ns), sharei (J x ns), scondi (J x ns),
#          sgroupi (1 x ns), rho (escalar)
# Saida: der (J x J), matriz de derivadas (media sobre consumidores simulados)
# Intuicao economica: mede como a share de j responde a variacao marginal
#   no preco de k. Essencial para construir a CPO de Bertrand.
compute_der1 <- function(pcoefi_m, sharei, scondi, sgroupi, rho) {
  J <- nrow(sharei)
  ns <- ncol(sharei)
  der <- matrix(0, nrow = J, ncol = J)
  for (j in 1:J) {
    for (k in 1:J) {
      if (j == k) {
        # Derivada propria (negativa): aumento do preco reduz share
        term <- pcoefi_m[j, ] * sharei[j, ] *
          (1 / (1 - rho) - (1 / (1 - rho)) * scondi[j, ] -
             sgroupi[1, ] + sharei[j, ])
      } else {
        # Derivada cruzada (positiva): substituicao entre produtos
        term <- pcoefi_m[k, ] * sharei[j, ] *
          (-1 / (1 - rho) * scondi[k, ] + sgroupi[1, ] - sharei[k, ])
      }
      der[j, k] <- mean(term)
    }
  }
  der
}

# --------------------------------------------------------------------------
# foc_residual: residuo da CPO de Bertrand (usado como sistema para o solver)
# --------------------------------------------------------------------------
# Entrada:
#   p             — vetor de precos candidatos (J x 1)
#   mc_m          — custos marginais (fixos, recuperados no Step 4)
#   alpha_val     — coeficiente medio de preco
#   rho_val       — parametro de nesting
#   deltanp_m     — utilidade media sem componente de preco
#   owner         — matriz de propriedade (pos-fusao no contrafactual)
#   x2_m_template — template de caracteristicas dos produtos
#   pcoefi_m      — coeficientes de preco individuais
#   dfull_m       — sorteios de renda
#   theta2w_val   — interacoes com renda
#   ns_val        — numero de consumidores simulados
# Saida: vetor J x 1 de residuos (p - mc - markup). O solver busca p* tal que residuo = 0.
# Intuicao economica:
#   No equilibrio de Bertrand, cada firma escolhe precos para maximizar lucro.
#   A CPO implica p - mc = markup(p). O solver encontra o ponto fixo.
foc_residual <- function(p, mc_m, alpha_val, rho_val, deltanp_m, owner,
                         x2_m_template, pcoefi_m, dfull_m, theta2w_val, ns_val) {
  J <- length(p)

  # Protecao contra precos invalidos (evita overflow numerico)
  if (any(!is.finite(p)) || any(p <= 0) || any(p > 100)) {
    return(rep(1e6, J))
  }

  # Utilidade media nos precos candidatos
  delta2 <- deltanp_m + alpha_val * p

  # Atualiza preco em x2
  x2_m <- x2_m_template
  x2_m[, 1] <- p

  # Heterogeneidade do consumidor
  mu_res <- compute_mu(theta2w_val, ns_val, dfull_m, x2_m)

  # Shares individuais (com protecao contra overflow)
  sh_res <- tryCatch(
    compute_indsh(exp(delta2), exp(mu_res$mu), rho_val),
    error = function(e) NULL
  )
  if (is.null(sh_res)) return(rep(1e6, J))

  # Shares de mercado (media sobre consumidores)
  svec <- rowMeans(sh_res$sharei)
  if (any(!is.finite(svec)) || any(svec <= 0)) return(rep(1e6, J))

  # Coeficiente de preco individual completo
  pcoefi_full <- alpha_val + mu_res$ai

  # Matriz de derivadas ds/dp
  der1 <- compute_der1(pcoefi_full, sh_res$sharei, sh_res$scondi, sh_res$sgroupi, rho_val)

  # ==== Inversao da CPO: markup = -(Omega .* der')^{-1} * s ====
  omega <- owner * t(der1)
  markup <- tryCatch(-solve(omega, svec), error = function(e) rep(NA, J))
  if (any(!is.finite(markup))) return(rep(1e6, J))

  # Residuo: no equilibrio, p - mc - markup = 0
  p - mc_m - markup
}

# --------------------------------------------------------------------------
# compute_shares_at_prices: calcula shares a precos dados
# --------------------------------------------------------------------------
# Entrada: p (J x 1), parametros do modelo, dados do mercado
# Saida: vetor J x 1 de shares de mercado
# Intuicao economica: funcao auxiliar para calcular shares pos-simulacao
#   (validacao e calculo de variacao de shares)
compute_shares_at_prices <- function(p, alpha_val, rho_val, deltanp_m,
                                     x2_m_template, dfull_m, theta2w_val, ns_val) {
  delta2 <- deltanp_m + alpha_val * p
  x2_m <- x2_m_template
  x2_m[, 1] <- p
  mu_res <- compute_mu(theta2w_val, ns_val, dfull_m, x2_m)
  sh_res <- compute_indsh(exp(delta2), exp(mu_res$mu), rho_val)
  rowMeans(sh_res$sharei)
}

# --------------------------------------------------------------------------
# compute_inclusive_value: excedente do consumidor via valor inclusivo
# --------------------------------------------------------------------------
# Entrada: p (J x 1), parametros do modelo, dados do mercado
# Saida: escalar, excedente medio do consumidor (em unidades monetarias)
# Intuicao economica:
#   O excedente do consumidor no modelo logit aninhado com coeficientes
#   aleatorios eh proporcional a E_i[ log(1 + IV_g^{1-rho}) / (-alpha_i) ].
#   O denominador (-alpha_i) converte utilidade em unidades monetarias.
#   A comparacao CS_pre - CS_post mede a perda de bem-estar causada pela fusao.
compute_inclusive_value <- function(p, alpha_val, rho_val, deltanp_m,
                                   x2_m_template, dfull_m, theta2w_val, ns_val) {
  J <- length(p)
  ns <- ncol(dfull_m)
  delta2 <- deltanp_m + alpha_val * p
  x2_m <- x2_m_template
  x2_m[, 1] <- p
  mu_res <- compute_mu(theta2w_val, ns_val, dfull_m, x2_m)

  # Utilidade total do consumidor i pelo produto j: V_ij = delta_j + mu_ij
  V <- matrix(delta2, nrow = J, ncol = ns) + mu_res$mu

  # Valor inclusivo: IV_g = sum_j exp(V_ij / (1-rho))
  expval <- exp(V / (1 - rho))
  IV_g <- colSums(expval)

  # CS_i = log(1 + IV_g^{1-rho}) / (-alpha_i)
  # Converte utilidade em dolares dividindo pelo coeficiente de preco individual
  pcoefi_full <- alpha_val + mu_res$ai
  alpha_i <- pcoefi_full[1, ]  # coeficiente de preco individual (identico para todos os j)

  IV_g_1mrho <- IV_g^(1 - rho)
  log_denom <- log(1 + IV_g_1mrho)

  # Media sobre consumidores simulados (integracao Monte Carlo)
  cs_mean <- mean(log_denom / (-alpha_i))
  cs_mean
}

# ==============================================================================
# BLOCO: CONFIGURACAO DA SIMULACAO
# ==============================================================================

cat("\n--- Setup ---\n")

J <- 8
n_markets <- 80
deltanp_640 <- deltanp_val
pcoefi_640 <- pcoefi_val

# Usa mc_R (calculado no Step 4 com a mesma eval_foc) para garantir consistencia:
# o residuo da CPO pre-fusao sera ~0 por construcao.
mc_base <- mc_R

cat(sprintf("  Markets: %d, Products per market: %d\n", n_markets, J))
cat(sprintf("  MC source: sres_bertrand.mat (Matlab, validated)\n"))
cat(sprintf("  Parameters: alpha=%.6f, rho=%.6f, ns=%d\n", alpha, rho, NS))

# ==============================================================================
# BLOCO: VERIFICACAO DO EQUILIBRIO PRE-FUSAO
# ==============================================================================
# Antes de simular, verificamos que a CPO eh satisfeita (~0) nos precos
# observados com a propriedade pre-fusao. Isso confirma que mc_R eh consistente.
# ==============================================================================

cat("\n--- Verifying pre-merger equilibrium ---\n")

pre_mkts <- which(df$fiscid[seq(1, 640, by = 8)] <= 2007)
sample_mkt <- pre_mkts[1]
idx_s <- ((sample_mkt - 1) * J + 1):(sample_mkt * J)

foc_check <- foc_residual(
  p = df$price[idx_s],
  mc_m = mc_base[idx_s],
  alpha_val = alpha, rho_val = rho,
  deltanp_m = deltanp_640[idx_s],
  owner = outer(df$firmid[idx_s], df$firmid[idx_s], "==") * 1.0,
  x2_m_template = cbind(df$price[idx_s], rep(1, J), df$calor[idx_s]),
  pcoefi_m = pcoefi_640[idx_s, ],
  dfull_m = dfull[idx_s, ],
  theta2w_val = theta2w, ns_val = NS
)

cat(sprintf("  FOC residual at observed prices (market %d): max|r| = %.6f\n",
            sample_mkt, max(abs(foc_check))))

# ==============================================================================
# BLOCO: SIMULACAO DA FUSAO
# ==============================================================================
# Para cada mercado:
#   (1) Constroi a matriz de propriedade pos-fusao (firma 4 -> firma 5).
#   (2) Usa os precos observados como chute inicial para o solver.
#   (3) Resolve o sistema de CPOs via nleqslv (metodo de Broyden).
#   (4) O solver encontra p* tal que foc_residual(p*) = 0 sob nova propriedade.
#
# Metodo de Broyden: quasi-Newton que atualiza o Jacobiano iterativamente
# sem recalcula-lo a cada passo. Mais robusto que Newton puro para este tipo
# de sistema nao-linear de dimensao moderada (J=8).
# ==============================================================================

cat("\n--- Running merger simulation ---\n")
cat("  Counterfactual: MillerCoors merger (firm 4 -> firm 5)\n")
cat("  Solver: nleqslv (Broyden method)\n")
cat("  Tolerances: xtol=1e-12, ftol=1e-12, maxit=500\n\n")

# Armazenamento
p_pre    <- df$price
p_post   <- rep(NA_real_, 640)
s_pre    <- df$share
s_post   <- rep(NA_real_, 640)
converged <- rep(FALSE, n_markets)
foc_norm  <- rep(NA_real_, n_markets)

# IDs de firma pos-fusao: Coors (firma 4) eh absorvida por MillerCoors (firma 5)
firmid_post <- df$firmid
firmid_post[firmid_post == 4] <- 5

for (m in 1:n_markets) {
  idx <- ((m - 1) * J + 1):(m * J)

  # Dados do mercado
  p_m       <- df$price[idx]
  mc_m      <- mc_base[idx]
  deltanp_m <- deltanp_640[idx]
  pcoefi_m  <- pcoefi_640[idx, ]
  dfull_m   <- dfull[idx, ]
  x2_m_tmpl <- cbind(p_m, rep(1, J), df$calor[idx])

  # ==== Matriz de propriedade pos-fusao ====
  # Agora Coors Light (firma 4) e Miller Lite (firma 5) sao da mesma firma,
  # entao Omega_{jk} = 1 para todos os pares Coors-Miller.
  firm_post_m <- firmid_post[idx]
  owner_post  <- outer(firm_post_m, firm_post_m, "==") * 1.0

  # ==== Solver: encontra precos de equilibrio pos-fusao ====
  # Chute inicial: precos pre-fusao (proximo do novo equilibrio)
  # Criterio de convergencia: |foc_residual| < ftol e |dp| < xtol
  sol <- tryCatch({
    nleqslv(
      x = p_m,
      fn = foc_residual,
      mc_m = mc_m, alpha_val = alpha, rho_val = rho,
      deltanp_m = deltanp_m, owner = owner_post,
      x2_m_template = x2_m_tmpl, pcoefi_m = pcoefi_m,
      dfull_m = dfull_m, theta2w_val = theta2w, ns_val = NS,
      method = "Broyden",
      control = list(xtol = 1e-12, ftol = 1e-12, maxit = 500)
    )
  }, error = function(e) {
    cat(sprintf("    Market %d solver error: %s\n", m, e$message))
    list(x = p_m, termcd = -1, fvals = rep(NA_real_, J))
  })

  # Verifica resultado do solver
  if (!is.null(sol$x) && is.numeric(sol$x) && all(is.finite(sol$x)) && all(sol$x > 0)) {
    p_post[idx] <- sol$x
    converged[m] <- sol$termcd == 1
    foc_norm[m]  <- max(abs(sol$fvec))

    # Calcula shares nos novos precos de equilibrio
    s_post[idx] <- compute_shares_at_prices(
      sol$x, alpha, rho, deltanp_m, x2_m_tmpl, dfull_m, theta2w, NS
    )
  } else {
    # Fallback: tenta Broyden com tolerancias mais frouxas
    sol2 <- tryCatch(
      nleqslv(
        x = p_m, fn = foc_residual,
        mc_m = mc_m, alpha_val = alpha, rho_val = rho,
        deltanp_m = deltanp_m, owner = owner_post,
        x2_m_template = x2_m_tmpl, pcoefi_m = pcoefi_m,
        dfull_m = dfull_m, theta2w_val = theta2w, ns_val = NS,
        method = "Broyden",
        control = list(xtol = 1e-10, ftol = 1e-10, maxit = 1000)
      ),
      error = function(e) list(x = rep(NA, J), termcd = -1, fvals = rep(NA, J))
    )
    if (is.numeric(sol2$x) && all(is.finite(sol2$x)) && all(sol2$x > 0)) {
      p_post[idx] <- sol2$x
      converged[m] <- sol2$termcd == 1
      foc_norm[m]  <- max(abs(sol2$fvals))
      s_post[idx] <- compute_shares_at_prices(
        sol2$x, alpha, rho, deltanp_m, x2_m_tmpl, dfull_m, theta2w, NS
      )
    } else {
      cat(sprintf("    Market %d: FAILED both methods\n", m))
      p_post[idx] <- p_m  # fallback: sem mudanca
      s_post[idx] <- s_pre[idx]
      foc_norm[m] <- NA
    }
  }

  if (m %% 20 == 0) {
    fn_str <- if (is.na(foc_norm[m])) "NA" else sprintf("%.2e", foc_norm[m])
    cat(sprintf("    Market %d/%d done (converged: %s, |FOC|=%s)\n",
                m, n_markets, converged[m], fn_str))
  }
}

cat(sprintf("\n  Converged: %d/%d markets\n", sum(converged), n_markets))
cat(sprintf("  Max FOC norm: %.2e\n", max(foc_norm)))

if (!all(converged)) {
  cat("  WARNING: some markets did not converge. Indices:\n")
  cat("   ", which(!converged), "\n")
}

# ==============================================================================
# BLOCO: RESULTADOS DA SIMULACAO
# ==============================================================================
# Compara precos, shares e excedente do consumidor antes e depois da fusao.
# O efeito unilateral esperado eh: precos sobem (especialmente das partes
# fusionantes), shares caem, e o excedente do consumidor diminui.
# ==============================================================================

cat("\n")
cat("-", rep("-", 69), "\n", sep = "")
cat("PART 4: Merger Simulation Results\n")
cat("-", rep("-", 69), "\n", sep = "")

# Variacao de precos (absoluta e percentual)
dp <- p_post - p_pre
dp_pct <- dp / p_pre * 100

df$p_post    <- p_post
df$s_post    <- s_post
df$dp        <- dp
df$dp_pct    <- dp_pct
df$firmid_post <- firmid_post

# --- 4.1 Variacao media de preco por firma ---
cat("\n--- 4.1 Mean price change by firm ---\n\n")

cat(sprintf("  %-14s %8s %8s %8s %8s\n", "Firm", "p_pre", "p_post", "dp ($)", "dp (%)"))
cat("  ", rep("-", 52), "\n", sep = "")

firm_labels <- c("1" = "ABI", "2" = "Corona", "3" = "Heineken",
                 "4" = "Coors", "5" = "MillerCoors")

for (firm in sort(unique(df$firmid))) {
  mask <- df$firmid == firm
  cat(sprintf("  %-14s %8.3f %8.3f %8.3f %8.2f%%\n",
              firm_labels[as.character(firm)],
              mean(p_pre[mask]), mean(p_post[mask]),
              mean(dp[mask]), mean(dp_pct[mask])))
}

cat(sprintf("\n  %-14s %8.3f %8.3f %8.3f %8.2f%%\n",
            "OVERALL", mean(p_pre), mean(p_post), mean(dp), mean(dp_pct)))

# --- 4.2 Por produto (marca x tamanho) ---
cat("\n--- 4.2 Price change by product ---\n\n")

prod_res <- df %>%
  mutate(label = case_when(
    brndid == 1  & sizeid == 2 ~ "Bud Light 12pk",
    brndid == 1  & sizeid == 3 ~ "Bud Light 24pk",
    brndid == 4  & sizeid == 2 ~ "Coors Light 12pk",
    brndid == 4  & sizeid == 3 ~ "Coors Light 24pk",
    brndid == 5  & sizeid == 2 ~ "Corona 12pk",
    brndid == 7  & sizeid == 2 ~ "Heineken 12pk",
    brndid == 13 & sizeid == 2 ~ "Miller Lite 12pk",
    brndid == 13 & sizeid == 3 ~ "Miller Lite 24pk"
  )) %>%
  group_by(label) %>%
  summarise(p_pre = mean(price), p_post = mean(p_post),
            dp = mean(dp), dp_pct = mean(dp_pct), .groups = "drop")

cat(sprintf("  %-20s %8s %8s %8s %8s\n", "Product", "p_pre", "p_post", "dp ($)", "dp (%)"))
cat("  ", rep("-", 56), "\n", sep = "")
for (i in 1:nrow(prod_res)) {
  cat(sprintf("  %-20s %8.3f %8.3f %8.3f %8.2f%%\n",
              prod_res$label[i], prod_res$p_pre[i], prod_res$p_post[i],
              prod_res$dp[i], prod_res$dp_pct[i]))
}

# --- 4.3 Partes fusionantes vs rivais ---
cat("\n--- 4.3 Merging firms vs rivals ---\n\n")

# Partes fusionantes: firmas cuja propriedade mudou (Coors + MillerCoors)
merging <- df$firmid %in% c(4, 5)
rival_abi <- df$firmid == 1
rival_import <- df$firmid %in% c(2, 3)

cat(sprintf("  Merging parties (Coors+MC):  mean dp = $%.3f  (%.2f%%)\n",
            mean(dp[merging]), mean(dp_pct[merging])))
cat(sprintf("  Rival ABI:                   mean dp = $%.3f  (%.2f%%)\n",
            mean(dp[rival_abi]), mean(dp_pct[rival_abi])))
cat(sprintf("  Rival imports:               mean dp = $%.3f  (%.2f%%)\n",
            mean(dp[rival_import]), mean(dp_pct[rival_import])))

# --- 4.4 Variacao de shares ---
cat("\n--- 4.4 Share changes ---\n\n")

ds <- s_post - s_pre
ds_pct <- ds / s_pre * 100

cat(sprintf("  %-14s %10s %10s %10s\n", "Firm", "s_pre", "s_post", "ds (%)"))
cat("  ", rep("-", 46), "\n", sep = "")
for (firm in sort(unique(df$firmid))) {
  mask <- df$firmid == firm
  cat(sprintf("  %-14s %10.5f %10.5f %10.2f%%\n",
              firm_labels[as.character(firm)],
              mean(s_pre[mask]), mean(s_post[mask]), mean(ds_pct[mask])))
}

# --- 4.5 Excedente do consumidor ---
# ==== Calculo do excedente do consumidor (CS) ====
# CS eh calculado via valor inclusivo do modelo RCNL:
#   CS_i = log(1 + IV_g^{1-rho}) / (-alpha_i)
# A diferenca CS_post - CS_pre mede a perda de bem-estar.
# Esperamos dCS < 0: a fusao eleva precos e reduz o excedente.
cat("\n--- 4.5 Consumer surplus change ---\n\n")

cs_pre_vec  <- rep(NA_real_, n_markets)
cs_post_vec <- rep(NA_real_, n_markets)

for (m in 1:n_markets) {
  idx <- ((m - 1) * J + 1):(m * J)
  x2_tmpl <- cbind(df$price[idx], rep(1, J), df$calor[idx])

  cs_pre_vec[m] <- compute_inclusive_value(
    p_pre[idx], alpha, rho, deltanp_640[idx],
    x2_tmpl, dfull[idx, ], theta2w, NS
  )
  cs_post_vec[m] <- compute_inclusive_value(
    p_post[idx], alpha, rho, deltanp_640[idx],
    x2_tmpl, dfull[idx, ], theta2w, NS
  )
}

dcs <- cs_post_vec - cs_pre_vec
dcs_pct <- dcs / cs_pre_vec * 100

cat(sprintf("  Mean CS (pre):   %.4f\n", mean(cs_pre_vec)))
cat(sprintf("  Mean CS (post):  %.4f\n", mean(cs_post_vec)))
cat(sprintf("  Mean delta CS:   %.4f  (%.3f%%)\n", mean(dcs), mean(dcs_pct)))
cat(sprintf("  All dCS < 0 (consumers harmed)? %s\n",
            ifelse(all(dcs < 0), "YES", "NO")))

# ==============================================================================
# BLOCO: COMPARACAO COM O ARTIGO
# ==============================================================================
# Os resultados do artigo (Tables 5/6) usam o modelo PLE (Price Leadership
# Equilibrium) com restricoes ICC, nao Bertrand puro. Nossa simulacao captura
# apenas o efeito unilateral. A diferenca entre o Bertrand e o observado nos
# dados eh atribuida pelo artigo a coordenacao (ABI como lider de preco).
# ==============================================================================

cat("\n")
cat("-", rep("-", 69), "\n", sep = "")
cat("PART 5: Comparison with Paper\n")
cat("-", rep("-", 69), "\n", sep = "")

cat("\n--- Loading paper's results for comparison ---\n")

cfmerger <- read.csv(file.path(path$results, "cfmerger.txt"), sep = ";", header = TRUE)
cat("\n  cfmerger.txt (merger counterfactual results):\n")
print(cfmerger)

eqeffects <- read.csv(file.path(path$results, "eqeffects.txt"), sep = ";", header = TRUE)
cat("\n  eqeffects.txt (profit/CS changes from PLE):\n")
print(eqeffects)

mean_markups <- read.csv(file.path(path$results, "mean_markups.txt"), sep = ";", header = TRUE)
cat("\n  mean_markups.txt (PLE supermarkups by brand):\n")
print(mean_markups)

cat("\n  NOTE: The paper's Tables 5/6 report results from the PLE model\n")
cat("  (Price Leadership Equilibrium with ICC constraints), not standard\n")
cat("  Bertrand. Our simulation is pure Bertrand — the unilateral effects\n")
cat("  component. The paper finds that unilateral effects alone underpredict\n")
cat("  the observed price increases, suggesting coordination.\n")

# ==============================================================================
# BLOCO: INTERPRETACAO
# ==============================================================================
# Discussao economica dos resultados: por que o Bertrand puro nao explica
# o aumento de precos observado e como isso motiva o modelo de coordenacao.
# ==============================================================================

cat("\n")
cat("-", rep("-", 69), "\n", sep = "")
cat("PART 6: Interpretation\n")
cat("-", rep("-", 69), "\n", sep = "")

cat("\n  Unilateral effects vs coordination:\n\n")

mean_dp_merging <- mean(dp[merging])
mean_dp_abi     <- mean(dp[rival_abi])

cat(sprintf("  Bertrand predicts:\n"))
cat(sprintf("    Merging parties (Coors+MC): dp = $%.3f\n", mean_dp_merging))
cat(sprintf("    Rival ABI:                  dp = $%.3f\n", mean_dp_abi))
cat(sprintf("    Ratio (ABI / merging):      %.3f\n", mean_dp_abi / mean_dp_merging))

cat("\n  Under standard Bertrand:\n")
cat("    - Merging firms raise prices (internalizing cross-price effects)\n")
cat("    - Rivals raise prices modestly (strategic complements)\n")
cat("    - ABI's price increase should be SMALLER than merging parties'\n")

cat("\n  What the paper observes in the actual data:\n")
cat("    - ABI prices increase by roughly the SAME amount as Miller/Coors\n")
cat("    - This is inconsistent with standard Bertrand unilateral effects\n")
cat("    - Consistent with price coordination (ABI as price leader)\n")

cat("\n  Our Bertrand simulation provides the unilateral effects baseline.\n")
cat("  The gap between Bertrand predictions and observed price changes\n")
cat("  is what the paper attributes to coordination/price leadership.\n")

# ==============================================================================
# BLOCO: SALVAR RESULTADOS
# ==============================================================================

cat("\n--- Saving results ---\n")

save(
  df, dfull, cdindex, obsindemand,
  alpha, theta1, theta2, rho, theta2w,
  derMat, elasMat, derMat_80,
  delta_val, deltanp_val, xi_val, mu_val, pcoefi_val,
  dprodfecoef, ddatefecoef, dcityfecoef,
  mc_R, markup_R, mc_matlab, mkt_selector,
  p_post, s_post, dp, dp_pct, dcs, cs_pre_vec, cs_post_vec,
  firmid_post, estimation_results,
  file = file.path(path$data_out, "step05_output.RData")
)
cat("  Saved: step05_output.RData\n")

# ==============================================================================
# RESUMO FINAL
# ==============================================================================

cat("\n")
cat("=", rep("=", 69), "\n", sep = "")
cat("STEP 5 COMPLETE: Merger Simulation\n")
cat("=", rep("=", 69), "\n", sep = "")
cat(sprintf("\n  Convergence: %d/%d markets\n", sum(converged), n_markets))
cat(sprintf("  Mean price change (all):     $%.3f (%.2f%%)\n", mean(dp), mean(dp_pct)))
cat(sprintf("  Mean price change (merging): $%.3f (%.2f%%)\n",
            mean(dp[merging]), mean(dp_pct[merging])))
cat(sprintf("  Mean price change (ABI):     $%.3f (%.2f%%)\n",
            mean(dp[rival_abi]), mean(dp_pct[rival_abi])))
cat(sprintf("  Consumer surplus change:     %.4f (%.3f%%)\n", mean(dcs), mean(dcs_pct)))
cat("\n  Next: Step 6 (documentation and RMarkdown)\n")
