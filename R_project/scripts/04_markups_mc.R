# ==============================================================================
# 04_markups_mc.R — Markups, marginal costs, ownership matrices
# Miller & Weinberg (2017) Replication in R
# Uses RCNL parameters from dresgmm2.mat as baseline
# ==============================================================================
#
# Este script recupera markups e custos marginais a partir da condicao de
# primeira ordem (CPO) do oligopolio de Bertrand multiproduto.
#
# A ideia central eh: dado que observamos precos e shares de equilibrio,
# e estimamos a funcao de demanda (elasticidades), podemos inverter a CPO
# para obter o markup implicito de cada produto em cada mercado:
#
#   markup = -(Omega .* der')^{-1} * s
#
# onde Omega eh a matriz de propriedade, der eh a matriz de derivadas ds/dp,
# e s eh o vetor de shares. O custo marginal segue como mc = p - markup.
# ==============================================================================

load("C:/Users/rafael.oliveira/Desktop/Problem Set/R_project/data/step03_output.RData")
source("C:/Users/rafael.oliveira/Desktop/Problem Set/R_project/scripts/00_setup.R")

cat("\n")
cat("=", rep("=", 69), "\n", sep = "")
cat("STEP 4: MARKUPS AND MARGINAL COSTS\n")
cat("=", rep("=", 69), "\n", sep = "")

# ==============================================================================
# BLOCO: FUNCOES RCNL
# ==============================================================================
# Estas funcoes traduzem o codigo Matlab do artigo (cf_foc_partial.m,
# rcnl_indsh.m, rcnl_der1.m) para R. Elas calculam:
#   (1) a heterogeneidade do consumidor (mu),
#   (2) as probabilidades de escolha individuais,
#   (3) as derivadas ds_j/dp_k,
#   (4) a avaliacao completa da CPO de Bertrand a precos dados.
# ==============================================================================

cat("\n--- Loading RCNL parameters ---\n")
cat(sprintf("  alpha = %.6f,  rho = %.6f,  ns = %d\n", alpha, rho, NS))

# --------------------------------------------------------------------------
# compute_mu: desvio de utilidade individual mu_ij
# --------------------------------------------------------------------------
# Entrada:
#   theta2w — matriz K x 2 de interacoes com renda (coluna 2 = renda)
#   ns      — numero de consumidores simulados
#   dfull_m — J x ns, sorteios de renda (mesmos para todos os produtos no mercado)
#   x2_m    — J x K, caracteristicas dos produtos (preco, constante, calorias)
# Saida:
#   mu — J x ns, desvio de utilidade por consumidor-produto
#   ai — J x ns, desvio do coeficiente de preco por consumidor
# Intuicao economica:
#   Captura como a sensibilidade ao preco e a preferencia por calorias variam
#   com a renda do consumidor. Consumidores de alta renda sao menos sensiveis ao preco.
compute_mu <- function(theta2w, ns, dfull_m, x2_m) {
  J <- nrow(x2_m)

  # Sorteios de renda (identicos para todos os J produtos dentro do mercado)
  income_draws <- dfull_m[1, ]  # ns-length vector
  theta2_col2 <- theta2w[, 2]   # K-length

  # mu[j,i] = (sum_k x2[j,k] * theta2_col2[k]) * income_draws[i]
  # Produto escalar das caracteristicas com os pesos de renda, escalado pela renda do consumidor
  product_loading <- as.numeric(x2_m %*% theta2_col2)  # J-length
  mu <- outer(product_loading, income_draws)             # J x ns

  # ai[j,i] = theta2w[1,2] * income_draws[i]
  # Desvio individual do coeficiente de preco (alpha_i = alpha + ai)
  ai <- matrix(theta2w[1, 2] * income_draws, nrow = J, ncol = ns, byrow = TRUE)

  list(mu = mu, ai = ai)
}

# --------------------------------------------------------------------------
# compute_indsh: probabilidades de escolha individuais no modelo RCNL
# --------------------------------------------------------------------------
# Entrada:
#   exp_delta — J x 1, exp(delta_j) onde delta_j = utilidade media do produto j
#   exp_mu    — J x ns, exp(mu_ij) desvios de utilidade individuais
#   rho       — escalar, parametro de aninhamento (nesting)
# Saida:
#   sharei  — J x ns, probabilidade incondicional de i escolher j
#   scondi  — J x ns, probabilidade condicional s_{j|g,i} (dentro do ninho)
#   sgroupi — 1 x ns, probabilidade de escolher o grupo "inside"
# Intuicao economica:
#   No nested logit, a escolha ocorre em dois estagios: primeiro o consumidor
#   decide entre comprar cerveja (inside) ou nao (outside), depois escolhe
#   qual cerveja. O parametro rho governa a correlacao de preferencias dentro
#   do ninho — rho alto implica alta substituicao entre cervejas.
compute_indsh <- function(exp_delta, exp_mu, rho) {
  J <- length(exp_delta)
  ns <- ncol(exp_mu)

  # exp(V_ij / (1-rho)): utilidade escalada pelo parametro de nesting
  # A divisao por (1-rho) amplifica diferencas dentro do grupo
  expval <- matrix(exp_delta^(1 / (1 - rho)), nrow = J, ncol = ns) *
    exp_mu^(1 / (1 - rho))

  # Valor inclusivo do grupo inside: IV_g = sum_j exp(V_ij/(1-rho))
  # Agrega a atratividade de todas as cervejas para cada consumidor
  IV_g <- colSums(expval)  # ns-length

  # Share condicional: s_{j|g,i} = expval_{j,i} / IV_g_i
  scondi <- expval / matrix(IV_g, nrow = J, ncol = ns, byrow = TRUE)

  # Probabilidade do grupo inside: P(inside) = IV_g^{1-rho} / (1 + IV_g^{1-rho})
  # Segue a estrutura do logit aninhado (outside good tem utilidade normalizada a zero)
  IV_g_1mrho <- IV_g^(1 - rho)
  sgroupi <- matrix(IV_g_1mrho / (1 + IV_g_1mrho), nrow = 1, ncol = ns)

  # Share incondicional: s_{j,i} = s_{j|g,i} * P(inside)_i
  sharei <- scondi * matrix(sgroupi, nrow = J, ncol = ns, byrow = TRUE)

  list(sharei = sharei, scondi = scondi, sgroupi = sgroupi)
}

# --------------------------------------------------------------------------
# compute_der1: derivadas de share em relacao ao preco, ds_j/dp_k
# --------------------------------------------------------------------------
# Entrada:
#   pcoefi_m — J x ns, coeficiente de preco individual (alpha + ai)
#   sharei   — J x ns, shares individuais (output de compute_indsh)
#   scondi   — J x ns, shares condicionais dentro do ninho
#   sgroupi  — 1 x ns, probabilidade do grupo inside
#   rho      — escalar, parametro de aninhamento
# Saida:
#   der — J x J, matriz de derivadas ds_j/dp_k (media sobre consumidores simulados)
# Intuicao economica:
#   Cada elemento der[j,k] mede como a share do produto j varia quando o preco
#   de k muda marginalmente. A diagonal (j=k) deve ser negativa (aumento do preco
#   proprio reduz share). Os elementos fora da diagonal sao positivos (substituicao).
#   A media sobre ns consumidores simulados integra a heterogeneidade.
compute_der1 <- function(pcoefi_m, sharei, scondi, sgroupi, rho) {
  J  <- nrow(sharei)
  ns <- ncol(sharei)

  der <- matrix(0, nrow = J, ncol = J)

  for (j in 1:J) {
    for (k in 1:J) {
      if (j == k) {
        # Derivada propria: ds_j/dp_j (RCNL, formula analitica)
        # Combina efeito direto (dentro do ninho) e efeito entre ninhos
        term <- pcoefi_m[j, ] * sharei[j, ] *
          (1 / (1 - rho) -
             (1 / (1 - rho)) * scondi[j, ] -
             sgroupi[1, ] + sharei[j, ])
      } else {
        # Derivada cruzada: ds_j/dp_k
        # Captura substituicao: aumento de p_k desloca demanda para j
        term <- pcoefi_m[k, ] * sharei[j, ] *
          (-1 / (1 - rho) * scondi[k, ] +
             sgroupi[1, ] - sharei[k, ])
      }
      # Media sobre consumidores simulados (integracao Monte Carlo)
      der[j, k] <- mean(term)
    }
  }
  der
}

# --------------------------------------------------------------------------
# eval_foc: avalia a CPO de Bertrand a precos dados
# --------------------------------------------------------------------------
# Entrada:
#   p          — vetor de precos dos produtos que otimizam
#   alpha_val  — coeficiente medio de preco
#   rho_val    — parametro de aninhamento
#   deltanp_m  — J x 1, utilidade media excluindo o componente de preco
#   mc_m       — J x 1, custos marginais (placeholder se estamos recuperando mc)
#   owner      — J x J, matriz de propriedade (1 se mesma firma, 0 caso contrario)
#   pfix       — precos fixados (de produtos que nao otimizam)
#   popt       — vetor binario indicando quais produtos otimizam
#   x2_m       — J x K, caracteristicas dos produtos
#   pcoefi_m   — J x ns, coeficientes de preco individuais
#   dfull_m    — J x ns, sorteios de renda
#   theta2w_val — K x 2, interacoes com renda
#   ns_val     — numero de consumidores simulados
# Saida:
#   zero    — residuo da CPO: p - mc - markup (deve ser ~0 no equilibrio)
#   svec    — J x 1, shares de mercado
#   der1    — J x J, derivadas ds/dp
#   bMarkup — markup implicito pela inversao da CPO
# Intuicao economica:
#   A CPO de Bertrand multiproduto eh: s + (Omega .* ds/dp') * (p - mc) = 0
#   Rearranjando: markup = p - mc = -(Omega .* ds/dp')^{-1} * s
#   Esta funcao calcula tudo e retorna o residuo para verificacao ou uso em solver.
eval_foc <- function(p, alpha_val, rho_val, deltanp_m, mc_m, owner,
                     pfix, popt, x2_m, pcoefi_m, dfull_m, theta2w_val, ns_val) {

  # Constroi vetor completo de precos (produtos otimizantes + fixos)
  pvec <- rep(0, length(popt))
  pvec[popt == 1] <- p
  pvec[popt == 0] <- pfix

  # Utilidade media: delta_j = deltanp_j + alpha * p_j
  # deltanp contem todos os componentes exceto preco (FEs, xi, etc.)
  delta2 <- deltanp_m + alpha_val * pvec

  # Atualiza x2 com precos correntes (necessario para mu que depende de preco)
  x2_m[, 1] <- pvec

  # Calcula heterogeneidade do consumidor
  mu_res <- compute_mu(theta2w_val, ns_val, dfull_m, x2_m)
  mu2 <- mu_res$mu
  ai  <- mu_res$ai

  # Shares individuais via modelo RCNL
  sh_res <- compute_indsh(exp(delta2), exp(mu2), rho_val)
  sharei  <- sh_res$sharei
  scondi  <- sh_res$scondi
  sgroupi <- sh_res$sgroupi

  # Shares de mercado: media sobre consumidores simulados
  svec <- rowMeans(sharei)

  # Coeficiente de preco completo: alpha_i = alpha + ai (individual)
  pcoefi_full <- alpha_val + ai

  # Matriz de derivadas ds_j/dp_k
  der1 <- compute_der1(pcoefi_full, sharei, scondi, sgroupi, rho_val)

  # ==== Inversao da CPO para obter markup ====
  # Omega = owner .* der1' (produto elemento-a-elemento da propriedade com a transposta das derivadas)
  # markup = -Omega^{-1} * s (inversao matricial)
  opt_idx <- which(popt == 1)
  omega <- owner * t(der1[opt_idx, opt_idx])
  bMarkup <- -solve(omega, svec[opt_idx])

  # Residuo da CPO: p - mc - markup (zero no equilibrio)
  zero <- p - mc_m[popt == 1] - bMarkup

  list(zero = zero, svec = svec, der1 = der1, bMarkup = bMarkup)
}

# ==============================================================================
# BLOCO: ESTRUTURAS POR MERCADO
# ==============================================================================
# Reconstroi o seletor de mercados (140 brutos -> 80 filtrados) e extrai as
# fatias correspondentes de derMat (matriz de derivadas pre-computada pelo Matlab).
# ==============================================================================

cat("\n--- Building market-level structures ---\n")

# Reconstroi seletor de mercados a partir dos IDs compostos do Matlab
mat_reload <- R.matlab::readMat(file.path(path$data_raw, "small_scanner.mat"))
id2_all <- mat_reload$small.scanner[, 1]
fid_all <- floor(id2_all / 1e10)
bid_all <- floor((id2_all - fid_all * 1e10) / 1e8)
sid_all <- floor((id2_all - fid_all * 1e10 - bid_all * 1e8) / 1e6)
cid_all <- floor((id2_all - fid_all * 1e10 - bid_all * 1e8 - sid_all * 1e6) / 1e4)
yid_all <- floor((id2_all - fid_all * 1e10 - bid_all * 1e8 - sid_all * 1e6 - cid_all * 1e4) / 1e2)
mid_all <- round(id2_all - fid_all * 1e10 - bid_all * 1e8 - sid_all * 1e6 - cid_all * 1e4 - yid_all * 1e2)

obsintemp_all <- (yid_all <= 3 | yid_all >= 6) | (yid_all == 4 & mid_all <= 2) | (yid_all == 5 & mid_all >= 3)
fiscid_all_vec <- yid_all + as.integer(mid_all >= 4) + 2004L
obsin_all_vec <- obsintemp_all & (fiscid_all_vec %in% c(2006, 2007, 2010, 2011))
raw_mkt_idx <- rep(1:140, each = 8)
mkt_selector <- tapply(obsin_all_vec, raw_mkt_idx, all)

derMat_80 <- derMat[, , mkt_selector]
cat(sprintf("  derMat_80: %s\n", paste(dim(derMat_80), collapse = " x ")))
stopifnot(dim(derMat_80)[3] == 80)

# ==============================================================================
# BLOCO: IMPUTACAO DE CUSTOS MARGINAIS VIA BERTRAND
# ==============================================================================
# Para cada mercado, avaliamos a CPO de Bertrand nos precos e shares observados.
# A inversao da CPO nos da o markup implicito; subtraindo-o do preco, obtemos mc.
#
# Formalmente: mc_j = p_j - markup_j, onde markup = -(Omega .* der')^{-1} * s.
#
# Isso assume que as firmas jogam um equilibrio de Nash-Bertrand em precos,
# cada uma maximizando lucro sobre seus proprios produtos dado os precos dos rivais.
# ==============================================================================

cat("\n--- Computing Bertrand markups and marginal costs ---\n")
cat("  Method: Full RCNL FOC evaluation (recomputing derivatives)\n")
cat("  This matches impute_bertrand.m -> f_impute_mc.m -> cf_foc_partial.m\n\n")

# Seleciona as 640 observacoes filtradas
deltanp_640 <- deltanp_val   # already selected in step 01
pcoefi_640  <- pcoefi_val    # 640 x 500

J <- 8
n_markets <- 80

mc_R <- rep(NA_real_, J * n_markets)
markup_R <- rep(NA_real_, J * n_markets)
svec_R <- rep(NA_real_, J * n_markets)

for (m in 1:n_markets) {
  idx <- ((m - 1) * J + 1):(m * J)

  firm_m    <- df$firmid[idx]
  p_m       <- df$price[idx]
  s_m       <- df$share[idx]
  deltanp_m <- deltanp_640[idx]
  pcoefi_m  <- pcoefi_640[idx, ]   # J x 500
  dfull_m   <- dfull[idx, ]         # J x 500

  # x2 para este mercado: [preco, constante, calorias]
  x2_m <- cbind(p_m, rep(1, J), df$calor[idx])

  # ==== Matriz de propriedade (Nash-Bertrand) ====
  # Omega_{jk} = 1 se j e k pertencem a mesma firma, 0 caso contrario.
  # Isso faz a firma internalizar a substituicao entre seus proprios produtos.
  owner_m <- outer(firm_m, firm_m, "==") * 1.0

  # Todos os produtos otimizam (sem split coalizao/franja)
  popt_m <- rep(1, J)
  pfix_m <- numeric(0)
  mc_dummy <- rep(0, J)  # placeholder, nao usado no calculo do markup

  # Avalia a CPO nos precos observados para obter markup e derivadas
  foc_res <- eval_foc(
    p = p_m, alpha_val = alpha, rho_val = rho,
    deltanp_m = deltanp_m, mc_m = mc_dummy,
    owner = owner_m, pfix = pfix_m, popt = popt_m,
    x2_m = x2_m, pcoefi_m = pcoefi_m,
    dfull_m = dfull_m, theta2w_val = theta2w, ns_val = NS
  )

  # Recupera markup e custo marginal: mc = p - markup
  markup_R[idx] <- foc_res$bMarkup
  mc_R[idx]     <- p_m - foc_res$bMarkup
  svec_R[idx]   <- foc_res$svec

  if (m %% 20 == 0) cat(sprintf("    Market %d/%d done\n", m, n_markets))
}

cat("  All 80 markets computed.\n")

# ==============================================================================
# BLOCO: VALIDACAO CONTRA MATLAB
# ==============================================================================
# Compara os custos marginais calculados em R com os do codigo Matlab original.
# Diferencas pequenas (< 1e-4) sao esperadas por questoes de precisao numerica.
# ==============================================================================

cat("\n--- Validation against sres_bertrand.mat ---\n")

mc_matlab <- R.matlab::readMat(
  file.path(path$data_ana, "sres_bertrand.mat")
)$mc[, 1]

diff_mc <- mc_R - mc_matlab
cat(sprintf("  max|mc_R - mc_matlab|  = %.6f\n", max(abs(diff_mc))))
cat(sprintf("  mean|mc_R - mc_matlab| = %.6f\n", mean(abs(diff_mc))))
cat(sprintf("  RMSE                   = %.6f\n", sqrt(mean(diff_mc^2))))
cat(sprintf("  correlation            = %.8f\n", cor(mc_R, mc_matlab)))
cat(sprintf("\n  mc_R:      mean=%.4f, range=[%.4f, %.4f]\n",
            mean(mc_R), min(mc_R), max(mc_R)))
cat(sprintf("  mc_matlab: mean=%.4f, range=[%.4f, %.4f]\n",
            mean(mc_matlab), min(mc_matlab), max(mc_matlab)))

# Custos marginais negativos sao economicamente suspeitos
cat(sprintf("\n  Negative MC: R=%d, Matlab=%d (of %d)\n",
            sum(mc_R < 0), sum(mc_matlab < 0), length(mc_R)))

# ==============================================================================
# BLOCO: RESUMO DE MARKUPS E CUSTOS MARGINAIS
# ==============================================================================

cat("\n")
cat("-", rep("-", 69), "\n", sep = "")
cat("PART 5: Markup and MC Summary\n")
cat("-", rep("-", 69), "\n", sep = "")

df$markup_bert <- markup_R
df$mc_bert     <- mc_R
df$mc_matlab   <- mc_matlab

cat("\n  By firm:\n")
cat(sprintf("  %-12s %8s %8s %8s %8s %8s\n",
            "Firm", "price", "markup", "mc_R", "mc_matl", "mc_diff"))
cat("  ", rep("-", 58), "\n", sep = "")

for (firm in sort(unique(df$firmid))) {
  mask <- df$firmid == firm
  label <- c("1" = "ABI", "2" = "Corona", "3" = "Heineken",
             "4" = "Coors", "5" = "MillerCoors")[as.character(firm)]
  cat(sprintf("  %-12s %8.3f %8.3f %8.3f %8.3f %8.4f\n",
              label,
              mean(df$price[mask]),
              mean(markup_R[mask]),
              mean(mc_R[mask]),
              mean(mc_matlab[mask]),
              mean(abs(diff_mc[mask]))))
}

cat("\n  By brand x size:\n")
cat(sprintf("  %-20s %8s %8s %8s %8s\n",
            "Product", "price", "markup", "mc_R", "mc_matlab"))
cat("  ", rep("-", 58), "\n", sep = "")

prod_summary <- df %>%
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
  summarise(
    price = mean(price),
    markup = mean(markup_bert),
    mc_R = mean(mc_bert),
    mc_matlab = mean(mc_matlab),
    .groups = "drop"
  )

for (i in 1:nrow(prod_summary)) {
  cat(sprintf("  %-20s %8.3f %8.3f %8.3f %8.3f\n",
              prod_summary$label[i],
              prod_summary$price[i],
              prod_summary$markup[i],
              prod_summary$mc_R[i],
              prod_summary$mc_matlab[i]))
}

# ==============================================================================
# BLOCO: MATRIZES DE PROPRIEDADE
# ==============================================================================
# A matriz de propriedade (ownership matrix) eh fundamental para o modelo de
# conducta. Ela define quais produtos sao co-otimizados pela mesma firma.
#
# Pre-fusao: Coors (firma 4) e MillerCoors (firma 5) sao independentes.
# Pos-fusao: Coors eh absorvida por MillerCoors — seus produtos passam a ser
# co-otimizados, gerando incentivo para aumentar precos (internalizacao de
# externalidades de substituicao entre Coors Light e Miller Lite).
# ==============================================================================

cat("\n")
cat("-", rep("-", 69), "\n", sep = "")
cat("PART 6: Ownership Matrices\n")
cat("-", rep("-", 69), "\n", sep = "")

# Pre-merger market (FY 2006 or 2007): firms 1, 2, 3, 4, 5
# Post-merger market (FY 2010 or 2011): firms 1, 2, 3, 5 (no firm 4)

# Exemplo: mercado 1 (pre-fusao)
m1_idx <- 1:8
firm_m1 <- df$firmid[m1_idx]
owner_pre <- outer(firm_m1, firm_m1, "==") * 1L

cat("\n  Pre-merger ownership (market 1, FY", df$fiscid[1], "):\n")
cat("  Products:", paste(firm_m1, collapse = " "), "\n")
cat("  Firm labels: ABI ABI Coors Coors Corona Heineken MC MC\n")
print(owner_pre)

# Pos-fusao: firma 4 (Coors) absorvida pela firma 5 (MillerCoors)
m_post_idx <- which(df$fiscid >= 2010)[1:8]
firm_mpost <- df$firmid[m_post_idx]
owner_post <- outer(firm_mpost, firm_mpost, "==") * 1L

cat("\n  Post-merger ownership (first FY 2010 market):\n")
cat("  Products:", paste(firm_mpost, collapse = " "), "\n")
cat("  Firm labels: ABI ABI MC MC Corona Heineken MC MC\n")
print(owner_post)

cat("\n  Key difference: Coors Light (products 3-4) now co-owned with\n")
cat("  Miller Lite (products 7-8) under MillerCoors. Ownership matrix\n")
cat("  has more 1s in the post-merger period.\n")

# Contagem de pares co-propriedade: mais 1s = mais internalizacao
n_pairs_pre  <- sum(owner_pre[upper.tri(owner_pre)])
n_pairs_post <- sum(owner_post[upper.tri(owner_post)])
cat(sprintf("\n  Co-owned product pairs: pre=%d, post=%d (of %d total)\n",
            n_pairs_pre, n_pairs_post, choose(8, 2)))

# ==============================================================================
# BLOCO: INDICE DE LERNER (markup / preco)
# ==============================================================================
# O indice de Lerner mede o poder de mercado como fracao do preco.
# Lerner = (p - mc) / p = markup / p.
# Valores mais altos indicam maior poder de mercado.
# ==============================================================================

cat("\n")
cat("-", rep("-", 69), "\n", sep = "")
cat("PART 7: Lerner Index (markup / price)\n")
cat("-", rep("-", 69), "\n", sep = "")

df$lerner <- markup_R / df$price

cat(sprintf("\n  Overall: mean=%.4f, range=[%.4f, %.4f]\n",
    mean(df$lerner), min(df$lerner), max(df$lerner)))

cat("\n  By firm:\n")
for (firm in sort(unique(df$firmid))) {
  mask <- df$firmid == firm
  label <- c("1" = "ABI", "2" = "Corona", "3" = "Heineken",
             "4" = "Coors", "5" = "MillerCoors")[as.character(firm)]
  cat(sprintf("    %-12s: Lerner = %.4f  (markup = $%.2f on $%.2f price)\n",
              label, mean(df$lerner[mask]),
              mean(markup_R[mask]), mean(df$price[mask])))
}

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
  estimation_results,
  file = file.path(path$data_out, "step04_output.RData")
)
cat("  Saved: step04_output.RData\n")

# ==============================================================================
# RESUMO FINAL
# ==============================================================================

cat("\n")
cat("=", rep("=", 69), "\n", sep = "")
cat("STEP 4 COMPLETE: Markups and Marginal Costs\n")
cat("=", rep("=", 69), "\n", sep = "")
cat(sprintf("\n  Bertrand MC: mean=%.3f, all positive: %s\n",
            mean(mc_R), ifelse(all(mc_R > 0), "YES", "NO")))
cat(sprintf("  Validation vs Matlab: corr=%.6f, mean|diff|=%.4f\n",
            cor(mc_R, mc_matlab), mean(abs(mc_R - mc_matlab))))
cat(sprintf("  Markup: mean=$%.2f (Lerner=%.1f%%)\n",
            mean(markup_R), mean(df$lerner) * 100))
cat("\n  Next: Step 5 (merger simulation)\n")
