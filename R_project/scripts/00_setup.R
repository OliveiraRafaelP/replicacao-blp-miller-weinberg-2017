# ==============================================================================
# 00_setup.R — Packages, paths, and helper functions
# Miller & Weinberg (2017) Replication in R
# ==============================================================================
#
# Este script inicializa o ambiente de trabalho para a replicacao de
# Miller & Weinberg (2017), "Understanding the Price Effects of the
# MillerCoors Joint Venture". Ele carrega pacotes, define caminhos,
# constantes globais e funcoes auxiliares reutilizadas nos demais scripts.
#
# Estrutura:
#   1. Instalacao/carregamento de pacotes
#   2. Definicao de caminhos do projeto
#   3. Constantes-chave (num. consumidores simulados, especificacao, fator desconto)
#   4. Funcoes auxiliares: decode_id, make_ownership, make_dummies,
#      contraction_map, compute_shares_rcnl, compute_derivatives
# ==============================================================================

# ==== BLOCO: Pacotes ====
# Instala (se necessario) e carrega todos os pacotes exigidos pelo projeto.
# Cada pacote tem papel especifico na estimacao BLP/RCNL.
required_pkgs <- c(
  "R.matlab",      # Load .mat files
  "tidyverse",     # Data manipulation + ggplot2
  "ivreg",         # 2SLS estimation (AER alternative)
  "nleqslv",       # Nonlinear equation solver (fsolve equivalent)
  "Matrix",        # Sparse matrices
  "sandwich",      # Clustered standard errors
  "lmtest",        # Coefficient tests with robust SEs
  "knitr",         # RMarkdown tables
  "kableExtra"     # Table formatting
)

options(repos = c(CRAN = "https://cran.r-project.org"))

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, quiet = TRUE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ==== BLOCO: Caminhos do projeto ====
# Lista nomeada com todos os diretorios relevantes.
# Facilita referencia cruzada entre scripts sem hardcodar caminhos repetidos.
path <- list(
  root     = normalizePath("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set"),
  data_raw = normalizePath("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/replication/data/raw"),
  data_ana = normalizePath("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/replication/data/analysis"),
  results  = normalizePath("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/replication/results"),
  rproject = normalizePath("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/R_project"),
  scripts  = normalizePath("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/R_project/scripts"),
  output   = normalizePath("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/R_project/output"),
  data_out = normalizePath("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/R_project/data")
)

# ==== BLOCO: Constantes globais ====
# Parametros fixos usados em toda a estimacao.
NS       <- 500L     # Number of simulated consumers
SPEC     <- "RCNL2"  # Demand specification (fiscal year starts Q4)
DF_MAIN  <- 0.26     # Main discount factor (from Table 3 model selection)

# Anos fiscais analisados: pre-fusao (2006, 2007) e pos-fusao (2010, 2011)
FISC_YEARS <- c(2006L, 2007L, 2010L, 2011L)

# ==== BLOCO: Funcao decode_id ====
# Objetivo: Extrair os componentes individuais (firma, marca, tamanho,
#   cidade, ano, mes) de um identificador composto numerico.
# Entrada: id2 — vetor numerico com IDs compostos de 11 digitos.
# Saida:   tibble com 6 colunas (firmid, brndid, sizeid, cityid, yearid, montid).
# Intuicao: Os dados originais do scanner (IMS) empacotam varias dimensoes
#   em um unico numero. Essa funcao reverte o processo para permitir
#   manipulacao por dimensao.
#' Decode the composite product ID used in small_scanner.mat
#' ID = firmid*1e10 + brndid*1e8 + sizeid*1e6 + cityid*1e4 + yearid*1e2 + montid
decode_id <- function(id2) {
  firmid <- floor(id2 / 1e10)
  brndid <- floor((id2 - firmid * 1e10) / 1e8)
  sizeid <- floor((id2 - firmid * 1e10 - brndid * 1e8) / 1e6)
  cityid <- floor((id2 - firmid * 1e10 - brndid * 1e8 - sizeid * 1e6) / 1e4)
  yearid <- floor((id2 - firmid * 1e10 - brndid * 1e8 - sizeid * 1e6 - cityid * 1e4) / 1e2)
  montid <- id2 - firmid * 1e10 - brndid * 1e8 - sizeid * 1e6 - cityid * 1e4 - yearid * 1e2

  tibble(firmid, brndid, sizeid, cityid, yearid, montid)
}

# ==== BLOCO: Funcao make_ownership ====
# Objetivo: Construir a matriz de propriedade (ownership matrix).
# Entrada: firmid_vec — vetor de IDs de firma para cada produto j no mercado.
# Saida:   Matriz logica J x J onde Owner[j,k] = TRUE se j e k pertencem
#          a mesma firma.
# Intuicao economica: Na teoria de oligopolio multi-produto, a firma
#   internaliza os efeitos cruzados de preco entre seus proprios produtos.
#   A matriz de propriedade determina quais termos de derivada cruzada
#   entram na condicao de primeira ordem (FOC) de precificacao.
#' Construct ownership matrix (f_ownMat.m equivalent)
#' Returns logical matrix: Owner[j,k] = TRUE if same firm
make_ownership <- function(firmid_vec) {
  n <- length(firmid_vec)
  outer(firmid_vec, firmid_vec, "==")
}

# ==== BLOCO: Funcao make_dummies ====
# Objetivo: Gerar uma matriz de variaveis dummy a partir de um vetor categorico.
# Entrada: x — vetor categorico (ex.: prodid, dateid).
# Saida:   Matriz 0/1 de dimensao length(x) x (num. categorias unicas).
# Equivale a cr_dum.m no codigo Matlab original.
#' Create dummy variable matrix from categorical vector
make_dummies <- function(x) {
  lvls <- sort(unique(x))
  mat <- matrix(0L, nrow = length(x), ncol = length(lvls))
  for (i in seq_along(lvls)) {
    mat[, i] <- as.integer(x == lvls[i])
  }
  colnames(mat) <- paste0("d_", lvls)
  mat
}

# ==== BLOCO: Funcao contraction_map ====
# Objetivo: Implementar o mapeamento de contracao de Berry, Levinsohn &
#   Pakes (1995) — BLP — para recuperar as utilidades medias (delta)
#   que racionalizam os market shares observados.
# Entrada:
#   s_obs   — market shares observados (vetor J x 1)
#   delta0  — valores iniciais de utilidade media
#   mu      — matriz de heterogeneidade dos consumidores (J x ns)
#   rho     — parametro de nesting (RCNL)
#   tol     — tolerancia de convergencia (padrao 1e-12)
#   max_iter — maximo de iteracoes
# Saida:   delta convergido (vetor J x 1)
# Intuicao: A cada iteracao, ajusta delta somando o log-erro entre
#   shares observados e preditos. Convergencia garante que o modelo
#   reproduz exatamente os shares de mercado observados.
#' BLP contraction mapping for RCNL
#' Iterates: delta^{r+1} = delta^r + log(s_obs) - log(s_pred(delta^r))
#' @param s_obs Observed market shares (J x 1)
#' @param delta0 Starting values for mean utility (J x 1)
#' @param mu Consumer heterogeneity matrix (J x ns)
#' @param rho Nesting parameter
#' @param tol Convergence tolerance
#' @param max_iter Maximum iterations
contraction_map <- function(s_obs, delta0, mu, rho, tol = 1e-12, max_iter = 1000) {
  delta <- delta0
  for (iter in seq_len(max_iter)) {
    s_pred <- compute_shares_rcnl(delta, mu, rho)
    delta_new <- delta + log(s_obs) - log(s_pred)
    if (max(abs(delta_new - delta)) < tol) {
      return(delta_new)
    }
    delta <- delta_new
  }
  warning("Contraction mapping did not converge in ", max_iter, " iterations")
  delta
}

# ==== BLOCO: Funcao compute_shares_rcnl ====
# Objetivo: Calcular market shares preditos pelo modelo RCNL
#   (Random Coefficients Nested Logit).
# Entrada:
#   delta — utilidade media (vetor J x 1)
#   mu    — heterogeneidade do consumidor (matriz J x ns)
#   rho   — parametro de nesting (0 = logit puro; proximo de 1 = forte correlacao)
# Saida:   Vetor J x 1 de market shares preditos, mediados sobre ns consumidores.
# Intuicao economica: O modelo RCNL generaliza o logit ao permitir:
#   (a) coeficientes aleatorios (via mu), capturando heterogeneidade de
#       preferencias entre consumidores; e
#   (b) correlacao intra-grupo (via rho), modelando substituicao mais
#       intensa entre produtos no mesmo ninho.
#   Shares sao calculados individualmente por consumidor simulado e depois
#   mediados (simulacao de Monte Carlo).
#' Compute predicted shares under RCNL
#' @param delta Mean utility vector (J x 1)
#' @param mu Consumer heterogeneity (J x ns)
#' @param rho Nesting parameter (0 = logit, >0 = nested)
#' @return Predicted market shares (J x 1), averaged over consumers
compute_shares_rcnl <- function(delta, mu, rho) {
  ns <- ncol(mu)
  J  <- length(delta)

  # Utilidade individual: V_ij = delta_j + mu_ij
  V <- matrix(delta, nrow = J, ncol = ns) + mu

  # Valor inclusivo dentro do grupo (todos produtos em um unico ninho)
  # IV_G = log(sum_j exp(V_ij / (1-rho)))
  V_scaled <- V / (1 - rho)
  log_sum <- apply(V_scaled, 2, function(col) {
    max_v <- max(col)
    max_v + log(sum(exp(col - max_v)))  # log-sum-exp for stability
  })

  # Probabilidade do grupo (inside): P(inside) = IV^(1-rho) / (1 + IV^(1-rho))
  log_IV_1mrho <- (1 - rho) * log_sum  # = log(IV^(1-rho))
  log_denom <- log(1 + exp(log_IV_1mrho))  # log(1 + IV^(1-rho))

  # Probabilidade condicional: P(j|inside) = exp(V_ij/(1-rho)) / sum_k exp(V_ik/(1-rho))
  log_cond <- V_scaled - matrix(log_sum, nrow = J, ncol = ns, byrow = TRUE)

  # Probabilidade incondicional: P(j) = P(j|inside) * P(inside)
  log_uncond <- log_cond + matrix(log_IV_1mrho - log_denom, nrow = J, ncol = ns, byrow = TRUE)

  shares_ind <- exp(log_uncond)
  rowMeans(shares_ind)  # Media sobre consumidores simulados
}

# ==== BLOCO: Funcao compute_derivatives ====
# Objetivo: Calcular a matriz J x J de derivadas de preco ds_j/dp_k.
# Entrada:
#   pcoefi  — coeficientes individuais de preco (J x ns): alpha + a_i
#   sharei  — shares individuais (J x ns)
#   scondi  — shares condicionais dentro do ninho (J x ns)
#   sgroupi — probabilidade de escolha do grupo/inside good (1 x ns ou J x ns)
#   rho     — parametro de nesting
# Saida:   Matriz J x J de derivadas ds_j/dp_k.
# Intuicao economica: Essa matriz entra diretamente na condicao de
#   primeira ordem (FOC) de Bertrand-Nash. Elementos diagonais (derivada
#   propria) medem a sensibilidade da demanda do produto j ao seu proprio
#   preco; elementos fora da diagonal (derivada cruzada) medem substituicao
#   entre produtos. Combinada com a ownership matrix, determina os markups.
#' Compute matrix of price derivatives ds_j/dp_k
#' @param pcoefi Individual price coefficients (J x ns): alpha + ai
#' @param sharei Individual shares (J x ns)
#' @param scondi Conditional shares within nest (J x ns)
#' @param sgroupi Group (inside) probability (1 x ns or J x ns)
#' @param rho Nesting parameter
#' @return J x J derivative matrix
compute_derivatives <- function(pcoefi, sharei, scondi, sgroupi, rho) {
  J  <- nrow(sharei)
  ns <- ncol(sharei)

  der <- matrix(0, nrow = J, ncol = J)

  for (j in seq_len(J)) {
    for (k in seq_len(J)) {
      if (j == k) {
        # Derivada propria: efeito do preco de j sobre sua propria demanda
        term <- pcoefi[j, ] * sharei[j, ] *
          (1 / (1 - rho) - (1 / (1 - rho)) * scondi[j, ] - sgroupi[1, ] + sharei[j, ])
        # Simplified: more accurate version below
        term <- pcoefi[j, ] * sharei[j, ] *
          (1 - (1 / (1 - rho)) * scondi[j, ] + (rho / (1 - rho)) * sgroupi[1, ] * scondi[j, ] /
             (sgroupi[1, ] + 1e-300) - sharei[j, ])
      } else {
        # Derivada cruzada: efeito do preco de k sobre a demanda de j
        term <- -pcoefi[k, ] * sharei[j, ] *
          ((1 / (1 - rho)) * scondi[k, ] - (rho / (1 - rho)) * sgroupi[1, ] * scondi[k, ] /
             (sgroupi[1, ] + 1e-300) + sharei[k, ])
      }
      # This is a simplified version; exact implementation follows rcnl_der1.m
      # Will be refined in Step 4 after validating against derMat_2
      der[j, k] <- mean(term)  # Media sobre consumidores simulados
    }
  }
  der
}

cat("Setup complete. Paths, constants, and helpers loaded.\n")
cat(sprintf("Working directory: %s\n", path$root))
cat(sprintf("Specification: %s | ns=%d | df=%.2f\n", SPEC, NS, DF_MAIN))
