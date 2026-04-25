# ==============================================================================
# 02_clean_shares.R — Instruments, fixed effects, validation, diagnostic plots
# Miller & Weinberg (2017) Replication in R
# ==============================================================================
# Este script constroi os instrumentos de demanda (Z1-Z12), monta as matrizes
# de efeitos fixos, valida os dados contra os resultados do Matlab original e
# gera graficos diagnosticos. E o passo preparatorio para a estimacao 2SLS/BLP.
# ==============================================================================

# ==== BLOCO: Carregamento de dados e configuracao ====
# Carrega os objetos do passo anterior (df, dfull, elasMat, delta, etc.)
load("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/R_project/data/step01_output.RData")
source("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/R_project/scripts/00_setup.R")

cat("\n")
cat("=", rep("=", 59), "\n", sep = "")
cat("STEP 02: Instruments, FE matrices, and validation\n")
cat("=", rep("=", 59), "\n", sep = "")

# ==============================================================================
# 1. CONSTRUCT DEMAND INSTRUMENTS (Z1–Z12)
#    Reconstructed from MW(2017) Section 4.2, pp. 1775-1776
#    NOT present in replication package code
# ==============================================================================

cat("\n--- Constructing demand instruments ---\n")

# ==== BLOCO: Instrumentos de preco (Set A: Z1-Z2) ====
# Z1 e Z2 sao deslocadores de custo (cost shifters), usados para instrumentar
# o preco na equacao de demanda. A logica e que custos afetam precos mas nao
# entram diretamente na utilidade do consumidor (restricao de exclusao).

# Z1: dist (already in df) — miles x diesel index, cost-side shifter
# Z2: coalpost — indicador para produtos ABI + MillerCoors apos a fusao
df <- df %>%
  mutate(coalpost = as.integer((firmid %in% c(1L, 5L)) & (yearid >= 5L)))

cat(sprintf("  Z1 (dist):     range [%.3f, %.3f], %d unique\n",
            min(df$dist), max(df$dist), n_distinct(df$dist)))
cat(sprintf("  Z2 (coalpost): %d of %d obs = 1\n", sum(df$coalpost), nrow(df)))

# ==== BLOCO: Instrumentos de nesting (Set B: Z3-Z8) ====
# Instrumentos estilo BLP: contagem de produtos e soma de caracteristicas dos
# rivais no mercado. Usados para identificar o parametro de nesting (sigma).
# ATENCAO: Z3, Z5 e Z6 sao degenerados nesta amostra pequena (J=8 fixo).

# Z3: num_products per market (DEGENERATE in small data: J=8 always)
df <- df %>%
  group_by(cdid) %>%
  mutate(
    num_products = n(),                              # Z3: num. de produtos no mercado (constante=8)
    sum_dist     = sum(dist)                          # Z4: soma de dist no mercado (variacao cross-section)
  ) %>%
  ungroup() %>%
  mutate(
    # Indicadores de grupo de firma para interacoes
    abi_ind = as.integer(firmid == 1L),
    mc_ind  = as.integer(firmid %in% c(4L, 5L)),     # Miller/Coors coalition
    # Z5-Z8: interacoes (produtos x firma)
    nj_abi     = num_products * abi_ind,              # Z5: degenerado (= 8*ABI sempre)
    nj_mc      = num_products * mc_ind,               # Z6: degenerado (= 8*MC sempre)
    sumdist_abi = sum_dist * abi_ind,                 # Z7: variacao util (dist agregado x ABI)
    sumdist_mc  = sum_dist * mc_ind                   # Z8: variacao util (dist agregado x MC)
  )

cat(sprintf("  Z3 (num_products): %s  %s\n",
            paste(unique(df$num_products), collapse = ","),
            ifelse(n_distinct(df$num_products) == 1, "[DEGENERATE]", "[OK]")))
cat(sprintf("  Z4 (sum_dist):  range [%.3f, %.3f], %d unique\n",
            min(df$sum_dist), max(df$sum_dist), n_distinct(df$sum_dist)))
cat(sprintf("  Z5 (nj_abi):    %s\n", ifelse(n_distinct(df$nj_abi) <= 2, "[DEGENERATE]", "[OK]")))
cat(sprintf("  Z6 (nj_mc):     %s\n", ifelse(n_distinct(df$nj_mc) <= 2, "[DEGENERATE]", "[OK]")))
cat(sprintf("  Z7 (sumdist_abi): range [%.3f, %.3f]\n", min(df$sumdist_abi), max(df$sumdist_abi)))
cat(sprintf("  Z8 (sumdist_mc):  range [%.3f, %.3f]\n", min(df$sumdist_mc), max(df$sumdist_mc)))

# ==== BLOCO: Instrumentos demograficos (Set C: Z9-Z12) ====
# Interacoes entre renda media do mercado e caracteristicas dos produtos.
# No modelo RCNL completo, esses instrumentos identificam os coeficientes
# aleatorios (random coefficients) que capturam heterogeneidade de preferencias.

# Mean income per observation (from dfull, which is 640 x 500 demeaned draws)
# dfull contem 500 draws de renda demeaned; a media das linhas recupera a renda media
df$mean_income <- rowMeans(dfull)

df <- df %>%
  mutate(
    z_inc_const  = mean_income,                       # Z9:  renda media (identifica interacao renda-preco)
    z_inc_calor  = mean_income * calor,               # Z10: renda x calorias
    z_inc_size   = mean_income * sizeid,              # Z11: renda x tamanho da embalagem
    z_inc_import = mean_income * import                # Z12: renda x dummy de importado
  )

cat(sprintf("  Z9  (inc x 1):     range [%.2f, %.2f]\n", min(df$z_inc_const), max(df$z_inc_const)))
cat(sprintf("  Z10 (inc x calor): range [%.2f, %.2f]\n", min(df$z_inc_calor), max(df$z_inc_calor)))
cat(sprintf("  Z11 (inc x size):  range [%.2f, %.2f]\n", min(df$z_inc_size), max(df$z_inc_size)))
cat(sprintf("  Z12 (inc x import):range [%.2f, %.2f]\n", min(df$z_inc_import), max(df$z_inc_import)))

# ==============================================================================
# 2. CONSTRUCT FIXED EFFECTS AND DEMAND/SUPPLY MATRICES
#    Source: main_data.m lines 157-188
# ==============================================================================

# ==== BLOCO: Matrizes de efeitos fixos ====
# Efeitos fixos de produto absorvem a utilidade media de cada marca x tamanho.
# Efeitos fixos de data absorvem choques comuns de tempo (sazonalidade, macro).
# Na oferta, adicionam-se EF de cidade para capturar custos locais fixos.

cat("\n--- Constructing FE matrices ---\n")

# Demand FE: fesd = [prodfe, datefe(:,2:end)]
# x1 = [p_jt, fesd]
# x2 = [p_jt, ones, calor]  (RCNL2 spec)
# Supply FE: fess = [prodfe, cityfe(:,2:end), datefe(:,2:end)]
# w = [mpost, cpost, dist, fess]

n_prods <- n_distinct(df$prodid)  # 8
n_dates <- n_distinct(df$dateid)  # should be ~16 (4 FY x 4 quarters)
n_cities <- n_distinct(df$cityid) # 5

cat(sprintf("  Products: %d, Dates: %d, Cities: %d\n", n_prods, n_dates, n_cities))

# For the cost regression (results_costregs.m):
# X = [apost, mpost, cpost, dist, prodFE, cityFE(:,2:end), dateFE(:,2:end)]
# Verified at line 41: X = [vars.apost vars.w]
# with w = [mpost cpost dist fess] and fess = [prodfe cityfe(:,2:end) datefe(:,2:end)]

# ==============================================================================
# 3. VALIDATION: CROSS-CHECK WITH MATLAB RESULTS
# ==============================================================================

# ==== BLOCO: Validacao cruzada com Matlab ====
# Compara estatisticas descritivas, elasticidades e delta calculados em R
# com os valores originais do Matlab para garantir que os dados foram
# carregados e transformados corretamente.

cat("\n--- Validation: cross-checking with Matlab ---\n")

# 3a. Compare summary statistics with sumstats_pq.txt
# sumstats_pq.txt has 5 rows (firms) x 7 cols
# Format: NaN;NaN;share_12pk;price_12pk;share_24pk;price_24pk;NaN
sumstats <- read.csv2(file.path(path$results, "sumstats_pq.txt"), header = TRUE)
cat("\n  Target (sumstats_pq.txt):\n")
print(sumstats)

cat("\n  Computed from R:\n")
firm_labels <- c("ABI", "Miller", "Corona", "Heineken", "Coors/MC")
# Map firms for summary: ABI=1, Miller (brand 13 under firm 5), Corona=2, Heineken=3, Coors (brand 4)
# The paper groups by brand identity, not firm
summary_r <- df %>%
  mutate(
    brand_label = case_when(
      brndid == 1  ~ "ABI (Bud Light)",
      brndid == 13 ~ "Miller Lite",
      brndid == 5  ~ "Corona",
      brndid == 7  ~ "Heineken",
      brndid == 4  ~ "Coors Light"
    ),
    size_label = ifelse(sizeid == 2, "12pk", "24pk")
  ) %>%
  group_by(brand_label, size_label) %>%
  summarise(
    mean_share = mean(share),
    mean_price = mean(price),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(brand_label, size_label)

print(as.data.frame(summary_r))

# 3b. logodds = ln(s_j) - ln(s_0): variavel dependente do modelo BLP
cat("\n  logodds = ln(s_j) - ln(s_0):\n")
cat(sprintf("    mean=%.4f, sd=%.4f, [%.4f, %.4f]\n",
            mean(df$logodds), sd(df$logodds), min(df$logodds), max(df$logodds)))

# 3c. Compare elasticities with elasMat_2
# First market: own-price elasticities (diagonal)
elas_m1 <- elasMat[,,1]
cat("\n  Elasticity matrix, market 1 (from Matlab):\n")
cat("    Own-price (diagonal):", round(diag(elas_m1), 3), "\n")

# Elasticidades proprias medias: esperamos valores negativos (lei da demanda)
mean_own_elas <- sapply(1:dim(elasMat)[3], function(m) diag(elasMat[,,m]))
cat("    Mean own-price across 140 mkts:", round(rowMeans(mean_own_elas), 3), "\n")

# 3d. Validate delta from daugfile
# delta = utilidade media estimada (mean utility) de cada produto-mercado
cat("\n  Delta validation (from daugfile):\n")
cat(sprintf("    delta_640: mean=%.4f, sd=%.4f, [%.4f, %.4f]\n",
            mean(delta_val), sd(delta_val), min(delta_val), max(delta_val)))

# 3e. Product FE coefficients from daugfile
cat("\n  Product FE coefficients (daugfile.dprodfecoef):\n")
cat("   ", round(dprodfecoef, 4), "\n")

# ==============================================================================
# 4. CONSTRUCT BLP-READY DATASET FOR pyblp
#    If we use pyblp for demand estimation, we need specific column names
# ==============================================================================

# ==== BLOCO: Exportacao do dataset em formato BLP ====
# Renomeia colunas para o padrao exigido pelo pacote pyblp (Python) e salva
# em CSV. Inclui todas as variaveis necessarias: shares, precos, instrumentos,
# caracteristicas e IDs de mercado/produto/firma.

cat("\n--- Preparing BLP-ready CSV ---\n")

df_blp <- df %>%
  mutate(
    market_ids  = cdid,
    product_ids = prodid,
    firm_ids    = firmid,
    shares      = share,
    prices      = price
  ) %>%
  select(
    market_ids, product_ids, firm_ids,
    shares, prices,
    # Product characteristics
    calor, sizeid, import,
    # Market structure
    msize, inshr, outshr, logcondshr, logodds,
    # IDs
    firmid, brndid, cityid, yearid, montid, fiscid, fisccity,
    coalid, leadid,
    # Cost shifters
    mpost, cpost, apost, mcpost, dist, distbutfor,
    # Instruments (demand-side)
    coalpost, num_products, sum_dist,
    sumdist_abi, sumdist_mc,
    z_inc_const, z_inc_calor, z_inc_size, z_inc_import,
    # Income
    mean_income
  )

write_csv(df_blp, file.path(path$data_out, "beer_blp_ready.csv"))
cat(sprintf("  Saved: beer_blp_ready.csv (%d obs, %d vars)\n", nrow(df_blp), ncol(df_blp)))

# Salva objetos R atualizados para uso nos proximos scripts
save(
  df, dfull, cdindex, obsindemand,
  alpha, theta1, theta2, rho, theta2w,
  derMat, elasMat,
  delta_val, deltanp_val, xi_val, mu_val, pcoefi_val,
  dprodfecoef, ddatefecoef, dcityfecoef,
  file = file.path(path$data_out, "step02_output.RData")
)
cat("  Saved: step02_output.RData\n")

# ==============================================================================
# 5. DIAGNOSTIC PLOTS
# ==============================================================================

# ==== BLOCO: Graficos diagnosticos ====
# Tres graficos exploram a estrutura basica dos dados antes da estimacao:
# (1) dispersao share x preco por firma, (2) distribuicao de logodds,
# (3) evolucao temporal das participacoes de mercado por marca.

cat("\n--- Generating diagnostic plots ---\n")

# 5a. Dispersao de shares vs precos — verifica relacao negativa esperada
p1 <- ggplot(df, aes(x = price, y = share, color = factor(firmid))) +
  geom_point(alpha = 0.5) +
  facet_wrap(~sizeid, labeller = labeller(sizeid = c("2" = "12-pack", "3" = "24/30-pack"))) +
  scale_color_manual(
    values = c("1" = "#E41A1C", "2" = "#FF7F00", "3" = "#4DAF4A",
               "4" = "#377EB8", "5" = "#984EA3"),
    labels = c("1" = "ABI", "2" = "Corona", "3" = "Heineken",
               "4" = "Coors", "5" = "MillerCoors")
  ) +
  labs(x = "Price ($/144oz)", y = "Market share", color = "Firm",
       title = "Shares vs. Prices by Firm and Pack Size") +
  theme_minimal()

ggsave(file.path(path$output, "fig_shares_vs_prices.png"), p1, width = 10, height = 5, dpi = 150)
cat("  Saved: fig_shares_vs_prices.png\n")

# 5b. Distribuicao de logodds — variavel dependente do modelo de demanda
p2 <- ggplot(df, aes(x = logodds)) +
  geom_histogram(bins = 25, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = mean(df$logodds), linetype = "dashed", color = "red") +
  labs(x = "log(s_j) - log(s_0)", y = "Count",
       title = "Distribution of BLP Dependent Variable") +
  theme_minimal()

ggsave(file.path(path$output, "fig_logodds_dist.png"), p2, width = 7, height = 4, dpi = 150)
cat("  Saved: fig_logodds_dist.png\n")

# 5c. Evolucao temporal das participacoes — detecta mudancas pos-fusao
p3 <- df %>%
  mutate(
    brand_label = case_when(
      brndid == 1  ~ "Bud Light",
      brndid == 13 ~ "Miller Lite",
      brndid == 5  ~ "Corona",
      brndid == 7  ~ "Heineken",
      brndid == 4  ~ "Coors Light"
    ),
    time_idx = dateid
  ) %>%
  group_by(brand_label, sizeid, time_idx) %>%
  summarise(mean_share = mean(share), .groups = "drop") %>%
  ggplot(aes(x = time_idx, y = mean_share, color = brand_label)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~sizeid, labeller = labeller(sizeid = c("2" = "12-pack", "3" = "24/30-pack")),
             scales = "free_y") +
  labs(x = "Time period", y = "Mean market share", color = "Brand",
       title = "Market Shares Over Time by Brand") +
  theme_minimal()

ggsave(file.path(path$output, "fig_shares_over_time.png"), p3, width = 10, height = 5, dpi = 150)
cat("  Saved: fig_shares_over_time.png\n")

# ==============================================================================
# 6. SUMMARY
# ==============================================================================

cat("\n")
cat("=", rep("=", 59), "\n", sep = "")
cat("STEP 02 COMPLETE\n")
cat("=", rep("=", 59), "\n", sep = "")
cat(sprintf("  Dataset: %d obs, %d variables\n", nrow(df), ncol(df)))
cat(sprintf("  Instruments: 12 total (9 functional in small data)\n"))
cat(sprintf("    Price (Z1-Z2):    2 cost-side shifters\n"))
cat(sprintf("    Nesting (Z3-Z8):  3 functional + 3 degenerate\n"))
cat(sprintf("    RC/demo (Z9-Z12): 4 income x characteristics\n"))
cat(sprintf("  Plots saved to: %s\n", path$output))
cat(sprintf("  Ready for Step 3 (estimation)\n"))
