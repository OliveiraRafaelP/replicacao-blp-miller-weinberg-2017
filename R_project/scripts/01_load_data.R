# ==============================================================================
# 01_load_data.R — Load .mat files, decode IDs, build raw tibble
# Miller & Weinberg (2017) Replication in R
# ==============================================================================
#
# Este script carrega os dados brutos do scanner (formato .mat do Matlab),
# decodifica os identificadores compostos, aplica os filtros temporais
# (exclusao do periodo ao redor da fusao), calcula market shares e
# carrega parametros de demanda previamente estimados.
#
# Fluxo principal:
#   1. Carregar small_scanner.mat (1120 obs)
#   2. Decodificar IDs compostos em firma, marca, tamanho, cidade, ano, mes
#   3. Montar tibble bruto (1120 obs)
#   4. Filtro obsintemp: excluir obs proximas a fusao -> 960 obs
#   5. Filtro fiscal: manter apenas anos fiscais de interesse -> 640 obs
#   6. Construir dataset final com IDs derivados
#   7. Calcular market shares e log-odds (variavel dependente do BLP)
#   8-13. Carregar parametros de demanda, demograficos e daugfile (validacao)
#   14. Salvar dados processados em CSV e RData
# ==============================================================================

source("C:/Users/rafael.oliveira/Desktop/Organização Industrial Empírica/Problem Set/R_project/scripts/00_setup.R")

# ==== BLOCO: Carregar dados brutos do scanner ====
# Arquivo small_scanner.mat contem 1120 observacoes (produto x mercado x mes)
# com 9 colunas. Esses dados vem do IMS Health scanner data para cerveja.

cat("\n--- Loading small_scanner.mat ---\n")
mat_scanner <- R.matlab::readMat(file.path(path$data_raw, "small_scanner.mat"))
scanner     <- mat_scanner$small.scanner  # 1120 x 9

cat(sprintf("  Raw scanner: %d obs x %d cols\n", nrow(scanner), ncol(scanner)))
stopifnot(nrow(scanner) == 1120, ncol(scanner) == 9)

# ==== BLOCO: Decodificar ID composto ====
# A coluna 1 do scanner contem um ID numerico de 11 digitos que empacota
# 6 dimensoes: firma, marca, tamanho, cidade, ano e mes.
# A decodificacao extrai cada dimensao por divisao inteira sucessiva.
# Referencia: main_data.m linhas 34-40 do codigo Matlab original.

cat("--- Decoding composite IDs ---\n")

id2    <- scanner[, 1]
firmid <- floor(id2 / 1e10)
brndid <- floor((id2 - firmid * 1e10) / 1e8)
sizeid <- floor((id2 - firmid * 1e10 - brndid * 1e8) / 1e6)
cityid <- floor((id2 - firmid * 1e10 - brndid * 1e8 - sizeid * 1e6) / 1e4)
yearid <- floor((id2 - firmid * 1e10 - brndid * 1e8 - sizeid * 1e6 - cityid * 1e4) / 1e2)
montid <- round(id2 - firmid * 1e10 - brndid * 1e8 - sizeid * 1e6 - cityid * 1e4 - yearid * 1e2)

# ==== BLOCO: Montar tibble bruto (1120 obs) ====
# Mapeia cada coluna do scanner para uma variavel nomeada.
# price: preco em dolares por equivalente de 144oz
# quantity: vendas em unidades de 144oz
# msize: tamanho do mercado potencial (consumidores)
# miles/dist: variaveis de distancia usadas como instrumentos de custo

df_raw <- tibble(
  id2       = id2,
  firmid    = as.integer(firmid),
  brndid    = as.integer(brndid),
  sizeid    = as.integer(sizeid),
  cityid    = as.integer(cityid),
  yearid    = as.integer(yearid),
  montid    = as.integer(montid),
  price     = scanner[, 2],          # p_jt: $/144oz equiv (line 133)
  quantity  = scanner[, 4],          # q_jt: unit sales, 144oz equiv (line 136)
  msize     = scanner[, 7],          # market size, potential consumers (line 139)
  miles     = scanner[, 5],          # miles to brewery, raw (line 148)
  calor_raw = scanner[, 6],          # calories, raw (line 153: /100 then demean)
  dist      = scanner[, 8],          # distance to ABI brewery (miles x diesel) (line 149)
  distbutfor= scanner[, 9]           # distance, but-for counterfactual (line 150)
)

# yearcityid: identificador unico de combinacao ano-cidade.
# Necessario para expandir os dados demograficos (35 combos = 7 anos x 5 cidades).
# Replica grp2idx(yearid*100+cityid) do Matlab — main_data.m linha 42.
df_raw <- df_raw %>%
  mutate(yc_combo = yearid * 100L + cityid) %>%
  mutate(yearcityid = as.integer(factor(yc_combo, levels = sort(unique(yc_combo))))) %>%
  select(-yc_combo)

cat(sprintf("  yearcityid: %d unique combos (expect 35 = 7yr x 5city)\n",
            n_distinct(df_raw$yearcityid)))
stopifnot(n_distinct(df_raw$yearcityid) == 35)

# ==== BLOCO: Filtro obsintemp (exclusao do periodo da fusao) ====
# A fusao MillerCoors foi anunciada/aprovada entre os anos 4 e 5 da amostra.
# O filtro RCNL2 exclui observacoes proximas a esse evento para evitar
# contaminacao por efeitos de transicao.
# Mantem: yearid<=3 (pre), yearid>=6 (pos), e meses nas bordas.
# Resultado: 960 de 1120 obs passam no filtro.

cat("--- Applying obsintemp filter (RCNL2) ---\n")

df_raw <- df_raw %>%
  mutate(
    obsintemp = (yearid <= 3L | yearid >= 6L) |
                (yearid == 4L & montid <= 2L) |
                (yearid == 5L & montid >= 3L)
  )

cat(sprintf("  obsintemp: %d of %d obs pass (expect 960)\n",
            sum(df_raw$obsintemp), nrow(df_raw)))
stopifnot(sum(df_raw$obsintemp) == 960)

# ==== BLOCO: Ano fiscal e filtro final ====
# O ano fiscal comeca no Q4 (mes >= 4 implica proximo ano fiscal).
# Formula: fiscid = yearid + 1*(montid >= 4) + 2004
# Depois, seleciona apenas os 4 anos fiscais de interesse:
#   2006, 2007 (pre-fusao) e 2010, 2011 (pos-fusao).
# Resultado final: 640 observacoes.

cat("--- Computing fiscal year and applying fiscal filter ---\n")

df_raw <- df_raw %>%
  mutate(
    fiscid = yearid + as.integer(montid >= 4L) + 2004L,
    obsin  = obsintemp & (fiscid %in% c(2006L, 2007L, 2010L, 2011L))
  )

cat(sprintf("  After fiscal filter: %d obs pass (expect 640)\n", sum(df_raw$obsin)))
stopifnot(sum(df_raw$obsin) == 640)

# obsindemand: vetor logico que mapeia as 960 obs (obsintemp) para as 640 finais.
# Usado para indexar arrays intermediarios do daugfile (que tem 960 linhas).
obsindemand <- df_raw$obsin[df_raw$obsintemp]
cat(sprintf("  obsindemand: %d of %d selected from daugfile\n",
            sum(obsindemand), length(obsindemand)))
stopifnot(sum(obsindemand) == 640, length(obsindemand) == 960)

# ==== BLOCO: Dataset filtrado (640 obs) ====
# Filtra e cria identificadores derivados para a estimacao.

cat("--- Building filtered dataset ---\n")

df <- df_raw %>%
  filter(obsin) %>%
  select(-obsintemp, -obsin)

# IDs derivados (main_data.m linhas 91-113):
# prodid   — identifica produto (combinacao marca x tamanho)
# dateid   — identifica periodo temporal (ano x mes)
# fisccity — efeito fixo cidade x ano fiscal
# coalid   — dummy para a coalizao pos-fusao (ABI + Coors + MillerCoors)
# leadid   — dummy para a firma lider (ABI / Anheuser-Busch)
# import   — dummy para marcas importadas
df <- df %>%
  mutate(
    prodid   = as.integer(factor(brndid * 100L + sizeid)),  # grp2idx(brndid*100+sizeid)
    dateid   = as.integer(factor(yearid * 100L + montid)),  # grp2idx(yearid*100+montid)
    fisccity = as.integer(factor(cityid * 100L + fiscid)),  # grp2idx(cityid*100+fiscid)
    coalid   = as.integer(firmid %in% c(1L, 4L, 5L)),      # Coalition: ABI + Coors + MC
    leadid   = as.integer(firmid == 1L),                     # Leader: ABI
    import   = as.integer(firmid %in% c(2L, 3L))            # Import brands
  )

# ==== BLOCO: Construcao do cdid (identificador sequencial de mercado) ====
# Cada mercado e definido pela combinacao unica de cidade x ano x mes.
# O cdid numera sequencialmente os mercados na ordem em que aparecem nos dados.
# Espera-se 80 mercados = 5 cidades x 16 periodos.

# Extrai os digitos 6 a 11 do ID composto (cidade + ano + mes) como string
id2_str <- sprintf("%.0f", df$id2)
mkt_code <- substr(id2_str, 6, 11)

# Incrementa cdid a cada mudanca de mercado
df$cdid <- 1L
for (i in 2:nrow(df)) {
  df$cdid[i] <- df$cdid[i - 1L] + as.integer(mkt_code[i] != mkt_code[i - 1L])
}

cat(sprintf("  cdid: %d markets (expect 80)\n", max(df$cdid)))
stopifnot(max(df$cdid) == 80)

# cdindex: indice da ultima linha de cada mercado (util para loops por mercado)
cdindex <- which(c(diff(df$cdid) != 0, TRUE))
stopifnot(length(cdindex) == 80)

# Verifica J=8 produtos em todos os mercados (painel balanceado)
prods_per_mkt <- table(df$cdid)
stopifnot(all(prods_per_mkt == 8))
cat("  J=8 in all markets: CONFIRMED\n")

# ==== BLOCO: Market shares e opcao externa (outside good) ====
# s_jt = q_jt / msize: participacao de mercado do produto j no periodo t.
# inshr: soma dos shares de todos os produtos dentro do mercado (inside share).
# outshr: 1 - inshr = share da opcao externa (nao comprar cerveja).
# logodds: ln(s_j) - ln(s_0) = variavel dependente na equacao de demanda BLP.
#   No caso do logit simples, logodds = delta (utilidade media).
# logcondshr: share condicional dentro do ninho (relevante para RCNL).

cat("--- Computing shares ---\n")

df <- df %>%
  mutate(share = quantity / msize) %>%
  group_by(cdid) %>%
  mutate(inshr = sum(share)) %>%
  ungroup() %>%
  mutate(
    outshr     = 1.0 - inshr,
    logcondshr = log(share) - log(inshr),  # main_data.m line 200
    logodds    = log(share) - log(outshr)   # main_data.m line 204 = ln(s_j) - ln(s_0)
  )

# Validacao: shares devem ser positivos e finitos
stopifnot(all(df$share > 0))
stopifnot(all(df$outshr > 0))
stopifnot(all(is.finite(df$logodds)))
cat(sprintf("  share:  [%.6f, %.6f]  all > 0: PASS\n", min(df$share), max(df$share)))
cat(sprintf("  inshr:  [%.4f, %.4f]\n", min(df$inshr), max(df$inshr)))
cat(sprintf("  outshr: [%.4f, %.4f]  all > 0: PASS\n", min(df$outshr), max(df$outshr)))
cat(sprintf("  logodds: [%.4f, %.4f]  all finite: PASS\n", min(df$logodds), max(df$logodds)))

# ==== BLOCO: Calorias centralizadas ====
# Divide por 100 e subtrai a media amostral (640 obs).
# A centralizacao (demeaning) e padrao em IO empirica para que o intercepto
# capture o efeito medio e o coeficiente de calorias capture desvios.

df <- df %>%
  mutate(calor = calor_raw / 100 - mean(calor_raw / 100))

# ==== BLOCO: Cost shifters e dummies pos-fusao ====
# Variaveis que capturam mudancas estruturais apos a fusao MillerCoors.
# mcpost: MillerCoors (firmid=5) no periodo pos-fusao
# mpost:  marcas Miller no periodo pos-fusao
# cpost:  marcas Coors no periodo pos-fusao
# apost:  ABI (Anheuser-Busch) no periodo pos-fusao
# Essas dummies permitem estimar mudancas de custo marginal pos-fusao.

df <- df %>%
  mutate(
    mcpost = as.integer(firmid == 5L & yearid >= 5L),                            # line 184
    mpost  = as.integer((brndid %in% c(11L, 12L, 13L)) & (yearid >= 5L)),       # line 185
    cpost  = as.integer((brndid %in% c(3L, 4L)) * (yearid >= 5L)),              # line 186
    apost  = as.integer(firmid == 1L & yearid >= 5L)                              # line 187
  )

# ==== BLOCO: Carregar parametros de demanda estimados ====
# Arquivo small_dresgmm2.mat contem os resultados da estimacao GMM do modelo
# de demanda RCNL, previamente obtidos no Matlab.
# alpha:   coeficiente medio de preco (espera-se negativo)
# theta1:  todos os parametros lineares (alpha + 62 efeitos fixos)
# theta2:  parametros de coeficientes aleatorios (interacao com renda)
# rho:     parametro de nesting (0 = logit, 1 = correlacao perfeita no ninho)
# derMat:  matrizes 8x8 de derivadas de preco por mercado (para calculo de markup)
# elasMat: matrizes 8x8 de elasticidades-preco por mercado

cat("\n--- Loading demand parameters ---\n")

mat_dres <- R.matlab::readMat(file.path(path$data_raw, "small_dresgmm2.mat"))

alpha   <- mat_dres$theta1.2[1, 1]       # Mean price coefficient
theta1  <- mat_dres$theta1.2[, 1]        # All linear params (alpha + 62 FEs)
theta2  <- mat_dres$theta2.2[, 1]        # RC params (3x1)
rho     <- mat_dres$rho.2[1, 1]          # Nesting parameter
derMat  <- mat_dres$derMat.2             # (8, 8, 140) derivative matrices
elasMat <- mat_dres$elasMat.2            # (8, 8, 140) elasticity matrices

# theta2w: matriz 3x2 que mapeia coeficientes aleatorios as demograficas.
# Coluna 2 contem interacoes com renda (income):
#   linha 1: preco x renda (consumidores mais ricos sao menos sensiveis a preco)
#   linha 2: constante x renda
#   linha 3: calorias x renda
theta2w <- matrix(0, nrow = 3, ncol = 2)
theta2w[1, 2] <- theta2[1]  # price x income
theta2w[2, 2] <- theta2[2]  # constant x income
theta2w[3, 2] <- theta2[3]  # calories x income

cat(sprintf("  alpha = %.6f\n", alpha))
cat(sprintf("  rho   = %.6f\n", rho))
cat(sprintf("  theta2 = [%.5f, %.5f, %.5f]\n", theta2[1], theta2[2], theta2[3]))
cat(sprintf("  derMat:  %s\n", paste(dim(derMat), collapse = " x ")))
cat(sprintf("  elasMat: %s\n", paste(dim(elasMat), collapse = " x ")))

# Validacao contra valores conhecidos da estimacao original
stopifnot(abs(alpha - (-0.10872)) < 1e-4)
stopifnot(abs(rho - 0.77788) < 1e-4)
stopifnot(all(dim(derMat) == c(8, 8, 140)))

# ==== BLOCO: Carregar demograficos (renda simulada) ====
# Arquivo small_demosE.mat contem 35 x 1500 draws de renda.
# As 35 linhas correspondem as 35 combinacoes ano-cidade (yearcityid).
# As primeiras 500 colunas sao os draws de renda usados na simulacao (NS=500).
# Expansao: cada observacao recebe os draws de renda correspondentes
# ao seu yearcityid (calculado sobre todas as 1120 obs antes de filtrar).
# Depois, centraliza-se (demean) a renda para interpretacao dos coeficientes.

cat("\n--- Loading demographics ---\n")

mat_demo <- R.matlab::readMat(file.path(path$data_raw, "small_demosE.mat"))
demosE   <- mat_demo$small.demosE  # 35 x 1500

cat(sprintf("  demosE: %d x %d\n", nrow(demosE), ncol(demosE)))
stopifnot(nrow(demosE) == 35, ncol(demosE) >= 1500)

# Expande para 1120 obs usando yearcityid, depois filtra para 640
demos_expanded <- demosE[df_raw$yearcityid, 1:NS]  # 1120 x 500
demos_filtered <- demos_expanded[df_raw$obsintemp & (df_raw$fiscid %in% FISC_YEARS), ]  # 640 x 500

# dfull: matriz 640 x 500 de renda centralizada (cada coluna = um consumidor simulado)
dfull <- demos_filtered - mean(demos_filtered)

cat(sprintf("  dfull: %d x %d (demeaned income draws)\n", nrow(dfull), ncol(dfull)))
cat(sprintf("  mean(dfull) = %.8f (expect ~0)\n", mean(dfull)))
stopifnot(nrow(dfull) == 640, ncol(dfull) == NS)
stopifnot(abs(mean(dfull)) < 1e-10)

# ==== BLOCO: Carregar daugfile (arquivo de validacao) ====
# O daugfile.mat e um arquivo intermediario gerado pelo Matlab contendo
# resultados da estimacao de demanda: delta (utilidade media), mu
# (heterogeneidade), pcoefi (coeficientes individuais de preco), etc.
# Usado aqui exclusivamente para validacao cruzada dos calculos em R.
# As 960 linhas correspondem ao filtro obsintemp; selecionamos 640 via obsindemand.

cat("\n--- Loading daugfile (validation target) ---\n")

mat_daug <- R.matlab::readMat(file.path(path$data_ana, "daugfile.mat"))
daug     <- mat_daug$daugfile[,,1]  # R.matlab reads struct as 3D array with named dims

# Extrair campos-chave (960 linhas cada)
delta_960   <- as.numeric(daug$delta)
mu_960      <- daug$mu             # 960 x 500
pcoefi_960  <- daug$pcoefi         # 960 x 500
deltanp_960 <- as.numeric(daug$deltanp)
xi_960      <- as.numeric(daug$xi)

# Selecionar 640 obs via obsindemand
delta_val   <- delta_960[obsindemand]
deltanp_val <- deltanp_960[obsindemand]
xi_val      <- xi_960[obsindemand]
mu_val      <- mu_960[obsindemand, ]
pcoefi_val  <- pcoefi_960[obsindemand, ]

cat(sprintf("  delta (640):   [%.4f, %.4f]\n", min(delta_val), max(delta_val)))
cat(sprintf("  deltanp (640): [%.4f, %.4f]\n", min(deltanp_val), max(deltanp_val)))

# Coeficientes de efeitos fixos para validacao posterior
dprodfecoef <- as.numeric(daug$dprodfecoef)  # 8 x 1
ddatefecoef <- as.numeric(daug$ddatefecoef)  # 24 x 1
dcityfecoef <- as.numeric(daug$dcityfecoef)  # 5 x 1

cat(sprintf("  dprodfecoef: %d values\n", length(dprodfecoef)))
cat(sprintf("  dcityfecoef: %s\n", paste(round(dcityfecoef, 4), collapse = ", ")))

# ==== BLOCO: Salvar dados processados ====
# Exporta o dataset final em dois formatos:
#   - CSV: para inspecao e uso externo
#   - RData: para carregamento rapido nos scripts subsequentes,
#     incluindo todos os objetos necessarios (dados, parametros, validacao)

cat("\n--- Saving processed data ---\n")

# Main dataset
write_csv(df, file.path(path$data_out, "beer_blp.csv"))
cat(sprintf("  Saved: beer_blp.csv (%d obs)\n", nrow(df)))

# Save R objects for downstream scripts
save(
  df, dfull, cdindex, obsindemand,
  alpha, theta1, theta2, rho, theta2w,
  derMat, elasMat,
  delta_val, deltanp_val, xi_val, mu_val, pcoefi_val,
  dprodfecoef, ddatefecoef, dcityfecoef,
  file = file.path(path$data_out, "step01_output.RData")
)
cat("  Saved: step01_output.RData\n")

# ==== BLOCO: Resumo final ====

cat("\n")
cat("=" , rep("=", 59), "\n", sep = "")
cat("STEP 01 COMPLETE: Data loaded and processed\n")
cat("=" , rep("=", 59), "\n", sep = "")
cat(sprintf("  Observations:  %d (8 products x 80 markets)\n", nrow(df)))
cat(sprintf("  Firms:         %s\n", paste(sort(unique(df$firmid)), collapse = ", ")))
cat(sprintf("  Brands:        %s\n", paste(sort(unique(df$brndid)), collapse = ", ")))
cat(sprintf("  Cities:        %d\n", n_distinct(df$cityid)))
cat(sprintf("  Fiscal years:  %s\n", paste(sort(unique(df$fiscid)), collapse = ", ")))
cat(sprintf("  alpha = %.5f, rho = %.5f\n", alpha, rho))
cat(sprintf("  Demographics:  %d x %d (demeaned income draws)\n", nrow(dfull), ncol(dfull)))
