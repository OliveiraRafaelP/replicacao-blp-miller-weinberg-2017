# Indice de Scripts
## Miller & Weinberg (2017) — Replicacao em R

### Pipeline de execucao

Cada script carrega os resultados do passo anterior via `.RData` e produz seus proprios outputs. A ordem de execucao é sequencial (00 → 05).

```
00_setup.R ──────────────────► Pacotes, caminhos, constantes e funcoes auxiliares
     │                          (decode_id, make_ownership, compute_shares_rcnl, etc.)
     │
01_load_data.R ──────────────► Carrega .mat, decodifica IDs, constroi dataset de 640 obs
     │                          Saidas: beer_blp.csv, step01_output.RData
     │
02_clean_shares.R ───────────► Constroi instrumentos (Z1-Z12), efeitos fixos, graficos
     │                          Saidas: beer_blp_ready.csv, step02_output.RData
     │
02b_instrument_diagnostics.R ► Degenerescencia, F do 1o estagio, correlacoes parciais
02c_iv_sensitivity.R ────────► Compara baseline (5 IVs) vs full (9 IVs)
02d_iv_diagnosis.R ──────────► Specs A-E, testes de Sargan, selecao de modelo
     │
03_demand_estimation.R ──────► OLS, 2SLS (logit/nested logit), parametros RCNL,
     │                          elasticidades proprias e cruzadas
     │                          Saidas: step03_output.RData
     │
04_markups_mc.R ─────────────► Inversao da CPO de Bertrand, markups, custos marginais
     │                          Matrizes de propriedade, indice de Lerner
     │                          Saidas: step04_output.RData
     │
05_merger_sim.R ─────────────► Contrafactual de fusao (Bertrand pos-fusao)
                                Variacoes de preco, shares e excedente do consumidor
                                Saidas: step05_output.RData
```

### Mapeamento scripts → questoes do problem set

| Secao do Problem Set | Script(s) | Outputs principais |
|---------------------|-----------|--------------------|
| Q1 (Teoria) | — | Respostas textuais no RMarkdown |
| Q2.1 (Carregar dados) | 01_load_data.R | 640 obs, 80 mercados, J=8 |
| Q2.2 (Shares) | 02_clean_shares.R | Shares, outside good, graficos |
| Q2.3 (Instrumentos) | 02_clean_shares.R, 02b | 9 IVs funcionais, justificativa economica |
| Q3.1 (2SLS linear) | 03_demand_estimation.R | Spec B: alpha=-0.059; Spec D: alpha=-0.027 |
| Q3.2 (BLP completo) | 03_demand_estimation.R | RCNL carregado: alpha=-0.109, rho=0.778 |
| Q3.3 (Markups) | 04_markups_mc.R | Markup medio $3.51, Lerner 28.7% |
| Q4.1 (Dataset pos-fusao) | 05_merger_sim.R | Propriedade: firma 4 → firma 5 |
| Q4.2 (Equilibrio Bertrand) | 05_merger_sim.R | 80/80 convergem, dp medio = +1.05% |
| Q4.3 (Tabelas 5/6) | 05_merger_sim.R | Fusao +1.9%, ABI +0.5%, dCS -0.85% |
| Q5 (Bonus) | — | Discussao qualitativa no RMarkdown |

### Documentacao de apoio

| Arquivo | Conteudo |
|---------|---------|
| step1_diagnosis.md | Engenharia reversa da pipeline Matlab, arquitetura de dados, riscos |
| step2_checkpoint.md | Validacao completa: dicionario de dados, instrumentos, sensibilidade IV, amostra pequena vs completa |
| step3_4_estimation.md | Resultados de estimacao, elasticidades, markups, limitacoes |
| step5_post_estimation.md | Elasticidades, markups, validacao contra Matlab |
| step6_counterfactual.md | Resultados da fusao, interpretacao, gap de coordenacao |

### Nota sobre os comentarios nos scripts

Todos os scripts contem comentarios didaticos em portugues, organizados em:

- **Blocos** (`# ==== BLOCO: [nome] ====`): explicam o que o bloco faz e por que e necessario na pipeline.
- **Linhas selecionadas**: comentam variaveis-chave, passos econometricos e formulas.
- **Funcoes**: documentam entradas, saidas e intuicao economica.