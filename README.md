# Replicação BLP — Mercado de Cervejas

Replicação em R do artigo **Miller & Weinberg (2017)**, *"Understanding the Price Effects of the MillerCoors Joint Venture"*, Econometrica 85(6): 1763–1791.

## Estrutura

```
Problem Set/
├── R_project/                   ← Implementação em R (scripts, dados, relatório)
│   ├── scripts/                 ← Scripts R numerados (00–05)
│   ├── data/                    ← Dados processados (.csv, .RData)
│   ├── output/                  ← Gráficos diagnósticos
│   ├── documentation markdown/  ← Documentação técnica por etapa
│   ├── problemset_pt.Rmd        ← Relatório final (RMarkdown)
│   └── problemset_pt.html       ← Relatório renderizado
├── replication/                 ← Replication package original (Matlab)
│   ├── data/raw/                ← Dados brutos (.mat)
│   ├── data/analysis/           ← Resultados intermediários (Matlab)
│   ├── code/                    ← Código Matlab original
│   └── results/                 ← Tabelas e figuras do Matlab
├── problemset.pdf               ← Enunciado da lista de exercícios
└── miller and Weinberg (2017).pdf ← Artigo original
```

## Como rodar

Veja o guia completo em [`R_project/documentation markdown/COMO_RODAR_DO_ZERO.md`](R_project/documentation%20markdown/COMO_RODAR_DO_ZERO.md).

Resumo: abra o RStudio, rode os scripts na ordem `00_setup.R` → `05_merger_sim.R`, depois knite `problemset_pt.Rmd`.

## Pipeline

| Script | Etapa | Output |
|--------|-------|--------|
| `00_setup.R` | Pacotes e configuração | — |
| `01_load_data.R` | Carregamento e filtros | 640 obs, 80 mercados |
| `02_clean_shares.R` | Instrumentos e validação | 9 IVs funcionais |
| `03_demand_estimation.R` | Demanda (2SLS + RCNL) | α = −0,027 (2SLS), −0,109 (RCNL) |
| `04_markups_mc.R` | Markups e custos marginais | Lerner médio ≈ 29% |
| `05_merger_sim.R` | Simulação de fusão | Δp fusão ≈ +1,9% |

## Referências

- Miller, N.H. & Weinberg, M.C. (2017). "Understanding the Price Effects of the MillerCoors Joint Venture." *Econometrica*, 85(6): 1763–1791.
- Berry, S., Levinsohn, J. & Pakes, A. (1995). "Automobile Prices in Market Equilibrium." *Econometrica*, 63(4): 841–890.
