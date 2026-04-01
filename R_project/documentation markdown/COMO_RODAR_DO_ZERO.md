# Como Rodar o Projeto do Zero

---

## 1. O que este projeto faz

Este projeto replica em R a análise empírica do artigo **Miller & Weinberg (2017)**, que estuda os efeitos da fusão MillerCoors sobre os preços no mercado de cervejas dos EUA. O trabalho segue cinco etapas: (1) carregar dados de scanner de um pacote de replicação em Matlab; (2) construir variáveis de mercado e instrumentos; (3) estimar a demanda por cerveja usando variáveis instrumentais; (4) recuperar custos marginais e markups; (5) simular o efeito da fusão sobre preços e bem-estar dos consumidores. Ao final, um documento RMarkdown (`problemset_pt.Rmd`) compila tudo em um relatório HTML com tabelas, gráficos e respostas às questões da lista de exercícios do curso "Organização Industrial Empírica", do Prof. Claudio Lucinda, ministrado em 2026.

---

## 2. O que você precisa ter antes de começar

### Arquivos obrigatórios

Você precisa da seguinte estrutura de pastas no seu computador:

```
Problem Set/
├── replication/
│   ├── data/
│   │   ├── raw/
│   │   │   ├── small_scanner.mat        ← dados de scanner (preços, quantidades)
│   │   │   ├── small_demosE.mat         ← sorteios de renda dos consumidores
│   │   │   └── small_dresgmm2.mat       ← parâmetros de demanda pré-estimados
│   │   └── analysis/
│   │       ├── daugfile.mat             ← dados intermediários (para validação)
│   │       └── sres_bertrand.mat        ← custos marginais do Matlab (para validação)
│   └── results/
│       ├── sumstats_pq.txt             ← estatísticas descritivas do Matlab
│       ├── cfmerger.txt                ← resultados de fusão do Matlab
│       └── (outros .txt)
└── R_project/
    ├── scripts/                         ← todos os scripts R (00 a 05)
    ├── data/                            ← dados processados (gerados pelos scripts)
    ├── output/                          ← gráficos (gerados pelos scripts)
    └── problemset_pt.Rmd               ← relatório final
```
A pasta Problem Set também contém o artigo **Miller & Weinberg (2017)** e a lista de exercícios do curso (problemset.pdf).

Os três arquivos `.mat` dentro de `replication/data/raw/` são **obrigatórios**. Eles vêm do zip do replication package (arquivo `140341-V2.zip`). Se você ainda não extraiu esse zip, faça isso primeiro.

### Arquivos gerados automaticamente

Os scripts criam os seguintes arquivos dentro de `R_project/data/`:

| Arquivo | Criado por | Conteúdo |
|---------|-----------|----------|
| `beer_blp.csv` | `01_load_data.R` | Dataset principal (640 obs) |
| `beer_blp_ready.csv` | `02_clean_shares.R` | Dataset com instrumentos |
| `step01_output.RData` | `01_load_data.R` | Objetos R do passo 1 |
| `step02_output.RData` | `02_clean_shares.R` | Objetos R do passo 2 |
| `step03_output.RData` | `03_demand_estimation.R` | Estimativas de demanda |
| `step04_output.RData` | `04_markups_mc.R` | Markups e custos marginais |
| `step05_output.RData` | `05_merger_sim.R` | Resultados da simulação de fusão |

Esses arquivos **não precisam existir antes** — eles são criados quando você roda os scripts.

---

## 3. O que abrir

1. **Abra o RStudio** (ou outro editor com suporte a R).
2. **Defina o diretório de trabalho** como a pasta `Problem Set/`. No RStudio: `Session > Set Working Directory > Choose Directory...` e selecione a pasta `Problem Set`.
3. **Olhe primeiro o arquivo** `R_project/SCRIPTS.md` — ele contém o mapa de todos os scripts e a ordem de execução.
4. **Para rodar cada script**, abra o arquivo `.R` correspondente no RStudio e clique em `Source` (ou pressione `Ctrl+Shift+S`).

**Importante**: os scripts devem ser rodados **na ordem** (00, 01, 02, 02b, 02c, 02d, 03, 04, 05). Cada um depende dos resultados do anterior.

---

## 4. Ordem de execução

### Passo 0: `00_setup.R`

- **O que faz**: instala e carrega todos os pacotes necessários (tidyverse, ivreg, nleqslv, etc.), define os caminhos das pastas e cria funções auxiliares usadas nos demais scripts.
- **O que cria**: nada em disco; apenas carrega objetos na memória do R.
- **O que verificar**: o console deve imprimir `"Setup complete."` sem erros. Se aparecer erro de pacote, veja a Seção 5.

### Passo 1: `01_load_data.R`

- **O que faz**: carrega os arquivos `.mat` do Matlab, decodifica os identificadores de produto, aplica os filtros temporais (exclusão do período da fusão), calcula market shares e carrega os parâmetros de demanda do artigo.
- **O que cria**: `data/beer_blp.csv` e `data/step01_output.RData`.
- **O que verificar**: o console deve imprimir:
  - `Raw scanner: 1120 obs x 9 cols`
  - `After fiscal filter: 640 obs pass`
  - `cdid: 80 markets`
  - `J=8 in all markets: CONFIRMED`
  - `STEP 01 COMPLETE`

### Passo 2: `02_clean_shares.R`

- **O que faz**: constrói os 12 instrumentos de demanda, monta matrizes de efeitos fixos, compara estatísticas com o Matlab original e gera gráficos diagnósticos.
- **O que cria**: `data/beer_blp_ready.csv`, `data/step02_output.RData`, e três gráficos em `output/`.
- **O que verificar**: o console deve listar os instrumentos Z1-Z12 e indicar quais são `[OK]` e quais são `[DEGENERATE]`. Deve terminar com `STEP 02 COMPLETE`.

### Passo 2b: `02b_instrument_diagnostics.R`

- **O que faz**: verifica a qualidade dos instrumentos — quais têm variação suficiente, qual a correlação com as variáveis endógenas, qual a estatística F do primeiro estágio.
- **O que cria**: nada em disco (apenas imprime diagnósticos no console).
- **O que verificar**: deve reportar 9 instrumentos funcionais e 3 degenerados. A estatística F do primeiro estágio deve ser > 10 para ambas as variáveis endógenas.

### Passo 2c: `02c_iv_sensitivity.R`

- **O que faz**: compara a estimação com 5 instrumentos (baseline) vs. 9 instrumentos (incluindo renda). Testa se os resultados são estáveis.
- **O que cria**: nada em disco.
- **O que verificar**: o coeficiente de preço deve ser negativo com o baseline (5 IVs) mas pode ficar positivo (sinal errado) com o conjunto completo (9 IVs). Isso é esperado — os dados são perturbados.

### Passo 2d: `02d_iv_diagnosis.R`

- **O que faz**: testa várias especificações (Specs A-E) para encontrar a mais confiável. Recomenda a Spec D (nested logit, 5 IVs).
- **O que cria**: nada em disco.
- **O que verificar**: a Spec D deve ser a única com alpha < 0 **e** sigma entre 0 e 1 simultaneamente.

### Passo 3: `03_demand_estimation.R`

- **O que faz**: estima a demanda por OLS e 2SLS (logit e nested logit), carrega os parâmetros RCNL do artigo, calcula elasticidades-preço próprias e cruzadas.
- **O que cria**: `data/step03_output.RData`.
- **O que verificar**: a tabela de comparação deve mostrar alpha negativo nas especificações IV. A elasticidade-preço própria média (RCNL) deve ser próxima de -4,7.

### Passo 4: `04_markups_mc.R`

- **O que faz**: inverte as condições de primeira ordem de Bertrand-Nash para recuperar markups e custos marginais, usando os parâmetros RCNL do artigo.
- **O que cria**: `data/step04_output.RData`.
- **O que verificar**:
  - Todos os custos marginais devem ser positivos (`Negative MC: R=0`).
  - A correlação com o Matlab deve ser > 0,95.
  - O markup médio deve ser aproximadamente $3,50.
- **Atenção**: este script pode levar **1-2 minutos** para rodar (80 mercados × simulação com 500 consumidores).

### Passo 5: `05_merger_sim.R`

- **O que faz**: simula o contrafactual da fusão MillerCoors. Altera a matriz de propriedade (Coors passa a pertencer à MillerCoors) e resolve o novo equilíbrio de Bertrand para os 80 mercados.
- **O que cria**: `data/step05_output.RData`.
- **O que verificar**:
  - Todos os 80 mercados devem convergir (`Converged: 80/80`).
  - As partes fusionadas devem ter aumento de preço positivo (~+1,9%).
  - A ABI deve ter aumento menor (~+0,5%).
  - Importados devem ter variação próxima de zero.
- **Atenção**: este script pode levar **2-5 minutos**.

### Passo final: `problemset_pt.Rmd`

- **O que faz**: compila todas as respostas, tabelas e gráficos em um relatório HTML.
- **Como rodar**: no RStudio, abra o arquivo e clique em `Knit` (botão com ícone de agulha de tricô no topo do editor). Ou rode no console: `rmarkdown::render("R_project/problemset_pt.Rmd")`.
- **O que cria**: `problemset_pt.html` (o arquivo final para entrega).
- **O que verificar**: o HTML deve abrir no navegador com todas as seções, tabelas e gráficos renderizados corretamente.

**Nota sobre o Pandoc**: se o RStudio reportar erro de Pandoc, o `Knit` pode não funcionar fora do RStudio. Dentro do RStudio, o Pandoc já está incluído e não há problema.

---

## 5. O que fazer se der erro

### "Pacote não encontrado" / "there is no package called 'X'"

Rode no console do R:
```r
install.packages("X")
```
Substitua `X` pelo nome do pacote que falta. Os pacotes necessários são: `R.matlab`, `tidyverse`, `ivreg`, `nleqslv`, `Matrix`, `sandwich`, `lmtest`, `knitr`, `kableExtra`.

### "Caminho não encontrado" / "cannot open file"

Os scripts usam caminhos absolutos que começam com `C:/Users/rafael.oliveira/Desktop/Problem Set/`. Se a sua pasta está em outro lugar, você precisa editar os caminhos no arquivo `00_setup.R` (na seção `path <- list(...)`) e no início do `problemset_pt.Rmd` (no bloco `setup`).

### "Arquivo .mat não encontrado"

Verifique que os três arquivos `.mat` estão dentro de `replication/data/raw/`. Se você não extraiu o zip `140341-V2.zip`, faça isso na pasta `Problem Set/`.

### "Objeto não encontrado" / "object 'df' not found"

Isso significa que você pulou um script anterior. Os scripts devem ser rodados **na ordem**: 00, 01, 02, ..., 05. Cada um cria objetos que o próximo precisa. Volte ao passo que faltou e rode novamente.

### "Script demora muito"

Os scripts `04_markups_mc.R` e `05_merger_sim.R` fazem simulações com 500 consumidores para 80 mercados. É normal levarem 1-5 minutos. Se travar por mais de 10 minutos, interrompa (`Esc` ou botão vermelho no RStudio) e verifique se há erro nos dados.

### "Pandoc not found" (ao knitar o Rmd)

Use o RStudio para knitar — ele inclui o Pandoc automaticamente. Se estiver usando R puro no terminal, instale o Pandoc separadamente: https://pandoc.org/installing.html

---

## 6. Checklist final

Após rodar todos os scripts e knitar o Rmd, verifique:

- [ ] A pasta `R_project/data/` contém 5 arquivos `.RData` (step01 a step05) e 2 `.csv`
- [ ] A pasta `R_project/output/` contém 3 arquivos `.png` (gráficos)
- [ ] O arquivo `problemset_pt.html` existe e abre no navegador
- [ ] No HTML, todas as tabelas mostram números (não `NA` ou `NaN`)
- [ ] O gráfico de shares vs. preços aparece corretamente
- [ ] A tabela de estimação mostra alpha negativo para Spec B e Spec D
- [ ] A tabela de markups mostra valores positivos para todas as firmas
- [ ] A tabela de fusão mostra aumento de preço para as partes fusionadas
- [ ] O relatório termina com a seção "Apêndice: Códigos Principais"

---

## 7. Glossário rápido

| Termo | O que significa |
|-------|----------------|
| **Share (market share)** | Participação de mercado de um produto. Calculada como quantidade vendida dividida pelo tamanho total do mercado. Exemplo: se Bud Light vende 10.000 unidades num mercado de 300.000, seu share é ~3,3%. |
| **Outside share** | A fração dos consumidores que **não compra nenhum** dos produtos modelados. Representa substitutos fora do modelo (outras bebidas, não consumir, etc.). É calculada como 1 menos a soma de todos os shares. |
| **Instrumento (IV)** | Variável que afeta o preço (por exemplo, custo de transporte) mas não afeta diretamente a preferência do consumidor. Usada para corrigir o problema de endogeneidade na estimação de demanda. |
| **Efeito fixo (FE)** | Variável dummy que absorve diferenças médias entre grupos. Efeito fixo de produto captura que Bud Light é diferente de Corona em características fixas; efeito fixo de tempo captura que um trimestre é diferente de outro. |
| **Markup** | Diferença entre o preço cobrado e o custo marginal de produção. Mede o poder de mercado da firma. Quanto maior o markup, maior a margem de lucro. |
| **Índice de Lerner** | Markup dividido pelo preço: `(preço - custo) / preço`. Varia entre 0 (concorrência perfeita) e 1 (monopólio). |
| **Custo marginal (MC)** | Custo de produzir uma unidade adicional do produto. Não é observado diretamente — é recuperado invertendo as condições de equilíbrio do modelo de demanda. |
| **Consumer surplus (CS)** | Excedente do consumidor: medida de bem-estar que captura o quanto os consumidores ganham por poderem comprar os produtos aos preços de mercado em vez de não comprar nada. |
| **FOC (first-order condition)** | Condição de primeira ordem. É a equação que descreve o preço ótimo de cada firma no equilíbrio de Nash. No modelo de Bertrand: `preço = custo marginal + markup`. |
| **Bertrand-Nash** | Modelo de concorrência em preços entre firmas que vendem produtos diferenciados. Cada firma escolhe seu preço supondo que os rivais mantêm os deles fixos. |
| **Nested logit / RCNL** | Modelos de demanda que permitem padrões de substituição mais realistas que o logit simples. O nested logit agrupa produtos em "ninhos"; o RCNL adiciona heterogeneidade de preferências entre consumidores. |
