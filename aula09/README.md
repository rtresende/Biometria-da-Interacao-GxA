# Aula 09 — Modelos Mistos IV: Regressão Aleatória e Normas de Reação

**Data:** 06 de maio de 2026  
**Carga horária:** 4h

## Conteúdo abordado

- Recapitulação da ANOVA conjunta, Finlay–Wilkinson, Eberhart–Russell e modelos mistos para G×A
- Construção do índice ambiental a partir da média fenotípica dos ambientes
- Transição da regressão de Finlay–Wilkinson para a regressão aleatória
- Intercepto e inclinação como efeitos aleatórios de genótipos
- Construção e interpretação da matriz `Z` na regressão aleatória
- Matriz de variâncias e covariância dos coeficientes aleatórios, `G_u`
- Variância genética ao longo do gradiente ambiental
- Covariância genética entre ambientes induzida pelo gradiente
- Normas de reação lineares e interpretação de responsividade genotípica
- Cuidados na interpretação do índice ambiental e da extrapolação

## Material

- 📄 `Biometria_da_Interação_GxA_GMP0164_aula9.pdf` — roteiro teórico da aula
- 💾 `script_GMP0164_aula9.R` — script utilizado na aula prática
- 📊 `maize_safrinha_MET.txt` — conjunto de dados utilizado na prática

## Atividade

Regressão aleatória e normas de reação em dados G×A desbalanceados

- Organizar os dados em formato longo, com as colunas `env`, `gen` e `y`
- Calcular as médias ambientais, a média geral e o índice ambiental centralizado ou padronizado
- Construir manualmente a matriz `Z_u` para intercepto e inclinação aleatórios por genótipo
- Ajustar no `lme4` o modelo `y ~ x + (1 + x | gen)`
- Extrair os BLUPs de intercepto e inclinação dos genótipos
- Interpretar a inclinação total de cada genótipo como efeito fixo médio mais desvio aleatório
- Extrair e interpretar a matriz `G_u`
- Calcular a variância genética ao longo do gradiente ambiental
- Construir gráficos das normas de reação
- Discutir o efeito do desbalanceamento em genótipos com pouca informação
- Discutir conceitualmente como o ajuste mudaria caso os genótipos fossem aparentados
