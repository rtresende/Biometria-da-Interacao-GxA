# Aula 03 — Regressão conjunta para adaptabilidade e estabilidade

**Data:** 25 de março de 2026  
**Carga horária:** 4h

## Conteúdo abordado

- Transição da ANOVA conjunta para modelagem da resposta individual dos genótipos  
- Modelo de regressão conjunta: \( \bar{y}_{ij} = \mu_i + b_i I_j + \varepsilon_{ij} \)  
- Construção e interpretação do índice ambiental (\(I_j\))  
- Conceitos biométricos de adaptabilidade e estabilidade (estática vs dinâmica)  
- Método de Finlay e Wilkinson: interpretação de \( \hat{b}_i \)  
- Método de Eberhart e Russell: inclusão dos desvios da regressão (\( \hat{\sigma}^2_{di} \))  
- Testes de hipótese para adaptabilidade (\( \beta_{1i} = 1 \)) e estabilidade (\( \sigma^2_{di} = 0 \))  
- Decomposição da soma de quadrados (parte linear e desvios)  
- Relação entre média, adaptabilidade e previsibilidade na recomendação  
- Limitações do índice ambiental e implicações práticas (dependência, linearidade, alavancagem)

## Material

- 📄 `Biometria_da_Interação_GxA_aula3.pdf` — roteiro completo da aula  
- 💾 `aula3_script.R` — script R utilizado durante a aula  

## Atividade

Leitura biométrica da regressão conjunta para genótipos selecionados:

- Construção de tabela-resumo com \( \bar{y}_{i\cdot} \), \( \hat{\beta}_{1i} \), \( \hat{\sigma}^2_{di} \) e \( R^2_i \)  
- Classificação quanto à adaptabilidade (geral, favorável, desfavorável)  
- Avaliação da previsibilidade (desvios da regressão)  
- Indicação de genótipos para recomendação ampla, ambientes restritivos e casos com cautela  
- Comparação conceitual entre Finlay–Wilkinson e Eberhart–Russell  
