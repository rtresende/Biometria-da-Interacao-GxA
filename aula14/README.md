# Aula 14 — Ambientômica I: Modelagem Integrativa

**Data:** 17 de junho de 2026

## Conteúdo abordado

- Diferenças entre SIG, ambientipagem, fenômica e Ambientômica
- Analogia entre genômica e Ambientômica: genótipo, genoma, ambiótipo, ambientoma e marcador ambientômico
- Premissas operacionais dos marcadores ambientômicos
- Estratégias para construção de marcadores: covariáveis brutas, transformações, índices ecofisiológicos e marcadores engenheirados por IA
- Normas de reação com gradientes ambientais construídos a partir de covariáveis
- Kernel ambientômico linear e não linear
- Regressão aleatória sobre marcadores ambientômicos
- Comparação entre kernel ambientômico e regressão aleatória
- Portabilidade entre GBLUP/RRBLUP e limitações da analogia em normas de reação
- Abordagem latente com FA, PLS e predição de ambientes não testados
- Validação preditiva em CV0, CV1, CV2 e diferentes definições de ambiente novo

## Material

- 📄 `Biometria_da_Interação_GxA_GMP0164_aula14.pdf` — roteiro teórico da aula e atividade prática
- 📊 Dados didáticos embutidos no roteiro — 14 observações, 4 genótipos, 6 locais e covariáveis `EC1` a `EC5`

## Atividade

Norma de reação ambientômica por regressão aleatória.

- Montar a matriz `Z` bloco-diagonal, com intercepto e covariáveis ambientais por genótipo
- Montar o kernel `K <- kronecker(A, diag(6))`, usando a matriz de parentesco dos quatro genótipos
- Ajustar o modelo misto `y = Xb + Zu + e` com `mixed.solve`
- Extrair os efeitos aleatórios `u_hat` por genótipo e covariável
- Predizer `y_hat` para um local não testado, usando valores plausíveis de `EC1` a `EC5`
- Comparar os genótipos no novo local
- Interpretar qual covariável mais desloca a resposta
- Discutir se a similaridade entre o local novo e os locais testados sustenta a predição
