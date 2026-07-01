# Aula 15 — Ambientômica II: Síntese e Integração Prática

**Data:** 24 de junho de 2026  
**Carga horária:** 1h40

## Conteúdo abordado

- População-alvo de ambientes (TPE) e domínio de recomendação
- Marcadores ambientômicos (`EA`) e matriz ambiental `W`
- Kernel ambientômico e PCA ambiental
- Modelos baseline para comparação
- Regressão aleatória com PC1–PC5
- Validação CV0 para ambientes não testados
- Capacidade preditiva (`pa`) dentro de ambientes ocultos
- Predição espacial para pixels da TPE
- Mapas `which-won-where` e desempenho predito do vencedor
- Normas de reação em gradiente ambiental amplo
- Ganho contra genótipo de referência
- Confiabilidade: interpolação versus extrapolação
- Herdabilidade espacial aproximada
- Zonas de melhoramento por perfis genotípicos preditos
- Diferença entre predição, recomendação e deployment

## Material

- 📄 `Biometria_da_Interação_GxA_GMP0164_aula15.pdf` — roteiro teórico da aula
- 💾 `script_GMP0164_aula15.R` — script utilizado na prática
- 📊 `data_GMP0164_aula15.rds` — dados preparados para a prática

## Atividade

Aplicação de um fluxo ambientômico completo, partindo dos marcadores ambientais até produtos de decisão para melhoramento.

- Carregar dados fenotípicos e ambientais preparados
- Visualizar a TPE e os ambientes testados
- Construir kernel ambientômico com todos os marcadores `EA`
- Sintetizar os marcadores ambientais por PCA
- Ajustar modelo baseline: `Yield ~ 1 + (1 | Genotype)`
- Ajustar regressão aleatória ambientômica com PC1–PC5
- Realizar validação CV0 ocultando ambientes inteiros
- Calcular capacidade preditiva (`pa`) geral e dentro de cada trial
- Predizer todos os genótipos em todos os pixels da TPE
- Gerar mapa `which-won-where`
- Construir normas de reação em gradiente ambiental amplo
- Calcular ganho contra o genótipo de referência `G001`
- Classificar pixels em interpolação ou extrapolação
- Mapear herdabilidade espacial aproximada
- Definir zonas de melhoramento por agrupamento Ward dos perfis preditos
