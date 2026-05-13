# Aula 10 — Predição de Desempenho em Múltiplos Ambientes

**Data:** 13 de maio de 2026  
**Carga horária:** 2h

## Conteúdo abordado

- Desbalanceamento em redes MET: cenários de dados faltantes internos e ambientes não observados
- Input para imputação: BLUEs, BLUPs e desregressão (problema do duplo encolhimento)
- Imputação por decomposição matricial: EM-AMMI e SVD iterativo (Yan, 2013)
- Imputação por modelos mistos (UN/FA) como preditores naturais da matriz G×A
- Imputação por aprendizado de máquina: MICE e missForest
- Acurácia preditiva e validação cruzada (CV1 e CV2)
- Estratégias de seleção com a matriz G×A completa: adaptação ampla, específica e índices multi-ambiente

## Material

- 📄 `Biometria_da_Interação_GxA_GMP0164_aula10.pdf` — roteiro teórico da aula
- 💾 `aula10_script.R` — script utilizado na aula prática
- 📊 `soy_MET_phenotypic.txt` — dados fenotípicos MET de soja simulados (40 genótipos × 6 ambientes × 3 blocos, DBC; gerados por `simu_soybean-data_balanced.R`)

## Atividade

Comparação de imputadores da matriz G×A sob validação cruzada, e avaliação do impacto do desbalanceamento no ranking genotípico.

- Construir a matriz de médias G×A a partir dos dados balanceados de soja
- Induzir ausências artificiais (CV2, 20% das células) para simular rede desbalanceada
- Imputar a matriz via SVD iterativo (testando sensibilidade ao número de componentes) e via missForest
- Comparar a acurácia preditiva (ρ_pred) dos dois métodos e discutir as diferenças à luz das suposições de cada abordagem
- Avaliar o impacto da imputação no ranking genotípico (correlação de Spearman entre rankings com e sem imputação)
- [Opcional] Executar esquema CV1 (remoção de ambiente inteiro) e comparar acurácia por local
