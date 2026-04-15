# Aula 06 — Modelos Mistos I: Introdução ao BLUP Multiambiente

**Data:** 15 de abril de 2026  
**Carga horária:** 4h

## Conteúdo abordado

- Limitações da abordagem baseada em efeitos fixos em dados desbalanceados
- Formulação do modelo linear misto (\( y = X\beta + Zu + e \))
- Distinção entre efeitos fixos e aleatórios no contexto de MET
- Componentes de variância (\( \sigma_g^2, \sigma_{ga}^2, \sigma_e^2 \))
- Estimação por REML
- Equações de modelos mistos de Henderson (MME)
- Definição de BLUE e BLUP
- Variância do erro de predição (PEV)
- Acurácia seletiva
- Interpretação do *shrinkage*
- Modelagem da interação G×A como efeito aleatório
- Extensão conceitual para incorporação de parentesco (\( \mathbf{A} \)) e GBLUP

## Material

- 📄 `Biometria_da_Interação_GxA_GMP0164_aula6.pdf` — roteiro completo da aula

## Atividade

Análise de dados multiambiente desbalanceados:

- Ajuste de modelo misto via `lme4`
- Estimação de componentes de variância por REML
- Cálculo de \( h^2 \), \( \hat{H}^2 \) e acurácia seletiva
- Extração de BLUPs e PEVs
- Comparação entre ranking por BLUP e por média aritmética
- Identificação de genótipos mais afetados pelo desbalanceamento
