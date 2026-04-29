# Aula 8 — WAASB e WAASBY (Índices de Adaptabilidade e Estabilidade via BLUP)

**Data:** [29 de abril de 2026]  
**Carga horária:** [4h]

**Responsável pela aula:** Prof. Dr. Tiago Olivoto (UFSC)

## Conteúdo abordado

- Motivação: quantificar e interpretar a interação G×A sob modelos mistos  
- Decomposição da interação G×E a partir de BLUPs (matriz GE)  
- SVD aplicada à matriz de efeitos da interação  
- Índice WAASB: definição, cálculo e interpretação  
- Relação WAASB × estabilidade  
- Biplot produtividade × estabilidade (WAASB × Y)  
- Índice WAASBY: integração entre produtividade e estabilidade  
- Reescalonamento (0–100) e ponderação (θ)  
- Interpretação por quadrantes e ranking genotípico  
- Comparação com AMMI e índices clássicos  

## Material

- 📄 `Biometria_da_Interação_GxA_GMP0164_Slides_aula8_Olivoto.pdf` — slides da aula (Prof. Tiago Olivoto)  
- 💾 `script_GMP0164_aula8.R` — script utilizado na aula  
- 📊 `soy_MET.txt` — conjunto de dados multiambiente utilizado  

## Atividade

Cálculo e interpretação de WAASB e WAASBY em dados MET

- Ajustar modelo misto para obter BLUPs genótipo × ambiente  
- Extrair a matriz de efeitos da interação (GE)  
- Aplicar SVD e obter escores principais  
- Calcular WAASB para cada genótipo  
- Construir gráfico WAASB × produtividade  
- Calcular WAASBY com diferentes pesos (θ)  
- Comparar ranking com média fenotípica e estabilidade clássa  
