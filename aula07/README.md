# Aula 07 — Modelos Mistos II: Estruturas de Variância-Covariância

Data: 21 de abril de 2026  
Carga horária: 4h

## Conteúdo abordado

- Reformulação da G×A em termos da matriz genética entre ambientes, G_e
- Produto de Kronecker na notação do modelo misto multiambiente
- Estruturas VCOV para G_e: ID, DIAG, CS, CSH, UN, AR(1) e FA
- Interpretação biométrica de variâncias, covariâncias e correlações genéticas entre ambientes
- Equivalência entre o modelo clássico g + ga e a estrutura CS
- Critérios de escolha de estrutura: plausibilidade biométrica, estabilidade numérica e ajuste
- Efeito da estrutura VCOV sobre os BLUPs e sobre o compartilhamento de informação entre ensaios
- Comparação prática de modelos via script em R

## Material

- 📄 `Biometria_da_Interação_GxA_GMP0164_aula7.pdf` — roteiro teórico da aula
- 💾 `script_GMP0164_aula7.R` — script utilizado na aula
- 📊 `maize_safrinha_MET.txt` — conjunto de dados utilizado na prática

## Atividade

Derivação da estrutura genética implícita no modelo g + ga

- Derivar algebricamente a variância de u_ij = g_i + (ga)_ij
- Derivar a covariância entre u_ij e u_ij' para j ≠ j'
- Escrever explicitamente a matriz G_e para J = 3 ambientes
- Identificar a estrutura VCOV correspondente
- Assumindo σ^2_g = 3,0 e σ^2_ga = 1,5, calcular a correlação genética uniforme implícita
- Interpretar o significado desse valor para a consistência do ranking dos genótipos entre ambientes
- Discutir em quais condições a hipótese de correlação uniforme pode ser razoável e em quais pode ser inadequada
