# Aula 06 — Modelos Mistos I: Introdução ao BLUP Multiambiente

Material da Aula 06 da disciplina GMP0164 — Biometria da Interação Genótipos × Ambientes, ofertada no Programa de Pós-Graduação em Genética e Melhoramento de Plantas da UFG.

📄 **Roteiro completo (PDF):** :contentReference[oaicite:0]{index=0}

Esta aula introduz o modelo linear misto no contexto de ensaios multiambiente (MET), estabelecendo a base para a formulação moderna da interação G×A. O foco está na passagem da análise baseada em efeitos fixos para uma abordagem em que parte da variação é modelada por componentes de variância e os valores genotípicos são obtidos por predição (BLUP).

---

## Objetivos da aula

Ao final desta aula, o estudante deve ser capaz de:

- compreender as limitações da abordagem baseada em efeitos fixos em dados desbalanceados;
- reconhecer a estrutura do modelo linear misto;
- distinguir efeitos fixos e aleatórios no contexto de MET;
- interpretar os componentes de variância associados a genótipos, interação G×A e erro;
- entender a lógica da estimação por REML;
- interpretar o BLUP como predição de valores genotípicos;
- compreender o papel do *shrinkage* e da PEV na qualidade da predição;
- relacionar acurácia seletiva e ganho esperado com base nos BLUPs.

---

## Escopo

A aula apresenta a formulação introdutória do BLUP multiambiente. A interação G×A é incorporada como efeito aleatório, mas sem detalhamento de estruturas de covariância entre ambientes.

---

## Conteúdo

O roteiro aborda:

- motivação para o uso de modelos mistos em MET;
- formulação matricial do modelo linear misto;
- interpretação dos termos \( \mathbf{X}\boldsymbol{\beta} \), \( \mathbf{Z}\mathbf{u} \) e \( \mathbf{e} \);
- distinção entre efeitos fixos e aleatórios;
- componentes de variância e sua leitura biométrica;
- estimação por REML;
- equações de modelos mistos de Henderson;
- definição de BLUE e BLUP;
- variância do erro de predição (PEV);
- acurácia seletiva;
- interpretação do *shrinkage*;
- modelagem da interação G×A no contexto misto;
- extensão conceitual para incorporação de parentesco (\( \mathbf{A} \)) e relação com GBLUP.

---

## Organização dos arquivos

Esta pasta pode conter:

- `aula6.pdf` — roteiro completo;
- arquivos `.tex`;
- scripts em R utilizados na atividade;
- dados auxiliares, quando aplicável.

---

## Software

As análises são conduzidas em R. Pacotes utilizados:

- `lme4`
- `metan`
- `dplyr`

---

## Atividade

A atividade utiliza o conjunto `data_ge` (pacote `metan`), com desbalanceamento artificial da tabela G×A.

O estudante deve:

1. ajustar o modelo misto via `lme4`;
2. estimar componentes de variância por REML;
3. calcular herdabilidade e acurácia seletiva;
4. extrair BLUPs e PEVs;
5. comparar o ranking por BLUP com o ranking por média.

---

## Leituras de apoio

- Henderson, C. R. (1975). Biometrics.
- Patterson, H. D., & Thompson, R. (1971). Biometrika.
- Searle, S. R. et al. (1992). Variance Components.
- Lynch, M., & Walsh, B. (1998). Genetics and Analysis of Quantitative Traits.
- Piepho, H.-P. et al. (2008). Euphytica.
- Resende, M. D. V. (2007). Embrapa.

---

## Observação

Este material estabelece a base conceitual para o uso de BLUP em dados multiambiente. A partir desta aula, a disciplina passa a tratar explicitamente a G×A como problema de modelagem de variâncias e predição, o que fundamenta as abordagens desenvolvidas nas aulas seguintes.
