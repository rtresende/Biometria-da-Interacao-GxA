# Biometria da Interação Genótipos × Ambientes

**Disciplina:** GMP0164 — Tópicos Especiais  
**Programa:** Pós-Graduação em Genética e Melhoramento de Plantas — PPGGMP/UFG  
**Professor:** Prof. Dr. Rafael Tassinari Resende  
**Semestre:** 2026/1 (início: março/2026)  
**Carga horária:** 64h | 4 créditos  

---

## Sobre a disciplina

Apresenta os fundamentos estatísticos da interação genótipos × ambientes (G×A)
e suas aplicações no melhoramento de plantas, integrando métodos clássicos e
modelos mistos multiambiente. Os tópicos cobrem desde a ANOVA conjunta até
modelagem ambientômica para integração G×A e recomendação de genótipos.

---

## Conteúdo programático

| Aula | Data | Tema |
|---|---|---|
| 1 | 11/03 | Introdução conceitual e histórica da G×A |
| 2 | 18/03 | ANOVA conjunta e componentes de variação |
| 3 | 25/03 | Regressão de Finlay-Wilkinson |
| 4 | 01/04 | Índices clássicos de adaptabilidade e estabilidade |
| 5 | 08/04 | AMMI e GGE Biplot |
| 6 | 15/04 | Modelos mistos I — BLUP multiambiente |
| 7 | 22/04 | Modelos mistos II — Estruturas de variância |
| 8 | 29/04 | Modelos mistos III — Índices via BLUP |
| 9 | 06/05 | Regressão aleatória e normas de reação |
| 10 | 13/05 | Predição e imputação em redes multiambiente |
| 11 | 20/05 | SIG para dados de melhoramento |
| 12 | 27/05 | Estratificação ambiental |
| 13 | 03/06 | Ambientipagem e covariáveis ambientais |
| 14 | 10/06 | Ambientômica I |
| 15 | 17/06 | Ambientômica II |
| 16 | 24/06 | Avaliação final |

---

## Estrutura do repositório
```
data/       # Datasets utilizados nas aulas
scripts/    # Scripts de simulação e suporte
lessons/    # Materiais e scripts por aula
docs/       # Plano de ensino e documentos
```

---

## Requisitos de software
```r
install.packages(c(
  "AlphaSimR",  # simulação genômica
  "lme4",       # modelos mistos
  "metan",      # AMMI, GGE, estabilidade
  "ggplot2"     # visualização
# entre vários outros!
))
```

---

## Licença

Material licenciado sob [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).  
Uso acadêmico livre com atribuição. Uso comercial não permitido.

© 2026 Rafael Tassinari Resende — Escola de Agronomia / UFG
