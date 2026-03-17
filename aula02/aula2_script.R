# =============================================================================
# GMP0164 - Biometria da Interacao Genotipos x Ambientes
# Aula 2: ANOVA Conjunta e Componentes de Variancia em G×A
# Dados: soy_MET.txt - 40 genotipos x 6 ambientes x 3 blocos (DBC)
#
# Modelo assumido (Cruz, Regazzi e Carneiro, 2012 - caso c):
#   Ambientes (A)                          : FIXOS
#   Genotipos (G), G×A, Blocos e Erro      : ALEATORIOS
#
# E(QM) sob esse modelo:
#   E(QM_Erro)   = s2
#   E(QM_GxA)    = s2 + K*l*s2_ga     onde l = J/(J-1)
#   E(QM_G)      = s2 + J*K*s2_g      [s2_ga NAO aparece aqui]
#   E(QM_B(A))   = s2 + I*s2_ba
#   E(QM_A)      = s2 + K*l*s2_ga + I*s2_ba + I*K*phi_a
#
# Testes F corretos:
#   F(G)   = QM_G   / QM_E      [denominador = QM_E, nao QM_GxA]
#   F(GxA) = QM_GxA / QM_E
#   F(A)   = QM_A   / QM_B(A)
# =============================================================================

dat <- read.table("soy_MET.txt", header = TRUE, sep = "\t") #mude o seu diretorio aqui.
dat$gen   <- factor(dat$gen)
dat$env   <- factor(dat$env)
dat$block <- factor(dat$block)

str(dat)
summary(dat)

# =============================================================================
# 1. TABELA DE MEDIAS E GRAFICO DE INTERACAO
# =============================================================================

# Matriz de medias genotipo x ambiente (media sobre os K blocos)
means_ge <- tapply(dat$y, list(dat$gen, dat$env), mean)

# Medias marginais por ambiente (ordenadas crescentemente)
cat("Medias por ambiente (kg/ha):\n")
print(round(sort(colMeans(means_ge)), 1))

# Top 10 genotipos pela media geral sobre todos os ambientes
cat("\nTop 10 genotipos pela media geral (kg/ha):\n")
geno_means <- sort(rowMeans(means_ge), decreasing = TRUE)
print(round(head(geno_means, 10), 1))

# Grafico de interacao: cada linha = um genotipo
# Ambientes ordenados pela media para facilitar leitura visual
# Linhas paralelas  -> interacao simples  (magnitude muda, ranking preservado)
# Cruzamento        -> interacao complexa (mudanca de ranking entre ambientes)
env_order <- order(colMeans(means_ge))

matplot(
  x    = seq_len(ncol(means_ge)),
  y    = t(means_ge[, env_order]),
  type = "b", pch = 16, lty = 1, cex = 0.7,
  xaxt = "n",
  xlab = "Ambiente (ordenado pela media)",
  ylab = "Produtividade media (kg/ha)",
  main = "Grafico de interacao G×A - Soja MET"
)
axis(1, at = seq_len(ncol(means_ge)),
     labels = colnames(means_ge)[env_order])

# =============================================================================
# 2. ANOVA CONJUNTA (DBC, ambientes fixos)
# =============================================================================

# O aov() calcula SQ corretamente para dados balanceados.
# Porem, seus testes F internos assumem modelo fixo.
# Os testes F corretos para o modelo misto (A fixo) sao calculados manualmente.

fit <- aov(y ~ env + env:block + gen + gen:env, data = dat)
tab <- anova(fit)

# Dimensoes do experimento
I <- nlevels(dat$gen)    # numero de genotipos
J <- nlevels(dat$env)    # numero de ambientes
K <- nlevels(dat$block)  # numero de blocos por ambiente
l <- J / (J - 1)         # coeficiente lambda (Cruz et al., 2012)

cat("\nGraus de liberdade:\n")
cat("  Ambiente       :", J - 1, "\n")
cat("  Bloco(Ambiente):", J * (K - 1), "\n")
cat("  Genotipo       :", I - 1, "\n")
cat("  G×A            :", (I - 1) * (J - 1), "\n")
cat("  Erro           :", J * (I - 1) * (K - 1), "\n")
cat("  Total          :", I * J * K - 1, "\n")

# Extrair quadrados medios
MS     <- setNames(tab[, "Mean Sq"], rownames(tab))
QM_E   <- as.numeric(MS["Residuals"])
QM_GxA <- as.numeric(MS["env:gen"])
QM_G   <- as.numeric(MS["gen"])
QM_BA  <- as.numeric(MS["env:block"])
QM_A   <- as.numeric(MS["env"])

# Graus de liberdade para os testes F
df_G   <- I - 1
df_GxA <- (I - 1) * (J - 1)
df_E   <- J * (I - 1) * (K - 1)
df_BA  <- J * (K - 1)
df_A   <- J - 1

# Estatisticas F e p-valores
F_G   <- QM_G   / QM_E
F_GxA <- QM_GxA / QM_E
F_A   <- QM_A   / QM_BA

p_G   <- pf(F_G,   df_G,   df_E,  lower.tail = FALSE)
p_GxA <- pf(F_GxA, df_GxA, df_E,  lower.tail = FALSE)
p_A   <- pf(F_A,   df_A,   df_BA, lower.tail = FALSE)

# Funcao auxiliar para formatar p-valores
fmt_p <- function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p))

cat("\n--- ANOVA conjunta (modelo misto: A fixo) ---\n")
cat(sprintf("  %-15s | %4s | %12s | %8s | %s\n",
            "Fonte", "GL", "QM", "F", "p"))
cat(sprintf("  %-15s | %4d | %12.1f | %8.3f | %s\n",
            "Ambiente",  df_A,   QM_A,   F_A,   fmt_p(p_A)))
cat(sprintf("  %-15s | %4d | %12.1f | %8s | %s\n",
            "Bloco(Amb)", df_BA, QM_BA,  "---", "---"))
cat(sprintf("  %-15s | %4d | %12.1f | %8.3f | %s\n",
            "Genotipo",  df_G,   QM_G,   F_G,   fmt_p(p_G)))
cat(sprintf("  %-15s | %4d | %12.1f | %8.3f | %s\n",
            "G×A",       df_GxA, QM_GxA, F_GxA, fmt_p(p_GxA)))
cat(sprintf("  %-15s | %4d | %12.1f | %8s | %s\n",
            "Erro",      df_E,   QM_E,   "---", "---"))

# =============================================================================
# 3. COMPONENTES DE VARIANCIA (modelo misto, A fixo)
#    Metodo dos momentos - Cruz, Regazzi e Carneiro (2012), caso c
# =============================================================================

# Estimadores:
#   s2    = QM_E
#   s2_ga = (QM_GxA - QM_E) * (J-1) / (K*J)    [denominador = K*l = K*J/(J-1)]
#   s2_g  = (QM_G   - QM_E) / (J*K)             [sem s2_ga no numerador]
#   s2_ba = (QM_BA  - QM_E) / I

s2    <- QM_E
s2_ga <- (QM_GxA - QM_E) * (J - 1) / (K * J)
s2_g  <- (QM_G   - QM_E) / (J * K)
s2_ba <- (QM_BA  - QM_E) / I

vc <- c(s2_g = s2_g, s2_ga = s2_ga, s2_ba = s2_ba, s2 = s2)

cat("\nComponentes de variancia (metodo dos momentos, A fixo):\n")
print(round(vc, 1))

# Proporcoes relativas
# Valores negativos (possivel em amostras finitas) sao truncados em zero
# apenas para o calculo das proporcoes; os valores brutos sao mantidos acima.
vc_pos <- pmax(vc, 0)
sP2    <- sum(vc_pos)

cat("\nProporcao relativa (%):\n")
print(round(100 * vc_pos / sP2, 1))

# Comparacao didatica: efeito do modelo sobre s2_ga
# Sob ambientes fixos, denominador = K*l > K => s2_ga estimada e MENOR
# pois parte da variacao de G×A e absorvida pelos efeitos fixos de ambiente
cat(sprintf(
  "\nComparacao s2_ga: A fixo = %.1f | A aleatorio = %.1f\n",
  (QM_GxA - QM_E) * (J - 1) / (K * J),   # modelo misto (correto aqui)
  (QM_GxA - QM_E) / K                      # modelo totalmente aleatorio
))

# =============================================================================
# 4. CORRELACAO GENETICA HOMOGENEA ENTRE AMBIENTES
# =============================================================================

# Sob o modelo aleatorio independente:
#   u_ij = g_i + (ga)_ij  (valor genetico especifico do ambiente j)
#   Var(u_ij)          = s2_g + s2_ga
#   Cov(u_ij, u_ij')   = s2_g        (j != j')
#
# => rg = s2_g / (s2_g + s2_ga)
#
# rg -> 1: ranking estavel entre ambientes -> adaptacao ampla viavel
# rg -> 0: ranking instavel               -> recomendacao regionalizada necessaria

rg <- vc_pos["s2_g"] / (vc_pos["s2_g"] + vc_pos["s2_ga"])

cat("\nCorrelacao genetica homogenea entre ambientes (rg):", round(rg, 3), "\n")
cat("Razao s2_ga / s2_g (indice de complexidade G×A)  :",
    round(vc_pos["s2_ga"] / vc_pos["s2_g"], 3), "\n")

# =============================================================================
# 5. CORRELACOES PAR A PAR ENTRE AMBIENTES
# =============================================================================

# r_jj' = correlacao de Pearson entre medias dos genotipos em j e j'
# Essa correlacao e usada diretamente na decomposicao S + C (Secao 6)
# e permite identificar pares de ambientes concordantes ou discordantes
# quanto ao ranking de genotipos.

r_pairwise <- cor(means_ge, method = "pearson")

cat("\nMatriz de correlacoes de Pearson entre ambientes:\n")
print(round(r_pairwise, 3))

# =============================================================================
# 6. DECOMPOSICAO DA G×A EM PARTES SIMPLES (S) E COMPLEXA (C)
#    Robertson (1959); Cruz e Castoldi (1991)
#
# Para cada par de ambientes (j, j'):
#   QMGxA_jj' = S + C
#   S  = (1/2) * (sqrt(Q1) - sqrt(Q2))^2
#   C  = sqrt[(1 - r)^3 * Q1 * Q2]       <- Cruz e Castoldi (1991)
#
# onde:
#   Q_j = QM entre genotipos na ANOVA individual do ambiente j
#   r   = correlacao de Pearson entre medias dos genotipos em j e j'
#
# S reflete diferenca de magnitude (ranking preservado)
# C reflete mudanca de ranking (%C > 100% indica r < 0)
# =============================================================================

# 6a. QM entre genotipos por ambiente (ANOVA individual)
Q      <- numeric(J)
names(Q) <- colnames(means_ge)

for (j in colnames(means_ge)) {
  dj    <- dat[dat$env == j, ]
  fit_j <- aov(y ~ block + gen, data = dj)
  Q[j]  <- anova(fit_j)["gen", "Mean Sq"]
}

cat("\nQM entre genotipos por ambiente (Q_j):\n")
print(round(Q, 1))

# 6b. Matrizes S, C e %C para todos os pares de ambientes
envs      <- colnames(means_ge)
mat_S     <- matrix(NA, J, J, dimnames = list(envs, envs))
mat_C     <- matrix(NA, J, J, dimnames = list(envs, envs))
mat_pC    <- matrix(NA, J, J, dimnames = list(envs, envs))
mat_QMGxA <- matrix(NA, J, J, dimnames = list(envs, envs))

for (j1 in 1:(J - 1)) {
  for (j2 in (j1 + 1):J) {
    e1 <- envs[j1];  e2 <- envs[j2]
    Q1 <- Q[e1];     Q2 <- Q[e2]
    r  <- r_pairwise[e1, e2]

    S         <- 0.5 * (sqrt(Q1) - sqrt(Q2))^2
    C         <- sqrt((1 - r)^3 * Q1 * Q2)
    QMGxA_par <- S + C

    mat_S[e1, e2]     <- S
    mat_C[e1, e2]     <- C
    mat_pC[e1, e2]    <- 100 * C / QMGxA_par
    mat_QMGxA[e1, e2] <- QMGxA_par
  }
}

cat("\nQM da G×A por par de ambientes (S + C):\n")
print(round(mat_QMGxA, 1))

cat("\nParte simples (S):\n")
print(round(mat_S, 1))

cat("\nParte complexa (C) - Cruz e Castoldi (1991):\n")
print(round(mat_C, 1))

cat("\n%C (proporcao complexa):\n")
print(round(mat_pC, 1))
cat("Nota: %C > 100% indica correlacao negativa entre os dois ambientes.\n")

# =============================================================================
# 7. ESTRATIFICACAO DE AMBIENTES - TESTES F PAR A PAR
#
# Para cada par (j, j'), testa se a G×A entre esses dois ambientes
# e significativa (Cruz et al., 2012):
#
#   theta_jj' = (1/2)*[sum_i(yb_ij - yb_ij')^2 - (T_j - T_j')^2 / I]
#   F_jj'     = [theta_jj' / (I-1)] / (QM_E / K)
#   GL        = (I-1) e J*(I-1)*(K-1)
#
# Pares nao significativos -> G×A nao detectada -> candidatos ao mesmo
# grupo de recomendacao (estratificacao ambiental; detalhado na Aula 12).
# =============================================================================

mat_F <- matrix(NA, J, J, dimnames = list(envs, envs))
mat_p <- matrix(NA, J, J, dimnames = list(envs, envs))

for (j1 in 1:(J - 1)) {
  for (j2 in (j1 + 1):J) {
    e1 <- envs[j1];  e2 <- envs[j2]
    y1 <- means_ge[, e1]
    y2 <- means_ge[, e2]
    T1 <- sum(y1);   T2 <- sum(y2)

    theta  <- 0.5 * (sum((y1 - y2)^2) - (T1 - T2)^2 / I)
    QM_par <- theta / (I - 1)
    F_par  <- QM_par / (QM_E / K)
    p_par  <- pf(F_par, I - 1, df_E, lower.tail = FALSE)

    mat_F[e1, e2] <- F_par
    mat_p[e1, e2] <- p_par
  }
}

cat("\nTeste F par a par entre ambientes:\n")
print(round(mat_F, 3))

cat("\nValores-p (teste F par a par):\n")
# Formatar p-valores: exibir "< 0.0001" quando muito pequeno
mat_p_fmt <- matrix(
  ifelse(is.na(mat_p), "---",
         ifelse(mat_p < 0.0001, "< 0.0001", sprintf("%.4f", mat_p))),
  nrow = J, dimnames = list(envs, envs)
)
print(mat_p_fmt, quote = FALSE)

cat("\nPares com G×A nao significativa (p > 0.05):\n")
nenhum <- TRUE
for (j1 in 1:(J - 1)) {
  for (j2 in (j1 + 1):J) {
    e1 <- envs[j1];  e2 <- envs[j2]
    if (!is.na(mat_p[e1, e2]) && mat_p[e1, e2] > 0.05) {
      cat(sprintf("  %s x %s : F = %.3f, p = %s\n",
                  e1, e2, mat_F[e1, e2], fmt_p(mat_p[e1, e2])))
      nenhum <- FALSE
    }
  }
}
if (nenhum) cat("  Nenhum par nao significativo encontrado (p > 0.05).\n")

# =============================================================================
# 8. HERDABILIDADE E IMPLICACOES NOS GANHOS DE SELECAO
#
# Herdabilidade para a media na rede MET (J ambientes, K blocos):
#   h2_MET = s2_g / (s2_g + s2_ga/J + s2/(J*K))
#
# Herdabilidade para a media em ambiente unico j (K blocos):
#   h2_j = (s2_g + s2_ga) / (s2_g + s2_ga + s2/K)
#
# ATENCAO: h2_MET e h2_j medem coisas distintas.
#   h2_MET mede a precisao da media sobre a rede -> base para DG_MET
#   h2_j   mede a precisao no ambiente j         -> base para DG_j e RC
#
# Ganho de selecao direto na rede (ganho genetico geral):
#   DG_MET = i * h_MET * sqrt(s2_g)         [unidade: kg/ha]
#
# Ganho de selecao em ambiente unico j (ganho no proprio ambiente):
#   DG_j   = i * h_j * sqrt(s2_g + s2_ga)   [unidade: kg/ha]
#
# DG_MET e DG_j NAO sao diretamente comparaveis:
#   DG_MET = ganho esperado em s2_g (parte estavel entre ambientes)
#   DG_j   = ganho esperado em s2_g + s2_ga (inclui parte especifica do ambiente)
#   DG_j > DG_MET nao significa que selecionar em 1 ambiente e melhor;
#   significa que o ganho e maior naquele ambiente especifico, mas pode
#   nao se repetir em outros ambientes.
#
# Resposta correlacionada (selecao em j, ganho esperado em j'):
#   RC(j'|j) = i * rg * h_j * sqrt(s2_g + s2_ga)
#
# Eficiencia relativa da selecao indireta vs. direta em j'
# (assumindo ambientes com mesma precisao, h_j = h_j'):
#   ER(j->j') = rg
#   ER < 1 sempre que rg < 1: selecionar em 1 ambiente e menos eficiente
#   para ganho em outro ambiente do que selecionar diretamente nele.
# =============================================================================

s2_g_pos  <- vc_pos["s2_g"]
s2_ga_pos <- vc_pos["s2_ga"]
s2_pos    <- vc_pos["s2"]

# Herdabilidades
h2_MET <- s2_g_pos / (s2_g_pos + s2_ga_pos / J + s2_pos / (J * K))
h2_j   <- (s2_g_pos + s2_ga_pos) / (s2_g_pos + s2_ga_pos + s2_pos / K)
h_MET  <- sqrt(h2_MET)
h_j    <- sqrt(h2_j)

cat("\n--- Herdabilidade ---\n")
cat(sprintf("  h2 na rede MET (J=%d ambientes, K=%d blocos) : %.3f\n", J, K, h2_MET))
cat(sprintf("  h2 em ambiente unico (K=%d blocos)           : %.3f\n", K, h2_j))
cat("  Nota: h2_MET e h2_j medem precisoes distintas (ver comentario acima).\n")

# Ganhos de selecao (i = 1 unidade padrao, para fins didaticos)
i_sel  <- 1
DG_MET <- i_sel * h_MET * sqrt(s2_g_pos)
DG_j   <- i_sel * h_j   * sqrt(s2_g_pos + s2_ga_pos)

cat("\n--- Ganho de selecao (i = 1 unidade padrao) ---\n")
cat(sprintf("  DG_MET: ganho genetico geral (rede)     = %.1f kg/ha\n", DG_MET))
cat(sprintf("  DG_j  : ganho no ambiente j             = %.1f kg/ha\n", DG_j))
cat("  DG_j > DG_MET reflete a inclusao de s2_ga em DG_j,\n")
cat("  nao que selecionar em 1 ambiente seja mais eficiente para a rede.\n")

# Resposta correlacionada
RC <- i_sel * rg * h_j * sqrt(s2_g_pos + s2_ga_pos)

cat("\n--- Resposta correlacionada (selecao em j, ganho em j') ---\n")
cat(sprintf("  RC(j'|j) = rg * h_j * sqrt(s2_g + s2_ga) = %.1f kg/ha\n", RC))

# Eficiencia relativa (h_j = h_j' assumido)
ER <- as.numeric(rg)
cat(sprintf("  ER(j->j') = rg = %.3f\n", ER))
cat(sprintf("  Interpretacao: para cada 1 kg/ha de ganho direto em j',\n"))
cat(sprintf("  a selecao indireta via j rende apenas %.3f kg/ha.\n", ER))
cat("  Quanto menor rg, menos eficiente e a selecao em ambiente unico\n")
cat("  para ganhos consistentes em toda a rede.\n")

# =============================================================================
# 9. RESUMO FINAL
# =============================================================================

cat("\n=============================================================\n")
cat("RESUMO FINAL\n")
cat("=============================================================\n")
cat(sprintf("Genotipos (I): %d | Ambientes (J): %d | Blocos (K): %d\n", I, J, K))
cat(sprintf("Teste F para G×A : F = %.2f, p = %s\n", F_GxA, fmt_p(p_GxA)))
cat(sprintf("Teste F para G   : F = %.2f, p = %s\n", F_G,   fmt_p(p_G)))
cat(sprintf("Componentes      : s2_g = %.1f | s2_ga = %.1f | s2 = %.1f\n",
            vc["s2_g"], vc["s2_ga"], vc["s2"]))
cat(sprintf("Correlacao genetica homogenea : rg     = %.3f\n", rg))
cat(sprintf("Herdabilidade na rede         : h2_MET = %.3f\n", h2_MET))
cat(sprintf("Herdabilidade por ambiente    : h2_j   = %.3f\n", h2_j))
cat(sprintf("Ganho genetico geral (i=1)    : DG_MET = %.1f kg/ha\n", DG_MET))
cat(sprintf("Resposta correlacionada (i=1) : RC     = %.1f kg/ha\n", RC))
cat("=============================================================\n")





# =============================================================================
# NOTA 1 - h2_MET vs h2_j
# Neste conjunto de dados, h2_MET (0.894) e h2_j (0.895) sao numericamente
# quase identicos. Isso ocorre por coincidencia dos valores estimados:
# s2_ga/J e s2/(J*K) sao aproximadamente iguais neste dataset, fazendo os
# denominadores das duas expressoes convergirem. Isso NAO significa que as
# formulas sao equivalentes. Em geral, h2_MET < h2_j, pois a media na rede
# diluiu s2_ga pelo fator J. Nao generalize esse resultado para outros datasets.
# =============================================================================

# =============================================================================
# NOTA 2 - %C elevado em todos os pares (87-100%)
# A decomposicao S + C revelou que a G×A nessa rede e predominantemente
# complexa em praticamente todos os pares de ambientes. Isso e coerente com
# rg = 0.654 (moderado) e com todos os testes F par a par sendo significativos.
# Do ponto de vista do melhoramento, isso indica que o ranking dos genotipos
# muda substancialmente entre os ambientes do Cerrado avaliados, tornando
# a recomendacao unica (adaptacao ampla) uma estrategia de risco. A
# regionalizacao da recomendacao e, portanto, justificada por esses dados.
# =============================================================================
