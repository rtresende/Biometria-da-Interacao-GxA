# =============================================================================
# GMP0164 - Biometria da Interacao Genotipos x Ambientes
# Aula 2: ANOVA Conjunta e Componentes de Variancia em G×A
# Dados: soy_MET.txt - 40 genotipos x 6 ambientes x 3 blocos (DBC)
# Prof. Dr. Rafael Tassinari Resende
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

dat <- read.table("G:/Meu Drive/UFG/PPGGMP/BiometriaDaInteracaoGxA/aulas/simu_data/soy_MET.txt", header = TRUE, sep = "\t")
dat$gen   <- factor(dat$gen)
dat$env   <- factor(dat$env)
dat$block <- factor(dat$block)

str(dat)
summary(dat)

# =============================================================================
# 1. TABELA DE MEDIAS E GRAFICO DE INTERACAO
# =============================================================================

means_ge  <- tapply(dat$y, list(dat$gen, dat$env), mean)
env_order <- order(colMeans(means_ge))
mat_plot  <- t(means_ge[, env_order])   # ambientes nas linhas, genotipos nas colunas
nx        <- nrow(mat_plot)             # numero de ambientes

cat("Medias por ambiente (kg/ha):\n")
print(round(sort(colMeans(means_ge)), 1))
cat("\nTop 10 genotipos pela media geral (kg/ha):\n")
print(round(head(sort(rowMeans(means_ge), decreasing = TRUE), 10), 1))

# Identificacao por prefixo
all_gens  <- colnames(mat_plot)
is_check  <- grepl("^CHK", all_gens)
chk_names <- all_gens[is_check]
lin_names <- all_gens[!is_check]

# Paleta colorblind-safe (Bang Wong, 2011) + formas distintas por check
chk_colors <- c(CHK01 = "#0072B2", CHK02 = "#D55E00",
                CHK03 = "#9400D3", CHK04 = "#000000")
chk_pch    <- c(CHK01 = 16, CHK02 = 17, CHK03 = 15, CHK04 = 18)
lin_color  <- "gray65"

# Grafico: linhagens no fundo, checks em destaque
plot(NA, xlim = c(1, nx + 0.6), ylim = range(means_ge) * c(0.97, 1.03),
     xaxt = "n", xlab = "Ambiente (ordenado pela media)",
     ylab = "Produtividade media (kg/ha)",
     main = "Grafico de interacao G×A — Soja MET")
axis(1, at = 1:nx, labels = rownames(mat_plot))

# Linhagens (cinza, fundo) — matlines/matpoints evitam loop
matlines( 1:nx, mat_plot[, lin_names], col = lin_color, lwd = 0.8, lty = 1)
matpoints(1:nx, mat_plot[, lin_names], col = lin_color, pch = 16,  cex = 0.5)

# Checks (cores + formas distintas, destaque) — rotulo na ponta direita
for (g in chk_names) {
  lines( 1:nx, mat_plot[, g], col = chk_colors[g], lwd = 2.8)
  points(1:nx, mat_plot[, g], col = chk_colors[g], pch = chk_pch[g], cex = 1.3)
  text(nx + 0.12, mat_plot[nx, g], labels = g,
       col = chk_colors[g], cex = 0.82, adj = 0, font = 2)
}

legend("topleft",
       legend = c(chk_names, sprintf("Linhagens (n=%d)", length(lin_names))),
       col    = c(chk_colors[chk_names], lin_color),
       lwd    = c(rep(2.8, length(chk_names)), 0.8),
       pch    = c(chk_pch[chk_names], 16),
       pt.cex = c(rep(1.3, length(chk_names)), 0.7),
       bty    = "n", cex = 0.82)

# Interpretacao do grafico de interacao com presenca de checks:
#
# (1) ESCALA DE REFERENCIA: a faixa CHK03--CHK04 delimita o intervalo de
#     cultivares comerciais estabelecidas. Linhagens que superam todos os
#     checks de forma consistente sao candidatas prioritarias ao avanco.
#
# (2) ATENCAO — os proprios checks se reordenam entre ambientes (CHK04
#     abaixo de CHK01 em E1-E2, liderando de E3 em diante): a G×A nao e
#     exclusividade das linhagens. %C elevado na decomposicao S+C reflete
#     discordancia estrutural da rede, nao apenas instabilidade das linhagens.
#
# (3) VIES NOS COMPONENTES DE VARIANCIA: a amplitude CHK03--CHK04
#     (~1.000-2.500 kg/ha) e comparavel a das linhagens e infla sigma2_g,
#     rg e h2_MET — pois checks sao genotipos fixos, nao amostras da
#     populacao em melhoramento. Os valores impressos a seguir descrevem
#     o conjunto heterogeneo (checks + linhagens), nao a populacao sob
#     selecao. Veja a Secao 3 do script para a analise restrita as linhagens.

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

# ATENCAO: sigma2_g = 51,6% inclui a amplitude entre checks (CHK03--CHK04).
# Restrito as 36 linhagens, sigma2_g tende a cair e sigma2_ga a subir
# proporcionalmente — a situacao real das linhagens e provavelmente
# menos favoravel do que esses numeros indicam. Ver analise restrita adiante.

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
# 5. CORRELACOES FENOTÍPICAS PAR A PAR ENTRE AMBIENTES
# =============================================================================

# r_jj' = correlacao de Pearson entre medias dos genotipos em j e j'
# Essa correlacao e usada diretamente na decomposicao S + C (Secao 6)
# e permite identificar pares de ambientes concordantes ou discordantes
# quanto ao ranking de genotipos.

r_pairwise <- cor(means_ge, method = "pearson")

cat("\nMatriz de correlacoes de Pearson entre ambientes:\n")
print(round(r_pairwise, 3))

# NOTA: r_jj' eh fenotipica (contem g_i, (ga)_ij e residuo) — tende a
# subestimar a concordancia genetica real. Para a populacao sob selecao,
# use rg (Secao 7); correlacoes geneticas par a par robustas requerem
# estrutura de covariancia nao-estruturada via REML (Aula 6).

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
# h2_MET = s2_g / (s2_g + s2_ga/J + s2/(J*K))
#   Numerador: s2_g apenas — parte genetica estavel entre ambientes.
#   Base para selecao de ampla adaptacao.
#
# h2_j (medio) = (s2_g + s2_ga) / (s2_g + s2_ga + s2/K)
#   Numerador: s2_g + s2_ga — em um unico ambiente as duas partes sao
#   indistinguiveis e ambas constituem variacao herdavel local.
#   Usa componentes medios da ANOVA conjunta — representa um ambiente
#   hipotetico medio da rede.
#
# h2_j (local) = (Q_j - QME_j/K) / Q_j
#   Estimado por ambiente via ANOVA individual. Nao assume homogeneidade
#   de variancias residuais entre ambientes.
#
# DG_MET = i.s. * h_MET * sqrt(s2_g)         [ganho genetico geral na rede]
# DG_j   = i.s. * h_j   * sqrt(s2_g + s2_ga) [ganho no ambiente medio]
#   DG_j > DG_MET porque inclui s2_ga no desvio-padrao genotipico.
#   Nao implica que selecionar em 1 ambiente seja melhor para a rede.
#
# RC(j'|j) = i * rg * h_j * sqrt(s2_g + s2_ga)
# ER(j->j') = rg * (h_j / h_j')
#   Com h_j ≈ h_j' entre ambientes, ER ≈ rg.
#   Selecionar em 1 ambiente desperdiça (1 - rg)*100% do potencial de ganho.
# =============================================================================

s2_g_pos  <- vc_pos["s2_g"]
s2_ga_pos <- vc_pos["s2_ga"]
s2_pos    <- vc_pos["s2"]

# --- Herdabilidade media (componentes da ANOVA conjunta) ---
h2_MET <- s2_g_pos / (s2_g_pos + s2_ga_pos / J + s2_pos / (J * K))
h2_j   <- (s2_g_pos + s2_ga_pos) / (s2_g_pos + s2_ga_pos + s2_pos / K)
h_MET  <- sqrt(h2_MET)
h_j    <- sqrt(h2_j)

cat("\n--- Herdabilidade ---\n")
cat(sprintf("  h2_MET (J=%d, K=%d) : %.3f\n", J, K, h2_MET))
cat(sprintf("  h2_j medio (K=%d)   : %.3f\n", K, h2_j))

# --- Herdabilidade local por ambiente (Q_j ja disponivel da Secao 6) ---
QME_local <- setNames(numeric(J), envs)
for (j in envs)
  QME_local[j] <- anova(aov(y ~ block + gen,
                             data = dat[dat$env == j, ]))["Residuals", "Mean Sq"]

h2_j_local <- (Q - QME_local / K) / Q
h_j_local  <- sqrt(h2_j_local)

cat("\n  h2_j por ambiente:\n")
print(round(h2_j_local, 3))
cat(sprintf("  Amplitude: %.3f — h2_j estatisticamente equivalentes entre ambientes.\n",
            diff(range(h2_j_local))))
cat("  O limitante da eficiencia de selecao e rg, nao h2_j.\n")

# --- Ganhos de selecao (i.s. = 1 unidade padrao) ---
DG_MET <- h_MET * sqrt(s2_g_pos)
DG_j   <- h_j   * sqrt(s2_g_pos + s2_ga_pos)

cat("\n--- Ganho de selecao (i.s. = 1) ---\n")
cat(sprintf("  DG_MET : %.1f kg/ha\n", DG_MET))
cat(sprintf("  DG_j   : %.1f kg/ha  (diferenca de %.1f kg/ha reflete s2_ga em DG_j)\n",
            DG_j, DG_j - DG_MET))

# --- Resposta correlacionada e eficiencia relativa ---
RC       <- as.numeric(rg) * h_j * sqrt(s2_g_pos + s2_ga_pos)
mat_ER   <- outer(h_j_local, h_j_local, function(a, b) as.numeric(rg) * a / b)
diag(mat_ER) <- NA

cat("\n--- Resposta correlacionada e eficiencia relativa ---\n")
cat(sprintf("  RC(j'|j) medio  : %.1f kg/ha\n", RC))
cat(sprintf("  rg              : %.3f\n", as.numeric(rg)))
cat("\n  Matriz ER(j->j') = rg * h_j / h_j':\n")
print(round(mat_ER, 3))
cat(sprintf("\n  ER ≈ rg = %.3f em todos os pares (h_j homogeneos).\n",
            as.numeric(rg)))
cat(sprintf("  Selecionar em 1 ambiente desperdiça %.1f%% do potencial de ganho na rede.\n",
            (1 - as.numeric(rg)) * 100))



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
cat(sprintf("Correlacao genetica homogenea : rg        = %.3f\n", as.numeric(rg)))
cat(sprintf("Herdabilidade na rede         : h2_MET    = %.3f\n", h2_MET))
cat(sprintf("Herdabilidade ambiente medio  : h2_j med  = %.3f\n", h2_j))
cat(sprintf("Herdabilidade por ambiente    : h2_j [%s] = %.3f  [%s] = %.3f  [%s] = %.3f\n",
            envs[1], h2_j_local[1], envs[2], h2_j_local[2], envs[3], h2_j_local[3]))
cat(sprintf("                                [%s] = %.3f  [%s] = %.3f  [%s] = %.3f\n",
            envs[4], h2_j_local[4], envs[5], h2_j_local[5], envs[6], h2_j_local[6]))
cat(sprintf("Ganho genetico geral (i.s.=1)    : DG_MET    = %.1f kg/ha\n", DG_MET))
cat(sprintf("Ganho ambiente medio (i.s.=1)    : DG_j      = %.1f kg/ha\n", DG_j))
cat(sprintf("Resposta correlacionada (i.s.=1) : RC        = %.1f kg/ha\n", RC))
cat(sprintf("Eficiencia relativa media     : ER        = %.3f\n", as.numeric(rg)))
cat(sprintf("Desperdicio selecao 1 amb.    : 1-ER      = %.1f%%\n",
            (1 - as.numeric(rg)) * 100))
cat("=============================================================\n")


# TO NOT RUN — analise restrita as linhagens (remove checks)
# Checks sao genotipos fixos (escolhidos deliberadamente) e nao pertencem
# a populacao sob selecao. Sua remocao tende a reduzir s2_g, rg e h2_MET,
# revelando a situacao real das linhagens experimentais.
# dat <- dat[!dat$gen %in% c("CHK01","CHK02","CHK03","CHK04"), ]
# dat$gen <- droplevels(dat$gen)
