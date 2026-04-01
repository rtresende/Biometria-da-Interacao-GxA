# =============================================================================
# GMP0164 - Biometria da Interacao Genotipos x Ambientes
# Aula 4: Indices Classicos de Estabilidade e Criterios de Recomendacao
# Dados: soy_MET.txt -- 40 genotipos, 6 ambientes, DBC (3 blocos/ambiente)
# Prof. Dr. Rafael Tassinari Resende
# =============================================================================

setwd("G:/Meu Drive/UFG/PPGGMP/BiometriaDaInteracaoGxA/aulas")

# -----------------------------------------------------------------------------
# 1. Leitura dos dados
# -----------------------------------------------------------------------------

dados <- read.table("simu_data/soy_MET.txt", header = TRUE, sep = "\t")

# -----------------------------------------------------------------------------
# 2. Pre-processamento basico
# -----------------------------------------------------------------------------

colunas_necessarias <- c("gen", "env", "block", "y")
if (!all(colunas_necessarias %in% names(dados))) {
  stop("O arquivo de entrada deve conter as colunas: gen, env, block, y.")
}

dados$gen   <- factor(dados$gen)
dados$env   <- factor(dados$env)
dados$block <- factor(dados$block)

dados <- dados[, c("gen", "env", "block", "y")]

cat("\nEstrutura dos dados de entrada:\n")
str(dados)

cat("\nResumo da variavel resposta (kg/ha):\n")
print(summary(dados$y))

# -----------------------------------------------------------------------------
# 3. Funcoes auxiliares
# -----------------------------------------------------------------------------

formata_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p)))
}

ajusta_faixa <- function(x, frac = 0.08) {
  r <- range(x, na.rm = TRUE)
  d <- diff(r)
  if (d == 0) d <- max(abs(r), 1) * frac
  r + c(-1, 1) * d * frac
}

indice_lin_binns <- function(mat_y, idx_env = seq_len(ncol(mat_y))) {
  if (length(idx_env) == 0) return(rep(NA_real_, nrow(mat_y)))
  mat_sub <- mat_y[, idx_env, drop = FALSE]
  max_j   <- apply(mat_sub, 2, max)
  d2      <- sweep(mat_sub, 2, max_j, "-")^2
  rowSums(d2) / (2 * ncol(mat_sub))
}

indice_annicchiarico <- function(mat_y, idx_env = seq_len(ncol(mat_y)), alpha = 0.25) {
  if (length(idx_env) == 0) return(rep(NA_real_, nrow(mat_y)))
  mat_sub  <- mat_y[, idx_env, drop = FALSE]
  media_j  <- colMeans(mat_sub)
  z_ij     <- sweep(mat_sub, 2, media_j, "/") * 100
  z_media  <- rowMeans(z_ij)
  z_dp     <- apply(z_ij, 1, sd)
  z_crit   <- qnorm(1 - alpha)
  z_media - z_crit * z_dp
}

# -----------------------------------------------------------------------------
# 4. Dimensoes experimentais e checagens do delineamento
# -----------------------------------------------------------------------------

I <- nlevels(dados$gen)   # numero de genotipos
J <- nlevels(dados$env)   # numero de ambientes

blocos_por_amb <- tapply(dados$block, dados$env, function(x) length(unique(x)))
if (length(unique(blocos_por_amb)) != 1) {
  stop("O numero de blocos nao e constante entre ambientes. Este script assume DBC balanceado.")
}
K <- unique(blocos_por_amb)[1]

n_celula <- xtabs(~ gen + env, data = dados)
if (any(n_celula == 0)) {
  stop("Existem celulas genotipo x ambiente ausentes. Os indices deste script assumem tabela G x A completa.")
}
if (length(unique(as.vector(n_celula))) != 1) {
  stop("O numero de observacoes por celula genotipo x ambiente nao e constante. Este script assume balanceamento.")
}

cat("\nDimensoes experimentais:\n")
cat("  Genotipos           :", I, "\n")
cat("  Ambientes           :", J, "\n")
cat("  Blocos por ambiente :", K, "\n")
cat("  Total de observacoes:", nrow(dados), "\n")

# -----------------------------------------------------------------------------
# 5. Tabela de medias G x A
# -----------------------------------------------------------------------------

# Tabela de medias fenotipicas:
# linhas  = genotipos
# colunas = ambientes
Y <- tapply(dados$y, list(dados$gen, dados$env), mean)

if (anyNA(Y)) {
  stop("A tabela de medias G x A contem valores NA.")
}

nomes_env <- colnames(Y)
nomes_gen <- rownames(Y)

cat("\nTabela de medias G x A (kg/ha):\n")
print(round(Y, 1))

# Medias marginais e media geral
media_gen   <- rowMeans(Y)   # media de cada genotipo
media_env   <- colMeans(Y)   # media de cada ambiente
media_geral <- mean(Y)       # media geral

cat("\nMedias por ambiente:\n")
print(round(media_env, 1))

cat("\nMedias por genotipo:\n")
print(round(media_gen, 1))

# -----------------------------------------------------------------------------
# 6. Grafico de interacao
# -----------------------------------------------------------------------------

# Cada linha representa um genotipo.
# Muitos cruzamentos sugerem interacao complexa.
mat_plot <- t(Y)

plot(
  NA,
  xlim = c(1, nrow(mat_plot)),
  ylim = range(mat_plot, na.rm = TRUE),
  xaxt = "n",
  xlab = "Ambiente",
  ylab = "Media (kg/ha)",
  main = "Grafico de interacao G x A"
)
axis(1, at = seq_len(nrow(mat_plot)), labels = rownames(mat_plot))
matlines(seq_len(nrow(mat_plot)), mat_plot, col = "gray70", lty = 1, lwd = 0.8)
matpoints(seq_len(nrow(mat_plot)), mat_plot, col = "gray40", pch = 16, cex = 0.5)
abline(h = mean(dados$y), lty = 2, col = "red", lwd = 1.2)

# -----------------------------------------------------------------------------
# 7. ANOVA conjunta (DBC, ambientes fixos)
#    Necessaria para obter QME -- denominador do teste de Shukla.
# -----------------------------------------------------------------------------

ajuste_aov <- aov(y ~ gen * env + block %in% env, data = dados)
tab_aov    <- summary(ajuste_aov)[[1]]

cat("\nANOVA conjunta:\n")
print(tab_aov)

QME    <- tab_aov["Residuals", "Mean Sq"]
gl_erro <- tab_aov["Residuals", "Df"]

cat(sprintf("\nQME (erro experimental) = %.4f  (GL = %d)\n", QME, gl_erro))

# -----------------------------------------------------------------------------
# 8. Quantidades basicas comuns a todos os indices
# -----------------------------------------------------------------------------

# Indice ambiental (mesma ideia usada na regressao conjunta):
# positivo = favoravel, negativo = desfavoravel, zero = neutro.
I_j <- media_env - media_geral

# Efeito de interacao aditivo-residual:
# media observada da celula - media do genotipo - media do ambiente + media geral
GA_hat <- Y -
          outer(media_gen, rep(1, J)) -
          outer(rep(1, I), media_env) + media_geral

cat("\nEfeitos de interacao aditivo-residual (GA_hat_ij):\n")
print(round(GA_hat, 2))

# Por construcao: soma em linha = 0 e soma em coluna = 0
cat("\nChecagem das restricoes de GA_hat:\n")
cat("  Max |soma por linha|  :", round(max(abs(rowSums(GA_hat))), 10), "\n")
cat("  Max |soma por coluna| :", round(max(abs(colSums(GA_hat))), 10), "\n")

# Subconjuntos ambientais
amb_fav  <- which(I_j > 0)
amb_desf <- which(I_j < 0)
amb_zero <- which(I_j == 0)

cat("\nIndice ambiental (I_j):\n")
print(round(I_j, 2))

cat("\nAmbientes favoraveis   :", if (length(amb_fav)  > 0) paste(nomes_env[amb_fav],  collapse = ", ") else "nenhum", "\n")
cat("Ambientes desfavoraveis:", if (length(amb_desf) > 0) paste(nomes_env[amb_desf], collapse = ", ") else "nenhum", "\n")
cat("Ambientes neutros      :", if (length(amb_zero) > 0) paste(nomes_env[amb_zero], collapse = ", ") else "nenhum", "\n")

# =============================================================================
# INDICES DE ESTABILIDADE
# =============================================================================
# Organizados em duas familias:
#   A) Indices centrados na contribuicao do genotipo a G x A
#      8.1 Metodo tradicional
#      8.2 Plaisted & Peterson
#      8.3 Wricke
#      8.4 Shukla
#   B) Indices que combinam media e estabilidade
#      8.5 Francis & Kannenberg
#      8.6 Lin & Binns
#      8.7 Annicchiarico
# =============================================================================

# -----------------------------------------------------------------------------
# 8.1 Metodo tradicional: QM(A/G_i)
# -----------------------------------------------------------------------------

QM_AGi <- apply(Y, 1, var)

tab_tradicional <- data.frame(
  Genotipo = nomes_gen,
  Media    = round(media_gen, 1),
  QM_AGi   = round(QM_AGi, 2),
  row.names = NULL
)

cat("\n--- Metodo tradicional: QM(A/G_i) ---\n")
print(tab_tradicional[order(tab_tradicional$QM_AGi), ], row.names = FALSE)

# -----------------------------------------------------------------------------
# 8.2 Plaisted & Peterson (theta_i)
# -----------------------------------------------------------------------------

# Para cada par de genotipos (i, i'), calcula-se o componente de variancia
# da interacao do par, e depois tira-se a media sobre todos os pares que
# envolvem o genotipo i.
theta_i <- numeric(I)
names(theta_i) <- nomes_gen

for (ii in seq_len(I)) {
  sigma_par <- numeric(I - 1)
  k_par <- 0

  for (ip in seq_len(I)) {
    if (ip == ii) next

    k_par <- k_par + 1
    d2_ii <- sum((Y[ii, ] - Y[ip, ])^2)

    # Usando medias por ambiente:
    # 1/J * (Y_i. - Y_i'.)^2  ==  J * (media_i - media_i')^2
    SQ_par <- (K / 2) * (d2_ii - J * (media_gen[ii] - media_gen[ip])^2)
    sigma_par[k_par] <- (SQ_par / (J - 1) - QME) / K
  }

  theta_i[ii] <- mean(sigma_par)
}

theta_i[theta_i < 0] <- 0

soma_theta <- sum(theta_i)
theta_pct  <- if (soma_theta > 0) 100 * theta_i / soma_theta else rep(NA_real_, I)

tab_pp <- data.frame(
  Genotipo = nomes_gen,
  Media    = round(media_gen, 1),
  theta_i  = round(theta_i, 4),
  theta_p  = round(theta_pct, 2),
  row.names = NULL
)

cat("\n--- Plaisted & Peterson (theta_i) ---\n")
print(tab_pp[order(tab_pp$theta_i), ], row.names = FALSE)

# -----------------------------------------------------------------------------
# 8.3 Ecovalencia de Wricke (omega_i)
# -----------------------------------------------------------------------------

omega_i <- rowSums(GA_hat^2)
SQ_GxA_medias <- sum(omega_i)

# Na escala da ANOVA, a SQ da interacao e multiplicada por K
SQ_GxA_aov_via_medias <- K * SQ_GxA_medias
SQ_GxA_aov_observada  <- tab_aov["gen:env", "Sum Sq"]

tab_wricke <- data.frame(
  Genotipo = nomes_gen,
  Media    = round(media_gen, 1),
  omega_i  = round(omega_i, 4),
  omega_p  = round(100 * omega_i / SQ_GxA_medias, 2),
  row.names = NULL
)

cat("\n--- Ecovalencia de Wricke (omega_i) ---\n")
print(tab_wricke[order(tab_wricke$omega_i), ], row.names = FALSE)

cat(sprintf("\nChecagem da identidade de Wricke:\n"))
cat(sprintf("  sum(omega_i) [escala da tabela de medias] = %.6f\n", SQ_GxA_medias))
cat(sprintf("  K * sum(omega_i) [escala da ANOVA]        = %.6f\n", SQ_GxA_aov_via_medias))
cat(sprintf("  SQ(gen:env) da ANOVA conjunta             = %.6f\n", SQ_GxA_aov_observada))
cat(sprintf("  Diferenca (ANOVA - K*sum(omega_i))        = %.10f\n",
            SQ_GxA_aov_observada - SQ_GxA_aov_via_medias))

# -----------------------------------------------------------------------------
# 8.4 Variancia de estabilidade de Shukla (sigma2_i)
# -----------------------------------------------------------------------------

sigma2_i <- (I * (I - 1) * omega_i - SQ_GxA_medias) /
            ((I - 1) * (I - 2) * (J - 1))
sigma2_i[sigma2_i < 0] <- 0

# Teste F aproximado usando QME da ANOVA conjunta
QM_int_i <- K * omega_i / (J - 1)
F_i      <- QM_int_i / QME
gl1_sh   <- J - 1
gl2_sh   <- (I - 1) * (I - 2) * (J - 1) / 2
p_i      <- pf(F_i, df1 = gl1_sh, df2 = gl2_sh, lower.tail = FALSE)

tab_shukla <- data.frame(
  Genotipo = nomes_gen,
  Media    = round(media_gen, 1),
  omega_i  = round(omega_i, 4),
  sigma2_i = round(sigma2_i, 4),
  F_stat   = round(F_i, 4),
  p_valor  = sapply(p_i, formata_p),
  row.names = NULL
)

cat("\n--- Variancia de estabilidade de Shukla (sigma2_i) ---\n")
print(tab_shukla[order(tab_shukla$sigma2_i), ], row.names = FALSE)

cat(sprintf("\nCorrelacoes de postos de Spearman:\n"))
cat(sprintf("  Plaisted & Peterson vs Wricke = %.4f\n",
            cor(theta_i, omega_i, method = "spearman")))
cat(sprintf("  Wricke vs Shukla              = %.4f\n",
            cor(omega_i, sigma2_i, method = "spearman")))

# -----------------------------------------------------------------------------
# 8.5 Francis & Kannenberg (CV_i)
# -----------------------------------------------------------------------------

CV_i <- apply(Y, 1, sd) / media_gen * 100

tab_fk <- data.frame(
  Genotipo = nomes_gen,
  Media    = round(media_gen, 1),
  CV_i     = round(CV_i, 2),
  row.names = NULL
)

plot(
  media_gen, CV_i,
  xlab = "Media geral do genotipo (kg/ha)",
  ylab = expression(CV[i] ~ "(%)"),
  main = "Francis & Kannenberg",
  pch  = 16,
  col  = "gray30",
  xlim = ajusta_faixa(media_gen),
  ylim = ajusta_faixa(CV_i)
)
abline(h = mean(CV_i), v = mean(media_gen), lty = 2, col = "red")
text(media_gen, CV_i, labels = nomes_gen, cex = 0.50, pos = 3)

cat("\n--- Francis & Kannenberg (CV_i) ---\n")
print(tab_fk[order(tab_fk$CV_i), ], row.names = FALSE)

# -----------------------------------------------------------------------------
# 8.6 Lin & Binns (P_i): geral, favoravel, desfavoravel
# -----------------------------------------------------------------------------

P_i    <- indice_lin_binns(Y)
P_fav  <- indice_lin_binns(Y, amb_fav)
P_desf <- indice_lin_binns(Y, amb_desf)

tab_lb <- data.frame(
  Genotipo = nomes_gen,
  Media    = round(media_gen, 1),
  P_i      = round(P_i, 2),
  P_fav    = round(P_fav, 2),
  P_desf   = round(P_desf, 2),
  row.names = NULL
)

cat("\n--- Lin & Binns (P_i) ---\n")
print(tab_lb[order(tab_lb$P_i), ], row.names = FALSE)

# -----------------------------------------------------------------------------
# 8.7 Annicchiarico (W_i): geral, favoravel, desfavoravel
# -----------------------------------------------------------------------------

alpha <- 0.25

W_i    <- indice_annicchiarico(Y, alpha = alpha)
W_fav  <- indice_annicchiarico(Y, amb_fav, alpha = alpha)
W_desf <- indice_annicchiarico(Y, amb_desf, alpha = alpha)

tab_ann <- data.frame(
  Genotipo = nomes_gen,
  Media    = round(media_gen, 1),
  W_i      = round(W_i, 2),
  W_fav    = round(W_fav, 2),
  W_desf   = round(W_desf, 2),
  row.names = NULL
)

cat("\n--- Annicchiarico (W_i, alpha = 0.25) ---\n")
print(tab_ann[order(tab_ann$W_i, decreasing = TRUE), ], row.names = FALSE)

# -----------------------------------------------------------------------------
# 9. Resumo final integrado
# -----------------------------------------------------------------------------

resultado_final <- data.frame(
  Genotipo = nomes_gen,
  Media    = round(media_gen, 1),
  QM_AGi   = round(QM_AGi, 2),
  theta_i  = round(theta_i, 4),
  omega_i  = round(omega_i, 4),
  sigma2_i = round(sigma2_i, 4),
  CV_i     = round(CV_i, 2),
  P_i      = round(P_i, 2),
  W_i      = round(W_i, 2),
  row.names = NULL
)

cat("\n=============================================================\n")
cat("RESUMO FINAL -- Aula 4: Indices Classicos de Estabilidade\n")
cat("=============================================================\n")
cat(sprintf("Genotipos              : %d\n", I))
cat(sprintf("Ambientes              : %d\n", J))
cat(sprintf("Blocos por ambiente    : %d\n", K))
cat(sprintf("QME (ANOVA conjunta)   : %.4f\n", QME))
cat(sprintf("SQ_GxA (escala medias) : %.6f\n", SQ_GxA_medias))
cat(sprintf("SQ_GxA (escala ANOVA)  : %.6f\n", SQ_GxA_aov_via_medias))
cat(sprintf("Amb. favoraveis        : %s\n",
            if (length(amb_fav) > 0) paste(nomes_env[amb_fav], collapse = ", ") else "nenhum"))
cat(sprintf("Amb. desfavoraveis     : %s\n",
            if (length(amb_desf) > 0) paste(nomes_env[amb_desf], collapse = ", ") else "nenhum"))
cat(sprintf("Amb. neutros           : %s\n",
            if (length(amb_zero) > 0) paste(nomes_env[amb_zero], collapse = ", ") else "nenhum"))

cat("\nTabela integrada (ordenada por media decrescente):\n")
print(resultado_final[order(resultado_final$Media, decreasing = TRUE), ], row.names = FALSE)

cat("\nTop 5 por Wricke (menor omega_i):\n")
print(tab_wricke[order(tab_wricke$omega_i), ][1:5, ], row.names = FALSE)

cat("\nTop 5 por Lin & Binns (menor P_i):\n")
print(tab_lb[order(tab_lb$P_i), ][1:5, ], row.names = FALSE)

cat("\nTop 5 por Annicchiarico (maior W_i):\n")
print(tab_ann[order(tab_ann$W_i, decreasing = TRUE), ][1:5, ], row.names = FALSE)

cat("=============================================================\n")
