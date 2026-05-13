# =============================================================================
# GMP0164 - Biometria da Interacao Genotipos x Ambientes
# Aula 10: Predicao de Desempenho em Multiplos Ambientes
# Dados: Soja MET simulado (40 genotipos x 6 ambientes x 3 blocos, DBC)
# Prof. Dr. Rafael Tassinari Resende
# =============================================================================
# Dependencias:
#   install.packages(c("missForest", "lme4breeding"))
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Carregar dados
# -----------------------------------------------------------------------------

# Gerar os dados (caso ainda nao tenha rodado o simulador):
# source("simu_soybean-data_balanced.R")

dat <- read.table("soy_MET_phenotypic.txt", header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE)

head(dat)
str(dat)


# -----------------------------------------------------------------------------
# 2. Preprocessamento
# -----------------------------------------------------------------------------

dat$gen   <- factor(dat$gen)
dat$env   <- factor(dat$env)
dat$block <- factor(dat$block)
dat$group <- factor(dat$group,
                    levels = c("check", "stable", "responsive",
                               "rustic", "unstable", "inferior"))

summary(dat[ , c("gen", "env", "block", "group", "y")])


# -----------------------------------------------------------------------------
# 3. Funcoes auxiliares
# -----------------------------------------------------------------------------

# Formatar p-valor
fmt_p <- function(p) {
  ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p))
}

# SVD iterativo (Yan, 2013)
svd_impute <- function(M, ncomp = 2, maxiter = 100, tol = 1e-4) {
  M_imp <- M
  col_means <- colMeans(M, na.rm = TRUE)
  for (j in seq_len(ncol(M_imp)))
    M_imp[is.na(M_imp[ , j]), j] <- col_means[j]

  for (iter in seq_len(maxiter)) {
    M_old <- M_imp
    sv    <- svd(M_imp)
    M_imp <- sv$u[ , 1:ncomp, drop = FALSE] %*%
             diag(sv$d[1:ncomp], ncomp)     %*%
             t(sv$v[ , 1:ncomp, drop = FALSE])
    dimnames(M_imp) <- dimnames(M)
    M_imp[!is.na(M)] <- M[!is.na(M)]
    if (max(abs(M_imp - M_old)) < tol) {
      cat(sprintf("  SVD convergiu em %d iteracoes.\n", iter))
      break
    }
  }
  M_imp
}

# Acuracia preditiva (correlacao de Pearson)
rho_pred <- function(pred, obs) cor(pred, obs, use = "complete.obs")


# -----------------------------------------------------------------------------
# 4. Dimensoes experimentais
# -----------------------------------------------------------------------------

I <- nlevels(dat$gen)
J <- nlevels(dat$env)
K <- nlevels(dat$block)

cat("\nDimensoes experimentais:\n")
cat("  Genotipos   :", I, "\n")
cat("  Ambientes   :", J, "\n")
cat("  Blocos      :", K, "\n")
cat("  Total obs.  :", nrow(dat), "\n")

cat("\nComposicao dos grupos:\n")
print(table(unique(dat[ , c("gen", "group")])$group))


# -----------------------------------------------------------------------------
# 5. Sumarios descritivos e matriz G x A
# -----------------------------------------------------------------------------

mat_obs <- tapply(dat$y, list(dat$gen, dat$env), mean)

cat("\nMatriz G x A (medias) — primeiras linhas:\n")
print(round(mat_obs[1:8, ], 0))

cat("\nMedias por ambiente:\n")
print(round(colMeans(mat_obs), 0))

cat("\nMedias por genotipo (ordenadas):\n")
print(round(sort(rowMeans(mat_obs), decreasing = TRUE), 0))

# --- Grafico de interacao ----------------------------------------------------

grp_lookup <- unique(dat[ , c("gen", "group")])
rownames(grp_lookup) <- grp_lookup$gen

cols_grp <- c(check = "black", stable = "royalblue", responsive = "forestgreen",
              rustic = "darkorange", unstable = "red3", inferior = "gray50")

mat_plot <- t(mat_obs)

plot(NA,
     xlim = c(1, nrow(mat_plot)), ylim = range(mat_plot, na.rm = TRUE),
     xaxt = "n", xlab = "Ambiente", ylab = "Produtividade media (kg/ha)",
     main = "Grafico de interacao G x A — dados completos")
axis(1, at = 1:nrow(mat_plot), labels = rownames(mat_plot), las = 2, cex.axis = 0.7)
abline(h = mean(mat_obs), lty = 2, col = "gray70")

for (i in seq_len(ncol(mat_plot))) {
  gname <- colnames(mat_plot)[i]
  grp   <- as.character(grp_lookup[gname, "group"])
  lines(1:nrow(mat_plot),  mat_plot[ , i], col = cols_grp[grp], lwd = 1)
  points(1:nrow(mat_plot), mat_plot[ , i], col = cols_grp[grp], pch = 16, cex = 0.5)
}

legend("topleft", legend = names(cols_grp), col = cols_grp,
       lwd = 1.5, pch = 16, cex = 0.75, bty = "n")


# =============================================================================
# 6. Simulacao de desbalanceamento (esquema CV2)
# =============================================================================

set.seed(42)
prop_miss <- 0.20
n_miss    <- round(prop_miss * length(mat_obs))
idx_miss  <- sample(seq_along(mat_obs), n_miss)
obs_true  <- mat_obs[idx_miss]   # valores verdadeiros — nao tocar!

mat_miss             <- mat_obs
mat_miss[idx_miss]   <- NA

cat(sprintf("\nCV2 — celulas removidas: %d (%.0f%% da matriz)\n",
            n_miss, 100 * prop_miss))
cat("Celulas faltantes por ambiente:\n")
print(colSums(is.na(mat_miss)))


# =============================================================================
# 7. Imputacao por SVD iterativo (Yan, 2013)
# =============================================================================

cat("\n--- SVD iterativo ---\n")

mat_svd  <- svd_impute(mat_miss, ncomp = 2)
pred_svd <- mat_svd[idx_miss]
rho_svd  <- rho_pred(pred_svd, obs_true)
cat(sprintf("  Acuracia preditiva (rho_pred) : %.3f\n", rho_svd))

cat("\n  Sensibilidade ao numero de componentes (ncomp):\n")
for (k in 1:5) {
  r_k <- rho_pred(svd_impute(mat_miss, ncomp = k)[idx_miss], obs_true)
  cat(sprintf("    ncomp = %d  ->  rho_pred = %.3f\n", k, r_k))
}


# =============================================================================
# 8. Imputacao por missForest
# =============================================================================

library(missForest)

cat("\n--- missForest ---\n")

imp_mf  <- missForest(mat_miss, verbose = FALSE)
pred_mf <- imp_mf$ximp[idx_miss]
rho_mf  <- rho_pred(pred_mf, obs_true)

cat(sprintf("  Acuracia preditiva (rho_pred) : %.3f\n", rho_mf))
cat(sprintf("  Erro OOB estimado (NRMSE)     : %.4f\n", imp_mf$OOBerror))


# =============================================================================
# 9. Imputacao por modelo misto com covariancia UN (lme4breeding)
# =============================================================================
# O modelo misto com estrutura de covariancia nao estruturada (UN) entre
# ambientes e o "imputador natural": a predicao de celulas faltantes emerge
# diretamente das equacoes de modelo misto (BLUPs), ponderada pela
# covariancia genetica entre os ambientes.
#
# Nota sobre dados desbalanceados e lme4breeding:
#   lmebreed() usa decomposicao espectral da matriz de relacao, que requer
#   dados balanceados. Para dados desbalanceados (nosso CV2), a abordagem
#   recomendada e a "imputation trick" (Covarrubias-Pazaran, 2024): preencher
#   as celulas faltantes com a media da coluna antes de ajustar o modelo.
#   Os BLUPs resultantes sao altamente correlacionados (>0.98) com os
#   obtidos pela formulacao densa completa — e ainda melhores que a media
#   sozinha, pois a estrutura UN redistribui a informacao entre ambientes.
# =============================================================================

library(lme4breeding)

cat("\n--- lme4breeding: modelo misto com covariancia UN ---\n")

# --- 9.1 Preparar dados: mean-imputation das celulas faltantes ---------------
# (necessario para que lmebreed utilize a decomposicao espectral)

mat_init <- mat_miss
for (j in seq_len(ncol(mat_init)))
  mat_init[is.na(mat_init[ , j]), j] <- mean(mat_init[ , j], na.rm = TRUE)

# Converter para long format (uma linha por genotipo x ambiente)
dat_un     <- as.data.frame(as.table(mat_init))
colnames(dat_un) <- c("gen", "env", "y")
dat_un$gen <- factor(dat_un$gen, levels = levels(dat$gen))
dat_un$env <- factor(dat_un$env, levels = levels(dat$env))

# --- 9.2 Ajustar modelo UN ---------------------------------------------------
# Formula: y ~ env + (0 + env | gen)
#   Efeitos fixos  : env — captura a media de cada ambiente
#   Efeitos aleatorios : (0 + env | gen) — cada genotipo tem um efeito
#     especifico por ambiente; a covariancia entre esses efeitos e modelada
#     por uma matriz UN (J x J), com J*(J+1)/2 parametros estimados via REML.
#     Isso e equivalente a assumir que cada par de ambientes tem sua propria
#     correlacao genetica, capturando a heterogeneidade da G x A.

fit_un <- lmebreed(
  y ~ env + (0 + env | gen),
  data    = dat_un,
  verbose = FALSE
)

cat("  Modelo ajustado. Componentes de variancia (matriz UN):\n")
print(VarCorr(fit_un), comp = "Variance")

# Correlacoes geneticas entre ambientes (estrutura UN estimada)
vc_un   <- as.matrix(VarCorr(fit_un)$gen)
corr_un <- cov2cor(vc_un)
cat("\n  Correlacoes geneticas entre ambientes:\n")
print(round(corr_un, 2))

# --- 9.3 Extrair BLUPs e montar matriz imputada ------------------------------

re_un <- ranef(fit_un)$gen    # genotipos (linhas) x env (colunas)
fe_un <- fixef(fit_un)        # efeitos fixos

# Efeito fixo por ambiente:
#   nivel de referencia = intercepto; demais = intercepto + desvio estimado
env_levels <- levels(dat_un$env)
fe_env_un  <- setNames(numeric(J), env_levels)
fe_env_un[env_levels[1]] <- fe_un["(Intercept)"]
for (j in 2:J) {
  key <- paste0("env", env_levels[j])
  fe_env_un[env_levels[j]] <- fe_un["(Intercept)"] + fe_un[key]
}

# BLUP(i,j) = efeito fixo(j) + efeito aleatorio genotipico(i,j)
re_cols <- colnames(re_un)   # nomes das colunas do ranef (conferir)
mat_un  <- matrix(NA, I, J,
                  dimnames = list(rownames(mat_obs), colnames(mat_obs)))

for (j in seq_len(J))
  mat_un[ , j] <- fe_env_un[env_levels[j]] +
                  re_un[rownames(mat_obs), re_cols[j]]

pred_un <- mat_un[idx_miss]
rho_un  <- rho_pred(pred_un, obs_true)
cat(sprintf("\n  Acuracia preditiva (rho_pred) : %.3f\n", rho_un))


# =============================================================================
# 10. Comparacao de acuracia preditiva (CV2)
# =============================================================================

# Baseline: imputacao simples pela media da coluna
baseline_pred <- colMeans(mat_miss, na.rm = TRUE)[col(mat_miss)[idx_miss]]
rho_baseline  <- rho_pred(baseline_pred, obs_true)

cat("\n=== Acuracia preditiva CV2 ===\n")
cat(sprintf("  Media da coluna (baseline)      : rho_pred = %.3f\n", rho_baseline))
cat(sprintf("  SVD iterativo (ncomp = 2)       : rho_pred = %.3f\n", rho_svd))
cat(sprintf("  missForest                      : rho_pred = %.3f\n", rho_mf))
cat(sprintf("  lme4breeding UN                 : rho_pred = %.3f\n", rho_un))

# Grafico: predito vs observado (3 metodos)
par(mfrow = c(1, 3))
lim <- range(c(obs_true, pred_svd, pred_mf, pred_un))

for (lst in list(
  list(p = pred_svd, r = rho_svd, col = "royalblue",   lab = "SVD iterativo"),
  list(p = pred_mf,  r = rho_mf,  col = "forestgreen", lab = "missForest"),
  list(p = pred_un,  r = rho_un,  col = "darkorange",  lab = "lme4breeding UN")
)) {
  plot(obs_true, lst$p,
       pch = 16, col = adjustcolor(lst$col, 0.6), cex = 0.8,
       xlim = lim, ylim = lim,
       xlab = "Observado (kg/ha)", ylab = "Predito (kg/ha)",
       main = sprintf("%s\nrho = %.3f", lst$lab, lst$r))
  abline(0, 1, col = "gray40", lty = 2)
}

par(mfrow = c(1, 1))


# =============================================================================
# 11. Uso estrategico: ranking com matriz completa vs incompleta
# =============================================================================
# Usa a matriz imputada pelo modelo UN como referencia, por ter o
# melhor fundamento probabilistico (covariancia genetica entre ambientes).
# =============================================================================

cat("\n=== Ranking: matriz incompleta vs matriz UN ===\n")

mean_incomplete <- rowMeans(mat_miss, na.rm = TRUE)
mean_un         <- rowMeans(mat_un)

rank_inc <- rank(-mean_incomplete)
rank_un  <- rank(-mean_un)
delta    <- rank_un - rank_inc

ranking_tab <- data.frame(
  genotipo   = rownames(mat_obs),
  grupo      = grp_lookup[rownames(mat_obs), "group"],
  media_obs  = round(mean_incomplete, 0),
  media_un   = round(mean_un, 0),
  rank_obs   = rank_inc,
  rank_un    = rank_un,
  delta_rank = delta,
  stringsAsFactors = FALSE
)
ranking_tab <- ranking_tab[order(ranking_tab$rank_un), ]

cat("\nTop-10 genotipos (ranking pela matriz imputada UN):\n")
print(ranking_tab[1:10, ], row.names = FALSE)

cat("\nGenotipos com maior mudanca de ranking (|delta| >= 5):\n")
big_shift <- ranking_tab[abs(ranking_tab$delta_rank) >= 5, ]
if (nrow(big_shift) > 0) print(big_shift, row.names = FALSE) else
  cat("  Nenhum genotipo com |delta| >= 5.\n")

rho_rank <- cor(rank_inc, rank_un, method = "spearman")
cat(sprintf("\nCorrelacao de Spearman (ranking incompleto vs UN): %.3f\n", rho_rank))

plot(rank_inc, rank_un,
     pch = 16, cex = 0.9,
     col = cols_grp[as.character(
       grp_lookup[rownames(mat_obs)[order(rownames(mat_obs))], "group"])],
     xlab = "Ranking — matriz incompleta",
     ylab = "Ranking — matriz imputada (UN)",
     main = sprintf("Mudanca de ranking apos imputacao UN  |  rho_S = %.3f", rho_rank))
abline(0, 1, lty = 2, col = "gray50")
legend("topleft", legend = names(cols_grp), col = cols_grp,
       pch = 16, cex = 0.75, bty = "n")


# =============================================================================
# [OPCIONAL] 12. Esquema CV1: predicao de ambiente nao observado (SVD)
# =============================================================================

cat("\n=== CV1: predicao de ambiente nao observado (SVD iterativo) ===\n")

env_levels_obs <- colnames(mat_obs)
rho_cv1        <- setNames(numeric(J), env_levels_obs)

for (env_j in env_levels_obs) {
  mat_cv1           <- mat_obs
  mat_cv1[ , env_j] <- NA
  mat_imp_j         <- svd_impute(mat_cv1, ncomp = 2)
  rho_cv1[env_j]    <- rho_pred(mat_imp_j[ , env_j], mat_obs[ , env_j])
}

cat("\nAcuracia preditiva CV1 por ambiente (SVD, ncomp = 2):\n")
for (env_j in env_levels_obs)
  cat(sprintf("  %-24s  rho_pred = %.3f\n", env_j, rho_cv1[env_j]))
cat(sprintf("  Media CV1: %.3f\n", mean(rho_cv1)))

barplot(rho_cv1,
        ylim = c(0, 1), col = "steelblue", las = 2, cex.names = 0.65,
        ylab = "Acuracia preditiva (rho_pred)",
        main = "CV1: acuracia por ambiente removido (SVD iterativo)")
abline(h = mean(rho_cv1), lty = 2, col = "red3")
legend("topright", legend = sprintf("Media = %.2f", mean(rho_cv1)),
       lty = 2, col = "red3", bty = "n", cex = 0.9)


# -----------------------------------------------------------------------------
# 13. Resumo final
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("RESUMO FINAL — Aula 10\n")
cat("=============================================================\n")
cat(sprintf("Genotipos              : %d\n", I))
cat(sprintf("Ambientes              : %d\n", J))
cat(sprintf("Blocos                 : %d\n", K))
cat(sprintf("Celulas faltantes (CV2): %d (%.0f%%)\n", n_miss, 100 * prop_miss))
cat("--- Acuracia preditiva (CV2) ---\n")
cat(sprintf("  Media da coluna (baseline)    : %.3f\n", rho_baseline))
cat(sprintf("  SVD iterativo (ncomp = 2)     : %.3f\n", rho_svd))
cat(sprintf("  missForest                    : %.3f\n", rho_mf))
cat(sprintf("  lme4breeding UN               : %.3f\n", rho_un))
cat("--- Acuracia preditiva (CV1, SVD) ---\n")
cat(sprintf("  Media por ambiente            : %.3f\n", mean(rho_cv1)))
cat("--- Ranking ---\n")
cat(sprintf("  Spearman (incompleto vs UN)   : %.3f\n", rho_rank))
cat(sprintf("  Genotipos com |delta| >= 5    : %d\n",
            sum(abs(ranking_tab$delta_rank) >= 5)))
cat("=============================================================\n")