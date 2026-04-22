# =============================================================================
# GMP0164 - Biometria da Interacao Genotipos x Ambientes
# Aula 7: Modelos Mistos II -- Estruturas de Variancia-Covariancia
# Dados: Milho safrinha MET 2025 -- 70 genotipos, 6 ensaios, 3 blocos/ensaio
#        Desbalanceado entre ensaios; balanceado dentro de cada ensaio
#        Variavel resposta: yield (kg/ha)
# Prof. Dr. Rafael Tassinari Resende
# =============================================================================
#
# Estrategia de software:
#
#   ID   (lme4)   : (1 | gen) + (1 | gen:trial)
#                   separa um componente comum de genotipo e um componente
#                   adicional de genotipo dentro de ensaio
#                   G_e = sigma^2_g J_J + sigma^2_ge I_J  (2 parametros VCOV)
#                   na pratica, induz uma estrutura homogenea do tipo CS
#
#   CS   (nlme)   : random = list(gen = pdCompSymm(~ trial - 1))
#                   mesma variancia genetica em todos os ensaios
#                   mesma covariancia genetica entre quaisquer dois ensaios
#                   G_e com diagonal constante e fora da diagonal constante
#                   (2 parametros VCOV)
#
#   UN   (lme4)   : (0 + trial | gen)
#                   G_e livre J x J  (J*(J+1)/2 = 21 parametros VCOV)
#
#   DIAG (sommer) : vsm(dsm(trial), ism(gen))
#                   G_e = diag(sigma^2_g1,...,sigma^2_gJ)  (J parametros VCOV)
#
#   UN   (sommer) : vsm(usm(trial), ism(gen))
#                   G_e livre J x J -- para comparacao com UN do lme4
#
# Observacoes:
#   1) Blocos dentro de ensaio entram como efeitos fixos em todos os modelos.
#   2) Os criterios de ajuste do sommer serao mostrados em tabela separada,
#      pois a escala reportada por mmes nao ficou diretamente comparavel com
#      lme4/nlme neste script.
#   3) O modelo FA foi deixado de fora.
# =============================================================================

library(lme4)
library(nlme)
library(sommer)
library(ggplot2)
library(reshape2)

# =============================================================================
# 1. Leitura dos dados
# =============================================================================

dat <- read.table(
  "G:/Meu Drive/UFG/PPGGMP/BiometriaDaInteracaoGxA/aulas/simu_data/maize_safrinha_MET.txt",
  header = TRUE, sep = "\t"
)

# =============================================================================
# 2. Pre-processamento
# =============================================================================

dat$gen   <- factor(dat$gen)
dat$trial <- factor(dat$trial)
dat$block <- factor(dat$block)
names(dat)[names(dat) == "yield"] <- "y"

dat <- dat[order(dat$trial), ]
rownames(dat) <- NULL

tnames <- levels(dat$trial)
gnames <- levels(dat$gen)
tlab   <- sub("TRL_", "", sub("_2025", "", tnames))

I <- nlevels(dat$gen)
J <- nlevels(dat$trial)

# =============================================================================
# 3. Funcoes auxiliares
# =============================================================================

vc_sommer_fit <- function(fit) c(
  logLik = as.numeric(fit$llik[1, ncol(fit$llik)]),
  AIC    = as.numeric(fit$AIC),
  BIC    = as.numeric(fit$BIC)
)

get_G_sommer <- function(fit, tnames) {
  vc <- summary(fit)$varcomp
  rn <- rownames(vc)[grepl("^trial:gen:", rownames(vc))]
  G  <- matrix(0, length(tnames), length(tnames), dimnames = list(tnames, tnames))
  for (i in rn) {
    ij <- strsplit(sub("^trial:gen:", "", i), ":", fixed = TRUE)[[1]]
    G[ij[1], ij[2]] <- vc[i, "VarComp"]
    G[ij[2], ij[1]] <- vc[i, "VarComp"]
  }
  G
}

plot_vcov <- function(Glist, tnames) {
  df <- do.call(rbind, lapply(names(Glist), function(nm) {
    x <- melt(Glist[[nm]])
    names(x) <- c("Env1", "Env2", "value")
    x$modelo <- nm
    x
  }))
  lev <- sub("TRL_", "", sub("_2025", "", tnames))
  df$Env1 <- factor(sub("TRL_", "", sub("_2025", "", df$Env1)), levels = rev(lev))
  df$Env2 <- factor(sub("TRL_", "", sub("_2025", "", df$Env2)), levels = lev)
  lim <- 0.6 * max(abs(df$value), na.rm = TRUE)
  ggplot(df, aes(Env2, Env1, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.0f", value), color = abs(value) > lim), size = 3) +
    scale_fill_gradient2(low = "#d73027", mid = "white", high = "#4575b4", midpoint = 0, name = "VCOV") +
    scale_color_manual(values = c("black", "white"), guide = "none") +
    facet_wrap(~ modelo, nrow = 2, ncol = 3) +
    labs(title = "Matrizes de variancia-covariancia genetica", x = NULL, y = NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(hjust = 0.5),
      panel.grid  = element_blank()
    )
}

# =============================================================================
# 4. Balanceamento e descritivas
# =============================================================================

tab_bal <- table(dat$gen, dat$trial)

print(colSums(tab_bal > 0))
print(table(rowSums(tab_bal > 0)))

means_ge <- tapply(dat$y, list(dat$gen, dat$trial), mean)

desc_env <- data.frame(
  ensaio = tlab,
  media  = round(colMeans(means_ge, na.rm = TRUE), 1),
  min    = round(apply(means_ge, 2, min, na.rm = TRUE), 1),
  max    = round(apply(means_ge, 2, max, na.rm = TRUE), 1),
  dp     = round(apply(means_ge, 2, sd,  na.rm = TRUE), 1)
)
print(desc_env, row.names = FALSE)

gm <- sort(rowMeans(means_ge, na.rm = TRUE), decreasing = TRUE)
print(round(head(gm, 5), 1))

# =============================================================================
# 5. Grafico de interacao
# =============================================================================

plot_long <- reshape(
  data.frame(gen = rownames(means_ge), means_ge, check.names = FALSE),
  varying   = colnames(means_ge),
  v.names   = "y",
  timevar   = "trial",
  times     = colnames(means_ge),
  direction = "long"
)

plot_long$trial_lab <- factor(
  sub("TRL_", "", sub("_2025", "", plot_long$trial)),
  levels = tlab[order(colMeans(means_ge, na.rm = TRUE))]
)

trial_means <- aggregate(y ~ trial_lab, data = plot_long, mean, na.rm = TRUE)
trial_means$trial_lab <- factor(trial_means$trial_lab, levels = levels(plot_long$trial_lab))

ggplot(plot_long, aes(trial_lab, y, group = gen)) +
  geom_line(color = "gray40", alpha = 0.28, linewidth = 0.40, na.rm = TRUE) +
  geom_point(color = "gray40", alpha = 0.38, size = 0.85, na.rm = TRUE) +
  geom_point(
    data = trial_means, aes(trial_lab, y),
    inherit.aes = FALSE, color = "firebrick", shape = 18, size = 3.5
  ) +
  labs(
    title = "Grafico de interacao G x A -- milho safrinha 2025",
    x = "Ensaio (local)", y = "Produtividade media (kg/ha)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )

# =============================================================================
# 6. Ajuste dos modelos
# =============================================================================

fit_id_lme4 <- lmer(
  y ~ trial + trial:block + (1 | gen) + (1 | gen:trial),
  data = dat, REML = TRUE
)

fit_cs_nlme <- lme(
  fixed = y ~ trial + trial:block,
  random = list(gen = pdCompSymm(~ trial - 1)),
  data = dat, method = "REML", na.action = na.omit,
  control = lmeControl(returnObject = TRUE)
)

fit_un_lme4 <- lmer(
  y ~ trial + trial:block + (0 + trial | gen),
  data = dat, REML = TRUE
)

fit_diag_sommer <- mmes(
  y ~ trial + trial:block,
  random = ~ vsm(dsm(trial), ism(gen)),
  rcov = ~ units,
  data = dat
)

fit_un_sommer <- mmes(
  y ~ trial + trial:block,
  random = ~ vsm(usm(trial), ism(gen)),
  rcov = ~ units,
  data = dat
)

# =============================================================================
# 7. Comparacao de ajuste
# =============================================================================

tab_fit_main <- rbind(
  data.frame(modelo = "ID_lme4", npar_vcov =  2, logLik = as.numeric(logLik(fit_id_lme4)),  AIC = AIC(fit_id_lme4),  BIC = BIC(fit_id_lme4)),
  data.frame(modelo = "CS_nlme", npar_vcov =  2, logLik = as.numeric(logLik(fit_cs_nlme)),  AIC = AIC(fit_cs_nlme),  BIC = BIC(fit_cs_nlme)),
  data.frame(modelo = "UN_lme4", npar_vcov = 21, logLik = as.numeric(logLik(fit_un_lme4)),  AIC = AIC(fit_un_lme4),  BIC = BIC(fit_un_lme4))
)
tab_fit_main[, -1] <- round(tab_fit_main[, -1], 2)
print(tab_fit_main, row.names = FALSE)

tab_fit_sommer <- rbind(
  data.frame(modelo = "DIAG_sommer", npar_vcov =  6, t(vc_sommer_fit(fit_diag_sommer))),
  data.frame(modelo = "UN_sommer",   npar_vcov = 21, t(vc_sommer_fit(fit_un_sommer)))
)
tab_fit_sommer[, -1] <- round(tab_fit_sommer[, -1], 2)
print(tab_fit_sommer, row.names = FALSE)

cat("\nNota: criterios de ajuste do sommer foram mantidos em tabela separada.\n")
cat("A escala reportada por mmes nao deve ser comparada diretamente com lme4/nlme.\n")

# =============================================================================
# 8. Matrizes G e correlacoes geneticas
# =============================================================================

vc_id <- as.data.frame(VarCorr(fit_id_lme4))
sigma2_g  <- vc_id$vcov[vc_id$grp == "gen"]
sigma2_ge <- vc_id$vcov[vc_id$grp == "gen:trial"]

G_id <- matrix(sigma2_g, J, J, dimnames = list(tnames, tnames))
diag(G_id) <- sigma2_g + sigma2_ge

G_cs <- as.matrix(pdMatrix(fit_cs_nlme$modelStruct$reStruct$gen)) * sigma(fit_cs_nlme)^2
rownames(G_cs) <- sub("^trial", "", rownames(G_cs))
colnames(G_cs) <- sub("^trial", "", colnames(G_cs))
G_cs <- G_cs[tnames, tnames, drop = FALSE]

G_un <- as.matrix(VarCorr(fit_un_lme4)$gen)
if (!all(dim(G_un) == c(J, J))) stop("G_un nao tem dimensao J x J.")
dimnames(G_un) <- list(tnames, tnames)

G_diag <- get_G_sommer(fit_diag_sommer, tnames)
G_uns  <- get_G_sommer(fit_un_sommer,  tnames)

R_id   <- cov2cor(G_id)
R_cs   <- cov2cor(G_cs)
R_un   <- cov2cor(G_un)
R_diag <- cov2cor(G_diag)
R_uns  <- cov2cor(G_uns)

cat("\nCorrelacao genetica - ID (lme4)\n");     print(round(R_id,   4))
cat("\nCorrelacao genetica - CS (nlme)\n");     print(round(R_cs,   4))
cat("\nCorrelacao genetica - UN (lme4)\n");     print(round(R_un,   4))
cat("\nCorrelacao genetica - DIAG (sommer)\n"); print(round(R_diag, 4))
cat("\nCorrelacao genetica - UN (sommer)\n");   print(round(R_uns,  4))

# =============================================================================
# 9. Heatmaps facetados das matrizes VCOV
# =============================================================================

Glist <- list(
  ID_lme4     = G_id,
  CS_nlme     = G_cs,
  UN_lme4     = G_un,
  DIAG_sommer = G_diag,
  UN_sommer   = G_uns
)

print(plot_vcov(Glist, tnames))

# =============================================================================
# 10. Comparacao direta entre UN lme4 e UN sommer
# =============================================================================

cmp_un <- data.frame(
  par = c("max_abs_diff_G", "max_abs_diff_R", "mean_abs_diff_R"),
  val = round(c(
    max(abs(G_un - G_uns), na.rm = TRUE),
    max(abs(R_un - R_uns), na.rm = TRUE),
    mean(abs(R_un - R_uns), na.rm = TRUE)
  ), 6)
)
print(cmp_un, row.names = FALSE)

# =============================================================================
# 11. BLUPs por ensaio, ranking medio e painel de comparacao
# =============================================================================

re_id_g   <- ranef(fit_id_lme4)$gen[, 1]
names(re_id_g) <- rownames(ranef(fit_id_lme4)$gen)

re_id_geM <- ranef(fit_id_lme4)$`gen:trial`
re_id_ge  <- re_id_geM[, 1]
rn_ge     <- rownames(re_id_geM)

lab_gt <- outer(gnames, tnames, paste, sep = ":")
lab_tg <- outer(gnames, tnames, function(g, t) paste(t, g, sep = ":"))

idx_gt <- match(c(lab_gt), rn_ge)
idx_tg <- match(c(lab_tg), rn_ge)

blup_id <- matrix(NA_real_, I, J, dimnames = list(gnames, tnames))
if (sum(!is.na(idx_gt)) >= sum(!is.na(idx_tg))) blup_id[] <- re_id_ge[idx_gt] else blup_id[] <- re_id_ge[idx_tg]
blup_id <- sweep(blup_id, 1, re_id_g[gnames], "+")

blup_cs <- as.matrix(ranef(fit_cs_nlme))
colnames(blup_cs) <- sub("^trial", "", colnames(blup_cs))
blup_cs <- blup_cs[gnames, tnames, drop = FALSE]

blup_un <- as.matrix(ranef(fit_un_lme4)$gen)
if (!all(dim(blup_un) == c(I, J))) stop("blup_un nao tem dimensao I x J.")
dimnames(blup_un) <- list(gnames, tnames)

blup_diag <- as.matrix(fit_diag_sommer$uList[[1]])
blup_diag <- blup_diag[gnames, tnames, drop = FALSE]

blup_uns <- as.matrix(fit_un_sommer$uList[[1]])
blup_uns <- blup_uns[gnames, tnames, drop = FALSE]

cat("\nNA counts in BLUP matrices\n")
print(c(
  ID_lme4     = sum(is.na(blup_id)),
  CS_nlme     = sum(is.na(blup_cs)),
  UN_lme4     = sum(is.na(blup_un)),
  DIAG_sommer = sum(is.na(blup_diag)),
  UN_sommer   = sum(is.na(blup_uns))
))

rank_id   <- sort(rowMeans(blup_id,   na.rm = TRUE), decreasing = TRUE)
rank_cs   <- sort(rowMeans(blup_cs,   na.rm = TRUE), decreasing = TRUE)
rank_un   <- sort(rowMeans(blup_un,   na.rm = TRUE), decreasing = TRUE)
rank_diag <- sort(rowMeans(blup_diag, na.rm = TRUE), decreasing = TRUE)
rank_uns  <- sort(rowMeans(blup_uns,  na.rm = TRUE), decreasing = TRUE)

cat("\nTop 10 genotipos - ID lme4\n");        print(round(head(rank_id,   10), 2))
cat("\nTop 10 genotipos - CS nlme\n");        print(round(head(rank_cs,   10), 2))
cat("\nTop 10 genotipos - UN lme4\n");        print(round(head(rank_un,   10), 2))
cat("\nTop 10 genotipos - DIAG sommer\n");    print(round(head(rank_diag, 10), 2))
cat("\nTop 10 genotipos - UN sommer\n");      print(round(head(rank_uns,  10), 2))

blup_long <- data.frame(
  gen         = rep(gnames, times = J),
  trial       = rep(tnames, each = I),
  ID_lme4     = c(blup_id),
  CS_nlme     = c(blup_cs),
  UN_lme4     = c(blup_un),
  DIAG_sommer = c(blup_diag),
  UN_sommer   = c(blup_uns)
)

blup_vars <- c("ID_lme4", "CS_nlme", "UN_lme4", "DIAG_sommer", "UN_sommer")
blup_mat  <- as.matrix(blup_long[, blup_vars])

trial_cols <- c(
  "TRL_CPS_2025" = "#D55E00",
  "TRL_DOU_2025" = "#0072B2",
  "TRL_JAT_2025" = "#009E73",
  "TRL_RVD_2025" = "#CC79A7",
  "TRL_SIN_2025" = "#E69F00",
  "TRL_SOR_2025" = "#56B4E9"
)
pt_col <- unname(trial_cols[blup_long$trial])

panel_scatter <- function(x, y, col_vec, ...) {
  usr <- par("usr"); on.exit(par(usr))
  points(x, y, pch = 16, cex = 0.8, col = col_vec)
  abline(0, 1, lty = 2, lwd = 1.2, col = "gray30")
}

panel_cor <- function(x, y, digits = 4, prefix = "r = ", cex.cor = 1.35, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use = "pairwise.complete.obs")
  text(0.5, 0.5, paste0(prefix, format(round(r, digits), nsmall = digits)), cex = cex.cor)
}

panel_hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  h <- hist(x, plot = FALSE, breaks = "Sturges")
  par(usr = c(usr[1:2], 0, max(h$counts) * 1.05))
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts, col = "gray80", border = "white")
}

op <- par(cex.axis = 1.15, cex.lab = 1.2, cex.main = 1.2, mar = c(4, 4, 3, 1))

pairs(
  blup_mat,
  labels = blup_vars,
  lower.panel = function(x, y, ...) panel_scatter(x, y, col_vec = pt_col, ...),
  upper.panel = panel_cor,
  diag.panel  = panel_hist,
  gap = 0.7,
  cex.labels = 1.2
)

legend(
  "topright", inset = 0.01, bty = "n", pch = 16, pt.cex = 1.1, cex = 1.05,
  col = unname(trial_cols),
  legend = sub("TRL_|_2025", "", names(trial_cols)),
  title = "Trial"
)

par(op)

cat("\nCorrelacoes entre modelos (todos os pontos gen x trial)\n")
print(round(cor(blup_mat, use = "pairwise.complete.obs"), 4))

# No modelo DIAG, a ausencia de covariancias geneticas entre ensaios faz com que
# varios BLUPs de genotipo dentro de ensaio encolham exatamente para zero; por
# isso, ao comparar esse modelo com estruturas como CS ou UN, surgem faixas de
# pontos alinhadas vertical ou horizontalmente nos paineis de dispersao.

# O modelo FA ficou de fora deste script.
