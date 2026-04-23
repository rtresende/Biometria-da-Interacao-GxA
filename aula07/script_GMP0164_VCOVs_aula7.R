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
#   GEN (lme4) :
#       y ~ trial + (1 | trial_block) + (1 | gen)
#       apenas um componente genetico comum a todos os ensaios
#       G_e = sigma^2_g J_J
#       correlacao genetica = 1 entre quaisquer ensaios (sem GxE)
#
#   ID (lme4) :
#       y ~ trial + (1 | trial_block) + (1 | gen:trial)
#       um componente independente por combinacao genotipo x ensaio
#       G_e = sigma^2_ge I_J
#       nenhuma covariancia entre ensaios
#
#   DG (lme4breeding) :
#       y ~ trial + (1 | trial_block) + (0 + trial || gen)
#       variancias especificas por ensaio
#       G_e = diag(sigma^2_g1,...,sigma^2_gJ)
#       covariancias iguais a zero
#
#   CS (lme4) :
#       y ~ trial + (1 | trial_block) + (1 | gen) + (1 | gen:trial)
#       decomposicao em efeito genetico principal + GxE
#       G_e = sigma^2_g J_J + sigma^2_ge I_J
#       variancia total constante e covariancia constante entre ensaios
#
#   CSH (lme4breeding) :
#       y ~ trial + (1 | trial_block) + (trial || gen)
#       combinacao de efeito genetico principal + variancias especificas
#       G_e = sigma^2_g J_J + diag(sigma^2_g1,...,sigma^2_gJ)
#       covariancia constante, variancias heterogeneas
#
#   US (lme4breeding) :
#       y ~ trial + (1 | trial_block) + (0 + trial | gen)
#       matriz G_e totalmente livre (J x J)
#       estima todas as variancias e covariancias entre ensaios
#
#   FA1 (lme4breeding) :
#       y ~ trial + (1 | trial_block) + (0 + PC1 | gen) + (0 + trial || gen)
#       modelo fatorial de posto 1
#       G_e = Lambda Lambda' + Psi
#       Lambda: loadings por ensaio; Psi: variancias especificas
#
#   FA2 (lme4breeding) :
#       y ~ trial + (1 | trial_block) + (0 + PC1 + PC2 | gen) + (0 + trial || gen)
#       modelo fatorial de posto 2
#       G_e = Lambda Lambda' + Psi
#       maior flexibilidade na representacao das covariancias
#
# Estrategia de analise:
#
#   1) Ajuste de todos os modelos sob REML.
#   2) Extracao de componentes de variancia via VarCorr().
#   3) Organizacao em tabela unica (variancias x modelos).
#   4) Reconstrucao explicita das matrizes G (6 x 6) para cada modelo.
#   5) Visualizacao via heatmap (ggplot2) com valores nas celulas.
#   6) Comparacao de ajuste via AIC e BIC (ranking de modelos).
#   7) Reconstrucao de pseudo-fenotipos:
#          y* = BLUE(trial) + BLUP(gen,trial)
#      com adaptacao conforme a estrutura de cada modelo.
#   8) Comparacao entre modelos via matriz de dispersao (GGally):
#       - lower: scatter geral com r global
#       - upper: correlacoes por ensaio
#       - diagonal: densidades por ensaio
#
# Observacoes:
#   1) Blocos dentro de ensaio entram como efeito aleatorio em todos os modelos.
#   2) Para FA, a matriz G_e nao e diretamente observada em VarCorr(),
#      sendo reconstruida como G = Lambda Lambda' + Psi.
#   3) Para CS e CSH, a estrutura e induzida pela combinacao de termos aleatorios.
#   4) A comparacao entre modelos considera tanto ajuste (AIC/BIC)
#      quanto comportamento das matrizes G e dos y*.
# =============================================================================


rm(list=ls());gc()

# pacotes
library(lme4)
library(lme4breeding)
library(Matrix)
library(ggplot2)

# dados
dat <- read.table(
  "G:/Meu Drive/UFG/PPGGMP/BiometriaDaInteracaoGxA/aulas/simu_data/maize_safrinha_MET.txt",
  header = TRUE, sep = "\t")

# preparo
names(dat)[names(dat) == "yield"] <- "y"
dat$gen <- factor(dat$gen)
dat$trial <- factor(dat$trial)
dat$block <- factor(dat$block)
dat$trial_block <- interaction(dat$trial, dat$block, drop = TRUE)

# matriz A identidade
A <- Diagonal(nlevels(dat$gen))
rownames(A) <- levels(dat$gen)
colnames(A) <- levels(dat$gen)

# modelos
m_gen <- lmer(
  y ~ trial + (1 | trial_block) + (1 | gen),
  data = dat, REML = TRUE)

m_id <- lmer(
  y ~ trial + (1 | trial_block) + (1 | gen:trial),
  data = dat, REML = TRUE)

m_dg <- lmeb(
  y ~ trial + (1 | trial_block) + (0 + trial || gen),
  relmat = list(gen = A),
  verbose = 0L, trace = 0L,
  data = dat, REML = TRUE)

m_cs <- lmer(
  y ~ trial + (1 | trial_block) + (1 | gen) + (1 | gen:trial),
  data = dat, REML = TRUE)

m_csh <- lmeb(
  y ~ trial + (1 | trial_block) + (trial || gen),
  relmat = list(gen = A),
  verbose = 0L, trace = 0L,
  data = dat, REML = TRUE)

m_us <- lmeb(
  y ~ trial + (1 | trial_block) + (0 + trial | gen),
  relmat = list(gen = A),
  verbose = 0L, trace = 0L,
  data = dat, REML = TRUE)

# FA
Z0 <- with(dat, smm(trial))
if (is.null(colnames(Z0))) colnames(Z0) <- levels(dat$trial)
for (j in seq_len(ncol(Z0))) dat[, colnames(Z0)[j]] <- Z0[, j]

m_fa_seed <- lmeb(
  as.formula(paste0(
    "y ~ trial + (1 | trial_block) + (0 + ",
    paste(colnames(Z0), collapse = " + "),
    " || gen)"
  )),
  relmat = list(gen = A),
  verbose = 0L, trace = 0L,
  data = dat, REML = TRUE)

H0 <- ranef(m_fa_seed)$gen

Z1 <- with(dat, rrm(trial, H = H0, nPC = 1))
if (is.null(dim(Z1))) Z1 <- matrix(Z1, ncol = 1)
colnames(Z1) <- "PC1"
dat$PC1 <- Z1[, 1]

m_fa1 <- lmeb(
  y ~ trial + (1 | trial_block) + (0 + PC1 | gen) + (0 + trial || gen),
  relmat = list(gen = A),
  verbose = 0L, trace = 0L,
  data = dat, REML = TRUE)

fa_loadings1 <- with(dat, rrm(trial, H = H0, nPC = 1, returnGamma = TRUE))$Gamma

Z2 <- with(dat, rrm(trial, H = H0, nPC = 2))
if (is.null(dim(Z2))) Z2 <- matrix(Z2, ncol = 2)
colnames(Z2) <- c("PC1", "PC2")
dat$PC1 <- Z2[, 1]
dat$PC2 <- Z2[, 2]

m_fa2 <- lmeb(
  y ~ trial + (1 | trial_block) + (0 + PC1 + PC2 | gen) + (0 + trial || gen),
  relmat = list(gen = A),
  verbose = 0L, trace = 0L,
  data = dat, REML = TRUE)

fa_loadings2 <- with(dat, rrm(trial, H = H0, nPC = 2, returnGamma = TRUE))$Gamma

# modelos em ordem
mods <- list(GEN = m_gen, ID = m_id, DG = m_dg, CS = m_cs,
             CSH = m_csh, US = m_us, FA1 = m_fa1, FA2 = m_fa2)


# ===================
# AIC/BIC e Variancias
# ===================

fit_stats <- data.frame(
  model = names(mods),
  AIC = sapply(mods, AIC),
  BIC = sapply(mods, BIC),
  row.names = NULL)

print(fit_stats, row.names = FALSE)

get_lab <- function(x) ifelse(
  x$grp == "trial_block", "Vbloco",
  ifelse(x$grp == "Residual", "Ve residual",
  ifelse(x$grp == "gen:trial", "Vgxe geral",
  ifelse(x$var1 == "(Intercept)" & grepl("^gen", x$grp), "Vg principal",
  ifelse(x$var2 %in% c(NA, "<NA>") & grepl("^TRL_", x$var1), paste0("Vg ", x$var1),
  ifelse(x$var1 %in% c("PC1", "PC2") & x$var2 %in% c(NA, "<NA>"), paste0("Vfator ", x$var1),
  NA))))))

tabs <- lapply(mods, function(m) {
  x <- as.data.frame(VarCorr(m))
  x$lab <- get_lab(x)
  tapply(x$vcov, x$lab, sum)})

ord <- c(
  "Vg principal", "Vgxe geral",
  paste0("Vg ", levels(dat$trial)),
  "Vfator PC1", "Vfator PC2",
  "Vbloco", "Ve residual", "AIC", "BIC")

var_tab <- sapply(names(mods), function(nm) {
  z <- tabs[[nm]]
  out <- rep(NA_real_, length(ord))
  names(out) <- ord
  out[names(z)] <- z
  out["AIC"] <- AIC(mods[[nm]])
  out["BIC"] <- BIC(mods[[nm]])
  out})

var_tab <- data.frame(Variancia = ord, var_tab, row.names = NULL)
print(var_tab)

# AIC/BIC por modelo
fit_stats <- data.frame(
  model = names(mods),
  AIC = sapply(mods, AIC),
  BIC = sapply(mods, BIC),
  row.names = NULL)

fit_stats <- fit_stats[order(fit_stats$AIC, fit_stats$BIC), ]
print(fit_stats, row.names = FALSE) # ordenado do melhor modelo pro pior

# ===================
# Variancias e Covariancias
# ===================

# ordem e rotulos dos ambientes
env <- levels(dat$trial)
env_lab <- sub("TRL_", "", sub("_2025", "", env))

# helper
vc <- function(m) as.data.frame(VarCorr(m))

# GEN
x <- vc(m_gen)
vg <- x$vcov[x$grp == "gen"]
G_gen <- vg * matrix(1, length(env), length(env), dimnames = list(env, env))

# ID
x <- vc(m_id)
vge <- x$vcov[x$grp == "gen:trial"]
G_id <- vge * diag(length(env))
dimnames(G_id) <- list(env, env)

# DG
x <- vc(m_dg)
v <- x$vcov[grepl("^TRL_", x$var1)]
nm <- x$var1[grepl("^TRL_", x$var1)]
G_dg <- diag(v[match(env, nm)])
dimnames(G_dg) <- list(env, env)

# CS
x <- vc(m_cs)
vg <- x$vcov[x$grp == "gen"]
vge <- x$vcov[x$grp == "gen:trial"]
G_cs <- vg * matrix(1, length(env), length(env)) + vge * diag(length(env))
dimnames(G_cs) <- list(env, env)

# CSH
x <- vc(m_csh)
vg <- x$vcov[x$var1 == "(Intercept)" & grepl("^gen", x$grp) & !is.na(x$vcov)]
v <- x$vcov[grepl("^TRL_", x$var1)]
nm <- x$var1[grepl("^TRL_", x$var1)]
G_csh <- vg * matrix(1, length(env), length(env)) + diag(v[match(env, nm)])
dimnames(G_csh) <- list(env, env)

# US
x <- vc(m_us)
G_us <- matrix(0, length(env), length(env), dimnames = list(env, env))
for (i in 1:nrow(x)) {
  if (x$grp[i] == "gen" && grepl("^TRL_", x$var1[i])) {
    if (is.na(x$var2[i]) || x$var2[i] == "<NA>") {
      G_us[x$var1[i], x$var1[i]] <- x$vcov[i]
    } else {
      G_us[x$var1[i], x$var2[i]] <- x$vcov[i]
      G_us[x$var2[i], x$var1[i]] <- x$vcov[i]
    }
  }
}

# FA1
x <- vc(m_fa1)
Lambda1 <- as.matrix(fa_loadings1[env, , drop = FALSE])
Phi1 <- matrix(
  x$vcov[x$grp == "gen.6" & x$var1 == "PC1" & is.na(x$var2)],
  1, 1, dimnames = list("PC1", "PC1"))
psi <- x$vcov[grepl("^TRL_", x$var1)]
nm <- x$var1[grepl("^TRL_", x$var1)]
Psi1 <- diag(psi[match(env, nm)])
dimnames(Psi1) <- list(env, env)
G_fa1 <- Lambda1 %*% Phi1 %*% t(Lambda1) + Psi1

# FA2
x <- vc(m_fa2)
Lambda2 <- as.matrix(fa_loadings2[env, , drop = FALSE])
Phi2 <- matrix(0, 2, 2, dimnames = list(c("PC1", "PC2"), c("PC1", "PC2")))
Phi2["PC1", "PC1"] <- x$vcov[x$grp == "gen.6" & x$var1 == "PC1" & is.na(x$var2)]
Phi2["PC2", "PC2"] <- x$vcov[x$grp == "gen.6" & x$var1 == "PC2" & is.na(x$var2)]
Phi2["PC1", "PC2"] <- x$vcov[x$grp == "gen.6" & x$var1 == "PC1" & !is.na(x$var2) & x$var2 == "PC2"]
Phi2["PC2", "PC1"] <- Phi2["PC1", "PC2"]
psi <- x$vcov[grepl("^TRL_", x$var1)]
nm <- x$var1[grepl("^TRL_", x$var1)]
Psi2 <- diag(psi[match(env, nm)])
dimnames(Psi2) <- list(env, env)
G_fa2 <- Lambda2 %*% Phi2 %*% t(Lambda2) + Psi2

# lista
G_list <- list(
  GEN = G_gen,
  ID = G_id,
  DG = G_dg,
  CS = G_cs,
  CSH = G_csh,
  US = G_us,
  FA1 = G_fa1,
  FA2 = G_fa2)

# formato longo
plot_dat <- do.call(rbind, lapply(names(G_list), function(nm) {
  g <- G_list[[nm]]
  d <- expand.grid(row = rownames(g), col = colnames(g), stringsAsFactors = FALSE)
  d$value <- as.vector(g)
  d$model <- nm
  d
}))

plot_dat$model <- factor(plot_dat$model, levels = names(G_list))
plot_dat$row <- factor(plot_dat$row, levels = rev(env), labels = rev(env_lab))
plot_dat$col <- factor(plot_dat$col, levels = env, labels = env_lab)
plot_dat$lab <- as.character(round(plot_dat$value, 0))

ggplot(plot_dat, aes(col, row, fill = value)) +
  geom_tile(color = "grey85") +
  geom_text(aes(label = lab), size = 2.7) +
  facet_wrap(~ model, nrow = 2) +
  coord_equal() +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0) +
  theme_bw() +
  labs(x = NULL, y = NULL, fill = "Value")


# ===================
# BLUPs / y_star
# ===================

env <- levels(dat$trial)

blue <- function(m) {
  b <- fixef(m)
  out <- setNames(rep(b["(Intercept)"], length(env)), env)
  out[names(b)[-1]] <- out[names(b)[-1]] + b[-1]
  names(out) <- sub("^trial", "", names(out))
  out[env]
}

make_grid <- function() expand.grid(
  gen = levels(dat$gen),
  trial = env,
  stringsAsFactors = FALSE)

get_gt <- function(x) {
  sp <- strsplit(rownames(x), ":", fixed = TRUE)
  data.frame(
    gen = sapply(sp, "[", 1),
    trial = sapply(sp, "[", 2),
    value = x[, 1],
    stringsAsFactors = FALSE)
}

make_blup <- function(model, m) {
  rr <- lme4::ranef(m)
  out <- make_grid()
  out$model <- model
  out$blue_trial <- blue(m)[out$trial]
  out$blup <- 0

  if (model == "GEN") {out$blup <- rr$gen[out$gen, 1]}

  if (model == "ID") {z <- get_gt(rr$`gen:trial`)
      out$blup <- z$value[match(paste(out$gen, out$trial), paste(z$gen, z$trial))]}

  if (model == "DG") {out$blup <- rr$gen[cbind(out$gen, out$trial)]}

  if (model == "CS") {z <- get_gt(rr$`gen:trial`)
      out$blup <- rr$gen[out$gen, 1] +
      z$value[match(paste(out$gen, out$trial), paste(z$gen, z$trial))]}

  if (model == "CSH") {out$blup <- rr$gen[out$gen, "(Intercept)"] +
      rr$gen[cbind(out$gen, out$trial)]}

  if (model == "US") {out$blup <- rr$gen[cbind(out$gen, out$trial)]}

  if (model == "FA1") {out$blup <- rr$gen[cbind(out$gen, out$trial)] +
      rr$gen[out$gen, "PC1"] * fa_loadings1[out$trial, "PC1"]}

  if (model == "FA2") {out$blup <- rr$gen[cbind(out$gen, out$trial)] +
      rr$gen[out$gen, "PC1"] * fa_loadings2[out$trial, "PC1"] +
      rr$gen[out$gen, "PC2"] * fa_loadings2[out$trial, "PC2"]}

  out$y_star <- out$blue_trial + out$blup
  out[, c("model", "gen", "trial", "blue_trial", "blup", "y_star")]
}

blup_tab <- do.call(rbind, lapply(names(mods), function(nm) make_blup(nm, mods[[nm]])))

# ===================
# ggpairs para y_star
# ===================

library(GGally)

wide <- reshape(
  blup_tab[, c("model", "gen", "trial", "y_star")],
  idvar = c("gen", "trial"),
  timevar = "model",
  direction = "wide")

names(wide) <- sub("^y_star\\.", "", names(wide))
wide$Trial <- sub("TRL_", "", sub("_2025", "", wide$trial))

mod_names <- names(mods)
wide <- wide[complete.cases(wide[, c(mod_names, "Trial")]), ]

lower_fun <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  r <- cor(x, y, use = "complete.obs")
  
  ggplot(data, mapping) +
    geom_point(aes(color = Trial), size = 1.5, alpha = 0.70) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
    annotate("text", x = -Inf, y = Inf, label = paste0("r = ", round(r, 3)),
             hjust = -0.05, vjust = 1.2, size = 3.2) +
    theme_bw()
}

upper_fun <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  tr <- data$Trial
  rr <- tapply(seq_along(x), tr, function(i) cor(x[i], y[i], use = "complete.obs"))
  lab <- paste(names(rr), sprintf("%.3f", rr), sep = ": ", collapse = "\n")
  
  ggplot(data, mapping) +
    annotate("text", x = 0.5, y = 0.5, label = lab, size = 2.9) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank())
}

diag_fun <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_density(aes(color = Trial, fill = Trial), alpha = 0.25, linewidth = 0.7) +
    theme_bw()
}

p_pairs <- ggpairs(
  wide,
  columns = mod_names,
  mapping = aes(color = Trial),
  lower = list(continuous = lower_fun),
  upper = list(continuous = upper_fun),
  diag = list(continuous = diag_fun),
  progress = FALSE)

p_pairs
