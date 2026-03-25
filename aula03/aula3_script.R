# =============================================================================
# GMP0164 - Biometria da Interacao Genotipos x Ambientes
# Aula 3: Regressao conjunta para adaptabilidade e estabilidade
# Dados: soy_MET.txt | 40 genotipos | 6 ambientes | 3 blocos | DBC | 720 obs
#        CHK01-CHK04 = testemunhas | L001-L036 = linhagens candidatas
#        Variavel resposta: produtividade (kg/ha)
# Prof. Dr. Rafael Tassinari Resende
# =============================================================================
rm(list = ls()); gc()
library(ggplot2)

# -----------------------------------------------------------------------------
# 0. DATA
# -----------------------------------------------------------------------------

dat <- read.table("soy_MET.txt", header = TRUE, sep = "\t")
dat$gen   <- factor(dat$gen)
dat$env   <- factor(dat$env)
dat$block <- factor(dat$block)

env_info <- data.frame(
  env      = levels(dat$env),
  location = c("Anapolis", "Luziania", "Jatai", "RioVerde", "Uberlandia", "Sorriso"),
  state    = c("GO", "GO", "GO", "GO", "MG", "MT"),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 1. DIMENSOES AND DESCRITORES
# -----------------------------------------------------------------------------

g <- nlevels(dat$gen)
a <- nlevels(dat$env)
r <- nlevels(dat$block)

gen_chk <- grep("^CHK", levels(dat$gen), value = TRUE)
gen_lin <- grep("^L",   levels(dat$gen), value = TRUE)

means_ge   <- with(dat, tapply(y, list(gen, env), mean))
env_means  <- colMeans(means_ge)
gen_means  <- rowMeans(means_ge)
grand_mean <- mean(means_ge)

cat(sprintf("\n=== DIMENSOES ===\n  Genotipos: %d  (%d CHK testemunhas + %d L linhagens candidatas)\n",
            g, length(gen_chk), length(gen_lin)))
cat(sprintf("  Ambientes: %d | Blocos: %d | Obs: %d\n", a, r, nrow(dat)))

cat("\nMedias por ambiente [kg/ha]:\n")
print(data.frame(
  Amb    = names(env_means),
  Local  = env_info$location[match(names(env_means), env_info$env)],
  Estado = env_info$state[match(names(env_means), env_info$env)],
  Media  = round(env_means, 1)
), row.names = FALSE)

cat("\nMedias por genotipo (top 10):\n")
print(round(sort(gen_means, decreasing = TRUE)[1:10], 1))

cat(sprintf("\nMedia geral: %.1f kg/ha | Amplitude: %.1f a %.1f\n",
            grand_mean, min(env_means), max(env_means)))

# Mean table in long format for plotting
mean_df <- as.data.frame(as.table(means_ge))
names(mean_df) <- c("gen", "env", "y_mean")
mean_df$gen   <- factor(mean_df$gen, levels = levels(dat$gen))
mean_df$env   <- factor(mean_df$env, levels = levels(dat$env))
mean_df$Tipo  <- ifelse(mean_df$gen %in% gen_chk, "CHK", "L")
mean_df$env_n <- match(mean_df$env, levels(dat$env))
mean_df$env_lab <- paste0(mean_df$env, "\n(", env_info$location[match(mean_df$env, env_info$env)], ")")

p_gxe <- ggplot(mean_df, aes(x = env_n, y = y_mean, group = gen, color = Tipo)) +
  geom_hline(yintercept = grand_mean, linetype = 2, linewidth = 0.7) +
  geom_line(linewidth = 0.55, alpha = 0.55) +
  scale_color_manual(values = c("CHK" = "tomato", "L" = "gray55")) +
  scale_x_continuous(breaks = 1:a, labels = paste0(env_info$env, "\n(", env_info$location, ")")) +
  labs(x = "Ambiente", y = "Produtividade media (kg/ha)",
       title = "Grafico de interacao G x A",
       subtitle = "Vermelho = testemunhas (CHK) | Cinza = linhagens (L)") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(size = 9),
        plot.title = element_text(face = "bold"))
print(p_gxe)

# =============================================================================
# SECAO A: INDICE AMBIENTAL
# I_j = ybar_.j - ybar_..
# Propriedades: sum(I_j) = 0; I_j < 0 = desfavoravel; I_j > 0 = favoravel
# Construido a partir das proprias respostas fenotipicas (endogeno)
# =============================================================================

cat("\n\n=== SECAO A: INDICE AMBIENTAL ===\n")
cat("I_j = media do ambiente j - media geral\n")
cat("Propriedade: sum(I_j) = 0  (centrado; escala empirica de favorabilidade)\n\n")

I_j <- env_means - grand_mean

idx_tab <- data.frame(
  Ambiente = names(I_j),
  Local    = env_info$location[match(names(I_j), env_info$env)],
  Media    = round(env_means, 1),
  I_j      = round(I_j, 1),
  Tipo     = ifelse(I_j >= 0, "Favoravel (I>0)", "Desfavoravel (I<0)")
)
print(idx_tab, row.names = FALSE)

cat(sprintf("\nVerificacao: sum(I_j) = %.6f  (deve ser 0)\n", sum(I_j)))
cat("\nObs.: I_j e construido a partir das PROPRIAS medias do ensaio (indice endogeno).\n")
cat("      Ele nao mede clima/solo diretamente; resume o efeito agregado do ambiente\n")
cat("      sobre todos os genotipos avaliados.\n")

idx_df <- data.frame(
  env      = names(I_j),
  location = env_info$location[match(names(I_j), env_info$env)],
  I_j      = as.numeric(I_j),
  stringsAsFactors = FALSE
)
idx_df$Tipo <- ifelse(idx_df$I_j >= 0, "Favoravel", "Desfavoravel")
idx_df$lab  <- paste0(idx_df$env, "\n(", idx_df$location, ")")
idx_df      <- idx_df[order(idx_df$I_j), ]

p_idx <- ggplot(idx_df, aes(x = reorder(lab, I_j), y = I_j, fill = Tipo)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.7) +
  geom_text(aes(label = round(I_j, 0)),
            vjust = ifelse(idx_df$I_j >= 0, -0.35, 1.15), size = 3.4) +
  scale_fill_manual(values = c("Favoravel" = "steelblue", "Desfavoravel" = "tomato")) +
  labs(x = NULL, y = "I_j (kg/ha)", title = "Indice ambiental por ambiente") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(face = "bold"))
print(p_idx)

cat("\n\n=== SECAO B: FINLAY E WILKINSON (1963) ===\n")
cat("Modelo: y_ij = mu_i + b_i*I_j + e_ij\n")
cat("Ajuste via lm(y ~ gen * I_j), no nivel de parcela.\n")
cat("Como sum(I_j)=0: o intercepto de cada genotipo coincide com sua media.\n\n")

dat$I_j <- I_j[as.character(dat$env)] #atribui o indice ambiental a cada parcela

fit_fw <- lm(y ~ gen * I_j, data = dat) #rodando o Finlay & Wilkinson (1963)
cf <- coef(fit_fw) #guardando os coeficientes

get_cf <- function(nm) if (nm %in% names(cf)) cf[nm] else 0 #funcao auxiliar para evitar NA

gen_levels <- levels(dat$gen) #lista de genotipos
base_gen   <- gen_levels[1]   #genotipo de referencia do lm()

mu_hat <- setNames(numeric(length(gen_levels)), gen_levels) #vetor de medias (interceptos)
b_hat  <- setNames(numeric(length(gen_levels)), gen_levels) #vetor de inclinacoes (b_i)

mu_hat[base_gen] <- get_cf("(Intercept)") #intercepto do genotipo base
b_hat[base_gen]  <- get_cf("I_j")         #inclinação do genotipo base

for (nm in gen_levels[-1]) {
  mu_hat[nm] <- get_cf("(Intercept)") + get_cf(paste0("gen", nm)) #ajuste do intercepto
  b_hat[nm]  <- get_cf("I_j") +                                  #inclinação base
                get_cf(paste0("gen", nm, ":I_j")) +              #efeito interação
                get_cf(paste0("I_j:gen", nm))                    #garante robustez ao nome
}

cat(sprintf("Verificacao: max|mu_i - media do genotipo| = %.10f\n\n",
            max(abs(mu_hat - gen_means[names(mu_hat)])))) #checa equivalencia com media

perfil_fw <- ifelse(b_hat > 1.10, "Alta responsividade (b>1)", #classificacao operacional
             ifelse(b_hat < 0.90, "Baixa responsividade (b<1)",
                                  "Adapt. geral (b~1)"))

tipo_gen <- ifelse(names(mu_hat) %in% gen_chk, "CHK", "L") #identifica CHK vs linhagem

fw_tab <- data.frame(
  Tipo     = tipo_gen,
  Genotipo = names(mu_hat),
  Media    = round(mu_hat, 1),
  b_i      = round(b_hat, 3),
  Perfil   = perfil_fw,
  stringsAsFactors = FALSE
)

fw_tab <- fw_tab[order(-fw_tab$Media), ] #ordena por media decrescente
rownames(fw_tab) <- NULL

cat("Tabela F&W (por media desc. | CHK=testemunha | L=linhagem candidata):\n")
print(fw_tab, row.names = FALSE)
cat("\nLimites operacionais: b>1.10 = alta resp. | b<0.90 = baixa resp. | outros = geral\n")
cat("A inferencia formal (H0: b_i=1) sera feita em Eberhart e Russell.\n")

# Data for plots
point_df <- mean_df #base de pontos observados (medias)
point_df$I_j <- I_j[as.character(point_df$env)] #adiciona I_j
point_df$Perfil <- ifelse(b_hat[as.character(point_df$gen)] > 1.10, "b>1",
                   ifelse(b_hat[as.character(point_df$gen)] < 0.90, "b<1", "b~1")) #perfil por ponto
point_df$Tipo <- ifelse(point_df$gen %in% gen_chk, "CHK", "L") #tipo grafico

I_seq <- seq(min(I_j) * 1.15, max(I_j) * 1.15, length.out = 200) #grade de I_j para retas

line_df <- expand.grid(
  gen = gen_levels,
  I_j = I_seq,
  stringsAsFactors = FALSE
)
line_df$gen <- factor(line_df$gen, levels = levels(dat$gen)) #garante mesma ordem
line_df$pred <- predict(fit_fw, newdata = line_df) #valores ajustados do lm()
line_df$Perfil <- ifelse(b_hat[as.character(line_df$gen)] > 1.10, "b>1",
                  ifelse(b_hat[as.character(line_df$gen)] < 0.90, "b<1", "b~1")) #perfil das retas

top3 <- names(sort(b_hat, decreasing = TRUE))[1:3] #maiores b_i
bot3 <- names(sort(b_hat))[1:3]                   #menores b_i

highlight_df <- line_df[line_df$gen %in% c(top3, bot3), ] #subset destaque
highlight_df$Grupo <- ifelse(as.character(highlight_df$gen) %in% top3, "Top b", "Bottom b")

label_x <- max(I_j) * 0.88 #posicao dos rótulos
label_df <- data.frame(
  gen   = c(top3, bot3),
  x     = label_x,
  y     = mu_hat[c(top3, bot3)] + b_hat[c(top3, bot3)] * label_x, #equacao da reta
  label = paste0(c(top3, bot3), "\n(b=", round(b_hat[c(top3, bot3)], 2), ")"),
  Grupo = c(rep("Top b", 3), rep("Bottom b", 3)),
  stringsAsFactors = FALSE
)

p_fw <- ggplot() +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.6, color = "gray40") + #eixo I_j=0
  geom_hline(yintercept = grand_mean, linetype = 2, linewidth = 0.6, color = "gray40") + #media geral
  geom_line(data = line_df,
            aes(x = I_j, y = pred, group = gen, color = Perfil),
            linewidth = 0.45, alpha = 0.35) + #todas as retas
  geom_point(data = point_df,
             aes(x = I_j, y = y_mean, color = Perfil, shape = Tipo),
             size = 1.8, alpha = 0.70) + #pontos observados
  geom_line(data = highlight_df,
            aes(x = I_j, y = pred, group = gen, color = Grupo),
            linewidth = 1.1, show.legend = FALSE) + #destaque extremos
  geom_text(data = label_df,
            aes(x = x, y = y, label = label, color = Grupo),
            hjust = 0, size = 3.1, show.legend = FALSE) + #rotulos
  scale_color_manual(values = c(
    "b>1"      = "tomato3",
    "b<1"      = "steelblue",
    "b~1"      = "gray60",
    "Top b"    = "darkred",
    "Bottom b" = "darkblue"
  )) +
  scale_shape_manual(values = c("CHK" = 17, "L" = 16)) + #triangulo vs ponto
  labs(x = "I_j (kg/ha)", y = "Produtividade media (kg/ha)",
       title = "F&W (1963) | Retas de regressao por genotipo",
       subtitle = "Triangulo = CHK | Ponto = linhagem | Cor = perfil de b_i") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_fw)

quad_df <- data.frame(
  Genotipo = names(mu_hat),
  Media    = as.numeric(mu_hat),
  b_i      = as.numeric(b_hat),
  Tipo     = ifelse(names(mu_hat) %in% gen_chk, "CHK", "L"),
  Perfil   = ifelse(b_hat > 1.10, "b>1",
             ifelse(b_hat < 0.90, "b<1", "b~1")),
  stringsAsFactors = FALSE
)

x_right <- max(quad_df$Media) - 0.01 * diff(range(quad_df$Media)) #limite direito
y_top   <- max(quad_df$b_i)   - 0.01 * diff(range(quad_df$b_i))    #topo
y_bot   <- min(quad_df$b_i)   + 0.01 * diff(range(quad_df$b_i))    #base

p_quad <- ggplot(quad_df, aes(x = Media, y = b_i, color = Perfil, shape = Tipo)) +
  geom_vline(xintercept = grand_mean, linetype = 2, linewidth = 0.6, color = "gray40") + #media geral
  geom_hline(yintercept = 1,          linetype = 2, linewidth = 0.6, color = "gray40") + #b=1
  geom_point(size = 2.2) +
  geom_text(aes(label = Genotipo), vjust = -0.6, size = 2.7) + #rotulos
  scale_color_manual(values = c("b>1" = "tomato3", "b<1" = "steelblue", "b~1" = "gray60")) +
  scale_shape_manual(values = c("CHK" = 17, "L" = 16)) +
  annotate("text", x = x_right, y = y_top, label = "Alta media\nAlta resp.",
           hjust = 1, vjust = 1, color = "gray50", size = 3.4) +
  annotate("text", x = x_right, y = y_bot, label = "Alta media\nBaixa resp.",
           hjust = 1, vjust = 0, color = "gray50", size = 3.4) +
  labs(x = "Media do genotipo (kg/ha)", y = "b_i",
       title = "F&W: Media x Adaptabilidade",
       subtitle = "Triangulo = CHK | Ponto = L") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_quad)


# =============================================================================
# SECAO C: EBERHART E RUSSELL (1966)
# Extensao de F&W: adiciona sigma2_di como medida de estabilidade
# Modelo (medias): Y_ij = beta0_i + beta1_i*I_j + delta_ij
# Modelo (parcelas): Y_ijk = beta0_i + beta1_i*I_j + delta_ij + epsilon_ijk
# TRES PARAMETROS:
#   beta0_i  = media do genotipo      (nivel de producao)
#   beta1_i  = coef. de regressao     (adaptabilidade; identico a b_i de F&W)
#   sigma2_di= var. dos desvios       (estabilidade; ~0 = previsivel)
# GENOTIPO IDEAL: alta media + beta1~1 + sigma2_di nao significativo
# =============================================================================

cat("\n\n=== SECAO C: EBERHART E RUSSELL (1966) ===\n")
cat("Extensao de F&W: inclui sigma2_di como parametro de estabilidade.\n")
cat("Tripe: beta0_i (nivel) | beta1_i (adaptabilidade) | sigma2_di (previsibilidade)\n\n")

fmt_p <- function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p)) #formata p-valor
sig_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(p < 0.10, ".", " ")))) #marcas de significancia

# -----------------------------------------------------------------------------
# C1. ANOVA conjunta para obter QMR
# -----------------------------------------------------------------------------

cat("--- C1. ANOVA conjunta (para obter QMR) ---\n")

fit_joint <- aov(y ~ gen * env + env:block, data = dat) #anova conjunta
tab_joint <- summary(fit_joint)[[1]] #tabela da anova

cat("\nANOVA conjunta:\n")
print(round(tab_joint, 3))

QMR <- tab_joint["Residuals", "Mean Sq"] #quadrado medio do residuo
nu  <- tab_joint["Residuals", "Df"] #gl do residuo

sum_Ij2  <- sum(I_j^2) #soma de quadrados do indice ambiental
sigma2_e <- QMR / r #variancia residual entre medias
V_beta1  <- sigma2_e / sum_Ij2 #variancia de beta1

cat(sprintf("\nQMR=%.2f | nu=%d | sigma2_e=QMR/r=%.2f | V(beta1)=%.6f\n",
            QMR, nu, sigma2_e, V_beta1))
cat(sprintf("GL esperado: a*(r-1)*(g-1) = %d*%d*%d = %d\n",
            a, r - 1, g - 1, a * (r - 1) * (g - 1)))

# -----------------------------------------------------------------------------
# C2. Estimacao dos parametros
# beta1 vem da SECAO B via lm()
# sigma2_di e R2 seguem a forma classica de E&R
# -----------------------------------------------------------------------------

cat("\n--- C2. Estimacao dos parametros ---\n")

gen_id <- rownames(means_ge) #ordem base dos genotipos
mu_vec <- as.numeric(mu_hat[gen_id]) #medias por genotipo
b_vec  <- as.numeric(b_hat[gen_id]) #inclinacoes por genotipo

SQ_AG  <- r * rowSums((means_ge - mu_vec)^2) #SQ de ambientes dentro de genotipo
SQ_reg <- r * (as.vector(means_ge %*% I_j)^2) / sum_Ij2 #SQ linear da regressao
QMD    <- (SQ_AG - SQ_reg) / (a - 2) #quadrado medio dos desvios

er <- data.frame(
  Tipo      = ifelse(gen_id %in% gen_chk, "CHK", "L"),
  Genotipo  = gen_id,
  Media     = round(mu_vec, 1),
  beta1     = round(b_vec, 3),
  sigma2_di = round((QMD - QMR) / r, 1),
  R2        = round(SQ_reg / SQ_AG * 100, 1),
  t_stat    = round((b_vec - 1) / sqrt(V_beta1), 3),
  p_t       = 2 * pt(-abs((b_vec - 1) / sqrt(V_beta1)), df = nu),
  F_stat    = round(QMD / QMR, 3),
  p_F       = pf(QMD / QMR, df1 = a - 2, df2 = nu, lower.tail = FALSE),
  stringsAsFactors = FALSE
)

er <- er[order(-er$Media), ] #ordena por media
rownames(er) <- NULL

# -----------------------------------------------------------------------------
# C3. Tabela de resultados
# -----------------------------------------------------------------------------

cat("\nTabela E&R (por media desc. | CHK=testemunha | L=linhagem candidata):\n\n")
cat(sprintf("%-4s %-8s %8s %7s %10s %5s  %7s %-9s  %7s %-9s\n",
            "Tipo", "Gen", "Media", "beta1", "sigma2_di", "R2", "t(b=1)", "p(t)", "F(s2=0)", "p(F)"))
cat(strrep("-", 90), "\n")

for (i in seq_len(nrow(er))) {
  cat(sprintf("%-4s %-8s %8.1f %7.3f %10.1f %5.1f  %7.3f %-9s  %7.3f %-9s\n",
              er$Tipo[i], er$Genotipo[i], er$Media[i], er$beta1[i], er$sigma2_di[i], er$R2[i],
              er$t_stat[i], paste0(fmt_p(er$p_t[i]), sig_stars(er$p_t[i])),
              er$F_stat[i], paste0(fmt_p(er$p_F[i]), sig_stars(er$p_F[i]))))
}

cat(strrep("-", 90), "\n")
cat("t(b=1): H0 beta1=1 (adapt. especifica se rejeitado)\n")
cat("F(s2=0): H0 sigma2_di=0 (comportamento imprevisto se rejeitado)\n")
cat("*** p<0.001  ** p<0.01  * p<0.05  . p<0.10\n")

# -----------------------------------------------------------------------------
# C4. Classificacao conjunta
# -----------------------------------------------------------------------------

er$Classe <- with(er, ifelse(
  Media < grand_mean, "Baixa media - pouco atrativo",
  ifelse(p_F < 0.05 & p_t < 0.05 & beta1 > 1, "Alta media | adapt. favoraveis | pouco previsivel",
  ifelse(p_F < 0.05 & p_t < 0.05 & beta1 < 1, "Alta media | adapt. desfavoraveis | pouco previsivel",
  ifelse(p_F < 0.05, "Alta media | adapt. geral | pouco previsivel",
  ifelse(p_t < 0.05 & beta1 > 1, "Alta media | adapt. favoraveis | previsivel",
  ifelse(p_t < 0.05 & beta1 < 1, "Alta media | adapt. desfavoraveis | previsivel",
         "Alta media | adapt. geral | previsivel"))))))
)

ideais <- er$Genotipo[er$Classe == "Alta media | adapt. geral | previsivel"] #genotipos ideais

cat("\nGenotipos: alta media + adapt. geral + comportamento previsivel:\n")
cat(" ", paste(ideais, collapse = ", "), "\n\n")
cat("Contagem por classe:\n")
print(table(er$Classe))

# -----------------------------------------------------------------------------
# C5. Retas coloridas pela classificacao
# -----------------------------------------------------------------------------

line_df_er <- expand.grid(
  gen = levels(dat$gen),
  I_j = seq(min(I_j) * 1.15, max(I_j) * 1.15, length.out = 200),
  stringsAsFactors = FALSE
)
line_df_er$gen <- factor(line_df_er$gen, levels = levels(dat$gen)) #mantem a ordem
line_df_er$pred <- predict(fit_fw, newdata = line_df_er) #predicoes do modelo da secao B
line_df_er$Classe <- er$Classe[match(as.character(line_df_er$gen), er$Genotipo)] #classe por genotipo

point_df_er <- mean_df
point_df_er$I_j <- I_j[as.character(point_df_er$env)] #indice ambiental nos pontos
point_df_er$Classe <- er$Classe[match(as.character(point_df_er$gen), er$Genotipo)] #classe dos pontos
point_df_er$Tipo   <- ifelse(point_df_er$gen %in% gen_chk, "CHK", "L") #tipo do genotipo

ideal_df <- line_df_er[line_df_er$gen %in% ideais, ] #subset dos ideais

class_cols <- c(
  "Alta media | adapt. geral | previsivel"               = "darkgreen",
  "Alta media | adapt. favoraveis | previsivel"          = "tomato3",
  "Alta media | adapt. desfavoraveis | previsivel"       = "steelblue",
  "Alta media | adapt. geral | pouco previsivel"         = "darkorange",
  "Alta media | adapt. favoraveis | pouco previsivel"    = "pink3",
  "Alta media | adapt. desfavoraveis | pouco previsivel" = "lightblue3",
  "Baixa media - pouco atrativo"                         = "gray75"
)

p_er <- ggplot() +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.6, color = "gray50") +
  geom_hline(yintercept = grand_mean, linetype = 2, linewidth = 0.6, color = "gray50") +
  geom_line(data = line_df_er, aes(x = I_j, y = pred, group = gen, color = Classe),
            linewidth = 0.45, alpha = 0.45) +
  geom_point(data = point_df_er, aes(x = I_j, y = y_mean, color = Classe, shape = Tipo),
             size = 1.7, alpha = 0.70) +
  geom_line(data = ideal_df, aes(x = I_j, y = pred, group = gen),
            color = "darkgreen", linewidth = 1.0, show.legend = FALSE) +
  scale_color_manual(values = class_cols) +
  scale_shape_manual(values = c("CHK" = 17, "L" = 16)) +
  labs(x = "I_j (kg/ha)", y = "Produtividade media (kg/ha)",
       title = "E&R (1966) | Classificacao conjunta",
       subtitle = "Triangulo = CHK | Ponto = linhagem | Cor = classe") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_er)

# -----------------------------------------------------------------------------
# C6. Inspecao individual: melhor e pior R2 entre genotipos de alta media
# -----------------------------------------------------------------------------

er_alta   <- er[er$Media >= grand_mean, ] #genotipos acima da media geral
nm_melhor <- er_alta$Genotipo[which.max(er_alta$R2)] #maior R2
nm_pior   <- er_alta$Genotipo[which.min(er_alta$R2)] #menor R2

cat(sprintf("\n--- C6. Inspecao visual individual ---\n"))
cat(sprintf("Maior R2 (mais previsivel): %s [%s]  R2=%.1f%%\n",
            nm_melhor, ifelse(nm_melhor %in% gen_chk, "CHK", "L"),
            er_alta$R2[er_alta$Genotipo == nm_melhor]))
cat(sprintf("Menor R2 (menos previsivel): %s [%s]  R2=%.1f%%\n",
            nm_pior, ifelse(nm_pior %in% gen_chk, "CHK", "L"),
            er_alta$R2[er_alta$Genotipo == nm_pior]))

sel_gen <- c(nm_melhor, nm_pior) #genotipos a inspecionar

inspect_pts <- data.frame(
  Genotipo = rep(sel_gen, each = a),
  I_j      = rep(as.numeric(I_j), times = 2),
  y_obs    = c(as.numeric(means_ge[nm_melhor, ]), as.numeric(means_ge[nm_pior, ])),
  y_fit    = c(mu_hat[nm_melhor] + b_hat[nm_melhor] * I_j,
               mu_hat[nm_pior]   + b_hat[nm_pior]   * I_j),
  env_lab  = rep(names(I_j), times = 2),
  stringsAsFactors = FALSE
)

inspect_pts$Titulo <- c(
  rep(paste0(nm_melhor, " [", ifelse(nm_melhor %in% gen_chk, "CHK", "L"), "]",
             "\n", "b=", round(er$beta1[er$Genotipo == nm_melhor], 2),
             " | s2d=", round(er$sigma2_di[er$Genotipo == nm_melhor], 0),
             " | R2=", round(er$R2[er$Genotipo == nm_melhor], 1), "%"), each = a),
  rep(paste0(nm_pior, " [", ifelse(nm_pior %in% gen_chk, "CHK", "L"), "]",
             "\n", "b=", round(er$beta1[er$Genotipo == nm_pior], 2),
             " | s2d=", round(er$sigma2_di[er$Genotipo == nm_pior], 0),
             " | R2=", round(er$R2[er$Genotipo == nm_pior], 1), "%"), each = a)
)

inspect_lines <- expand.grid(
  Genotipo = sel_gen,
  I_j      = seq(min(I_j), max(I_j), length.out = 200),
  stringsAsFactors = FALSE
)
inspect_lines$y_fit <- ifelse(inspect_lines$Genotipo == nm_melhor,
                              mu_hat[nm_melhor] + b_hat[nm_melhor] * inspect_lines$I_j,
                              mu_hat[nm_pior]   + b_hat[nm_pior]   * inspect_lines$I_j)

inspect_lines$Titulo <- ifelse(
  inspect_lines$Genotipo == nm_melhor,
  paste0(nm_melhor, " [", ifelse(nm_melhor %in% gen_chk, "CHK", "L"), "]",
         "\n", "b=", round(er$beta1[er$Genotipo == nm_melhor], 2),
         " | s2d=", round(er$sigma2_di[er$Genotipo == nm_melhor], 0),
         " | R2=", round(er$R2[er$Genotipo == nm_melhor], 1), "%"),
  paste0(nm_pior, " [", ifelse(nm_pior %in% gen_chk, "CHK", "L"), "]",
         "\n", "b=", round(er$beta1[er$Genotipo == nm_pior], 2),
         " | s2d=", round(er$sigma2_di[er$Genotipo == nm_pior], 0),
         " | R2=", round(er$R2[er$Genotipo == nm_pior], 1), "%")
)

p_inspect <- ggplot(inspect_pts, aes(x = I_j, y = y_obs)) +
  geom_segment(aes(x = I_j, xend = I_j, y = y_fit, yend = y_obs),
               color = "gray60", linewidth = 0.6) +
  geom_point(color = "steelblue", size = 2.2) +
  geom_line(data = inspect_lines, aes(x = I_j, y = y_fit),
            color = "tomato", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, color = "gray50") +
  geom_hline(data = data.frame(Titulo = unique(inspect_pts$Titulo),
                               y = c(mu_hat[nm_melhor], mu_hat[nm_pior])),
             aes(yintercept = y), linetype = 3, linewidth = 0.5, color = "gray50") +
  geom_text(aes(label = env_lab), vjust = -0.7, size = 3.0, color = "gray30") +
  facet_wrap(~ Titulo, scales = "free_y") +
  labs(x = "I_j (kg/ha)", y = "Produtividade (kg/ha)",
       title = "Inspecao individual de previsibilidade") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))
print(p_inspect)

cat("\nR2 alto: reta captura bem a resposta | R2 baixo: inspecionar grafico.\n")
cat("sigma2_di pode ser significativo mesmo com R2 alto (desvio sistematico).\n")
cat("A combinacao R2 + sigma2_di + inspecao grafica e mais informativa que cada um isolado.\n")

# -----------------------------------------------------------------------------
# C7. R^2 por genotipo
# -----------------------------------------------------------------------------

r2_df <- er[, c("Genotipo", "R2")]
r2_df$Grupo <- ifelse(r2_df$Genotipo %in% gen_chk, "CHK",
               ifelse(r2_df$R2 >= 70, "L R2>=70%", "L R2<70%"))
r2_df$Genotipo <- factor(r2_df$Genotipo, levels = r2_df$Genotipo[order(r2_df$R2)])

p_r2 <- ggplot(r2_df, aes(x = Genotipo, y = R2, fill = Grupo)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 70, linetype = 2, linewidth = 0.7, color = "gray40") +
  coord_flip() +
  scale_fill_manual(values = c("CHK" = "tomato3", "L R2>=70%" = "steelblue", "L R2<70%" = "gray65")) +
  labs(x = NULL, y = "R2 (%)",
       title = "R2_i por genotipo (E&R)",
       subtitle = "Vermelho = CHK | Azul = L R2>=70% | Cinza = L R2<70%") +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold"))
print(p_r2)

# =============================================================================
# SECAO D: DECOMPOSICAO DA SQ (E&R)
# SQ(A/G) = SQ(A linear) + SQ(GA linear) + SQDC
# =============================================================================

cat("\n\n=== SECAO D: DECOMPOSICAO DA SQ (E&R) ===\n")

SQ_AG_tot  <- r * sum((means_ge - rowMeans(means_ge))^2) #SQ total de ambientes dentro de genotipos
SQ_reg_all <- sum(r * (means_ge %*% I_j)^2 / sum_Ij2) #SQ total explicada pela regressao
SQ_A_lin   <- g * r * sum_Ij2 #componente linear medio dos ambientes
SQ_GA_lin  <- SQ_reg_all - SQ_A_lin #heterogeneidade das inclinacoes entre genotipos
SQDC       <- SQ_AG_tot - SQ_reg_all #desvio combinado nao explicado pela reta

QM_GA_lin  <- SQ_GA_lin / (g - 1) #QM da parte linear de GxA
QM_DC      <- SQDC / (g * (a - 2)) #QM do desvio combinado

F_GA <- QM_GA_lin / QMR #teste para heterogeneidade das inclinacoes
p_GA <- pf(F_GA, df1 = g - 1, df2 = nu, lower.tail = FALSE) #p correto da cauda superior

F_DC <- QM_DC / QMR #teste para desvio nao linear
p_DC <- pf(F_DC, df1 = g * (a - 2), df2 = nu, lower.tail = FALSE) #p correto da cauda superior

cat(sprintf("SQ(A/G) total ..............: %.1f\n", SQ_AG_tot))
cat(sprintf("SQ regressao linear ........: %.1f\n", SQ_reg_all))
cat(sprintf("SQ A linear .................: %.1f\n", SQ_A_lin))
cat(sprintf("SQ GA linear ................: %.1f\n", SQ_GA_lin))
cat(sprintf("SQ desvio combinado .........: %.1f\n", SQDC))
cat(sprintf("Verificacao SQ(A/G)=SQ_reg+SQDC: %.2f\n", SQ_AG_tot - SQ_reg_all - SQDC))

cat("\nLeitura inferencial:\n")
cat(sprintf("GA linear ........ F=%.2f | p=%s%s\n", F_GA, fmt_p(p_GA), sig_stars(p_GA)))
cat(sprintf("Desvio combinado . F=%.2f | p=%s%s\n", F_DC, fmt_p(p_DC), sig_stars(p_DC)))
cat("Interpretacao: GA linear significativo indica diferencas em adaptabilidade.\n")
cat("Interpretacao: desvio significativo indica componente nao linear em torno da reta.\n")

sq_df <- data.frame(
  Componente = factor(c("A linear", "GA linear", "Desvio combinado"),
                      levels = c("A linear", "GA linear", "Desvio combinado")),
  SQ = c(SQ_A_lin, SQ_GA_lin, SQDC)
)

p_sq <- ggplot(sq_df, aes(x = Componente, y = SQ, fill = Componente)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = round(SQ, 0)), vjust = -0.35, size = 3.4) +
  scale_fill_manual(values = c("gray70", "steelblue", "tomato")) +
  labs(x = NULL, y = "Soma de quadrados",
       title = "Decomposicao da SQ(A/G) em E&R") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))
print(p_sq)

# =============================================================================
# SECAO E: CAUTELAS
# =============================================================================

cat("\n\n=== SECAO E: CAUTELAS ===\n")

R2_vec <- er$R2; names(R2_vec) <- er$Genotipo #vetor nomeado de R2

cat(sprintf("\nLinearidade da resposta:\n"))
cat(sprintf("R2 medio = %.1f%%\n", mean(er$R2)))
cat(sprintf("Menor R2 = %.1f%% (%s)\n", min(er$R2), er$Genotipo[which.min(er$R2)]))
cat(sprintf("Maior R2 = %.1f%% (%s)\n", max(er$R2), er$Genotipo[which.max(er$R2)]))

bad_r2 <- sort(R2_vec[R2_vec < 70]) #genotipos com ajuste linear fraco
cat("\nGenotipos com R2 < 70%:\n")
if (length(bad_r2) > 0) print(round(bad_r2, 1)) else cat("Nenhum.\n")

# L027: b_i proximo de zero nao significa estabilidade; a reta nao explica quase nada.
# L028: b_i alto com R2 baixo e sigma2_di alto exige cautela na interpretacao.

env_rem <- names(which.max(abs(I_j))) #ambiente mais extremo no gradiente
means_red <- means_ge[, colnames(means_ge) != env_rem] #remove esse ambiente
I_j_red <- colMeans(means_red) - mean(means_red) #recalcula o indice apos remocao
b_red <- as.vector((means_red %*% I_j_red) / sum(I_j_red^2)) #reestima b_i sem o ambiente extremo
names(b_red) <- rownames(means_red) #recoloca nomes dos genotipos
delta_b <- b_red - b_hat[names(b_red)] #impacto da remocao sobre b_i

cat(sprintf("\nAlavancagem ambiental:\n"))
cat(sprintf("Ambiente mais extremo = %s\n", env_rem))
cat(sprintf("I_j desse ambiente = %.1f\n", I_j[env_rem]))
cat(sprintf("Maior mudanca em b_i apos remocao = %.3f (%s)\n",
            max(abs(delta_b)), names(which.max(abs(delta_b)))))
cat(sprintf("Mudanca media absoluta em b_i = %.3f\n", mean(abs(delta_b))))

cat(sprintf("\nDependencia do indice:\n"))
cat(sprintf("Cada genotipo contribui %.3f ao proprio I_j\n", 1 / g))
cat("Leitura: b_i e relativo ao painel avaliado, nao uma propriedade absoluta.\n")

cat(sprintf("\nNumero e amplitude:\n"))
cat(sprintf("Ambientes = %d | gl do desvio por genotipo = %d\n", a, a - 2))
cat(sprintf("Amplitude observada de I_j = %.1f kg/ha\n", diff(range(I_j))))
cat("Leitura: poucos ambientes ou pouca amplitude dificultam a discriminacao das inclinacoes.\n")

cat("\nLeitura conjunta dos parametros:\n")
cat("1. Media = nivel de producao\n")
cat("2. sigma2_di = previsibilidade da resposta\n")
cat("3. beta1 = perfil de adaptacao\n")
cat("Ordem pratica: media > previsibilidade > perfil de adaptacao\n")

caut_df <- data.frame(
  Genotipo = factor(names(bad_r2), levels = names(bad_r2)),
  R2 = as.numeric(bad_r2)
)

if (nrow(caut_df) > 0) {
  p_caut <- ggplot(caut_df, aes(x = Genotipo, y = R2)) +
    geom_col(fill = "gray65", width = 0.7) +
    geom_hline(yintercept = 70, linetype = 2, linewidth = 0.7, color = "gray40") +
    coord_flip() +
    labs(x = NULL, y = "R2 (%)",
         title = "Genotipos com ajuste linear questionavel") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"))
  
  print(p_caut)
}

# =============================================================================
# SECAO F: FREEMAN E PERKINS (1971)
# =============================================================================

cat("\n\n=== SECAO F: FREEMAN & PERKINS (1971) ===\n")
cat("Comparacao entre indice endogeno e indice exogeno baseado apenas em CHK.\n\n")

means_chk <- means_ge[gen_chk, ] #medias das testemunhas
I_j_exo <- colMeans(means_chk) - mean(means_chk) #indice ambiental exogeno

idx_cmp <- data.frame(
  Ambiente = names(I_j),
  Local    = env_info$location[match(names(I_j), env_info$env)],
  I_endo   = round(I_j, 1),
  I_exo    = round(I_j_exo, 1),
  Dif      = round(I_j_exo - I_j, 1)
)

cat("Comparacao dos indices:\n")
print(idx_cmp, row.names = FALSE)
cat(sprintf("\nCor(I_endo, I_exo) = %.4f\n", cor(I_j, I_j_exo)))

b_exo <- as.vector((means_ge[gen_lin, ] %*% I_j_exo) / sum(I_j_exo^2)) #b_i das linhagens com indice exogeno
names(b_exo) <- gen_lin
b_endo_lin <- b_hat[gen_lin] #b_i das linhagens com indice endogeno

cat(sprintf("Cor(b_endo, b_exo) = %.4f\n", cor(b_endo_lin, b_exo)))
cat("Leitura: alta correlacao indica que o vies do indice endogeno e pequeno neste ensaio.\n")

df_fp <- data.frame(b_endo = b_endo_lin, b_exo = b_exo, gen = gen_lin)

p_fp <- ggplot(df_fp, aes(b_endo, b_exo, label = gen)) +
  geom_point(color = "steelblue", size = 2.0) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.7, color = "gray40") +
  geom_hline(yintercept = 1, linetype = 3, linewidth = 0.6, color = "tomato") +
  geom_vline(xintercept = 1, linetype = 3, linewidth = 0.6, color = "tomato") +
  geom_text(size = 2.5, vjust = -0.5, color = "gray30") +
  theme_bw() +
  labs(x = "b_i endogeno", y = "b_i exogeno",
       title = "Freeman & Perkins (1971)",
       subtitle = "Linhagens: indice endogeno vs indice exogeno")
print(p_fp)

cat("\nMensagem final desta secao:\n")
cat("O indice exogeno e conceitualmente melhor.\n")
cat("Na pratica, com muitos genotipos, os resultados tendem a ser muito proximos.\n")
cat("CHK fixas ajudam na comparabilidade entre redes e safras.\n")

# =============================================================================
# RESUMO FINAL
# =============================================================================

cat("\n\n=============================================================\n")
cat("RESUMO FINAL - AULA 3\n")
cat("=============================================================\n")

cat(sprintf("Estrutura do ensaio: %d genotipos (%d CHK + %d L) | %d ambientes | %d blocos\n",
            g, length(gen_chk), length(gen_lin), a, r))
cat(sprintf("Media geral: %.1f kg/ha\n", grand_mean))
cat(sprintf("QMR da anova conjunta: %.2f (gl=%d)\n", QMR, nu))
cat(sprintf("Faixa do indice ambiental: %.1f a %.1f kg/ha\n", min(I_j), max(I_j)))
cat(sprintf("Amplitude do gradiente: %.1f kg/ha\n", diff(range(I_j))))
cat(sprintf("Menor b_i: %.2f (%s)\n", min(b_hat), names(which.min(b_hat))))
cat(sprintf("Maior b_i: %.2f (%s)\n", max(b_hat), names(which.max(b_hat))))

cat("\nGenotipos ideais segundo E&R:\n")
cat(" ", paste(ideais, collapse = ", "), "\n")

cat("\nContagem por classe:\n")
print(table(er$Classe))

cat(sprintf("\nGA linear ........ F=%.2f | p=%s%s\n", F_GA, fmt_p(p_GA), sig_stars(p_GA)))
cat(sprintf("Desvio combinado . F=%.2f | p=%s%s\n", F_DC, fmt_p(p_DC), sig_stars(p_DC)))
cat(sprintf("Cor(I_endo, I_exo) = %.4f\n", cor(I_j, I_j_exo)))
cat(sprintf("Cor(b_endo, b_exo) = %.4f\n", cor(b_endo_lin, b_exo)))

cat("\nCautelas principais:\n")
cat("1. Ler media antes de adaptabilidade.\n")
cat("2. Ler sigma2_di antes de concluir sobre recomendacao.\n")
cat("3. R2 ajuda a detectar quando a reta resume mal a resposta.\n")
cat("4. Ambientes extremos podem alterar b_i de forma relevante.\n")
cat("5. b_i e relativo ao conjunto de genotipos e ambientes analisado.\n")
cat("=============================================================\n")
