# =============================================================================
# GMP0164 - Biometria da Interacao Genotipos x Ambientes
# Aula 6: Modelos Mistos I - Introducao ao BLUP Multiambiente
# Dados: maize_safrinha_MET.txt
#   70 hibridos de milho safrinha (NC Design II + 8 checks comerciais)
#   6 ambientes do Cerrado | DBC, 3 blocos | desbalanceamento estrategico
#   traits: yield = "Grain yield (kg/ha)" | maturity = "Days to physiological maturity"
# Prof. Dr. Rafael Tassinari Resende
# =============================================================================

library(pedigreemm)
library(lme4)
library(lme4breeding)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Carregar dados
# -----------------------------------------------------------------------------

dat <- read.table("G:/Meu Drive/UFG/PPGGMP/BiometriaDaInteracaoGxA/aulas/simu_data/maize_safrinha_MET.txt", header=TRUE, sep="\t",
                  stringsAsFactors=FALSE)
dat$gen   <- factor(dat$gen)
dat$trial <- factor(dat$trial)
dat$block <- factor(dat$block)

# -----------------------------------------------------------------------------
# 2. Dimensoes e descritivas
# -----------------------------------------------------------------------------

I <- nlevels(dat$gen); J <- nlevels(dat$trial); K <- nlevels(dat$block)
n_exp <- length(unique(dat$gen[dat$gen_type=="experimental"]))
n_chk <- length(unique(dat$gen[dat$gen_type=="check"]))
 
cat(sprintf("Genotipos: %d (%d exp + %d checks) | Ambientes: %d | Blocos: %d\n",
            I, n_exp, n_chk, J, K))
cat(sprintf("Obs: %d / %d balanceadas (%.1f%% missing)\n",
            nrow(dat), I*J*K, 100*(1-nrow(dat)/(I*J*K))))
 
# Ambientes por genotipo — experimentais e checks
amb_por_gen <- tapply(dat$trial[dat$gen_type=="experimental"],
                      droplevels(dat$gen[dat$gen_type=="experimental"]),
                      function(x) length(unique(x)))
amb_por_chk <- tapply(dat$trial[dat$gen_type=="check"],
                      droplevels(dat$gen[dat$gen_type=="check"]),
                      function(x) length(unique(x)))
cat("Distribuicao de ambientes por genotipo experimental:\n")
print(amb_por_gen)
cat("Ambientes por testemunha (checks):\n")
print(amb_por_chk)
 
# Medias por ambiente — ordenar crescente (estressante -> favoravel)
med_amb  <- sort(tapply(dat$yield, dat$trial, mean))
env_order <- names(med_amb)   # ordem usada nos graficos
cat("\nMedias por ambiente (kg/ha):\n"); print(round(med_amb))
 
# Medias por genotipo — todos (exp + checks)
med_gen <- tapply(dat$yield, dat$gen, mean)
cat("\nTop 15 genotipos — media aritmetica (kg/ha):\n")
print(round(sort(med_gen, decreasing=TRUE)[1:15]))
 
# Tabela GxA (medias por genotipo x ambiente)
tab_ga <- tapply(dat$yield, list(dat$gen, dat$trial), mean)[ , env_order]
cat("\nTabela GxA — primeiras 6 linhas:\n")
print(round(tab_ga[1:6, ], 0))

# -----------------------------------------------------------------------------
# 3. Grafico de interacao GxA — ggplot2
# X: do menos favoravel (esquerda) para o mais favoravel (direita)
# -----------------------------------------------------------------------------

exp_ids    <- sort(unique(as.character(dat$gen[dat$gen_type=="experimental"])))
check_ids  <- sort(unique(as.character(dat$gen[dat$gen_type=="check"])))
trial_city <- unique(dat[, c("trial","city")])
city_order <- trial_city$city[match(env_order, trial_city$trial)]

chk_colors <- c("#E63946","#F4A261","#2A9D8F","#457B9D",
                "#9B2226","#BB3E03","#005F73","#0A9396")
chk_names  <- unique(dat[dat$gen_type=="check", c("gen","cultivar_name")])
chk_names  <- chk_names[order(chk_names$gen), ]

# Formato longo para ggplot
ga_long        <- as.data.frame(as.table(tab_ga))
names(ga_long) <- c("gen","trial","yield")
ga_long        <- ga_long[!is.na(ga_long$yield), ]
ga_long$gen    <- as.character(ga_long$gen)
ga_long$trial  <- as.character(ga_long$trial)
ga_long$city   <- city_order[match(ga_long$trial, env_order)]
ga_long$city   <- factor(ga_long$city, levels=city_order)  # estressante->favoravel
ga_long$tipo   <- ifelse(ga_long$gen %in% check_ids, "check", "experimental")
ga_long$label  <- chk_names$cultivar_name[match(ga_long$gen, chk_names$gen)]

chk_long <- ga_long[ga_long$tipo == "check", ]
chk_long$leg <- paste(chk_long$gen, chk_long$label, sep=" - ")
color_map  <- setNames(chk_colors, sort(unique(chk_long$leg)))

ggplot() +
  geom_hline(yintercept=mean(dat$yield), linetype="dashed",
             color="gray50", linewidth=0.4) +
  geom_line(data=ga_long[ga_long$tipo=="experimental", ],
            aes(x=city, y=yield, group=gen),
            color="gray75", linewidth=0.35, alpha=0.6) +
  geom_line(data=chk_long,
            aes(x=city, y=yield, group=leg, color=leg),
            linewidth=1.1) +
  geom_point(data=chk_long,
             aes(x=city, y=yield, color=leg),
             size=2.2) +
  scale_color_manual(values=color_map, name="Testemunha") +
  scale_y_continuous(labels=scales::comma, breaks=seq(2000,12000,by=1000)) +
  labs(title="Interacao GxA — Milho safrinha MET 2025",
       subtitle="Cinza: 62 hibridos experimentais | Colorido: 8 testemunhas comerciais",
       x="Ambiente (estressante -> favoravel)", y="Yield medio (kg/ha)") +
  theme_bw(base_size=12) +
  theme(axis.text.x   = element_text(angle=25, hjust=1),
        legend.text   = element_text(size=8),
        legend.title  = element_text(size=9),
        panel.grid.minor = element_blank())

# -----------------------------------------------------------------------------
# 4. Modelo misto I — todos os genotipos aleatorios (G = I)
# -----------------------------------------------------------------------------

mod1 <- lmer(yield ~ -1 + trial + (1|gen) + (1|gen:trial) + (1|trial:block),
             data=dat, REML=TRUE)
cat("\n--- Componentes de variancia (mod1: todos aleatorios, G = I) ---\n")
print(VarCorr(mod1))

vc   <- as.data.frame(VarCorr(mod1))
s2g  <- vc$vcov[vc$grp=="gen"]
s2ga <- vc$vcov[vc$grp=="gen:trial"]
s2b  <- vc$vcov[vc$grp=="trial:block"]
s2e  <- vc$vcov[vc$grp=="Residual"]
s2P  <- s2g + s2ga + s2b + s2e

cat(sprintf("V_G=%.0f(%.1f%%) | V_GxA=%.0f(%.1f%%) | V_Blk=%.0f(%.1f%%) | V_Err=%.0f(%.1f%%)\n",
            s2g,100*s2g/s2P, s2ga,100*s2ga/s2P, s2b,100*s2b/s2P, s2e,100*s2e/s2P))
cat(sprintf("Razao V_GxA / V_G = %.2f\n", s2ga/s2g))

# H2 em base de medias
a_med <- mean(amb_por_gen)
H2    <- s2g / (s2g + s2ga/a_med + s2e/(a_med*K))
cat(sprintf("H2 (base de medias, %.1f amb x %d blocos): %.3f\n", a_med, K, H2))

# -----------------------------------------------------------------------------
# 5. Modelo misto II — checks como efeito fixo, experimentais aleatorios
# -----------------------------------------------------------------------------

# chk_fix: CHKxx para testemunhas, "EXP" para experimentais (nivel de referencia)
# gen_r  : ID do genotipo so para experimentais (NA para checks)
dat$chk_fix <- factor(ifelse(dat$gen_type=="check", as.character(dat$gen), "EXP"))
dat$gen_ran   <- factor(ifelse(dat$gen_type=="experimental", as.character(dat$gen), 1))

mod2 <- lmer(yield ~ -1 + trial + chk_fix + (1|gen_ran) + (1|gen:trial) + (1|trial:block),
             data=dat, REML=TRUE)

cat("\n--- Componentes de variancia (mod2: checks fixos) ---\n")
print(VarCorr(mod2))

# BLUEs dos checks (fixos) — desvio em relacao ao nivel "EXP"
fe_chk <- fixef(mod2)[grep("chk_fix", names(fixef(mod2)))]
names(fe_chk) <- gsub("chk_fix","", names(fe_chk))
cat("\n--- BLUEs dos checks (desvio relativo ao grupo EXP) ---\n")
print(round(sort(fe_chk, decreasing=TRUE), 1))

# -----------------------------------------------------------------------------
# 6. BLUPs e PEVs (mod1)
# -----------------------------------------------------------------------------

re1 <- ranef(mod1, condVar=TRUE)

# --- Genotipo ---
blup_g           <- re1$gen[ ,1]
pev_g            <- drop(attr(re1$gen, "postVar"))
acc_g            <- sqrt(1 - pev_g/s2g)   # acuracia seletiva
blup_gen_df      <- data.frame(gen=rownames(re1$gen), blup=blup_g,
                               pev=pev_g, acc=acc_g, row.names=NULL)
blup_gen_df      <- blup_gen_df[order(-blup_gen_df$blup), ]
cat("\n--- BLUPs de genotipo + PEV + acuracia (top 15) ---\n")
print(head(blup_gen_df, 15), digits=3)

# --- Interacao GxA ---
blup_ga_val  <- re1$`gen:trial`[ ,1]
pev_ga_val   <- drop(attr(re1$`gen:trial`, "postVar"))
acc_ga_val   <- sqrt(1 - pev_ga_val/s2ga)
blup_ga_df   <- data.frame(gen_trial=rownames(re1$`gen:trial`),
                            blup_ga=blup_ga_val, pev_ga=pev_ga_val,
                            acc_ga=acc_ga_val, row.names=NULL)
blup_ga_df   <- blup_ga_df[order(-abs(blup_ga_df$blup_ga)), ]
cat("\n--- BLUPs GxA + PEV + acuracia (maiores interacoes, top 15) ---\n")
print(head(blup_ga_df, 15), digits=3)

cat("\n--- Resumo PEV genotipo ---\n")
print(round(quantile(blup_gen_df$pev, c(0,.25,.5,.75,1))))
cat("--- Resumo PEV GxA ---\n")
print(round(quantile(blup_ga_df$pev_ga, c(0,.25,.5,.75,1))))

# -----------------------------------------------------------------------------
# 7. BLUP vs. media aritmetica — comparacao central (experimentais)
# -----------------------------------------------------------------------------

blup_exp <- blup_gen_df[blup_gen_df$gen %in% exp_ids, ]
comp     <- merge(blup_exp,
                  data.frame(gen=names(med_gen), media=med_gen),
                  by="gen")
comp$n_amb      <- amb_por_gen[comp$gen]
comp$rank_blup  <- rank(-comp$blup)
comp$rank_media <- rank(-comp$media)
comp$delta_rank <- abs(comp$rank_blup - comp$rank_media)

cat(sprintf("\nCorrelacao de Spearman BLUP x media: %.3f\n",
            cor(comp$blup, comp$media, method="spearman")))
cat("Maior mudanca de ranking (desbalanceamento):\n")
print(head(comp[order(-comp$delta_rank), ], 10))

# -----------------------------------------------------------------------------
# 8. Grafico de encolhimento — facet por ambiente
# Media (raw, por ambiente) vs BLUP predito (fixef + ranef) por ambiente
# X: "Media" | "BLUP" | Y: yield (kg/ha) | linhas = genotipos experimentais
# -----------------------------------------------------------------------------

# Media raw por genotipo x ambiente
med_ge <- aggregate(yield ~ gen + trial, dat[dat$gen_type=="experimental", ], mean)

# BLUP predito = fixef[trial] + ranef[gen] + ranef[gen:trial]
fe      <- fixef(mod1)
blup_g  <- setNames(re1$gen[,1],          rownames(re1$gen))
blup_ga <- setNames(re1$`gen:trial`[,1],  rownames(re1$`gen:trial`))

med_ge$blup_pred <- mapply(function(g, tr) {
  env_fe  <- fe[paste0("trial", tr)]
  ga_key  <- paste(g, tr, sep=":")
  env_fe + blup_g[g] + ifelse(ga_key %in% names(blup_ga), blup_ga[ga_key], 0)
}, as.character(med_ge$gen), as.character(med_ge$trial))

# Formato longo
long_sh <- rbind(
  data.frame(gen=med_ge$gen, trial=med_ge$trial, tipo="Media",    valor=med_ge$yield),
  data.frame(gen=med_ge$gen, trial=med_ge$trial, tipo="BLUP pred",valor=med_ge$blup_pred)
)
long_sh$tipo  <- factor(long_sh$tipo, levels=c("Media","BLUP pred"))
long_sh$trial <- factor(long_sh$trial, levels=env_order)  # estressante->favoravel

# Labeller: trial -> cidade
city_lab <- setNames(city_order, env_order)

ggplot(long_sh, aes(x=tipo, y=valor, group=gen)) +
  geom_line(color="gray60", alpha=0.45, linewidth=0.4) +
  geom_point(aes(color=tipo), size=1, alpha=0.7) +
  facet_wrap(~trial, nrow=2,
             labeller=labeller(trial=as_labeller(city_lab))) +
  scale_color_manual(values=c("Media"="#457B9D","BLUP pred"="#E63946"), guide="none") +
  labs(title="Encolhimento (shrinkage) — Media raw vs BLUP predito por ambiente",
       subtitle="Cada linha = 1 genotipo experimental | GxA visivel no cruzamento das linhas",
       x=NULL, y="Yield (kg/ha)") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(size=8))

# -----------------------------------------------------------------------------
# 9. Extensao: incorporando parentesco (G = A) via lme4breeding
# -----------------------------------------------------------------------------

# Pedigree NC Design II
ped_exp  <- unique(dat[dat$gen_type=="experimental", c("gen","female","male")])
ped_fund <- data.frame(gen=c(unique(ped_exp$female), unique(ped_exp$male)), female=0, male=0)
ped_chk  <- data.frame(gen=unique(dat$gen[dat$gen_type=="check"]),          female=0, male=0)
ped      <- rbind(ped_fund, ped_exp, ped_chk)

A <- as.matrix(getA(pedigree(sire=ped$male, dam=ped$female, label=ped$gen)))
heatmap(A)
cat(sprintf("Matriz A: %dx%d | diag: %.2f-%.2f | off-diag: %.3f-%.3f\n",
            nrow(A), ncol(A), min(diag(A)), max(diag(A)),
            min(A[lower.tri(A)]), max(A[lower.tri(A)])))

# Matriz A ⊗ I_J para o termo gen:trial
# Logica: dois niveis gen:trial so sao correlacionados se estiverem no mesmo ambiente
dat$gen_trial <- paste(dat$gen, dat$trial, sep=":")
gxt           <- levels(factor(dat$gen_trial))
g_idx         <- sub(":.*", "", gxt)
t_idx         <- sub(".*:", "", gxt)

A_GxA <- matrix(0, length(gxt), length(gxt), dimnames=list(gxt, gxt))
for(i in seq_along(gxt))
  for(j in seq_along(gxt))
    if(t_idx[i]==t_idx[j]) A_GxA[i,j] <- A[g_idx[i], g_idx[j]]

# Modelos
mod_A  <- lmebreed(yield ~ -1 + trial + (1|gen) + (1|gen:trial) + (1|trial:block),
                   relmat=list(gen=A),                  data=dat, REML=TRUE)
mod_AK <- lmebreed(yield ~ -1 + trial + (1|gen) + (1|gen_trial) + (1|trial:block),
                   relmat=list(gen=A, gen_trial=A_GxA), data=dat, REML=TRUE)

cat("\n--- VarCorr: G=A, GxA=I ---\n");    print(VarCorr(mod_A))
cat("\n--- VarCorr: G=A, GxA=A⊗I_J ---\n"); print(VarCorr(mod_AK))

# Comparacao de BLUPs (genotipos)
blup_A  <- data.frame(gen=rownames(ranef(mod_A)$gen),  blup_A =ranef(mod_A)$gen[,1])
blup_AK <- data.frame(gen=rownames(ranef(mod_AK)$gen), blup_AK=ranef(mod_AK)$gen[,1])
comp2   <- merge(merge(blup_gen_df[,c("gen","blup")], blup_A, by="gen"), blup_AK, by="gen")

cat(sprintf("Corr BLUP(I)  x BLUP(A)  : %.4f\n", cor(comp2$blup,   comp2$blup_A)))
cat(sprintf("Corr BLUP(I)  x BLUP(AK) : %.4f\n", cor(comp2$blup,   comp2$blup_AK)))
cat(sprintf("Corr BLUP(A)  x BLUP(AK) : %.4f\n", cor(comp2$blup_A, comp2$blup_AK)))

# -----------------------------------------------------------------------------
# 10. Grafico comparativo — 4 paineis: medias vs modelos
# -----------------------------------------------------------------------------

marg <- function(mod, dat, gen_term="gen", gxa_term="gen:trial") {
  fe  <- fixef(mod);  re <- ranef(mod)
  bg  <- setNames(re[[gen_term]][,1],  rownames(re[[gen_term]]))
  bga <- setNames(re[[gxa_term]][,1],  rownames(re[[gxa_term]]))
  mapply(function(g, tr) {
    fe[paste0("trial",tr)] +
      ifelse(g %in% names(bg),  bg[g], 0) +
      ifelse(paste(g,tr,sep=":") %in% names(bga), bga[paste(g,tr,sep=":")], 0)
  }, as.character(dat$gen), as.character(dat$trial))
}

marg2 <- function(mod, dat) {
  fe  <- fixef(mod);  re <- ranef(mod)
  bg  <- setNames(re[["gen_ran"]][,1],   rownames(re[["gen_ran"]]))
  bga <- setNames(re[["gen:trial"]][,1], rownames(re[["gen:trial"]]))
  mapply(function(g, tr, gtype) {
    fe[paste0("trial",tr)] +
      ifelse(gtype=="check", fe[paste0("chk_fix",g)], 0) +
      bg[ifelse(gtype=="experimental", g, "1")] +
      ifelse(paste(g,tr,sep=":") %in% names(bga), bga[paste(g,tr,sep=":")], 0)
  }, as.character(dat$gen), as.character(dat$trial), dat$gen_type)
}

dat$fit1  <- marg(mod1,  dat)
dat$fit2  <- marg2(mod2, dat)
dat$fitAK <- marg(mod_AK, dat, gxa_term="gen_trial")
fit_ge    <- aggregate(cbind(fit1,fit2,fitAK) ~ gen+trial, dat, mean)

# Formato longo
med_long        <- as.data.frame(as.table(tab_ga))
names(med_long) <- c("gen","trial","valor")
med_long        <- med_long[!is.na(med_long$valor), ]

mk <- function(df, col, lab)
  data.frame(gen=df$gen, trial=df$trial, valor=df[[col]], modelo=lab)

all_long <- rbind(
  mk(med_long, "valor",  "1. Medias raw"),
  mk(fit_ge,   "fit1",   "2. Mod1 — G = I"),
  mk(fit_ge,   "fit2",   "3. Mod2 — Checks fixos"),
  mk(fit_ge,   "fitAK",  "4. Mod_AK — G = A")
)
all_long$city   <- factor(city_order[match(all_long$trial, env_order)], levels=city_order)
all_long$tipo   <- ifelse(as.character(all_long$gen) %in% check_ids, "check", "experimental")
all_long$modelo <- factor(all_long$modelo,
                          levels=c("1. Medias raw","2. Mod1 — G = I",
                                   "3. Mod2 — Checks fixos","4. Mod_AK — G = A"))

ggplot(all_long, aes(x=city, y=valor, group=gen)) +
  geom_line(data=all_long[all_long$tipo=="experimental",],
            color="gray75", linewidth=0.3, alpha=0.6) +
  geom_line(data=all_long[all_long$tipo=="check",],
            aes(color=gen), linewidth=1.0) +
  facet_wrap(~modelo, nrow=1) +
  scale_color_manual(values=setNames(chk_colors, check_ids), name="Check") +
  scale_y_continuous(labels=scales::comma) +
  labs(title="Comparacao de ajustes por modelo — GxA Milho safrinha MET 2025",
       subtitle="Cinza: experimentais | Colorido: testemunhas comerciais",
       x="Ambiente (estressante -> favoravel)", y="Yield (kg/ha)") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=7),
        legend.text=element_text(size=7), strip.text=element_text(face="bold"),
        legend.position="bottom", panel.grid.minor=element_blank())

# Pairs panel
wide <- merge(
  data.frame(gen=as.character(med_long$gen), trial=as.character(med_long$trial),
             tipo=ifelse(med_long$gen %in% check_ids, "check", "experimental"),
             Medias=med_long$valor),
  fit_ge[, c("gen","trial","fit1","fit2","fitAK")],
  by=c("gen","trial"), all.x=TRUE
)
names(wide)[5:7] <- c("Mod1 G=I","Mod2 Chk-fix","Mod_AK G=A")

pares <- list(c("Medias","Mod1 G=I"), c("Medias","Mod2 Chk-fix"),
              c("Medias","Mod_AK G=A"),
              c("Mod1 G=I","Mod2 Chk-fix"),
              c("Mod1 G=I","Mod_AK G=A"),
              c("Mod2 Chk-fix","Mod_AK G=A"))

pl <- do.call(rbind, lapply(pares, function(p)
  data.frame(xlab=p[1], ylab=p[2], x=wide[[p[1]]], y=wide[[p[2]]],
             tipo=wide$tipo, gen=wide$gen)))
pl$panel <- factor(paste(pl$xlab,"vs",pl$ylab), levels=unique(paste(pl$xlab,"vs",pl$ylab)))

ggplot(pl, aes(x=x, y=y)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray40", linewidth=0.5) +
  geom_point(data=pl[pl$tipo=="experimental",], color="gray75", size=1.2, alpha=0.5) +
  geom_point(data=pl[pl$tipo=="check",], aes(color=gen), size=2.8, alpha=0.9) +
  scale_color_manual(values=setNames(chk_colors, check_ids), name="Check") +
  facet_wrap(~panel, nrow=2, scales="free") +
  labs(title="Comparacao par-a-par entre modelos — GxA Milho safrinha MET 2025",
       subtitle="1 ponto = 1 genotipo x ambiente | diagonal = concordancia perfeita",
       x=NULL, y=NULL) +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(size=8, face="bold"),
        legend.position="bottom", legend.direction="horizontal")

# -----------------------------------------------------------------------------
# 11. Sumario final
# -----------------------------------------------------------------------------

cat("\n=============================================================\n")
cat("SUMARIO FINAL — Aula 6: Modelos Mistos I\n")
cat("=============================================================\n")
cat(sprintf("Genotipos: %d | Ambientes: %d | Obs: %d (%.1f%% missing)\n",
            I, J, nrow(dat), 100*(1-nrow(dat)/(I*J*K))))
cat(sprintf("V_GxA / V_G = %.2f | H2 (medias) = %.3f\n", s2ga/s2g, H2))
cat(sprintf("Melhor BLUP(I)  : %s (%.0f kg/ha) | acc = %.3f\n",
            blup_gen_df$gen[1], blup_gen_df$blup[1], blup_gen_df$acc[1]))
cat(sprintf("Melhor BLUP(AK) : %s (%.0f kg/ha)\n",
            cmp_pat$gen[which.max(cmp_pat$blup_AK)],
            max(cmp_pat$blup_AK)))
cat(sprintf("Corr BLUP(I) x BLUP(AK) : %.4f\n",
            cor(cmp_pat$blup_1, cmp_pat$blup_AK)))
cat("=============================================================\n")



# -----------------------------------------------------------------------------
# 12. [PLUS] Parentesco (tende) a puxar individuos para a media da familia paterna (Mod1 vs Mod_AK)
# -----------------------------------------------------------------------------

# BLUPs experimentais — mod1 e mod_AK
re1   <- ranef(mod1);   reAK  <- ranef(mod_AK)
b1    <- setNames(re1$gen[,1],   rownames(re1$gen))
bAK   <- setNames(reAK$gen[,1], rownames(reAK$gen))

fam_pat <- unique(dat[dat$gen_type=="experimental", c("gen","male")])

cmp_pat <- data.frame(
  gen    = fam_pat$gen,
  male   = fam_pat$male,
  blup_1 = b1[as.character(fam_pat$gen)],
  blup_AK= bAK[as.character(fam_pat$gen)],
  n_amb  = amb_por_gen[as.character(fam_pat$gen)]
)

fam1_pat  <- aggregate(blup_1  ~ male, cmp_pat, mean); names(fam1_pat)[2]  <- "fam1"
famAK_pat <- aggregate(blup_AK ~ male, cmp_pat, mean); names(famAK_pat)[2] <- "famAK"

bar_df <- rbind(
  data.frame(male=fam1_pat$male,  blup=fam1_pat$fam1,   modelo="Mod1 — G=I"),
  data.frame(male=famAK_pat$male, blup=famAK_pat$famAK, modelo="Mod_AK — G=A")
)
bar_df$modelo <- factor(bar_df$modelo, levels=c("Mod1 — G=I","Mod_AK — G=A"))

pt_df <- rbind(
  data.frame(male=cmp_pat$male, blup=cmp_pat$blup_1,  n_amb=cmp_pat$n_amb, modelo="Mod1 — G=I"),
  data.frame(male=cmp_pat$male, blup=cmp_pat$blup_AK, n_amb=cmp_pat$n_amb, modelo="Mod_AK — G=A")
)
pt_df$modelo <- factor(pt_df$modelo, levels=c("Mod1 — G=I","Mod_AK — G=A"))

ggplot(bar_df, aes(x=male)) +
  geom_hline(yintercept=0, linewidth=0.4, color="gray40") +
  geom_col(aes(y=blup, fill=modelo),
           position=position_dodge(0.75), width=0.7, alpha=0.3) +
  geom_point(data=pt_df, aes(y=blup, color=n_amb, group=modelo),
             position=position_dodge(0.75), size=2.5, alpha=0.85) +
  scale_fill_manual(values=c("Mod1 — G=I"="#457B9D","Mod_AK — G=A"="#2A9D8F"), name=NULL) +
  scale_color_gradientn(colors=c("#d73027","#fee08b","#1a9850"),
                        name="N ambientes", breaks=2:6) +
  labs(title="BLUPs por familia paterna — Mod1 vs Mod_AK",
       subtitle="Barra = media familiar | Pontos = individuos | Cor = n ambientes observados",
       x="Familia paterna", y="BLUP (kg/ha)") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=40, hjust=1, size=8),
        legend.position="right", panel.grid.minor=element_blank())
