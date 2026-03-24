# =============================================================================
# GMP0164 - Biometria da Interacao Genotipos x Ambientes
# Aula 3: Regressao conjunta para adaptabilidade e estabilidade
# Dados: soy_MET.txt | 40 genotipos | 6 ambientes | 3 blocos | DBC | 720 obs
#        CHK01-CHK04 = testemunhas | L001-L036 = linhagens candidatas
#        Variavel resposta: produtividade (kg/ha)
# Prof. Dr. Rafael Tassinari Resende
# =============================================================================

# install.packages(c("ggplot2", "reshape2", "RColorBrewer"))
library(ggplot2); library(reshape2); library(RColorBrewer)

dat <- read.table("soy_MET.txt", header = TRUE, sep = "\t")
dat$gen <- factor(dat$gen); dat$env <- factor(dat$env); dat$block <- factor(dat$block)

env_info <- data.frame(
  env      = levels(dat$env),   # usa os niveis reais do fator
  location = c("Anapolis","Luziania","Jatai","RioVerde","Uberlandia","Sorriso"),
  state    = c("GO","GO","GO","GO","MG","MT"), stringsAsFactors = FALSE)

fmt_p     <- function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p))
sig_stars <- function(p) ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.10,"."," "))))

# -----------------------------------------------------------------------------
# 1. DIMENSOES E DESCRITIVAS
# -----------------------------------------------------------------------------

g <- nlevels(dat$gen); a <- nlevels(dat$env); r <- nlevels(dat$block)
gen_chk <- grep("^CHK", levels(dat$gen), value = TRUE)  # testemunhas
gen_lin <- grep("^L",   levels(dat$gen), value = TRUE)  # linhagens candidatas

cat(sprintf("\n=== DIMENSOES ===\n  Genotipos: %d  (%d CHK testemunhas + %d L linhagens candidatas)\n",
    g, length(gen_chk), length(gen_lin)))
cat(sprintf("  Ambientes: %d | Blocos: %d | Obs: %d\n", a, r, nrow(dat)))

means_ge   <- tapply(dat$y, list(dat$gen, dat$env), mean)
env_means  <- colMeans(means_ge)
gen_means  <- rowMeans(means_ge)
grand_mean <- mean(means_ge)

cat("\nMedias por ambiente [kg/ha]:\n")
print(data.frame(Amb=names(env_means), Local=env_info$location,
                 Estado=env_info$state, Media=round(env_means,1)), row.names=FALSE)
cat(sprintf("\nMedias por genotipo (top 10):\n"))
print(round(sort(gen_means, decreasing=TRUE)[1:10], 1))
cat(sprintf("\nMedia geral: %.1f kg/ha | Amplitude: %.1f a %.1f\n",
    grand_mean, min(env_means), max(env_means)))

# Grafico de interacao G x A — CHK em vermelho, linhagens em cinza
mat_plot <- t(means_ge)
col_tipo <- ifelse(rownames(means_ge) %in% gen_chk, "tomato", "gray55")
lwd_tipo <- ifelse(rownames(means_ge) %in% gen_chk, 1.8, 0.7)

par(mar=c(5,5,4,2))
plot(NA, xlim=c(1,a), ylim=range(mat_plot)+c(-200,200), xaxt="n",
     xlab="Ambiente", ylab="Produtividade media (kg/ha)",
     main="Grafico de interacao G x A\nVermelho = testemunhas (CHK) | Cinza = linhagens (L)")
axis(1, at=1:a, labels=paste0(env_info$env,"\n(",env_info$location,")"), cex.axis=0.75)
abline(h=grand_mean, lty=2, col="black", lwd=1.5)
for(i in seq_len(nrow(means_ge)))
  lines(1:a, mat_plot[, rownames(means_ge)[i]],
        col=adjustcolor(col_tipo[i],0.55), lwd=lwd_tipo[i])
legend("topleft", bty="n", cex=0.8,
  legend=c("Media geral","Testemunhas (CHK)","Linhagens (L)"),
  lty=c(2,1,1), col=c("black","tomato","gray55"), lwd=c(1.5,1.8,0.7))

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

idx_tab <- data.frame(Ambiente=names(I_j), Local=env_info$location,
  Media=round(env_means,1), I_j=round(I_j,1),
  Tipo=ifelse(I_j>=0,"Favoravel (I>0)","Desfavoravel (I<0)"))
print(idx_tab, row.names=FALSE)
cat(sprintf("\nVerificacao: sum(I_j) = %.6f  (deve ser 0)\n", sum(I_j)))
cat("\nObs.: I_j e construido a partir das PROPRIAS medias do ensaio (indice endogeno).\n")
cat("      Ele nao mede clima/solo diretamente; resume o efeito agregado do ambiente\n")
cat("      sobre todos os genotipos avaliados. Ver Freeman & Perkins (Secao F).\n")

par(mar=c(6,5,4,2))
bp <- barplot(sort(I_j), col=ifelse(sort(I_j)>=0,"steelblue","tomato"),
  ylab="I_j (kg/ha)", main="Indice ambiental por ambiente",
  names.arg=paste0(names(sort(I_j)),"\n(",env_info$location[match(names(sort(I_j)),env_info$env)],")"),
  las=2, cex.names=0.7)
abline(h=0, lwd=1.5)
text(bp, sort(I_j)+sign(sort(I_j))*50, labels=round(sort(I_j),0), cex=0.8, font=2)

# =============================================================================
# SECAO B: FINLAY E WILKINSON (1963)
# Modelo: y_ij = mu_i + b_i * I_j + e_ij
# mu_i = media do genotipo (intercepto, pois sum(I_j) = 0)
# b_i  = coef. de regressao (responsividade ao gradiente ambiental)
# b_i ~ 1: adaptabilidade geral
# b_i > 1: alta responsividade (melhor em ambientes favoraveis)
# b_i < 1: baixa responsividade (menor perda em ambientes desfavoraveis)
# CHK e L diferenciados por pch (triangulo vs ponto) e cor nas tabelas/graficos
# =============================================================================

cat("\n\n=== SECAO B: FINLAY E WILKINSON (1963) ===\n")
cat("Modelo: y_ij = mu_i + b_i*I_j + e_ij\n")
cat("b_i = sum_j(y_ij*I_j) / sum_j(I_j^2)  [estimador de MQO]\n")
cat("Como sum(I_j)=0: intercepto = media do genotipo (mu_i = ybar_i.)\n\n")

sum_Ij2 <- sum(I_j^2)
b_hat   <- apply(means_ge, 1, function(y_i) sum(y_i*I_j)/sum_Ij2)
mu_hat  <- gen_means

# Perfil operacional (inferencia formal via teste t no metodo de E&R, Secao C)
perfil_fw <- ifelse(b_hat>1.10,"Alta responsividade (b>1)",
             ifelse(b_hat<0.90,"Baixa responsividade (b<1)","Adapt. geral (b~1)"))
tipo_gen  <- ifelse(names(b_hat) %in% gen_chk, "CHK", "L")

fw_tab <- data.frame(Tipo=tipo_gen, Genotipo=names(mu_hat),
  Media=round(mu_hat,1), b_i=round(b_hat,3), Perfil=perfil_fw)
fw_tab <- fw_tab[order(-fw_tab$Media),]; rownames(fw_tab) <- NULL
cat("Tabela F&W (por media desc. | CHK=testemunha | L=linhagem candidata):\n")
print(fw_tab, row.names=FALSE)
cat("\nLimites operacionais: b>1.10 = alta resp. | b<0.90 = baixa resp. | outros = geral\n")
cat("A inferencia formal (H0: b_i=1) e feita via teste t no metodo de Eberhart e Russell.\n")

# Retas de regressao
I_seq  <- seq(min(I_j)*1.15, max(I_j)*1.15, length.out=200)
col_fw <- ifelse(b_hat>1.10,"tomato3",ifelse(b_hat<0.90,"steelblue","gray60"))
pch_fw <- ifelse(names(b_hat)%in%gen_chk, 17, 16)
names(col_fw) <- names(b_hat)

par(mar=c(5,5,4,2))
plot(NA, xlim=range(I_seq), ylim=range(means_ge)+c(-300,300),
     xlab="I_j (kg/ha)", ylab="Produtividade media (kg/ha)",
     main="F&W (1963) | Retas de regressao por genotipo\nTriangulo=CHK | Ponto=linhagem | Cor=perfil de b_i")
abline(v=0, lty=2, col="gray40", lwd=1); abline(h=grand_mean, lty=2, col="gray40", lwd=1)
for(nm in names(b_hat)) {
  abline(a=mu_hat[nm], b=b_hat[nm], col=adjustcolor(col_fw[nm],0.4), lwd=0.8)
  points(I_j, means_ge[nm,], pch=pch_fw[nm],
         col=adjustcolor(col_fw[nm],0.6), cex=ifelse(nm%in%gen_chk,0.9,0.4))
}
top3 <- names(sort(b_hat,dec=TRUE))[1:3]; bot3 <- names(sort(b_hat))[1:3]
for(nm in c(top3,bot3)) {
  cor2 <- ifelse(nm%in%top3,"darkred","darkblue")
  abline(a=mu_hat[nm], b=b_hat[nm], col=cor2, lwd=2.5)
  text(max(I_j)*0.88, mu_hat[nm]+b_hat[nm]*max(I_j)*0.88,
       paste0(nm,"\n(b=",round(b_hat[nm],2),")"), col=cor2, cex=0.62, pos=4)
}
rug(I_j, side=1, col="gray30", lwd=2)
legend("topleft", bty="n", cex=0.72,
  legend=c("b>1","b<1","b~1","CHK (triangulo)","L (ponto)"),
  col=c("tomato3","steelblue","gray60","black","black"),
  lty=c(1,1,1,NA,NA), pch=c(NA,NA,NA,17,16), lwd=c(2,2,1,NA,NA))

# Quadrantes: media x b_i
par(mar=c(5,5,4,2))
plot(mu_hat, b_hat, pch=pch_fw, col=col_fw, cex=1.2,
     xlab="Media do genotipo (kg/ha)", ylab="b_i",
     main="F&W: Media x Adaptabilidade | Triangulo=CHK | Ponto=L")
abline(h=1, lty=2, col="gray40", lwd=1.5); abline(v=grand_mean, lty=2, col="gray40", lwd=1.5)
text(mu_hat, b_hat, labels=names(mu_hat), cex=0.57, pos=3, col=col_fw)
usr <- par("usr")
text(usr[2]-0.01*(usr[2]-usr[1]),usr[4]-0.01*(usr[4]-usr[3]),
     "Alta media\nAlta resp.", cex=0.65, adj=c(1,1), col="gray50")
text(usr[2]-0.01*(usr[2]-usr[1]),usr[3]+0.01*(usr[4]-usr[3]),
     "Alta media\nBaixa resp.", cex=0.65, adj=c(1,0), col="gray50")

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

# C1. ANOVA conjunta -> QMR
# Modelo: y ~ gen*env + env:block  (DBC dentro de cada ambiente)
# GL residuo = a*(r-1)*(g-1) = 6*2*39 = 468
cat("--- C1. ANOVA conjunta (para obter QMR) ---\n")
fit_joint <- aov(y ~ gen*env + env:block, data=dat)
tab_joint <- summary(fit_joint)[[1]]
cat("\nANOVA conjunta:\n"); print(round(tab_joint,3))

QMR <- tab_joint["Residuals","Mean Sq"]; nu <- tab_joint["Residuals","Df"]
sigma2_e <- QMR/r;  V_beta1 <- sigma2_e/sum_Ij2
cat(sprintf("\nQMR=%.2f | nu=%d | sigma2_e=QMR/r=%.2f | V(beta1)=%.6f\n",
    QMR, nu, sigma2_e, V_beta1))
cat(sprintf("GL esperado: a*(r-1)*(g-1) = %d*%d*%d = %d\n", a,r-1,g-1,a*(r-1)*(g-1)))

# C2. Estimacao dos parametros de E&R
# SQ(A/G)_i = r*sum_j(y_ij - ybar_i)^2
# SQ(reg)_i = r*(sum_j y_ij*I_j)^2 / sum_j I_j^2
# QMD_i     = (SQ(A/G)_i - SQ(reg)_i) / (a-2)
# sigma2_di = (QMD_i - QMR) / r          [valores negativos: sem efeito extra]
# H0: beta1_i=1  -> t_i = (beta1_i-1)/sqrt(V_beta1), gl=nu
# H0: sigma2_di=0 -> F_i = QMD_i/QMR, gl=(a-2, nu)
cat("\n--- C2. Estimacao dos parametros ---\n")

er <- data.frame(Tipo=ifelse(rownames(means_ge)%in%gen_chk,"CHK","L"),
  Genotipo=rownames(means_ge), Media=NA_real_, beta1=NA_real_,
  sigma2_di=NA_real_, R2=NA_real_, t_stat=NA_real_,
  p_t=NA_real_, F_stat=NA_real_, p_F=NA_real_, stringsAsFactors=FALSE)

for(i in seq_len(g)) {
  nm <- rownames(means_ge)[i]; y_i <- means_ge[i,]; mu_i <- mean(y_i)
  SQ_AG  <- r*sum((y_i-mu_i)^2)
  SQ_reg <- r*(sum(y_i*I_j))^2/sum_Ij2
  QMD    <- (SQ_AG-SQ_reg)/(a-2)
  t_i    <- (b_hat[nm]-1)/sqrt(V_beta1)
  er[i,"Media"]     <- round(mu_i,1)
  er[i,"beta1"]     <- round(b_hat[nm],3)
  er[i,"sigma2_di"] <- round((QMD-QMR)/r,1)
  er[i,"R2"]        <- round(SQ_reg/SQ_AG*100,1)
  er[i,"t_stat"]    <- round(t_i,3)
  er[i,"p_t"]       <- 2*pt(-abs(t_i), df=nu)
  er[i,"F_stat"]    <- round(QMD/QMR,3)
  er[i,"p_F"]       <- pf(QMD/QMR, df1=a-2, df2=nu, lower.tail=FALSE)
}
er <- er[order(-er$Media),]; rownames(er) <- NULL

# C3. Tabela de resultados
cat("\nTabela E&R (por media desc. | CHK=testemunha | L=linhagem candidata):\n\n")
cat(sprintf("%-4s %-8s %8s %7s %10s %5s  %7s %-9s  %7s %-9s\n",
    "Tipo","Gen","Media","beta1","sigma2_di","R2","t(b=1)","p(t)","F(s2=0)","p(F)"))
cat(strrep("-",90),"\n")
for(i in seq_len(nrow(er)))
  cat(sprintf("%-4s %-8s %8.1f %7.3f %10.1f %5.1f  %7.3f %-9s  %7.3f %-9s\n",
    er$Tipo[i],er$Genotipo[i],er$Media[i],er$beta1[i],er$sigma2_di[i],er$R2[i],
    er$t_stat[i], paste0(fmt_p(er$p_t[i]),sig_stars(er$p_t[i])),
    er$F_stat[i], paste0(fmt_p(er$p_F[i]),sig_stars(er$p_F[i]))))
cat(strrep("-",90),"\n")
cat("t(b=1): H0 beta1=1 (adapt. especifica se rejeitado)\n")
cat("F(s2=0): H0 sigma2_di=0 (comportamento imprevisto se rejeitado)\n")
cat("*** p<0.001  ** p<0.01  * p<0.05  . p<0.10\n")

# C4. Classificacao conjunta
er$Classe <- with(er, sapply(seq_len(nrow(er)), function(i) {
  if(Media[i]<grand_mean) return("Baixa media - pouco atrativo")
  prev <- p_F[i]>=0.05; fav <- p_t[i]<0.05 & beta1[i]>1; desf <- p_t[i]<0.05 & beta1[i]<1
  if(!prev & fav)  return("Alta media | adapt. favoraveis | pouco previsivel")
  if(!prev & desf) return("Alta media | adapt. desfavoraveis | pouco previsivel")
  if(!prev)        return("Alta media | adapt. geral | pouco previsivel")
  if(fav)          return("Alta media | adapt. favoraveis | previsivel")
  if(desf)         return("Alta media | adapt. desfavoraveis | previsivel")
  return("Alta media | adapt. geral | previsivel")
}))
ideais <- er$Genotipo[er$Classe == "Alta media | adapt. geral | previsivel"]
cat("\nGenotipos: alta media + adapt. geral + comportamento previsivel:\n")
cat(" ", paste(ideais, collapse=", "), "\n\n")
cat("Contagem por classe:\n"); print(table(er$Classe))

# C5. Retas coloridas pela classificacao
cls_col <- c("Alta media | adapt. geral | previsivel"="darkgreen",
  "Alta media | adapt. favoraveis | previsivel"="tomato3",
  "Alta media | adapt. desfavoraveis | previsivel"="steelblue",
  "Alta media | adapt. geral | pouco previsivel"="darkorange",
  "Alta media | adapt. favoraveis | pouco previsivel"="pink3",
  "Alta media | adapt. desfavoraveis | pouco previsivel"="lightblue3",
  "Baixa media - pouco atrativo"="gray75")
col_er <- cls_col[er$Classe]; names(col_er) <- er$Genotipo
pch_er <- ifelse(er$Genotipo%in%gen_chk, 17, 16); names(pch_er) <- er$Genotipo

par(mar=c(5,5,4,2))
plot(NA, xlim=range(I_seq), ylim=range(means_ge)+c(-300,300),
     xlab="I_j (kg/ha)", ylab="Produtividade media (kg/ha)",
     main="E&R (1966) | Classificacao conjunta\nTriangulo=CHK | Ponto=linhagem | Cor=classe")
abline(v=0,lty=2,col="gray50"); abline(h=grand_mean,lty=2,col="gray50")
for(nm in er$Genotipo) {
  abline(a=mu_hat[nm],b=b_hat[nm],col=adjustcolor(col_er[nm],0.45),lwd=0.9)
  points(I_j,means_ge[nm,],pch=pch_er[nm],
         col=adjustcolor(col_er[nm],0.6),cex=ifelse(nm%in%gen_chk,0.9,0.35))
}
for(nm in ideais) abline(a=mu_hat[nm],b=b_hat[nm],col="darkgreen",lwd=2.5)
rug(I_j,side=1,col="gray30",lwd=2)
legend("topleft",bty="n",cex=0.66,
  legend=c("geral|prev (IDEAL)","fav|prev","desf|prev","*|imprev","Baixa media","CHK","L"),
  col=c("darkgreen","tomato3","steelblue","darkorange","gray75","black","black"),
  lty=c(1,1,1,1,1,NA,NA),pch=c(NA,NA,NA,NA,NA,17,16),lwd=2)

# C6. Inspecao individual: melhor e pior R2 entre genotipos de alta media
er_alta <- er[er$Media>=grand_mean,]
nm_melhor <- er_alta$Genotipo[which.max(er_alta$R2)]
nm_pior   <- er_alta$Genotipo[which.min(er_alta$R2)]
cat(sprintf("\n--- C6. Inspecao visual individual ---\n"))
cat(sprintf("Maior R2 (mais previsivel): %s [%s]  R2=%.1f%%\n",
    nm_melhor, ifelse(nm_melhor%in%gen_chk,"CHK","L"),
    er_alta$R2[er_alta$Genotipo==nm_melhor]))
cat(sprintf("Menor R2 (menos previsivel): %s [%s]  R2=%.1f%%\n",
    nm_pior, ifelse(nm_pior%in%gen_chk,"CHK","L"),
    er_alta$R2[er_alta$Genotipo==nm_pior]))

par(mfrow=c(1,2), mar=c(5,5,4,2))
for(nm in c(nm_melhor,nm_pior)) {
  rw <- er[er$Genotipo==nm,]; y_o <- means_ge[nm,]; y_f <- mu_hat[nm]+b_hat[nm]*I_j
  plot(I_j,y_o, xlim=range(I_j)*1.1, ylim=range(c(y_o,y_f))+c(-200,200),
       xlab="I_j (kg/ha)", ylab="Produtividade (kg/ha)", pch=16, col="steelblue", cex=1.4,
       main=paste0(nm," [",ifelse(nm%in%gen_chk,"CHK","L"),"]\n",
         "b=",round(rw$beta1,2)," | s2d=",round(rw$sigma2_di,0)," | R2=",round(rw$R2,1),"%"))
  abline(a=mu_hat[nm],b=b_hat[nm],col="tomato",lwd=2)
  abline(v=0,lty=2,col="gray50"); abline(h=mu_hat[nm],lty=3,col="gray50")
  for(j in seq_along(I_j)) segments(I_j[j],y_f[j],I_j[j],y_o[j],col="gray60",lwd=1)
  text(I_j,y_o,labels=names(I_j),pos=3,cex=0.75,col="gray30")
  legend("topleft",bty="n",cex=0.72,legend=c("Obs.","Reta","Desvio"),
    pch=c(16,NA,NA),lty=c(NA,1,1),col=c("steelblue","tomato","gray60"),lwd=c(NA,2,1))
}
par(mfrow=c(1,1))
cat("\nR2 alto: reta captura bem a resposta | R2 baixo: inspecionar grafico.\n")
cat("sigma2_di pode ser significativo mesmo com R2 alto (desvio sistematico).\n")
cat("A combinacao R2 + sigma2_di + inspecao grafica e mais informativa que cada um isolado.\n")

# C7. R2 por genotipo — CHK destacados em cor diferente
R2_vec <- er$R2; names(R2_vec) <- er$Genotipo
par(mar=c(5,7,4,2))
barplot(sort(R2_vec), horiz=TRUE, las=1, xlab="R2 (%)",
  main="R2_i por genotipo (E&R)\nVermelho=CHK | Azul=L R2>=70% | Cinza=L R2<70%",
  col=ifelse(names(sort(R2_vec))%in%gen_chk,"tomato3",
      ifelse(sort(R2_vec)>=70,"steelblue","gray65")),
  cex.names=0.62, xlim=c(0,105))
abline(v=70,lty=2,col="gray40",lwd=1.5)
text(72,1,"R2=70%",cex=0.7,col="gray40",adj=0)

# =============================================================================
# SECAO D: DECOMPOSICAO DA SQ (Eberhart & Russell)
# SQ(A/G) = SQ(A linear) + SQ(GA linear) + SQDC
# A linear : 1 gl       — gradiente ambiental medio (b_A = 1 sempre)
# GA linear: g-1 gl     — heterogeneidade das inclinacoes entre genotipos
# SQDC     : g*(a-2) gl — desvios nao explicados pela componente linear
# =============================================================================

cat("\n\n=== SECAO D: DECOMPOSICAO DA SQ (E&R) ===\n")

SQ_AG_tot  <- r*sum(apply(means_ge,1,function(y) sum((y-mean(y))^2)))
SQ_reg_all <- sum(sapply(seq_len(g),function(i) r*(sum(means_ge[i,]*I_j))^2/sum_Ij2))
SQ_A_lin   <- g*r*sum_Ij2                    # b_A = 1 sempre (propriedade do indice)
SQ_GA_lin  <- SQ_reg_all - SQ_A_lin
SQDC       <- SQ_AG_tot  - SQ_reg_all
QM_GA_lin  <- SQ_GA_lin/(g-1); QM_DC <- SQDC/(g*(a-2))
F_GA <- QM_GA_lin/QMR; p_GA <- pf(F_GA,df1=g-1,   df2=nu,lower.tail=FALSE)
F_DC <- QM_DC/QMR;     p_DC <- pf(F_DC,df1=g*(a-2),df2=nu,lower.tail=FALSE)

cat(sprintf("Verificacao: SQ(A/G) = SQ_reg + SQDC -> %.1f = %.1f + %.1f (dif=%.2f)\n",
    SQ_AG_tot, SQ_reg_all, SQDC, SQ_AG_tot-SQ_reg_all-SQDC))
cat(sprintf("\n%-28s %6s %15s %12s %8s %10s\n","FV","GL","SQ","QM","F","p"))
cat(strrep("-",82),"\n")
linhas <- list(
  list("A linear",           1,       SQ_A_lin,  SQ_A_lin,   NA,   NA),
  list("GA linear",          g-1,     SQ_GA_lin, QM_GA_lin,  F_GA, p_GA),
  list("Desvio comb. (A/G)", g*(a-2), SQDC,      QM_DC,      F_DC, p_DC),
  list("  Desvio/Gi (media)",a-2,     SQDC/g,    SQDC/g/(a-2),NA,  NA),
  list("Residuo (ANOVA)",    nu,      QMR*nu,    QMR,         NA,   NA))
for(ln in linhas) {
  f_s <- if(is.na(ln[[5]])) "       --" else sprintf("%8.3f",ln[[5]])
  p_s <- if(is.na(ln[[6]])) "         --" else
         paste0(sprintf("%10s",fmt_p(ln[[6]])),sig_stars(ln[[6]]))
  cat(sprintf("%-28s %6d %15.1f %12.1f %s %s\n",
      ln[[1]],ln[[2]],ln[[3]],ln[[4]],f_s,p_s))
}
cat(strrep("-",82),"\n")
cat("GA linear significativo: genotipos diferem em adaptabilidade (inclinacoes heterogeneas).\n")
cat("Desvio combinado signif.: ha componente nao linear na interacao nao capturada pela reta.\n")

# =============================================================================
# SECAO E: QUESTOES PRATICAS E CAUTELAS
# =============================================================================

cat("\n\n=== SECAO E: QUESTOES PRATICAS E CAUTELAS ===\n")

# E1. Linearidade
cat(sprintf("\n[E1] Linearidade (R2_i):\n"))
cat(sprintf("     R2 medio=%.1f%% | min=%.1f%% (%s) | max=%.1f%% (%s)\n",
    mean(er$R2),min(er$R2),er$Genotipo[which.min(er$R2)],
    max(er$R2),er$Genotipo[which.max(er$R2)]))
bad_r2 <- sort(R2_vec[R2_vec<70])
cat("     Genotipos com R2<70% (ajuste linear questionavel):\n")
if(length(bad_r2)>0) print(bad_r2) else cat("     Nenhum.\n")

# E2. Alavancagem ambiental (I_j recalculado apos remocao)
env_rem   <- names(which.max(abs(I_j)))
means_red <- means_ge[,colnames(means_ge)!=env_rem]
I_j_red   <- colMeans(means_red)-mean(means_red)  # recalculo obrigatorio
b_red     <- apply(means_red,1,function(y) sum(y*I_j_red)/sum(I_j_red^2))
delta_b   <- b_red - b_hat[names(b_red)]
cat(sprintf("\n[E2] Alavancagem: amb. mais extremo = %s (%s) I_j=%.1f\n",
    env_rem, env_info$location[env_info$env==env_rem], I_j[env_rem]))
cat("     I_j recalculado apos remocao (escala corrigida) antes de reestimar b_i.\n")
cat(sprintf("     Delta max: %.3f (%s) | Delta medio: %.3f\n",
    max(abs(delta_b)),names(which.max(abs(delta_b))),mean(abs(delta_b))))
cat("     Inspecao grafica individual recomendada para genotipos com R2<70%.\n")

# E3. Dependencia do indice
cat(sprintf("\n[E3] Dependencia do I_j: cada genotipo contribui 1/%d=%.3f do proprio preditor.\n",g,1/g))
cat("     b_i e propriedade RELATIVA ao experimento (painel + ambientes + escala).\n")
cat("     Comparar b_i entre redes/anos distintos exige cautela.\n")

# E4. Numero e amplitude
cat(sprintf("\n[E4] a=%d ambientes -> %d gl para desvio por genotipo.\n",a,a-2))
cat(sprintf("     Amplitude de I_j: %.1f kg/ha. Quanto maior, mais informativo.\n",diff(range(I_j))))
cat("     Situacoes com a<5 geram estimativas de b_i e sigma2_di pouco estaveis.\n")

# E5. Leitura conjunta dos tres parametros
cat("\n[E5] Os tres parametros de E&R formam um TRIPE; nenhum substitui os demais:\n")
cat("     Media  -> nivel de producao (criterio prioritario de selecao)\n")
cat("     sigma2_di -> previsibilidade (criterio de descarte por inconsistencia)\n")
cat("     beta1  -> perfil de adaptacao (orienta zona de recomendacao)\n")
cat("     Ordem pratica: media > previsibilidade > perfil de adaptacao.\n")

# =============================================================================
# SECAO F: PALHINHA — FREEMAN E PERKINS (1971)
# Problema: I_j endogeno nao e independente das respostas que se pretende modelar.
# Solucao proposta: usar genotipos-testemunha (CHK) para construir I_j exogeno.
# Resultado pratico aqui: com g=40 e 4 CHK fixas em todos os ambientes, a
# correlacao entre os dois indices e entre os dois conjuntos de b_i revela
# quao critico (ou nao) o problema e neste ensaio.
# Vantagem adicional das CHK fixas: permitem comparar b_i entre redes/anos
# na mesma escala de referencia ambiental.
# =============================================================================

cat("\n\n=== SECAO F: PALHINHA — FREEMAN & PERKINS (1971) ===\n")
cat("Problema: I_j classico e endogeno — construido a partir das proprias\n")
cat("respostas fenotipicas, incluindo as das linhagens que se quer avaliar.\n")
cat("Solucao: I_j exogeno, calculado apenas com as testemunhas (CHK),\n")
cat("cujas respostas sao independentes das linhagens candidatas.\n\n")

means_chk <- means_ge[gen_chk,]
I_j_exo   <- colMeans(means_chk) - mean(means_chk)

cat("Comparacao dos indices ambientais (endogeno vs. exogeno):\n")
print(data.frame(Amb=names(I_j), Local=env_info$location,
  I_endo=round(I_j,1), I_exo=round(I_j_exo,1),
  Dif=round(I_j_exo-I_j,1)), row.names=FALSE)
cat(sprintf("\nCorrelacao I_endo x I_exo: r = %.4f\n", cor(I_j,I_j_exo)))
cat("Alta correlacao: testemunhas e linhagens respondem de forma similar ao gradiente.\n\n")

# b_i das linhagens sobre o indice exogeno
b_exo      <- apply(means_ge[gen_lin,],1,function(y) sum(y*I_j_exo)/sum(I_j_exo^2))
b_endo_lin <- b_hat[gen_lin]
cat(sprintf("Correlacao b_endo x b_exo (linhagens): r = %.4f\n", cor(b_endo_lin,b_exo)))
cat("Alta correlacao confirma que o vies de atenuacao do indice endogeno e\n")
cat(sprintf("negligivel com g=%d genotipos (cada um contribui apenas 1/%d=%.3f).\n",g,g,1/g))

par(mar=c(5,5,4,2))
plot(b_endo_lin, b_exo, pch=16, col="steelblue", cex=1.2,
     xlab="b_i endogeno (todos os genotipos)",
     ylab="b_i exogeno (apenas testemunhas CHK)",
     main="Freeman & Perkins (1971) | Linhagens candidatas\nb_i endogeno vs. exogeno")
abline(0,1,lty=2,col="gray40",lwd=1.5)
abline(h=1,lty=3,col="tomato",lwd=1); abline(v=1,lty=3,col="tomato",lwd=1)
text(b_endo_lin,b_exo,labels=gen_lin,cex=0.58,pos=3,col="gray30")
legend("topleft",bty="n",cex=0.8,
  legend=c("Linhagens","Igualdade perfeita","b=1"),
  pch=c(16,NA,NA),lty=c(NA,2,3),col=c("steelblue","gray40","tomato"))

cat("\nMensagem final desta secao:\n")
cat("  O indice exogeno (CHK) e conceitualmente superior ao endogeno.\n")
cat("  Na pratica, com muitos genotipos, os dois tendem a coincidir.\n")
cat("  A vantagem operacional das CHK fixas vai alem da teoria:\n")
cat("  permite comparar b_i entre redes de ensaios e safras diferentes\n")
cat("  na mesma escala de referencia ambiental — algo impossivel com\n")
cat("  o indice endogeno, que varia com o painel avaliado.\n")

# =============================================================================
# RESUMO FINAL
# =============================================================================

cat("\n\n=============================================================\n")
cat("RESUMO FINAL - AULA 3\n")
cat("=============================================================\n")
cat(sprintf("Genotipos: %d (%d CHK + %d L) | Ambientes: %d | Blocos: %d\n",
    g,length(gen_chk),length(gen_lin),a,r))
cat(sprintf("Media geral: %.1f kg/ha | QMR: %.2f (gl=%d)\n", grand_mean,QMR,nu))
cat(sprintf("I_j: %.1f a +%.1f kg/ha (amplitude=%.1f)\n",
    min(I_j),max(I_j),diff(range(I_j))))
cat(sprintf("b_i: min=%.2f (%s) | max=%.2f (%s) | media=%.2f\n",
    min(b_hat),names(which.min(b_hat)),max(b_hat),names(which.max(b_hat)),mean(b_hat)))
cat("\nGenotipos ideais (alta media | adapt. geral | previsivel):\n")
cat(" ", paste(ideais, collapse=", "), "\n")
cat("\nContagem por classe:\n"); print(table(er$Classe))
cat(sprintf("\nGA linear: F=%.2f p=%s | Desvio comb: F=%.2f p=%s\n",
    F_GA,fmt_p(p_GA),F_DC,fmt_p(p_DC)))
cat(sprintf("I_j endo x exo (CHK): r=%.4f | b_i endo x exo (L): r=%.4f\n",
    cor(I_j,I_j_exo), cor(b_endo_lin,b_exo)))
cat("\nCautelas:\n")
cat("  1. Tripe E&R: media > previsibilidade > perfil de adaptacao.\n")
cat("  2. b_i relativo ao painel — comparar entre redes/anos com cautela.\n")
cat("  3. R2_i: primeiro diagnostico de nao-linearidade.\n")
cat("  4. Ambientes extremos exercem alta alavancagem — inspecionar graficos.\n")
cat("  5. CHK fixas: referencia exogena para I_j e comparabilidade entre safras.\n")
cat("=============================================================\n")
