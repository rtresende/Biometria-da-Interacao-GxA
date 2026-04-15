# =============================================================================
#  MAIZE SAFRINHA MULTI-ENVIRONMENT TRIAL SIMULATOR
#  Discipline: Biometria da Interação Genótipos × Ambientes — GMP0164/UFG 2026
#  Prof. Rafael Tassinari Resende
# =============================================================================
#  Crossing design : NC Design II (Comstock & Robinson, 1952)
#                    20 females (F01–F20) × 5 males (M01–M05) → 62 hybrids
#  Checks          : 8 commercial cultivars (CHK01–CHK08)
#  Environments    : 6 safrinha locations — Cerrado biome, 2025
#  Unbalancing     : strategic, between environments only
#                    checks always complete (6/6); experimentals stratified
#  Response        : grain yield (kg/ha)
#  Lessons         : 6 (Mixed Models I: BLUP), 7, 8, 9, 10
# =============================================================================
#  Dependencies: dplyr
# =============================================================================

setwd("G:/Meu Drive/UFG/PPGGMP/BiometriaDaInteracaoGxA/aulas/simu_data")

library(dplyr)
set.seed(2026)


# =============================================================================
# SECTION 1 — CONFIGURATION
# =============================================================================

MU       <- 6000   # grand mean (kg/ha)
SD_ENV   <- 1500   # environment main effect SD
SD_GEN   <-  500   # genotype main effect SD
SD_GXE   <-  700   # GxE interaction SD
SD_BLOCK <-  200   # block-within-environment SD
SD_ERROR <-  650   # residual SD

N_BLOCKS <- 3

ENV <- data.frame(
  trial   = c("TRL_CPS_2025", "TRL_RVD_2025", "TRL_JAT_2025",
              "TRL_SOR_2025", "TRL_SIN_2025", "TRL_DOU_2025"),
  city    = c("Chapadão do Sul", "Rio Verde", "Jataí",
              "Sorriso", "Sinop", "Dourados"),
  state   = c("MS", "GO", "GO", "MT", "MT", "MS"),
  env_eff = c(1400, 900, 650, 50, -600, -1200),
  stringsAsFactors = FALSE
)

n_envs  <- nrow(ENV)
env_ids <- ENV$trial


# =============================================================================
# SECTION 2 — GENOTYPE METADATA
# =============================================================================

# --- NC Design II crossing plan ---
# Elite testers (M01, M02): more crosses; M05: fewest
cross_plan <- rbind(
  data.frame(female = sprintf("F%02d", 1:15),  male = "M01"),
  data.frame(female = sprintf("F%02d", 4:18),  male = "M02"),
  data.frame(female = sprintf("F%02d", 6:18),  male = "M03"),
  data.frame(female = sprintf("F%02d", 8:18),  male = "M04"),
  data.frame(female = sprintf("F%02d", 13:20), male = "M05")
)
cross_plan$gen <- sprintf("H%03d", seq_len(nrow(cross_plan)))  # H001-H062

# GxE adaptation profile (internal -- not exported in final dataset)
gxe_profiles <- c(rep("stable",15), rep("responsive",15), rep("rustic",12),
                  rep("unstable",12), rep("inferior",8))
cross_plan$gxe_profile <- sample(gxe_profiles)

# --- Parental maturity (days to physiological maturity) ---
# Each parental line has ONE fixed maturity value -- intrinsic genetic characteristic.
# Females (linhagens candidatas): broader range, more diversity among candidates
# Males (testadores): narrower range -- elite, well-characterized lines
dtm_f <- setNames(sample(98:122, 20, replace = TRUE), sprintf("F%02d", 1:20))
dtm_m <- setNames(sample(103:118,  5, replace = TRUE), sprintf("M%02d", 1:5))

# Hybrid maturity = midparent + small SCA deviation (sd ~ 2 days)
# Clamped to biologically plausible range for safrinha hybrids [100, 125]
cross_plan$maturity <- pmin(pmax(
  round((dtm_f[cross_plan$female] + dtm_m[cross_plan$male]) / 2 +
          rnorm(nrow(cross_plan), 0, 2)),
  100L), 125L)

cross_plan$gen_type       <- "experimental"
cross_plan$cultivar_name  <- NA_character_

# --- Commercial checks (pedigree unknown) ---
checks <- data.frame(
  gen           = sprintf("CHK%02d", 1:8),
  female        = NA_character_,
  male          = NA_character_,
  gxe_profile   = "check",
  maturity      = c(112L, 118L, 107L, 120L, 115L, 110L, 113L, 116L),
  gen_type      = "check",
  cultivar_name = c("BRS 4510 VIP3", "DKX 389 PRO4", "AGX 917 VIP3", "BRS 4071 PWR",
                    "NVX 632 PRO4",  "SYX 284 VIP3", "PIX 753 PRO4", "CRX 128 VIP3"),
  stringsAsFactors = FALSE
)

all_gen <- bind_rows(cross_plan, checks)
gen_ids <- all_gen$gen
n_gen   <- nrow(all_gen)   # 70

cat(sprintf("\n--- Genotypes: %d (%d experimentals + %d checks) ---\n",
            n_gen, nrow(cross_plan), nrow(checks)))
cat("\n--- Parental maturity (days) ---\n")
cat("  Females:", paste(sprintf("F%02d=%d", 1:20, dtm_f), collapse = ", "), "\n")
cat("  Males  :", paste(sprintf("M%02d=%d", 1:5,  dtm_m), collapse = ", "), "\n")


# =============================================================================
# SECTION 3 — GENOTYPIC EFFECTS (NC Design II: GCA + SCA)
# =============================================================================

# GCA components -- additive effects within each pool
gca_f <- setNames(rnorm(20, 0, SD_GEN * 0.60), sprintf("F%02d", 1:20))
gca_m <- setNames(rnorm(5,  0, SD_GEN * 0.40), sprintf("M%02d", 1:5))

gen_eff <- setNames(numeric(n_gen), gen_ids)

for (i in seq_len(n_gen)) {
  gid <- gen_ids[i]
  if (all_gen$gen_type[i] == "experimental") {
    sca          <- rnorm(1, 0, SD_GEN * 0.45)
    gen_eff[gid] <- gca_f[all_gen$female[i]] + gca_m[all_gen$male[i]] + sca
  } else {
    gen_eff[gid] <- rnorm(1, 300, 100)    # checks: slightly above-average
  }
}

# Rescale experimentals to target SD_GEN
exp_mask <- all_gen$gen_type == "experimental"
ge       <- gen_eff[exp_mask]
gen_eff[exp_mask] <- (ge - mean(ge)) / sd(ge) * SD_GEN

# Pull inferior group below average
inf_ids <- all_gen$gen[all_gen$gxe_profile == "inferior"]
gen_eff[inf_ids] <- gen_eff[inf_ids] - 1.2 * SD_GEN


# =============================================================================
# SECTION 4 — GxE INTERACTION MATRIX
# =============================================================================

env_scaled <- as.numeric(scale(ENV$env_eff))

GXE_MULT <- c(check = 0.40, stable = 0.50, responsive = 1.20,
              rustic = 1.00, unstable = 2.20, inferior = 0.60)

gxe_mat <- matrix(0, nrow = n_gen, ncol = n_envs,
                  dimnames = list(gen_ids, env_ids))

for (i in seq_len(n_gen)) {
  gid <- gen_ids[i]
  grp <- all_gen$gxe_profile[i]
  sig <- SD_GXE * GXE_MULT[grp]

  if (grp == "responsive") {
    gxe_mat[gid, ] <- env_scaled * runif(1, 0.40, 0.80) * sig +
                      rnorm(n_envs, 0, sig * 0.35)
  } else if (grp == "rustic") {
    gxe_mat[gid, ] <- -env_scaled * runif(1, 0.30, 0.65) * sig +
                      rnorm(n_envs, 0, sig * 0.35)
  } else {
    gxe_mat[gid, ] <- rnorm(n_envs, 0, sig)
  }
}

# Double-center (AMMI convention)
gxe_mat <- sweep(gxe_mat, 1, rowMeans(gxe_mat), "-")
gxe_mat <- sweep(gxe_mat, 2, colMeans(gxe_mat), "-")


# =============================================================================
# SECTION 5 — STRATEGIC UNBALANCING (between environments only)
# =============================================================================

# Checks      : all 6 environments (always complete)
# Elite   (15): 6/6 environments -- most advanced materials
# Advanced(20): 4-5 environments
# Intermediate(17): 3-4 environments
# Novice  (10): 2-3 environments -- first cycle in network

exp_ids    <- all_gen$gen[all_gen$gen_type == "experimental"]
adv_labels <- c(rep("elite", 15), rep("advanced", 20),
                rep("intermediate", 17), rep("novice", 10))
adv_of     <- setNames(sample(adv_labels), exp_ids)

env_present <- list()

for (g in checks$gen) env_present[[g]] <- env_ids

for (g in exp_ids) {
  n <- switch(adv_of[g],
    elite        = 6L,
    advanced     = sample(4:5, 1L),
    intermediate = sample(3:4, 1L),
    novice       = sample(2:3, 1L)
  )
  env_present[[g]] <- sort(sample(env_ids, n))
}


# =============================================================================
# SECTION 6 — ASSEMBLE DATASET
# =============================================================================

blk_keys <- paste(rep(env_ids, each = N_BLOCKS),
                  paste0("B", 1:N_BLOCKS), sep = ":")
blk_eff  <- setNames(rnorm(length(blk_keys), 0, SD_BLOCK), blk_keys)

rows <- list()

for (i in seq_len(n_gen)) {
  gid   <- gen_ids[i]
  ginfo <- all_gen[i, ]

  for (eid in env_present[[gid]]) {
    erow <- ENV[ENV$trial == eid, ]

    for (b in paste0("B", 1:N_BLOCKS)) {
      y <- MU +
           erow$env_eff +
           gen_eff[gid] +
           gxe_mat[gid, eid] +
           blk_eff[paste(eid, b, sep = ":")] +
           rnorm(1, 0, SD_ERROR)

      rows[[length(rows) + 1L]] <- data.frame(
        gen           = gid,
        female        = ginfo$female,
        male          = ginfo$male,
        gen_type      = ginfo$gen_type,
        cultivar_name = ginfo$cultivar_name,
        maturity      = ginfo$maturity,
        city          = erow$city,
        state         = erow$state,
        trial         = eid,
        block         = b,
        yield         = round(y, 1),
        stringsAsFactors = FALSE
      )
    }
  }
}

dat <- bind_rows(rows)


# =============================================================================
# SECTION 7 — VARIANCE COMPONENTS & SUMMARY
# =============================================================================

vG <- SD_GEN^2; vGxE <- SD_GXE^2; vB <- SD_BLOCK^2; vE <- SD_ERROR^2
vP <- vG + vGxE + vB + vE
mean_e <- mean(sapply(exp_ids, function(g) length(env_present[[g]])))

cat("\n--- Variance components ---\n")
cat(sprintf("  V_G   = %7.0f  (%.1f%%)\n", vG,   100 * vG   / vP))
cat(sprintf("  V_GxE = %7.0f  (%.1f%%)\n", vGxE, 100 * vGxE / vP))
cat(sprintf("  V_Blk = %7.0f  (%.1f%%)\n", vB,   100 * vB   / vP))
cat(sprintf("  V_Err = %7.0f  (%.1f%%)\n", vE,   100 * vE   / vP))
cat(sprintf("  h2 (plot level)              : %.2f\n", vG / vP))
cat(sprintf("  H2 (means, %.1f envs x %dr) : %.2f\n", mean_e, N_BLOCKS,
            vG / (vG + vGxE / mean_e + vE / (mean_e * N_BLOCKS))))

n_balanced <- n_gen * n_envs * N_BLOCKS
cat("\n--- Dataset ---\n")
cat(sprintf("  Balanced obs (if complete) : %d\n", n_balanced))
cat(sprintf("  Final obs                  : %d\n", nrow(dat)))
cat(sprintf("  Unbalancing                : %.1f%% missing\n",
            100 * (1 - nrow(dat) / n_balanced)))
cat(sprintf("  yield range   : %.0f - %.0f kg/ha\n", min(dat$yield), max(dat$yield)))
cat(sprintf("  yield mean    : %.0f kg/ha\n",         mean(dat$yield)))
cat(sprintf("  maturity range: %d - %d days\n",       min(dat$maturity), max(dat$maturity)))

cat("\n--- Environments per genotype (experimentals) ---\n")
print(table(sapply(exp_ids, function(g) length(env_present[[g]]))))


# =============================================================================
# SECTION 8 — EXPORT
# =============================================================================

OUT_FILE <- "maize_safrinha_MET.txt"
write.table(dat, file = OUT_FILE, row.names = FALSE, quote = FALSE, sep = "\t")
cat(sprintf("\n  Exported -> %s\n\n", OUT_FILE))


# =============================================================================
# SECTION 9 — SANITY CHECK (optional -- uncomment to run)
# =============================================================================
 library(lme4)
 mod_yield <- lmer(yield ~ -1 + trial + (1|gen) + (1|gen:trial) + (1|trial:block),
                  data = dat, REML = TRUE)
 print(VarCorr(mod_yield))

 mod_matu <- lmer(maturity ~ 1 + (1|gen), data = dat, REML = TRUE)
 print(VarCorr(mod_matu))   # only gen-level variance expected; no GxE, no residual

 plot(dat$maturity, dat$yield)
 abline(lm(yield ~ maturity, dat), col = "blue", lty = 2)
 summary(lm(yield ~ maturity, dat))
# =============================================================================











# =============================================================================
#  DIAGNOSTIC SCRIPT — maize_safrinha_MET.txt
#  Discipline: GMP0164/UFG 2026 — Prof. Rafael Tassinari Resende
# =============================================================================

library(dplyr)
library(ggplot2)

dat <- read.table("maize_safrinha_MET.txt", header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE)


# =============================================================================
# SECTION 1 — STRUCTURE
# =============================================================================

cat("\n========== SECTION 1: STRUCTURE ==========\n")

cat(sprintf("  Rows              : %d\n", nrow(dat)))
cat(sprintf("  Columns           : %d\n", ncol(dat)))
cat(sprintf("  Genotypes (total) : %d\n", n_distinct(dat$gen)))
cat(sprintf("    experimentals   : %d\n", n_distinct(dat$gen[dat$gen_type == "experimental"])))
cat(sprintf("    checks          : %d\n", n_distinct(dat$gen[dat$gen_type == "check"])))
cat(sprintf("  Environments      : %d\n", n_distinct(dat$trial)))
cat(sprintf("  Blocks per env    : %s\n", paste(sort(unique(dat$block)), collapse = ", ")))
cat(sprintf("  Missing values    : %d\n", sum(is.na(dat$yield_kg_ha))))

cat("\n  Column classes:\n")
for (col in names(dat)) cat(sprintf("    %-20s %s\n", col, class(dat[[col]])))


# =============================================================================
# SECTION 2 — CROSSING DESIGN INTEGRITY (NC Design II)
# =============================================================================

cat("\n========== SECTION 2: NC DESIGN II INTEGRITY ==========\n")

exp <- dat %>% filter(gen_type == "experimental") %>%
  select(gen, female, male) %>% distinct()

cat(sprintf("  Unique hybrids          : %d  (expected 62)\n", nrow(exp)))
cat(sprintf("  Unique females          : %d  (expected 20)\n", n_distinct(exp$female)))
cat(sprintf("  Unique males            : %d  (expected  5)\n", n_distinct(exp$male)))

# No female used as male and vice-versa
overlap <- intersect(unique(exp$female), unique(exp$male))
cat(sprintf("  Female/male pool overlap: %d  (expected  0)\n", length(overlap)))

# Crosses per tester
cat("\n  Crosses per tester:\n")
exp %>% count(male, name = "n_crosses") %>%
  arrange(male) %>%
  print(row.names = FALSE)

# Half-sib counts
cat("\n  Half-sibs per female (maternal):\n")
hs_f <- exp %>% count(female, name = "n_testers") %>% arrange(female)
print(summary(hs_f$n_testers))

cat("\n  Half-sibs per male (paternal):\n")
hs_m <- exp %>% count(male, name = "n_females") %>% arrange(male)
print(hs_m, row.names = FALSE)

# Checks: pedigree NA
chk <- dat %>% filter(gen_type == "check") %>%
  select(gen, female, male, cultivar_name) %>% distinct()
cat(sprintf("\n  Checks with NA pedigree : %d  (expected 8)\n",
            sum(is.na(chk$female) & is.na(chk$male))))
cat("  Check names:\n")
print(chk[, c("gen","cultivar_name")], row.names = FALSE)


# =============================================================================
# SECTION 3 — UNBALANCING DIAGNOSTICS
# =============================================================================

cat("\n========== SECTION 3: UNBALANCING ==========\n")

n_balanced <- 70 * 6 * 3
cat(sprintf("  Balanced obs (expected)  : %d\n", n_balanced))
cat(sprintf("  Observed obs             : %d\n", nrow(dat)))
cat(sprintf("  Missing                  : %d  (%.1f%%)\n",
            n_balanced - nrow(dat), 100*(1 - nrow(dat)/n_balanced)))

# Environments per genotype
env_per_gen <- dat %>%
  group_by(gen, gen_type) %>%
  summarise(n_envs = n_distinct(trial), .groups = "drop")

cat("\n  Environments per genotype — experimentals:\n")
exp_envs <- env_per_gen %>% filter(gen_type == "experimental")
print(table(exp_envs$n_envs))

cat("\n  Environments per genotype — checks (all should be 6):\n")
chk_envs <- env_per_gen %>% filter(gen_type == "check")
print(table(chk_envs$n_envs))

# Genotypes per environment
gen_per_env <- dat %>%
  group_by(trial, gen_type) %>%
  summarise(n_gen = n_distinct(gen), .groups = "drop")

cat("\n  Genotypes per environment:\n")
print(as.data.frame(gen_per_env), row.names = FALSE)


# =============================================================================
# SECTION 4 — YIELD DIAGNOSTICS
# =============================================================================

cat("\n========== SECTION 4: YIELD (kg/ha) ==========\n")

cat("\n  Overall:\n")
cat(sprintf("    Mean   : %.0f kg/ha\n", mean(dat$yield_kg_ha)))
cat(sprintf("    SD     : %.0f kg/ha\n", sd(dat$yield_kg_ha)))
cat(sprintf("    Min    : %.0f kg/ha\n", min(dat$yield_kg_ha)))
cat(sprintf("    Max    : %.0f kg/ha\n", max(dat$yield_kg_ha)))
cat(sprintf("    CV     : %.1f%%\n",    100 * sd(dat$yield_kg_ha) / mean(dat$yield_kg_ha)))

cat("\n  By environment (mean ± SD):\n")
dat %>%
  group_by(city, state, trial) %>%
  summarise(mean = round(mean(yield_kg_ha)), sd = round(sd(yield_kg_ha)),
            n = n(), .groups = "drop") %>%
  arrange(desc(mean)) %>%
  print(row.names = FALSE)

cat("\n  By gen_type (mean ± SD):\n")
dat %>%
  group_by(gen_type) %>%
  summarise(mean = round(mean(yield_kg_ha)), sd = round(sd(yield_kg_ha)),
            n = n(), .groups = "drop") %>%
  print(row.names = FALSE)

# Plausibility checks
n_neg    <- sum(dat$yield_kg_ha < 0)
n_low    <- sum(dat$yield_kg_ha < 2000)
n_high   <- sum(dat$yield_kg_ha > 13000)
cat(sprintf("\n  Implausible values:\n"))
cat(sprintf("    Negative      : %d\n", n_neg))
cat(sprintf("    < 2000 kg/ha  : %d\n", n_low))
cat(sprintf("    > 13000 kg/ha : %d\n", n_high))

# Environmental gradient check (should decrease CPS > RVD > JAT > SOR > SIN > DOU)
cat("\n  Environment ranking (descending mean yield):\n")
env_means <- dat %>%
  group_by(trial, city) %>%
  summarise(mean_yield = round(mean(yield_kg_ha)), .groups = "drop") %>%
  arrange(desc(mean_yield))
print(env_means, row.names = FALSE)


# =============================================================================
# SECTION 5 — DAYS TO MATURITY DIAGNOSTICS
# =============================================================================

cat("\n========== SECTION 5: DAYS TO MATURITY ==========\n")

dtm <- dat %>% select(gen, gen_type, days_to_maturity) %>% distinct()

cat(sprintf("  Range (experimentals): %d – %d days\n",
            min(dtm$days_to_maturity[dtm$gen_type == "experimental"]),
            max(dtm$days_to_maturity[dtm$gen_type == "experimental"])))
cat(sprintf("  Range (checks)       : %d – %d days\n",
            min(dtm$days_to_maturity[dtm$gen_type == "check"]),
            max(dtm$days_to_maturity[dtm$gen_type == "check"])))

cat("\n  Distribution (experimentals):\n")
cat(sprintf("    Superprecoce (100–112 d) : %d hybrids\n",
            sum(dtm$days_to_maturity[dtm$gen_type == "experimental"] <= 112)))
cat(sprintf("    Precoce      (113–125 d) : %d hybrids\n",
            sum(dtm$days_to_maturity[dtm$gen_type == "experimental"] >= 113)))

cat("\n  Checks — days to maturity:\n")
dtm %>% filter(gen_type == "check") %>%
  left_join(distinct(dat[, c("gen","cultivar_name")]), by = "gen") %>%
  select(gen, cultivar_name, days_to_maturity) %>%
  arrange(gen) %>%
  print(row.names = FALSE)


# =============================================================================
# SECTION 6 — GxE INTERACTION FINGERPRINT PLOT
# =============================================================================

cat("\n========== SECTION 6: GxE INTERACTION PLOT ==========\n")

gen_env_means <- dat %>%
  group_by(gen, gen_type, cultivar_name, trial, city) %>%
  summarise(yield = mean(yield_kg_ha), .groups = "drop") %>%
  mutate(city = factor(city, levels = dat %>% group_by(city) %>% summarise(m = mean(yield_kg_ha), .groups="drop") %>% arrange(m) %>% pull(city)),
         label = ifelse(gen_type == "check", paste0(gen, "\n", cultivar_name), gen))

check_colors <- c("#D62828","#F77F00","#FCBF49","#2A9D8F","#1D3557","#6A4C93","#118AB2","#073B4C")
check_ids    <- sort(unique(dat$gen[dat$gen_type == "check"]))
color_map    <- setNames(check_colors, check_ids)

ggplot() +
  geom_line(data = gen_env_means %>% filter(gen_type=="experimental"),
            aes(city, yield, group=gen), color="gray70", alpha=0.5, linewidth=0.35) +
  geom_line(data = gen_env_means %>% filter(gen_type=="check"),
            aes(city, yield, group=gen, color=gen), linewidth=1.1) +
  geom_point(data = gen_env_means %>% filter(gen_type=="check"),
             aes(city, yield, color=gen), size=2.5) +
  scale_color_manual(
    values = color_map,
    labels = setNames(
      paste0(check_ids, " — ", dat %>% filter(gen_type=="check") %>% distinct(gen, cultivar_name) %>% arrange(gen) %>% pull(cultivar_name)),
      check_ids
    ),
    name = "Commercial checks"
  ) +
  scale_y_continuous(labels = scales::comma, breaks = seq(2000,12000,1000)) +
  labs(title="GxE interaction fingerprint — Maize safrinha MET 2025",
       subtitle="Gray: 62 experimental hybrids  |  Colored: 8 commercial checks",
       x="Environment (stressful → favorable)", y="Mean yield (kg/ha)") +
  theme_bw(base_size=12) +
  theme(legend.position="right",
        legend.text=element_text(size=8),
        axis.text.x=element_text(angle=25,hjust=1),
        plot.title=element_text(face="bold"),
        plot.subtitle=element_text(color="gray40"),
        panel.grid.minor=element_blank())

cat("  Plot saved → gxe_fingerprint_maize.png\n\n")
cat("========== DIAGNOSTICS COMPLETE ==========\n\n")