# =============================================================================
#  SOYBEAN MULTI-ENVIRONMENT TRIAL SIMULATOR
#  Discipline: Biometria da Interação Genótipos × Ambientes — GMP0164/UFG 2026
#  Prof. Rafael Tassinari Resende
# =============================================================================
#  Design:   Balanced MET — all genotypes × all environments × all blocks
#  Crop:     Soybean (Glycine max) — advanced inbred lines F7
#  Response: Grain yield (kg/ha)
#  Lessons:  2 (ANOVA), 3 (Finlay-Wilkinson), 4 (stability indices),
#            5 (AMMI / GGE biplot)
# =============================================================================
#  Dependencies:
#    install.packages(c("AlphaSimR", "MASS", "dplyr"))
# =============================================================================


library(AlphaSimR)
library(MASS)
library(dplyr)


# =============================================================================
# SECTION 1 — CONFIGURATION
# Edit parameters here. Do not change code below unless necessary.
# =============================================================================

# --- Seed ---
SEED <- 2026

# --- Experimental structure ---
N_LINES   <- 36          # F7 inbred lines under evaluation
N_CHECKS  <- 4           # commercial check cultivars (present in all envs)
N_ENVS    <- 6           # number of environments (locations)
N_BLOCKS  <- 3           # complete blocks per environment (RCBD)
# Total observations: (N_LINES + N_CHECKS) x N_ENVS x N_BLOCKS = 720

# --- Phenotypic parameters (kg/ha) ---
MU        <- 3600        # grand mean
SD_ENV    <- 750         # environment main effect (largest — typical in MET)
SD_GEN    <- 380         # genotype main effect
SD_GXE    <- 520         # GxE interaction
SD_BLOCK  <- 160         # block within environment
SD_ERROR  <- 420         # plot-level residual error

# Resulting heritabilities (computed in Section 5):
#   h2 at plot level  ≈ 0.23   [V_G / (V_G + V_GxE + V_Block + V_Error)]
#   H2 at means level ≈ 0.72   [V_G / (V_G + V_GxE/e + V_Error/(e*r))]

# --- Genotype groups (must sum to N_LINES) ---
# Groups define the GxE pattern each genotype will express:
#   stable     : high mean, beta ≈ 1.0, low deviation  → broad adaptation
#   responsive : high mean in favorable envs, beta > 1  → specific adaptation
#   rustic     : better in unfavorable envs, beta < 1   → specific adaptation
#   unstable   : any mean, high GxE deviation           → unpredictable
#   inferior   : low mean, stable by low performance    → discard candidates
N_STABLE     <- 8
N_RESPONSIVE <- 8
N_RUSTIC     <- 6
N_UNSTABLE   <- 8
N_INFERIOR   <- 6
# Check: N_STABLE + N_RESPONSIVE + N_RUSTIC + N_UNSTABLE + N_INFERIOR == N_LINES

# --- GxE magnitude multiplier per group (relative to SD_GXE) ---
# Controls how strongly each group deviates across environments
GXE_MULT <- list(
  check      = 0.40,   # checks: well-known, predictable
  stable     = 0.50,   # small GxE — the ideal genotype profile
  responsive = 1.20,   # large GxE with positive env gradient covariance
  rustic     = 1.00,   # moderate GxE with negative env gradient covariance
  unstable   = 2.20,   # large random GxE — pedagogically important contrast
  inferior   = 0.60    # small GxE but low mean
)

# --- Environments: Cerrado biome, Brazil ---
# env_index: deviation from MU (kg/ha) — used as Finlay-Wilkinson index
# Gradient spans ~2000 kg/ha from most stressful to most productive
ENV <- data.frame(
 	id    = 	c("E1_Anapolis_GO", "E2_Luziania_GO", "E3_Jatai_GO",
            		     "E4_RioVerde_GO", "E5_Uberlandia_MG", "E6_Sorriso_MT"),
	index = c(-920, -480, -60, 390, 680, 1150),
	lat   = c(-16.33, -16.25, -17.88, -17.80, -18.91, -12.54),  # Jatai e Uberlandia
  	lon   = c(-48.95, -47.95, -51.72, -51.93, -48.27, -55.71)   # Rio Verde
)

# --- AlphaSimR genomic architecture ---
N_FOUNDERS    <- 20     # founder individuals (genetic diversity base)
N_CHR         <- 20     # soybean chromosomes
SEG_SITES_CHR <- 300    # segregating sites per chromosome
N_QTL_CHR     <- 10     # QTLs per chromosome → 200 total (polygenic)
N_SNP_CHR     <- 250    # SNPs per chromosome → 5000 total
H2_PLOT       <- 0.20   # heritability at plot level (AlphaSimR parameter)
N_SELFING     <- 6      # selfing generations F2 → F7
SELECT_PROP   <- 0.40   # proportion selected per generation (truncation)

# --- Output files ---
OUT_PHENO  <- "soy_MET_phenotypic.txt"
OUT_GRM    <- "soy_MET_GRM.csv"
OUT_SNP    <- "soy_MET_SNP.csv"
OUT_META   <- "soy_MET_genotype_metadata.csv"


# =============================================================================
# SECTION 2 — GENOMIC SIMULATION (AlphaSimR)
# =============================================================================
# Delivers:
#   bv_scaled : breeding values (kg/ha) used as genotype effects in Section 4
#   GRM       : genomic relationship matrix — 40 x 40  (for Lessons 6–8)
#   SNP_mat   : SNP marker matrix — 40 x 5000          (for Lessons 13–15)
# =============================================================================

set.seed(SEED)
cat("\n--- Section 2: Genomic simulation (AlphaSimR) ---\n")

n_total <- N_LINES + N_CHECKS

# Founder haplotypes
founder_hap <- quickHaplo(
  nInd     = N_FOUNDERS,
  nChr     = N_CHR,
  segSites = SEG_SITES_CHR
)

# Simulation parameters
SP <- SimParam$new(founder_hap)
SP$addTraitA(nQtlPerChr = N_QTL_CHR, mean = MU, var = SD_GEN^2)
SP$addSnpChip(nSnpPerChr = N_SNP_CHR)
SP$setVarE(h2 = H2_PLOT)

# Founder population
pop <- newPop(founder_hap, simParam = SP)

# Selfing cycles with truncation selection
for (cycle in seq_len(N_SELFING)) {
  pop <- self(pop, nProgeny = 10, simParam = SP)
  n_keep <- min(n_total * 4, nInd(pop))
  pop <- selectInd(pop, nInd = n_keep, trait = 1, simParam = SP)
}

# Final selection of n_total genotypes
pop_final <- selectInd(pop, nInd = n_total, trait = 1, simParam = SP)

# Extract outputs
bv_raw  <- gv(pop_final, simParam = SP)[, 1]
SNP_mat <- pullSnpGeno(pop_final, simParam = SP)       # n_total x 5000
GRM     <- A.mat(SNP_mat - 1)                          # centered GRM

# Assign IDs
gen_ids <- c(sprintf("CHK%02d", seq_len(N_CHECKS)),
             sprintf("L%03d",   seq_len(N_LINES)))
names(bv_raw)      <- gen_ids
rownames(GRM)      <- colnames(GRM) <- gen_ids
rownames(SNP_mat)  <- gen_ids

# Rescale BVs: center = 0, sd = SD_GEN (kg/ha units)
bv_scaled <- (bv_raw - mean(bv_raw)) / sd(bv_raw) * SD_GEN
names(bv_scaled) <- gen_ids

cat(sprintf("  Genotypes  : %d (%d checks + %d lines)\n",
            n_total, N_CHECKS, N_LINES))
cat(sprintf("  SNP markers: %d\n", ncol(SNP_mat)))
cat(sprintf("  BV range   : %.0f to %.0f kg/ha\n",
            min(bv_scaled), max(bv_scaled)))


# =============================================================================
# SECTION 3 — GENOTYPE GROUP ASSIGNMENT
# =============================================================================

checks_id     <- sprintf("CHK%02d", seq_len(N_CHECKS))
stable_id     <- sprintf("L%03d", 1:N_STABLE)
responsive_id <- sprintf("L%03d", (N_STABLE + 1):(N_STABLE + N_RESPONSIVE))
rustic_id     <- sprintf("L%03d", (N_STABLE + N_RESPONSIVE + 1):
                            (N_STABLE + N_RESPONSIVE + N_RUSTIC))
unstable_id   <- sprintf("L%03d", (N_STABLE + N_RESPONSIVE + N_RUSTIC + 1):
                            (N_STABLE + N_RESPONSIVE + N_RUSTIC + N_UNSTABLE))
inferior_id   <- sprintf("L%03d", (N_STABLE + N_RESPONSIVE + N_RUSTIC +
                                     N_UNSTABLE + 1):N_LINES)

group_of <- c(
  setNames(rep("check",      N_CHECKS),      checks_id),
  setNames(rep("stable",     N_STABLE),      stable_id),
  setNames(rep("responsive", N_RESPONSIVE),  responsive_id),
  setNames(rep("rustic",     N_RUSTIC),      rustic_id),
  setNames(rep("unstable",   N_UNSTABLE),    unstable_id),
  setNames(rep("inferior",   N_INFERIOR),    inferior_id)
)


# =============================================================================
# SECTION 4 — PHENOTYPIC SIMULATION
# =============================================================================
# Model:
#   y_ijk = mu + env_j + gen_i + gxe_ij + block_k(j) + epsilon_ijk
#
# GxE profiles are structured by group:
#   responsive → positive covariance with env_index  (beta > 1)
#   rustic     → negative covariance with env_index  (beta < 1)
#   unstable   → large random variance               (high s2d)
#   stable     → small random variance               (low s2d, beta ≈ 1)
#   check      → smallest variance                   (predictable)
#   inferior   → small variance, pulled down by low BV
# =============================================================================

cat("\n--- Section 4: Phenotypic simulation ---\n")

n_envs  <- nrow(ENV)
env_ids <- ENV$id
env_idx <- setNames(ENV$index, env_ids)   # named Finlay-Wilkinson index

# ---- GxE effect matrix (genotypes x environments) --------------------------

gxe_mat <- matrix(0, nrow = n_total, ncol = n_envs,
                  dimnames = list(gen_ids, env_ids))

env_scaled <- as.numeric(scale(ENV$index))   # standardized gradient

for (gid in gen_ids) {

  grp   <- group_of[gid]
  sigma <- SD_GXE * GXE_MULT[[grp]]

  if (grp == "responsive") {
    slope          <- runif(1, 0.40, 0.80)
    noise          <- rnorm(n_envs, 0, sigma * 0.35)
    gxe_mat[gid, ] <- env_scaled * slope * sigma + noise

  } else if (grp == "rustic") {
    slope          <- runif(1, 0.30, 0.65)
    noise          <- rnorm(n_envs, 0, sigma * 0.35)
    gxe_mat[gid, ] <- -env_scaled * slope * sigma + noise

  } else {
    # stable, unstable, inferior, check: random across environments
    gxe_mat[gid, ] <- rnorm(n_envs, 0, sigma)
  }
}

# Center GxE matrix (rows and columns → AMMI convention)
gxe_mat <- sweep(gxe_mat, 1, rowMeans(gxe_mat), "-")
gxe_mat <- sweep(gxe_mat, 2, colMeans(gxe_mat), "-")

# ---- Build fully balanced design -------------------------------------------

design <- expand.grid(
  gen   = gen_ids,
  env   = env_ids,
  block = paste0("B", seq_len(N_BLOCKS)),
  stringsAsFactors = FALSE
)

# ---- Block effects (nested within environment) -----------------------------

blk_key <- paste(design$env, design$block, sep = ":")
blk_ids <- unique(blk_key)
blk_eff <- setNames(rnorm(length(blk_ids), 0, SD_BLOCK), blk_ids)

# ---- Simulate response variable --------------------------------------------

design$y <- with(design,
  MU +
  env_idx[env] +
  bv_scaled[gen] +
  gxe_mat[cbind(gen, env)] +
  blk_eff[paste(env, block, sep = ":")] +
  rnorm(nrow(design), 0, SD_ERROR)
)

design$y <- round(design$y, 1)

# ---- Attach metadata -------------------------------------------------------

design <- design %>%
  left_join(
    data.frame(env       = ENV$id,
               env_index = ENV$index,
               lat       = ENV$lat,
               lon       = ENV$lon,
               stringsAsFactors = FALSE),
    by = "env"
  ) %>%
  mutate(group = unname(group_of[gen])) %>%
  select(gen, group, env, env_index, lat, lon, block, y)

cat(sprintf("  Total observations : %d\n", nrow(design)))
cat(sprintf("  y range (kg/ha)    : %.0f – %.0f\n",
            min(design$y), max(design$y)))
cat(sprintf("  y mean             : %.0f kg/ha\n", mean(design$y)))


# =============================================================================
# SECTION 5 — VARIANCE COMPONENTS & HERITABILITY (analytical)
# =============================================================================

vG   <- SD_GEN^2
vGxE <- SD_GXE^2
vB   <- SD_BLOCK^2
vE   <- SD_ERROR^2
vP   <- vG + vGxE + vB + vE

h2_plot  <- vG / vP
H2_means <- vG / (vG + vGxE / N_ENVS + vE / (N_ENVS * N_BLOCKS))

cat("\n--- Section 5: Variance components ---\n")
cat(sprintf("  V_G   = %7.0f  (%.1f%%)\n", vG,   100 * vG   / vP))
cat(sprintf("  V_GxE = %7.0f  (%.1f%%)\n", vGxE, 100 * vGxE / vP))
cat(sprintf("  V_Blk = %7.0f  (%.1f%%)\n", vB,   100 * vB   / vP))
cat(sprintf("  V_Err = %7.0f  (%.1f%%)\n", vE,   100 * vE   / vP))
cat(sprintf("  h2 (plot level)         : %.2f\n", h2_plot))
cat(sprintf("  H2 (means, %de x %dr)  : %.2f\n", N_ENVS, N_BLOCKS, H2_means))


# =============================================================================
# SECTION 6 — EXPORT
# =============================================================================

cat("\n--- Section 6: Exporting files ---\n")

# Phenotypic data (main dataset for lessons 2–5)
write.table(design, file = OUT_PHENO,
            row.names = FALSE, quote = FALSE, sep = "\t")
cat("  Phenotypic data   →", OUT_PHENO, "\n")

# GRM — genomic relationship matrix (lessons 6–8)
write.csv(as.data.frame(GRM), file = OUT_GRM)
cat("  GRM matrix        →", OUT_GRM, "\n")

# SNP matrix (lessons 13–15)
snp_df <- cbind(data.frame(gen = rownames(SNP_mat)), as.data.frame(SNP_mat))
write.csv(snp_df, file = OUT_SNP, row.names = FALSE)
cat(sprintf("  SNP matrix        → %s  (%d genotypes x %d markers)\n",
            OUT_SNP, nrow(SNP_mat), ncol(SNP_mat)))

# Genotype metadata
meta <- data.frame(
  gen      = names(group_of),
  group    = unname(group_of),
  bv_kg_ha = round(bv_scaled[names(group_of)], 1)
)
write.csv(meta, file = OUT_META, row.names = FALSE)
cat("  Genotype metadata →", OUT_META, "\n")

cat("\nDone.\n\n")


# =============================================================================
# SECTION 7 — QUICK SANITY CHECK  (optional)
# =============================================================================
# Uncomment to verify variance components via lme4.
# Expected: VarCorr should recover approximate SD values from CONFIG.
#
# library(lme4)
# mod <- lmer(y ~ -1 + env + (1|gen) + (1|gen:env) + (1|env:block),
#             data = design, REML = TRUE)
# print(VarCorr(mod))
# =============================================================================
