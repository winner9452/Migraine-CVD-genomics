# =============================================================================
# 03 | Cross-trait meta-analysis -- CPASSOC (SHet)
# =============================================================================
# Detects pleiotropic loci shared by a migraine phenotype and a CVD trait
# using the SHet statistic of CPASSOC (Zhu et al., Am J Hum Genet 2015).
# SHet up-weights traits with larger effects and is robust to effect-size
# heterogeneity across phenotypes.
#
# CPASSOC is distributed by the authors as plain R source code (functions
# SHom / SHet), NOT as an installable package. Download the official source
# and place it at  tools/CPASSOC/CPASSOC.R :
#   https://github.com/Xuexia/CPASSOC
#
# IMPORTANT: the input-construction and output-filtering logic below is what
# the downstream scripts rely on. The single call to SHet() depends on the
# version of CPASSOC you download -- verify the argument names (X / SigmaO)
# and the names of the returned statistic / p-value against that file.
# =============================================================================

source("config.R")
library(data.table)

source(file.path("tools", "CPASSOC", "CPASSOC.R"))   # provides SHom() / SHet()

out_dir <- file.path(RESULTS_DIR, "cpassoc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_cpassoc <- function(trait1, trait2) {
  d1 <- data.table::fread(file.path(DATA_DIR, paste0(trait1, ".txt")))
  d2 <- data.table::fread(file.path(DATA_DIR, paste0(trait2, ".txt")))
  m  <- merge(d1, d2, by = COL$SNP, suffixes = c(".x", ".y"))

  # per-trait Z-scores
  Z <- cbind(m[[paste0(COL$BETA, ".x")]] / m[[paste0(COL$SE, ".x")]],
             m[[paste0(COL$BETA, ".y")]] / m[[paste0(COL$SE, ".y")]])

  # trait correlation matrix, estimated from variants null for both traits
  null   <- m[[paste0(COL$P, ".x")]] > 0.05 & m[[paste0(COL$P, ".y")]] > 0.05
  SigmaO <- cor(Z[null, ], use = "complete.obs")

  # SHet statistic -- see header note re: argument / output names
  sh <- SHet(X = Z, SigmaO = SigmaO)

  data.table::data.table(
    SNP    = m[[COL$SNP]],
    CHR    = m[[paste0(COL$CHR, ".x")]],
    BP     = m[[paste0(COL$BP, ".x")]],
    pval.x = m[[paste0(COL$P, ".x")]],
    pval.y = m[[paste0(COL$P, ".y")]],
    Shet   = sh$TS,
    p_Shet = sh$pval)
}

for (m in MIGRAINE_TRAITS) {
  for (cv in CVD_TRAITS) {
    res <- run_cpassoc(m, cv)
    data.table::fwrite(res,
      file.path(out_dir, paste0(m, "_", cv, "_cpassoc_result.txt")), sep = "\t")
    # genome-wide significant SHet SNPs
    data.table::fwrite(res[p_Shet < 5e-8],
      file.path(out_dir, paste0(m, "_", cv, "_cpassoc_result_sig.txt")),
      sep = "\t")
  }
}

cat("CPASSOC cross-trait meta-analysis complete. Results in:", out_dir, "\n")
