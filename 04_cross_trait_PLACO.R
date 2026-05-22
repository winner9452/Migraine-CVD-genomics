# =============================================================================
# 04 | Cross-trait meta-analysis -- PLACO
# =============================================================================
# Detects pleiotropic loci under a composite null hypothesis using PLACO
# (Ray & Chatterjee, PLoS Genet 2020).
#
# PLACO is distributed as a single R source file. Download the official
# version and place it at  tools/PLACO/PLACO.R :
#   https://github.com/RayDebashree/PLACO
#
# IMPORTANT: follow the official README. Verify the argument names of
# var.placo() and placo(), and the names of the returned statistics, against
# the version you download -- the workflow below reflects the documented
# interface but the authors may revise it.
# =============================================================================

source("config.R")
library(data.table)

source(file.path("tools", "PLACO", "PLACO.R"))   # provides var.placo() / placo()

out_dir <- file.path(RESULTS_DIR, "placo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_placo <- function(trait1, trait2) {
  d1 <- data.table::fread(file.path(DATA_DIR, paste0(trait1, "_MTAG.txt")))
  d2 <- data.table::fread(file.path(DATA_DIR, paste0(trait2, "_MTAG.txt")))
  m  <- merge(d1, d2, by = COL$SNP, suffixes = c(".x", ".y"))

  # Z-score matrix (compute from beta/se if a Z column is absent)
  if (all(c("Z.x", "Z.y") %in% names(m))) {
    Z <- cbind(m$Z.x, m$Z.y)
  } else {
    Z <- cbind(m[[paste0(COL$BETA, ".x")]] / m[[paste0(COL$SE, ".x")]],
               m[[paste0(COL$BETA, ".y")]] / m[[paste0(COL$SE, ".y")]])
  }
  # corresponding p-value matrix (derive from Z if absent)
  P <- 2 * pnorm(-abs(Z))

  # remove variants with extreme squared Z (PLACO recommendation)
  keep <- rowSums(Z^2 > 80) == 0
  m <- m[keep, ]; Z <- Z[keep, , drop = FALSE]; P <- P[keep, , drop = FALSE]

  # estimate variance parameters, then test every variant
  VarZ <- var.placo(Z, P, p.threshold = 1e-4)
  res  <- t(vapply(seq_len(nrow(Z)),
                   function(i) unlist(placo(Z = Z[i, ], VarZ = VarZ)),
                   numeric(2)))

  data.table::data.table(
    SNP     = m[[COL$SNP]],
    CHR     = m[[paste0(COL$CHR, ".x")]],
    BP      = m[[paste0(COL$BP, ".x")]],
    T.placo = res[, "T.placo"],
    p.placo = res[, "p.placo"])
}

for (m in MIGRAINE_TRAITS) {
  for (cv in CVD_TRAITS) {
    res <- run_placo(m, cv)
    data.table::fwrite(res,
      file.path(out_dir, paste0(m, "_", cv, "_PLACO_result.txt")), sep = "\t")
  }
}

cat("PLACO cross-trait meta-analysis complete. Results in:", out_dir, "\n")
