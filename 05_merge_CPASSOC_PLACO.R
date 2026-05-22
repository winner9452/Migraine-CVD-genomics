# =============================================================================
# 05 | Merge CPASSOC + PLACO significant loci
# =============================================================================
# A locus is retained as pleiotropic only if it is genome-wide significant in
# BOTH CPASSOC and PLACO. Significant SNPs from the two methods are
# intersected for every migraine x CVD pair, and all pairs are then stacked
# into a single table that feeds the colocalization step (script 06).
# =============================================================================

source("config.R")
library(data.table)

cpassoc_dir <- file.path(RESULTS_DIR, "cpassoc")
placo_dir   <- file.path(RESULTS_DIR, "placo")
out_dir     <- file.path(RESULTS_DIR, "merge_cpassoc_placo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

merged_all <- list()
for (m in MIGRAINE_TRAITS) {
  for (cv in CVD_TRAITS) {
    cp_file <- file.path(cpassoc_dir, paste0(m, "_", cv, "_cpassoc_result.txt"))
    pl_file <- file.path(placo_dir,   paste0(m, "_", cv, "_PLACO_result.txt"))
    if (!file.exists(cp_file) || !file.exists(pl_file)) next

    cp <- data.table::fread(cp_file)
    pl <- data.table::fread(pl_file)

    # CPASSOC: genome-wide significant SHet, both single-trait p < 1e-3
    cp <- cp[p_Shet < 5e-8 & pval.x < 1e-3 & pval.y < 1e-3]
    # PLACO: genome-wide significant composite p
    pl <- pl[p.placo < 5e-8]

    shared <- merge(cp, pl, by = "SNP", suffixes = c(".x", ".y"))
    if (nrow(shared) == 0) next

    shared[, prefix := paste0(m, "&", cv)]
    data.table::fwrite(shared,
      file.path(out_dir, paste0(m, "_", cv, "_merge_cpassoc_placo_sig.txt")),
      sep = "\t")
    merged_all[[paste0(m, "_", cv)]] <- shared
  }
}

all_loci <- data.table::rbindlist(merged_all, fill = TRUE)
data.table::fwrite(all_loci,
  file.path(out_dir, "merged_all_placo_sig.txt"), sep = "\t")

cat("Merged", nrow(all_loci), "pleiotropic loci across",
    length(merged_all), "trait pairs.\n")
cat("Combined table:",
    file.path(out_dir, "merged_all_placo_sig.txt"), "\n")
