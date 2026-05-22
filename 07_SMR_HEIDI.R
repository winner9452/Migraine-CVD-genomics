# =============================================================================
# 07 | SMR + HEIDI (multi-SNP)
# =============================================================================
# Summary-data-based Mendelian Randomization (Zhu et al., Nat Genet 2016)
# integrating GWAS summary statistics with GTEx v8 cis-eQTL data to identify
# genes whose expression is pleiotropically associated with each trait. The
# HEIDI test distinguishes genuine pleiotropy from linkage.
#
# SMR is an open-source command-line binary -- no in-house wrapper required:
#   https://yanglab.westlake.edu.cn/software/smr/
#
# Inputs (place under reference/ and data/):
#   - data/gwas/ma/<TRAIT>.ma          GWAS in SMR .ma format
#                                      (cols: SNP A1 A2 freq b se p n)
#   - reference/eQTL_besd/             GTEx v8 cis-eQTL BESD files, one set
#                                      (.besd/.esi/.epi) per tissue
#   - reference/1000G_EUR/EUR          1000 Genomes EUR PLINK bfile (LD ref)
# =============================================================================

source("config.R")
library(data.table)

GWAS_MA_DIR <- file.path(DATA_DIR, "ma")            # <TRAIT>.ma files
EQTL_DIR    <- file.path(REF_DIR, "eQTL_besd")      # GTEx v8 BESD files
LD_BFILE    <- file.path(REF_DIR, "1000G_EUR", "EUR")

all_traits <- c(MIGRAINE_TRAITS, CVD_TRAITS)

# each BESD set is a prefix shared by .besd / .esi / .epi files
eqtl_prefixes <- unique(tools::file_path_sans_ext(
  list.files(EQTL_DIR, pattern = "\\.(besd|esi|epi)$", full.names = TRUE)))

out_root <- file.path(RESULTS_DIR, "smr")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

## ---- Step 1: run SMR-multi + HEIDI for every trait x tissue -----------------
for (trait in all_traits) {
  dir.create(file.path(out_root, trait), showWarnings = FALSE)
  for (eqtl in eqtl_prefixes) {
    tissue <- basename(eqtl)
    system2(SMR_BIN, c(
      "--bfile",          LD_BFILE,
      "--gwas-summary",   file.path(GWAS_MA_DIR, paste0(trait, ".ma")),
      "--beqtl-summary",  eqtl,
      "--smr-multi",                     # multi-SNP-based SMR test
      "--cis-wnd",        2000,          # cis window, kb
      "--peqtl-smr",      "5e-8",
      "--ld-multi-snp",   0.1,
      "--heidi-mtd",      1,
      "--diff-freq",      0.2,
      "--diff-freq-prop", 0.05,
      "--maf",            0.01,
      "--thread-num",     8,
      "--out", file.path(out_root, trait, tissue)
    ))
  }
}

## ---- Step 2: merge .msmr output across tissues, per trait -------------------
for (trait in all_traits) {
  files <- list.files(file.path(out_root, trait),
                       pattern = "\\.msmr$", full.names = TRUE)
  if (length(files) == 0) next
  merged <- data.table::rbindlist(lapply(files, function(f) {
    d <- data.table::fread(f)
    d[, tissue := tools::file_path_sans_ext(basename(f))]
    d
  }), fill = TRUE)
  data.table::fwrite(merged,
    file.path(out_root, paste0(trait, "_merged.msmr")), sep = "\t")
}

## ---- Step 3: significant pleiotropic genes ----------------------------------
# Retain: significant multi-SNP SMR (p_SMR_multi < 0.05) AND
#         non-significant HEIDI (p_HEIDI > 0.05, i.e. no evidence of linkage).
msmr_files <- list.files(out_root, pattern = "_merged\\.msmr$",
                         full.names = TRUE)
all_smr <- data.table::rbindlist(lapply(msmr_files, function(f) {
  d <- data.table::fread(f)
  d[, trait := sub("_merged\\.msmr$", "", basename(f))]
  d
}), fill = TRUE)

sig <- all_smr[p_SMR_multi < 0.05 & p_HEIDI > 0.05]
data.table::fwrite(sig, file.path(out_root, "SMR_significant_genes.csv"))

cat("SMR + HEIDI complete.",
    nrow(sig), "significant gene-trait associations.\n")
cat("Result:", file.path(out_root, "SMR_significant_genes.csv"), "\n")
