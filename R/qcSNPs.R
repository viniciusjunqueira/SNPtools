#' Quality Control for SNPDataLong with optional criteria
#'
#' Applies flexible quality control filters on an object of class \code{SNPDataLong}. 
#' Supports call rate filtering, minor allele frequency (MAF), Hardy-Weinberg equilibrium (HWE),
#' removal of monomorphic SNPs, exclusion of specific chromosomes, optionally removing SNPs without positions,
#' and optionally removing SNPs at the same genomic position.
#'
#' @param x An object of class SNPDataLong.
#' @param missing_ind Maximum allowed proportion of missing data per individual (optional). *[Currently not implemented]*
#' @param missing_snp Maximum allowed proportion of missing data per SNP (optional). *[Currently not implemented]*
#' @param min_snp_cr Minimum acceptable call rate for SNPs (e.g., 0.95).
#' @param min_maf Minimum minor allele frequency allowed for SNPs (e.g., 0.05).
#' @param hwe p-value threshold for Hardy-Weinberg equilibrium test (e.g., 1e-6).
#' @param snp_position Logical. If TRUE, removes SNPs mapped to the same position, keeping the one with highest MAF.
#' @param no_position Logical. If TRUE, removes SNPs without defined map positions.
#' @param snp_mono Logical. If TRUE, removes monomorphic SNPs.
#' @param remove_chr Character vector of chromosomes to exclude (e.g., c("X", "Y")).
#' @param action One of "report" (list of removed SNPs), "filter" (returns filtered SNPDataLong), or "both" (both results).
#'
#' @return Depending on the action argument:
#' - "report": list of SNPs removed by each filter and SNPs retained.
#' - "filter": filtered SNPDataLong object.
#' - "both": list containing the filtered object and detailed report.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' mat <- matrix(sample(c(0, 1, 2, NA), 100, replace = TRUE, prob = c(0.4, 0.4, 0.15, 0.05)),
#'               nrow = 10, ncol = 10)
#' colnames(mat) <- paste0("snp", 1:10)
#' rownames(mat) <- paste0("ind", 1:10)
#' map <- data.frame(Name = colnames(mat), Chromosome = 1, Position = 1:10)
#' x <- new("SNPDataLong", geno = mat, map = map, path = "dummy_path", xref_path = rep("chip1", 10))
#'
#' # Example with multiple filters
#' qcSNPs(x, min_snp_cr = 0.8, min_maf = 0.05, snp_mono = TRUE, no_position = TRUE, snp_position = TRUE, action = "filter")
#' }
#'
#' @importFrom reshape2 acast
#' @import data.table
#' @export
setGeneric("qcSNPs", function(x, ...) standardGeneric("qcSNPs"))

setMethod("qcSNPs", "SNPDataLong", function(x, 
                                            missing_ind = NULL,
                                            missing_snp = NULL,
                                            min_snp_cr = NULL,
                                            min_maf = NULL,
                                            hwe = NULL,
                                            snp_position = NULL,
                                            no_position = NULL,
                                            snp_mono = FALSE,
                                            remove_chr = NULL,
                                            action = c("report", "filter", "both")) {
  action <- match.arg(action)

  geno <- x@geno
  map <- x@map
  keep_snps <- colnames(geno)

  qc_header("🚦 Quality Control on SNPs")

  message("ℹ️ Initial number of SNPs: ", length(keep_snps))
  message("🔬 Applying quality control filters...")

  snpsum <- col.summary(geno)

  # Initialize vectors to store SNPs to be removed
  low_callrate_snps <- low_maf <- dev.hwe <- mono <- discard_chr <- snpstoremove <- no_pos <- character()

  # Call rate filter
  if (!is.null(min_snp_cr)) {
    low_callrate_snps <- fQC::check.call.rate(snpsum, min.call.rate = min_snp_cr)
    keep_snps <- setdiff(keep_snps, low_callrate_snps)
    message(sprintf("  • Call rate filter: %d SNP(s) removed; %d retained.",
                    length(low_callrate_snps), length(keep_snps)))
  }

  # MAF filter
  if (!is.null(min_maf)) {
    low_maf <- fQC::check.snp.maf(snpsum, min.maf = min_maf)
    keep_snps <- setdiff(keep_snps, low_maf)
    message(sprintf("  • MAF filter: %d SNP(s) removed; %d retained.",
                    length(low_maf), length(keep_snps)))
  }

  # Hardy-Weinberg filter
  if (!is.null(hwe)) {
    dev.hwe <- fQC::check.snp.hwe.chi2(snpsum, hwe)
    keep_snps <- setdiff(keep_snps, dev.hwe)
    message(sprintf("  • HWE filter: %d SNP(s) removed; %d retained.",
                    length(dev.hwe), length(keep_snps)))
  }

  # No position filter
  if (!is.null(no_position) && no_position) {
    no_pos <- fQC::check.snp.no.position(map)
    keep_snps <- setdiff(keep_snps, no_pos)
    message(sprintf("  • No-position filter: %d SNP(s) removed; %d retained.",
                    length(no_pos), length(keep_snps)))
  }

  # Same genomic position filter
  if (!is.null(snp_position) && snp_position) {
    snp_same <- fQC::check.snp.same.position(map)
    message("  • Positions with overlapping SNPs: ", length(snp_same))
    n <- length(snp_same)
    if (n > 0) {
      for (i in seq_len(n)) {
        snpsum1 <- snpsum[snp_same[[i]], , drop = FALSE]
        snp.high.maf <- rownames(snpsum1[snpsum1[, "MAF"] == max(snpsum1[, "MAF"]), , drop = FALSE])[1]
        snpstoremove <- union(snpstoremove, setdiff(rownames(snpsum1), snp.high.maf))
      }
      keep_snps <- setdiff(keep_snps, snpstoremove)
      message(sprintf("  • Same-position filter: %d SNP(s) removed; %d retained.",
                      length(snpstoremove), length(keep_snps)))
    }
  }

  # Monomorphic SNP filter
  if (snp_mono) {
    mono <- fQC::check.snp.monomorf(snpsum)
    keep_snps <- setdiff(keep_snps, mono)
    message(sprintf("  • Monomorphic filter: %d SNP(s) removed; %d retained.",
                    length(mono), length(keep_snps)))
  }

  # Chromosome filter
  if (!is.null(remove_chr)) {
    discard_chr <- fQC::check.snp.chromo(map, remove_chr)
    keep_snps <- setdiff(keep_snps, discard_chr)
    message(sprintf("  • Chromosome filter: %d SNP(s) removed; %d retained.",
                    length(discard_chr), length(keep_snps)))
  }

  # Report only
  if (action == "report") {
    return(list(
      removed_by_callrate = low_callrate_snps,
      removed_by_maf = low_maf,
      removed_by_hwe = dev.hwe,
      removed_no_position = no_pos,
      removed_same_position = snpstoremove,
      removed_monomorphic = mono,
      removed_by_chr = discard_chr,
      kept_snps = keep_snps
    ))
  }

  # Apply filter
  filtered_geno <- geno[, keep_snps, drop = FALSE]
  filtered_map  <- map[map$Name %in% keep_snps, , drop = FALSE]
  filtered_obj  <- new("SNPDataLong", 
                       geno = filtered_geno, 
                       map = filtered_map,
                       path = x@path,
                       xref_path = x@xref_path)

  if (action == "filter") {
    return(filtered_obj)
  }

  if (action == "both") {
    return(list(
      filtered = filtered_obj,
      report = list(
        removed_by_callrate = low_callrate_snps,
        removed_by_maf = low_maf,
        removed_by_hwe = dev.hwe,
        removed_no_position = no_pos,
        removed_same_position = snpstoremove,
        removed_monomorphic = mono,
        removed_by_chr = discard_chr,
        kept_snps = keep_snps
      )
    ))
  }
})
