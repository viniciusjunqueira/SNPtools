setMethod("qcSNPs", "SNPDataLong", function(x, 
                                            missing_ind = NULL,
                                            missing_snp = NULL,
                                            min_snp_cr = NULL,
                                            min_maf = NULL,
                                            hwe = NULL,
                                            snp_mono = FALSE,
                                            remove_chr = NULL,
                                            action = c("report", "filter", "both")) {
  action <- match.arg(action)

  geno <- x@geno
  map <- x@map
  keep_snps <- colnames(geno)

  # Cabeçalho informativo
  qc_header("Quality Control on SNPs")

  message("Initial number of SNPs: ", length(keep_snps))
  message("Applying quality control filters...")

  snpsum <- col.summary(geno)

  # Inicializa vetores de SNPs removidos por critério
  low_callrate_snps <- low_maf <- dev.hwe <- mono <- discard_chr <- character()

  if (!is.null(min_snp_cr)) {
    low_callrate_snps <- fQC::check.call.rate(snpsum, min.call.rate = min_snp_cr)
    keep_snps <- setdiff(keep_snps, low_callrate_snps)
    message(sprintf("  • Call rate filter: %d SNP(s) removed; %d retained.",
                    length(low_callrate_snps), length(keep_snps)))
  }

  if (!is.null(min_maf)) {
    low_maf <- fQC::check.snp.maf(snpsum, min.maf = min_maf)
    keep_snps <- setdiff(keep_snps, low_maf)
    message(sprintf("  • MAF filter: %d SNP(s) removed; %d retained.",
                    length(low_maf), length(keep_snps)))
  }

  if (!is.null(hwe)) {
    dev.hwe <- fQC::check.snp.hwe.chi2(snpsum, hwe)
    keep_snps <- setdiff(keep_snps, dev.hwe)
    message(sprintf("  • HWE filter: %d SNP(s) removed; %d retained.",
                    length(dev.hwe), length(keep_snps)))
  }

  if (snp_mono) {
    mono <- fQC::check.snp.monomorf(snpsum)
    keep_snps <- setdiff(keep_snps, mono)
    message(sprintf("  • Monomorphic filter: %d SNP(s) removed; %d retained.",
                    length(mono), length(keep_snps)))
  }

  if (!is.null(remove_chr)) {
    discard_chr <- fQC::check.snp.chromo(map, remove_chr)
    keep_snps <- setdiff(keep_snps, discard_chr)
    message(sprintf("  • Chromosome filter: %d SNP(s) removed; %d retained.",
                    length(discard_chr), length(keep_snps)))
  }

  # Se action == "report", retorna apenas as listas de removidos e mantidos
  if (action == "report") {
    return(list(
      removed_by_callrate = low_callrate_snps,
      removed_by_maf = low_maf,
      removed_by_hwe = dev.hwe,
      removed_monomorphic = mono,
      removed_by_chr = discard_chr,
      kept_snps = keep_snps
    ))
  }

  # Criação do novo objeto filtrado
  filtered_geno <- geno[, keep_snps, drop = FALSE]
  filtered_map  <- map[map$Name %in% keep_snps, , drop = FALSE]
  filtered_obj  <- new("SNPDataLong", geno = filtered_geno, map = filtered_map)

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
        removed_monomorphic = mono,
        removed_by_chr = discard_chr,
        kept_snps = keep_snps
      )
    ))
  }
})
