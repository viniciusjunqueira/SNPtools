#' Controle de Qualidade para SNPDataLong com critérios opcionais
#'
#' Permite análise de qualidade genotípica com critérios definidos pelo usuário.
#'
#' @param x Objeto da classe SNPDataLong
#' @param missing_ind Proporção máxima de dados faltantes permitida por indivíduo (opcional)
#' @param missing_snp Proporção máxima de dados faltantes permitida por SNP (opcional)
#' @param min_snp_cr Mínimo call rate aceitável para SNPs 
#' @param min_maf Frequência alélica mínima permitida para SNPs (opcional)
#' @param hwe Identifica SNPs com desvio de HWE
#' @param snp_mono Identifica SNPs monomórficos
#' @param remove_chr Descarta SNPs localizados nos cromossomos listados
#' @param action "report", "filter" ou "both"
#'
#' @return Dependendo do argumento `action`, retorna:
#' - "report": lista com SNPs removidos por critério;
#' - "filter": objeto `SNPDataLong` com SNPs filtrados;
#' - "both": lista contendo o objeto filtrado e o relatório.
#'
#' @examples
#' ## Not run:
#' set.seed(123)
#' mat <- matrix(sample(c(0, 1, 2, NA), 100, replace = TRUE, prob = c(0.4, 0.4, 0.15, 0.05)),
#'               nrow = 10, ncol = 10)
#' colnames(mat) <- paste0("snp", 1:10)
#' rownames(mat) <- paste0("ind", 1:10)
#' map <- data.frame(Name = colnames(mat), Chrom = 1, Position = 1:10)
#' x <- new("SNPDataLong", geno = mat, map = map)
#'
#' qcSNPs(x, min_snp_cr = 0.8, min_maf = 0.05, snp_mono = TRUE, action = "report")
#' ## End(Not run)
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
                                            snp_mono = FALSE,
                                            remove_chr = NULL,
                                            action = c("report", "filter", "both")) {
  action <- match.arg(action)

  geno <- x@geno
  map <- x@map
  keep_snps <- colnames(geno)
  message("Number of SNPs before QC: ", length(keep_snps))

  snpsum <- col.summary(geno)

  # Armazena os removidos por critério
  low_callrate_snps <- low_maf <- dev.hwe <- mono <- discard_chr <- character()

  if (!is.null(min_snp_cr)) {
    low_callrate_snps <- fQC::check.call.rate(snpsum, min.call.rate = min_snp_cr)
    keep_snps <- setdiff(keep_snps, low_callrate_snps)
    message("After call rate filter: ", length(keep_snps))
  }

  if (!is.null(min_maf)) {
    low_maf <- fQC::check.snp.maf(snpsum, min.maf = min_maf)
    keep_snps <- setdiff(keep_snps, low_maf)
    message("After MAF filter: ", length(keep_snps))
  }

  if (!is.null(hwe)) {
    dev.hwe <- fQC::check.snp.hwe.chi2(snpsum, hwe)
    keep_snps <- setdiff(keep_snps, dev.hwe)
    message("After HWE filter: ", length(keep_snps))
  }

  if (snp_mono) {
    mono <- fQC::check.snp.monomorf(snpsum)
    keep_snps <- setdiff(keep_snps, mono)
    message("After monomorphic SNP filter: ", length(keep_snps))
  }

  if (!is.null(remove_chr)) {
    discard_chr <- fQC::check.snp.chromo(map, remove_chr)
    keep_snps <- setdiff(keep_snps, discard_chr)
    message("After chromosome filter: ", length(keep_snps))
  }

  # Report apenas
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

  # Filtro
  filtered_geno <- geno[, keep_snps, drop = FALSE]
  filtered_map <- map[map$Name %in% keep_snps, , drop = FALSE]
  filtered_obj <- new("SNPDataLong", geno = filtered_geno, map = filtered_map)

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
