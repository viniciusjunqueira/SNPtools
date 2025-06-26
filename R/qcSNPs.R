#' Controle de Qualidade para SNPDataLong com critérios opcionais
#'
#' Permite análise de qualidade genotípica com critérios definidos pelo usuário.
#'
#' @param x Objeto da classe SNPDataLong
#' @param missing_ind Proporção máxima de dados faltantes permitida por indivíduo (opcional)
#' @param missing_snp Proporção máxima de dados faltantes permitida por SNP (opcional)
#' @param min_snp_cr
#' @param min_maf Frequência alélica mínima permitida para SNPs (opcional)
#' @param hwe
#' @param snp_mono
#' @param action "report", "filter" ou "both"
# ' @param plot Se TRUE, exibe gráficos de qualidade via fQC
#' @return Lista com sumário dos filtros e (opcionalmente) os dados filtrados
# ' @importFrom fQC missing_summary maf_summary plotMissing
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
                                            action = c("report", "filter", "both")
                                            ) {
  action <- match.arg(action)

  # Acessa a matriz de genótipos no formato longo
  geno <- x@geno
#   mat <- reshape2::acast(geno, id ~ snp, value.var = "genotype")

  # Inicializa vetores de filtro como TRUE para todos
#   ids_ok <- rep(TRUE, nrow(geno))
  keep_snps <- colnames(geno) ; print(length(keep_snps))

  # Estatísticas de qualidade
  snpsum <- col.summary(geno)
  sample.qc <- row.summary(geno)

  snp_cr <- miss_ind <- miss_snp <- maf_stat <- NULL

  if(!is.null(min_snp_cr)){
    low_callrate_snps <- fQC::check.call.rate(snpsum, min.call.rate = min_snp_cr)
    keep_snps <- setdiff(keep_snps, low_callrate_snps)
  }
  print(length(keep_snps))

  if (!is.null(min_maf)) {
    low_maf <- fQC::check.snp.maf(snp.summary = snpsum, min.maf = min_maf)
    keep_snps <- setdiff(keep_snps, low_maf)
  }
  print(length(keep_snps))

  if (!is.null(hwe)) {
    dev.hwe <- fQC::check.snp.hwe.chi2(snpsum, hwe)
    keep_snps <- setdiff(keep_snps, dev.hwe)
  }

  if(snp_mono){
    mono <- fQC::check.snp.monomorf(snpsum)
    keep_snps <- setdiff(keep_snps, mono)
  }

#   # Aplicação dos filtros, se for o caso
  filtered_object <- NULL
  if (action %in% c("filter", "both")) {
    x <- geno[,keep_snps]
    filtered_object <- x

    return(list(
      # summary = list(individual = miss_ind, snp = miss_snp, maf = maf_stat),
      geno = filtered_object
    ))
  }  
})
