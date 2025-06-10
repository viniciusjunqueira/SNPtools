#' Summary para objetos da classe SNPDataLong
#'
#' @param object Objeto do tipo SNPDataLong
#' @return Imprime resumo na tela
#' @export
setMethod("summary", "SNPDataLong", function(object, ...) {
  cat("Resumo do objeto SNPDataLong\n")
  cat("----------------------------\n")

  ## Validações básicas
  if (!inherits(object@geno, "SnpMatrix")) {
    cat("Aviso: Slot 'geno' não é um SnpMatrix válido.\n")
    return(invisible(NULL))
  }

  if (!is.data.frame(object@map)) {
    cat("Aviso: Slot 'map' não é um data.frame.\n")
    return(invisible(NULL))
  }

  if (nrow(object@geno) == 0 || ncol(object@geno) == 0) {
    cat("Objeto vazio: sem indivíduos ou SNPs.\n")
    return(invisible(NULL))
  }

  ## Informações gerais
  n_ind <- nrow(object@geno)
  n_snp <- ncol(object@geno)
  cat("Indivíduos :", n_ind, "\n")
  cat("SNPs       :", n_snp, "\n\n")

  ## Dados ausentes
  n_total <- n_ind * n_snp
  n_na <- sum(is.na(object@geno))
  pct_na <- round(100 * n_na / n_total, 2)
  cat("Dados ausentes (NA):\n")
  cat(" - Total     :", n_na, "de", n_total, "\n")
  cat(" - Proporção :", pct_na, "%\n\n")

  ## Análise por cromossomo (se aplicável)
  if ("Name" %in% colnames(object@map) && "Chromosome" %in% colnames(object@map)) {
    idx <- match(colnames(object@geno), object@map$Name)
    chr_info <- object@map$Chromosome[idx]

    if (any(is.na(chr_info))) {
      cat("Aviso: Alguns SNPs não foram encontrados no mapa e serão ignorados na contagem por cromossomo.\n")
      chr_info <- chr_info[!is.na(chr_info)]
    }

    if (length(chr_info) > 0) {
      cat("Distribuição de SNPs por cromossomo:\n")
      print(table(chr_info))
      cat("\n")

      ## SNPs com NA por cromossomo
      snp_na_count <- colSums(is.na(object@geno))
      chr_na <- chr_info
      names(chr_na) <- colnames(object@geno)[idx]
      chr_na_tab <- tapply(snp_na_count > 0, chr_na, sum)
      cat("SNPs com dados ausentes por cromossomo:\n")
      print(chr_na_tab)
    } else {
      cat("Nenhum cromossomo encontrado após alinhamento com o mapa.\n")
    }
  } else {
    cat("Mapa sem colunas esperadas ('Name', 'Chromosome'); resumo por cromossomo omitido.\n")
  }

  invisible(NULL)
})
