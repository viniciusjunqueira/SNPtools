#' Combina múltiplos objetos SNPDataLong (painéis diferentes)
#'
#' @param lista Lista de objetos da classe SNPDataLong
#'
#' @return Objeto SNPDataLong unificado
#' @export
combinarSNPData <- function(lista) {
  stopifnot(length(lista) > 0)

  # União de todos os SNPs
  snps_todos <- Reduce(union, lapply(lista, function(x) colnames(x@geno)))

  # Preenche cada geno com todos os SNPs (NA onde não existir)
  geno_list <- lapply(lista, function(x) {
    geno <- x@geno
    ausentes <- setdiff(snps_todos, colnames(geno))
    if (length(ausentes) > 0) {
      # Cria SnpMatrix de NAs para os SNPs ausentes
      na_block <- new("SnpMatrix", matrix(as.raw(0), nrow = nrow(geno), ncol = length(ausentes)))
      colnames(na_block) <- ausentes
      geno <- cbind(geno, na_block)
    }
    # Reordena colunas para manter a mesma ordem
    geno[, snps_todos, drop = FALSE]
  })

  # Junta os genótipos
  geno_comb <- do.call(rbind, geno_list)

  # Junta informações FAM
#   fam_comb <- do.call(rbind, lapply(lista, function(x) x@fam))

  # Junta mapas e remove SNPs duplicados mantendo primeiro
  map_all <- do.call(rbind, lapply(lista, function(x) x@map))
  map_all <- map_all[!duplicated(map_all$Name), , drop = FALSE]
  map_final <- map_all[match(snps_todos, map_all$Name), , drop = FALSE]

  # Validação final
  if (!inherits(geno_comb, "SnpMatrix")) {
    geno_comb <- methods::as(geno_comb, "SnpMatrix")
  }

  new("SNPDataLong",
      geno = geno_comb,
    #   fam  = fam_comb,
      map  = map_final,
      path = paste(sapply(lista, function(x) x@path), collapse = ";"))
}
