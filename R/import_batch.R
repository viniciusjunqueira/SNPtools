#' Importa múltiplos conjuntos genotípicos a partir de uma lista de configurações
#'
#' @param config_list Lista com configurações (cada elemento é uma lista com path, fields, sep, etc.)
#'
#' @return Objeto SNPDataLong unificado
#' @export
import_geno_list <- function(config_list) {
  stopifnot(is.list(config_list))

  resultados <- lapply(config_list, function(cfg) {
    getGeno(
      path = cfg$path,
      fields = cfg$fields,
      codes = if (!is.null(cfg$codes)) cfg$codes else c("A", "B"),
      threshold = if (!is.null(cfg$threshold)) cfg$threshold else 0.15,
      sep = if (!is.null(cfg$sep)) cfg$sep else "\t",
      skip = if (!is.null(cfg$skip)) cfg$skip else 0,
      verbose = if (!is.null(cfg$verbose)) cfg$verbose else TRUE
    )
  })

  combinarSNPData(resultados)
}
