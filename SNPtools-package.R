#' SNPtools
#'
#' S4 tools for reading and organizing genetic data.
#'
#' @docType package
#' @name SNPtools
#' @author
#' Vinícius Junqueira \email{junqueiravinicius@hotmail.com}
#'
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics hist par text
#' @importFrom methods new as
#' @importFrom stats dist hclust pchisq prcomp sd
#' @importFrom utils read.table write.table
#' @importFrom magrittr %>%
#' @importFrom snpStats col.summary row.summary snp.pre.multiply snp.post.multiply
#' @importFrom MASS isoMDS
#' @importFrom reshape2 acast
#' @importFrom anticlust fast_anticlustering
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal theme element_rect ggsave
#' @importClassesFrom snpStats SnpMatrix
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("Name", "PC1", "PC2", "Group"))
}

dummy_imports <- function() {
  MASS::isoMDS
  anticlust::fast_anticlustering
  dplyr::select
  ggplot2::ggplot
  ggplot2::aes
  ggplot2::geom_point
  ggplot2::labs
  ggplot2::theme_minimal
  ggplot2::theme
  ggplot2::element_rect
  ggplot2::ggsave
  grDevices::dev.off
  graphics::hist
  graphics::par
  graphics::text
  magrittr::`%>%`
  reshape2::acast
  snpStats::col.summary
  snpStats::row.summary
  snpStats::snp.pre.multiply
  snpStats::snp.post.multiply
  stats::dist
  stats::hclust
  stats::pchisq
  stats::prcomp
  stats::sd
  utils::read.table
  utils::write.table
  Rcpp::evalCpp
  invisible(NULL)
}

if (FALSE) {
  dummy_imports()
}
