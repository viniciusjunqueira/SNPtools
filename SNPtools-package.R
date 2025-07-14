#' SNPtools
#'
#' S4 tools for reading and organizing genetic data.
#'
#' @docType package
#' @name SNPtools
#' @author
#' Vinícius Junqueira \email{junqueiravinicius@hotmail.com}
"_PACKAGE"

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
