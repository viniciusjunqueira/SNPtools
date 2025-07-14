if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("PC1", "PC2", "Group"))
}


#' Convert geno slot from SNPDataLong to a data.frame
#'
#' Converts the genotype matrix (geno slot) of a SNPDataLong object to a data.frame,
#' with optional centering and scaling per SNP (column).
#'
#' @param object An object of class SNPDataLong.
#' @param center Logical or numeric. If TRUE (default FALSE), center columns to mean zero.
#' @param scale Logical or numeric. If TRUE (default FALSE), scale columns to standard deviation one.
#'
#' @return A data.frame with individuals as rows and SNPs as columns (numeric 0/1/2, or centered/scaled values).
#'
#' @examples
#' \dontrun{
#' df <- genoToDF(nelore_imputed, center = TRUE, scale = TRUE)
#' head(df[, 1:5])
#' }
#' @export
genoToDF <- function(object, center = FALSE, scale = FALSE) {
  if (!inherits(object, "SNPDataLong")) {
    stop("Input object must be of class SNPDataLong.")
  }

  snpsum <- col.summary(object = object@geno)
  mono <- check.snp.monomorf(snpsum)
  object <- Subset(object = object, index = mono, margin = 2, keep = FALSE)

  geno_matrix <- as(object@geno, "numeric")
  geno_df <- as.data.frame(geno_matrix)

  rownames(geno_df) <- rownames(object@geno)
  colnames(geno_df) <- colnames(object@geno)

  if (isTRUE(center) || isTRUE(scale) || is.numeric(center) || is.numeric(scale)) {
    cat("Applying centering and/or scaling to SNP columns...\n")
    geno_df <- as.data.frame(scale(geno_df, center = center, scale = scale))
  }

  cat("Genotype data converted to data.frame with dimensions:",
      nrow(geno_df), "x", ncol(geno_df), "\n")

  return(geno_df)
}

#' Run PCA and Anticlustering on SNPDataLong
#'
#' Converts a SNPDataLong object to a data.frame, runs PCA, and performs anticlustering grouping.
#'
#' @param object An object of class SNPDataLong.
#' @param K Number of groups for anticlustering.
#' @param n_pcs Number of top principal components to use (default: 20).
#' @param center Logical or numeric. Center columns before PCA (default: TRUE).
#' @param scale Logical or numeric. Scale columns before PCA (default: TRUE).
#'
#' @return A list with:
#' - groups: vector with group assignments.
#' - pca: the PCA result object (prcomp).
#' - pcs: matrix of top PCs used in anticlustering.
#'
#' @examples
#' \dontrun{
#' res <- runAnticlusteringPCA(nelore_imputed, K = 2, n_pcs = 20)
#' table(res$groups)
#' }
#' @export
runAnticlusteringPCA <- function(object, K = 2, n_pcs = 20, center = TRUE, scale = TRUE) {
  if (!inherits(object, "SNPDataLong")) {
    stop("Input object must be of class SNPDataLong.")
  }

  geno_df <- genoToDF(object, center = center, scale = scale)
  cat("Genotype data frame created for PCA.\n")

  cat("Running PCA...\n")
  pca_res <- stats::prcomp(geno_df, center = FALSE, scale. = FALSE)

  top_pcs <- pca_res$x[, seq_len(n_pcs)]
  cat("Top", n_pcs, "PCs extracted.\n")

  cat("Running anticlustering with", K, "groups...\n")
  groups <- anticlust::fast_anticlustering(top_pcs, K = K)
  cat("Anticlustering completed. Groups assigned.\n")

  return(list(
    groups = groups,
    pca = pca_res,
    pcs = top_pcs
  ))
}

#' Plot PCA groups from anticlustering result
#'
#' @param pca_res A prcomp object.
#' @param groups A factor or vector of group assignments.
#' @param pcs Vector of length 2 indicating which PCs to plot (default: c(1, 2)).
#' @param filename Optional. If provided, saves plot to this file (e.g., "antic.png").
#'
#' @return A ggplot object (also prints to screen).
#'
#' @examples
#' \dontrun{
#' res <- runAnticlusteringPCA(nelore_imputed, K = 2, n_pcs = 20)
#' plotPCAgroups(res$pca, res$groups)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal theme element_rect ggsave
#' @export
plotPCAgroups <- function(pca_res, groups, pcs = c(1, 2), filename = NULL) {
  explained_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  pc1_var <- round(100 * explained_var[pcs[1]], 2)
  pc2_var <- round(100 * explained_var[pcs[2]], 2)

  pc_df <- data.frame(
    PC1 = pca_res$x[, pcs[1]],
    PC2 = pca_res$x[, pcs[2]],
    Group = as.factor(groups)
  )

  p <- ggplot2::ggplot(pc_df, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::labs(
      title = "PCA plot colored by Anticlustering Group",
      x = paste0("PC", pcs[1], " (", pc1_var, "%)"),
      y = paste0("PC", pcs[2], " (", pc2_var, "%)")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))

  if (!is.null(filename)) {
    ggplot2::ggsave(filename, p, width = 7, height = 5, dpi = 300)
    cat("Plot saved to:", filename, "\n")
  } else {
    print(p)
  }

  return(p)
}
