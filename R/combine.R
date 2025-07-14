#' Safe cbind for SnpMatrix preserving dimnames
#'
#' This function performs a column-wise binding of multiple \code{SnpMatrix} objects,
#' explicitly preserving row names and column names, avoiding unexpected "object has no names" warnings.
#'
#' @param ... SnpMatrix objects to combine (must have identical row names).
#'
#' @return A single combined \code{SnpMatrix} with preserved row and column names.
#'
#' @examples
#' \dontrun{
#' cbind_SnpMatrix(matrix1, matrix2)
#' }
#' @export
cbind_SnpMatrix <- function(...) {
  mats <- list(...)

  # Check that all matrices have identical row names
  row_names_list <- lapply(mats, rownames)
  if (!all(sapply(row_names_list, function(x) all(x == row_names_list[[1]])))) {
    stop("All matrices must have identical row names to cbind safely.")
  }

  message("Performing safe cbind on SnpMatrix objects...")

  # Perform cbind using base
  res <- do.call(base::cbind, mats)

  # Reassign dimnames explicitly
  rownames(res) <- row_names_list[[1]]
  colnames(res) <- unlist(lapply(mats, colnames), use.names = FALSE)

  message("Row and column names preserved after cbind.")

  res
}

#' Safe rbind for SnpMatrix preserving dimnames
#'
#' This function performs a row-wise binding of multiple \code{SnpMatrix} objects,
#' explicitly preserving row names and column names, avoiding unexpected "object has no names" warnings.
#'
#' @param ... SnpMatrix objects to combine (must have identical column names).
#'
#' @return A single combined \code{SnpMatrix} with preserved row and column names.
#'
#' @examples
#' \dontrun{
#' rbind_SnpMatrix(matrix1, matrix2)
#' }
#' @export
rbind_SnpMatrix <- function(...) {
  mats <- list(...)

  # Check that all matrices have identical column names
  col_names_list <- lapply(mats, colnames)
  if (!all(sapply(col_names_list, function(x) all(x == col_names_list[[1]])))) {
    stop("All matrices must have identical column names to rbind safely.")
  }

  message("Performing safe rbind on SnpMatrix objects...")

  # Perform rbind using base
  res <- do.call(base::rbind, mats)

  # Reassign dimnames explicitly
  rownames(res) <- unlist(lapply(mats, rownames), use.names = FALSE)
  colnames(res) <- col_names_list[[1]]

  message("Row and column names preserved after rbind.")

  res
}

#' Combine multiple SNPDataLong objects
#'
#' This function merges a list of \code{SNPDataLong} objects, typically representing different SNP panels
#' or datasets, into a single unified \code{SNPDataLong} object. It ensures that all genotype matrices
#' have the same set of SNPs (filling missing SNPs with NA), and merges the marker map information while
#' removing duplicate SNP entries.
#'
#' @param lista A list of \code{SNPDataLong} objects to be combined.
#'
#' @return A single \code{SNPDataLong} object containing the combined genotype matrix, merged map,
#' and a concatenated path string.
#'
#' @examples
#' \dontrun{
#' combined <- combinarSNPData(list(snp_obj1, snp_obj2, snp_obj3))
#' }
#'
#' @export
combinarSNPData <- function(lista) {
  stopifnot(length(lista) > 0)

  message("Starting SNPDataLong combination...")

  # Get the union of all SNP names across objects
  snps_all <- Reduce(union, lapply(lista, function(x) colnames(x@geno)))
  message("Unified SNP panel with ", length(snps_all), " SNPs.")

  # Ensure each genotype matrix has all SNPs (fill missing SNPs with NAs)
  geno_list <- lapply(lista, function(x) {
    geno <- x@geno
    missing_snps <- setdiff(snps_all, colnames(geno))

    if (length(missing_snps) > 0) {
      message("Adding ", length(missing_snps), " missing SNPs filled with NA for one matrix...")
      na_block <- new("SnpMatrix", matrix(as.raw(0), nrow = nrow(geno), ncol = length(missing_snps)))
      colnames(na_block) <- missing_snps
      rownames(na_block) <- rownames(geno)
      geno <- cbind_SnpMatrix(geno, na_block)
    }

    # Reorder columns
    geno <- geno[, snps_all, drop = FALSE]

    # Ensure row names exist
    if (is.null(rownames(geno)) || any(rownames(geno) == "")) {
      rownames(geno) <- sprintf("Sample_%d", seq_len(nrow(geno)))
    }

    # Ensure column names
    if (is.null(colnames(geno)) || any(colnames(geno) == "")) {
      colnames(geno) <- snps_all
    }

    geno
  })

  # Combine matrices using safe rbind
  geno_comb <- do.call(rbind_SnpMatrix, geno_list)

  # Combine marker maps and remove duplicates
  map_all <- do.call(rbind, lapply(lista, function(x) x@map))
  map_all <- map_all[!duplicated(map_all$Name), , drop = FALSE]
  map_final <- map_all[match(snps_all, map_all$Name), , drop = FALSE]

  # Validate and coerce if needed
  if (!inherits(geno_comb, "SnpMatrix")) {
    geno_comb <- methods::as(geno_comb, "SnpMatrix")
  }

  message("Combination complete. Final matrix: ", nrow(geno_comb), " samples x ", ncol(geno_comb), " SNPs.")

  # Return new object
  new("SNPDataLong",
      geno = geno_comb,
      map  = map_final,
      path = paste(sapply(lista, function(x) x@path), collapse = ";"))
}
