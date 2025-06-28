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
#' @details
#' The function performs the following steps internally:
#' \enumerate{
#'   \item Computes the union of all SNPs across input objects.
#'   \item Fills missing SNP columns in each genotype matrix with NA-coded columns.
#'   \item Combines genotype matrices by rows (individuals).
#'   \item Merges marker maps, removing duplicates (retaining the first occurrence).
#'   \item Creates and returns a new \code{SNPDataLong} object with the combined data.
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' combined_data <- combinarSNPData(list(snp_data1, snp_data2, snp_data3))
#' }
#'
#' @export
combinarSNPData <- function(lista) {
  stopifnot(length(lista) > 0)

  # Get the union of all SNP names across all objects
  snps_all <- Reduce(union, lapply(lista, function(x) colnames(x@geno)))

  # Ensure each genotype matrix has all SNPs (fill with NAs where missing)
  geno_list <- lapply(lista, function(x) {
    geno <- x@geno
    missing_snps <- setdiff(snps_all, colnames(geno))
    if (length(missing_snps) > 0) {
      # Create a SnpMatrix of NAs for missing SNPs
      na_block <- new("SnpMatrix", matrix(as.raw(0), nrow = nrow(geno), ncol = length(missing_snps)))
      colnames(na_block) <- missing_snps
      geno <- cbind(geno, na_block)
    }
    # Reorder columns to maintain consistent SNP order
    geno[, snps_all, drop = FALSE]
  })

  # Combine genotype matrices row-wise
  geno_comb <- do.call(rbind, geno_list)

  # Combine marker maps and remove duplicate SNPs, keeping the first occurrence
  map_all <- do.call(rbind, lapply(lista, function(x) x@map))
  map_all <- map_all[!duplicated(map_all$Name), , drop = FALSE]
  map_final <- map_all[match(snps_all, map_all$Name), , drop = FALSE]

  # Final validation and type coercion if needed
  if (!inherits(geno_comb, "SnpMatrix")) {
    geno_comb <- methods::as(geno_comb, "SnpMatrix")
  }

  # Create and return the combined SNPDataLong object
  new("SNPDataLong",
      geno = geno_comb,
      map  = map_final,
      path = paste(sapply(lista, function(x) x@path), collapse = ";"))
}
