#' Faster row-bind for SnpMatrix objects with differing columns
#'
#' Combines multiple SnpMatrix objects by rows, automatically handling differing SNP columns, optimized for large matrices.
#'
#' @param ... One or more SnpMatrix objects.
#'
#' @return A single SnpMatrix object with all rows combined.
#' @examples
#' \dontrun{
#' combined <- rbindSnpFlexible(brangus_geno, batch_BM@geno)
#' }
#' @importClassesFrom snpStats SnpMatrix
#' @export
rbindSnpFlexible <- function(...) {
  matrices <- list(...)
  
  if (length(matrices) < 2) {
    stop("❌ At least two SnpMatrix objects are required.")
  }
  
  # Check all are SnpMatrix
  are_snp <- vapply(matrices, function(x) inherits(x, "SnpMatrix"), logical(1))
  if (!all(are_snp)) {
    stop("❌ All inputs must be SnpMatrix objects.")
  }
  
  # Union of columns
  all_snps <- unique(unlist(lapply(matrices, colnames)))
  n_total_rows <- sum(sapply(matrices, nrow))
  n_total_cols <- length(all_snps)
  
  # Create final big raw matrix
  big_raw <- matrix(as.raw(0), nrow = n_total_rows, ncol = n_total_cols)
  colnames(big_raw) <- all_snps
  
  all_row_names <- character(n_total_rows)
  row_ptr <- 1
  
  for (mat in matrices) {
    n_rows <- nrow(mat)
    rows_idx <- row_ptr:(row_ptr + n_rows - 1)
    
    common_snps <- intersect(colnames(mat), all_snps)
    cols_idx <- match(common_snps, all_snps)
    
    # Fill directly
    big_raw[rows_idx, cols_idx] <- mat[, common_snps, drop = FALSE]@.Data
    
    all_row_names[rows_idx] <- rownames(mat)
    row_ptr <- row_ptr + n_rows
  }
  
  # Set row names after loop
  dimnames(big_raw) <- list(all_row_names, all_snps)
  
  # Create SnpMatrix
  combined <- new("SnpMatrix", big_raw)
  
  return(combined)
}

#' Subset method for SNPDataLong
#'
#' @param object A SNPDataLong object.
#' @param index Character vector with row (individual) or column (SNP) names to filter.
#' @param margin Integer: 1 = rows (individuals), 2 = columns (SNPs).
#' @param keep Logical: TRUE to keep specified levels, FALSE to discard them.
#' @return A new SNPDataLong object, subsetted accordingly.
#' @export
setGeneric("Subset", function(object, index, margin = 1, keep = TRUE) standardGeneric("Subset"))

setMethod("Subset", "SNPDataLong",
          function(object, index, margin = 1, keep = TRUE) {
            
            if (!inherits(object@geno, "SnpMatrix")) {
              stop("The 'geno' slot must be a valid SnpMatrix.")
            }
            if (!is.character(index)) {
              stop("'index' must be a character vector with row or column names.")
            }
            if (!(margin %in% c(1, 2))) {
              stop("'margin' must be 1 (rows) or 2 (columns).")
            }
            if (!is.logical(keep) || length(keep) != 1) {
              stop("'keep' must be a single logical value (TRUE or FALSE).")
            }
            
            rn <- rownames(object@geno)
            cn <- colnames(object@geno)
            
            if (margin == 1) {
              # Subset rows (individuals)
              to_select <- rn %in% index
              if (!keep) {
                to_select <- !to_select
              }
              
              new_geno <- object@geno[to_select, , drop = FALSE]
              new_xref <- object@xref_path[to_select]
              
              # Atualiza path para conter apenas paths únicos restantes
              new_path <- paste(unique(new_xref), collapse = ";")
              
              # Atualiza slots
              object@geno <- new_geno
              object@xref_path <- new_xref
              object@path <- new_path
              
            } else {
              # Subset columns (SNPs)
              to_select <- cn %in% index
              if (!keep) {
                to_select <- !to_select
              }
              
              new_geno <- object@geno[, to_select, drop = FALSE]
              new_map <- object@map[object@map$Name %in% colnames(new_geno), , drop = FALSE]
              
              # Atualiza slots
              object@geno <- new_geno
              object@map <- new_map
            }
            
            return(object)
          })
