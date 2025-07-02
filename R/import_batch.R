#' Import multiple genotype datasets from a list of configurations
#'
#' @param config_list List of configurations (see documentation).
#' @return Combined SNPDataLong object with xref_path slot correctly filled.
#'
#' @export
import_geno_list <- function(config_list) {
  stopifnot(is.list(config_list))
  
  xref_path_list <- list()  # To store paths per individual

  # Import each genotype dataset using configuration parameters
  results <- lapply(config_list, function(cfg) {
    geno_obj <- getGeno(
      path      = cfg$path,
      fields    = cfg$fields,
      codes     = if (!is.null(cfg$codes)) cfg$codes else c("A", "B"),
      threshold = if (!is.null(cfg$threshold)) cfg$threshold else 0.15,
      sep       = if (!is.null(cfg$sep)) cfg$sep else "\t",
      skip      = if (!is.null(cfg$skip)) cfg$skip else 0,
      verbose   = if (!is.null(cfg$verbose)) cfg$verbose else TRUE
    )
    
    if (is.null(geno_obj)) {
      return(NULL)
    }
    
    # Apply subset if present
    if (!is.null(cfg$subset)) {
      subset_ids <- cfg$subset
      geno_obj@geno <- geno_obj@geno[rownames(geno_obj@geno) %in% subset_ids, , drop = FALSE]
    }
    
    # Update path slot
    geno_obj@path <- cfg$path

    # Create xref path vector
    n_ind <- nrow(geno_obj@geno)
    if (n_ind > 0) {
      geno_obj@xref_path <- rep(cfg$path, n_ind)
    }

    return(geno_obj)
  })
  
  # Remove NULL results
  results <- Filter(Negate(is.null), results)
  
  # Check if at least one valid object remains
  if (length(results) == 0) {
    stop("❌ No valid genotype datasets were imported. All configurations returned NULL.")
  }
  
  # Combine imported genotype datasets
  combined <- combinarSNPData(results)
  
  # After combination, force combined xref_path
  all_xref_paths <- unlist(lapply(results, function(x) x@xref_path), use.names = FALSE)
  if (is.null(all_xref_paths) || length(all_xref_paths) == 0) {
    all_xref_paths <- character(0)
  }
  combined@xref_path <- all_xref_paths

  return(combined)
}
