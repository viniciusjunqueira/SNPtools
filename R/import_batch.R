#' Import multiple genotype datasets from a list of configurations
#'
#' @description
#' Reads and imports multiple genotype datasets specified in a list of configurations.
#' Each configuration must include the path to the genotype data and information on field mapping.
#' Optionally, you can also specify codes, quality threshold, separator, lines to skip, and a subset of IDs to retain.
#' The function automatically fills the `xref_path` slot per individual and combines maps into a single data.frame,
#' adding a `SourcePath` column indicating their origin and removing duplicated SNP rows (by Name).
#' Prints progress messages indicating the current path being loaded (with counter).
#'
#' @param config_list A list of configuration lists. Each element should contain:
#'   - `path` (character): Path to the genotype file or folder.
#'   - `fields` (list): Named list defining the columns (e.g., SNP ID, sample ID, alleles, confidence).
#'   - `codes` (character vector, optional): Allele codes (default is c("A", "B")).
#'   - `threshold` (numeric, optional): Maximum allowed missingness or confidence threshold (default 0.15).
#'   - `sep` (character, optional): Field separator in the input file (default tab "\t").
#'   - `skip` (integer, optional): Number of lines to skip at the beginning of the file (default 0).
#'   - `verbose` (logical, optional): Whether to print detailed messages (default TRUE).
#'   - `subset` (character vector, optional): Vector of sample IDs to retain after import.
#'
#' @return An object of class `SNPDataLong` containing:
#'   - Combined genotype matrix (`geno`).
#'   - Combined map (`map`) as a single data.frame with `SourcePath` column and without duplicated rows.
#'   - Combined `xref_path` vector (one entry per individual).
#'   - `path` slot as a semicolon-separated string of all input dataset paths.
#'
#' @export
import_geno_list <- function(config_list) {
  # Validate input
  stopifnot(is.list(config_list))
  
  # Total number of paths to process
  total_paths <- length(config_list)
  counter <- 1
  
  # Process each configuration
  results <- lapply(config_list, function(cfg) {
    # Print progress message with counter
    message(sprintf("📄 Loading path %d/%d: %s", counter, total_paths, cfg$path))
    
    # Call getGeno with parameters from config
    geno_obj <- getGeno(
      path      = cfg$path,
      fields    = cfg$fields,
      codes     = if (!is.null(cfg$codes)) cfg$codes else c("A", "B"),
      threshold = if (!is.null(cfg$threshold)) cfg$threshold else 0.15,
      sep       = if (!is.null(cfg$sep)) cfg$sep else "\t",
      skip      = if (!is.null(cfg$skip)) cfg$skip else 0,
      verbose   = if (!is.null(cfg$verbose)) cfg$verbose else TRUE
    )
    
    # If failed, warn and continue
    if (is.null(geno_obj)) {
      warning("⚠️ No data returned for path: ", cfg$path)
      counter <<- counter + 1
      return(NULL)
    }
    
    # Apply subset filter if provided
    if (!is.null(cfg$subset)) {
      subset_ids <- cfg$subset
      geno_obj@geno <- geno_obj@geno[rownames(geno_obj@geno) %in% subset_ids, , drop = FALSE]
    }
    
    # Update path slot
    geno_obj@path <- cfg$path
    
    # Fill xref_path with path per individual
    n_ind <- nrow(geno_obj@geno)
    if (n_ind > 0) {
      geno_obj@xref_path <- rep(cfg$path, n_ind)
    }
    
    # Check map validity
    if (is.null(geno_obj@map) || is.character(geno_obj@map)) {
      stop("❌ Invalid map in configuration with path: ", cfg$path, 
           ". Check that getGeno returns a valid data.frame for map.")
    }
    
    # Add source path column in map
    geno_obj@map$SourcePath <- cfg$path
    
    counter <<- counter + 1
    return(geno_obj)
  })
  
  # Remove failed or NULL results
  results <- Filter(Negate(is.null), results)
  
  # Check if at least one valid object remains
  if (length(results) == 0) {
    stop("❌ No valid genotype datasets were imported. All configurations returned NULL.")
  }
  
  # Combine genotype matrices (custom function you already have)
  combined <- combinarSNPData(results)
  
  # Combine xref paths
  all_xref_paths <- unlist(lapply(results, function(x) x@xref_path), use.names = FALSE)
  combined@xref_path <- if (is.null(all_xref_paths) || length(all_xref_paths) == 0) character(0) else all_xref_paths
  
  # Combine all maps
  all_maps_df <- do.call(rbind, lapply(results, function(x) x@map))
  
  # Remove duplicated rows by Name
  dup_idx <- duplicated(all_maps_df[, c("Name")])
  all_maps_df <- all_maps_df[!dup_idx, , drop = FALSE]
  
  # Assign merged map to combined object
  combined@map <- all_maps_df
  
  # Update combined path slot with semicolon-separated paths
  paths <- vapply(results, function(x) x@path, character(1))
  combined@path <- paste(paths, collapse = ";")
  
  return(combined)
}
