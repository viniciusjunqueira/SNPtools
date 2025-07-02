#' Import multiple genotype datasets from a list of configurations
#'
#' @description
#' Reads and imports multiple genotype datasets specified in a list of configurations.
#' Each configuration must include the path to the genotype data and information on field mapping.
#' Optionally, you can also specify codes, quality threshold, separator, lines to skip, and a subset of IDs to retain.
#' The function automatically fills the `xref_path` slot per individual and combines maps into a single data.frame,
#' adding a `SourcePath` column indicating their origin and removing duplicated SNP rows (by Name, Chromosome, Position).
#'
#' @param config_list A list of configuration lists. Each element should contain:
#'   - `path` (character): Path to the genotype file or folder.
#'   - `fields` (list): Named list defining the columns (e.g., SNP ID, sample ID, alleles).
#'   - `codes` (character, optional): Allele codes, default is c("A", "B").
#'   - `threshold` (numeric, optional): Maximum allowed missingness or confidence threshold (default 0.15).
#'   - `sep` (character, optional): Field separator in the input file, default tab ("\t").
#'   - `skip` (integer, optional): Number of lines to skip at the beginning of the file.
#'   - `verbose` (logical, optional): Whether to print detailed messages (default TRUE).
#'   - `subset` (character vector, optional): Vector of sample IDs to retain after import.
#'
#' @return An object of class `SNPDataLong` containing:
#'   - Combined genotype matrix (`geno`).
#'   - Combined map (`map`) as a single data.frame with `SourcePath` column and without duplicated rows.
#'   - Combined `xref_path` vector (one entry per individual).
#'   - Path string with semicolon-separated paths of all input datasets.
#'
#' @export
import_geno_list <- function(config_list) {
  stopifnot(is.list(config_list))
  
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

    # Create xref_path vector
    n_ind <- nrow(geno_obj@geno)
    if (n_ind > 0) {
      geno_obj@xref_path <- rep(cfg$path, n_ind)
    }

    # Check map validity
    if (is.null(geno_obj@map) || is.character(geno_obj@map)) {
      stop("❌ Invalid map in configuration with path: ", cfg$path, 
           ". Check that getGeno returns a valid data.frame for map.")
    }

    # Add SourcePath column to map
    geno_obj@map$SourcePath <- cfg$path

    return(geno_obj)
  })
  
  # Remove NULL results
  results <- Filter(Negate(is.null), results)
  
  # Check if at least one valid object remains
  if (length(results) == 0) {
    stop("❌ No valid genotype datasets were imported. All configurations returned NULL.")
  }
  
  # Combine genotype matrices
  combined <- combinarSNPData(results)
  
  # Combine xref paths
  all_xref_paths <- unlist(lapply(results, function(x) x@xref_path), use.names = FALSE)
  if (is.null(all_xref_paths) || length(all_xref_paths) == 0) {
    all_xref_paths <- character(0)
  }
  combined@xref_path <- all_xref_paths

  # Combine maps into a single data.frame
  all_maps_df <- do.call(rbind, lapply(results, function(x) x@map))

  # Remove duplicated rows by (Name, Chromosome, Position)
  # dup_idx <- duplicated(all_maps_df[, c("Name", "Chromosome", "Position")])
  dup_idx <- duplicated(all_maps_df[, c("Name")])
  all_maps_df <- all_maps_df[!dup_idx, , drop = FALSE]

  # Assign merged map to slot
  combined@map <- all_maps_df

  # Update path slot: concatenation of all paths
  paths <- vapply(results, function(x) x@path, character(1))
  combined@path <- paste(paths, collapse = ";")

  return(combined)
}
