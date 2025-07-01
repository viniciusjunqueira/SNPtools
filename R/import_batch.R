#' Import multiple genotype datasets from a list of configurations
#'
#' This function iterates over a list of configuration lists (each specifying parameters such as path, fields, separators, etc.), 
#' imports each genotype dataset using \code{getGeno()}, and then combines them into a single \code{SNPDataLong} object.
#'
#' @param config_list A list of configuration lists. Each element must include at least \code{path} and \code{fields}. 
#' Optional elements include \code{codes}, \code{threshold}, \code{sep}, \code{skip}, and \code{verbose}.
#'
#' @return A unified \code{SNPDataLong} object containing combined genotype data from all configurations.
#'
#' @examples
#' \dontrun{
#' configs <- list(
#'   list(path = "panel1", fields = list(sample = 2, snp = 1, allele1 = 7, allele2 = 8, confidence = 9)),
#'   list(path = "panel2", fields = list(sample = 2, snp = 1, allele1 = 7, allele2 = 8, confidence = 9), threshold = 0.10)
#' )
#' combined_data <- import_geno_list(configs)
#' }
#'
#' @export
import_geno_list <- function(config_list) {
  stopifnot(is.list(config_list))
  
  # Import each genotype dataset using configuration parameters
  results <- lapply(config_list, function(cfg) {
    message("📄 Importing genotypes from: ", cfg$path)
    
    geno_obj <- getGeno(
      path = cfg$path,
      fields = cfg$fields,
      codes = if (!is.null(cfg$codes)) cfg$codes else c("A", "B"),
      threshold = if (!is.null(cfg$threshold)) cfg$threshold else 0.15,
      sep = if (!is.null(cfg$sep)) cfg$sep else "\t",
      skip = if (!is.null(cfg$skip)) cfg$skip else 0,
      verbose = if (!is.null(cfg$verbose)) cfg$verbose else TRUE
    )
    
    if (is.null(geno_obj)) {
      message("⚠️ Skipping file '", cfg$path, "' because it returned NULL.")
      return(NULL)
    }
    
    return(geno_obj)
  })
  
  # Remove any NULL results (skip)
  results <- Filter(Negate(is.null), results)
  
  # Check if at least one valid object remains
  if (length(results) == 0) {
    stop("❌ No valid genotype datasets were imported. All configurations returned NULL.")
  }
  
  # Combine all imported genotype datasets into one SNPDataLong object
  combinarSNPData(results)
}
