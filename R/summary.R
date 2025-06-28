#' Summary for SNPDataLong objects
#'
#' Provides a detailed summary of an \code{SNPDataLong} object, including sample and SNP counts,
#' proportion of missing data, and SNP distribution by chromosome if mapping information is available.
#'
#' @param object An object of class \code{SNPDataLong}.
#' @return Prints a summary to the console. Returns \code{NULL} (invisible).
#' @export
setMethod("summary", "SNPDataLong", function(object, ...) {
  cat("Summary of SNPDataLong object\n")
  cat("-----------------------------\n")

  # Basic validations
  if (!inherits(object@geno, "SnpMatrix")) {
    cat("Warning: Slot 'geno' is not a valid SnpMatrix.\n")
    return(invisible(NULL))
  }

  if (!is.data.frame(object@map)) {
    cat("Warning: Slot 'map' is not a data.frame.\n")
    return(invisible(NULL))
  }

  if (nrow(object@geno) == 0 || ncol(object@geno) == 0) {
    cat("Empty object: no individuals or SNPs.\n")
    return(invisible(NULL))
  }

  # General information
  n_ind <- nrow(object@geno)
  n_snp <- ncol(object@geno)
  cat("Individuals :", n_ind, "\n")
  cat("SNPs        :", n_snp, "\n\n")

  # Missing data summary
  n_total <- n_ind * n_snp
  n_na <- sum(is.na(object@geno))
  pct_na <- round(100 * n_na / n_total, 2)
  cat("Missing data (NA):\n")
  cat(" - Total     :", n_na, "of", n_total, "\n")
  cat(" - Proportion:", pct_na, "%\n\n")

  # Chromosome analysis (if map columns are available)
  if ("Name" %in% colnames(object@map) && "Chromosome" %in% colnames(object@map)) {
    idx <- match(colnames(object@geno), object@map$Name)
    chr_info <- object@map$Chromosome[idx]

    if (any(is.na(chr_info))) {
      cat("Warning: Some SNPs were not found in the map and will be ignored in chromosome counts.\n")
      chr_info <- chr_info[!is.na(chr_info)]
    }

    if (length(chr_info) > 0) {
      cat("Distribution of SNPs by chromosome:\n")
      print(table(chr_info))
      cat("\n")

      # Count SNPs with missing data by chromosome
      snp_na_count <- colSums(is.na(object@geno))
      chr_na <- chr_info
      names(chr_na) <- colnames(object@geno)[idx]
      chr_na_tab <- tapply(snp_na_count > 0, chr_na, sum)
      cat("SNPs with missing data by chromosome:\n")
      print(chr_na_tab)
    } else {
      cat("No chromosomes found after aligning with the map.\n")
    }
  } else {
    cat("Map does not contain expected columns ('Name', 'Chromosome'); chromosome summary omitted.\n")
  }

  invisible(NULL)
})
