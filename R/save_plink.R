#' Save SNPDataLong object to PLINK format
#'
#' Saves genotype and map data from an SNPDataLong object in PLINK format (.ped/.map and optionally binary files).
#'
#' @param object An object of class SNPDataLong.
#' @param path Character. Directory where files will be saved.
#' @param name Character. Base name for PLINK output files.
#' @param run_plink Logical. If TRUE (default), runs PLINK1 to convert to binary files. If FALSE, only .ped and .map files are saved.
#' @param chunk_size Integer. Number of individuals per chunk for writing .ped file (default: 1000).
#'
#' @return No return value. Files are saved to disk.
#'
#' @examples
#' \dontrun{
#' savePlink(genotypes_qc, path = "plink_out", name = "nelore_qc", run_plink = TRUE, chunk_size = 2000)
#' }
#' @importFrom utils write.table
#' @export
savePlink <- function(object, path = "plink_out", name = "plink_data", run_plink = TRUE, chunk_size = 1000) {
  if (!inherits(object, "SNPDataLong")) {
    stop("Input object must be of class SNPDataLong.")
  }

  qc_header("Saving Files in Plink Format")

  geno <- object@geno
  map <- object@map
  n_ind <- nrow(geno)

  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    cat("Created output directory:", path, "\n")
  }

  ## ----- PED file -----
  ped_file <- file.path(path, paste0(name, ".ped"))
  smp <- rownames(geno)

  cat("Writing .ped file in chunks...\n")
  con <- file(ped_file, "wt")

  for (start in seq(1, n_ind, by = chunk_size)) {
    end <- min(start + chunk_size - 1, n_ind)
    idx <- start:end

    # Convert block of individuals
    geno_chr <- as(geno[idx, , drop = FALSE], "character")
    geno_chr <- gsub("/", " ", geno_chr)
    geno_chr[is.na(geno_chr)] <- "0 0"

    # Build lines for this chunk
    lines <- vapply(seq_along(idx), function(i) {
      paste("NA", smp[idx[i]], "0 0 -9 -9", paste(geno_chr[i, ], collapse = " "))
    }, character(1L))

    writeLines(lines, con)
    cat(sprintf(" Wrote individuals %d to %d\n", start, end))
  }

  close(con)
  cat(".ped file written:", ped_file, "\n")

  ## ----- MAP file -----
  cat("Writing .map file...\n")
  map_file <- file.path(path, paste0(name, ".map"))
  map_out <- data.frame(
    Chromosome = map$Chromosome,
    SNP_ID = map$Name,
    # Genetic_distance = 0,
    Position = map$Position,
    stringsAsFactors = FALSE
  )
  utils::write.table(map_out, map_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
  cat(".map file written:", map_file, "\n")

  ## ----- Optionally run PLINK -----
  if (run_plink) {
    cat("Running PLINK to generate binary files...\n")
    cmd <- paste("cd", shQuote(path), "&& plink1 --file", shQuote(name), "--map3 --out", shQuote(name), "--make-bed --noweb")
    status <- system(cmd)

    if (status == 0) {
      cat("PLINK binary files created successfully.\n")
    } else {
      cat("PLINK execution failed. Please check your installation and logs.\n")
    }
  } else {
    cat("Skipping PLINK binary conversion as requested.\n")
  }

  cat("All done!\n")
}
