#' @title Save genotype and map files in FImpute format
#'
#' @description
#' S4 method to export genotype (`.gen`), map (`.map`), and parameter (`data.par`) files compatible with the [FImpute](https://www.aps.uoguelph.ca/~msargol/fimpute/) software.
#'
#' @param object An object of class `FImputeExport` or `SNPDataLong`.
#' @param path Character. Output directory where files will be written (only for `SNPDataLong` method; default = `"fimpute_run"`).
#' @param ... Further arguments passed to methods.
#'
#' @return No return value. Files are saved to disk.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("snpStats", quietly = TRUE)) {
#'   mat <- matrix(sample(c(0L, 1L, 2L), 50, replace = TRUE), nrow = 5)
#'   colnames(mat) <- paste0("snp", 1:10)
#'   rownames(mat) <- paste0("ind", 1:5)
#'
#'   sm <- new("SnpMatrix", data = as.raw(mat))
#'   map <- data.frame(Name = colnames(mat), Chromosome = 1, Position = 1:10)
#'   x <- new("SNPDataLong", geno = sm, map = map)
#'
#'   saveFImpute(x, path = tempdir())
#' }
#' }
#'
#' @export
setGeneric("saveFImpute", function(object, ...) standardGeneric("saveFImpute"))

#' @rdname saveFImpute
#' @export
setMethod("saveFImpute", "FImputeExport", function(object, ...) {
  save_fimpute_raw(object@geno, object@map, object@path)
})

#' @rdname saveFImpute
#' @export
setMethod("saveFImpute", "SNPDataLong", function(object, path = "fimpute_run") {
  if (!requireNamespace("snpStats", quietly = TRUE)) {
    stop("The 'snpStats' package is required. Please install it with install.packages('snpStats').")
  }

  if (!inherits(object@geno, "SnpMatrix")) {
    stop("The 'geno' slot of the object must be of class 'SnpMatrix'.")
  }

  if (!is.data.frame(object@map)) {
    stop("The 'map' slot of the object must be a data.frame.")
  }

  fimpute_export <- new("FImputeExport",
                        geno = object@geno,
                        map = object@map,
                        path = path)

  saveFImpute(fimpute_export)
})

#' @title Export genotypes and map using basic arguments
#'
#' @description
#' Convenience function to export FImpute files directly from a `SnpMatrix` and map `data.frame`.
#'
#' @param geno A `SnpMatrix` object (from the `snpStats` package).
#' @param map A `data.frame` with columns 'Name', 'Chromosome', and 'Position'.
#' @param path Path where the files will be saved.
#'
#' @return No return value. Files are saved to disk.
#' @export
saveFImputeRaw <- function(geno, map, path) {
  export <- new("FImputeExport", geno = geno, map = map, path = path)
  saveFImpute(export)
}

#' Internal function: writes files in FImpute format (.gen, .map, data.par)
#'
#' @noRd
save_fimpute_raw <- function(genotype, map, caminho) {
  if (!requireNamespace("snpStats", quietly = TRUE)) {
    stop("The 'snpStats' package is required. Please install it with install.packages('snpStats').")
  }

  if (!inherits(genotype, "SnpMatrix")) {
    stop("The 'genotype' object must be of class 'SnpMatrix'.")
  }

  if (!is.data.frame(map)) {
    stop("The 'map' argument must be a data.frame.")
  }

  required_cols <- c("Name", "Chromosome", "Position")
  missing_cols <- setdiff(required_cols, colnames(map))
  if (length(missing_cols) > 0) {
    stop("The map is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.character(caminho) || length(caminho) != 1) {
    stop("'caminho' must be a single character string indicating the output directory.")
  }

  if (!dir.exists(caminho)) {
    message("Creating output directory: ", caminho)
    dir.create(caminho, recursive = TRUE)
  } else {
    existing_files <- list.files(caminho, pattern = "data\\.(gen|map|par)$", full.names = TRUE)
    if (length(existing_files) > 0) {
      warning("The following files will be overwritten:\n  ",
              paste(basename(existing_files), collapse = "\n  "))
    }
  }

  ## Write .gen file
  gen_file <- file.path(caminho, "data.gen")
  con <- file(gen_file, "wt")

  smp <- rownames(genotype)
  if (is.null(smp)) stop("The 'genotype' object must have row names (individual IDs).")

  nC <- max(nchar(smp), na.rm = TRUE)
  writeLines("ID    Chip                   Call...", con)

  for (i in seq_len(nrow(genotype))) {
    temp <- as(genotype[i, ], "character")
    temp[temp == "A/A"] <- "0"
    temp[temp == "A/B"] <- "1"
    temp[temp == "B/B"] <- "2"
    temp[is.na(temp) | temp == "NA"] <- "5"
    linha <- paste(format(smp[i], width = nC, justify = "left"), "1", paste(temp, collapse = ""))
    writeLines(linha, con)
  }

  close(con)
  message("✔ File written: ", gen_file)

  ## Write .map file
  map_file <- file.path(caminho, "data.map")
  write("SNP_ID           Chr      Pos   Chip1", file = map_file)

  map <- map[map$Name %in% colnames(genotype), ]
  if (nrow(map) == 0) {
    warning("No SNPs in the map matched the columns in the genotype object.")
  }

  map_out <- data.frame(map[, required_cols, drop = FALSE], 1:nrow(map))
  write.table(
    map_out,
    file = map_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = " ",
    append = TRUE
  )
  message("✔ File written: ", map_file)

  ## Write data.par file
  par_file <- file.path(caminho, "fimpute.par")
  writeLines(c(
    'title="FImpute imputation";',
    'genotype_file="data.gen";',
    'snp_info_file="data.map";',
    'parentage_test /ert_mm=0.02 /remove_conflict /find_match_cnflt /find_match_mp /find_match_ugp;',
    'output_folder="output_fimpute";',
    'save_genotype;',
    'njob=24;'
  ), con = par_file)
  message("✔ File written: ", par_file)

  invisible(NULL)
}
