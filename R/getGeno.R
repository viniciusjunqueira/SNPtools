#' Flexible and efficient genotype file reading with autodetection using fread
#'
#' This generic and method allow flexible import of SNP genotype data from Illumina FinalReport files,
#' supporting fast initial column detection using \code{data.table::fread}, followed by full genotype
#' matrix construction via \code{snpStats::read.snps.long}.
#'
#' @param path Path to the directory containing \code{FinalReport.txt}
#' @param fields A list specifying column indices for sample, SNP, allele1, allele2, and confidence
#' @param codes A character vector with allele codes (e.g., \code{c("A", "B")})
#' @param threshold Confidence threshold for genotype calling
#' @param sep Field separator used in the files
#' @param skip Number of lines to skip at the start of the file
#' @param verbose Logical; if \code{TRUE}, displays progress messages
#' @param every Frequency of progress update (number of SNPs)
#'
#' @return An \code{SNPDataLong} object containing the genotype matrix and map, or \code{NULL} if an error occurs
#' @export
setGeneric("getGeno", function(...) standardGeneric("getGeno"))

setMethod("getGeno", signature(),
  function(path,
           fields = list(sample = 2, snp = 1, allele1 = 7, allele2 = 8, confidence = 9),
           codes = c("A", "B"),
           threshold = 0.15,
           sep = "\t",
           skip = 0,
           verbose = TRUE,
           every = NULL) {

    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("The 'data.table' package is required. Install it using install.packages('data.table').")
    }

    if (!file.exists(file.path(path, "FinalReport.txt"))) {
      warning("File FinalReport.txt not found at: ", path)
      return(NULL)
    }

    # Fast initial read using fread to detect samples and SNPs
    cols_to_read <- unique(unlist(fields[c("sample", "snp")]))
    col_select <- as.integer(cols_to_read)

    fread_result <- tryCatch({
      data.table::fread(
        file = file.path(path, "FinalReport.txt"),
        sep = sep,
        skip = skip,
        header = TRUE,
        select = col_select,
        data.table = FALSE,
        colClasses = "character"
      )
    }, error = function(e) {
      warning("Failed to read FinalReport.txt with fread: ", e$message)
      return(NULL)
    })

    if (is.null(fread_result)) return(NULL)

    sample.col <- match(fields$sample, col_select)
    snp.col    <- match(fields$snp, col_select)

    sample.id <- unique(fread_result[[sample.col]])
    snp.id    <- unique(fread_result[[snp.col]])

    if (length(sample.id) == 0) {
      warning("Sample IDs could not be determined correctly; setting default numeric row names.")
    }
    if (length(snp.id) == 0) {
      warning("SNP IDs could not be determined correctly; setting default numeric column names.")
    }

    if (is.null(every)) every <- length(snp.id)

    # Full genotype matrix reading with snpStats
    data <- tryCatch({
      snpStats::read.snps.long(
        file = file.path(path, "FinalReport.txt"),
        sample.id = sample.id,
        snp.id = snp.id,
        fields = fields,
        codes = codes,
        threshold = threshold,
        sep = sep,
        skip = skip,
        verbose = verbose,
        every = every
      )
    }, error = function(e) {
      warning("Error while running read.snps.long: ", e$message)
      return(NULL)
    })

    if (is.null(data)) return(NULL)

    # Force row and column names
    if (is.null(rownames(data))) {
      rownames(data) <- sample.id
    }
    if (is.null(colnames(data))) {
      colnames(data) <- snp.id
    }

    # Read map file
    map_file <- file.path(path, "SNP_Map.txt")
    if (!file.exists(map_file)) {
      warning("SNP_Map.txt file not found at: ", path)
      return(NULL)
    }

    map <- tryCatch({
      read.table(map_file, colClasses = "character", sep = sep, header = TRUE)
    }, error = function(e) {
      warning("Error reading SNP_Map.txt: ", e$message)
      return(NULL)
    })

    if (is.null(map)) return(NULL)

    # Select columns of interest from map
    map <- map[, c("Name", "Chromosome", "Position")]#, "GenTrain.Score")]

    new("SNPDataLong",
        geno = data,
        map  = map,
        path = path)
  }
)

#' Import and combine multiple genotype configurations
#'
#' This generic and method import genotype data from multiple configurations defined in an
#' \code{SNPImportList} object, then combine them into a single unified \code{SNPDataLong} object.
#'
#' @param object An object of class \code{SNPImportList} containing import configurations
#'
#' @return A single combined \code{SNPDataLong} object
#' @export
setGeneric("importAllGenos", function(object) standardGeneric("importAllGenos"))

setMethod("importAllGenos", "SNPImportList", function(object) {
  all_genos <- lapply(object@configs, function(cfg) {
    tryCatch(getGeno(cfg), error = function(e) {
      warning("Error in getGeno(): ", e$message)
      NULL
    })
  })

  all_genos <- Filter(Negate(is.null), all_genos)
  if (length(all_genos) == 0) stop("No genotype data was successfully imported.")
  combinarSNPData(all_genos)
})
