#' @title Run FImpute from a FImputeRunner object
#'
#' @description
#' This function runs the external FImpute software using a `FImputeRunner` object,
#' ensuring that all required input files are present and the results are imported.
#'
#' @param object An object of class `FImputeRunner`.
#' @param verbose Logical. If TRUE (default), FImpute output will be printed to the console.
#'
#' @return An updated `FImputeRunner` object with the `results` slot populated (SnpMatrix).
#' @examples
#' \dontrun{
#' # Example: Running FImpute from a FImputeRunner object
#'
#' path_fimpute <- "fimpute_run_example"
#' param_file <- file.path(path_fimpute, "fimpute.par")
#' fimpute_exec <- "FImpute3"  # assuming it is in PATH
#'
#' export_obj <- new("FImputeExport",
#'                   geno = geno_obj@geno,
#'                   map = geno_obj@map,
#'                   path = path_fimpute)
#'
#' runner <- new("FImputeRunner",
#'               export = export_obj,
#'               par_file = param_file,
#'               exec_path = fimpute_exec)
#'
#' runner <- runFImpute(runner, verbose = TRUE)
#' head(runner@results)
#' }
#' @export
setGeneric("runFImpute", function(object, verbose = TRUE) standardGeneric("runFImpute"))

#' @rdname runFImpute
#' @export
setMethod("runFImpute", "FImputeRunner", function(object, verbose = TRUE) {
  dir <- object@export@path
  par_file <- object@par_file
  exec <- object@exec_path

  if (!dir.exists(dir)) stop("Working directory does not exist: ", dir)
  if (!file.exists(par_file)) stop("Parameter file not found: ", par_file)

  resolved_exec <- Sys.which(exec)
  if (resolved_exec == "") {
    stop("FImpute executable not found in PATH: ", exec)
  }

  if (.Platform$OS.type == "unix" && file.access(resolved_exec, 1) != 0) {
    stop("FImpute executable exists but is not executable: ", resolved_exec)
  }

  destino_par <- file.path(dir, "fimpute.par")
  if (normalizePath(par_file) != normalizePath(destino_par)) {
    file.copy(par_file, destino_par, overwrite = TRUE)
    message("✓ Parameter file copied to working directory: ", destino_par)
  } else {
    message("✓ Parameter file already in working directory, no copy needed.")
  }

  output_dir <- file.path(dir, "output_fimpute")
  if (dir.exists(output_dir)) {
    message("⚠️ Removing previous output directory: ", output_dir)
    unlink(output_dir, recursive = TRUE, force = TRUE)
  }

  message("🧬 Running FImpute...")

  command <- paste("cd", shQuote(dir), "&&", shQuote(resolved_exec), "fimpute.par")
  status <- system(
    command,
    intern = FALSE,
    ignore.stdout = !verbose,
    ignore.stderr = !verbose
  )

  if (status != 0) {
    stop("❌ FImpute execution failed. Please check the executable and input files.")
  }

  if (!dir.exists(output_dir)) {
    stop("❌ Output directory 'output_fimpute' was not created. FImpute may have failed silently.")
  }

  if (!exists("read.fimpute", mode = "function")) {
    stop("The function 'read.fimpute()' must be defined and available in the current environment.")
  }

  res <- read.fimpute(file = output_dir)

  # if (!inherits(res, "SnpMatrix")) {
    # stop("The function 'read.fimpute()' did not return a SnpMatrix object. Please verify its behavior.")
  # }

  message("✔ Results successfully read and stored in 'results' slot of FImputeRunner.")
  object <- res
  # object@results <- res

  return(object)
})

#' Read imputed genotypes from FImpute output and return SNPDataLong object
#'
#' Reads imputed genotypes and SNP information from FImpute output,
#' builds a SnpMatrix and a corresponding map, and returns an SNPDataLong object.
#'
#' @param file Character. Path to the FImpute output directory (usually "output_fimpute").
#'
#' @return An object of class \code{SNPDataLong} containing the imputed genotypes and SNP map.
#'
#' @examples
#' \dontrun{
#' snp_long <- read.fimpute("output_fimpute")
#' }
#' @importClassesFrom snpStats SnpMatrix
#' @export
read.fimpute <- function(file) {
  # Define file paths
  genotype_file <- file.path(file, "genotypes_imp.txt")
  snp_info_file <- file.path(file, "snp_info.txt")
  
  # Read SNP map
  map <- read.table(snp_info_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Check required columns
  required_cols <- c("SNPID", "Chr", "BPPos")
  if (!all(required_cols %in% colnames(map))) {
    stop("snp_info.txt must contain columns: SNP_ID, Chromosome, Position")
  }
  
  # Extract SNP names
  snp_names <- map$SNPID
  
  # Read genotype data
  temp <- read.table(genotype_file, colClasses = "character", header = TRUE, stringsAsFactors = FALSE)
  
  nr <- nrow(temp)
  nc <- nchar(temp[1, 3])
  
  cat("Number of individuals: ", nr, "\n")
  cat("Number of SNPs: ", nc, "\n")
  
  # Split genotype strings into matrix
  nTemp <- matrix(0L, nr, nc)
  for (i in seq_len(nr)) {
    if (i %% 500 == 0) cat("  ", i, " lines processed\n")
    nTemp[i, ] <- as.integer(strsplit(temp[i, 3], "")[[1]])
  }
  
  # Convert codes
  nTemp[nTemp == 2] <- 3L
  nTemp[nTemp == 1] <- 2L
  nTemp[nTemp == 0] <- 1L
  
  # Set row and column names
  rownames(nTemp) <- temp[, 1]
  colnames(nTemp) <- snp_names
  
  # Create SnpMatrix object
  snp_matrix <- new("SnpMatrix", nTemp)
  
  # Build SNP map data.frame to match SNPDataLong
  map_out <- data.frame(
    Name = map$SNPID,
    Chromosome = map$Chr,
    Position = map$BPPos,
    stringsAsFactors = FALSE
  )
  
  # Create SNPDataLong object, including path
  snp_long <- new("SNPDataLong", geno = snp_matrix, map = map_out, path = file)
  
  return(snp_long)
}
