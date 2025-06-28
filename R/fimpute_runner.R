#' @title Run FImpute from a FImputeRunner object
#'
#' @description
#' This function runs the external FImpute software using a `FImputeRunner` object,
#' ensuring that all required input files are present and the results are imported.
#'
#' @param object An object of class `FImputeRunner`
#'
#' @return An updated `FImputeRunner` object with the `results` slot populated.
#' @export
setGeneric("runFImpute", function(object) standardGeneric("runFImpute"))

#' @rdname runFImpute
#' @export
setMethod("runFImpute", "FImputeRunner", function(object) {
  dir <- object@export@path
  par_file <- object@par_file
  exec <- object@exec_path

  ## --- Initial validation ---
  if (!dir.exists(dir)) stop("Working directory does not exist: ", dir)
  if (!file.exists(par_file)) stop("Parameter file not found: ", par_file)
  if (!file.exists(exec)) stop("FImpute executable not found at: ", exec)

  ## --- Copy parameter file into working directory ---
  destino_par <- file.path(dir, "fimpute.par")
  file.copy(par_file, destino_par, overwrite = TRUE)
  message("✓ Parameter file copied to working directory: ", destino_par)

  ## --- Clean previous output, if exists ---
  output_dir <- file.path(dir, "output_fimpute")
  if (dir.exists(output_dir)) {
    message("⚠️ Removing previous output directory: ", output_dir)
    unlink(output_dir, recursive = TRUE, force = TRUE)
  }

  ## --- Run FImpute ---
  command <- paste("cd", shQuote(dir), "&&", shQuote(exec), "fimpute.par")
  message("🧬 Running FImpute...")
  message("🖥️  Command: ", command)

  status <- system(command, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)

  if (status != 0) {
    stop("❌ FImpute execution failed. Please check the executable and input files.")
  }

  ## --- Check for output folder ---
  if (!dir.exists(output_dir)) {
    stop("❌ Output directory 'output_fimpute' was not created. FImpute may have failed silently.")
  }

  ## --- Read results ---
  if (!exists("read.fimpute", mode = "function")) {
    stop("The function 'read.fimpute()' must be defined and available in the current environment.")
  }

  res <- read.fimpute(file = output_dir)

  if (!is.data.frame(res)) {
    stop("The function 'read.fimpute()' did not return a data.frame. Please verify its behavior.")
  }

  message("✔ Results successfully read from: ", output_dir)
  object@results <- res

  return(object)
})
