#' Build FImputeRunner object
#'
#' A convenience function to construct a `FImputeRunner` object from a `SNPDataLong` object.
#'
#' @param object An object of class `SNPDataLong`, from which `geno` and `map` slots will be extracted.
#' @param path A character string indicating the directory to save FImpute files.
#' @param exec_path Path to the FImpute executable (default = "FImpute3").
#' @param name Name for the dataset (used internally, default = "data").
#'
#' @return An object of class `FImputeRunner`.
#' @export
FImputeRunner <- function(object, path, exec_path = "FImpute3", name = "data") {

  geno <- object@geno
  map <- object@map

  export <- new("FImputeExport", geno = geno, map = map, path = path, name = name)
  par_file <- file.path(path, "fimpute.par")
  runner <- new("FImputeRunner", export = export, par_file = par_file, exec_path = exec_path, results = data.frame())
  return(runner)
}
