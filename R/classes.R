# Ensure that the SnpMatrix class is registered
if (requireNamespace("snpStats", quietly = TRUE)) {
  getClass("SnpMatrix", where = asNamespace("snpStats"))
  setClassUnion("SnpMatrixOrNULL", c("SnpMatrix", "NULL"))
} else {
  warning("Package 'snpStats' not found. Using fallback 'NULL' for SnpMatrixOrNULL.")
  setClassUnion("SnpMatrixOrNULL", "NULL")
}

# Define union for map slot before using it
setClassUnion("MapDataFrameOrList", c("data.frame", "list"))


#' SNPDataLong Class
#'
#' A class to store SNP genotype data in long format, including genotype matrix, marker map, and file paths.
#'
#' @slot geno A SnpMatrix containing genotype data.
#' @slot map A data.frame or list containing marker information.
#' @slot path A character string with the file path or identifier.
#' @slot xref_path A character string with per-individual paths or identifiers.
#'
#' @export
setClass("SNPDataLong",
         slots = c(
           geno = "SnpMatrix",
           map = "MapDataFrameOrList",
           path = "character",
           xref_path = "character"
         ))


#' SNPFileConfig Class
#'
#' A class for configuring SNP file import options.
#'
#' @slot path Path to the SNP file.
#' @slot fields A list specifying column mappings or field configurations.
#' @slot codes Character vector for genotype or allele codes.
#' @slot threshold Numeric value for filtering or quality control.
#' @slot sep Character specifying the field separator.
#' @slot skip Number of lines to skip at the top of the file.
#'
#' @export
setClass("SNPFileConfig",
         slots = c(
           path = "character",
           fields = "list",
           codes = "character",
           threshold = "numeric",
           sep = "character",
           skip = "numeric"
         ))


#' SNPImportList Class
#'
#' A class for managing a list of SNP file import configurations.
#'
#' @slot configs A list of SNPFileConfig objects.
#'
#' @export
setClass("SNPImportList",
         slots = c(
           configs = "list"
         ))


#' FImputeExport Class
#'
#' A class to handle export preparation for FImpute.
#'
#' @slot geno A SnpMatrix or NULL containing genotype data.
#' @slot map A data.frame containing marker information.
#' @slot path Output file path.
#' @slot name Project or file name.
#'
#' @export
setClass("FImputeExport",
         slots = c(
           geno = "SnpMatrixOrNULL",
           map = "data.frame",
           path = "character",
           name = "character"
         ))


#' FImputeRunner Class
#'
#' A class to manage FImpute execution and results.
#'
#' @slot export An FImputeExport object.
#' @slot par_file Path to parameter file.
#' @slot exec_path Path to FImpute executable.
#' @slot results A data.frame containing results or summary information.
#'
#' @export
setClass("FImputeRunner",
         slots = c(
           export = "FImputeExport",
           par_file = "character",
           exec_path = "character",
           results = "data.frame"
         ))
