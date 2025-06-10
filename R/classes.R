# Assegura que a classe SnpMatrix esteja registrada
if (requireNamespace("snpStats", quietly = TRUE)) {
  getClass("SnpMatrix", where = asNamespace("snpStats"))
  setClassUnion("SnpMatrixOrNULL", c("SnpMatrix", "NULL"))
} else {
  warning("Pacote 'snpStats' não encontrado. Usando fallback 'NULL' para SnpMatrixOrNULL.")
  setClassUnion("SnpMatrixOrNULL", "NULL")
}

# Classes S4
setClass("SNPDataLong",
         slots = c(
           geno = "SnpMatrix",   # sem o prefixo snpStats::
           map = "data.frame",
           path = "character"
         ))

setClass("SNPFileConfig",
         slots = c(
           path = "character",
           fields = "list",
           codes = "character",
           threshold = "numeric",
           sep = "character",
           skip = "numeric"
         ))

setClass("SNPImportList",
         slots = c(configs = "list"))

setClass("FImputeExport",
         slots = c(
           geno = "SnpMatrixOrNULL",
           map = "data.frame",
           path = "character",
           name = "character"
         ))

setClass("FImputeRunner",
         slots = c(
           export = "FImputeExport",
           par_file = "character",
           exec_path = "character",
           results = "data.frame"
         ))
