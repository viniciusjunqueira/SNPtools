# Ensure that the SnpMatrix class is registered
if (requireNamespace("snpStats", quietly = TRUE)) {
  getClass("SnpMatrix", where = asNamespace("snpStats"))
  setClassUnion("SnpMatrixOrNULL", c("SnpMatrix", "NULL"))
} else {
  warning("Package 'snpStats' not found. Using fallback 'NULL' for SnpMatrixOrNULL.")
  setClassUnion("SnpMatrixOrNULL", "NULL")
}

# Define S4 classes for SNP data and related configurations

# Class for storing SNP data in long format
setClass("SNPDataLong",
         slots = c(
           geno = "SnpMatrix",    # Genotype matrix (without snpStats:: prefix)
           map = "data.frame",    # Marker map information
           path = "character"     # File path or identifier
         ))

# Class for configuration of SNP file import
setClass("SNPFileConfig",
         slots = c(
           path = "character",    # File path
           fields = "list",       # Field mapping or column specifications
           codes = "character",   # Codes for genotype or allele representation
           threshold = "numeric", # Quality or filtering threshold
           sep = "character",     # Field separator
           skip = "numeric"       # Number of lines to skip at the start
         ))

# Class for storing a list of SNP file configurations
setClass("SNPImportList",
         slots = c(
           configs = "list"       # List of SNPFileConfig objects
         ))

# Class for preparing data export for FImpute
setClass("FImputeExport",
         slots = c(
           geno = "SnpMatrixOrNULL", # Genotype matrix or NULL
           map = "data.frame",       # Marker map
           path = "character",       # Output file path
           name = "character"        # File or project name
         ))

# Class for managing FImpute execution and results
setClass("FImputeRunner",
         slots = c(
           export = "FImputeExport", # FImpute export object
           par_file = "character",   # Path to FImpute parameter file
           exec_path = "character",  # Path to FImpute executable
           results = "data.frame"    # Imputation results or summary
         ))
