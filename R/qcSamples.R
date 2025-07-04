#' Quality control on samples
#'
#' Applies quality control (QC) procedures to samples in a `SNPDataLong` object,
#' based on heterozygosity and call rate thresholds.
#'
#' @param x An object of class `SNPDataLong`.
#' @param heterozygosity A numeric threshold or range for heterozygosity. Samples outside this threshold are removed.
#' @param smp_cr Minimum acceptable sample call rate (between 0 and 1). Samples below this value are removed.
#' @param action Character string indicating the action to perform. One of:
#'   - `"report"`: only returns a list of samples to remove and those kept;
#'   - `"filter"`: returns a filtered object without reporting;
#'   - `"both"`: performs filtering and returns the filtered object.
#'
#' @return Depending on the `action` argument:
#'   - `"report"`: returns a list with removed and kept samples;
#'   - `"filter"`: returns a new `SNPDataLong` object with filtered genotypes;
#'   - `"both"`: returns a list with:
#'       - `filtered`: the filtered `SNPDataLong` object;
#'       - `report`: a list of removed and kept samples.
#'
#' @importFrom fQC check.sample.heterozygosity
#' @export
setGeneric("qcSamples", function(x, ...) standardGeneric("qcSamples"))

#' @rdname qcSamples
setMethod("qcSamples", "SNPDataLong", function(x, 
                                               heterozygosity = NULL,
                                               smp_cr = NULL,
                                               action = c("report", "filter", "both")) {
  action <- match.arg(action)

  geno <- x@geno
  map <- x@map
  xref_path <- x@xref_path
  path <- x@path

  # Verifica e trata amostras duplicadas
  dups_logical <- duplicated(rownames(geno))
  if (any(dups_logical)) {
    message("\n⚠️  Duplicated sample identifiers detected in the SnpMatrix object.")
    message("   Only the first occurrence of each duplicated sample will be retained.")
    geno <- geno[!dups_logical, , drop = FALSE]
    if (!is.null(xref_path) && length(xref_path) == length(dups_logical)) {
      xref_path <- xref_path[!dups_logical]
    }
  }

  qc_header("🧬 Quality Control on Samples")
  message("Initial number of samples: ", nrow(geno))
  message("Applying quality control filters:")

  keep_samples <- rownames(geno)

  # Estatísticas por indivíduo
  sample.qc <- row.summary(geno)

  removed_hetero <- removed_cr <- character()

  if (!is.null(heterozygosity)) {
    removed_hetero <- fQC::check.sample.heterozygosity(sample.qc, heterozygosity)
    keep_samples <- setdiff(keep_samples, removed_hetero)
    message(sprintf("  • Heterozygosity filter: %d sample(s) removed, %d remaining.",
                    length(removed_hetero), length(keep_samples)))
  }

  if (!is.null(smp_cr)) {
    removed_cr <- check.sample.call.rate(sample.qc, smp_cr)
    keep_samples <- setdiff(keep_samples, removed_cr)
    message(sprintf("  • Call rate filter: %d sample(s) removed, %d remaining.",
                    length(removed_cr), length(keep_samples)))
  }

  if (length(keep_samples) == 0) {
    stop("No samples passed the quality control filters.")
  }

  if (action == "report") {
    return(list(
      removed_heterozygosity = removed_hetero,
      removed_call_rate = removed_cr,
      kept_samples = keep_samples
    ))
  }

  # Aplicar filtro
  filtered_geno <- geno[rownames(geno) %in% keep_samples, , drop = FALSE]
  if (!is.null(xref_path) && length(xref_path) == nrow(geno)) {
    filtered_xref <- xref_path[rownames(geno) %in% keep_samples]
  } else {
    filtered_xref <- xref_path
  }

  filtered_obj <- new("SNPDataLong",
                      geno = filtered_geno,
                      map = map,
                      path = path,
                      xref_path = filtered_xref)

  if (action == "filter") {
    return(filtered_obj)
  }

  if (action == "both") {
    return(list(
      filtered = filtered_obj,
      report = list(
        removed_heterozygosity = removed_hetero,
        removed_call_rate = removed_cr,
        kept_samples = keep_samples
      )
    ))
  }
})

#' Check Sample Call Rate
#'
#' Identifies samples with call rate below a given threshold.
#'
#' @param sample.summary A data frame with a "Call.rate" column for each sample.
#' @param min.call.rate Minimum acceptable call rate (between 0 and 1).
#'
#' @return A character vector with the names of samples to remove.
#' @export
check.sample.call.rate <- function(sample.summary, min.call.rate) {
  result <- sample.summary[, "Call.rate"] < min.call.rate
  result[is.na(result)] <- FALSE
  smps <- NULL
  if (sum(result) > 0) {
    smps <- rownames(sample.summary[result, ])
  }
  return(smps)
}
