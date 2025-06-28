#' Formatted header message
#'
#' Prints a formatted message with a border for section titles in the console.
#'
#' @param title Character string to be printed inside the header box.
#'
#' @return No return value. Used for side effects (message).
#'
#' @examples
#' qc_header("Quality Control on Samples")
#'
#' @export
qc_header <- function(title) {
  bar <- strrep("═", nchar(title) + 8)
  message("\n╔", bar, "╗")
  message("║    ", title, "    ║")
  message("╚", bar, "╝\n")
}
