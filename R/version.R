#' Bundled SparseDiffEngine version
#'
#' Returns the version string of the SparseDiffEngine C library bundled with
#' this package.
#'
#' @return A character scalar, e.g. \code{"0.3.0"}.
#' @export
#' @examples
#' engine_version()
engine_version <- function() {
  sd_engine_version()
}
