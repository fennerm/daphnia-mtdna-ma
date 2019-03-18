#' Remove NA values from a vector
#'
#' @export
remove_na <- function(vec) {
  vec[!is.na(vec)]
}

#' Count occurences of a value in a possibly nested list
#'
#' @export
count_occurences <- function(a_list, value) {
  length(which(unlist(a_list) == "transversion"))
}
