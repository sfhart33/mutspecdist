#' Calculate Cosine Distance
#'
#' This function calculates the cosine distance between two vectors
#'
#' @param A A vector (e.g. mutation spectrum: 96 3mer contexts)
#' @param B A second vector
#' @return Cosine distance between the vectors (0-1, 0 = identical)
#' @export
cosine_dist <- function(A, B){
  cosine_similarity <- sum(A * B) / (sqrt(sum(A^2)) * sqrt(sum(B^2)))
  cosine_distance <- 1 - cosine_similarity
  return(cosine_distance)
}