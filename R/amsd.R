#' Run Aggragate Mutation Spectrum Distance Method
#'
#' AMSD calculates significance of a difference between two groups of
#' mutation spectra by comparing the distance between the aggregate spectra
#' to a null distribution produced by random sampling, outputting a p-value.
#'
#' @param set1 A matrix or data.frame of mutation spectra for group 1
#' @param set2 A matrix or data.frame of mutation spectra for group 2
#' @param mean_or_sum Either "mean" (default) or "sum" to specify whether to
#'	aggregate spectrum is an average of the group's spectra ("mean": when
#'	comparing mutation spectra fractions that add up to 1) or sum ("sum"
#'	when comapring mutation counts)
#' @param n_sim Number of null resampling permutations to run (default = 1000)
#' @param seed Specify seed for reproducibility (default = NULL)
#' @return A list of three results: $cosine = observed cosine distance between
#'	the aggragate group spectra, $p  = p-value denoting the fraction of the
#'	random sampling permutations that are greater than or equal to the
#'	observed cosine distance, $sims = cosine distance for each of the random
#'	sampling permutations.
#' @export
amsd <- function(set1,
                 set2,
                 mean_or_sum = "mean",  # or "sum"
                 n_sim = 1000,
                 seed = NULL) {   
  
  # Validate inputs
  if (!is.matrix(set1) && !is.data.frame(set1)) stop("set1 must be a matrix or data.frame")
  if (!is.matrix(set2) && !is.data.frame(set2)) stop("set2 must be a matrix or data.frame")
  
  if (!mean_or_sum %in% c("mean", "sum")) {
    stop("Argument 'mean_or_sum' must be either 'mean' or 'sum'")
  }	

  if (mean_or_sum == "mean" & max(c(rowSums(set1), rowSums(set2)) > 1)) {
    stop("Spectra fractions should add up to 1 when running 'mean_or_sum' = 'mean'")
  }

  if (mean_or_sum == "sum" & max(c(rowSums(set1), rowSums(set2)) <= 1)) {
    stop("Run 'mean_or_sum' = 'mean' for fractional mutation spectra that add up to 1")
  }
  
  # Define aggregation function
  aggragate_spectra <- if (mean_or_sum == "mean") {
    function(data) colMeans(data, na.rm = TRUE)
  } else {
    function(data) colSums(data, na.rm = TRUE)
  }
  
  
  # Compute observed distance
  spectra1 <- aggragate_spectra(set1)
  spectra2 <- aggragate_spectra(set2)
  observed_distance <- cosine_dist(spectra1, spectra2)[[1]]
  
  # Prepare permutation dataset
  combined_set <- rbind(set1, set2)

  # Run permutations
  set.seed(seed)
  dist_rands <- numeric(n_sim)
  group_size <- nrow(set1)
  
  for (k in seq_len(n_sim)) {
    indices <- sample(seq_len(nrow(combined_set)), group_size)
    spectra_group1 <- aggragate_spectra(combined_set[indices, , drop = FALSE])
    spectra_group2 <- aggragate_spectra(combined_set[-indices, , drop = FALSE])
    dist_rands[k] <- cosine_dist(spectra_group1, spectra_group2)[[1]]
  }
  
  # Calculate p-value
  p_value <- max(c(mean(dist_rands >= observed_distance), 1 / n_sim))
  
  # Return output
  return(list(cosine = observed_distance, p = p_value, sims = dist_rands))
}