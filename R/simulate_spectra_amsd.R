#' Simulates Spectra for AMSD Input and Runs AMSD
#'
#' Combines the amsd and simluate_spectra functions. See ?amsd and ?simluate_spectra for details
#'
#' @param n_samples Number of samples to simulate 
#' @param n_mutations Number of total mutations per sample
#' @param sig_probs Specify signatures to sample from and the fraction of
#'	signatures to sample from each
#' @param signatures Matrix of signature spectra, with each row a signature 
#'	specified by rowname. Must incude names from sig_probs
#' @param additional_sig Additional signature to add as an exposure
#' @param n_extra Number of mutations to sample from additional signature
#' @param n_sim Number of null resampling permutations to run
#' @param seed Specify seed for reproducibility
#' @return A list of three results: $cosine = observed cosine distance between
#'	the aggragate group spectra, $p  = p-value denoting the fraction of the
#'	random sampling permutations that are greater than or equal to the
#'	observed cosine distance, $sims = cosine distance for each of the random
#'	sampling permutations.

#' @export
simulate_spectra_amsd <- function(n_samples, n_mutations, sig_probs, signatures, additional_sig, n_extra, n_sim, seed) {
  
  # Run on a base set, then with exposures
  no_exposure_test <- simulate_spectra(n_samples = n_samples,
                                       n_mutations = n_mutations,
                                       sig_probs = sig_probs,
                                       signatures = signatures)
  no_exposure_test <- as.data.frame(do.call(rbind, no_exposure_test))
  no_exposure_test <- no_exposure_test/rowSums(no_exposure_test) # convert to mutaion fractions
  with_exposure_test <- simulate_spectra(n_samples = n_samples,
                                         n_mutations = n_mutations,
                                         sig_probs = sig_probs,
                                         signatures = signatures,
                                         additional_sig = additional_sig,
                                         n_extra = n_extra)
  with_exposure_test <- as.data.frame(do.call(rbind, with_exposure_test))
  with_exposure_test <- with_exposure_test/rowSums(with_exposure_test) # convert to mutaion fractions
  amsd(no_exposure_test, with_exposure_test, n_sim = n_sim, seed = seed) %>%
    return()
}
