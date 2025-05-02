#' Simulates Mutation Spectra for AMSD Input
#'
#' Generates simulated mutation spectra 
#'
#' @param n_samples Number of samples to simulate 
#' @param n_mutations Number of total mutations per sample
#' @param sig_probs Specify signatures to sample from and the fraction of
#'	signatures to sample from each
#' @param signatures Matrix of signature spectra, with each row a signature 
#'	specified by rowname. Must incude names from sig_probs
#' @param additional_sig Additional signature to add as an exposure
#' @param n_extra Number of mutations to sample from additional signature
#' @param seed Specify seed for reproducibility (default = NULL)
#' @return List n_samples long of simluated spectra mutation counts
#' @export
simulate_spectra <- function(n_samples = 1,
                             n_mutations = 100,
                             sig_probs = c(SBS1 = 0.3, SBS5 = 0.6, SBS18 = 0.1),
                             signatures,
                             additional_sig = NULL,
                             n_extra = 0,
                             seed = NULL) {
  set.seed(seed)
  
# loop through samples and save 3mer counts
  spectras <- list()
  for (i in 1:n_samples) {
    
  # Loop through each signature and sample 3mer counts
    spectra <- c()
    for(s in 1:length(sig_probs)){
      
    # Which signature
      sig <- names(sig_probs)[s]
      sig_spectra <- signatures[sig,]
      
    # Sample mutations from given signatures
      mutations <- sample(names(sig_spectra),
                          size = (n_mutations) * sig_probs[sig],
                          replace = TRUE,
                          prob = sig_spectra)
      
      spectra[[s]] <- table(factor(mutations, levels = names(sig_spectra)))
    }
    
  # Add additional mutations if applicable
    if (!is.null(additional_sig) && n_extra[i] > 0) {
      extra_mutations <- sample(names(signatures[additional_sig,]),
                                size = n_extra[i],
                                replace = TRUE,
                                prob = signatures[additional_sig,])
      spectra[[length(sig_probs)+1]] <- table(factor(extra_mutations, levels = names(sig_spectra)))
    } 
    
    spectras[[i]] <- colSums(do.call(rbind, spectra))
  }
  return(spectras)
}