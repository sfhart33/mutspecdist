% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_spectra_amsd.R
\name{simulate_spectra_amsd}
\alias{simulate_spectra_amsd}
\title{Simulates Spectra for AMSD Input and Runs AMSD}
\usage{
simulate_spectra_amsd(
  n_samples,
  n_mutations,
  sig_probs,
  signatures,
  additional_sig,
  n_extra,
  n_sim,
  seed
)
}
\arguments{
\item{n_samples}{Number of samples to simulate}

\item{n_mutations}{Number of total mutations per sample}

\item{sig_probs}{Specify signatures to sample from and the fraction of
signatures to sample from each}

\item{signatures}{Matrix of signature spectra, with each row a signature
specified by rowname. Must incude same names from sig_probs}

\item{additional_sig}{Additional signature to add as an exposure}

\item{n_extra}{Number of mutations to sample from additional signature}

\item{n_sim}{Number of null resampling permutations to run}

\item{seed}{Specify seed for reproducibility}
}
\value{
A list of three results: $cosine = observed cosine distance between
the aggragate group spectra, $p  = p-value denoting the fraction of the
random sampling permutations that are greater than or equal to the
observed cosine distance, $sims = cosine distance for each of the random
sampling permutations.
}
\description{
Combines the amsd and simluate_spectra functions. See ?amsd and ?simluate_spectra for details
}
