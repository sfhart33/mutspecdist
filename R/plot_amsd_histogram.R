#' Plot Histogram of AMSD Output
#'
#' Generates a ggplot histogram of the AMSD permutations, with the observed
#' cosine distance marked by a line and the p-value in the title.
#'
#' @param output The output from amsd (list of three: $cosine, $p, $sims)
#' @return A ggplot histogram
#' @export
plot_amsd_histogram <- function(output){
  plot1 <- ggplot2::ggplot(data.frame(sims = output$sims), aes(sims)) +
    geom_histogram()+
    geom_vline(xintercept = output$cosine, linetype = "dashed")+
    xlab("Cosine distance")+
    ggtitle(paste0("p = ", output$p)) +
    theme_classic()
  return(plot1)
}