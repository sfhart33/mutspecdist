#' Plot Histogram of AMSD Output
#'
#' Generates a ggplot histogram of the AMSD permutations, with the observed
#' cosine distance marked by a line and the p-value in the title.
#'
#' @param output The output from amsd (list of three: $cosine, $p, $sims)
#' @return A ggplot histogram
#' @export
plot_amsd_histogram <- function(output){
  plot1 <- ggplot2::ggplot(data.frame(sims = output$sims), ggplot2::aes(sims)) +
    ggplot2::geom_histogram()+
    ggplot2::geom_vline(xintercept = output$cosine, linetype = "dashed")+
    ggplot2::xlab("Cosine distance")+
    ggplot2::ggtitle(paste0("p = ", output$p)) +
    ggplot2::theme_classic()
  return(plot1)
}