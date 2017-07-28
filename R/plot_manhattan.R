##################
# plot_manhattan #
##################

#' Manhattan plot of Association Test results
#'
#' Plot the results of an association test using manhattan plot using the
#' \code{manhattan} function of the qqman package.
#'
#' @param res Four column \code{data.frame} with: marker identifier, chromosome,
#' position (cM or bp) and -log10(p-value). For example result of function
#' \code{AssTest_SNP}.
#'
#' @param ... Arguments passed on to other plot/points functions
#'
#' @author Vincent Garin
#'
#' @export
#'


plot_manhattan <- function(res, ...){

  # rownames(res) <- res[, 1]
  colnames(res) <- c("SNP","CHR", "BP", "P")
  qqman::manhattan(res, logp = FALSE, ...)

}
