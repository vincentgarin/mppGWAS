##################
# plot_manhattan #
##################

#' Manhattan plot of Association Test results
#'
#' Plot the results of an association test using manhattan plot using the
#' \code{manhattan} function of the qqman package.
#'
#' @param res Object of class \code{G_res} obtained with function
#' \code{\link{GWAS_SNP}}, \code{\link{GWAS_haplo}}, \code{\link{GWAS_GSM}}.
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
