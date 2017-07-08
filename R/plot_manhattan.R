##################
# plot_manhattan #
##################

#' Manhattan plot of Association Test results
#'
#' Plot the results of an association test using manhattan plot using the
#' \code{manhattan} function of the qqman package.
#'
#' @param res result of function \code{SNP_AssTest}.
#'
#' @author Vincent Garin
#'
#' @export
#'


plot_manhattan <- function(res){

  # rownames(res) <- res[, 1]
  colnames(res) <- c("SNP","CHR", "BP", "P")
  qqman::manhattan(res, logp = FALSE)

}
