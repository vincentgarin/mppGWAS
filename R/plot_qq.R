###########
# plot_qq #
###########

#' QQ plot of Association Test results
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


plot_qq <- function(res){

  p.val <- 10^(-res[, 4])
  qq(p.val)

}
