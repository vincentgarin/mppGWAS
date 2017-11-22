###########
# plot_qq #
###########

#' QQ plot of Association Test results
#'
#' Plot the results of an association test using manhattan plot using the
#' \code{manhattan} function of the qqman package.
#'
#' @param res Object of class \code{G_res} obtained with function
#' \code{\link{GWAS_SNP}}, \code{\link{GWAS_haplo}}, \code{\link{GWAS_GSM}}.
#'
#' @author Vincent Garin
#'
#' @export
#'


plot_qq <- function(res){

  p.val <- 10^(-res[, 4])
  qqman::qq(p.val)

}
