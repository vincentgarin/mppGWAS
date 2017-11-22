##############
# plot.G_res #
##############

#' Plot Association test results
#'
#' S3 \code{plot} method for object of class \code{G_res}. The function
#' return either a manhattan plot (\code{type = "man"}) or a QQ plot
#' (\code{type = "qq"})
#'
#' The function call functions from the package qqman (Turner, 2014).
#'
#' @param x Object of class \code{G_res} obtained with function
#' \code{\link{GWAS_SNP}}, \code{\link{GWAS_haplo}}, \code{\link{GWAS_GSM}}.
#'
#' @param pl.type \code{Character} expression. Either "man" for a manhattan plot
#' or "qq" for a QQ plot. Default = "man".
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @author Vincent Garin
#'
#' @references
#'
#' Turner, S.D. qqman: an R package for visualizing GWAS results using Q-Q and
#' manhattan plots. biorXiv DOI: 10.1101/005165 (2014).
#'
#' # @export: not functionning for the moment
#'

plot.G_res <- function(x, pl.type = "man", ...){

  stopifnot(inherits(x, "G_res"))

  if(pl.type == "man"){

    colnames(x) <- c("SNP","CHR", "BP", "P")
    qqman::manhattan(x, logp = FALSE, ...)


  } else if (pl.type == "qq"){

    p.val <- 10^(-x[, 4])
    qqman::qq(p.val, ...)

  }

}
