###############
# mpp_kinship #
###############

#' Compute linkage desequilibrium adjusted kinship
#'
#' Compute linkage disequilibrium adjusted kinship using the weights computed by
#' the program LDAK. (http://dougspeed.com/ldak/).
#'
#' The formula used to compute the kinship is the following:
#'
#' (1/sum(wj)) [sqrt(wj)(Xj - 2pj)][sqrt(wj)(Xj - 2pj)]'/[(2pj(1-pj))^alpha]
#'
#' The results return by the function are two times the results given by LDAK
#' using --calc-kins-direct.
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2.
#'
#' @param weights Two columns \code{data.frame} containing the marker
#' identifiers and the weights marker weights. For example the weights obtained
#' with the LDAK program (see \code{\link{LDAK_weights}}). By default, the
#' kinship matrix is computed unweighted or all weights are equal to 1, which
#' correspond to the Astle and Balding kinship matrix
#'
#' @param power \code{Numerical} value specifying the value of the
#' parameter for marker scores standardization. The column of the marker matrix
#' (X.j) are multiplied by var(X.j)^(power/2). It correspond to alpha in the
#' formula .Default = -1.
#'
#' @param mk.sel \code{Character vector} specifying a list of marker to use
#' for the kinship matrix computation. By default, the function use all markers
#' of the \code{gpData} object.
#'
#' @return Return:
#'
#' \item{K}{Kinship matrix computed with the LDAK weights.}
#'
#' @author Vincent Garin
#'
#' @references
#'
#' Astle, W., & Balding, D. J. (2009). Population structure and cryptic
#' relatedness in genetic association studies. Statistical Science, 451-471.
#'
#' Speed, D., Hemani, G., Johnson, M. R., & Balding, D. J. (2012).
#' Improved heritability estimation from genome-wide SNPs. The American Journal
#' of Human Genetics, 91(6), 1011-1021.
#'
#' @export
#'


# gp <- gp.imp
#
# # weights file
# # order the weights (same as map file)
#
# weights <- read.table("/home/vincent/Haplo_GRM/EUNAM/data/geno/LDAK_weights/EUNAM_LDAK_weights",
#                     h = TRUE)
#
# weights <- weights[, 1:2]
# power <- -1
# pz.ind <- substr(rownames(map), 1, 2)
# map_pz <- map[pz.ind == "PZ", ]
# mk.sel <- map_pz[map_pz[, 2] != 1, 1]

mpp_kinship <- function(gp, weights = NULL, power = -1, mk.sel = NULL){

  X <- gp$geno

  if(is.null(weights)){

    weights <- data.frame(colnames(X), rep(1, dim(X)[2]),
                          stringsAsFactors = FALSE)

  }

  # check that all markers have a weights

  if(sum(colnames(X) %in% weights[, 1]) != dim(X)[2]){

    stop(paste("Some marker of the genotype matrix of the gp object do not",
               "have a weight."))

  }

  # subset the X matrix nd the weight vector according to the mk selection

  if(!is.null(mk.sel)){

  X <- X[, colnames(X) %in% mk.sel]
  weights <- weights[weights[, 1] %in% mk.sel, ]

  }

  # order the weights information

  rownames(weights) <- weights[, 1]
  weights <- weights[colnames(X), ]

  wj <- weights[, 2]

  pj <- colMeans(X)/2
  varj <- 2 * pj * (1-pj)
  facj <- (varj^(-power/2)) * (1/sqrt(wj))

  X.w <- scale(x = X, center = 2*pj, scale = facj)
  den <- sum(wj)

  K <- tcrossprod(X.w) / den

  return(K)

}
