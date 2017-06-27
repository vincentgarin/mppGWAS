#################
# LDAK_kinship2 #
#################

#' Compute linkage desequilibrium adjusted kinship
#'
#' Compute linkage disequilibrium adjusted kinship using the weights computed by
#' the program LDAK. (http://dougspeed.com/ldak/).
#'
#' The results return by the function are two times the results given by LDAK
#' using --calc-kins-direct.
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2.
#'
#' @param weights Four columns \code{data.frame} containing LDAK weights
#' information. The columns must be: 1) marker identifier, chromosome,
#' position (in cM or bp), LDAK weights (see \code{\link{LDAK_weights}}).
#'
#' @param power \code{Numerical} value specifying the value of the
#' parameter for marker scores standardization. The column of the marker matrix
#' (X.j) are multiplied by var(X.j)^(power/2) .Default = -1.
#'
#' @param K_i \code{Numerical} value specifying a unique chromosome number that
#' should be removed from the kinship computation. By default \code{K_i = NULL},
#' which means that the kinship is computed using all markers.
#'
#' @return Return:
#'
#' \item{K}{Kinship matrix computed with the LDAK weights.}
#'
#' @author Vincent Garin
#'
#' @references
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
# wghts2 <- wghts[map[, 1], ]
# weights <- data.frame(map[, c(1, 2, 4)], wghts2$Weight, stringsAsFactors = FALSE)
# K_i <- 1
# power <- -1

LDAK_kinship2 <- function(gp, weights, power = -1, K_i = NULL){

  X <- gp$geno

  # check that all markers have a weights

  if(sum(colnames(X) %in% weights[, 1]) != dim(X)[2]){

    stop(paste("Some marker of the genotype matrix of the gp object do not",
               "have a weight."))

  }

  # order (and subset) the weights information

  rownames(weights) <- weights[, 1]
  weights <- weights[colnames(X), ]

  if(!is.null(K_i)){

  weights[weights[, 2] == K_i, 4] <- 0

  }

  wj <- weights[, 4]

  pj <- colMeans(X)/2
  varj <- 2 * pj * (1-pj)
  facj <- (varj^(-power/2)) * (1/sqrt(wj))

  X.w <- scale(x = X, center = 2*pj, scale = facj)
  den.wj <- sum(wj)

  K <- (X.w %*% t(X.w))/den.wj

  return(K)

}
