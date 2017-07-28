###############
# AssTest_SNP #
###############

#' Single SNP Association Test
#'
#' Perform a single SNP genome wide association test using a kinship matrix
#' to correct for the genetic background.
#'
#' The function is a wrapper for the \code{mmer} function from the \code{sommer}
#' package (Covarrubias-Pazaran, 2016). The model is fitted using the EMMA
#' algorithm proposed by Kang et al. (2008).
#'
#' By default, the kinship matrix is computed using the method of Astle and
#' Balding (2009). It is possible to compute a linkage disequilibrium adjusted
#' kinship (LDAK) kinship matrix using the method of Speed et al. (2012) by
#' introducing the weights computed with the function \code{\link{LDAK_weights}}.
#' The model can be fitted using the kinship containing all markers or removing
#' the markers of the scanned chromosome (\code{K_i = TRUE}).
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2 and family.
#'
#' @param trait \code{Numerical} vector of phenotypic trait values.
#'
#' @param map Three columns \code{data.frame} with marker identifier,
#' chromosome, and marker position (cM or bp).
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
#' of the \code{gp}.
#'
#' @param K_i \code{Logical} specifying if the kinship correction should be done
#' by removing the markers of the scanned chromosome. Default = TRUE.
#'
#' @param n.cores Default = 1.
#'
#' @return Return:
#'
#' \item{results}{Four column data.frame with marker identifier, chromosome,
#' position and -log10(p-value).}
#'
#' @author Vincent Garin
#'
#' @references
#'
#' Astle, W., & Balding, D. J. (2009). Population structure and cryptic
#' relatedness in genetic association studies. Statistical Science, 451-471.
#'
#' Covarrubias-Pazaran G. 2016. Genome assisted prediction of quantitative traits
#' using the R package sommer. PLoSONE 11(6):1-15.
#'
#' Speed, D., Hemani, G., Johnson, M. R., & Balding, D. J. (2012).
#' Improved heritability estimation from genome-wide SNPs. The American Journal
#' of Human Genetics, 91(6), 1011-1021.
#'
#' @export
#'


AssTest_SNP <- function(gp, trait, map, weights = NULL, power = -1,
                        mk.sel = NULL, K_i = TRUE, n.cores = 1){

  # check that the list of markers in the map is the same as the list of the
  # genotype matrix

  if(!identical(colnames(gp$geno), map[, 1])){

  stop(paste("The list of marker in the genotype matrix of the gpData object",
             "and in the map are not stricly equivalent",
             "(same content, same order)."))

  }

  # check that there is a weight for each marker if weights are given

  if(!is.null(weights)){

    if(!(sum(colnames(gp$geno) %in% weights[, 1]) == dim(gp$geno)[2])){

     stop(paste("The weight argument does not contain a value for each marker",
                "in the gpData object (gp)."))

    }
  }

  # check that the selection of marker is present in the map

  if( sum(!(mk.sel %in% map[, 1])) !=0 ){

    prob.mk <- paste(mk.sel[!(mk.sel %in% map[, 1])], collapse = ", ")

    stop(paste("the following markers:", prob.mk, "are present in mk.sel but",
         "not in the map."))

  }

  ############ end checks

  X <- IncMat_cross(gp$covar$family)

  if(!K_i){ # Use the whole genome for K.

    K <- mpp_kinship(gp = gp, weights = weights, power = power, mk.sel = mk.sel)

    Z1 <- diag(length(trait))
    ETA <- list( list(Z=Z1, K=K))
    ans <- mmer(Y = trait, X = X, Z = ETA, W = gp$geno, method = "EMMA",
                n.cores = n.cores)

    results <- data.frame(map, ans$W.scores$additive[, 1])
    colnames(results) <- c("mk.id", "Chrom", "Position", "p.val")

  } else { # Remove the kth chromosome for the computation of K.

    # K all markers - ith chromosome

    p.val <- c()
    n.chr <- length(unique(map[, 2]))
    chr.id <- unique(map[, 2])

    if(!is.null(mk.sel)){

      mk.sel_temp <- map[map[, 1] %in% mk.sel, ]

    } else { mk.sel_temp <- map }

    for(i in 1:n.chr){

      mk_chr_i <- map[map[, 2] == chr.id[i], 1]

      # adapt the list of selected markers

    mk.sel_i <- mk.sel_temp[mk.sel_temp[, 2] != chr.id[i], 1]

     # obtain a reduced genotype matrix

      geno_i <- gp$geno[, colnames(gp$geno) %in% mk_chr_i]

      K <- mpp_kinship(gp = gp, weights = weights, power = power,
                       mk.sel = mk.sel_i)

      Z1 <- diag(length(trait))
      ETA <- list( list(Z=Z1, K=K))
      ans <- mmer(Y = trait, X = X, Z = ETA, W = geno_i, method = "EMMA",
                  n.cores = n.cores)
      p.val <- c(p.val, ans$W.scores$additive[, 1])

    }

    results <- data.frame(map, p.val)
    colnames(results) <- c("mk.id", "Chrom", "Position", "p.val")

  }

  return(results)

}
