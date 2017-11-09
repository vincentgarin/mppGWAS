############
# GWAS_SNP #
############

#' Single SNP GWAS
#'
#' Perform a single SNP genome wide association test using a kinship matrix
#' to correct for the genetic background.
#'
#' The function is a wrapper for the \code{GWAS} function from the \code{sommer}
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
#' @param gp \code{gpData} object with elements geno coded 0 1 2, map with
#' marker position in cM, phenotype and family indicator.
#'
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the gp object should be used. Default = 1.
#'
#' @param weights object of class \code{LD_wgh} obtained by the function
#' LDAK_weights() representing a data.frame with two columns: the marker
#' identifier and the LD adjusted weights. These weight will be used to compute
#' a LDAK as defined by Speed et al. (2012) for the genetic background
#' adjustement. Default = NULL.
#'
#' @param power \code{Numerical} value specifying the value of the
#' parameter for marker scores standardization. The column of the marker matrix
#' (X.j) are multiplied by var(X.j)^(power/2). It correspond to alpha in the
#' formula. Default = -1.
#'
#' @param mk.sel \code{Character vector} specifying a list of marker to use
#' for the kinship matrix computation. By default, the function use all markers
#' of the \code{gp}.
#'
#' @param K_i \code{Logical} specifying if the kinship correction should be done
#' by removing the markers of the scanned chromosome. Default = TRUE.
#'
#' @param n.cores \code{Numeric} value indicating the number of core to be used
#' if the user want to run the K-i (\code{K_i = TRUE}) scan in parallel.
#' Default = NULL.
#'
#' @param verbose \code{Logical} indicating if function outputs should be printed.
#' Default = FALSE.
#'
#' @return Return:
#'
#' \item{G.res}{Object of class \code{G_res} representing a data.frame with four
#' columns: marker identifier, chromosome, position in cM and -log10(p-value).}
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
#' Kang, H. M., Zaitlen, N. A., Wade, C. M., Kirby, A., Heckerman, D., Daly,
#' M. J., & Eskin, E. (2008). Efficient control of population structure in model
#' organism association mapping. Genetics, 178(3), 1709-1723.
#'
#' Speed, D., Hemani, G., Johnson, M. R., & Balding, D. J. (2012).
#' Improved heritability estimation from genome-wide SNPs. The American Journal
#' of Human Genetics, 91(6), 1011-1021.
#'
#' @export
#'

# arguments

# source('~/Haplo_GRM/mppGWAS/R/test_class.R')
# source('~/Haplo_GRM/mppGWAS/R/test_gpData_content.R')
#
# data("EUNAM_gp")
# data("EUNAM_LD_weights")
#
# gp <- EUNAM_gp
# trait <- 1
# weights <- EUNAM_LD_weights
# power = -1
# mk.sel = NULL
# K_i = TRUE
# verbose = FALSE


GWAS_SNP <- function(gp, trait = 1, weights = NULL, power = -1, mk.sel = NULL,
                     K_i = TRUE, n.cores = NULL, verbose = FALSE){

  # test gpData and his content

  test_gpData_content(gp)

  if(!is.null(weights)){

    if(!is_LD_wgh(weights)){

      stop("The LD weights are not of class LD_wgh obtained with LDAK_weights().")

    }

  }

  # check that there is a weight for each marker if weights are given

  if(!is.null(weights)){

    if(!(sum(colnames(gp$geno) %in% weights[, 1]) == dim(gp$geno)[2])){

     stop(paste("The weight argument does not contain a value for each marker",
                "in the gpData object (gp)."))

    }
  }

  # reconstruct the map

  map <- data.frame(rownames(gp$map), gp$map, stringsAsFactors = FALSE)
  colnames(map) <- c("mk.id", "chr", "cM")

  # check that the selection of marker is present in the map

  if( sum(!(mk.sel %in% map[, 1])) !=0 ){

    prob.mk <- paste(mk.sel[!(mk.sel %in% map[, 1])], collapse = ", ")

    stop(paste("the following markers:", prob.mk, "are present in mk.sel but",
         "not in the map."))

  }

  # check the value of trait

  if(!(is.character(trait) || is.numeric(trait))){

    stop("The trait object must be either numeric or character")

  }

  if(is.numeric(trait)){

  nb.trait <- dim(gp$pheno)[3]

  if(!((0 < trait) && (trait <=nb.trait))){

    stop(paste("The trait indicator should be between 1 and", nb.trait,
               "the total number of traits in the gp object"))

  }

  }

  if(is.character(trait)){

    trait.names <- attr(gp$pheno, "dimnames")[[2]]

    if (!(trait %in% trait.names)){

    stop(paste("trait must be one of:", trait.names))

    }

  }

  ############ end checks

  # select the trait

  if(is.numeric(trait)){

    pheno <- gp$pheno[, 1, trait]

  } else {

    pheno <- gp$pheno[, 1, which(trait %in% trait.names)]

  }

  X <- IncMat_cross(gp$covar$family)

  if(!K_i){ # Use the whole genome for K.

    K <- mpp_kinship(gp = gp, weights = weights, power = power, mk.sel = mk.sel)

    Z1 <- diag(length(pheno))
    ETA <- list( list(Z=Z1, K=K))
    ans <- sommer::GWAS(Y = pheno, X = X, Z = ETA, W = gp$geno, method = "EMMA",
                silent = !verbose, gwas.plots = FALSE)

    G.res <- data.frame(map, ans$W.scores$score[1, ], stringsAsFactors = FALSE)
    colnames(G.res) <- c("mk.id", "Chrom", "Position", "p.val")
    class(G.res) <- c("data.frame", "G_res")

  } else { # Remove the kth chromosome for the computation of K.

    # K all markers - ith chromosome

    p.val <- c()
    n.chr <- length(unique(map[, 2]))
    chr.id <- unique(map[, 2])

    if(!is.null(mk.sel)){

      mk.sel_temp <- map[map[, 1] %in% mk.sel, ]

    } else { mk.sel_temp <- map }

    # Execution of the K-i scan in parallel

    if(!is.null(n.cores)){

      cl <- makeCluster(n.cores)
      registerDoParallel(cl)

      res <- foreach(i=1:n.chr) %dopar% {

        GWAS_SNP_i(i = i, chr.id = chr.id, gp = gp, X = X, pheno = pheno, map = map,
                   mk.sel_temp = mk.sel_temp, weights = weights, power = power)

      }

      stopCluster(cl)

      p.val <- unlist(res)

    } else { # or execution of the chromosome scan in a regular for loop

      for(i in 1:n.chr){

        mk_chr_i <- map[map[, 2] == chr.id[i], 1]

        # adapt the list of selected markers

        mk.sel_i <- mk.sel_temp[mk.sel_temp[, 2] != chr.id[i], 1]

        # obtain a reduced genotype matrix

        geno_i <- gp$geno[, colnames(gp$geno) %in% mk_chr_i]

        K <- mpp_kinship(gp = gp, weights = weights, power = power,
                         mk.sel = mk.sel_i)

        Z1 <- diag(length(pheno))
        ETA <- list( list(Z=Z1, K=K))
        ans <- sommer::GWAS(Y = pheno, X = X, Z = ETA, W = geno_i, method = "EMMA",
                    silent = !verbose, gwas.plots = FALSE)
        p.val <- c(p.val, ans$W.scores$score[1, ])

      }

    }

    G.res <- data.frame(map, p.val)
    colnames(G.res) <- c("mk.id", "Chrom", "Position", "p.val")

  }

  return(G.res)

}
