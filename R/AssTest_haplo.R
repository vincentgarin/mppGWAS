#################
# AssTest_haplo #
#################

#' Haplotype Association Test
#'
#' Perform an haplotype genome wide association test using a kinship matrix
#' to correct for the genetic background.
#'
#' The function is an adaptation of the Fixed-A and Random-A models proposed
#' by Wei and Xu (2016). We modified the code of the function \code{magicScan}
#' from the \code{magicQTL} package
#' (\url{http://www.genetics.org/content/202/2/471.supplemental})
#' (Wei and Xu, 2016). We replace the founder QTL term by an haplotype QTL term.
#' As in \code{magicQTL}, the QTL term can be fitted as fixed or random. The
#' model is fitted using the EMMA algorithm proposed by Kang et al. (2008).
#'
#' A correction for the genetic background is applied using a kinship matrix.
#' By default, the kinship is computed using the method of Astle and Balding
#' (2009). It is possible to compute a linkage disequilibrium adjusted kinship
#' (LDAK) kinship matrix using the method of Speed et al. (2012) by introducing
#' the weights computed with the function \code{\link{LDAK_weights}}. The model
#' can be fitted using the kinship containing all markers or removing
#' the markers of the scanned chromosome (\code{K_i = TRUE}).
#'
#' @param haplo.block \code{List} containing the required data to compute an
#' haplotype model. This object can be obtained using the function
#' \code{\link{haplo_blocks}}.
#'
#' @param haplo.term \code{Character} variable indicating if the haplotype term
#' should be computed as fixed ("fixed") or random ("random").
#' Default = "fixed".
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2 and family
#' \strong{containing the list of markers that will be used to compute the
#' kinship matrix}
#'
#' @param map Three columns \code{data.frame} corresponding to the list of
#' markers present in \code{gp}. The three columns represent the marker
#' identifier, the chromosome, and the marker position (cM or bp).
#'
#' @param trait \code{Numerical} vector of phenotypic trait values.
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
#' @param K_i \code{Logical} specifying if the kinship correction should be done
#' by removing the markers of the scanned chromosome. Default = TRUE.
#'
#'
#' @return Return:
#'
#' \item{res}{Data.frame with in row the haplotype blocks and in column the
#' LRT and Wald statistics of each haplotype blocks.}
#'
#' @author Vincent Garin, Julong Wei and Shizhong Xu
#' (original underlying functions)
#'
#' @references
#'
#' Astle, W., & Balding, D. J. (2009). Population structure and cryptic
#' relatedness in genetic association studies. Statistical Science, 451-471.
#'
#' Kang, H. M., Zaitlen, N. A., Wade, C. M., Kirby, A., Heckerman, D., Daly,
#' M. J., & Eskin, E. (2008). Efficient control of population structure in
#' model organism association mapping. Genetics, 178(3), 1709-1723.
#'
#' Speed, D., Hemani, G., Johnson, M. R., & Balding, D. J. (2012).
#' Improved heritability estimation from genome-wide SNPs. The American Journal
#' of Human Genetics, 91(6), 1011-1021.
#'
#' Wei, J., & Xu, S. (2016). A random-model approach to QTL mapping in
#' multiparent advanced generation intercross (MAGIC) populations.
#' Genetics, 202(2), 471-486.
#'
#' @import data.table
#' @import SKAT
#' @import zoo
#' @importFrom utils read.table setTxtProgressBar tail txtProgressBar write.table
#' @importFrom sommer mmer
#' @importFrom stats as.formula model.matrix optim pchisq
#' @importFrom qqman manhattan qq
#'
#' @export
#'

############# arguments

# # haplo.block production
#
# path <- "/home/vincent/Haplo_GRM/EUNAM"
# setwd(path)
#
# # genotypic data
#
# load("./data/geno/geno_gp_imp.RData")
#
# # map
#
# map <- read.table("./data/map/map_imp.txt", h = TRUE, stringsAsFactors = FALSE)
# map <- map[, c(1, 3, 4, 6)]
# colnames(map) <- c("mk.id", "chr", "cM", "bp")
#
# pheno <- read.csv("./data/pheno/Adj_means.csv", row.names = 1)
#
# hp.block <- haplo_blocks(gp = gp.imp, map = map, hap = 3, hap.unit = 1)
#
# haplo.term <- "fixed"
# trait <- pheno[, 1]
# gp <- gp.imp
# map <- map
# weights = NULL
# power = -1
# K_i = FALSE
#
# source('~/Haplo_GRM/mppAssTest/R/magicScan_mod.R')
# source('~/Haplo_GRM/mppAssTest/R/mixedPar.R')
# source('~/Haplo_GRM/mppAssTest/R/cal.xx.R')
# source('~/Haplo_GRM/mppAssTest/R/Fpoly_mod.R')
# source('~/Haplo_GRM/mppAssTest/R/Rpoly_mod.R')

AssTest_haplo <- function(haplo.block, haplo.term = "fixed", gp, map, trait,
                          weights = NULL, power = -1, K_i = TRUE){

  # 1. Checks
  ###########

  # check that the list of haplotype in  haplo.map is the same as the list of the
  # haplotypes in haplo.geno

  # if(!identical(rownames(haplo.geno), haplo.map[, 1])){
  #
  #   stop(paste("The list of marker in the haplotype matrix (haplo.geno) ",
  #              "and in the haplotype map (haplo.map) are not stricly equivalent",
  #              "(same content, same order)."))
  #
  # }

  # check that the list of markers in the map is the same as the list of the
  # genotype matrix of the gp object

  if(!identical(colnames(gp$geno), map[, 1])){

    stop(paste("The list of marker in the gp object ",
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

  ######################### end checks

  # 2. Process the data for the magicScan_mod function
  ####################################################

  if(haplo.term == "random"){model <- "Random-A"
  } else if (haplo.term == "fixed")  {model <- "Fixed-A"}

  # 2.1 datafram

  d <- data.frame(trait, IncMat_cross(gp$covar$family))
  colnames(d) <- c("trait", paste("cr", 1:(dim(d)[2] - 1)) )


  # 3. Model computation
  ######################

  if(!K_i){ # Use the whole genome for K.

    ### 3.1 Preparation of the kinship matrix

    K <- mpp_kinship(gp = gp, weights = weights, power = power)

    kk.eigen <- list()
    kk.eigen[[1]] <- K
    kk.eigen[[2]] <- eigen(K)

    scan <- magicScan_mod(dataframe = d, gen = haplo.block[[1]],
                          map = haplo.block[[2]], kk.eigen = kk.eigen,
                          nalleles = haplo.block[[3]], model = model)

    res <- lapply(1:length(haplo.block[[2]]), function(i){ return(scan[[i]]) })
    res <- do.call(rbind, res)


  } else { # Remove the kth chromosome for the computation of K.

    # K all markers - ith chromosome

    res <- c()
    n.chr <- length(haplo.block[[2]])
    chr.id <- unique(map[, 2])

    mk.sel_temp <- map

    for(i in 1:n.chr){

      mk_chr_i <- map[map[, 2] == chr.id[i], 1]

      # adapt the list of selected markers

      mk.sel_i <- mk.sel_temp[mk.sel_temp[, 2] != chr.id[i], 1]

      K <- mpp_kinship(gp = gp, weights = weights, power = power,
                       mk.sel = mk.sel_i)

      kk.eigen <- list()
      kk.eigen[[1]] <- K
      kk.eigen[[2]] <- eigen(K)

      scan <- magicScan_mod(dataframe = d, gen = haplo.block[[1]][i],
                            map = haplo.block[[2]][i], kk.eigen = kk.eigen,
                            nalleles = haplo.block[[3]][i],
                            model = model)


      res <- rbind(res, scan[[1]])

    }

  }

  return(res)

}
