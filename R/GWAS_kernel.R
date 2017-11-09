############
# GWAS_GSM #
############

#' Kernel Association Test
#'
#' Perform a kernel genome wide association test using a kinship matrix
#' to correct for the genetic background.
#'
#' The first step of the analysis is to partition the genome into a certain
#' number of markers sets. Several strategies are possibles. \code{hap}
#' determine the size of the block or window. It can be defined in number of
#' markers (\code{hap.unit = 1}) or in cM (\code{hap.unit = 2}). The blocks
#' can be adjacent (no overlapping, e.g. b1 = 123, b2 = 456, etc.) or can
#' be deterimed by sliding window (\code{sliding.window = TRUE}). When user
#' choose sliding window it is also possible to define the number of marker
#' between the two centers of a window using \code{gap}.
#'
#' The genome scan is based on the sequence kernel association test (SKAT)
#' proposed by Wu et al. (2011) calling the SKAT package. This means that the
#' function determines the genetic relatedness between lines using a set of
#' markers and a kernel function. The variance covariance structure determined
#' by the kernel function is then used to test for association between the
#' phenotype and a tested region. The marker VCOV is assoiciated to a random
#' QTL term. The procedure is repeated at different positions along the genome
#' to perform a full genome scan.
#'
#' Three models are possible. The first one does not correct for the genetic
#' background using a kinship matrix (\code{kin.correct = FALSE}). The second
#' and the third models correct for the genetic background using a
#' kinship matrix (\code{kin.correct = TRUE}) that can include all markers
#' (\code{K_i = FALSE}) or remove the markers on the currently scanned
#' chromsome (\code{K_i = TRUE}). All models include a fixed cross-specific
#' intercept term. The model using the kinship correction are fitted using
#' the EMMA algorithm proposed by Kang et al. (2008).
#'
#' By default, the kinship is computed using the method of Astle and Balding
#' (2009). It is possible to compute a linkage disequilibrium adjusted kinship
#' (LDAK) kinship matrix using the method of Speed et al. (2012) by introducing
#' the weights computed with the function \code{\link{LDAK_weights}}. The model
#' can be fitted using the kinship containing all markers or removing
#' the markers of the scanned chromosome (\code{K_i = TRUE}).
#'
#' @param SSD_file file containing the genome partition in SNP set obtained
#' with the function \code{\link{generate_SSD_file}}.
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2 and family
#' \strong{containing the list of markers that will be used to compute the
#' kinship matrix}.
#'
#' @param map Four columns \code{data.frame} with marker id, chromosome,
#' getic position in cM and physical position in bp.
#'
#' @param trait \code{Numerical} vector of phenotypic trait values.
#'
#' @param kernel \code{Character} element indicating the type of kernel function
#' used to calculate the variance covariance structure using marker
#' set information. can e one of: "linear", "linear.weighted", "IBS",
#' "IBS.weighted", "quadratic" and "2wayIX". For details see details of the
#' function SKAT (pacakge SKAT). Default = "linear".
#'
#' @param weights.beta a numeric vector of parameters for the beta weights for
#' the weighted kernels. Default = c(1, 25).
#'
#' @param weights.kernel NOT AVAILABLE NOW. Vector of own provided weights for
#' for the kernel function. Default = NULL.
#'
#' @param kin.correct \code{Logical} value. If \code{kin.correct = TRUE}, the
#' model will use a kinship matrix to correct for the genetic background.
#' Default = TRUE.
#'
#' @param K_i \code{Logical} specifying if the kinship correction should be done
#' by removing the markers of the scanned chromosome. Default = TRUE.
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
#' formula. Default = -1.
#'
#'
#' @return Return:
#'
#' \item{res}{Data.frame with in row the haplotype blocks and in column the
#' LRT and Wald statistics of each haplotype blocks.}
#'
#' @author Vincent Garin
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
#' Wu, M. C., Lee, S., Cai, T., Li, Y., Boehnke, M., & Lin, X. (2011).
#' Rare-variant association testing for sequencing data with the sequence kernel
#' association test. The American Journal of Human Genetics, 89(1), 82-93.
#'
#' @export
#'

############# arguments


# gp <- gp.imp
# map <- map
# trait <- pheno[, 1]
# kernel <- "linear"
# weights.beta = c(1, 25)
# weights.kernel = NULL
# weights = weights
# power = -1
# K_i = TRUE

GWAS_GSM <- function(SSD_file, gp, map, trait, kernel = "linear",
                           weights.beta = c(1, 25), weights.kernel = NULL,
                           kin.correct = TRUE, K_i = TRUE, weights = NULL,
                           power = -1){

  # 1. Checks
  ###########

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


  ######################### end checks

  # 2. Open the SSD file
  ######################

  SSD.INFO <- Open_SSD(SSD_file[[1]], SSD_file[[2]])

  # 3. Computation of the genome scan
  ###################################

  # 3.1 cross-specific intercept term (used in all three models)

  cr.mat <- IncMat_cross(gp$covar$family)
  dataset <- data.frame(trait, cr.mat)
  colnames(dataset) <- c("trait", paste0("cr", 1:dim(cr.mat)[2]))
  null.mod.form <- paste("trait ~", paste(paste0("cr", 1:dim(cr.mat)[2]),
                                          collapse = "+"))


  if(!kin.correct){ # No kinship correction

    obj <- SKAT_Null_Model(formula = as.formula(null.mod.form), data = dataset,
                           out_type="C")

    out <- SKAT.SSD.All(SSD.INFO, obj, kernel = kernel,
                        weights.beta = weights.beta,
                        obj.SNPWeight = weights.kernel)

    res <- data.frame(SSD_file[[3]], -log10(out$results[, 2]),
                      stringsAsFactors = FALSE)
    colnames(res) <- c("mk.id", "Chrom", "Position", "p.val")

  } else {

    if(!K_i){ # Use the whole genome for K.

      # Preparation of the kinship matrix

      K <- mpp_kinship(gp = gp, weights = weights, power = power)

      # compute the NULL model cr + kinship (EMMA_x)

      obj <- SKAT_NULL_emmaX(formula = as.formula(null.mod.form), data = dataset,
                             K = K, llim = -10, ulim = 10)

      out <- SKAT.SSD.All(SSD.INFO, obj, kernel = kernel,
                          weights.beta = weights.beta,
                          obj.SNPWeight = weights.kernel)

      res <- data.frame(SSD_file[[3]], -log10(out$results[, 2]),
                        stringsAsFactors = FALSE)
      colnames(res) <- c("set.id", "chr", "cM", "bp", "p.val")


    } else { # Remove the kth chromosome for the computation of K.

      # make a temporary map to remove the markers on the considered chromsome.

      mk.sel_temp <- map

      res <- c()
      chr.id <- unique(mk.sel_temp[, 2])
      n.chr <- length(chr.id)

      for(i in 1:n.chr){

        # modify the SSD.INFO

        SSD.INFO_i <- SSD.INFO
        SSD.INFO_i$SetInfo <- SSD.INFO_i$SetInfo[SSD_file[[3]][, 2] == chr.id[i], ]
        SSD.INFO_i$nSets <- sum(SSD_file[[3]][, 2] == chr.id[i])

        # adapt the list of selected markers

        mk.sel_i <- mk.sel_temp[mk.sel_temp[, 2] != chr.id[i], 1]

        K <- mpp_kinship(gp = gp, weights = weights, power = power,
                         mk.sel = mk.sel_i)


        # compute the NULL model cr + kinship (EMMA_x)

        obj <- SKAT_NULL_emmaX(formula = as.formula(null.mod.form), data = dataset,
                               K = K, llim = -10, ulim = 10)

        out_i <- SKAT.SSD.All_Ki(SSD.INFO = SSD.INFO_i, obj, kernel = kernel,
                                 weights.beta = weights.beta,
                                 obj.SNPWeight = weights.kernel)


        res_i <- data.frame(SSD_file[[3]][SSD_file[[3]][, 2] == chr.id[i], ],
                            -log10(out_i$results[, 2]), stringsAsFactors = FALSE)

        colnames(res_i) <- c("set.id", "chr", "cM", "bp", "p.val")

        res <- rbind(res, res_i)

      }

    } # end K_i models

  }

  Close_SSD()

  # system(paste("rm -rf", temp.dir))


  return(res)

}
