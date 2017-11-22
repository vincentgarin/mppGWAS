############
# GWAS_GSM #
############

#' GSM Association Test
#'
#' Perform a genetic similarity matrix (GSM) genome wide association test using
#' a kinship matrix to correct for the genetic background.
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
#' @param SSD_file object of class \code{SSD_file} obtained with the function
#' \code{\link{generate_SSD_file}}.
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2, map with
#' marker position in cM, phenotype and family. \strong{containing the list of
#' markers that will be used to compute the kinship matrix}.
#'
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the gp object should be used. Default = 1.
#'
#' @param kernel \code{Character} element indicating the type of kernel function
#' used to calculate the variance covariance structure using marker
#' set information. can e one of: "linear", "linear.weighted", "IBS",
#' "IBS.weighted", "quadratic" and "2wayIX". For details see details of the
#' function SKAT (pacakge SKAT). Default = "linear.weighted".
#'
#' @param weights.beta a numeric vector of parameters for the beta weights for
#' the weighted kernels. Default = c(0.5, 0.5).
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
#' @param weights object of class \code{LD_wgh} obtained by the function
#' LDAK_weights() representing a data.frame with two columns: the marker
#' identifier and the LD adjusted weights. These weight will be used to compute
#' a LDAK as defined by Speed et al. (2012) for the genetic background
#' adjustement. Default = NULL. By default, the kinship matrix is computed
#' unweighted or all weights are equal to 1, which correspond to the Astle and
#' Balding kinship matrix.
#'
#' @param power \code{Numerical} value specifying the value of the
#' parameter for marker scores standardization. The column of the marker matrix
#' (X.j) are multiplied by var(X.j)^(power/2). It correspond to alpha in the
#' formula. Default = -1.
#'
#' @param n.cores \code{Numeric} value indicating the number of core to be used
#' if the user want to run the K-i (\code{K_i = TRUE}) scan in parallel.
#' Default = NULL.
#'
#' @param verbose \code{Logical} indicating if function outputs should be printed.
#' This argument is only used for the K-i model (\code{K_i = TRUE}).
#' Default = FALSE.
#'
#'
#' @return Return:
#'
#' \item{G_res}{Object of class \code{G_res} representing a data.frame with four
#' columns: marker identifier, chromosome, position in cM and -log10(p-value).}
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
#' @examples
#'
#' data("EUNAM_gp")
#' data("EUNAM_LD_weights")
#'
#' # for gp object construction and weights computation see examples of
#' # LDAK_weights()
#'
#' \dontrun{
#'
#'   plink.dir <- "/home/.../PLINK"
#'
#'   SSD_file <- generate_SSD_file(gp = EUNAM_gp, out.dir = getwd(), prefix = "Test",
#'                                 plink.dir = plink.dir, hap = 1, hap.unit = 2,
#'                                 sliding.window = TRUE)
#'
#'   res <- GWAS_GSM(SSD_file = SSD_file, gp = EUNAM_gp, kin.correct = TRUE,
#'                   weights = EUNAM_LD_weights, K_i = TRUE)
#'
#'   plot_manhattan(res)
#'   plot_qq(res)
#'
#' }
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

GWAS_GSM <- function(SSD_file, gp, trait = 1, kernel = "linear.weighted",
                     weights.beta = c(0.5, 0.5), weights.kernel = NULL,
                     kin.correct = TRUE, K_i = TRUE, weights = NULL,
                     power = -1, n.cores = NULL, verbose = FALSE){

  # 1. Checks
  ###########

  # test haplo.block

  if(!is_SSD_file(SSD_file)){

    stop(paste("The SSD_file object is not of class SSD_file. Such an object",
               "can be obtained with the function generate_SSD_file()"))

  }

  # check gp, trait and weights

  check_gpData(gp)

  check_weights(weights = weights, gp = gp)

  check_trait(trait = trait, gp = gp)

  ######################### end checks

  # 2. Open the SSD file
  ######################

  SSD.INFO <- Open_SSD(SSD_file[[1]], SSD_file[[2]])

  map <- data.frame(rownames(gp$map), gp$map, stringsAsFactors = FALSE)
  colnames(map) <- c("mk.id", "chr", "cM")

  if(is.numeric(trait)){

    pheno <- gp$pheno[, 1, trait]

  } else {
    trait.names <- attr(gp$pheno, "dimnames")[[2]]
    pheno <- gp$pheno[, 1, which(trait %in% trait.names)]

  }

  # 3. Computation of the genome scan
  ###################################

  # 3.1 cross-specific intercept term (used in all three models)

  cr.mat <- IncMat_cross(gp$covar$family)
  dataset <- data.frame(pheno, cr.mat)
  colnames(dataset) <- c("trait", paste0("cr", 1:dim(cr.mat)[2]))
  null.mod.form <- paste("trait ~", paste(paste0("cr", 1:dim(cr.mat)[2]),
                                          collapse = "+"))


  if(!kin.correct){ # No kinship correction

    obj <- SKAT_Null_Model(formula = as.formula(null.mod.form), data = dataset,
                           out_type="C")

    out <- SKAT.SSD.All(SSD.INFO, obj, kernel = kernel,
                        weights.beta = weights.beta,
                        obj.SNPWeight = weights.kernel)

    G_res <- data.frame(SSD_file[[3]], -log10(out$results[, 2]),
                        stringsAsFactors = FALSE)
    colnames(G_res) <- c("mk.id", "Chrom", "Position", "p.val")

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

      G_res <- data.frame(SSD_file[[3]], -log10(out$results[, 2]),
                          stringsAsFactors = FALSE)
      colnames(G_res) <- c("mk.id", "Chrom", "Position", "p.val")


    } else { # Remove the kth chromosome for the computation of K.

      # make a temporary map to remove the markers on the considered chromsome.

      mk.sel_temp <- map

      G_res <- c()
      chr.id <- unique(mk.sel_temp[, 2])
      n.chr <- length(chr.id)

      ######### multiple cores version

      if(!is.null(n.cores)){

        cl <- parallel::makeCluster(n.cores)
        doParallel::registerDoParallel(cl)

        res <- foreach::foreach(i=1:n.chr) %dopar% {

          GWAS_GSM_i(i = i, SSD_file = SSD_file,
                     chr.id = chr.id, mk.sel_temp = mk.sel_temp, gp = gp,
                     weights = weights, power = power,
                     null.mod.form = null.mod.form, dataset = dataset,
                     kernel = kernel, weights.beta = weights.beta,
                     weights.kernel = weights.kernel)

        }

        parallel::stopCluster(cl)

        G_res <- do.call(what = rbind, res)

      } else { # without multiple cores.

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
                                   obj.SNPWeight = weights.kernel,
                                   verbose = verbose)


          res_i <- data.frame(SSD_file[[3]][SSD_file[[3]][, 2] == chr.id[i], ],
                              -log10(out_i$results[, 2]), stringsAsFactors = FALSE)

          colnames(res_i) <- c("mk.id", "Chrom", "Position", "p.val")


          G_res <- rbind(G_res, res_i)

        }

      }



    } # end K_i models

  }

  Close_SSD()

  class(G_res) <- c("data.frame", "G_res")

  return(G_res)

}
