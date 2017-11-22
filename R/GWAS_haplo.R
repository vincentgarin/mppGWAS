##############
# GWAS_haplo #
##############

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
#' the markers of the scanned chromosome (\code{K_i = TRUE}). In the sitation,
#' where the markers of the scanned chromosomes are removed (\code{K_i = TRUE}),
#' if the computation failed for one or several chromosomes, their results
#' are set to 0 and a warning message is produced at the end of the scan.
#'
#' @param haplo.block Object of class \code{hap_bl} obtained with the function
#' \code{\link{haplo_blocks}} containing the haplotype matrices necessary to
#' compute the haplotype model.
#'
#' @param haplo.term \code{Character} variable indicating if the haplotype term
#' should be computed as fixed ("fixed") or random ("random").
#' Default = "fixed".
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2, map with
#' marker position in cM, phenotype and family. \strong{containing the list of
#' markers that will be used to compute the kinship matrix}.
#'
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the gp object should be used. Default = 1.
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
#' formula .Default = -1.
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
#'
#' @return Return:
#'
#' \item{G_res}{Object of class \code{G_res} representing a data.frame with four
#' columns: marker identifier, chromosome, position in cM and -log10(p-value).}
#'
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
#' @examples
#'
#' # GWAS haplotype
#'
#' data("EUNAM_gp")
#' data("EUNAM_LD_weights")
#'
#' # for gp object construction and weights computation see examples of
#' # LDAK_weights()
#'
#' hap_bl <- haplo_blocks(gp = EUNAM_gp, hap = 3, hap.unit = 1,
#'                        sliding.window = FALSE, gap = 1)
#'
#' res <- GWAS_haplo(haplo.block = hap_bl, haplo.term = "fixed", gp = EUNAM_gp,
#'                   weights = EUNAM_LD_weights, K_i = TRUE, verbose = FALSE)
#'
#' plot_manhattan(res)
#' plot_qq(res)
#'
#' @import data.table
#' @import doParallel
#' @import foreach
#' @import parallel
#' @import SKAT
#' @import stats
#' @import synbreed
#' @import utils
#' @import zoo
#' @importFrom utils read.table setTxtProgressBar tail txtProgressBar write.table
#' @importFrom sommer GWAS
#' @importFrom stats as.formula model.matrix optim pchisq
#' @importFrom qqman manhattan qq
#' @importFrom foreach %dopar%
#'
#' @export
#'

############# arguments

# haplo.block production

# haplo.term <- "fixed"
# trait <- pheno[, 1]
# gp <- gp.imp
# map <- map
# weights = NULL
# power = -1
# K_i = FALSE
#
# source('~/Haplo_GRM/mppGWAS/R/magicScan_mod.R')
# source('~/Haplo_GRM/mppGWAS/R/mixedPar.R')
# source('~/Haplo_GRM/mppGWAS/R/cal.xx.R')
# source('~/Haplo_GRM/mppGWAS/R/Fpoly_mod.R')
# source('~/Haplo_GRM/mppGWAS/R/Rpoly_mod.R')

GWAS_haplo <- function(haplo.block, haplo.term = "fixed", gp, trait = 1,
                       weights = NULL, power = -1, K_i = TRUE, n.cores = NULL,
                       verbose = FALSE){

  # 1. Checks
  ###########

  # test haplo.block

  if(!is_hap_bl(haplo.block)){

    stop(paste("The haplo.block object is not of class hap_bl. Such an object",
               "can be obtained with the function haplo_blocks()"))

  }

  # test haplo.term

  if (!(haplo.term %in% c("fixed", "random"))){

    stop("haplo.term must be 'fixed' or 'random'.")

  }

  # check gp, trait and weights

  check_gpData(gp)

  check_weights(weights = weights, gp = gp)

  check_trait(trait = trait, gp = gp)

  ######## end checks

  # 2. Process the data for the magicScan_mod function
  ####################################################

  # map

  map <- data.frame(rownames(gp$map), gp$map, stringsAsFactors = FALSE)
  colnames(map) <- c("mk.id", "chr", "cM")

  # pheno

  if(is.numeric(trait)){

    pheno <- gp$pheno[, 1, trait]

  } else {

    trait.names <- attr(gp$pheno, "dimnames")[[2]]
    pheno <- gp$pheno[, 1, which(trait %in% trait.names)]

  }

  if(haplo.term == "random"){model <- "Random-A"
  } else if (haplo.term == "fixed")  {model <- "Fixed-A"}

  # 2.1 datafram

  d <- data.frame(pheno, IncMat_cross(gp$covar$family))
  colnames(d) <- c("trait", paste("cr", 1:(dim(d)[2] - 1)) )

  ################### stop there: parallel, classed object

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
                          nalleles = haplo.block[[3]], model = model,
                          verbose = verbose)

    res <- lapply(1:length(haplo.block[[2]]), function(i){ return(scan[[i]]) })
    res <- do.call(rbind, res)


  } else { # Remove the kth chromosome for the computation of K.

    # K all markers - ith chromosome

    res <- c()
    n.chr <- length(haplo.block[[2]])
    chr.id <- unique(map[, 2])
    mk.sel_temp <- map
    ind.failed <- c()


    if(!is.null(n.cores)){

      cl <- parallel::makeCluster(n.cores)
      doParallel::registerDoParallel(cl)

      res <- foreach::foreach(i=1:n.chr) %dopar% {

        GWAS_haplo_i(i = i, chr.id = chr.id, map = map, mk.sel_temp = mk.sel_temp,
                     gp = gp, weights = weights, power = power, d = d,
                     haplo.block = haplo.block, haplo.term = haplo.term,
                     model = model, verbose = verbose)

      }

      parallel::stopCluster(cl)

      ind.failed <- unlist(lapply(X = res, FUN = function(x) x[[2]]))
      res <- lapply(X = res, FUN = function(x) x[[1]])
      res <- do.call(what = rbind, res)

    } else {

      # regular execution of K-i

      for(i in 1:n.chr){

        mk_chr_i <- map[map[, 2] == chr.id[i], 1]

        # adapt the list of selected markers

        mk.sel_i <- mk.sel_temp[mk.sel_temp[, 2] != chr.id[i], 1]

        K <- mpp_kinship(gp = gp, weights = weights, power = power,
                         mk.sel = mk.sel_i)

        kk.eigen <- list()
        kk.eigen[[1]] <- K
        kk.eigen[[2]] <- eigen(K)

        scan <- tryCatch(magicScan_mod(dataframe = d, gen = haplo.block[[1]][i],
                                       map = haplo.block[[2]][i], kk.eigen = kk.eigen,
                                       nalleles = haplo.block[[3]][i],
                                       model = model, verbose = verbose),
                         error = function(e) NULL)

        if(!is.null(scan)){

          res <- rbind(res, scan[[1]])
          ind.failed <- c(ind.failed, NULL)

        } else {

          ind.failed <- c(ind.failed, i)
          map.info <- haplo.block[[2]][[i]][, c(2, 3)]
          n.pos <- dim(map.info)[1]
          map.info <- data.frame(1:n.pos, map.info)

          if(haplo.term == "fixed"){

            res_i <- data.frame(map.info, matrix(0, n.pos, 7))
            colnames(res_i) <- c("Num", "chr", "ccM", "lrt", "lrt.p", "lrt.logp",
                                 "wald", "wald.p", "wald.logp", "sigma2")

          } else if (haplo.term == "random"){

            res_i <- data.frame(map.info, matrix(0, n.pos, 10))
            colnames(res_i) <- c("Num", "chr", "ccM", "lrt", "lrt.p", "lrt.logp",
                                 "wald", "wald.p", "wald.logp", "tau_k", "sigma2",
                                 "lam_k", "conv")

          }

          res <- rbind(res, res_i)

        }

      }

    }

    if(!is.null(ind.failed)){

      mess <- paste("The model computation failed for chromosome(s):",
                    paste(ind.failed, collapse = ", "))

      warning(mess)

    }

  }



  bl_nm <- paste0("chr", res[, 2], "_", res[, 1])

  if(haplo.term == "fixed"){

    G_res <- data.frame(bl_nm, res[, c(2, 3, 9)])

  } else if (haplo.term == "random") {

    G_res <- data.frame(bl_nm, res[, c(2, 3, 6)])

  }

  colnames(G_res) <- c("mk.id", "Chrom", "Position", "p.val")
  class(G_res) <- c("data.frame", "G_res")

  return(G_res)

}
