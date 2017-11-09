################
# LDAK_weights #
################

#' Compute Linkage disequilibrium markers weights
#'
#' Compute linkage disequilibrium weights using the program LDAK.
#' The function is a wrapper for the LDAK software.
#'
#' This function use the LDAK 5 beta executable file that can be download here:
#' http://dougspeed.com/ldak/ and the plink executable file that can be downloaded
#' here: http://zzz.bwh.harvard.edu/plink/download.shtml.
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2, map with
#' marker position in cM, phenotype and family indicator.
#'
#' @param bp.pos \code{Numeric vector} of marker physical base pair positions.
#'
#' @param out.dir output directory where the ldak weights will be stored.
#' Defaut = working directory.
#'
#' @param plink.dir directory where the plink executable program is located.
#'
#' @param ldak.dir directory where the ldak executable program is located
#'
#' @param verbose \code{Logical} indicating if function outputs should be printed.
#' Default = FALSE
#'
#' @return Return:
#'
#' An object of class \code{LD_wgh} which is a \code{data.frame} with two columns
#' containing: the marker identifier and the computed LDAK weight.
#'
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


LDAK_weights <- function(gp, bp.pos, out.dir, plink.dir, ldak.dir,
                         verbose = FALSE){

  # test gpData and his content

  test_gpData_content(gp)

  # reform the map

  map <- data.frame(rownames(gp$map), gp$map, bp.pos, stringsAsFactors = FALSE)
  colnames(map) <- c("mk.id", "chr", "cM", "bp")

  # create a temporary output directory to store all the intermediary files

  temp.dir <- file.path(out.dir, "temp_dir")
  system(paste("mkdir", temp.dir))

  prefix <- "Pop"

  # compute the plink .bed file

  write_plink_bed(gp = gp, map = map, out.dir = temp.dir, prefix = prefix,
                  plink.dir = plink.dir, verbose = verbose)


 # compute the LDAK

 bed.file <- file.path(temp.dir, prefix)
 ldak.loc <- file.path(ldak.dir, "ldak5.beta")

 cmd1 <- paste(ldak.loc, "--cut-weights", temp.dir, "--bfile", bed.file)
 system(cmd1, intern = !verbose) # determine the sections

 cmd2 <- paste(ldak.loc, "--calc-weights-all", temp.dir, "--bfile", bed.file)
 system(cmd2, intern = !verbose) # compute the weights


 # get the results

 wgh.file <- file.path(temp.dir, "weights.all")

 weights <- read.table(wgh.file, header = TRUE, stringsAsFactors = FALSE)
 weights <- weights[, 1:2]

 class(weights) <- c("data.frame", "LD_wgh")

 # remove the temporary directory

 system(paste("rm -rf", temp.dir))

 return(weights)


}

