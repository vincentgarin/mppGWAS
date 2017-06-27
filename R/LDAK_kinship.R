################
# LDAK_kinship #
################

#' Compute linkage desequilibrium adjusted kinship
#'
#' Compute linkage disequilibrium adjusted kinship using the program LDAK.
#' (http://dougspeed.com/ldak/). The function is a wrapper for the LDAK
#' software.
#'
#' @param ldak.dir directory where the ldak executable program is located
#'
#' @param weights.loc path to the LDAK weights files (output of function
#' \code{\link{LDAK_weights}}). For example, /home/.../weights.
#'
#' @param bed.file.loc path to the .bed files (output of function
#' \code{\link{write_plink_bed}}) without .bed extension. For example,
#' /home/.../my_file.
#'
#' @param out.dir output directory where temporary file will be saved.
#' These files will be removed. Default = getwd().
#'
#' @param power \code{Numerical} value specifying the value of the
#' parameter for marker scores standardization. The column of the marker matrix
#' (X.j) are multiplied by var(X.j)^(power/2) .Default = -1.
#'
#' @param K_i \code{Numerical} value specifying a unique chromosome number that
#' should be removed from the kinship computation. By default \code{K_i = NULL},
#' which means that the kinship is computed using all markers.
#'
#' @param map If \code{K_i} is not NULL, \code{data.frame} map information with
#' at least a colum for marker identifier labeled \code{'mk.id'}, and one column
#' for chromosome indicator labeled \code{'chr'}. \strong{The marker identificer
#' must be the same as the one of the weight file.}. Default = NULL.
#'
#' @return Return:
#'
#' \item{K}{kinship matrix computed with the LDAK weights.}
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

# ldak.dir <- "/home/vincent/Haplo_GRM/software/LDAK"
# weights.loc <- "/home/vincent/Haplo_GRM/EUNAM/data/geno/LDAK_test/Test_LDAK_weights"
# bed.file.loc <- "/home/vincent/Haplo_GRM/EUNAM/data/geno/plink_files/Test"
# out.dir <- "/home/vincent/Haplo_GRM/EUNAM/data/geno/LDAK_test"
# power <- -1
# K_i <- 1
# map <- map

LDAK_kinship <- function(ldak.dir, weights.loc, bed.file.loc, out.dir = getwd(),
                         power = -1, K_i = NULL, map = NULL){

  # create a temporary output directory to store all the intermediary files

  temp.dir <- file.path(out.dir, "temp_dir")
  system(paste("mkdir", temp.dir))

  # save a copy of the weights in the temporary directory

  wgh.file2 <- file.path(temp.dir, "wgh")

  system(paste("cp", weights.loc, wgh.file2))

  if (!is.null(K_i)){

    wgh <- read.table(wgh.file2, header = TRUE, stringsAsFactors = FALSE)

    # Check the name of the marker in the map and weight file

    if(dim(map)[1] != dim(wgh)[1]){

      stop(paste("The list of marker in the map and in the weight file do not",
                 "have the same length."))

    }

    if((sum(map$mk.id %in% wgh[, 1]) != dim(map)[1])){

      stop(paste("The list of marker in the map and in the weight file are",
                 "different."))

    }

    # list of markers to set to zero

    mk.list <- map$mk.id[map$chr == K_i]
    wgh[wgh[, 1] %in% mk.list, 2] <- 0

    # save the modified weights

    write.table(x = wgh, file = wgh.file2, row.names = FALSE, quote = FALSE)

  }

  # kinship computation

  out.kin <- file.path(temp.dir, "kin")
  ldak.loc <- file.path(ldak.dir, "ldak5.beta")

  cmd3 <- paste(ldak.loc, "--calc-kins-direct", out.kin, "--weights",
                wgh.file2, "--bfile", bed.file.loc, paste("--power", power),
                "--kinship-raw YES")

  system(cmd3)

  # load kinship

  kin.loc <- file.path(temp.dir, "kin.grm.raw")
  K <- read.table( kin.loc)

  # delete temp directory

  system(paste("rm -rf", temp.dir))

  return(K)


}
