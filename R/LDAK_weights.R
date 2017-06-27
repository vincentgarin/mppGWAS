################
# LDAK_weights #
################

#' Compute Linkage disequilibrium markers weights
#'
#' Compute linkage disequilibrium weights using the program LDAK.
#' (http://dougspeed.com/ldak/). The function is a wrapper for the LDAK
#' software.
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2 and family
#'
#' @param map Four columns \code{data.frame} with chromosome, marker id,
#' getic position in cM and physical position in bp. with column names being
#' "mk.id", "chr", "cM", "bp".
#'
#' @param out.dir output directory where the ldak weights will be stored
#'
#' @param prefix \code{Character}. Prefix for the plink input file
#'
#' @param plink.dir directory where the plink executable program is located.
#'
#' @param ldak.dir directory where the ldak executable program is located
#'
#' @param kinship \code{Logical} value specifying if the kinshi matrix should be
#' returned. Default = FALSE.
#'
#' @return Return:
#'
#' the function save the LDAK weights in the directory specified in argument
#' \code{out.dir}.
#'
#' If \code{kinship = TRUE}, kinship matrix computed with the LDAK weights.
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
# out.dir <- "/home/vincent/Haplo_GRM/EUNAM/data/geno/LDAK_test"
# gp <- gp.imp
# prefix = "Test"
# plink.dir = "/home/vincent/Haplo_GRM/software/PLINK"

LDAK_weights <- function(gp, map, out.dir, prefix, plink.dir, ldak.dir,
                         kinship = FALSE){

  # create a temporary output directory to store all the intermediary files

  temp.dir <- file.path(out.dir, "temp_dir")
  system(paste("mkdir", temp.dir))

  # compute the plink .bed file

  write_plink_bed(gp = gp, map = map, out.dir = temp.dir, prefix = prefix,
                  plink.dir = plink.dir)


 # compute the LDAK

 bed.file <- file.path(temp.dir, prefix)
 ldak.loc <- file.path(ldak.dir, "ldak5.beta")

 cmd1 <- paste(ldak.loc, "--cut-weights", temp.dir, "--bfile", bed.file)
 system(cmd1) # determine the sections

 cmd2 <- paste(ldak.loc, "--calc-weights-all", temp.dir, "--bfile", bed.file)
 system(cmd2) # compute the weights

 ##################### computation of the kinship K matrix

 if(kinship){

   out.kin <- file.path(temp.dir, paste(prefix, "kin", sep = "_"))
   wght.loc <- file.path(temp.dir, "weights.all")

   cmd3 <- paste(ldak.loc, "--calc-kins-direct", out.kin, "--weights", wght.loc,
                 "--bfile", bed.file, "--power -1 --kinship-raw YES")

   system(cmd3)

   # move the kinship to the raw directory

   kin.file <- file.path(temp.dir, paste(prefix, "kin.grm.raw", sep = "_"))
   kin.file2 <- file.path(out.dir, paste(prefix, "kinship.grm.raw", sep = "_"))

   system(paste("cp", kin.file, kin.file2))

 }

 ################################# end computation of the K matrix

 # re-organise the results

 wgh.file <- file.path(temp.dir, "weights.all")
 wgh.file2 <- file.path(out.dir, paste(prefix, "LDAK", "weights", sep = "_"))

 system(paste("cp", wgh.file, wgh.file2))

 # remove the useless files

 # cmd3 <- paste("cd", out.dir, "\n", paste("find . ! -name", paste0("'",
 # paste(prefix, "LDAK", "weights", sep = "_"), "'"),"-type f -exec rm -f {} +"))

 # remove the temporary directory

 system(paste("rm -rf", temp.dir))


}

