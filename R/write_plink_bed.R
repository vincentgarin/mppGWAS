###################
# write_plink_bed #
###################

#' Produce plink bed data
#'
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2 and family
#'
#' @param map Four columns \code{data.frame} with chromosome, marker id,
#' getic position in cM and physical position in bp.
#'
#' @param out.dir output directory where the .bed files will be saved.
#'
#' @param prefix \code{Character}. Prefix for the plink input file
#'
#' @param plink.dir directory where the plink executable program is located.
#'
#' @return Return:
#'
#' the function save the plink .bed files (.bed, .bim and .fam) in the directory
#' specified in argument \code{out.dir}
#'
#' @author Vincent Garin
#'
#' @export
#'


write_plink_bed <- function(gp, map, out.dir, prefix, plink.dir){

  geno <- gp$geno

  # checks

  if(!identical(colnames(geno), map[, 1])){

    stop("The markers in the map and the genotype matrix are not identical.")

  }

  # .ped file

  fam.id <- as.character(gp.imp$covar$family)
  gen.id <- as.character(gp.imp$covar$id)
  pat.id <- mat.id <- rep(0, dim(geno)[1])
  sex <- rep(9, dim(geno)[1])
  pheno <- rep(-9, dim(geno)[1])

  geno[geno == 0] <- "A A"
  geno[geno == 1] <- "A B"
  geno[geno == 2] <- "B B"

  ped_file <- data.frame(fam.id, gen.id, pat.id, mat.id, sex, pheno, geno,
                         stringsAsFactors = FALSE)

  # save .ped file

  write.table(ped_file, file = file.path(out.dir, paste0(prefix, ".ped")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  # .map file

  map_file <- data.frame(map[, 2], map[, 1], map[, 3], map[, 4],
                         stringsAsFactors = FALSE)

  write.table(map_file, file = file.path(out.dir, paste0(prefix, ".map")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  # convert the .ped and .map into bed files

  file.loc <- file.path(out.dir, prefix)
  soft.cmd <- file.path(plink.dir, "plink")
  cmd1 <- paste(soft.cmd, "--file", file.loc, "--make-bed --out", file.loc,
                "--noweb")

  system(cmd1)

  # remove the .ped amd .map files

  system(paste("rm", file.path(out.dir, paste0(prefix, ".ped"))))

  system(paste("rm", file.path(out.dir, paste0(prefix, ".map"))))

}
