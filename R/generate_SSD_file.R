#####################
# generate_SSD_file #
#####################

#' Kernel association test SSD files
#'
#' Generate the SSD files for the kernel association test. The function
#' partitions the genome in adjacent or sliding window blocks based on
#' a fixed number of marker or cM distance.
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2.
#'
#' @param map Four columns \code{data.frame} with marker id, chromosome,
#' getic position in cM and physical position in bp.
#'
#' @param out.dir output directory where temporary files will be stored.
#'
#' @param prefix \code{Character} string. The SSD files will be saved as
#' prefix.SetID, prefix.SSD and prefix.SSD.info in \code{out.dir}.
#'
#' @param plink.dir directory where the plink executable program is located.
#'
#' @param hap Defines either a number of markers (if \code{hap.unit=1}) or a
#' genetic distance in cM (if \code{hap.unit=2}) that should be used for building
#' haplo blocks. Default = 1.
#'
#' @param hap.unit Set to 1 for building haplo blocks based on number of
#' markers, set to 2 for building haplo blocks based on markers lying in a
#' specified stretch of the chromosome (in cM). Default = 1.
#'
#'
#' @param sliding.window \code{logical} value specifying if
#' the haplotype blocks should be construct using a sliding window. If not,
#' the haplotype blocks are adjacent and not overlapping
#' (e.g. b1 = 123, b2 = 456, etc.). Default = TRUE.
#'
#' @param gap \code{Numeric} value indicating the number of markers between the
#' center of two consecutive haplotype blocks. Default = 1.
#'
#' @param LD.based NOT AVAILABLE NOW. \code{logical} value specifying if the
#' haplotype blocks should be constructed based on LD. Default = FALSE.
#'
#' @param ld.threshold NOT AVAILABLE NOW.
#'
#'
#'
#' @return
#'
#' Create the SSD files in the location specified in \code{out.dir}. return
#' also \code{map_plot} a two column data.frame with the chromosome and the
#' center bp position of each marker set.
#'
#' @author Vincent Garin
#'
#' @export
#'

############ arguments

# gp <- gp.imp
# map <- map[, c(1, 3, 4, 6)]
# out.dir <- "/home/vincent/Haplo_GRM/software/SKAT/data"
# prefix <- "Test"
# plink.dir <-  "/home/vincent/Haplo_GRM/software/PLINK"
#
# hap <- 1
# hap.unit <- 1
#
#
# sliding.window = TRUE
# gap <- 1
#
# LD.based = FALSE
# ld.threshold = 0.9


generate_SSD_file <- function(gp, map, out.dir, prefix, plink.dir, hap = 1,
                              hap.unit = 1, sliding.window = FALSE, gap = 1,
                              LD.based = FALSE, ld.threshold = 0.9){

  temp.dir <- file.path(out.dir, "temp_dir")
  system(paste("mkdir", temp.dir))

  # 1. Make a marker partition
  ############################

  if(LD.based){


  } else {

    if(sliding.window){

      if(hap.unit == 1){# mk number

        partition <- sliding.window.mk(map = map, hap = hap, gap = gap)

        SetID.file <- partition[[1]]
        map_plot <- partition[[2]]

      } else if (hap.unit == 2){ # cM distance

        partition <- sliding.window.cM(map = map, hap = hap, gap = gap)

        SetID.file <- partition[[1]]
        map_plot <- partition[[2]]

      }

    } else { # no sliding window: adjacent blocks

      if(hap.unit == 1){

        partition <- adj.block.mk(map = map, hap = hap)
        SetID.file <- partition[[1]]
        map_plot <- partition[[2]]

      } else if (hap.unit == 2){

        partition <- adj.block.cM(map = map, hap = hap)
        SetID.file <- partition[[1]]
        map_plot <- partition[[2]]

      }

    }

  }

  # save the SetID file

  File.SetID <- file.path(out.dir, paste0(prefix, ".SetID"))
  write.table(x = SetID.file, file = File.SetID, quote = FALSE, row.names = FALSE,
              col.names = FALSE)


  # 2. Produce the Plink files
  ###########################

  write_plink_bed(gp = gp, map = map, out.dir = temp.dir, prefix = prefix,
                  plink.dir = plink.dir)

  # 3. Generate the SSID files
  ############################

  # setwd("/home/vincent/Haplo_GRM/")


  # Create the MW File
  File.Bed <- file.path(temp.dir, paste0(prefix, ".bed"))
  File.Fam <- file.path(temp.dir, paste0(prefix, ".fam"))
  File.Bim <- file.path(temp.dir, paste0(prefix, ".bim"))

  File.SetID <- file.path(out.dir, paste0(prefix, ".SetID"))
  File.SSD <- file.path(out.dir, paste0(prefix, ".SSD"))
  File.Info <- file.path(out.dir, paste0(prefix, ".SSD.info"))

  Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD,
                     File.Info)


  # delete the temporary directory

  system(paste("rm -rf", temp.dir))


  return(map_plot)

}
