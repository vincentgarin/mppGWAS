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

    stop("LD.based option is not available for the moment")

  } else {

    if(sliding.window){

      if(hap.unit == 1){# mk number

        SetID.file <- sliding.window.mk(map = map, hap = hap, gap = gap)

      } else if (hap.unit == 2){ # cM distance

        SetID.file <- sliding.window.cM(map = map, hap = hap, gap = gap)

      }

    } else { # no sliding window: adjacent blocks

      if(hap.unit == 1){

        SetID.file <- adj.block.mk(map = map, hap = hap)


      } else if (hap.unit == 2){

        SetID.file <- adj.block.cM(map = map, hap = hap)

      }

    }

  }

  # produce the corresponding map

  set_id <- unique(SetID.file$set.id)

  chr.ind <- unlist(lapply(X = set_id,
                           FUN = function(x) strsplit(x = x, split = "_")[[1]][1]))

  chr.id <- unique(chr.ind)

  chr <- rep(1:length(chr.id), time = table(factor(chr.ind, levels = chr.id)))

  cM <- rep(0, length(set_id))
  bp <- rep(0, length(set_id))

  for(j in 1:length(set_id)){

    mk.j <- SetID.file[SetID.file[, 1] == set_id[j], 2]
    cM[j] <- mean(map[map[, 1] %in% mk.j, 3])
    bp[j] <- mean(map[map[, 1] %in% mk.j, 4])

  }

  map.hp <- data.frame(set_id, chr, cM, bp, stringsAsFactors = FALSE)


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


  return(list(File.SSD = File.SSD, File.Info = File.Info, map.hp = map.hp))

}
