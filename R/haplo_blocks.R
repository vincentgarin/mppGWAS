################
# haplo_blocks #
################

#' Determine haplotype blocks and haplotypes
#'
#' Determine haplotype block and haplotypes. The haplotype blocks can be
#' determined based on several criteria: number of markers, distance in cM or
#' LD based (not for the moment). The blocks can be adjacent (non-overlapping)
#' or they can be determined using a sliding window.
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2.
#'
#' @param map Four columns \code{data.frame} with marker id, chromosome,
#' getic position in cM and physical position in bp.
#'
#' @param hap Defines either a number of markers (if \code{hap.unit=1}) or a
#' genetic distance in cM (if \code{hap.unit=2}) that should be used for building
#' haplo blocks. Default = 1.
#'
#' @param hap.unit Set to 1 for building haplo blocks based on number of
#' markers, set to 2 for building haplo blocks based on markers lying in a
#' specified stretch of the chromosome (in cM). Default = 1.
#'
#' @param sliding.window \code{logical} value specifying if
#' the haplotype blocks should be construct using a sliding window.
#' Default = FALSE.
#'
#' @param gap \code{Numeric} value indicating the number of markers between the
#' center of two consecutive haplotype blocks. Default = 1.
#'
#' @param LD.based NOT AVAILABLE NOW. \code{logical} value specifying if the
#' haplotype blocks should be constructed based on LD. Default = FALSE.
#'
#' @param ld.threshold NOT AVAILABLE NOW.
#'
#' @return Return a list of the three following lists:
#'
#' \item{gen}{List of haplotype incidence matrix per chromosome. For each
#' chromosome, the haplotype incidences matrices of all haplotype block are
#' agregated in a single numeric matrix. The column represent the genotype
#' and the row the haplotype of the different blocks. The total number of row
#' is equal to the sum of the number of different haplotype per block. For
#' example if the chromosome contain two haplotype block with respectively 3 and
#' 6 haplotype, then the number or row will be 9 (3 + 6). For each haplotype
#' block, the incidence matrix specifies the haplotype owned by the different
#' genotypes using 0 1 number.}
#'
#' \item{map}{List of haplotype block map per chromosome. For each chromosome,
#' the list contain a four column data frame with: the haplotype block
#' identifier, the chromosome, the position in cM and in bp. The positions
#' are the average cM and bp values of the marker present in the haplotype
#' block.}
#'
#' \item{nalleles}{List of vector per chrsomosome. For each chromosome, the
#' length of the vector correspond to the number of haplotype blocks contained
#' on the current chromosome. For each haplotype block it specify the number of
#' different haplotypes.}
#'
#'
#' @author Vincent Garin
#'
#' @export
#'

############ arguments

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
# gp <- gp.imp
# # out.dir <- "/home/vincent/Haplo_GRM/software/SelectionTools"
# hap <- 1
# hap.unit <- 2
#
# sliding.window = TRUE
# gap = 1
# LD.based = FALSE
# ld.threshold = 0.9
#
# # sub-functions
#
# source('~/Haplo_GRM/mppAssTest/R/sliding.window.cM.R')
# source('~/Haplo_GRM/mppAssTest/R/sliding.window.mk.R')
# source('~/Haplo_GRM/mppAssTest/R/adj.block.cM.R')
# source('~/Haplo_GRM/mppAssTest/R/adj.block.mk.R')


haplo_blocks <- function(gp, map, hap = 1, hap.unit = 1, sliding.window = FALSE,
                         gap = 1, LD.based = FALSE, ld.threshold = 0.9){

  # 1. Make a marker partition
  ############################

  if(LD.based){

    stop("LD.based option is not available for the moment")

    # ld <- st.calc.ld ( ld.measure="r2",
    #                    data.set="geno", auxfiles = FALSE)
    #
    # h <- st.def.hblocks ( ld.threshold = 0.8,
    #                       tolerance = 3,
    #                       ld.criterion = "flanking",
    #                       data.set="geno" )

  } else {

    if(sliding.window){

      if(hap.unit == 1){# mk number

        SNP_set <- sliding.window.mk(map = map, hap = hap, gap = gap)

      } else if (hap.unit == 2){ # cM distance

        SNP_set <- sliding.window.cM(map = map, hap = hap, gap = gap)

      }

    } else { # no sliding window: adjacent blocks

      if(hap.unit == 1){

        SNP_set <- adj.block.mk(map = map, hap = hap)

      } else if (hap.unit == 2){

        SNP_set <- adj.block.cM(map = map, hap = hap)

      }

    }

  }

  # 2. Convert the markers block into haplotypes
  ##############################################

  chr.ind <- unlist(lapply(X = SNP_set[, 1],
                    FUN = function(x) strsplit(x = x, split = "_")[[1]][1]))

  chr.id <- unique(chr.ind)

  # Initialise empty lists to store the haplotype information

  gen <- vector(mode = "list", length = length(chr.id))
  map.hp <- vector(mode = "list", length = length(chr.id))
  nalleles <- vector(mode = "list", length = length(chr.id))
  geno.names <- paste0("g", 1:dim(gp$geno)[1])

  for(i in 1:length(chr.id)){

  # selec data of the ith chromosome

    SNP_set_i <- SNP_set[chr.ind  == chr.id[i], ]
    block.id <- unique(SNP_set_i[, 1])

    gen_i <- vector(mode = "list", length = length(block.id))
    n.allele_i <- rep(0, length(block.id))
    cM_i <- rep(0, length(block.id))
    bp_i <- rep(0, length(block.id))

    for(j in 1:length(block.id)){

      mk.j <- SNP_set_i[SNP_set_i[, 1] == block.id[j], 2]
      geno.j <- gp$geno[, mk.j, drop = FALSE]
      cM_i[j] <- mean(map[map[, 1] %in% mk.j, 3])
      bp_i[j] <- mean(map[map[, 1] %in% mk.j, 4])

      haplo.list <- apply(X = geno.j, MARGIN = 1,
                          FUN = function(x) paste(x, collapse = ""))

      hap.mat.j <- t(model.matrix(~ -1 + haplo.list))
      colnames(hap.mat.j) <- geno.names
      n.al <- dim(hap.mat.j)[1]
      rownames(hap.mat.j) <- paste0("h", j, "_", 1:n.al)
      n.allele_i[j] <- n.al

      gen_i[[j]] <- data.frame(hap.mat.j)


    }

    gen[[i]] <- rbindlist(gen_i)

    map.hp_i <- data.frame(block.id, rep(i, length(block.id)),
                           cM_i, bp_i, stringsAsFactors = FALSE)
    colnames(map.hp_i) <- c("markers", "chr", "cm", "bp")

    map.hp[[i]] <- map.hp_i
    nalleles[[i]] <- n.allele_i

  }

  return(list(gen, map.hp, nalleles))

}
