################
# haplo_blocks #
################

#' Determine haplotype blocks and haplotypes
#'
#' Determine haplotype block and haplotypes using the SelectionTools package
#' (\url{http://www.population-genetics.de/~frisch-m/}).
#' The haplotype blocks can be determined based on several criteria:
#' number of markers, distance in cM or LD based.
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2.
#'
#' @param map Four columns \code{data.frame} with marker id, chromosome,
#' getic position in cM and physical position in bp.
#'
#' @param out.dir output directory where the ldak weights will be stored
#'
#' @param hap Defines either a number of markers (if \code{hap.unit=1}) or a
#' genetic distance in cM (if \code{hap.unit=2}) that should be used for building
#' haplo blocks. Default = 1.
#'
#' @param hap.unit Set to 1 for building haplo blocks based on number of
#' markers, set to 2 for building haplo blocks based on markers lying in a
#' specified stretch of the chromosome (in cM). Default = 1.
#'
#' @param sliding.window NOT AVAILABLE NOW. \code{logical} value specifying if
#' the haplotype blocks should be construct using a sliding window. ... .
#' Default = FALSE.
#'
#' [...]
#'
#' @param LD.based NOT AVAILABLE NOW. \code{logical} value specifying if the
#' haplotype blocks should be constructed based on LD. Default = FALSE.
#'
#' @param ld.threshold NOT AVAILABLE NOW.
#'
#' @return Return a list containing the following elements:
#'
#' \item{haplo.geno}{Data frame with in row the haplotype block and in column
#' the genotypes. For each haplotype, the haplotype alleles are coded with
#' a different haplotype score.}
#'
#' \item{haplo.map}{Four column map data frame with the haplotpype block id
#' the chromsome, the cM and the bp position of the haplotype block.}
#'
#' @author Vincent Garin
#'
#' @export
#'

############ arguments

# gp <- gp.imp
# map <- map[, c(1, 3, 4, 6)]
# out.dir <- "/home/vincent/Haplo_GRM/software/SelectionTools"
# hap <- 3
# hap.unit <- 1


haplo_blocks <- function(gp, map, out.dir, hap = 1, hap.unit = 1,
                         sliding.window = FALSE, LD.based = FALSE,
                         ld.threshold = 0.9){

  old.wd <- getwd()

  # create a temporary directories to store all the intermediary files

  temp.dir <- file.path(out.dir, "temp_dir")
  temp.dir.in <- file.path(out.dir, "temp_dir/input")
  temp.dir.out <- file.path(out.dir, "temp_dir/output")

  system(paste("mkdir", temp.dir))
  system(paste("mkdir", temp.dir.in))
  system(paste("mkdir", temp.dir.out))

 # Format the genotype information

  geno <- gp.imp$geno
  geno <- as.matrix(geno)
  geno <- t(geno)
  geno[geno == 0] <- "1/1"
  geno[geno == 1] <- "1/2"
  geno[geno == 2] <- "2/2"
  geno[is.na(geno)] <- "-1/-1"

  # save geno to the temporary directory

  file.pop <- file.path(temp.dir.in, "geno.pop")

  write.table(x = geno, file = file.pop, col.names = TRUE, row.names = TRUE,
              quote = FALSE)

  # Format the map information

  map.hp <- map[, 1:3]
  colnames(map.hp) <- c("name", "chrom", "pos")

  file.map <- file.path(temp.dir.in, "Map.map")

  write.table(x = map.hp, file = file.map, col.names = TRUE, row.names = FALSE,
              quote = FALSE)

  # Load SelectionTools data

  setwd(temp.dir)
  # st.input.dir <- "input"
  assign(x = "st.input.dir" , value = "input", envir = .GlobalEnv)
  # st.output.dir <- "output"
  assign(x = "st.output.dir" , value = "output", envir = .GlobalEnv)
  st.read.marker.data ("geno.pop", format = "m")
  st.read.map ("Map.map", format = "mcp", skip = 1)
  st.copy.marker.data ("geno", "default" )

  # make the haplotype partition

  if(LD.based){ # LD based haplotype definition

    stop("LD.based option is not available for the moment")

    # ld <- st.calc.ld ( ld.measure="r2",
    #                    data.set="geno", auxfiles = FALSE)
    #
    # h <- st.def.hblocks ( ld.threshold = 0.8,
    #                       tolerance = 3,
    #                       ld.criterion = "flanking",
    #                       data.set="geno" )

  } else {

    if(sliding.window){ # sliding window fixed mk nb or cM distance

      stop("sliding.window option is not available for the moment")

      # # sliding window partition based on marker numbers
      #
      # h.list <- haplo.list.mk(map = map, n.mk = hap)
      # h <- st.set.hblocks(haplotype.list = h.list, hap.symbol = "c",
      #                     data.set = "geno")
      #
      # # sliding windo partition based on cM distance
      #
      # h.list <- ?
      # h <- st.set.hblocks(haplotype.list = h.list, hap.symbol = "c",
      #                     data.set = "geno")

    } else { # non-overlapping mk or cM blocks

      h <- st.def.hblocks(hap = hap, hap.unit = hap.unit, data.set = "geno")

    }

  }

  st.recode.hil (data.set = "geno") # Haplotype alleles
  x <- st.marker.data.statistics("geno")

  haplo.alleles <- x$marker.list
  haplo.geno <- x$genotypes
  rownames(haplo.geno) <- haplo.geno[, 1]
  haplo.geno <- haplo.geno[, -1]
  haplo.geno <- as.matrix(haplo.geno)

  # recode the map

  # get the bp position

  s <- strsplit(as.character(h$Markers), ";")
  bp.pos <- rep(0, length(s))
  for(i in 1:length(s)){ ############################ problem there with the indices of s[[i]]

    bp.pos[i] <- mean(map[map[, 1] %in% s[[i]], 4])

  }

  haplo.map <- data.frame(as.character(h$Name), h$Chrom, h$Pos, bp.pos,
                          stringsAsFactors = FALSE)
  colnames(haplo.map) <- c("haplo.id", "chr", "cM", "bp")


  # re-establish the old wd

  setwd(old.wd)

  # delete the temporary directory

  system(paste("rm -rf", temp.dir))

  haplo <- list(haplo.geno = haplo.geno, haplo.map = haplo.map)

  return(haplo)

}
