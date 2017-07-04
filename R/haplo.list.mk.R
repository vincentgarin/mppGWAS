#################
# haplo.list.mk #
#################

# Sub-function to partition the the genome using a sliding window and a fixed
# number of markers.


haplo.list.mk <- function(map, n.mk){

  haplo.list <- c()
  chr.id <- unique(map[, 2])
  n.chr <- length(chr.id)

  for(k in 1:n.chr){

    map_i <- map[map[, 2] == chr.id[k], ]
    nb_mk <- dim(map_i)[1]
    win.part <- rollapply(data = 1:nb_mk, width = n.mk, FUN = function(x) x)

    mk.vect_i <- rep("", dim(win.part)[1])
    pos_i <- rep(0, dim(win.part)[1])

    for(l in 1:dim(win.part)[1]){

      mk.vect_i[l] <- paste0(paste(map_i[win.part[l, ], 1], collapse = ";"), ";")
      pos_i[l] <- mean(map_i[win.part[l, ], 3])

    }

    chr_i <- rep(chr.id[k], dim(win.part)[1])
    name_i <- paste0("h", k, "_", 1:dim(win.part)[1])

    haplo.list_i <- data.frame(chr_i, pos_i, name_i,
                               rep("b", dim(win.part)[1]), mk.vect_i,
                               stringsAsFactors = FALSE)
    colnames(haplo.list_i) <- c("Chrom", "Pos", "Name", "Class", "Markers")
    haplo.list <- rbind(haplo.list, haplo.list_i)

  }


  return(haplo.list)

}
