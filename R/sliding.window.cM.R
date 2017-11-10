#####################
# sliding.window.cM #
#####################

# function to map a partition of a map using regular sliding blocks
# with a fixed width in cM.


sliding.window.cM <- function(map, hap, gap){

  # Partition per chromosome

  set.id <- c()
  snp.id <- c()

  chr.id <- unique(map[, 2])

  for(i in 1:length(chr.id)){

    map_i <- map[map[, 2] == chr.id[i], ]
    dist.cM <- map_i[, 3]

    # determination of the point to stop

    dist.vect <- (dist.cM) - max(dist.cM)

    mk.stop <- map_i[( (dist.cM) - max(dist.cM) ) > -hap, 1][1]
    p.stop <- which(map_i[, 1] == mk.stop) -1

    win.cent <- seq(from = 1, to = p.stop, by = gap)
    list.mk_i <- vector(mode = "list", length = length(win.cent))
    set.id_i <- c()

    for(j in 1:length(win.cent)){

      dist.j <- dist.cM - dist.cM[win.cent[j]]
      list.mk_i[[j]] <- map_i[((dist.j >= 0) & (dist.j <= hap)), 1]
      set.id_i <- c(set.id_i, paste0("chr", i, "_",
                                     rep(j, length(list.mk_i[[j]]))))

    }

    set.id <- c(set.id, set.id_i)
    snp.id <- c(snp.id, unlist(list.mk_i))


  }

  SetID.file <- data.frame(set.id, snp.id, stringsAsFactors = FALSE)

  return(SetID.file)

}
