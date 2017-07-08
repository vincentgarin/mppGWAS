#####################
# sliding.window.cM #
#####################

# function to map a partition of a map using regular sliding blocks
# with a fixed width in cM.


sliding.window.cM <- function(map, hap, gap){

  # Partition per chromosome

  set.id <- c()
  snp.id <- c()
  map_plot <- c() # position to use for the plotting

  for(i in 1:10){

    map_i <- map[map[, 2] == i, ]
    dist.cM <- map_i[, 3]

    # determination of the point to stop

    dist.vect <- (dist.cM) - max(dist.cM)

    mk.stop <- map_i[( (dist.cM) - max(dist.cM) ) > -hap, 1][1]
    p.stop <- which(map_i[, 1] == mk.stop) -1

    win.cent <- seq(from = 1, to = p.stop, by = gap)
    list.mk_i <- vector(mode = "list", length = length(win.cent))
    set.id_i <- c()
    map_plot_i <- c()

    for(j in 1:length(win.cent)){

      dist.j <- dist.cM - dist.cM[win.cent[j]]
      list.mk_i[[j]] <- map_i[((dist.j >= 0) & (dist.j <= hap)), 1]
      set.id_i <- c(set.id_i, paste0("chr", i, "_",
                                     rep(j, length(list.mk_i[[j]]))))
      map_plot_i <- c(map_plot_i, mean(map_i[list.mk_i[[j]], 4]))

    }

    set.id <- c(set.id, set.id_i)
    snp.id <- c(snp.id, unlist(list.mk_i))
    map_plot <- rbind.data.frame(map_plot, cbind(unique(set.id_i),
                                      rep(i, length(map_plot_i)), map_plot_i),
                                 stringsAsFactors = FALSE)

  }

  SetID.file <- data.frame(set.id, snp.id, stringsAsFactors = FALSE)

  return(list(SetID.file, map_plot))

}
