################
# adj.block.mk #
################

adj.block.mk <- function(map, hap){

  # Partition per chromosome

  set.id <- c()
  snp.id <- c()
  map_plot <- c() # position to use for the plotting

  chr.id <- unique(map[, 2])

  for(i in 1:length(chr.id)){

    map_i <- map[map[, 2] == chr.id[i], ]
    nb_mk <- dim(map_i)[1]

    partition <- rollapply(data = 1:nb_mk, width = hap, by = hap, FUN = function(x) x)
    map_plot <- rbind.data.frame(map_plot, map_i[partition[, round(hap/2, 0)],
                                                 c(1, 2, 4)],
                                 stringsAsFactors = FALSE)

    set.id_i <- rep(1:dim(partition)[1], each = hap)
    set.id_i <- paste0("chr", i, "_", set.id_i)
    set.id <- c(set.id, set.id_i)

    snp.id_i <- apply(X = partition, MARGIN = 1, FUN = function(x, names) names[x],
                      names = map_i[, 1])

    snp.id <- c(snp.id, c(snp.id_i))

  }

  SetID.file <- data.frame(set.id, snp.id, stringsAsFactors = FALSE)

  return(list(SetID.file, map_plot))

}
