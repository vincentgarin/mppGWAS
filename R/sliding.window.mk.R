#####################
# sliding.window.mk #
#####################

# function to map a partition of a map using regular sliding blocks
# with a fixed width in marker eventually spaced by gap.


sliding.window.mk <- function(map, hap, gap){

  # Partition per chromosome

  set.id <- c()
  snp.id <- c()
  map_plot <- c() # position to use for the plotting

  for(i in 1:10){

    map_i <- map[map[, 2] == i, ]
    nb_mk <- dim(map_i)[1]

    partition <- rollapply(data = 1:nb_mk, width = hap, by = gap, FUN = function(x) x)
    map_plot <- rbind(map_plot, map_i[partition[, round(hap/2, 0)], c(2, 4)])

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