################
# adj.block.cM #
################

adj.block.cM <- function(map, hap){

  # Partition per chromosome

  set.id <- c()
  snp.id <- c()

  chr.id <- unique(map[, 2])

  for(i in 1:length(chr.id)){

    map_i <- map[map[, 2] == chr.id[i], ]
    nb_mk <- dim(map_i)[1]
    mk.names <- map_i[, 1]

    list.mk_i <- c()
    set.id_i <- c()
    map_plot_i <- c()

    dist.cM <- map_i[, 3]
    ref.dist <- map_i[, 3][1]
    block.const = TRUE
    j <- 1


    while(block.const){

      dist.cM_j <- dist.cM - ref.dist

      sel.vect <- (dist.cM_j >= 0) & (dist.cM_j <= hap)

      if(sum(sel.vect) > 0){

        list.mk_i[[j]] <- map_i[sel.vect, 1]
        set.id_i <- c(set.id_i, paste0("chr", i, "_",
                                       rep(j, length(list.mk_i[[j]]))))

        # update the variables

        last.mk <- tail(list.mk_i[[j]], n=1)
        ind.mk <- which(mk.names == last.mk)

        if(!(ind.mk >= dim(map_i)[1])){

          ref.dist <- map_i[ind.mk + 1, 3]
          j <- j + 1

        } else {block.const <- FALSE}


      } else {block.const <- FALSE}

    }

    set.id <- c(set.id, set.id_i)
    snp.id <- c(snp.id, unlist(list.mk_i))

  }

  SetID.file <- data.frame(set.id, snp.id, stringsAsFactors = FALSE)

  return(SetID.file)

}
