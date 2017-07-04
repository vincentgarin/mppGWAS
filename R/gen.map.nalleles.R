####################
# gen.map.nalleles #
####################

# Function to process the data gen, map and nalleles data for the function
# magicScan_mod


gen.map.nalleles <- function(haplo.geno, haplo.map){

  chr.id <- unique(haplo.map[, 2])

  gen <- vector(mode = "list", length = length(chr.id))
  map.hp <- vector(mode = "list", length = length(chr.id))
  nalleles <- vector(mode = "list", length = length(chr.id))
  geno.names <- paste0("g", 1:dim(haplo.geno)[2])

  for(i in 1:length(chr.id)){

    map_i <- haplo.map[haplo.map[, 2] == chr.id[i], ]
    colnames(map_i) <- c("markers", "chr", "cm", "bp")
    map.hp[[i]] <- map_i

    # subset the genotypes

    geno_i <- haplo.geno[map_i$markers, ]
    gen_i <- c()
    n.allele_i <- rep(0, dim(geno_i)[1])

    for(j in 1:dim(geno_i)[1]){

      hap.mat.j <- t(model.matrix(~ -1 + geno_i[j, ]))
      colnames(hap.mat.j) <- geno.names
      n.al <- dim(hap.mat.j)[1]
      rownames(hap.mat.j) <- paste0("h", j, "_", 1:n.al)
      n.allele_i[j] <- n.al
      gen_i <- rbind(gen_i, hap.mat.j)

    }

    gen[[i]] <- gen_i
    nalleles[[i]] <- n.allele_i

  }

  return(list(gen, map.hp, nalleles))

}
