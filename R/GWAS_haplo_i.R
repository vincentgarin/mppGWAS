################
# GWAS_haplo_i #
################

# argument


GWAS_haplo_i <- function(i, chr.id, map, mk.sel_temp, gp, weights, power, d,
                         haplo.block, model, verbose){

  mk_chr_i <- map[map[, 2] == chr.id[i], 1]

  # adapt the list of selected markers

  mk.sel_i <- mk.sel_temp[mk.sel_temp[, 2] != chr.id[i], 1]

  K <- mpp_kinship(gp = gp, weights = weights, power = power,
                   mk.sel = mk.sel_i)

  kk.eigen <- list()
  kk.eigen[[1]] <- K
  kk.eigen[[2]] <- eigen(K)

  scan <- tryCatch(magicScan_mod(dataframe = d, gen = haplo.block[[1]][i],
                                 map = haplo.block[[2]][i], kk.eigen = kk.eigen,
                                 nalleles = haplo.block[[3]][i],
                                 model = model, verbose = verbose),
                   error = function(e) NULL)

  if(!is.null(scan)){

    list(scan = scan[[1]], ind.failed = NULL)

  } else {

    map.info <- haplo.block[[2]][[i]][, c(2, 3)]
    n.pos <- dim(map.info)[1]
    map.info <- data.frame(1:n.pos, map.info)

    if(haplo.term == "fixed"){

      res_i <- data.frame(map.info, matrix(0, n.pos, 7))
      colnames(res_i) <- c("Num", "chr", "ccM", "lrt", "lrt.p", "lrt.logp",
                           "wald", "wald.p", "wald.logp", "sigma2")

    } else if (haplo.term == "random"){

      res_i <- data.frame(map.info, matrix(0, n.pos, 10))
      colnames(res_i) <- c("Num", "chr", "ccM", "lrt", "lrt.p", "lrt.logp",
                           "wald", "wald.p", "wald.logp", "tau_k", "sigma2",
                           "lam_k", "conv")

    }

    list(scan = res_i, ind.failed = i)

  }


}
