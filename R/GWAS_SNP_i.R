##############
# GWAS_SNP_i #
##############

# function to execute the GWAS 1SNP scan of a single chromosome. This function
# can be used to be called in parallel.

GWAS_SNP_i <- function(i, chr.id, gp, X, pheno, map, mk.sel_temp, weights,
                       power){

  mk_chr_i <- map[map[, 2] == chr.id[i], 1]

  mk.sel_i <- mk.sel_temp[mk.sel_temp[, 2] != chr.id[i], 1]

  geno_i <- gp$geno[, colnames(gp$geno) %in% mk_chr_i]

  K_i <- mpp_kinship(gp = gp, weights = weights, power = power,
                     mk.sel = mk.sel_i)

  Z1 <- diag(length(pheno))
  ETA <- list(list(Z=Z1, K=K_i))
  ans <- sommer::GWAS(Y = pheno, X = X, Z = ETA, W = geno_i, method = "EMMA",
                      silent = TRUE, gwas.plots = FALSE)
  ans$W.scores$score[1, ]


}
