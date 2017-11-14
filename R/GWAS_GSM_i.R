##############
# GWAS_GSM_i #
##############

GWAS_GSM_i <- function(i, SSD_file, chr.id, mk.sel_temp, gp, weights,
                       power, null.mod.form, dataset, kernel, weights.beta,
                       weights.kernel){

  SSD.INFO <- SKAT::Open_SSD(SSD_file[[1]], SSD_file[[2]])

  SSD.INFO_i <- SSD.INFO
  SSD.INFO_i$SetInfo <- SSD.INFO_i$SetInfo[SSD_file[[3]][, 2] == chr.id[i], ]
  SSD.INFO_i$nSets <- sum(SSD_file[[3]][, 2] == chr.id[i])

  # adapt the list of selected markers

  mk.sel_i <- mk.sel_temp[mk.sel_temp[, 2] != chr.id[i], 1]

  K <- mpp_kinship(gp = gp, weights = weights, power = power,
                   mk.sel = mk.sel_i)


  # compute the NULL model cr + kinship (EMMA_x)

  obj <- SKAT::SKAT_NULL_emmaX(formula = as.formula(null.mod.form), data = dataset,
                         K = K, llim = -10, ulim = 10)

  out_i <- SKAT.SSD.All_Ki(SSD.INFO = SSD.INFO_i, obj, kernel = kernel,
                           weights.beta = weights.beta,
                           obj.SNPWeight = weights.kernel, verbose = FALSE)


  res_i <- data.frame(SSD_file[[3]][SSD_file[[3]][, 2] == chr.id[i], ],
                      -log10(out_i$results[, 2]), stringsAsFactors = FALSE)

  colnames(res_i) <- c("mk.id", "Chrom", "Position", "p.val")

  SKAT::Close_SSD()


  res_i

}
