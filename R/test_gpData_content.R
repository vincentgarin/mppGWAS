#######################
# test_gpData_content #
#######################


test_gpData_content <- function(gp){

  if(!is_gpData(gp)){

    stop("The object given in gp is not a gpData object.")

  }

  # test geno

  if(is.null(gp$geno)){

    stop("The gpData object do not contain genotype data.")

  }

  # test map

  if(is.null(gp$map)){

    stop("The gpData object do not contain map data.")

  }

  # test pheno

  if(is.null(gp$pheno)){

    stop("The gpData object do not contain phenotype data.")

  }

  # test cross.ind

  if(any(is.na(gp$covar$family))){

    stop("The cross/family indicator of the gpData object is not valid.")

  }

}
