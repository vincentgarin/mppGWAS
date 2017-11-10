###############################
# check format of the objects #
###############################

# functions to check that the passed object have the correct format

# check gpData

check_gpData <- function(gp){

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

# check trait

check_trait <- function(trait, gp){

  if(!(is.character(trait) || is.numeric(trait))){

    stop("The trait object must be either numeric or character")

  }

  if(is.numeric(trait)){

    nb.trait <- dim(gp$pheno)[3]

    if(!((0 < trait) && (trait <=nb.trait))){

      stop(paste("The trait indicator should be between 1 and", nb.trait,
                 "the total number of traits in the gp object"))

    }

  }

  if(is.character(trait)){

    trait.names <- attr(gp$pheno, "dimnames")[[2]]

    if (!(trait %in% trait.names)){

      stop(paste("trait must be one of:", trait.names))

    }

  }

}


# check weights

check_weights <- function(weights, gp){

  if(!is.null(weights)){

    if(!is_LD_wgh(weights)){

      stop("The LD weights are not of class LD_wgh obtained with LDAK_weights().")

    }

  }

  # check that there is a weight for each marker if weights are given

  if(!is.null(weights)){

    if(!(sum(colnames(gp$geno) %in% weights[, 1]) == dim(gp$geno)[2])){

      stop(paste("The weight argument does not contain a value for each marker",
                 "in the gpData object (gp)."))

    }
  }

}
