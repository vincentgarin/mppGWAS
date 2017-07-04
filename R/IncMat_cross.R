################
# IncMat_cross #
################

# Form the cross effect incidence matrix for the cross-specific intercept term

IncMat_cross <- function(cross.ind){

  if(length(unique(cross.ind)) == 1){

    cross.mat <- matrix(1, length(cross.ind), 1)
    colnames(cross.mat) <- paste0("Cr", cross.ind[1])

  } else {

    Cr <- factor(x = cross.ind, levels = unique(cross.ind))
    cross.mat <- model.matrix( ~ Cr - 1, data = Cr)

  }

  return(cross.mat)

}
