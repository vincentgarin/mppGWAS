########################
# Test class functions #
########################

# gpData

is_gpData <- function(x){

  inherits(x = x, what = "gpData")

}

is_LD_wgh <- function(x){

  inherits(x = x, what = "LD_wgh")

}
