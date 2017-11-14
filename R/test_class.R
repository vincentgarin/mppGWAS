########################
# Test class functions #
########################

# gpData

is_gpData <- function(x){

  inherits(x = x, what = "gpData")

}

# LD_wgh

is_LD_wgh <- function(x){

  inherits(x = x, what = "LD_wgh")

}

# hap_bl

is_hap_bl <- function(x){

  inherits(x = x, what = "hap_bl")

}

# SSD_file

is_SSD_file <- function(x){

  inherits(x = x, what = "SSD_file")

}
