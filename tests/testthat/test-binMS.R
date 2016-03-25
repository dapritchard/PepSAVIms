
# `````````````````````````````````````````````````````` #
#  Construct some baseline values used to generate data  #
# ...................................................... #

## Set of values 1 (numbered 1-4)

# Create a baseline set of values
info_1 <- list(mass = 5000,
               chg  = 5,
               time = 20)
info_1$mtoz <- info_1$mass / info_1$chg

# Similar to 1 but charges are differenct
info_2 <- list(mass = 5000,
               chg  = 4,
               time = 20)
info_2$mtoz <- info_2$mass / info_2$chg

# Similar to 1 but time is different
info_3 <- list(mass = 5000,
               chg  = 5,
               time = 25)
info_3$mtoz <- info_3$mass / info_3$chg

# Similar to 1 but mass is different
info_4 <- list(mass = 6000,
               chg  = 5,
               time = 20)
info_4$mtoz <- info_4$mass / info_4$chg


## Set of values 2 (numbered 5-8)

# Create a baseline set of values
info_5 <- list(mass = 9000,
               chg  = 5,
               time = 20)
info_5$mtoz <- info_5$mass / info_5$chg

# Similar to 5 but charges are differenct
info_6 <- list(mass = 9000,
               chg  = 4,
               time = 20)
info_6$mtoz <- info_6$mass / info_6$chg

# Similar to 5 but time is different
info_7 <- list(mass = 9000,
               chg  = 5,
               time = 25)
info_7$mtoz <- info_7$mass / info_7$chg

#Similar to 5 but mass is different
info_8 <- list(mass = 8000,
               chg  = 5,
               time = 20)
info_8$mtoz <- info_8$mass / info_8$chg


## Set of values 3 (numbered 9-12)

# Charges out-of-bounds
info_9 <- list(mass = 9000,
               chg  = 1,
               time = 20)
info_9$mtoz <- info_9$mass / info_9$chg

# Time out-of-bounds
info_10 <- list(mass = 9000,
                chg  = 5,
                time = 0)
info_10$mtoz <- info_10$mass / info_10$chg

# Mass out-of-bounds
info_11 <- list(mass = 0,
                chg  = 5,
                time = 20)
info_11$mtoz <- info_11$mass / info_11$chg




# ``````````````````````````````````````````` #
#  Construct function to sample observations  #
# ........................................... #

# The relationship b/w mtoz and mass is: mass = (m/z - 1.007825) * charge

sampObs <- function(base, sd_mtoz=0.01, sd_time=0.2) {

  mtoz <- base$mtoz + rnorm(1, sd=sd_mtoz)
  chg <- base$chg
  time <- base$time + rnorm(1, sd=sd_time)
  mass <- (mtoz - 1.007825) * chg

  list(mtoz=mtoz, chg=chg, time=time, mass=mass)
}




# ````````````````````````` #
#  Sample data for testing  #
# ......................... #

# 1:11 corresponds to the 11 base values made earlier
sampleList <- lapply(1:11, function(x) {
  base_nm <- paste0("info_", x)
  replicate(5, sampObs(get(base_nm)))
})
testMS <- matrix(unlist( t( Reduce(cbind, sampleList) ) ), ncol=4)
testMS <- as.matrix( cbind(testMS, rep(1, nrow(testMS))) )
colnames(testMS) <- c("mtoz", "chg", "time", "mass", "ms")




# ```````````````````````````````````````` #
#  Calculate the binned data for test set  #
# ........................................ #

idx <- lapply(1:8, function(j) seq(5 * (j - 1) + 1, 5 * j))
binDat <- t( sapply(idx, function(j) colMeans(testMS[j, ])) )
binDat <- binDat[order(binDat[, "mtoz"], binDat[, "time"], binDat[, "chg"]), ]

msObj <- msDat(binDat[, "ms", drop=FALSE], binDat[, "mtoz"], binDat[, "chg"])

out <- binMS(testMS, "mtoz", "chg", "mass", "time", "ms", c(14, 45), c(2e3, 15e3), c(2, 10),
             0.05, 1)
