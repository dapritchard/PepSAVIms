
# ``````````` #
#  Load data  #
# ........... #

# Contains only the object testMS.  See Construct-Data-BinMS.R for the script
# used to create the data.

load("tests/testthat/Data-BinMS.RData")




# ```````````````````````````````````````` #
#  Calculate the binned data for test set  #
# ........................................ #

# Consider the following settings for consolidation:
#   peak retention time:  14 - 45
#   mass:  2000 - 15000
#   charge: 2 - 10
#   mass-to-charge difference:  0.05 Da
#   peak retention time difference:  1 minute
#   
# Under these settings, the data was constructed so that there are 8 bins, each
# consisting of a consecutive block of 5 observations.  I.e. the first bin
# consists of observations 1,...,5, the second bin consists of observations
# 6,...,10, etc. for the first 40 observations.  Observations 41,...,55 do not
# pass the inclusion criteria.

# Indices in data for each bin
idx <- lapply(1:8, function(j) seq(5 * (j - 1) + 1, 5 * j))
# Consolidate into bins
binDat <- t( sapply(idx, function(j) colMeans(testMS[j, ])) )
# Reorder the binned data
binDat <- binDat[order(binDat[, "mtoz"], binDat[, "time"], binDat[, "chg"]), ]
binDat[, "ms"] <- 5

# Construct msDat object
# msObj <- msDat(binDat[, "ms", drop=FALSE], binDat[, "mtoz"], binDat[, "chg"])
msObj <- msDat(binDat, "mtoz", "chg", "ms")
# Construct summary function information.  See Construct-Data-BinMS.R for where
# these values come from.
summ_info <- list(n_tot = 55L,
                    n_time_pr = 50L,
                    n_mass = 50L,
                    n_charge = 50L,
                    n_tiMaCh = 40L,
                    n_binned = 8L,
                    time_range = c(14, 45),
                    mass_range = c(2000, 15000),
                    charge_range = c(2, 10),
                    mtoz_diff = 0.05,
                    time_diff = 1)

trueBin <- list(msObj     = msObj,
                summ_info = summ_info)
class(trueBin) <- "binMS"







# ``````````````````` #
#    Begin testing    #
# ................... #

context("binMS method")

test_that("binMS: compare results with toy dataset", {
  
  out1 <- binMS(testMS, "mtoz", "chg", "mass", "time", "ms", c(14, 45), c(2e3, 15e3),
                c(2, 10), 0.05, 1)
  expect_equal( out1, trueBin )

  out2 <- binMS(testMS, 1, 2, "mass", "time", "ms", c(14, 45), c(2e3, 15e3),
             c(2, 10), 0.05, 1)
  expect_equal( out2, trueBin )

  out3 <- binMS(testMS, "mtoz", "chg", "mass", "time", NULL, c(14, 45), c(2e3, 15e3),
                c(2, 10), 0.05, 1)
  expect_equal( out3, trueBin )

  out4 <- binMS(data.frame(testMS), 1, 2, "mass", "time", NULL, c(14, 45), c(2e3, 15e3),
                c(2, 10), 0.05, 1)
  expect_equal( out4, trueBin )

  expect_warning( out5 <- binMS(data.frame(testMS), 1, 2, "mass", "time", NULL,
                                c(14, 45), c(2e3, 15e3), c(20, 100), 0.05, 1),
                  "No observations satisfied all of the inclusion criteria")
  expect_identical( out5$msObj, NULL )
})


# test_that("binMS: with invalid input", {
  

# })
