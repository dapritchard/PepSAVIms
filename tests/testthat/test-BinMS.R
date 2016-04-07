
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

# Randomly permute the testMS data; the binning algorithm should be invariant
# to this
nr <- nrow(testMS)
testMS <- testMS[sample.int(nr, nr), ]

# Add some unused columns to testMS
testMS <- cbind(testMS, matrix(rnorm(2 * nr), ncol=2, dimnames=list(NULL, c("extr1", "extr2"))))

# Create a data.frame
testdf <- data.frame(testMS)
testdf$extr3 <- "a"




# ``````````````````` #
#    Begin testing    #
# ................... #

context("binMS method")

test_that("binMS: compare results using toy dataset", {
  
  out1 <- binMS(testMS, "mtoz", "chg", "mass", "time", "ms", c(14, 45), c(2e3, 15e3),
                c(2, 10), 0.05, 1)
  expect_equal( out1, trueBin )

  out2 <- binMS(testMS, 1, 2, "mass", "time", "ms", c(14, 45), c(2e3, 15e3),
             c(2, 10), 0.05, 1)
  expect_equal( out2, trueBin )

  out3 <- binMS(testMS[, 1:5], "mtoz", "chg", "mass", "time", NULL, c(14, 45), c(2e3, 15e3),
                c(2, 10), 0.05, 1)
  expect_equal( out3, trueBin )

  out4 <- binMS(testdf[, 1:5], 1, 2, "mass", "time", NULL, c(14, 45), c(2e3, 15e3),
                c(2, 10), 0.05, 1)
  expect_equal( out4, trueBin )

  expect_warning(out5 <- binMS(data.frame(testMS), 1, 2, "mass", "time", NULL,
                               c(14, 45), c(2e3, 15e3), c(20, 100), 0.05, 1),
                 "No observations satisfied all of the inclusion criteria")
  expect_identical( out5$msObj, NULL )

  out6 <- binMS(testdf, 1, 2, "mass", "time", "ms", c(14, 45), c(2e3, 15e3),
                c(2, 10), 0.05, 1)
  expect_equal( out6, trueBin )

  out7 <- binMS(testdf, testMS[, "mtoz"], 2, "mass", "time", "ms", c(14, 45), c(2e3, 15e3),
                c(2, 10), 0.05, 1)
  expect_equal( out7, trueBin )
})


test_that("binMS: with missing input", {
  
  expect_error( binMS(mtoz="mtoz", charge="chg", mass="mass", time_peak_reten="time",
                      ms_inten="ms", time_range=c(14, 45), mass_range=c(2e3, 15e3),
                      charge_range=c(2, 10), mtoz_diff=0.05, time_diff=1),
                "Must provide an argument for mass_spec" )
  expect_error( binMS(mass_spec=testMS, charge="chg", mass="mass", time_peak_reten="time",
                      ms_inten="ms", time_range=c(14, 45), mass_range=c(2e3, 15e3),
                      charge_range=c(2, 10), mtoz_diff=0.05, time_diff=1),
                "Must provide an argument for mtoz" )
  expect_error( binMS(mass_spec=testMS, mtoz="mtoz", mass="mass", time_peak_reten="time",
                      ms_inten="ms", time_range=c(14, 45), mass_range=c(2e3, 15e3),
                      charge_range=c(2, 10), mtoz_diff=0.05, time_diff=1),
               "Must provide an argument for charge" )
  expect_error( binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                      ms_inten="ms", time_range=c(14, 45), mass_range=c(2e3, 15e3),
                      charge_range=c(2, 10), mtoz_diff=0.05, time_diff=1),
                "Must provide an argument for time_peak_reten" )
  expect_error( binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                      time_peak_reten="time", ms_inten="ms", mass_range=c(2e3, 15e3),
                      charge_range=c(2, 10), mtoz_diff=0.05, time_diff=1),
                "Must provide an argument for time_range" )
  expect_error( binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                      time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                      charge_range=c(2, 10), mtoz_diff=0.05, time_diff=1),
                "Must provide an argument for mass_range" )
  expect_error( binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                      time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                      mass_range=c(2e3, 15e3), mtoz_diff=0.05, time_diff=1),
                "Must provide an argument for charge_range" )
  expect_error( binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                      time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                      mass_range=c(2e3, 15e3), charge_range=c(2, 10), time_diff=1),
                "Must provide an argument for mtoz_diff" )
  expect_error( binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                      time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                      mass_range=c(2e3, 15e3), charge_range=c(2, 10), mtoz_diff=0.05),
                "Must provide an argument for time_diff" )
})


test_that("binMS: wrong types of arguments", {
  
  expect_error(binMS(mass_spec=list(), mtoz="mtoz", charge="chg", mass="mass",
                     time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                     mass_range=c(2e3, 15e3), charge_range=c(2, 10), mtoz_diff=0.05,
                     time_diff=1),
               "mass_spec must be either a matrix or data.frame")
  
  expect_error(binMS(mass_spec=testMS, mtoz=list(), charge="chg", mass="mass",
                     time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                     mass_range=c(2e3, 15e3), charge_range=c(2, 10), mtoz_diff=0.05,
                     time_diff=1),
               "mtoz must be either of mode numeric or character")

  expect_error(binMS(mass_spec=testMS, mtoz="mtoz", charge=list(), mass="mass",
                     time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                     mass_range=c(2e3, 15e3), charge_range=c(2, 10), mtoz_diff=0.05,
                     time_diff=1),
               "charge must be either of mode numeric or character")

  expect_error(binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass=list(),
                     time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                     mass_range=c(2e3, 15e3), charge_range=c(2, 10), mtoz_diff=0.05,
                     time_diff=1),
               "mass must be either NULL or of mode numeric or character")

  expect_error(binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                     time_peak_reten=list(), ms_inten="ms", time_range=c(14, 45),
                     mass_range=c(2e3, 15e3), charge_range=c(2, 10), mtoz_diff=0.05,
                     time_diff=1),
               "time_peak_reten must be either of mode numeric or character")
  
  expect_error(binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                     time_peak_reten="time", ms_inten=list(), time_range=c(14, 45),
                     mass_range=c(2e3, 15e3), charge_range=c(2, 10), mtoz_diff=0.05,
                     time_diff=1),
               "ms_inten must be either NULL or of mode numeric or character")

  expect_error(binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                     time_peak_reten="time", ms_inten="ms", time_range=list(),
                     mass_range=c(2e3, 15e3), charge_range=c(2, 10), mtoz_diff=0.05,
                     time_diff=1),
               "time_range must be of mode numeric")

  expect_error(binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                     time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                     mass_range=list(), charge_range=c(2, 10), mtoz_diff=0.05,
                     time_diff=1),
               "mass_range must be of mode numeric")

  expect_error(binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                     time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                     mass_range=c(2e3, 15e3), charge_range=list(), mtoz_diff=0.05,
                     time_diff=1),
               "charge_range must be of mode numeric")

  expect_error(binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                     time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                     mass_range=c(2e3, 15e3), charge_range=c(2, 10), mtoz_diff=list(),
                     time_diff=1),
               "mtoz_diff must be of mode numeric")

  expect_error(binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                     time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                     mass_range=c(2e3, 15e3), charge_range=c(2, 10), mtoz_diff=0.05,
                     time_diff=list()),
               "time_diff must be of mode numeric")
})




test_that("binMS: arguments of the right type but illegal values", {
  
  expect_error(binMS(mass_spec=testMS, mtoz="nocolumn", charge="chg", mass="mass",
                     time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                     mass_range=c(2e3, 15e3), charge_range=c(2, 10), mtoz_diff=1,
                     time_diff=1),
               "name provided not in data - nocolumn element in mtoz")

  expect_error(binMS(mass_spec=testMS, mtoz="mtoz", charge=543, mass="mass",
                     time_peak_reten="time", ms_inten="ms", time_range=c(14, 45),
                     mass_range=c(2e3, 15e3), charge_range=c(2, 10), mtoz_diff=1,
                     time_diff=1),
               "out of bounds value provided for charge")

  expect_error(binMS(mass_spec=testMS, mtoz="mtoz", charge="chg", mass="mass",
                     time_peak_reten="time", ms_inten="ms", time_range=c(45, 45),
                     mass_range=c(2e3, 15e3), charge_range=c(2, 10), mtoz_diff=0.05,
                     time_diff=1),
               "The values of time_range must be in increasing order")
})
