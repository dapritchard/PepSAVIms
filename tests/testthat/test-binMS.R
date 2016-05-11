
# ``````````` #
#  Load data  #
# ........... #

load("tests/testthat/Data-BinMS.RData")




# ```````````````````````````````````````` #
#  Valid input with various types of data  #
# ........................................ #


context("binMS method")

test_that("binMS: compare results with various methods of providing / specifying data", {

  # Specify data by name
  out <- binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
               c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1)
  expect_equal( out, trueBin )

  # Specify data by index
  out <- binMS(testMS, 1, 2, 4, 3, c(5, 6), c(14, 45),
               c(2e3, 15e3), c(2, 10), 0.05, 1)
  expect_equal( out, trueBin )

  # Provide data as vectors
  out <- binMS(testMS, testMS[, "mtoz"], testMS[, "chg"], testMS[, "mass"], testMS[, "time"],
               c("ms1", "ms2"), c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1)
  expect_equal( out, trueBin )

  # Specify data by a mix of indices and names
  out <- binMS(testMS, 1, 2, "mass", "time", c("ms1", "ms2"), c(14, 45),
               c(2e3, 15e3), c(2, 10), 0.05, 1)
  expect_equal( out, trueBin )

  # Specify / provide data by a mix of indices and names and data vectors
  out <- binMS(testMS, 1, 2, "mass", testMS[, "time"], c("ms1", "ms2"), c(14, 45),
               c(2e3, 15e3), c(2, 10), 0.05, 1)
  expect_equal( out, trueBin )

  # Specify ms data by NULL
  out <- binMS(testMS[, 1:n_datacols], "mtoz", "chg", "mass", "time", NULL,
               c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1)
  expect_equal( out, trueBin )

  # Specify data by a mix of indices and names and ms data by NULL
  out <- binMS(testMS[, 1:n_datacols], 1, 2, "mass", "time", NULL, c(14, 45),
               c(2e3, 15e3), c(2, 10), 0.05, 1)
  expect_equal( out, trueBin )
  
  # Test data as a data.frame and specify data by name
  out <- binMS(testdf, 1, 2, "mass", "time", c("ms1", "ms2"), c(14, 45),
               c(2e3, 15e3), c(2, 10), 0.05, 1)
  expect_equal( out, trueBin )

  # Test data as a data.frame and specify data by index
  out <- binMS(testdf, 1, 2, 4, 3, c(5, 6), c(14, 45),
               c(2e3, 15e3), c(2, 10), 0.05, 1)
  expect_equal( out, trueBin )

  # Test data as a data.frame and provide data as vectors
  out <- binMS(testdf, testdf[, "mtoz"], testdf[, "chg"], testdf[, "mass"], testdf[, "time"],
               c("ms1", "ms2"), c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1)
  expect_equal( out, trueBin )

  # Test data as a data.frame and specify data by a mix of indices and names
  out <- binMS(testdf, 1, 2, "mass", "time", c("ms1", "ms2"), c(14, 45),
               c(2e3, 15e3), c(2, 10), 0.05, 1)
  expect_equal( out, trueBin )

  # Test data as a data.frame and specify / provide data by a mix of indices and
  # names and data vectors
  out <- binMS(testdf, 1, 2, "mass", testdf[, "time"], c("ms1", "ms2"), c(14, 45),
               c(2e3, 15e3), c(2, 10), 0.05, 1)
  expect_equal( out, trueBin )

  # Specify minimum value for charge such that no observations meets criterion
  expect_warning(out <- binMS(testMS, 1, 2, "mass", "time", NULL,
                              c(14, 45), c(2e3, 15e3), c(20, 100), 0.05, 1),
                 "No observations satisfied all of the inclusion criteria")
  expect_identical( out$msObj, NULL )
})




# ``````````````````````````` #
#  Missing param for formals  #
# ........................... #


# Note that the formal arguments mass and ms_inten have defaults

test_that("binMS: with missing params", {

  # Missing mass_spec
  expect_error( binMS( , "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                      c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
                "Must provide an argument for mass_spec" )
  
  # Missing mtoz
  expect_error( binMS(testMS, , "chg", "mass", "time", c("ms1", "ms2"),
                      c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
                "Must provide an argument for mtoz" )
  
  # Missing charge
  expect_error( binMS(testMS, "mtoz", , "mass", "time", c("ms1", "ms2"),
                      c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "Must provide an argument for charge" )
  
  # Missing time_peak_reten
  expect_error( binMS(testMS, "mtoz", "chg", "mass", , c("ms1", "ms2"),
                      c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
                "Must provide an argument for time_peak_reten" )

  # Missing time_range
  expect_error( binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                      , c(2e3, 15e3), c(2, 10), 0.05, 1),
                "Must provide an argument for time_range" )

  # Missing mass_range
  expect_error( binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                      c(14, 45), , c(2, 10), 0.05, 1),
                "Must provide an argument for mass_range" )

  # Missing charge_range
  expect_error( binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                      c(14, 45), c(2e3, 15e3), , 0.05, 1),
                "Must provide an argument for charge_range" )

  # Missing mtoz_diff
  expect_error( binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                      c(14, 45), c(2e3, 15e3), c(2, 10), , 1),
                "Must provide an argument for mtoz_diff" )

  # Missing time_diff  
  expect_error( binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                      c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, ),
                "Must provide an argument for time_diff" )
})




# ````````````````````` #
#  Wrong argument type  #
# ..................... #


test_that("binMS: wrong types of arguments", {

  # Wrong type for mass_spec
  expect_error(binMS(list(), "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "mass_spec must be either a matrix or data.frame")

  # Wrong type for mtoz
  expect_error(binMS(testMS, mtoz=list(), "chg", "mass", "time",
                     c("ms1", "ms2"), c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "mtoz must be either of mode numeric or character")

  # Wrong type for charge
  expect_error(binMS(testMS, "mtoz", list(), "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "charge must be either of mode numeric or character")

  # Wrong type for mass
  expect_error(binMS(testMS, "mtoz", "chg", list(), "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "mass must be either NULL or of mode numeric or character")

  # Wrong type for time_peak_reten
  expect_error(binMS(testMS, "mtoz", "chg", "mass", list(), c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "time_peak_reten must be either of mode numeric or character")

  # Wrong type for ms_inten
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", list(),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "ms_inten must be either NULL or of mode numeric or character")

  # Wrong type for time_range
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                     list(), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "time_range must be of mode numeric")

  # Wrong type for mass_range
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), list(), c(2, 10), 0.05, 1),
               "mass_range must be of mode numeric")

  # Wrong type for charge_range
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), list(), 0.05, 1),
               "charge_range must be of mode numeric")

  # Wrong type for mtoz_diff
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), list(), 1),
               "mtoz_diff must be of mode numeric")

  # Wrong type for time_diff
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, list()),
               "time_diff must be of mode numeric")
})




# ``````````````````````````````` #
#  Right type but illegal values  #
# ............................... #


test_that("binMS: arguments of the right type but illegal values", {

  # Provide name specifying the column for a column that doesn't exist
  expect_error(binMS(testMS, "asdf", "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "column names in mass_spec do not contain asdf element in mtoz")

  # What if there are no colnames
  
  # Provide index that is too small
  expect_error(binMS(testMS, 0, "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "out of bounds index 0 provided for mtoz relative to mass_spec")

  # Provide index that is too big
  expect_error(binMS(testMS, 999, "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "out of bounds index 999 provided for mtoz relative to mass_spec")

  # Provide empty numeric vector
  expect_error(binMS(testMS, numeric(0), "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "If non-NULL, then mtoz must have length no less than 1")

  # Provide empty character vector
  expect_error(binMS(testMS, character(0), "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "If non-NULL, then mtoz must have length no less than 1")

  # Provide numeric vector of wrong length
  expect_error(binMS(testMS, c(1, 2), "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "mtoz must have length 1 or length equal to the number of observations")

  # Provide character vector of wrong length
  expect_error(binMS(testMS, c("mtoz", "chg"), "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "mtoz must have length 1 or length equal to the number of observations")

  # Provide numeric vector with duplicates
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", c(5, 5L),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "ms_inten cannot have any duplicate values")

  # Provide character vector with duplicates
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms1"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "ms_inten cannot have any duplicate values")

  # Provide ambiguous character vector (i.e. multiple matches in colnames)
  expect_error(binMS(testMS, "m", "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "name provided had multiple matches in data - m element in mtoz")

  # Provide nonincreasing values for time_range
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                     c(45, 45L), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "The values of time_range must be in increasing order")

  # Provide missing numeric value for mtoz
  expect_error(binMS(testMS, numer_NA, "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "mtoz cannot contain any missing")

  # Provide missing character value for mtoz
  expect_error(binMS(testMS, char_NA, "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "mtoz cannot contain any missing")
  
  # Provide missing value for time_range
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, numer_NA), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "time_range cannot have any missing")

  # Provide missing value for mass_range
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, numer_NA), c(2, 10), 0.05, 1),
               "mass_range cannot have any missing")

  # Provide missing value for charge_range
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, numer_NA), 0.05, 1),
               "charge_range cannot have any missing")

  # Provide missing value for mtoz_diff
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), numer_NA, 1),
               "mtoz_diff cannot have any missing")

  # Provide missing value for time_diff
  expect_error(binMS(testMS, "mtoz", "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, numer_NA),
               "mtoz_diff cannot have any missing")

  # ********* CHeck if can pass nonexistent name of arg
  expect_error(binMS(testMS, asdf, "chg", "mass", "time", c("ms1", "ms2"),
                     c(14, 45), c(2e3, 15e3), c(2, 10), 0.05, 1),
               "")
})

f <- function(x, y) {
  for (var in c("x", "y")) {
    tryCatch(get(var), error=function(err) 1)
  }
}
