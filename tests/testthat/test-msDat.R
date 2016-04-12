
# `````````````````````````````` #
#  Create data used for testing  #
# .............................. #


# Create raw data to use as msDat arguments ------------------------------------

# This is the raw data in several forms, which is to be coerced into class
# msDat.  Columns 1 and 2 of msMat are the mass-to-charge ratio and charge
# values, respectively.

msMat <- matrix(1:25, nrow=5, dimnames=list(paste0(1:5), c("mtoz", "chg", "ms1", "ms2", "ms3")))
msDf <- data.frame(msMat)
mtoz <- msMat[, 1]
chg <- msMat[, 2]


# Manually construct "true" msDat object ---------------------------------------

# What the data above should be coerced to by msDat function

ms_true <- matrix(11:25, nrow=5, dimnames=list(paste0(1:5), c("ms1", "ms2", "ms3")))
trueDatOut <- structure(list(ms   = ms_true,
                             mtoz = setNames(1:5, paste(1:5)),
                             chg  = setNames(6:10, paste(1:5))),
                        class="msDat")

trueOneCol <- structure(list(ms   = ms_true[, "ms1", drop=FALSE],
                             mtoz = setNames(1:5, paste(1:5)),
                             chg  = setNames(6:10, paste(1:5))),
                        class="msDat" )




# ``````````````````` #
#    Begin testing    #
# ................... #

context("msDat constructor")


# Testing valid inputs ---------------------------------------------------------


# Check that function works as intended when supplied a matrix for msDat

test_that("msDat: use matrix as input", {
  
  # Check various combinations of supplying columns, column names, indices, NULL
  expect_identical( msDat(msMat, 1, 2),                              trueDatOut )
  expect_identical( msDat(msMat[, -1], mtoz, 1),                     trueDatOut )
  expect_identical( msDat(msMat, mtoz, 2, c("ms1", "ms2", "ms3")),   trueDatOut )
  expect_identical( msDat(msMat[, -2], 1, chg),                      trueDatOut )
  expect_identical( msDat(msMat, 1, chg, c("ms1", "ms2", "ms3")),    trueDatOut )
  expect_identical( msDat(msMat[, -c(1, 2)], mtoz, chg),             trueDatOut )
  expect_identical( msDat(msMat, mtoz, chg, c("ms1", "ms2", "ms3")), trueDatOut )
  expect_identical( msDat(msMat, mtoz, chg, 3:5),                    trueDatOut )
  expect_identical( msDat(msMat, "mtoz", 2),                         trueDatOut )
  expect_identical( msDat(msMat, 1, "chg"),                          trueDatOut )
  expect_identical( msDat(msMat, "mtoz", "chg"),                     trueDatOut )
  
  # Intensity matrix has 1 dimension (should still be a matrix in msDat obj)
  expect_identical( msDat(msMat, "mtoz", "chg", "ms1"), trueOneCol )
})




# Check that function works as intended when supplied a matrix for msDat

test_that("msDat: use data frame as input", {
  
  # Check various combinations of supplying columns, column names, indices, NULL
  expect_identical( msDat(msDf, 1, 2),                              trueDatOut )
  expect_identical( msDat(msDf[, -1], mtoz, 1),                     trueDatOut )
  expect_identical( msDat(msDf, mtoz, 2, c("ms1", "ms2", "ms3")),   trueDatOut )
  expect_identical( msDat(msDf[, -2], 1, chg),                      trueDatOut )
  expect_identical( msDat(msDf, 1, chg, c("ms1", "ms2", "ms3")),    trueDatOut )
  expect_identical( msDat(msDf[, -c(1, 2)], mtoz, chg),             trueDatOut )
  expect_identical( msDat(msDf, mtoz, chg, c("ms1", "ms2", "ms3")), trueDatOut )
  expect_identical( msDat(msDf, mtoz, chg, 3:5),                    trueDatOut )
  expect_identical( msDat(msDf, "mtoz", 2),                         trueDatOut )
  expect_identical( msDat(msDf, 1, "chg"),                          trueDatOut )
  expect_identical( msDat(msDf, "mtoz", "chg"),                     trueDatOut )
  
  # Intensity matrix has 1 dimension (should still be a matrix in msDat obj)
  expect_identical( msDat(msDf, "mtoz", "chg", "ms1"), trueOneCol )
})




# Check that function stops when input is missing or is of the wrong type

test_that("msDat: missing or wrong type of input", {
  
  # Check various combinations of supplying columns, column names, indices, NULL
  expect_error( msDat(mtoz=mtoz, charge=chg), "Must provide an argument for mass_spec" )
  expect_error( msDat(mass_spec=msMat, charge=1), "Must provide an argument for mtoz" )
  expect_error( msDat(mass_spec=msMat, mtoz=1), "Must provide an argument for charge" )

  # Check if wrong type of input supplied
  expect_error(msDat(mass_spec=list(), "mtoz", "chg"),
               "mass_spec must be either a matrix or data.frame")
  expect_error(msDat(mass_spec=list(), "mtoz", "chg"),
               "mass_spec must be either a matrix or data.frame")
  expect_error(msDat(msMat, mtoz=logical(nrow(msMat)), "chg"),
               "mtoz must be either of mode numeric or character")
  expect_error(msDat(msMat, "mtoz", charge=logical(nrow(msMat))),
               "charge must be either of mode numeric or character")
  expect_error(msDat(msMat, "mtoz", "chg", ms_inten=list()),
               "ms_inten must be either NULL or of mode numeric or character")
})




# Provide input of the right type but values that don't make sense

test_that("msDat: input that doesn't make sense", {

  expect_error(msDat(ms_true, 1:100, chg),
               "mtoz must have length 1 or length equal to the number of observations")
  expect_error(msDat(ms_true, integer(0), chg),
               "mtoz must have length 1 or length equal to the number of observations")
  expect_error(msDat(ms_true, mtoz, 1:100),
               "charge must have length 1 or length equal to the number of observations")
  expect_error(msDat(ms_true, mtoz, integer(0)),
               "charge must have length 1 or length equal to the number of observations")
  expect_error(msDat(ms_true[, integer(0)], "mtoz", "chg"),
               "mass_spec cannot have 0 columns")
  expect_error(msDat(ms_true[, integer(0)], mtoz, chg),
               "mass_spec cannot have 0 columns")
  expect_error(msDat(msMat[integer(0), ], "mtoz", "chg"),
               "mass_spec must have 2 or more rows")
  expect_error(msDat(msDf[integer(1), ], "mtoz", "chg"),
               "mass_spec must have 2 or more rows")
  expect_error(msDat(ms_true[, integer(0)], mtoz, chg),
               "mass_spec cannot have 0 columns")
  expect_error(msDat(msMat, 0, 2),
               "out of bounds value 0 provided for mtoz")
  expect_error(msDat(msMat, 1000, 2),
               "out of bounds value 1000 provided for mtoz")
  expect_error(msDat(msMat, 1, -100),
               "out of bounds value -100 provided for charge")
  expect_error(msDat(msMat, 1, 2, 3:6),
               "out of bounds value 6 provided for ms_inten")
  expect_error(msDat(msMat, "nomatch", "chg"),
               "name provided not in data - nomatch element in mtoz")
  expect_error(msDat(msMat, "", "chg"),
               "name provided had multiple matches in data -  element in mtoz")
  expect_error(msDat(msMat, "mtoz", "chg", c("ms1", "ms2", "ms1")),
               "ms_inten cannot have any duplicate values")
  expect_error(msDat(msMat, "mtoz", "chg", character(0)),
               "If non-NULL, then ms_inten must have length >= 1")
  expect_error(msDat(msMat, character(0), "chg"),
               "mtoz must have length 1 or length equal to the number of observations")
  expect_error(msDat(msMat[, c("mtoz", "chg")], "mtoz", "chg"),
               "There cannot be 0 columns left for ms_inten after removing data for other variables")
})
