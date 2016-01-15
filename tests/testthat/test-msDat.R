
context("msDat constructor")

# Create raw data to use as msDat arguments ------------------------------------

# This is the raw data in several forms, which is to be coerced into class
# msDat.  Columns 1 and 2 of msMat are the mass-to-charge ratio and charge
# values, respectively.

# Colnames chosen to match those assi
msMat <- matrix(1:25, nrow=5)
msDf <- data.frame(msMat)
colnames(msMat) <- colnames(msDf)
mtoz <- msMat[, 1]
chg <- msMat[, 2]


# Manually construct "true" msDat object ---------------------------------------

# What the data above should be coerced to by msDat function

ms_true <- matrix(11:25, nrow=5)
colnames(ms_true) <- paste0("X", 3:nrow(msMat))
trueDatOut <- structure( list( ms   = ms_true,
                               mtoz = 1:5,
                               chg  = 6:10 ),
                         class="msDat" )


# Testing valid inputs ---------------------------------------------------------

test_that("use matrix as input", {
  expect_identical( msDat(msMat, 1, 2),                  trueDatOut )
  expect_identical( msDat(msMat[, -1], mtoz, 1),         trueDatOut )
  expect_identical( msDat(msMat[, -2], 1, chg),          trueDatOut )
  expect_identical( msDat(msMat[, -c(1, 2)], mtoz, chg), trueDatOut )
})

test_that("use data frame as input", {
  expect_identical( msDat(msDf, 1, 2),                  trueDatOut )
  expect_identical( msDat(msDf[, -1], mtoz, 1),         trueDatOut )
  expect_identical( msDat(msDf[, -2], 1, chg),          trueDatOut )
  expect_identical( msDat(msDf[, -c(1, 2)], mtoz, chg), trueDatOut )
  expect_identical( msDat(msDf, "X1", 2),               trueDatOut )
  expect_identical( msDat(msDf, 1, "X2"),               trueDatOut )
  expect_identical( msDat(msDf, "X1", "X2"),            trueDatOut )
})


# Test invalid input -----------------------------------------------------------

# TODO



