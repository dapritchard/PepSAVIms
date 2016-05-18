

# ``````````````````` #
#  Load testing data  #
# ................... #

# See construct-data/cd-filterMS.R for the file used to create the data
load("../data/data-rankEN_test.RData")




# `````````````````````````````````````````````````` #
#  Valid input with various argument specifications  #
# .................................................. #

context("rankEN method")

test_that("filterMS: msObj with valid input", {

  # msObj: filterMS object
  out <- rankEN(filterMS_obj, bioact, reg_idx, reg_idx, lambda)
  expect_identical(true_rankEN$pos, out)

  # bioact: matrix
  out <- rankEN(msDatObj, bioact, reg_idx, reg_idx, lambda)
  expect_identical(true_rankEN$pos, out)

  # bioact: data.frame
  out <- rankEN(msDatObj, bio_df, reg_idx, reg_idx, lambda)
  expect_identical(true_rankEN$pos, out)

  # bioact: vector
  out <- rankEN(msDatObj, bio_vec, reg_idx, reg_idx, lambda)
  expect_identical(true_rankEN_bio_ave, out)

  # bioact: (1 x n) matrix
  out <- rankEN(msDatObj, bio_matr_ave, reg_idx, reg_idx, lambda)
  expect_identical(true_rankEN_bio_ave, out)

  # region_ms, region_bio: character, numeric
  out <- rankEN(msDatObj, bioact, paste0(reg_idx), reg_idx, lambda)
  expect_identical(true_rankEN$pos, out)

  # region_ms, region_bio: character, characterr
  out <- rankEN(msDatObj, bioact, paste0(reg_idx), paste0(reg_idx), lambda)
  expect_identical(true_rankEN$pos, out)

  # region_ms, region_bio: NULL, numeric
  out <- rankEN(ms_reg_only, bioact, NULL, reg_idx, lambda)
  expect_identical(true_rankEN$pos, out)

  # region_ms, region_bio: numeric, NULL
  out <- rankEN(msDatObj, bio_reg_only, reg_idx, NULL, lambda)
  expect_identical(true_rankEN$pos, out)

  # region_ms, region_bio: numeric, numeric (with bioact a vector)
  out <- rankEN(msDatObj, bio_vec_reg_only, reg_idx, NULL, lambda)
  expect_identical(true_rankEN_bio_ave, out)

  # pos_only, ncomp: TRUE, 10
  out <- rankEN(msDatObj, bioact, reg_idx, reg_idx, lambda, TRUE, 10L)
  expect_identical(true_rankEN$pos_10, out)

  # pos_only, ncomp: FALSE, NULL
  out <- rankEN(msDatObj, bioact, reg_idx, reg_idx, lambda, FALSE, NULL)
  expect_identical(true_rankEN$all, out)

  # pos_only, ncomp: FALSE, 10
  out <- rankEN(msDatObj, bioact, reg_idx, reg_idx, lambda, FALSE, 10L)
  expect_identical(true_rankEN$all_10, out)
})



                 

# # A few variations of pos_only and ncomp; matrix for bioact and indices for
# # region specifiers
# test_that("rankEN: various choices of ncomp, pos_only", {
#   expect_identical( rankEN_default,         true_default )
#   expect_identical( rankEN_keep_10,         true_keep_10 )
#   expect_identical( rankEN_allcomp,         true_allcomp )
#   expect_identical( rankEN_keep_10_allcomp, true_keep_10_allcomp )
# })

# # bioact as a data.frame or vector
# # TODO: NA outside region

# test_that("rankEN: bioact as a data.frame or vector", {
#   expect_identical( rankEN_biodf_default,  true_default )
#   expect_identical( rankEN_biovec_default, true_biovec_default )
# })

# # region specifiers as character vectors
# test_that("rankEN: region specifiers as character vectors", {
#   expect_identical( rankEN_ms_reg_char_default,     true_default )
#   expect_identical( rankEN_bio_reg_char_default,    true_default )
#   expect_identical( rankEN_biovec_reg_char_default, true_biovec_default )
# })

# # region specifiers as NULL
# test_that("rankEN: region specifiers as NULL", {
#   expect_identical( rankEN_ms_region_null_allcomp,     true_allcomp )
#   expect_identical( rankEN_bio_region_null_keep10,     true_keep_10 )
#   expect_identical( rankEN_biodf_region_null_default,  true_default )
#   expect_identical( rankEN_biovec_region_null_default, true_biovec_default )
# })

# # bioact as data.frame / vector + character region specifier
# test_that("rankEN: bioact as data.frame / vector + character region specifier", {
#   expect_identical( rankEN_biodf_region_char_default,  true_default )
#   expect_identical( rankEN_biovec_region_char_default, true_biovec_default )
# })

# # bioact with NA or character data outside of region of interest
# test_that("rankEN: bioact as data.frame / vector + character region specifier", {
#   expect_identical( rankEN_char_outside_default,  true_default )
#   expect_identical( rankEN_char_outside_default, true_default )
# })

# # region as a matrix
# test_that("rankEN: region as a matrix", {
#   expect_identical( rankEN_ms_region_matr_default,     true_default )
#   expect_identical( rankEN_biovec_region_matr_default, true_biovec_default )
# })




# ``````````````````````````` #
#  Missing param for formals  #
# ........................... #

test_that("rankEN: missing input", {
  
  # Note that the remaining formal args have defaults
  
  # msObj missing
  expect_error(rankEN( , bioact, reg_idx, reg_idx, lambda),
               "Must provide an argument for msObj" )
  
  # bioact missing
  expect_error(rankEN(msDatObj, , reg_idx, reg_idx, lambda),
               "Must provide an argument for bioact" )
  
  # lambda missing
  expect_error(rankEN(msDatObj, bioact, reg_idx, reg_idx, ),
               "Must provide an argument for lambda" )
})




# ```````````````````````````````` #
#  Argument is nonexistent object  #
# ................................ #

test_that("rankEN: nonexistent object", {
  
  # msObj arg a nonexistent object
  expect_error(rankEN(asdf, bioact, reg_idx, reg_idx, lambda),
               "object 'asdf' not found for msObj")

  # bioact arg a nonexistent object
  expect_error(rankEN(filterMS_obj, asdf, reg_idx, reg_idx, lambda),
               "object 'asdf' not found for bioact")

  # region_ms arg a nonexistent object
  expect_error(rankEN(filterMS_obj, bioact, asdf, reg_idx, lambda),
               "object 'asdf' not found for region_ms")

  # region_bio arg a nonexistent object
  expect_error(rankEN(filterMS_obj, bioact, reg_idx, asdf, lambda),
               "object 'asdf' not found for region_bio")

  # lambda arg a nonexistent object
  expect_error(rankEN(filterMS_obj, bioact, reg_idx, reg_idx, asdf),
               "object 'asdf' not found for lambda")

  # pos_only arg a nonexistent object
  expect_error(rankEN(filterMS_obj, bioact, reg_idx, reg_idx, lambda, asdf),
               "object 'asdf' not found for pos_only")

  # ncomp arg a nonexistent object
  expect_error(rankEN(filterMS_obj, bioact, reg_idx, reg_idx, lambda, , asdf),
               "object 'asdf' not found for ncomp")
})




# ``````````````````````````````` #
#  Argument is of the wrong type  #
# ............................... #




## Invalid input

# rankEN <- function(msObj, bioact, region_ms=NULL, region_bio=NULL, lambda,
#                    pos_only=TRUE, ncomp=NULL)

# Check for missing arguments
expect_error( rankEN(, bioact, reg_idx, reg_idx, lambda),
              "Must provide an argument for msObj")
expect_error( rankEN(msDatObj, , reg_idx, reg_idx, lambda),
              "Must provide an argument for bioact")
expect_error( rankEN(msDatObj, bioact, reg_idx, reg_idx, ),
              "Must provide an argument for lambda")

expect_error( rankEN(list(), bioact, reg_idx, reg_idx, lambda),
              "msObj must be an object of class \"msDat\"")
expect_error( rankEN(msDatObj, list(), reg_idx, reg_idx, lambda),
              "bioact must be either a data.frame or of mode numeric")
expect_error( rankEN(msDatObj, bioact, list(), reg_idx, lambda),
              "region_ms must be either NULL or either of mode numeric or character")
expect_error( rankEN(msDatObj, bioact, reg_idx, logical(10), lambda),
              "region_bio must be either NULL or either of mode numeric or character")
expect_error( rankEN(msDatObj, bioact, reg_idx, reg_idx, list()),
              "lambda must be a numeric value")
expect_error( rankEN(msDatObj, bioact, reg_idx, reg_idx, -11),
              "lambda cannot be smaller than 0")
expect_error( rankEN(msDatObj, bioact, reg_idx, reg_idx, 1:3),
              "lambda must be an atomic vector of length 1")
expect_error( rankEN(msDatObj, bioact, reg_idx, reg_idx, 1, "asdf"),
              "pos_only must be either TRUE or FALSE")
expect_error( rankEN(msDatObj, bioact, reg_idx, reg_idx, 1, , list()),
              "ncomp must be either NULL or a numeric value")
expect_error( rankEN(msDatObj, bioact, reg_idx, reg_idx, 1, , 0.5),
              "If non-NULL then ncomp must be >= 1")
expect_error( rankEN(msDatObj, bioact, reg_idx, reg_idx, 1, , 2:4),
              "If non-NULL then ncomp must be an atomic vector of length 1")

expect_error( rankEN(msDatObj, biodf_char_in_region, reg_idx, reg_idx, lambda),
              "bioact cannot have any non-numeric in the region provided by region_bio" )
expect_error( rankEN(msDatObj, biodf_NA_in_region, reg_idx, reg_idx, lambda),
              "bioact cannot have any missing in the region provided by region_bio" )
expect_error( rankEN(msDatObj, bioact, reg_idx_NA, reg_idx, lambda),
              "region_ms cannot contain any missing" )
expect_error( rankEN(msDatObj, bioact, reg_idx, reg_idx_NA, lambda),
              "region_bio cannot contain any missing" )
expect_error( rankEN(msDatObj, bioact, reg_idx, reg_idx, na_num),
              "lambda cannot contain any missing" )
expect_error( rankEN(msDatObj, bioact, reg_idx, reg_idx, lambda, na_log),
              "pos_only must be either TRUE or FALSE" )
expect_error( rankEN(msDatObj, bioact, reg_idx, reg_idx, lambda, , na_num),
              "If non-NULL then ncomp cannot contain any missing" )
expect_error( rankEN(msDatObj, bioact, reg_idx, region_idx[-1L], lambda),
              paste0("Number of fractions for mass spectrometry (10) does not ",
                     "match the number of fractions for bioactivity (9)"), fixed=TRUE )





