



# ``````````````````` #
#    Begin testing    #
# ................... #

context("rankEN method")

# # Test list elements one-at-a-time for easier location of differences
# test_that("compare rankEN list elements one-by-one", {
#   expect_identical( rankEN_default$mtoz,      true_default$mtoz )
#   expect_identical( rankEN_default$charge,    true_default$charge )
#   expect_identical( rankEN_default$comp_cor,  true_default$comp_cor )
#   expect_identical( rankEN_default$enet_fit,  true_default$enet_fit )
#   expect_identical( rankEN_default$summ_info, true_default$summ_info )
# })

bio_df <- data.frame(bioact)
bio_df[, "bio1"] <- "asdf"

bio_vec <- colMeans(bioact)
bio_matr_ave <- matrix(bio_vec, nrow=1, dimnames=list(NULL, colnames(bioact)))

true_rankEN_bio_ave <- true_rankEN$pos
true_rankEN_bio_ave$summ_info$data_dim$repl <- 1L






# bioact args

out <- rankEN(msDatObj, bioact, region_idx, region_idx, lambda)
expect_identical(true_rankEN$pos, out)

out <- rankEN(msDatObj, bio_df, region_idx, region_idx, lambda)
expect_identical(true_rankEN$pos, out)

out <- rankEN(msDatObj, bio_vec, region_idx, region_idx, lambda)
expect_identical(true_rankEN_bio_ave, out)

out <- rankEN(msDatObj, bio_matr_ave, region_idx, region_idx, lambda)
expect_identical(true_rankEN_bio_ave, out)

# region args






                 

# A few variations of pos_only and ncomp; matrix for bioact and indices for
# region specifiers
test_that("rankEN: various choices of ncomp, pos_only", {
  expect_identical( rankEN_default,         true_default )
  expect_identical( rankEN_keep_10,         true_keep_10 )
  expect_identical( rankEN_allcomp,         true_allcomp )
  expect_identical( rankEN_keep_10_allcomp, true_keep_10_allcomp )
})

# bioact as a data.frame or vector
# TODO: NA outside region

test_that("rankEN: bioact as a data.frame or vector", {
  expect_identical( rankEN_biodf_default,  true_default )
  expect_identical( rankEN_biovec_default, true_biovec_default )
})

# region specifiers as character vectors
test_that("rankEN: region specifiers as character vectors", {
  expect_identical( rankEN_ms_reg_char_default,     true_default )
  expect_identical( rankEN_bio_reg_char_default,    true_default )
  expect_identical( rankEN_biovec_reg_char_default, true_biovec_default )
})

# region specifiers as NULL
test_that("rankEN: region specifiers as NULL", {
  expect_identical( rankEN_ms_region_null_allcomp,     true_allcomp )
  expect_identical( rankEN_bio_region_null_keep10,     true_keep_10 )
  expect_identical( rankEN_biodf_region_null_default,  true_default )
  expect_identical( rankEN_biovec_region_null_default, true_biovec_default )
})

# bioact as data.frame / vector + character region specifier
test_that("rankEN: bioact as data.frame / vector + character region specifier", {
  expect_identical( rankEN_biodf_region_char_default,  true_default )
  expect_identical( rankEN_biovec_region_char_default, true_biovec_default )
})

# bioact with NA or character data outside of region of interest
test_that("rankEN: bioact as data.frame / vector + character region specifier", {
  expect_identical( rankEN_char_outside_default,  true_default )
  expect_identical( rankEN_char_outside_default, true_default )
})

# region as a matrix
test_that("rankEN: region as a matrix", {
  expect_identical( rankEN_ms_region_matr_default,     true_default )
  expect_identical( rankEN_biovec_region_matr_default, true_biovec_default )
})





## Invalid input

# rankEN <- function(msObj, bioact, region_ms=NULL, region_bio=NULL, lambda,
#                    pos_only=TRUE, ncomp=NULL)

# Check for missing arguments
expect_error( rankEN(, bioact, region_idx, region_idx, lambda),
              "Must provide an argument for msObj")
expect_error( rankEN(msDatObj, , region_idx, region_idx, lambda),
              "Must provide an argument for bioact")
expect_error( rankEN(msDatObj, bioact, region_idx, region_idx, ),
              "Must provide an argument for lambda")

expect_error( rankEN(list(), bioact, region_idx, region_idx, lambda),
              "msObj must be an object of class \"msDat\"")
expect_error( rankEN(msDatObj, list(), region_idx, region_idx, lambda),
              "bioact must be either a data.frame or of mode numeric")
expect_error( rankEN(msDatObj, bioact, list(), region_idx, lambda),
              "region_ms must be either NULL or either of mode numeric or character")
expect_error( rankEN(msDatObj, bioact, region_idx, logical(10), lambda),
              "region_bio must be either NULL or either of mode numeric or character")
expect_error( rankEN(msDatObj, bioact, region_idx, region_idx, list()),
              "lambda must be a numeric value")
expect_error( rankEN(msDatObj, bioact, region_idx, region_idx, -11),
              "lambda cannot be smaller than 0")
expect_error( rankEN(msDatObj, bioact, region_idx, region_idx, 1:3),
              "lambda must be an atomic vector of length 1")
expect_error( rankEN(msDatObj, bioact, region_idx, region_idx, 1, "asdf"),
              "pos_only must be either TRUE or FALSE")
expect_error( rankEN(msDatObj, bioact, region_idx, region_idx, 1, , list()),
              "ncomp must be either NULL or a numeric value")
expect_error( rankEN(msDatObj, bioact, region_idx, region_idx, 1, , 0.5),
              "If non-NULL then ncomp must be >= 1")
expect_error( rankEN(msDatObj, bioact, region_idx, region_idx, 1, , 2:4),
              "If non-NULL then ncomp must be an atomic vector of length 1")

expect_error( rankEN(msDatObj, biodf_char_in_region, region_idx, region_idx, lambda),
              "bioact cannot have any non-numeric in the region provided by region_bio" )
expect_error( rankEN(msDatObj, biodf_NA_in_region, region_idx, region_idx, lambda),
              "bioact cannot have any missing in the region provided by region_bio" )
expect_error( rankEN(msDatObj, bioact, region_idx_NA, region_idx, lambda),
              "region_ms cannot contain any missing" )
expect_error( rankEN(msDatObj, bioact, region_idx, region_idx_NA, lambda),
              "region_bio cannot contain any missing" )
expect_error( rankEN(msDatObj, bioact, region_idx, region_idx, na_num),
              "lambda cannot contain any missing" )
expect_error( rankEN(msDatObj, bioact, region_idx, region_idx, lambda, na_log),
              "pos_only must be either TRUE or FALSE" )
expect_error( rankEN(msDatObj, bioact, region_idx, region_idx, lambda, , na_num),
              "If non-NULL then ncomp cannot contain any missing" )
expect_error( rankEN(msDatObj, bioact, region_idx, region_idx[-1L], lambda),
              paste0("Number of fractions for mass spectrometry (10) does not ",
                     "match the number of fractions for bioactivity (9)"), fixed=TRUE )


# Assume that region_ms and region_bio will provide appropriate errors when they
# are passed the right type



# checkwhich <- function(curr, targ) {
#   cat(sep="",
#       identical(curr$mtoz, targ$mtoz), "\n",
#       identical(curr$charge, targ$charge), "\n",
#       identical(curr$comp_cor, targ$comp_cor), "\n",
#       identical(curr$enet_fit, targ$enet_fit), "\n",
#       identical(curr$summ_info, targ$summ_info), "\n")
# }

# checkwhich(rankEN_allcomp, true_allcomp)








