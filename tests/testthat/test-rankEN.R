
# `````````````````````````````````` #
#  Load a saved dataset for testing  #
# .................................. #

# Load saved simulated data ----------------------------------------------------

load("tests/data/data-rankEN.RData")
# load("tests/Sim_Ms_Bio.RData")
#load("../Sim_Ms_Bio.RData")

# Create links for convenience for data from tests/Sim_Ms_Bio.RData
msDatObj <- testDat$msDat
bioact <- testDat$bioact

region_idx <- sim_args$regIdx
region_idx_matr <- matrix(region_idx, nrow=2)
region_idx_NA <- replace(region_idx, 1, NA)



biodf_NA_out_region <- biodf_NA_in_region <-
  biodf_char_out_region <- biodf_char_in_region <- data.frame(bioact)
biodf_NA_out_region[1, 1] <- NA
biodf_char_out_region$char_col <- rep("a", nrow(bioact))
biodf_NA_in_region[1, region_idx[1]] <- NA
biodf_char_in_region[, region_idx[1]] <- rep("a", nrow(bioact))

# Create a single numeric value with missing
na_num <- numeric(1)
na_num[1] <- NA

na_log <- logical(1)
na_log[1] <- NA

# Create data restriced to region of interest
msDat_region_only <- msDatObj[, region_idx]
bioact_region_only <- bioact[, region_idx]

# colnamesMS(testDat$msDat) <- paste0("mass_spec", 1:50)
# colnames(testDat$bioact) <- paste0("bioact", 1:50)



# `````````````````````````````````` #
#  Create a 'true' rankLasso object  #
# .................................. #

# Calculate model fits for comparison ------------------------------------------

# Explanetory data
ms_regr <- t( msDatObj$ms[, region_idx] )

# Outcome values
bio_regr <- colMeans(bioact[, region_idx])

# Fit model
lambda <- 0.1
enet_default <- elasticnet::enet(x=ms_regr, y=bio_regr, lambda=lambda)
# enet_all_cor <- elasticnet::enet(x=ms_regr, y=bio_regr, lambda=lambda, pos_only=FALSE)
# enet_keep_10 <- elasticnet::enet(x=ms_regr, y=bio_regr, lambda=lambda, ncomp=10)

# Obtain action index.  This is a vector of the compound indices as the enter /
# leave the model - a positive value means it entered, a negative value means it
# exited.
actions <- unlist(enet_default$actions)
actions <- actions[-length(actions)]

# Indices in the order that the compounds first entered the model
enter_idx <- unique( actions[actions > 0] )

# Obtain correlation values of proposed compounds
comp_cor <- apply(ms_regr, 2, function(x) cor(x, bio_regr))

# Indices in the order that the compounds first entered the model, after
# removing compounds that had a nonpositive correlation with the within-fraction
# average bioactivity levels
pos_idx <- enter_idx[comp_cor[enter_idx] > 0]
pos_10_idx <- pos_idx[ 1:min(length(pos_idx), 10L) ]
allcomp_10_idx <- enter_idx[ 1:min(length(pos_idx), 10L) ]

## Begin constructing rankEN object

# Dimensions of the data used for analysis
data_dim  = list(reg  = nrow(ms_regr),
                 comp = ncol(ms_regr),
                 repl = nrow(bioact))

# Column names for region of interest
region_nm <- list(ms  = colnames(msDatObj$ms)[region_idx],
                  bio = colnames(bioact)[region_idx])

# Create info for the summary function
summ_info_default <- list(
  data_dim  = data_dim,
  region_nm = region_nm,
  lambda    = 0.1,
  ncomp     = NULL,
  pos_only  = TRUE
)
summ_info_pos_10 <- summ_info_allcomp <- summ_info_allcomp_10 <- summ_info_default
summ_info_pos_10$ncomp <- summ_info_allcomp_10$ncomp <- 10L
summ_info_allcomp$pos_only <- summ_info_allcomp_10$pos_only <- FALSE

# Construct rankEN object using default settings
true_default <- list(mtoz      = msDatObj$mtoz[pos_idx],
                     charge    = msDatObj$chg[pos_idx],
                     comp_cor  = comp_cor[pos_idx],
                     enet_fit  = enet_default,
                     summ_info = summ_info_default)
class(true_default) <- "rankEN"

true_biovec_default <- true_default
true_biovec_default$summ_info$data_dim$repl <- 1L


true_keep_10 <- list(mtoz      = msDatObj$mtoz[pos_10_idx],
                     charge    = msDatObj$chg[pos_10_idx],
                     comp_cor  = comp_cor[pos_10_idx],
                     enet_fit  = enet_default,
                     summ_info = summ_info_pos_10)
true_allcomp <- list(mtoz      = msDatObj$mtoz[enter_idx],
                     charge    = msDatObj$chg[enter_idx],
                     comp_cor  = comp_cor[enter_idx],
                     enet_fit  = enet_default,
                     summ_info = summ_info_allcomp)
true_keep_10_allcomp <- list(mtoz      = msDatObj$mtoz[allcomp_10_idx],
                             charge    = msDatObj$chg[allcomp_10_idx],
                             comp_cor  = comp_cor[allcomp_10_idx],
                             enet_fit  = enet_default,
                             summ_info = summ_info_allcomp_10)
class(true_keep_10) <- class(true_allcomp) <- class(true_keep_10_allcomp) <- "rankEN"



# Create rankLasso object using function ---------------------------------------

# A few variations of pos_only and ncomp; matrix for bioact and indices for
# region specifiers
rankEN_default <- rankEN(msDatObj, bioact, region_idx, region_idx, lambda)
rankEN_keep_10 <- rankEN(msDatObj, bioact, region_idx, region_idx, lambda, , 10L)
rankEN_allcomp <- rankEN(msDatObj, bioact, region_idx, region_idx, lambda, FALSE)
rankEN_keep_10_allcomp <- rankEN(msDatObj, bioact, region_idx, region_idx, lambda, FALSE, 10L)

# bioact as a data.frame or vector
rankEN_biodf_default <- rankEN(msDatObj, data.frame(bioact), region_idx, region_idx, lambda)
rankEN_biovec_default <- rankEN(msDatObj, colMeans(bioact), region_idx, region_idx, lambda)

# region specifiers as character vectors
rankEN_ms_reg_char_default <-
  rankEN(msDatObj, bioact, paste0("ms", 21:30), region_idx, lambda)
rankEN_bio_reg_char_default <-
  rankEN(msDatObj, bioact, region_idx, paste0("bio", 21:30), lambda)
rankEN_biovec_reg_char_default <-
  rankEN(msDatObj, colMeans(bioact), region_idx, paste0("bio", 21:30), lambda)

# region specifiers as NULL
rankEN_ms_region_null_allcomp <-
  rankEN(msDat_region_only, bioact, NULL, region_idx, lambda, FALSE)
rankEN_bio_region_null_keep10 <-
  rankEN(msDatObj, bioact_region_only, region_idx, NULL, lambda, , 10L)
rankEN_biodf_region_null_default <-
  rankEN(msDatObj, data.frame(bioact_region_only), region_idx, NULL, lambda)
rankEN_biovec_region_null_default <-
  rankEN(msDatObj, colMeans(bioact_region_only), region_idx, NULL, lambda)

# bioact as data.frame / vector + character region specifier
rankEN_biodf_region_char_default <-
  rankEN(msDatObj, data.frame(bioact), region_idx, paste0("bio", 21:30), lambda)
rankEN_biovec_region_char_default <-
  rankEN(msDatObj, colMeans(bioact), region_idx, paste0("bio", 21:30), lambda)

# bioact with missing or character outside region of interest
rankEN_NA_outside_default <-
  rankEN(msDatObj, biodf_NA_out_region, region_idx, region_idx, lambda)
rankEN_char_outside_default <-
  rankEN(msDatObj, biodf_char_out_region, region_idx, region_idx, lambda)

# region as a matrix
rankEN_ms_region_matr_default <-
  rankEN(msDatObj, bioact, region_idx_matr, region_idx, lambda)
rankEN_biovec_region_matr_default <-
  rankEN(msDatObj, colMeans(bioact), region_idx, region_idx_matr, lambda)



# ``````````````````` #
#    Begin testing    #
# ................... #

context("rankEN method")

# Test list elements one-at-a-time for easier location of differences
test_that("compare rankEN list elements one-by-one", {
  expect_identical( rankEN_default$mtoz,      true_default$mtoz )
  expect_identical( rankEN_default$charge,    true_default$charge )
  expect_identical( rankEN_default$comp_cor,  true_default$comp_cor )
  expect_identical( rankEN_default$enet_fit,  true_default$enet_fit )
  expect_identical( rankEN_default$summ_info, true_default$summ_info )
})

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








