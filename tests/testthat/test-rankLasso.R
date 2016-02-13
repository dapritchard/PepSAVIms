
# `````````````````````````````````` #
#  Load a saved dataset for testing  #
# .................................. #

# Load saved simulated data ----------------------------------------------------

# The following are the commands used to simulate the data used for testing

# sim_args <- list(
#   nCmp   = 200L,   # Number of compounds
#   nFrac  = 50L,    # Number of fractions
#   nRepl  = 4L,     # Number of bioactivity replicates
#   regIdx = 21:30,  # Indices for the nonzero predictors
#   nPred  = 3L,     # Number of nonzero predictors
#   sigma  = 0.1     # Additive variance for bioactivity data
# )
#
# testDat <- with(sim_args, simData(nCmp, nFrac, nRepl, nPred, regIdx, sigma))
#
# save(testDat, sim_args, file="tests/Sim_Ms_Bio.RData")
#
# load("tests/Sim_Ms_Bio.RData")

load("../Sim_Ms_Bio.RData")




# `````````````````````````````````` #
#  Create a 'true' rankLasso object  #
# .................................. #

# Calculate model fits for comparison ------------------------------------------

# Rename for convenience.  testDat and sim_args are from Sim_Ms_Bio.RData.
msDat <- testDat$msDat
bioact <- testDat$bioact
regIdx <- sim_args$regIdx

# Explanetory data
ms_regr <- t( msDat$ms[, regIdx] )

# Outcome values
bio_regr <- colMeans(bioact[, regIdx])

# Fit model
lars_fit <- lars::lars(x=ms_regr, y=bio_regr)

# Obtain compound indices
actions <- unlist(lars_fit$actions)
cmpIdx <- unique( actions[actions > 0] )

# Obtain correlation values of proposed compounds
cmp_cor <- apply(ms_regr[, cmpIdx], 2, function(x) cor(x, bio_regr))

# Column names for region of interest
regionNm <- list( ms  = colnames(msDat$ms)[regIdx],
                  bio = colnames(bioact)[regIdx] )

# Indices for region of interest
regionIdx <- list( ms  = 21:30,
                   bio = 21:30 )

# Data statistics included in rankLasso output
data_desc <- list( msDim     = c(200L, 50L),
                   bioDim    = c(4L, 50L),
                   regionNm  = regionNm,
                   regionIdx = regionIdx,
                   cmpIdx    = cmpIdx )

# Object to compare to rankLasso output
true_out <- list( mtoz      = msDat$mtoz[cmpIdx],
                  charge    = msDat$chg[cmpIdx],
                  cmp_cor   = cmp_cor,
                  data_desc = data_desc,
                  lars_fit  = lars_fit )
class(true_out) <- "rankLasso"


# Create rankLasso object using function ---------------------------------------

rl_out <- rankLasso(msDat, bioact, regIdx)




# ``````````````````` #
#    Begin testing    #
# ................... #

context("rankLasso method")

test_that("use matrix as input", {
  expect_identical( rl_out$mtoz,      true_out$mtoz )
  expect_identical( rl_out$charge,    true_out$charge )
  expect_identical( rl_out$cmp_cor,   true_out$cmp_cor )
  expect_identical( rl_out$data_desc, true_out$data_desc )
  expect_identical( rl_out$lars_fit,  true_out$lars_fit )
  expect_identical( rl_out,           true_out )
})










