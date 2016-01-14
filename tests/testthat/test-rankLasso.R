
context("rankLasso method")

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
# save(testDat, sim_args, file="tests/Sim_Ms_Bio_v2.RData")


load("tests/Sim_Ms_Bio_v2.RData")

msDat <- testDat$msDat
bioact <- testDat$bioact
regIdx <- sim_args$regIdx

# Explanetory variables
ms <- t( msDat$ms[, regIdx] )

# Outcome vars
bio <- colMeans(bioact[, regIdx])

lars_fit <- lars::lars(x=ms, y=bio)

actions <- unlist(lars_fit$actions)
cmpIdx <- unique( actions[actions > 0] )

true_out <- list( mtoz     = msDat$mtoz[cmpIdx],
                  charge   = msDat$chg[cmpIdx],
                  nRegions = length(regIdx),
                  nRepl    = sim_args$nRepl,
                  useAve   = TRUE,
                  lars_fit = lars_fit )
class(true_out) <- "rankCmp"


# Some testing here ------------------------------------------------------------

rl_out <- rankLasso(msDat, bioact, regIdx, TRUE)

expect_identical(true_out$mtoz, rl_out$mtoz)
expect_identical(true_out$charge, rl_out$charge)
expect_identical(true_out$nRegions, rl_out$nRegions)
expect_identical(4L, rl_out$nRepl)
expect_identical(true_out$useAve, rl_out$useAve)
expect_identical(true_out$lars_fit, rl_out$lars_fit)












