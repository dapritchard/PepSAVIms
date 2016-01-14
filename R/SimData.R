#
#
# nCmp <- 200                       # Number of compounds
# nFrac <- 50                       # Number of fractions
# nRepl <- 4                        # Number of bioactivity replicates
# regIdx <- 21:25                   # Indices for the nonzero predictors
# nPred <- 3       # Number of nonzero predictors

nCmp=3
nFrac=50
nRepl=1
regIdx=11:40
nPred=3

simData <- function(nCmp, nFrac, nRepl, nPred, regIdx, sigma) {

  # Mass spec data
  ms <- matrix(rnorm(nCmp * nFrac), nrow=nCmp, ncol=nFrac)

  # Define true beta (just the nonzero part of beta)
  beta <- list( idx = seq_len(nPred),
                val = seq(1, 1 / nPred, length.out=nPred) )

  # True mean for response, across the fractions (fractions with with
  # association that is, as specified by regIdx).  Arbitrarily chooses the first
  # nPred compounds to be the significant compounds.
  repl_mean <- drop( crossprod( ms[beta$idx, regIdx], beta$val ) )

  # Sample the response.  sapply, and hence in turn replicate, writes each
  # replicate as a column to a matrix.  Thus t yields a matrix with rows as
  # replicates and columns as fractions
  nfrac_assoc <- length(regIdx)
  bio_resp <- t( replicate(nRepl, repl_mean + rnorm(nfrac_assoc, sd=sigma)) )

  # Column indexes for the non-associative fractions before and after the
  # associative region
  non_assoc_idx <- list( bef = seq(1, head(regIdx, 1) - 1),
                         aft = seq(tail(regIdx, 1) + 1, nFrac) )

  # Sample bioactivity matrix
  bio <- matrix( c( rnorm(nRepl * length(non_assoc_idx$bef)),
                    bio_resp,
                    rnorm(nRepl * length(non_assoc_idx$aft) ) ),
                 nrow=nRepl,
                 ncol=nFrac )

  # Construct msDat object with arbitrary mtoz and charge values
  msDat_out <- msDat(ms, seq_len(nCmp), rep(1, nCmp))

  out <- list( msDat   = msDat_out,
               bioact  = bio,
               region  = regIdx,
               beta    = beta )

  return (out)
}




