
# Simulate data for testing purposes
#
# Very, very, naive / unrealistic function used for simulating data.  It creates
# an msDat object and bioactivity matrix based on a given region and beta vector
# of regressors.  The data is unrealistic in that the mass spectrometry data is
# sampled from N(0, 1) distributions independent across both compounds and
# fractions, and hence does not have the structure characteristic of mass
# spectrometry intensities.


simData <- function(nCmp, nFrac, nRepl, nPred, regIdx, sigma) {

  # Mass spec data
  ms <- matrix(rnorm(nCmp * nFrac), nrow=nCmp, ncol=nFrac)
  colnames(ms) <- paste0("ms", seq_len(nFrac))

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
  colnames(bio) <- paste0("bio", seq_len(nFrac))

  # Construct msDat object with arbitrary mtoz and charge values
  msDat_out <- msDat(ms, seq_len(nCmp), rep(1, nCmp))

  out <- list( msDat   = msDat_out,
               bioact  = bio,
               region  = regIdx,
               beta    = beta )

  return (out)
}




