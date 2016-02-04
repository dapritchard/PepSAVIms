
#' Ranks compounds using the Lasso path
#'
#' Returns identifying information for the compounds in the order in which they
#' first enter the Lasso model
#'
#' @param msDat An object of class \code{\link{msDat}} containing the mass
#'   spectrometry data and identifying information
#'
#' @param bioact Either a non-list vector, a matrix, or a data frame providing
#'   bioactivity data.  If a non-list vector, then it is assumed that each entry
#'   corresponds to a particular fraction.  If the data is 2-dimensional, then
#'   it is assumed that each column corresponds to a particular fraction, and
#'   that each row corresponds to a particular bioactivity replicate.
#'
#' @param region Either \code{NULL}, a non-list vector, a matrix, or a list
#'   containing exactly two atomic vectors, providing information specifying
#'   which fractions are to be included in the Lasso model.  Note that a data
#'   frame satisfies these requirements.
#'
#'   If \code{NULL}, then it is assumed that all fractions included in the data
#'   are to be used in the model.  This requires that the number of fractions in
#'   the data for the parameter passed to \code{bioact} be the same as the
#'   number of fractions in the mass spectrometry data for the paramter passed
#'   to \code{msDat}.  It further assumes that the \code{k}-th column must refer
#'   to the same fraction for both the mass spectrometry data and the
#'   bioactivity data for every \code{k}.
#'
#'   If a non-list vector then it must be either a numeric or charcter vector,
#'   such that the vector specifies which columns (and hence which fractions) to
#'   include in the model.  If the vector is numeric, then the desired columns
#'   are specified by number, and if the vector is character, then the desired
#'   columns are specified by name (partial matching is allowed).  Note that
#'   this assumes that the corresponding columns in the mass spectrometry data
#'   and the bioactivity data refer to the same fractions.
#'
#'   If a matrix then it must be either numeric or character with exactly two
#'   columns - one is to be named \code{ms} and the other is to be named
#'   \code{bio}.  The \code{ms} column specifies the columns (and hence
#'   fractions) to include in the model from the mass spectrometry data, either
#'   as a vector of the column numbers or the column names.  The \code{bio}
#'   vector specifies the columns (and hence fractions) to include in the model
#'   from the bioactivity data, either as a vector of the column numbers or the
#'   column names.  It is assumed that two entries from a given row refer to the
#'   same fraction.
#'
#'   If a list, then it must be a list with two named non-list vectors of equal
#'   length - one is to be named \code{ms} and the other is to be named
#'   \code{bio} (note that a \code{n x 2} data frame satisfies this
#'   requirement).  The \code{ms} vector specifies the columns (and hence
#'   fractions) to include in the model from the mass spectrometry data, either
#'   as a vector of the column numbers or the column names.  The \code{bio}
#'   vector specifies the columns (and hence fractions) to include in the model
#'   from the bioactivity data, either as a vector of the column numbers or the
#'   column names.  It is assumed that the column from the mass spectrometry
#'   data specified by the \code{k}-th value in the \code{ms} vector corresponds
#'   to the same fraction as the column specified by the \code{k}-th value in
#'   the \code{bio} vector, for each \code{k}.
#'
#' @details Note that in the current incarnation of rankLasso, the solution to
#'   the Lasso path is the same when using either the average of bioactivity
#'   replicates or individual replicates.  A parameter, useAve, used to be
#'   offered - but since the result is the same either way, now it is just set
#'   \code{TRUE}.  The rest of the code is left unchanged, in the event another
#'   way is found to use individual replicates.
#'
#' @export

# TODO: document function output


rankLasso <- function(msDat, bioact, region=NULL) {

  # Set useAve to true since it is cheaper and result is invariant to choice.
  # See details in documentation for more information.
  useAve <- TRUE

  # Check that data arguments are of the right type
  checkValInp_rankLasso(msDat, bioact, region, useAve)

  # Convert (if necessary) bioact to matrix form.  Data is guaranteed to be
  # numeric with no missing.
  bioMat <- format_bio(bioact)

  # Convert (if necessary) region to either a list containing exactly two
  # vectors, or NULL.  If applicable, each vector is guaranteed to be numeric or
  # character with no missing.
  regList <- format_reg(region)

  # Creates region indices for the mass spec and bioactivity data based on the
  # region input.  Also checks to see if the input matches the data.
  regionIdx <- getRegionIdx(msDat, bioMat, regList)

  # Puts the data into a form where the rows are the fractions of interest and
  # columns are the compounds.  This is done so as to cast the problem as a
  # regression problem where the fractions are the observations and the
  # compounds are the potential predictor variables.  In addition, if useAve is
  # FALSE then replicates of the observations are created a la an ANOVA setting.
  ms_regr <- conv_ms(msDat$ms, regionIdx$ms, nrow(bioact), useAve)

  # Creates a numeric vector for the bioactivity response
  bio_regr <- conv_bioact(bioMat, regionIdx$bio, useAve)

  # Calculate the Lasso path
  lars_fit <- lars::lars(x=ms_regr, y=bio_regr)

  # Obtain indices for compounds as they first enter the model
  cmpIdx <- getCmpIdx(lars_fit)

  # Correlation of chosen compounds
  cmp_cor <- getCmpCor(msDat, bioMat, regionIdx, cmpIdx)

  # Record some summary statistics for the data for use by bioact.summary
  data_desc <- get_data_desc(msDat$ms, bioMat, regionIdx, cmpIdx)

  # Construct output object
  outDat <- list( mtoz      = msDat$mtoz[cmpIdx],
                  charge    = msDat$chg[cmpIdx],
                  cmp_cor   = cmp_cor,
                  data_desc = data_desc,
                  lars_fit  = lars_fit )

  structure(outDat, class="rankLasso")
}




summary.rankLasso <- function(rl_out) {

  # Create a link for convenience
  data_desc <- rl_out$data_desc

  # Print mass spectrometry and bioactivity data dimensions
  cat(sep="",
      "\n",
      "Mass spectrometry data:\n",
      "-----------------------\n",
      "  ", format(data_desc$msDim[1], width=6, big.mark=","), " compounds\n",
      "  ", format(data_desc$msDim[2], width=6, big.mark=","), " fractions\n",
      "\n",
      "Bioactivity data:\n",
      "-----------------\n",
      "  ", format(data_desc$bioDim[1], width=6, big.mark=","), " replicates\n",
      "  ", format(data_desc$bioDim[2], width=6, big.mark=","), " fractions\n",
      "\n")


  # TODO: have to consider case where (data_desc$regionNm$ms == NULL) and similalry for bio

  blkCol <- rep("    ", length(data_desc$regionIdx$ms))
  regDat <- setNames( data.frame( data_desc$regionNm$ms,
                                  data_desc$regionIdx$ms,
                                  blkCol,
                                  data_desc$regionNm$bio,
                                  data_desc$regionIdx$bio ),
                      c("M.S. Names", "  M.S. Index", "", "Bio. Names", "  Bio. Index") )

  # Print a table with the names and indices for mass spectrometry and bioactivity data
  cat(sep="",
      "Fractions included in region of interest:\n",
      "-----------------------------------------\n")
  print(regDat, row.names=FALSE)
  cat("\n")

  # Create a table with the mass-to-charge, charge, and correlation values for
  # chosen compounds
  sigCmp <- setNames( data.frame(rl_out$mtoz,
                                 rl_out$charge,
                                 rl_out$cmp_cor),
                      c("Mass-to-charge", "  Charge", "  Correlation") )
  sig_dig <- nchar( length( data_desc$cmpIdx ) )

  # Print a table with the mass-to-charge, charge, and correlation values for
  # chosen compounds
  cat(sep="",
      "Compounds chosen by the model (", length(data_desc$cmpIdx), " total):\n",
      "-------------------------------",        rep("-", sig_dig), "--------\n")
  print(sigCmp, row.names=FALSE, digits=4)
  cat("\n")

}





