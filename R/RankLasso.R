
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
#'   and the bioactivity data refer to the same fractions.  In other words if
#'   the \code{k}-th column is selected, then the \code{k}-th column must refer
#'   to the same fraction for both the mass spectrometry data and the
#'   bioactivity data.
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


rankLasso <- function(msDat, bioact, region=NULL) {

  # Set useAve to true since it is cheaper and result is invariant to choice
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

  structure(outDat, class="rankCmp")
}




#' Calculates indices for the fractions of interest
#'
#' See description of the \code{region} parameter for the behavior of
#' this function
#'
#' @inheritParams rankLasso
#'
#' @return A list providing the indices for the fractions of interest
#'   for the mass spectrometry and bioactivity data.  The list contains
#'   the elements described below.  Note that the \code{k}-th index
#'   for the mass spectrometry index corresponds to the \code{k}-th
#'   index for the bioactivity index.
#'
#'   \describe{
#'
#'   \item{\code{ms}}{An integer vector providing the indices for the
#'   fractions of interest in the mass spectrometry data }
#'
#'   \item{\code{bio}}{An integer vector providing the indices for the
#'   fractions of interest in the bioactivity data }
#'
#'   }


getRegionIdx <- function(msDat, bioact, region) {

  # case: null
  if (is.null(region)) {
    regionIdx <- null_to_idx(msDat, bioact)
  }

  # case: list
  else {
    regionIdx <- list( ms  = reg_to_idx(msDat, bioact, region$ms, "ms"),
                       bio = reg_to_idx(msDat, bioact, region$bio, "bio") )
  }

  return (regionIdx)
}




#' Convert mass spectrometry data
#'
#' Creates the explanetory data matrix for use by the \code{lars} function. This
#' involves transposing the data and reducing the data to the fractions of
#' interest.  Furthermore, if there are replicates for the bioactivity data,
#' then multiple copies of each fraction of interest are made, one for each
#' replicate.
#'
#' It is important to define how the ordering of the replicates is performed so
#' that this behavior matches that of \code{\link{conv_bioact}}.  The current
#' implementation places replicates within a fraction consecutively, so that for
#' example, the data may look like (going down the rows): fraction 1, fraction
#' 1, fraction 1, fraction 2, fraction 2, fraction 2, fraction 3, ...
#'
#' @inheritParams rankLasso
#'
#' @param ms A \code{matrix} containing mass spectrometry intensity readings.
#'   Each column provides the mass spectrometry values for a given fraction, and
#'   each row provides the mass spectrometry values for a given mass-to-charge
#'   ratio value across the fractions.
#'
#' @param msIdx An integer vector providing the indices (and hence fractions)
#'   for the mass spectrometry data which are to be used in the Lasso model
#'
#' @param nRepl An integer giving the number of replications provided in the
#'   bioactivity data.  Note that this paramter is not needed when \code{useAve}
#'   is \code{TRUE}.
#'
#' @return A matrix providing the mass spectrometry (and hence the explanetory)
#'   data for the Lasso model.  The rows are the observations (possibly expanded
#'   from the original data) for the model, and the columns are the compounds
#'   (and hence are the predictor variables).


conv_ms <- function(ms, msIdx, nRepl, useAve) {

  # If there are no bioactivity replicates or we are averaging the replicates,
  # then this is all that needs to be done to the data
  ms <- t( ms[, msIdx] )

  # If we have bioactivity replicates, then we need to copy each fraction of the
  # data (now the rows after transposition), one time for each replicate
  if (!useAve && (nRepl >= 2L)) {
    ms <- ms[rep(seq_len(nrow(ms)), each=nRepl), ]
  }

  return (ms)
}




#' Convert bioactivity data
#'
#' Creates the response vector for use by the \code{lars} function. This
#' involves transposing the data and reducing the data to the fractions of
#' interest.  If replicates are to be used, then the fractions containing the
#' replicates are concatenated together.  If the replicates are to be averaged,
#' then the data within a fraction is averaged before concatenating the values.
#'
#' It is important to define how the ordering of the replicates is performed so
#' that this behavior matches that of \code{\link{conv_ms}}.  The current
#' implementation places replicates within a fraction consecutively, so that for
#' example, the data may look like: fraction 1 replicate 1, fraction 1 replicate
#' 2, fraction 1 replicate 3, fraction 2 replicate 1, fraction 2 replicate 2,
#' fraction 2 replicate 3, fraction 3 replicate 1, ...
#'
#' @param bioMat The object returned by the format_bio function
#'
#' @param bioIdx An integer vector providing the indices (and hence fractions of
#'   interest) for the bioactivity data which are to be used in the Lasso model


conv_bioact <- function(bioMat, bioIdx, useAve) {

  if (useAve) {
    # Want to average the replicates within a fraction.  Note that the input is
    # provided as a column for each fraction.
    return ( colMeans(bioact[, bioIdx]) )
  }
  else {
    # Want to convert matrix to a non-list vector; uses c() to do this.  Note
    # that the behavior is to "stack" the columns.  Choice of doing this is
    # arbitrary, *but it must match the behavior of ms_conv()*.
    return ( c( bioact[, bioIdx] ) )
  }
}




#' Extract compound ordering from \code{lars} fit
#'
#' Returns an integer vector indexing (some of the) compounds by the order in
#' which they first enter the Lasso model, from first to last.
#'
#' @param lars_fit An object of class \code{lars}
#'
#' @details   A \code{lars} object is a list containing (among other things) an
#'   obect '\code{actions}'. This object is a list of all the actions taken in
#'   the model Lasso path, where actions are taken to be either adding or
#'   removing a predictor variable from model.  If the variable is added to the
#'   model, then the column index is inserted into the next position in the
#'   list.  If the variable is removed from the model, then -1 times the column
#'   index of the variable is inserted into the next position in the list.
#'
#'   \code{getCmpIdx} works by reducing '\code{actions}' list to the first
#'   positive time that a number is in the list.
#'
#' @return An integer vector indexing (some of the) compounds by the order in
#'   which they first enter the Lasso model, from first to last.


getCmpIdx <- function(lars_fit) {

  # See 'details' in function documentation to see why this function works

  actions <- unlist(lars_fit$actions)
  unique(actions[actions > 0])
}




getCmpCor <- function(msDat, bioMat, regionIdx, cmpIdx) {

  cmpMat <- msDat$ms[cmpIdx, regionIdx$ms]
  bioAve <- colMeans( bioMat[, regionIdx$bio] )

  apply(cmpMat, 1, function(x) cor(x, bioAve))
}




get_data_desc <- function(ms, bioMat, regionIdx, cmpIdx) {

  msDim <- dim(ms)
  bioDim <- dim(bioMat)

  regionNm <- list( ms  = colnames(ms)[regionIdx$ms],
                    bio = colnames(bioMat)[regionIdx$bio] )

  data_desc <- list( msDim     = msDim,
                     bioDim    = bioDim,
                     regionNm  = regionNm,
                     regionIdx = regionIdx,
                     cmpIdx    = cmpIdx )
  return (data_desc)
}



