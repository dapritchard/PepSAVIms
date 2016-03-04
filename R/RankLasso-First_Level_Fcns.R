
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

  if (useAve && (nrow(bioMat) > 1L)) {
    # Want to average the replicates within a fraction.  Note that the input is
    # provided as a column for each fraction.
    return ( colMeans(bioMat[, bioIdx]) )
  }
  # case: either (i) we're not using the average values and we stack the
  # replicates into a single column, or (ii) we have a matrix with 1 row and
  # hence using c() strips the vector of it's dim attributes
  else {
    # Want to convert matrix to a non-list vector; uses c() to do this.  Note
    # that the behavior is to "stack" the columns.  Choice of doing this is
    # arbitrary, *but it must match the behavior of ms_conv()*.
    return ( c( bioMat[, bioIdx] ) )
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

  # calculate bioAve: the average for the bioactivity over the replicates
  if ( identical(nrow(bioMat), 1L) ) {
    bioAve <- c( bioMat[regionIdx$bio] )
  } else {
    bioAve <- colMeans( bioMat[, regionIdx$bio] )
  }

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



