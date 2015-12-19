
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
#' @param region Either \code{NULL}, a non-list vector, or a list containing
#'   exactly two atomic vectors, providing information specifying which
#'   fractions are to be included in the Lasso model.
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
#'   If a list, then it must be a list with two named non-list vectors of equal
#'   length - one is to be named \code{ms} and the other is to be named
#'   \code{bio}.  The \code{ms} vector specifies the columns (and hence
#'   fractions) to include in the model from the mass spectrometry data, either
#'   as a vector of the column numbers or the column names.  The \code{bio}
#'   vector specifies the columns (and hence fractions) to include in the model
#'   from the bioactivity data, either as a vector of the column numbers or the
#'   column names.  It is assumed that the column from the mass spectrometry
#'   data specified by the \code{k}-th value in the \code{ms} vector corresponds
#'   to the same fraction as the column specified by the \code{k}-th value in
#'   the \code{bio} vector, for each \code{k}.
#'
#' @param useAve A logical value specifying whether or not to average replicate
#'   bioactivity observations.  Ignored if only one bioactivity observation is
#'   provided.


rankLasso <- function(msDat, bioact, region=NULL, useAve=TRUE) {

  # Checks the form of the input.  Note that some checks for the correctness of
  # the input are performed later.
  checkValInp_rankLasso(msDat, bioact, region, useAve)

  # Creates region indices for the mass spec and bioactivity data based on the
  # region input.  Also checks to see if the input matches the data.
  regionIdx <- getRegionIdx(msDat, bioact, region)

  # Puts the data into a form where the rows are the fractions and cols are the
  # compounds.  This is done so as to cast the problem as a regression problem
  # where the fractions are the observations and the compounds are the potential
  # predictor variables.  In addition, if useAve is FALSE then replicates of the
  # observations are created a la an ANOVA setting.
  ms <- conv_ms(msDat$ms, regionIdx$ms, ncol(bioact), useAve)

  # Creates a non-list vector for the bioactivity response
  bio <- conv_bioact(bioact, regionIdx$bio, useAve)

  # Clean up unneeded data for garbage collection
  mtoz <- msDat$mtoz
  chg <- msDat$chg
  rm(msDat, bioact)

  # Calculate the Lasso path
  fit <- lars(x=ms, y=bio)

  # Indices for compounds as they first enter the model
  cmpIdx <- getCmpIdx(fit)

  # Construct output object
  outDat <- list( mtoz     = mtoz[cmpIdx],
                  charge   = chg[cmpIdx],
                  nRegions = length(regionIdx$ms),
                  nRepl    = 1,
                  useAve   = useAve )

  structure(outDat, class="rankCmp")
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
#' example, the data may look like (going down the rows): fraction 1,
#' fraction 1, fraction 1, fraction 2, fraction 2, fraction 2, fraction 3, ...
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
  if (!useAve) {
    ms <- ms[rep(1:nrow(ms), each=nRepl), ]
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
#' @inheritParams rankLasso
#'
#' @param msIdx An integer vector providing the indices (and hence fractions)
#'   for the bioactivity data which are to be used in the Lasso model


conv_bioact <- function(bioact, bioIdx, useAve) {

  # case: bioact is non-list vector
  if (is.strictVec(bioact)) {
    # useAve ignored in this case as there are no replicates.   Simply return
    # one bioactivity observation for each fraction of interest.
    return ( bioact[bioIdx] )
  }

  # case: bioact is matrix or data frame
  bioact <- as.matrix(bioact)

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
#' @return An integer vector indexing (some of the) compounds by the order in
#'   which they first enter the Lasso model, from first to last.


getCmpIdx <- function(lars_fit) {

  # A lars object is a list containing (among other things) an obect 'actions'.
  # This object is a list of all the actions taken in the model Lasso path,
  # where actions are taken to be either adding or removing a predictor variable
  # from model.  If the variable is added to the model, then the column index is
  # inserted into the next position in the list.  If the variable is removed
  # from the model, then -1 times the column index of the variable is inserted
  # into the next position in the list.

  actions <- unlist(lars_fit$actions)
  unique(actions[actions > 0])
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

    # Check that ms and bioact dimensions match
    nfrac <- ifelse(is.vector(bioact), length(bioact), ncol(bioact))
    if ( !identical(ncol(msDat$ms), nfrac) ) {
      stop("If region is NULL then the number of fractions must be the same for
           the mass spectrometry data and the bioactivity data\n")
    }

    # Simply make an index entry for every fraction
    regionIdx <- list( ms  = seq_len(nfrac),
                       bio = seq_len(nfrac) )
  }

  # case: non-list vector
  else if (is.strictVec(region)) {

    # case: numeric vector
    if (is.numeric(region)) {
      check_numeric_region(msDat, bioact, region, "both")
      regionIdx <- list( ms  = as.integer(region),
                         bio = as.integer(region) )
    }
    # case: character vector
    else {
      check_character_region(msDat, bioact, region, "both")
      regionIdx <- list( ms  = char_to_idx(colnames(msDat$ms), region),
                         bio = char_to_idx(colnames(bioact), region) )
    }
  }

  # case: list
  else {
    # case: numeric region$ms
    if (is.numeric(region$ms)) {
      check_numeric_region(msDat, bioact, region$ms, "ms")
      regionIdx <- list(ms=region)
    }
    # case: character region$ms
    else {
      check_character_region(msDat, bioact, region, "ms")
      regionIdx <- list(ms=char_to_idx(colnames(msDat$ms), region$ms))
    }

    # case: numeric region$bio
    if (is.numeric(region$bio)) {
      check_numeric_region(msDat, bioact, region$bio, "bio")
      regionIdx$bio <- region$bio
    }
    # case: character region$bio
    else {
      check_character_region(msDat, bioact, region$bio, "bio")
      regionIdx$bio <- char_to_idx(colnames(bioact), region$bio)
    }
  }

  return (regionIdx)
}




conv_regionIdx <- function(msDat, bioact, regvec, whichReg) {

  # case: numeric region provided
  if (is.numeric(regvec)) {
    check_numeric_region(msDat, bioact, regvec, whichReg)
    return (regvec)
  }

  # case: character region provided
  else {
    check_character_region(msDat, bioact, regvec, whichReg)
    char_to_idx(1, 2)
  }
}




check_numeric_region <- function(msDat, bioact, regionIdx, whichRegion) {

  # Get applicable largest index value
  if (identical(whichRegion, "ms")) {
    max_region_colnum <- ncol(msDat$ms)
  }
  else if (identical(whichRegion, "bio")) {
    max_region_colnum <- ifelse(is.strictVec(bioact), length(bioact), ncol(bioact))
  }
  # case: "both"
  else {
    max_region_colnum <- max( msDat$ms,
                              ifelse(is.strictVec(bioact),
                                     length(bioact),
                                     ncol(bioact)) )
  }

  # Check valid input values
  if ((min(regionIdx) < 1) || (max(regionIdx) > max_region_colnum)) {
    stop("out of bounds region value provided\n")
  }
  else if (length(unique(as.integer(regionIdx))) < length(regionIdx)) {
    stop("no duplicate values allowed in region input\n")
  }
}




check_character_region <- function(msDat, bioact, regionNames, whichRegion) {

  if (length(unique(regionNames)) < length(regionNames)) {
    stop("no duplicate values allowed in region input\n")
  }


  if (identical(whichRegion, "ms")) {
    if (FALSE %in% (regionNames %in% colnames(msDat$ms))) {
      stop("names in provided region not in data\n")
    }
  }
  else if (identical(whichRegion, "bio")) {
    if (FALSE %in% (regionNames %in% colnames(bioact))) {
      stop("names in provided region not in data\n")
    }
  }
  # case: both
  else {
    if ( FALSE %in% (regionNames %in% colnames(msDat$ms))
         || FALSE %in% (regionNames %in% colnames(bioact)) ) {
      stop("names in provided region not in data\n")
    }
  }
}




char_to_idx <- function(datNames, charVec) {
  sapply(charVec, function(x) which(charVec == datNames))
}



checkValInp_rankLasso <- function(msDat, bioact, region, useAve) {

  # Check if msDat of class msDat
  if (class(msDat) != "msDat") {
    stop("msDat must be of class msDat\n")
  }

  # Check that bioact is the right form
  if ( !(is.strictVec(drop(bioact)) || is.matrix(drop(bioact)) || is.data.frame(bioact)) ) {
    stop("bioact must be either a non-list vector or a matrix or a data frame\n")
  }
  else if ( !is.numeric( as.matrix(bioact) ) ) {
    stop("bioact must be populated with only numeric values\n")
  }
  else if ( TRUE %in% is.na(bioact) ) {
    stop("bioact cannot have any missing values\n")
  }

  # Check if region is the right form.  Note that this doesn't yet check that the values
  # provided by region make sense (i.e. the numbers provided are not too big / small
  # or the names don't match).

  if ( !(is.null(region) || is.vector(region)) ){
    stop("region must either be NULL or a non-list vector or a list\n")
  }
  # Non-list vector case
  else if ( is.strictVec(region) ) {
    if ( !(is.numeric(region) || is.character(region)) ) {
      stop("region must be either numeric or character\n")
    }
    if ( TRUE %in% is.na(region) ) {
      stop("region cannot have any missing values\n")
    }
  }
  # List case
  else {
    if ( !( identical(names(region), c("ms", "bio"))
                 || identical(names(region), c("bio", "ms")) ) ) {
      stop("if region is a list then the objects must have the names ms and bio\n")
    }
    else if ( !(is.strictVec(bioact$ms) && is.strictVec(bioact$bio)) ) {
      stop("if region is a list then it must contain exactly two non-list vectors\n")
    }
    else if ( length(bioact$ms) != length(bioact$bio) ) {
      stop("if region is a list then it must contain
           exactly two non-list vectors of the same length\n")
    }
    else if ( (TRUE %in% is.na(bioact$ms))
              || (TRUE %in% is.na(bioact$bio)) ) {
      stop("region cannot have any missing values\n")
    }
  }
}




