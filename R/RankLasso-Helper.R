
# Delegates generatation of region indices to num_to_idx or char_to_idx,
# dependent on type of data.  Duplicates are also checked for.

reg_to_idx <- function(msDat, bioact, regVec, whichDat) {

  # Check that there are no duplicates
  if ( length(unique(regVec)) != length(regVec) ) {
    stop("region cannot have any duplicate values\n")
  }

  # case: numeric region provided
  if (is.numeric(regVec)) {
    return ( num_to_idx(msDat, bioact, regVec, whichDat) )
  }
  # case: character region provided
  else {
    return ( char_to_idx(msDat, bioact, regVec, whichDat) )
  }
}




null_to_idx <- function(msDat, bioact) {

  # Check that ms and bioact dimensions match
  if ( !identical(ncol(msDat$ms), ncol(bioact)) ) {
    stop("If region is NULL then the number of fractions must be the same for
           the mass spectrometry data and the bioactivity data\n")
  }
  nfrac <- ncol(msDat$ms)

  # Make an index entry for every fraction
  regionIdx <- list( ms  = seq_len(nfrac),
                     bio = seq_len(nfrac) )
  return (regionIdx)
}




num_to_idx <- function(msDat, bioact, regVec, whichDat) {

  # whichDat is either ms or bio, specifying which data the regions are
  # referring to
  nfrac <- ifelse(identical(whichDat, "ms"), ncol(msDat$ms), ncol(bioact))

  # Check valid input values
  if ((min(regVec) < 1) || (max(regVec) > nfrac)) {
    stop("out of bounds region value provided\n")
  }

  return ( as.integer(regVec) )
}




char_to_idx <- function(msDat, bioact, regVec, whichReg) {

  # nmFrac: a vector of the fraction names
  nmFrac <- get_data_nm(msDat, bioact, whichReg)

  # regionIdx: container for the indexes corresponding to provided regVec
  regionIdx <- integer( length(regVec) )

  # Loop iterates over provided names for region and checks each one to see
  # that it has exactly one match
  for (i in seq_along(regVec)) {

    # Number of matches for current element of regVec in fraction names
    matchIdx <- grep(regVec[i], nmFrac)
    nMatch <- length(matchIdx)

    # Check that current name has exactly one match
    if (identical(nMatch, 0L)) {
      stop("names in provided region not in data\n")
    }
    if (nMatch >= 2L) {
      stop("names in provided region had multiple matches in data\n")
    }

    regionIdx[i] <- matchIdx
  }

  return (regionIdx)
}




get_data_nm <- function(msDat, bioact, whichDat) {

  # case: regVec provided for mass spec data
  if (identical(whichDat, "ms")) {
    nmFrac <- colnames(msDat$ms)
  }
  # case: regVec provided for bioactivity data
  else {
    nmFrac <- colnames(bioact)
  }

  return (nmFrac)
}




is.strictvec <- function(x) {
  is.atomic(x) && !is.array(x)
}
