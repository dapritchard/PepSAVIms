
is.strictvec <- function(x) {
  is.atomic(x) && !is.array(x)
}



check_form_msDat <- function(msDat) {

  if ( missing(msDat) ) {
    stop("msDat must have an argument provided\n")
  }
  if (class(msDat) != "msDat") {
    stop("msDat must be of class msDat\n")
  }
}




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

  # Make an index entry for every fraction
  regionIdx <- list( ms  = seq_len(nfrac),
                     bio = seq_len(nfrac) )
  return (regionIdx)
}




num_to_idx <- function(msDat, bioact, regVec, whichDat) {

  nfrac <- ifelse(identical(whichDat, "ms"), ncol(msDat$ms), ncol(bioact))

  # Check valid input values
  if ((min(regVec) < 1) || (max(regVec) > nfrac)) {
    stop("out of bounds region value provided\n")
  }

  return (regVec)
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
    matchBool <- grepl(regVec[i], nmFrac)
    nMatch <- sum(matchBool)

    # Check that current name has exactly one match
    if (identical(nMatch, 0L)) {
      stop("names in provided region not in data\n")
    }
    if (nMatch >= 2L) {
      stop("names in provided region had multiple matches in data\n")
    }

    regionIdx[i] <- which(matchBool)
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




checkValInp_rankLasso <- function(msDat, bioact, region, useAve) {

  # Check if msDat of class msDat
  if (class(msDat) != "msDat") {
    stop("msDat must be of class msDat\n")
  }

  # Check if region is the right form.  Note that this doesn't yet check that the values
  # provided by region make sense (i.e. the numbers provided are not too big / small
  # or the names don't match).  These checks are performed in reg_to_idx().

  if ( !(is.null(region) || is.strictvec(region) || is.list(region)) ) {
    stop("region must either be NULL or a non-array atomic vector or a list\n")
  }

  # case: non-array atomic vector
  if ( is.strictvec(region) ) {
    if ( !(is.numeric(region) || is.character(region)) ) {
      stop("region must be either numeric or character\n")
    }
    if ( TRUE %in% is.na(region) ) {
      stop("region cannot have any missing values\n")
    }
  }

  # case: list
  else if ( is.list(region) ) {

    # Check that (exactly two objects) with names ms and bio contained in region
    if ( !( identical(names(region), c("ms", "bio"))
            || identical(names(region), c("bio", "ms")) ) ) {
      stop("if region is a list then the objects must have the names ms and bio\n")
    }

    # Check that objects are non-array atomic vectors
    if ( !(is.strictvec(region$ms) || is.strictvec(region$bio)) ){
      stop("if region is a list then the objects must be non-array atomic vectors\n")
    }

    # Check that objects are numeric or character
    if ( !(is.numeric(bioact$ms) || is.character(bioact$ms))
         || !(is.numeric(bioact$ms) || is.character(bioact$ms)) ) {
      stop("if region is a list then the objects must numeric or character\n")
    }
  } # end list case

  # Check that useAve is the right form
  if ( !(identical(useAve, TRUE) || identical(useAve, FALSE)) ) {
    stop("useAve must either be TRUE or FALSE\n")
  }
}


