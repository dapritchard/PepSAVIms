
reg_to_idx <- function(msDat, bioact, regVec, whichReg) {

  # case: numeric region provided
  if (is.numeric(regVec)) {
    check_num_reg(msDat, bioact, regVec, whichReg)
    return (regVec)
  }

  # case: character region provided
  else {
    check_char_reg(msDat, bioact, regVec, whichReg)
    char_to_idx(1, 2)
  }
}




check_num_reg <- function(msDat, bioact, regVec, whichReg) {

  # Get applicable largest index value
  if (identical(whichReg, "ms")) {
    max_region_colnum <- ncol(msDat$ms)
  }
  else if (identical(whichReg, "bio")) {
    max_region_colnum <- ifelse(is.strictVec(bioact), length(bioact), ncol(bioact))
  }

#   # case: "both"
#   else {
#     max_region_colnum <- max( msDat$ms,
#                               ifelse(is.strictVec(bioact),
#                                      length(bioact),
#                                      ncol(bioact)) )
#   }

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


