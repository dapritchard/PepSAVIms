
checkValInp_rankLasso <- function(msDat, bioact, region, useAve) {

  check_form_msDat(msDat)
  check_form_bioact(bioact)
  check_form_region(bioact, region)
  check_form_useAve(useAve)
}




check_form_msDat <- function(msDat) {

  if ( missing(msDat) ) {
    stop("msDat must have an argument provided\n")
  }
  if (class(msDat) != "msDat") {
    stop("msDat must be of class msDat\n")
  }
}




check_form_bioact <- function(bioact) {

  if ( missing(bioact) ) {
    stop("bioact must have an argument provided\n")
  }

  # Check that bioact either a non-array atomic vector or a matrix or a data frame
  if ( !( is.strictvec(bioact)
          || is.matrix(bioact)
          || is.data.frame(bioact) ) ) {
    stop("bioact must be either a non-array atomic vector or a matrix or a data frame\n")
  }
}




# Check if region is the right form.  Note that this doesn't yet check that the values
# provided by region make sense (i.e. the numbers provided are not too big / small
# or the names don't match).  These checks are performed in reg_to_idx().

check_form_region <- function(bioact, region) {

  # case: no argument provided
  # noop; a default is provided

  # case: region is a non-array atomic vector
  if ( is.strictvec(region) ) {
    if ( !(is.numeric(region) || is.character(region)) ) {
      stop("region must be either numeric or character\n")
    }
    if ( TRUE %in% is.na(region) ) {
      stop("region cannot have any missing values\n")
    }
  }

  # case: region is a list
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

  # case: region is NULL
  else if ( is.null(region) ) {
    # noop
  }

  # case: region neither NULL or a non-array atomic vector or a list
  else {
    stop("region must either be NULL or a non-array atomic vector or a list\n")
  }

}




check_form_useAve <- function(useAve) {

  if ( missing(useAve) ) {
    stop("useAve must have an argument provided\n")
  }
  if ( !(identical(useAve, TRUE) || identical(useAve, FALSE)) ) {
    stop("useAve must either be TRUE or FALSE\n")
  }
}
