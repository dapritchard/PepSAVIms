
format_bio <- function(bioact) {

  # Delete the dimensions of an array which have only one level
  bioact <- drop(bioact)

  # Check that bioact either a non-array atomic vector or a matrix or a data
  # frame.  If so then convert to a matrix.
  # case: vector
  if ( is.strictvec(bioact) ) {
    bioMat <- matrix(bioact, nrow=1)
  }
  # case: matrix
  else if ( is.matrix(bioact) ) {
    # noop; this is the desired form
  }
  # case: data frame
  else if ( is.data.frame(bioact) ) {
    bioact <- as.matrix(bioact)
  }
  # case: invalid input
  else {
    stop("bioact must be either a non-array atomic vector or a matrix or a data frame\n")
  }

  # Check that bioact is numeric w/o missing
  if ( !is.numeric(bioact) ) {
    stop("bioact must be populated with only numeric values\n")
  }
  if ( any( is.na(bioact) ) ) {
    stop("region cannot have missing values\n")
  }

  return (bioact)
}




format_reg <- function(region) {

  # Delete the dimensions of an array which have only one level
  region <- drop(region)

  # Call specialized function depending on type of input
  if (is.null(region)) {
    return ( NULL )
  }
  if (is.strictvec(region)) {
    return ( format_reg_vec(region) )
  }
  if (is.matrix(region)) {
    return ( format_reg_mat(region) )
  }
  if (is.list(region)) {
    return ( format_reg_list(region) )
  }

  # If we've made it here, then region doesn't satisfy function req's
  stop("region must be either a non-array vector, a matrix, or a list\n")
}




format_reg_vec <- function(region) {

  if ( !(is.numeric(region) || is.character(region)) ) {
    stop("region must be either numeric or character\n")
  }
  if ( any( is.na(region) ) ) {
    stop("region cannot have missing values\n")
  }

  regList <- list( ms  = region,
                   bio = region )
  return (regList)
}




format_reg_mat <- function() {

  if (ncol(region) != 2) {
    stop("if a matrix then region must have exactly two columns\n")
  }
  if ( !( identical(colnames(region), c("ms", "bio"))
          || identical(colnames(region), c("bio", "ms")) ) ) {
    stop("region columns must have the names \"ms\" and \"bio\"\n")
  }
  if ( !(is.numeric(region) || is.character(region)) ) {
    stop("if a matrix then region must be either numeric or character\n")
  }
  if ( any( is.na(region) ) ) {
    stop("region cannot have missing values\n")
  }

  regList <- list( ms  = region[, "ms"],
                   bio = region[, "bio"] )
  return (regList)
}




format_reg_list <- function() {

  # Check that length and names are correct
  if (length(region) != 2) {
    stop("if a list then region must have exactly two elements")
  }
  if ( !( identical(names(region), c("ms", "bio"))
          || identical(names(region), c("bio", "ms")) ) ) {
    stop("region columns must have the names \"ms\" and \"bio\"\n")
  }

  # Delete the dimensions of an array which have only one level for list elements
  ms <- drop(region$ms)
  bio <- drop(region$bio)

  # Check element form, length, and types
  if ( !(is.strictvec(ms) && is.strictvec(bio)) ) {
    stop("the elements in region must be non-array vectors\n")
  }
  if (length(ms) != length(bio)) {
    stop("the elements in region must have the same length\n")
  }
  if ( !(is.numeric(ms) || is.character(ms))
       || !(is.numeric(bio) || is.character(bio)) ) {
    stop("the elements in region must each be either numeric or character\n")
  }
  if ( any( is.na(ms) ) || any( is.na(bio) ) ) {
    stop("region cannot have missing values\n")
  }

  regList <- list( ms  = ms,
                   bio = bio )
  return (regList)
}




