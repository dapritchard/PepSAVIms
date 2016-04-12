
#' Extract embedded mass spectrometry data
#'
#' Extract mass spectrometry data from an object with class \code{binMS}, class
#' \code{filterMS}, or class \code{msDat}.
#'
#' @param msObj An an object with class \code{binMS}, class \code{filterMS}, or
#'   class \code{msDat}.
#'
#' @param type A character string with value either "matrix", or "msDat".  If
#'   "matrix" is provided as the argument, then the mass-to-charge values,
#'   charge values, and mass spectrometry data are combined into a single matrix
#'   and returned.  If "msDat" is provided as the argument, then an \code{msDat}
#'   object containing this data is returned.
#'
#' @details A convenience function for extracting and inspecting the mass
#'   spectrometry data in a \code{binMS}, \code{filterMS}, or \code{msDat}
#'   object.  \code{binMS} and \code{filterMS} objects are lists that contain an
#'   msDat object, and specifying \code{"msDat"} for \code{type} merely returns
#'   the \code{msDat} element from the list for these classes of object.
#'   specifying \code{"msDat"} for an object with class "msDat" merely returns
#'   the argument, i.e. is the identity function.  When \code{"matrix"} is
#'   specified, then the elements in the embedded \code{msDat} object are
#'   combined into a single matrix using \code{cbind} and returned.
#'
#' @return Returns either a matrix containing the mass spectrometry data if
#'   \code{"matrix"} is specified as the argument to \code{type}, a an object
#'   with class \code{msDat} if \code{"msDat"} is specified as the argument to
#'   \code{type}.  See \code{Details} for more detail regarding the return
#'   objects.
#'
#' @export


extractMS <- function(msObj, type="matrix") {

  # Check args are of the right type
  if (missing(msObj)) {
    stop("Must provide an argument for msObj")
  }
  else if (!is.character(type)) {
    stop("type must have mode character")
  }

  # class() returns a character vector regardless of input
  class_nm <- class(msObj)
  msDatObj <- switch(class_nm,
                     binMS    = msObj$msObj,
                     filterMS = msObj$msObj,
                     msDat    = msObj,
                     stop("msObj must be an object of class ",
                          "\"binMS\", \"filterMS\", or \"msDat\""))

  # Return data in the desired form
  switch(type,
         matrix = with(msDatObj, cbind(mtoz, chg, ms)),
         msDat  = msDatObj,
         stop("type must have a value of \"matrix\" or \"msDat\""))
}




#' Mass spectrometry abundances fraction names
#'
#' Retrieve or set the mass spectrometry abundances fraction names for objects
#' of class \code{binMS}, \code{filterMS}, or \code{msDat}.
#'
#' @param msObj An object of class \code{binMS}, \code{filterMS}, or
#'   \code{msDat}
#'
#' @param value Either \code{NULL} or a character vector or an object that can
#'   be coerced to a character vector.  If non-\code{NULL}, then must be of the same
#'   length as the number of fractions in the mass spectrometry data.
#'
#' @export

colnamesMS <- function(msObj) {
  # Error-checking performed in extractMS
  msDatObj <- extractMS(msObj, "msDat")
  # returns NULL when msDatObj is NULL
  colnames(msDatObj$ms)
}


#' @rdname colnamesMS
#'
#' @export

`colnamesMS<-` <- function(msObj, value) {
  
  # Note: don't need to check if msObj, value are missing.  R won't invoke this
  # function unless both args are present

  # If non-NULL, ensure that value is a character vector
  if (!is.null(value) && !is.character(value)) {
    # Will throw an error if cannot coerce; use default error message
    value <- as.character(value)
  }
  
  # Ensure that msObj of a valid class.  note: class() returns a character
  # vector regardless of input
  class_nm <- class(msObj)
  if (!identical(class_nm, "binMS")
      && !identical(class_nm, "filterMS")
      && !identical(class_nm, "msDat")) {
    stop("msObj must be an object of class \"binMS\", \"filterMS\", or \"msDat\"")
  }

  # Extract dimnames from mass spectrometry abundances matrix
  if (identical(class_nm, "msDat")) {
    dn <- dimnames(msObj$ms)
    nfrac <- NCOL(msObj$ms)
  }
  else {
    dn <- dimnames(msObj$msObj$ms)
    nfrac <- NCOL(msObj$msObj$ms)
  }

  # Ensure that length of replacement is correct
  if (!is.null(value) && !identical(length(value), nfrac)) {
    stop("Length of replacement vector (", length(value), ") ",
         "not equal to the number of fractions in mass spectrometry data (", nfrac, ")")
  }
  
  # case: dimnames from ms matrix are null; if value is null then just return,
  # otherwise create container to store new values
  if (is.null(dn)) {
    if (is.null(value)) {
      return(msObj)
    }
    dn <- vector("list", 2)
  }
  
  # Store value in dn
  if (is.null(value)) {
    dn[2L] <- list(NULL)
  } else {
    dn[[2L]] <- value
  }

  # Write updated dimnames to ms matrix
  if (identical(class_nm, "msDat")) {
    dimnames(msObj$ms) <- dn
  } else {
    dimnames(msObj$msObj$ms) <- dn
  }
  
  return (msObj)
}
