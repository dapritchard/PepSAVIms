
#' Extract embedded mass spectrometry data
#'
#' Extract mass spectrometry data from an object with class \code{binMS}, class
#' \code{filterMS}, or class \code{msDat}.
#'
#' @param msObj An an object with class \code{binMS}, class \code{filterMS}, or
#'     class \code{msDat}.
#'
#' @param type A character string with value either "matrix", or "msDat".  If
#'     "matrix" is provided as the argument, then the mass-to-charge values,
#'     charge values, and mass spectrometry data are combined into a single
#'     matrix and returned.  If "msDat" is provided as the argument, then an
#'     \code{msDat} object containing this data is returned.
#'
#' @details A convenience function for extracting and inspecting the mass
#'     spectrometry data in a \code{binMS}, \code{filterMS}, or \code{msDat}
#'     object.  \code{binMS} and \code{filterMS} objects are lists that contain
#'     an msDat object, and specifying \code{"msDat"} for \code{type} merely
#'     returns the \code{msDat} element from the list for these classes of
#'     object.  specifying \code{"msDat"} for an object with class "msDat"
#'     merely returns the argument, i.e. is the identity function.  When
#'     \code{"matrix"} is specified, then the elements in the embedded
#'     \code{msDat} object are combined into a single matrix using \code{cbind}
#'     and returned.
#'
#' @return Returns either a matrix containing the mass spectrometry data if
#'     \code{"matrix"} is specified as the argument to \code{type}, or an object
#'     with class \code{msDat} if \code{"msDat"} is specified as the argument to
#'     \code{type}.  See \code{Details} for more detail regarding the return
#'     objects.
#'
#' @export


extractMS <- function(msObj, type="matrix") {

    # Check msObj is of the right type
    if (missing(msObj)) {
        stop("Must provide an argument for msObj", call.=FALSE)
    }
    else if (!inherits(msObj, "msDat")) {
        stop("msObj must be of class \"msDat\"", call.=FALSE)
    }

    # Check that type is of the right type
    if (!identical(type, "matrix") && !identical(type, "msDat")) {
        stop("type must have a value of \"matrix\" or \"msDat\"", call.=FALSE)
    }

    # Point msObj to the location of the msDat object
    if (! identical(class(msObj), "msDat")) {
        msObj <- msObj$msDatObj
    }

    # Return data in the desired form
    switch(type,
           msDat  = msObj,
           matrix = cbind(mtoz   = msObj$mtoz,
                          charge = msObj$chg,
                          inten  = msObj$ms))
}


#' @export
dimnames.msDat <- function(x) {

    x <- extractMS(x, "msDat")
    dimnames(x$ms)
}


#' @export
`dimnames<-.msDat` <- function(x, value) {

    # case: one of the classes that decorates msDat object.  Recursively call
    # with inner msDat object.
    if (!identical(class(x), "msDat")) {
        dimnames(x$msDatObj) <- value
        return (x)
    }

    dimnames(x$ms) <- value

    x
}


#' @export
`[.msDat` <- function(msObj, i, j) {

    # case: one of the classes that decorates msDat object.  Recursively call
    # with inner msDat object.
    if (!identical(class(msObj), "msDat")) {
        msObj$msDatObj <- msObj$msDatObj[i, j]
        return (msObj)
    }

    msObj$ms <- msObj$ms[i, j, drop=FALSE]
    msObj$mtoz <- msObj$mtoz[i]
    msObj$chg <- msObj$chg[i]

    msObj
}


#' @export
`[<-.msDat` <- function(msObj, i, j, value) {

    # case: one of the classes that decorates msDat object.  Recursively call
    # with inner msDat object.
    if (!identical(class(msObj), "msDat")) {
        msObj$msDatObj[i, j] <- value
        return (msObj)
    }

    msObj$ms[i, j] <- value

    msObj
}


#' @export
dim.msDat <- function(x) {
    msDatObj <- extractMS(x, type="msDat")
    dim(msDatObj$ms)
}
