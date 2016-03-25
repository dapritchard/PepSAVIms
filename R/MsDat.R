
#' Constructor for class \code{msDat}
#'
#' Creates a data structure encapsulating the mass spectrometry intensity
#' readings as well as identifying information
#'
#' @param mass_spec Either a \code{matrix} or \code{data.frame}.  This object
#'   must contain mass spectrometry abundances, and may optionally contain
#'   mass-to-charge values and/or charge state information.  The mass
#'   spectrometry data is expected to be in a form such that a given column
#'   provides the mass spectrometry intensity values for a particular fraction.
#'   Then either 0, 1, or 2 additional columns may be included in the parameter
#'   input containing a possible column for the mass-to-charge values and a
#'   possible column for the charge information. Any ordering of the columns is
#'   allowed.
#'
#'   For example, suppose that a collection of mass spectrometry readings has
#'   provided data for 50 fractions across 20,000 mass-to-charge values.  Then
#'   the input for \code{mass_spec} should be a \code{matrix} or
#'   \code{data.frame} with 20,000 rows and one of 50, 51, or 52 columns.  Each
#'   row in the 50 columns containing the mass spectrometry readings should be
#'   for the same mass-to-charge value and charge information.  If columns are
#'   included for either the mass-to-charge values or charge information then
#'   the entries in these columns should provide the identifying information for
#'   the mass spectrometry data in the same row.
#'
#' @param mtoz A vector of either length 1 or length equal to the number of
#'   mass-to-charge values for which mass spectrometry data was collected, and
#'   which helps identify the mass-to-charge values for this data in one of
#'   several ways.
#'
#'   One way to provide the information is to provide a vector where each entry
#'   provides the mass-to-charge value for a corresponding row of mass
#'   spectrometry data.  Then the \code{k}-th entry of the vector would provide
#'   the mass-to-charge value for the \code{k}-th row of the mass spectrometry
#'   data.
#'
#'   A second way is to provide a single number which provides the column number
#'   in the input to the \code{mass_spec} parameter for a column which contains
#'   this information.
#'
#'   A third way is provide a single character string which provides the column
#'   name in the input to the \code{mass_spec} parameter for a column which
#'   contains this information.  Partial matching is allowed.
#'
#' @param charge The information for the \code{charge} parameter can be provided
#'   in the same manner as for the mass-to-charge values.
#'
#' @return Returns an object of class \code{msDat}.  This class is a \code{list}
#'   with elements described below.  The class is equipped with a summary
#'   function.
#'
#'   \describe{
#'
#'   \item{\code{ms}}{ A \code{matrix} containing mass spectrometry intensity
#'   readings. Each column provides the mass spectrometry values for a given
#'   fraction, and each row provides the mass spectrometry values for a given
#'   mass-to-charge ratio value across the fractions. }
#'
#'   \item{\code{mtoz}}{ A vector with length equal to the number of
#'   mass-to-charge values provided in the mass spectrometry data, such that the
#'   \code{k}-th entry in the vector provides the mass-to-charge value for the
#'   \code{k}-th row of mass spectrometry data }
#'
#'   \item{\code{chg}}{ A vector with length equal to the number of
#'   mass-to-charge values provided in the mass spectrometry data, such that the
#'   \code{k}-th entry in the vector provides the charge information for the
#'   \code{k}-th row of mass spectrometry data }
#'
#'   }
#'
#' @export


msDat <- function(mass_spec, mtoz, charge) {

  # Perform some checks for validity of input.  Some forms of invalid input may
  # still exist that are checked for as the function progresses.
  checkValInp_msDat(mass_spec, mtoz, charge)

  # Delete the dimensions of an array which have only one level (if necessary)
  #mass_spec <- drop(mass_spec)
  #mtoz <- drop(mtoz)
  #charge <- drop(charge)

  # cmpInfo: a list with seperate vectors containing the mass-to-charge values
  # and charge information, as well as integer values providing (if applicable)
  # the column number in mass_spec containing this information
  cmpInfo <- getCmpInfo(mass_spec, mtoz, charge)

  # keepIdx: indexes the columns in mass_spec that contain the mass spectrometry
  # intensity data
  keepIdx <- setdiff(seq_len(ncol(mass_spec)),
                     c(cmpInfo$mtoz$loc, cmpInfo$chg$loc))

  outDat <- list( ms   = as.matrix( mass_spec[, keepIdx] ),
                  mtoz = cmpInfo$mtoz$val,
                  chg  = cmpInfo$chg$val )

  if ( !is.numeric(outDat$ms) ) {
    stop("mass spectrometry data must be numeric\n")
  }

  structure(outDat, class="msDat")
}




#' Basic statistics for mass spectrometry data
#'
#' Summary function for class \code{msDat}

summary.msDat <- function(msDat) {

  cat("\nThe mass spectrometry data has",
      format(ncol(msDat$ms), big.mark=","),
      "fractions across",
      format(nrow(msDat$ms), big.mark=","),
      "mass-to-charge values\n\n")
}




#' Check for valid msDat arguments
#'
#' Check that arguments \code{mass_spec}, \code{mtoz}, and \code{charge} are
#' of the right data type
#'
#' @inheritParams msDat

checkValInp_msDat <- function(mass_spec, mtoz, charge) {

  # Check if parameters are one of matrices / data frames / non-list vectors,
  # as appropriate

  if ( !(is.matrix(mass_spec) || is.data.frame(mass_spec)) ) {
    stop("mass_spec must be a matrix or data frame\n")
  }
  else if ( !( is.strictVec(drop(mtoz)) ) ) {
    stop("mtoz must be a non-list vector\n")
  }
  else if ( !( is.strictVec(drop(charge)) ) ) {
    stop("charge must be a non-list vector\n")
  }

  # TODO: check for missing.  Is this not performed elsewhere?
}




# TODO: there is a is.strictvec in RankLasso_Helper.R.  Should merge the two.

is.strictVec <- function(x) {
  return( is.vector(x) && !is.list(x) )
}




# TODO: need documentation for this

getCmpInfo <- function(mass_spec, mtoz, charge) {

  # Calculate number of compounds
  nCmp <- nrow(mass_spec)

  # Check if mtoz is numeric or character.  Pre: mtoz is a non-list vector
  if ( !(is.numeric(mtoz) || is.character(mtoz)) ) {
    stop("mtoz must be numeric or character\n")
  }
  else if ( !( identical(length(mtoz), 1L) || identical(length(mtoz), nCmp) ) ) {
    stop("mtoz must have either have a length of 1 or length equal to the ",
         "number of compounds", call.=FALSE)
  }

  # Check if charge is numeric or character.  Pre: charge is a non-list vector
  if ( !(is.numeric(charge) || is.character(charge)) ) {
    stop("charge must be numeric or character\n")
  }
  else if ( !( identical(length(mtoz), 1L) || identical(length(mtoz), nCmp) ) ) {
    stop("charge must have either have a length of 1 or length equal to the ",
         "number of compounds", call.=FALSE)
  }

  # Initialize a list to save mass spectrometry identifying information (and
  # possibly location of this information is mass_spec) into
  outDat <- list()
  # Create iterable object for for loop
  varList <- list(mtoz=mtoz, chg=charge)


  # Each iteration in loop creates an entry in outDat where the entry is a list
  # containing elements loc and val.  loc is an index for the column in mass_spec
  # containing the mass-to-charge or charge information (which may be numeric(0) if
  # not in mass_spec).  val a vector with the actual information.

  for (i in 1:length(varList)) {
    outDat[[i]] <- list()
    thisVar <- varList[[i]]

    # nCmp length atomic vector
    if ( identical(length(thisVar), nCmp) ) {
      outDat[[i]]$loc <- numeric(0)
      outDat[[i]]$val <- thisVar
    }
    # length 1 numeric
    else if (is.numeric(thisVar)) {
      if ((thisVar < 1) || (thisVar > ncol(mass_spec))) {
        stop(paste("not a valid column number for", names(varList)[i], "\n"))
      }
      outDat[[i]]$loc <- as.integer(thisVar)
      outDat[[i]]$val <- mass_spec[, thisVar]
    }
    # length 1 character
    else if (is.character(thisVar)) {
      if ( !(thisVar %in% colnames(mass_spec)) ) {
        stop(paste("invalid column name given for", names(varList)[i], "\n"))
      }
      outDat[[i]]$loc <- which(thisVar == colnames(mass_spec))
      outDat[[i]]$val <- mass_spec[, thisVar]
    }
    else {
      # Should not reach here.  All cases should be covered either in the error
      # checking or conditionals.
    }
  } # End loop entering data into outDat

  names(outDat) <- names(varList)

  return (outDat)
}









