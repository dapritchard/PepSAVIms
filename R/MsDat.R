
#' Constructor for class \code{msDat}
#'
#' Creates a data structure encapsulating the mass spectrometry intensity
#' readings as well as identifying information
#'
#' @param mass_spec Either a \code{matrix} or \code{data.frame}.  This object
#'   must contain mass spectrometry abundances, and may optionally contain
#'   mass-to-charge values, charge state information, or additional extraneous
#'   variables.  The mass spectrometry data is expected to be in a form with
#'   each column corresponding to a variable and each row corresponding to a
#'   mass-to-charge level.
#'
#'   For example, suppose that a collection of mass spectrometry intensity
#'   observations has provided data for 50 fractions across 20,000
#'   mass-to-charge values.  Then the input for \code{mass_spec} should be a
#'   \code{matrix} or \code{data.frame} with 20,000 rows and 50 or more columns.
#'   The additional columns beyond the 50 containing the mass spectrometry
#'   intensities can be the mass-to-charge data, the charge data, or other
#'   extraneous variables (the extraneous variables will be discarded when
#'   constructing the \code{msDat} object).
#'
#' @param mtoz A vector of either length 1 or length equal to the number of
#'   mass-to-charge values for which mass spectrometry data was collected, and
#'   which helps identify the mass-to-charge values for this data in one of
#'   several ways.
#'
#'   One way to provide the information is to provide a numeric vector where
#'   each entry provides the mass-to-charge value for a corresponding row of
#'   mass spectrometry data.  Then the \code{k}-th entry of the vector would
#'   provide the mass-to-charge value for the \code{k}-th row of the mass
#'   spectrometry data.
#'
#'   A second way is to provide a single number which specifies the column index
#'   in the \code{matrix} or \code{data.frame} provided as the argument for the
#'   \code{mass_spec} parameter, such that this column contains the
#'   mass-to-charge information.
#'
#'   A third way is provide a single character string which provides the column
#'   name in the \code{matrix} or \code{data.frame} provided as the argument for
#'   the \code{mass_spec} parameter, such that this column contains the
#'   mass-to-charge information.  Partial matching is supported.
#'
#' @param charge The information for the \code{charge} parameter can be provided
#'   in the same manner as for the mass-to-charge values.
#'
#' @param ms_inten Either \code{NULL} or a vector either of mode character or
#'   mode numeric specifying which of the variables in the argument to
#'   \code{mass_spec} are to be retained as the mass spectrometry intensity
#'   data.  If \code{NULL}, then it is taken to mean that the entirety of the
#'   data in \code{mass_spec}, after removing variables in the data that are
#'   specified as arguments, is the mass spectrometry intensity data.  If it is
#'   a numeric vector, then the entries should provide the indices for the
#'   region of interest in the mass spectrometry data in the argument for
#'   \code{msObj}. If it is a character vector, then the entries should uniquely
#'   specify the region of interest through partial string matching.
#'
#' @details Since the mass spectrometry data could conceivably be available to
#'   the researcher in a variety forms, this function attempts to provide a
#'   uniform data structure for encapsulating this information.  It is the
#'   fundamental data structure containing the mass spectrometry data used
#'   internally by the \code{filterMS} and \code{rankEN} routines.  The external
#'   interface for \code{msDat} is provided to the user so that specifying the
#'   mass spectrometry information can be made in a distinct step from
#'   performing statistical analyses, which it is hoped makes interfaces for the
#'   downstream analysis routines simpler and more intuitive to use.
#'
#' @return Returns an object of class \code{msDat}.  This class is a \code{list}
#'   with elements described below.  The class is equipped with a \code{print}
#'   and \code{extractMS} function.
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


msDat <- function(mass_spec, mtoz, charge, ms_inten=NULL) {

  # Check that input is of the right form
  msDat_check_valid_input(mass_spec, mtoz, charge, ms_inten)

  # Obtain mass-to-charge, charge, and mass spectrometry abundances variables
  dmtoz <- extract_var(mass_spec, mtoz)
  dcharge <- extract_var(mass_spec, charge)
  ms_inten <- extract_var(mass_spec, ms_inten, TRUE, mtoz, charge)

  # Construct msDat object
  outDat <- list(ms   = ms_inten,
                 mtoz = dmtoz,
                 chg  = dcharge)
  
  structure(outDat, class="msDat")
}




# Class methods ----------------------------------------------------------------


#' Basic information for class \code{msDat}
#'
#' Displays the number of candidate compounds left in the data
#'
#' @param x An object of class \code{\link{msDat}}
#'
#' @param ... Arguments passed to dot-dot-dot are ignored
#'
#' @export

print.msDat <- function(x, ...) {

  print.default(x)
  
  cat("\nAn object of class \"msDat\" with ", format(NROW(x$ms), big.mark=","),
      " compounds and ", NCOL(x$ms), " fractions.\n", sep="")
  cat("Use extractMS to column-bind the data together into a single matrix.\n\n")

}




# Helper functions -------------------------------------------------------------


msDat_check_valid_input <- function(mass_spec, mtoz, charge, ms_inten) {

  ## Check for missing arguments
  
  all_var_nm <- c("mass_spec", "mtoz", "charge", "ms_inten")
  for (var_nm in all_var_nm) {
    if (!eval(substitute(hasArg(var_nm)))) {
      stop("Must provide an argument for ", var_nm, call.=FALSE)
    }
    # Check that an object exists for provided argument 
    tryCatch(get(var_nm), error = function(err) {
      err <- as.character(err)
      obj_nm <- regmatches(err, gregexpr("(?<=\')(.*?)(?=\')", err, perl=TRUE))[[1L]]
      stop("object \'", obj_nm, "\' not found for ", var_nm, call.=FALSE)
    })
  }

  ## Check mass_spec

  if (!is.matrix(mass_spec) && !is.data.frame(mass_spec)) {
    stop("mass_spec must be either a matrix or data.frame", call.=FALSE)
  }
  # The decision to not allow 0 rows is a design decision, and the reason that
  # we can't have 1 row is because that if a number is supplied then we don't
  # know if it is the value of the data or an index being supplied
  else if (NROW(mass_spec) <= 1L) {
    stop("mass_spec must have 2 or more rows", call.=FALSE)
  }
  # Don't allow 0 columns in data as a design decision
  else if (identical(NCOL(mass_spec), 0L)) {
    stop("mass_spec cannot have 0 columns", call.=FALSE)
  }
  
  ## Check mtoz, charge

  for (var_nm in c("mtoz", "charge")) {
    x <- get(var_nm)
    if (!is.numeric(x) && !is.character(x)) {
      stop(var_nm, " must be either of mode numeric or character", call.=FALSE)
    }
  }

  ## Check ms_inten

  if (!is.null(ms_inten) && !is.numeric(ms_inten) && !is.character(ms_inten)) {
    stop("ms_inten must be either NULL or of mode numeric or character", call.=FALSE)
  }
}







