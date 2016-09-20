
#' Constructor for class \code{msDat}
#'
#' Creates a data structure encapsulating the mass spectrometry intensity
#' readings as well as identifying information
#'
#' @param mass_spec Either a \code{matrix} or \code{data.frame}.  This object
#'     must contain mass spectrometry abundances, and may optionally contain
#'     mass-to-charge values, charge state information, or additional extraneous
#'     variables.  The mass spectrometry data is expected to be in a form with
#'     each column corresponding to a variable and each row corresponding to a
#'     mass-to-charge level.
#'
#'     For example, suppose that a collection of mass spectrometry intensity
#'     observations has provided data for 50 fractions across 20,000
#'     mass-to-charge values.  Then the input for \code{mass_spec} should be a
#'     \code{matrix} or \code{data.frame} with 20,000 rows and 50 or more
#'     columns.  The additional columns beyond the 50 containing the mass
#'     spectrometry intensities can be the mass-to-charge data, the charge data,
#'     or other extraneous variables (the extraneous variables will be discarded
#'     when constructing the \code{msDat} object).
#'
#' @param mtoz A vector of either length 1 or length equal to the number of
#'     mass-to-charge values for which mass spectrometry data was collected, and
#'     which helps identify the mass-to-charge values for this data in one of
#'     several ways.
#'
#'     One way to provide the information is to provide a numeric vector where
#'     each entry provides the mass-to-charge value for a corresponding row of
#'     mass spectrometry data.  Then the \code{k}-th entry of the vector would
#'     provide the mass-to-charge value for the \code{k}-th row of the mass
#'     spectrometry data.
#'
#'     A second way is to provide a single number which specifies the column
#'     index in the \code{matrix} or \code{data.frame} provided as the argument
#'     for the \code{mass_spec} parameter, such that this column contains the
#'     mass-to-charge information.
#'
#'     A third way is provide a single character string which provides the
#'     column name in the \code{matrix} or \code{data.frame} provided as the
#'     argument for the \code{mass_spec} parameter, such that this column
#'     contains the mass-to-charge information.  Partial matching is supported.
#'
#' @param charge The information for the \code{charge} parameter can be provided
#'     in the same manner as for the mass-to-charge values.
#'
#' @param ms_inten Either \code{NULL} or a vector either of mode character or
#'     mode numeric specifying which of the variables in the argument to
#'     \code{mass_spec} are to be retained as the mass spectrometry intensity
#'     data.  If \code{NULL}, then it is taken to mean that the entirety of the
#'     data in \code{mass_spec}, after removing variables in the data that are
#'     specified as arguments, is the mass spectrometry intensity data.  If it
#'     is a numeric vector, then the entries should provide the indices for the
#'     region of interest in the mass spectrometry data in the argument for
#'     \code{msObj}. If it is a character vector, then the entries should
#'     uniquely specify the region of interest through partial string matching.
#'
#' @details Since the mass spectrometry data could conceivably be available to
#'     the researcher in a variety forms, this function attempts to provide a
#'     uniform data structure for encapsulating this information.  It is the
#'     fundamental data structure containing the mass spectrometry data used
#'     internally by the \code{filterMS} and \code{rankEN} routines.  The
#'     external interface for \code{msDat} is provided to the user so that
#'     specifying the mass spectrometry information can be made in a distinct
#'     step from performing statistical analyses, which it is hoped makes
#'     interfaces for the downstream analysis routines simpler and more
#'     intuitive to use.
#'
#' @return Returns an object of class \code{msDat}.  This class is a \code{list}
#'     with elements described below.  The class is equipped with a \code{print}
#'     and \code{extractMS} function.
#'
#'     \describe{
#'
#'     \item{\code{ms}}{ A \code{matrix} containing mass spectrometry intensity
#'         readings. Each column provides the mass spectrometry values for a
#'         given fraction, and each row provides the mass spectrometry values
#'         for a given mass-to-charge ratio value across the fractions. }
#'
#'     \item{\code{mtoz}}{ A vector with length equal to the number of
#'         mass-to-charge values provided in the mass spectrometry data, such
#'         that the \code{k}-th entry in the vector provides the mass-to-charge
#'         value for the \code{k}-th row of mass spectrometry data }
#'
#'     \item{\code{chg}}{ A vector with length equal to the number of
#'         mass-to-charge values provided in the mass spectrometry data, such
#'         that the \code{k}-th entry in the vector provides the charge
#'         information for the \code{k}-th row of mass spectrometry data }
#'
#'     }
#'
#' @examples
#'
#' # Load mass spectrometry data
#' data(mass_spec)
#'
#' # Convert mass_spec from a data.frame to an msDat object
#' ms <- msDat(mass_spec = mass_spec,
#'             mtoz = "m/z",
#'             charge = "Charge",
#'             ms_inten = c(paste0("_", 11:43), "_47"))
#'
#' # Dimension of the data
#' dim(ms)
#'
#' # Print the first few rows and columns
#' ms[1:5, 1:2]
#'
#' # Let's change the fraction names to something more concise
#' colnames(ms) <- c(paste0("frac", 11:43), "frac47")
#'
#' # Print the first few rows and columns with the new fraction names
#' ms[1:5, 1:8]
#'
#' # Suppose there are some m/z levels that we wish to remove
#' ms <- ms[-c(2, 4), ]
#' # Print the first few rows and columns after removing rows 2 and 4
#' ms[1:5, 1:8]
#'
#' # Suppose that there was an instrumentation error and that we need to change
#' # some values
#' ms[1, paste0("frac", 12:17)] <- c(55, 57, 62, 66, 71, 79)
#' # Print the first few rows and columns after changing some of the values in
#' # the first row
#' ms[1:5, 1:10]
#'
#' @export


msDat <- function(mass_spec, mtoz, charge, ms_inten=NULL) {

    # Check that input is of the right form
    msDat_check_valid_input(mass_spec, mtoz, charge, ms_inten)

    # Obtain mass-to-charge, charge, and mass spectrometry abundances variables.
    # Note that we use the original version of the variables in each of the
    # extract_* calls.
    mtoz_ <- extract_var(mass_spec, mtoz)
    charge_ <- extract_var(mass_spec, charge)
    ms_inten <- extract_var(mass_spec, ms_inten, TRUE, mtoz, charge)

    # Row names for mass spectrometry data are their mtoz and charge states
    row.names(ms_inten) <- paste0(format(round(mtoz_, 4), nsmall=4),
                                  "/",
                                  charge_)

    # Construct msDat object
    outDat <- list(ms   = ms_inten,
                   mtoz = mtoz_,
                   chg  = charge_)

    structure(outDat, class="msDat")
}




# Class methods ----------------------------------------------------------------


#' Print method for class \code{msDat}
#'
#' Prints the mass spectrometry data encapsulated by the \code{msDat} object
#'
#' @param x An object of class \code{\link{msDat}}
#'
#' @param ... Arguments passed to dot-dot-dot are ignored
#'
#' @export

print.msDat <- function(x, ...) {

    print(x$ms)
}
