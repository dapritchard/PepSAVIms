
#' Rank compounds using the Elastic Net path
#'
#' Returns identifying information for the compounds in the order in which they
#' first enter the Elastic Net model
#'
#' @param msObj An object of class \code{{msDat} containing mass spectrometry
#'   abundances data and identifying information.  Note that this includes
#'   objects created by the functions \code{binMS}, \code{filterMS}, and
#'   \code{msDat}.
#'
#' @param bioact Either a numeric vector or matrix, or a data frame providing
#'   bioactivity data.  If a numeric vector, then it is assumed that each entry
#'   corresponds to a particular fraction.  If the data is 2-dimensional, then
#'   it is assumed that each column corresponds to a particular fraction, and
#'   that each row corresponds to a particular bioactivity replicate.
#'
#' @param region_ms Either \code{NULL}, or a vector either of mode character or
#'   mode numeric providing information specifying which fractions from the mass
#'   spectrometry abundances data are to be included in the data analysis.  If
#'   \code{NULL}, then it is assumed that the entirety of the mass spectrometry
#'   abundances data encapsulated in the argument to \code{msObj} is to be
#'   included in the analysis.  If numeric then the entries should provide the
#'   indices for the region of interest in the mass spectrometry data (i.e. the
#'   indices of the columns corresponding to the appropriate fractions in the
#'   data).  If character then the entries should uniquely specify the region of
#'   interest through partial string matching (i.e. the names of the columns
#'   corresponding to the appropriate fractions in the data).  The methods
#'   \code{dim}, \code{dimnames}, and \code{colnamesMS} can be used as
#'   interfaces to the mass spectrometry data encapsulated in \code{msObj}.
#'
#' @param region_bio Either \code{NULL}, or a vector either of mode character or
#'   mode numeric providing information specifying which fractions from the
#'   bioactivity data are to be included in the data analysis.  If \code{NULL},
#'   then it is assumed that the entirety of bioactivity data provided as the
#'   argument to \code{bioact} is to be included in the analysis.  If numeric
#'   then the entries should provide the indices for the region of interest in
#'   the bioactivity data (i.e. the indices of the columns corresponding to the
#'   appropriate fractions in the data).  If character then the entries should
#'   uniquely specify the region of interest through partial string matching
#'   (i.e. the names of the columns corresponding to the appropriate fractions
#'   in the data).
#'
#' @param lambda A single nonnegative numeric value providing the quadratic
#'   penalty mixture parameter argument for the elastic net model.  The elastic
#'   net fits the least squares model with penalty function
#'   \deqn{\gamma|\beta|_1 + \lambda|\beta|^2} where \eqn{\beta} is the vector
#'   of regression coefficients and \eqn{\gamma, \lambda \ge 0}.  \code{rankEN}
#'   constructs a list of candidate compounds by tracking the entrance of
#'   compounds into the elastic net model as \eqn{\gamma} is decreased from
#'   \eqn{\infty} to \eqn{0}.
#'
#' @param pos_only Either \code{TRUE} or \code{FALSE}; specifies whether the
#'   list of candidate compounds that the algorithm produces should include only
#'   those compounds that are positively correlated with bioactivity levels, or
#'   conversely should include all compounds.  The correlation is calculated
#'   using only observations from the region of interest, and when bioactivity
#'   replicates are present, the within-fraction replicates are averaged prior
#'   to calculation.
#'
#' @param ncomp Either \code{NULL}, or a numeric value no less than 1 specifying
#'   the maximum number of candidate compounds that the function should report.
#'   When \code{NULL}, this is taken to mean that all compounds that enter the
#'   model should be reported, possibly after removing compounds nonpositively
#'   correlated with bioactivity levels, as specified by \code{pos_only}.
#'
#' @details \code{rankEN} prepares the data by extracting the region of interest
#'   from the mass spectrometry abundance data and from the bioactivity data.
#'   If bioactivity replicates are present, then the within-fraction replicates
#'   are averaged.  Once the data has been converted into the appropriate form,
#'   then an elastic net model is fitted by invoking the \code{enet} function
#'   from the \code{elasticnet} package, and an ordered list of candidate
#'   compounds is constructed such that compounds are ranked by the order in
#'   which they first enter the model.  The list may be filtered and / or pruned
#'   before being returned to the user, as determined by the arguments to
#'   \code{pos_only} and \code{ncomp}.
#'
#' @return Returns an object of class \code{rankEN}.  This object is a
#'   \code{list} with elements described below.  The class is equipped with a
#'   \code{print}, \code{summary}, and \code{extract_candidates} function.
#'
#'   \describe{
#'
#'   \item{\code{mtoz}}{ A vector providing the mass-to-charge values of the
#'   candidate compounds, such that the \code{k}-th element of the vector
#'   provides the mass-to-charge value of the \code{k}-th compound to enter the
#'   elastic net model, possibly after removing compounds nonpositively
#'   correlated with bioactivity levels. }
#'
#'   \item{\code{charge}}{ A vector providing the charge state of the candidate
#'   compounds, such that the \code{k}-th element of the vector provides the
#'   charge state of the \code{k}-th compound to enter the elastic net model,
#'   possibly after removing compounds nonpositively correlated with bioactivity
#'   levels. }
#'
#'   \item{\code{comp_cor}}{ A vector providing the correlation between each of
#'   the candidate compounds and the bioactivity levels, such that the
#'   \code{k}-th element of the vector provides the correlation between the
#'   \code{k}-th compound to enter the elastic net model and the bioactivity
#'   levels, possibly after removing compounds nonpositively correlated with
#'   bioactivity levels. }
#'
#'   \item{\code{enet_fit}}{ The fitted model object produced by \code{rankEN}'s
#'   internal invokation of the \code{enet} function from the \code{elasticnet}
#'   package.}
#'
#'   \item{\code{summ_info}}{ A list containing information related to the data
#'   used to fit the elastic net model; used by the summary function. }
#'
#'   }
#'
#' @export


rankEN <- function(msObj, bioact, region_ms=NULL, region_bio=NULL, lambda,
                   pos_only=TRUE, ncomp=NULL) {

  # Ensure that arguments are of the right type
  rankEN_check_valid_input(msObj, bioact, region_ms, region_bio, lambda, pos_only, ncomp)


  # Mung data into the right form ----------------------------------------------
  
  # Obtain msDat obj
  msDatObj <- extractMS(msObj, type="msDat")
  if (is.null(msDatObj)) {
    stop("mass spec object encapsulated by msObj cannot be NULL", call.=FALSE)
  }
  ms <- msDatObj$ms
  
  # If we have a vector convert to a 1-row matrix.  Leave unchanged otherwise.
  bioact <- rankEN_vector_to_matrix(bioact)
  
  # Extract the region of interest for the ms data
  ms <- extract_var(ms, region_ms, TRUE)
  bio <- extract_var(bioact, region_bio, TRUE)
  
  # Check for missing and that dimensions match
  rankEN_check_regr_args(ms, bio)

  # Convert ms to form where rows are an observation (i.e. fraction) and cols
  # arms_te a variable (i.e. a compound)
  ms_regr <- t(ms)
  
  # Obtain the mean of the bioactivity replicates and convert to a vector
  bio_regr <- colMeans(bio)
  

  # Fit model ------------------------------------------------------------------

  # Calculate the elastic net path
  enet_fit <- tryCatch({
    elasticnet::enet(ms_regr, bio_regr, lambda)
  }, warning = function(war) {
    warning("message produced by call to enet from package elasticnet ==>\n", war, call.=FALSE)
  }, error = function(err) {
    stop("message produced by call to enet from package elasticnet ==>\n", err, call.=FALSE)
  })
  

  # Extract compound entrance results ------------------------------------------

  # Obtain indices for compounds as they first enter the model
  comp_idx <- rankEN_comp_entrance(enet_fit)

  # Correlation between chosen compounds and mean bioact levels in region
  comp_cor <- rankEN_comp_cor(ms_regr, bio_regr)

  # Filter compounds by correlation and by how many we wish to keep
  comp_idx_out <- rankEN_filter_compIdx(comp_idx, comp_cor, ncomp, pos_only)


  # Construct return object ----------------------------------------------------

  # Extract mass spec and bioactivity fraction names  
  ms_nm <- colnames(ms)
  if (is.null(ms_nm)) {
    ms_nm <- paste0("ms", 1:ncol(ms))
  }
  bio_nm <- colnames(bio)
  if (is.null(bio_nm)) {
    bio_nm <- paste0("bio", 1:ncol(bio))
  }  

  # Create info for the summary function
  summ_info <- list(
    data_dim  = list(reg  = ncol(ms),
                     comp = nrow(ms),
                     repl = nrow(bio)),
    region_nm = list(ms  = ms_nm,
                     bio = bio_nm),
    lambda    = lambda,
    pos_only  = pos_only,
    ncomp     = ncomp
  )

  # Construct output object
  outDat <- list( mtoz      = msDatObj$mtoz[comp_idx_out],
                  charge    = msDatObj$chg[comp_idx_out],
                  comp_cor  = comp_cor[comp_idx_out],
                  enet_fit  = enet_fit,
                  summ_info = summ_info )

  structure(outDat, class="rankEN")
}




#' Extract candidate compounds
#'
#' Extract an ordered list of candidate compounds from a \code{rankEN} object.
#' The list is presented in the form of a \code{data.frame}, such that each row
#' provides the identifying information for a particular candidate compound, and
#' with the rows arranged in the order that the compounds entered the elastic
#' net model (i.e. row 1 is the earliest, row 2 the 2nd earliest, etc.).  The
#' columns of the \code{data.frame} provide the mass-to-charge information,
#' charge information, and possibly the correlation between the compound and the
#' within-fraction average of the bioactivity replicates in the region of
#' interest.
#'
#' @param rankEN_obj An object of class \code{rankEN}.
#'
#' @param include_cor Either \code{TRUE} or \code{FALSE}, specifying whether a
#'   column should be included in the returning \code{data.frame} providing the
#'   correlation between the compound and the within-fraction average of the
#'   bioactivity replicates in the region of interest.
#'
#' @export


extract_candidates <- function(rankEN_obj, include_cor=TRUE) {

  if (!identical(class(rankEN_obj), "rankEN")) {
    stop("rankEN_obj must be of class rankEN")
  }
  else if (!identical(include_cor, TRUE) && !identical(include_cor, FALSE)) {
    stop("include_cor must be either TRUE or FALSE")
  }

  out <- data.frame(rankEN_obj$mtoz, rankEN_obj$charge)
  if (include_cor) {
    out$comp_cor <- rankEN_obj$comp_cor
  }

  return (out)
}




#' Basic information for class \code{rankEN}
#'
#' Displays the data dimensions used to fit the elastic net model
#'
#' @param x An object of class \code{rankEN}
#'
#' @param ... Arguments passed to dot-dot-dot are ignored
#'
#' @export


print.rankEN <- function(x, ...) {

  cat(sep="",
      "An object of class rankEN.\n",
      "Use summary to print a list of the compounds entering the model.\n",
      "Use extract_candidates to extract the compound info as a data.frame.\n\n")

  dd_char <- format(unlist(x$summ_info$data_dim), big.mark=",", justify="right")
  cat(sep="",
      "Data dimensions:\n",
      "----------------\n",
      "    region of interest:     ", dd_char[1], "\n",
      "    candidate compounds:    ", dd_char[2], "\n",
      "    bioactivity replicates: ", dd_char[3], "\n\n")
}




#' Overview of the elastic net selection process
#'
#' Prints a description of the elastic net variable selection process.  Includes
#' the dimensions used to fit the elastic net model, the fraction names for the
#' mass spectrometry and the bioactivity data in the region of interest, the
#' parameter specifications for the model, and a table with the identifying
#' information of the candidate compounds produced by the model fit.
#'
#' @param object An object of class \code{rankEN}.
#'
#' @param max_comp_print A numeric value >= 1 specifying the maximum number of
#'   compounds to print
#'
#' @param ... Arguments passed to dot-dot-dot are ignored
#'
#' @export

summary.rankEN <- function(object, max_comp_print=20L, ...) {

  # Check argument to max_comp_print
  if (!is.numeric(max_comp_print)) {
    stop("max_comp_print mus be of mode numeric")
  }
  else if ( !identical(length(max_comp_print), 1L) ) {
    stop("max_comp_print must have length 1")
  }
  else if (max_comp_print < 1L) {
    stop("max_comp_print cannot have value less than 1")
  }

  # Create links for convenience
  summ_info <- object$summ_info
  mtoz      <- object$mtoz
  chg       <- object$charge
  ccor      <- object$comp_cor
  
  # Print restricted mass spectrometry and bioactivity data dimensions
  dd_char <- format(unlist(summ_info$data_dim), big.mark=",", justify="right")
  cat(sep="",
      "\n",
      "Data dimensions:\n",
      "----------------\n",
      "    region of interest:     ", dd_char[1], "\n",
      "    candidate compounds:    ", dd_char[2], "\n",
      "    bioactivity replicates: ", dd_char[3], "\n\n")
  
  # Print a table with the fraction names used for the mass spectrometry and
  # bioactivity data
  region_nm_df <- data.frame(summ_info$region_nm)
  colnames(region_nm_df) <- c("Mass spec", "Bioactivity")
  region_table <- capture.output( print(region_nm_df, row.names=FALSE, print.gap=2)  )
  lead_blanks <- rep(" ", min(getOption("width") - max(nchar(region_table)), 2L))
  cat(sep="",
      "Fractions included in region of interest:\n",
      "-----------------------------------------\n")
  for (row_entry in region_table) {
    cat(lead_blanks, row_entry, "\n", sep="")
  }
  cat("\n")

  # Print arguments to rankEN
  parm_char <- format(c(format(summ_info$lambda, digits=4),
                        ifelse(summ_info$pos_only, "yes", "no"),
                        ifelse(is.null(summ_info$ncomp),
                               "all",
                               as.character(summ_info$ncomp))),
                      justify="right")
  cat(sep="",
      "Parameter arguments provided to rankEN:\n",
      "---------------------------------------\n",
      "    Quadratic penalty parameter:          ", parm_char[1], "\n",
      "    Consider only positive correlations:  ", parm_char[2], "\n",
      "    Max number of candidate compounds:    ", parm_char[3], "\n\n")

  # Print a table of the m/z, charge, and correlation values for selected comps
  maxpr_char <- format(as.integer(max_comp_print), big.mark=",")
  if (length(mtoz) <= max_comp_print) {
    cat(sep="",
        "Compounds in order of entrance (all compounds, earliest at top):\n",
        "----------------------------------------------------------------\n")
  }
  else {
    cat(sep="",
        "Compounds in order of entrance (first ", maxpr_char, " compounds, earliest at top):\n",
        rep("-", 67 + nchar(maxpr_char)), "\n")
  }
  comp_df <- data.frame(mtoz, chg, format(round(ccor, 4), nsmall=4))
  colnames(comp_df) <- c("Mass spec", "Charge", "Correlation")
  comp_table <- capture.output( print(comp_df, row.names=FALSE, print.gap=2) )
  lead_blanks <- rep(" ", min(getOption("width") - max(nchar(comp_table)), 2L))
  for (i in 1:(min(max_comp_print, length(mtoz)) + 1)) {
    cat(lead_blanks, comp_table[i], "\n", sep="")
  }
  cat("\n")
}





