
#' Rank compounds via the Elastic Net path
#'
#' Returns identifying information for the compounds in the order in which the
#' corresponding regression coefficient for a given compound first becomes
#' nonzero as part of the Elastic Net path
#'
#' @param msObj An object of class \code{msDat} containing mass spectrometry
#'     abundances data and identifying information.  Note that this includes
#'     objects created by the functions \code{binMS}, \code{filterMS}, and
#'     \code{msDat}.
#'
#' @param bioact Either a numeric vector or matrix, or a data frame providing
#'     bioactivity data.  If a numeric vector, then it is assumed that each
#'     entry corresponds to a particular fraction.  If the data is
#'     2-dimensional, then it is assumed that each column corresponds to a
#'     particular fraction, and that each row corresponds to a particular
#'     bioactivity replicate.
#'
#' @param region_ms Either \code{NULL}, or a vector either of mode character or
#'     mode numeric providing information specifying which fractions from the
#'     mass spectrometry abundances data are to be included in the data
#'     analysis.  If \code{NULL}, then it is assumed that the entirety of the
#'     mass spectrometry abundances data encapsulated in the argument to
#'     \code{msObj} is to be included in the analysis.  If numeric then the
#'     entries should provide the indices for the region of interest in the mass
#'     spectrometry data (i.e. the indices of the columns corresponding to the
#'     appropriate fractions in the data).  If character then the entries should
#'     uniquely specify the region of interest through partial string matching
#'     (i.e. the names of the columns corresponding to the appropriate fractions
#'     in the data).  The methods \code{dim}, \code{dimnames}, and
#'     \code{colnamesMS} can be used as interfaces to the mass spectrometry data
#'     encapsulated in \code{msObj}.
#'
#' @param region_bio Either \code{NULL}, or a vector either of mode character or
#'     mode numeric providing information specifying which fractions from the
#'     bioactivity data are to be included in the data analysis.  If
#'     \code{NULL}, then it is assumed that the entirety of bioactivity data
#'     provided as the argument to \code{bioact} is to be included in the
#'     analysis.  If numeric then the entries should provide the indices for the
#'     region of interest in the bioactivity data (i.e. the indices of the
#'     columns corresponding to the appropriate fractions in the data).  If
#'     character then the entries should uniquely specify the region of interest
#'     through partial string matching (i.e. the names of the columns
#'     corresponding to the appropriate fractions in the data).
#'
#' @param lambda A single nonnegative numeric value providing the quadratic
#'     penalty mixture parameter argument for the elastic net model.  The
#'     elastic net fits the least squares model with penalty function
#'     \deqn{\gamma|\beta|_1 + \lambda|\beta|^2} where \eqn{\beta} is the vector
#'     of regression coefficients and \eqn{\gamma, \lambda \ge 0}.
#'     \code{rankEN} constructs a list of candidate compounds by tracking the
#'     entrance of compounds into the elastic net model as \eqn{\gamma} is
#'     decreased from \eqn{\infty} to \eqn{0}.
#'
#' @param pos_only Either \code{TRUE} or \code{FALSE}; specifies whether the
#'     list of candidate compounds that the algorithm produces should include
#'     only those compounds that are positively correlated with bioactivity
#'     levels, or conversely should include all compounds.  The correlation is
#'     calculated using only observations from the region of interest, and when
#'     bioactivity replicates are present, the within-fraction replicates are
#'     averaged prior to calculation.
#'
#' @param ncomp Either \code{NULL}, or a numeric value no less than 1 specifying
#'     the maximum number of candidate compounds that the function should
#'     report.  When \code{NULL}, this is taken to mean that all compounds that
#'     enter the model should be reported, possibly after removing compounds
#'     nonpositively correlated with bioactivity levels, as specified by
#'     \code{pos_only}.
#'
#' @details \code{rankEN} prepares the data by extracting the region of interest
#'     from the mass spectrometry abundance data and from the bioactivity data.
#'     If bioactivity replicates are present, then the within-fraction
#'     replicates are averaged.  Once the data has been converted into the
#'     appropriate form, then an elastic net model is fitted by invoking the
#'     \code{enet} function from the \code{elasticnet} package, and an ordered
#'     list of candidate compounds is constructed such that compounds are ranked
#'     by the order in which they first enter the model.  The list may be
#'     filtered and / or pruned before being returned to the user, as determined
#'     by the arguments to \code{pos_only} and \code{ncomp}.
#'
#' @return Returns an object of class \code{rankEN}.  This object is a
#'     \code{list} with elements described below.  The class is equipped with a
#'     \code{print}, \code{summary}, and \code{extract_candidates} function.
#'
#'     \describe{
#'
#'     \item{\code{mtoz}}{ A vector providing the mass-to-charge values of the
#'     candidate compounds, such that the \code{k}-th element of the vector
#'     provides the mass-to-charge value of the \code{k}-th compound to enter
#'     the elastic net model, possibly after removing compounds nonpositively
#'     correlated with bioactivity levels. }
#'
#'     \item{\code{charge}}{ A vector providing the charge state of the
#'         candidate compounds, such that the \code{k}-th element of the vector
#'         provides the charge state of the \code{k}-th compound to enter the
#'         elastic net model, possibly after removing compounds nonpositively
#'         correlated with bioactivity levels. }
#'
#'     \item{\code{comp_cor}}{ A vector providing the correlation between each
#'         of the candidate compounds and the bioactivity levels, such that the
#'         \code{k}-th element of the vector provides the correlation between
#'         the \code{k}-th compound to enter the elastic net model and the
#'         bioactivity levels, possibly after removing compounds nonpositively
#'         correlated with bioactivity levels. }
#'
#'     \item{\code{enet_fit}}{ The fitted model object produced by
#'         \code{rankEN}'s internal invokation of the \code{enet} function from
#'         the \code{elasticnet} package.}
#'
#'     \item{\code{summ_info}}{ A list containing information related to the
#'         data used to fit the elastic net model; used by the summary
#'         function. }
#'
#'     }
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

    # Extract mtoz and charge info
    mtoz <- msDatObj$mtoz
    chg <- msDatObj$chg

    # Check for missing and that dimensions match
    rankEN_check_regr_args(ms, bio)

    # Check that compounds are nonconstant.  If any are constant then remove
    # them from the data.
    is_const <- (apply(ms, 1, stats::var) < .Machine$double.eps)
    if (sum(is_const) > 0) {
        rem_idx <- which(is_const)
        ms <- ms[-rem_idx, ]
        mtoz <- mtoz[-rem_idx]
        chg <- chg[-rem_idx]
    }

    # Convert ms to form where rows are an observation (i.e. fraction) and cols
    # are a variable (i.e. a compound)
    ms_regr <- t(ms)

    # Obtain the mean of the bioactivity replicates and convert to a vector
    bio_regr <- colMeans(bio)


    # Fit model ------------------------------------------------------------------

    # Calculate the elastic net path
    enet_fit <- tryCatch({
        elasticnet::enet(ms_regr, bio_regr, lambda)
    }, warning = function(war) {
        warning("message produced by call to enet from package elasticnet ==>\n",
                war[[1L]], call.=FALSE)
    }, error = function(err) {
        stop("message produced by call to elasticnet::enet ==>\n",
             err[[1L]], call.=FALSE)
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
        data_dim  = list(reg  = ncol(bio),
                         comp = nrow(msDatObj),
                         repl = nrow(bio)),
        region_nm = list(ms  = ms_nm,
                         bio = bio_nm),
        lambda    = lambda,
        pos_only  = pos_only,
        ncomp     = ncomp
    )
    if (sum(is_const) > 0) {
        summ_info$cmp_rm <- list(mtoz = msDatObj$mtoz[rem_idx],
                                 chg  = msDatObj$chg[rem_idx])
    }

    # Construct output object
    outDat <- list(mtoz      = mtoz[comp_idx_out],
                   charge    = chg[comp_idx_out],
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
#'     column should be included in the returning \code{data.frame} providing
#'     the correlation between the compound and the within-fraction average of
#'     the bioactivity replicates in the region of interest.
#'
#' @export


extract_candidates <- function(rankEN_obj, include_cor=TRUE) {

    if (!identical(class(rankEN_obj), "rankEN")) {
        stop("rankEN_obj must be of class rankEN")
    }
    else if (!identical(include_cor, TRUE) && !identical(include_cor, FALSE)) {
        stop("include_cor must be either TRUE or FALSE")
    }

    out <- data.frame(mtoz=rankEN_obj$mtoz, charge=rankEN_obj$charge)
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
        "\n",
        "An object of class \"rankEN\".\n",
        "Use summary to print a list of the compounds coefs along the elastic net path.\n",
        "Use extract_candidates to extract the compound info as a data.frame.\n\n")

    dd_char <- format_int( unlist(x$summ_info$data_dim) )
    cat(sep="",
        "Data dimensions:\n",
        "----------------\n",
        "    region of interest:     ", dd_char[1], "\n",
        "    candidate compounds:    ", dd_char[2], "\n",
        "    bioactivity replicates: ", dd_char[3], "\n",
        "\n")
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
    cat(format(object, max_comp_print), sep="")
}




format.rankEN <- function(x, max_comp_print, ...) {


    # Check argument to max_comp_print
    if (!is.numeric(max_comp_print)) {
        stop("max_comp_print must be of mode numeric")
    }
    else if ( !identical(length(max_comp_print), 1L) ) {
        stop("max_comp_print must have length 1")
    }
    else if (max_comp_print < 1L) {
        stop("max_comp_print cannot have value less than 1")
    }

    # Create links for convenience
    summ_info <- x$summ_info
    mtoz      <- x$mtoz
    charge    <- x$charge
    comp_cor  <- x$comp_cor

    # Number of compounds we should print
    print_len <- min(as.integer(max_comp_print), length(mtoz))


    # Begin variable construction --------------------------

    # Compound removal notification
    if (! is.null(summ_info$cmp_rm)) {
        rm_notif     <- " (some compounds were removed)"
        rm_notif_bar <- "------------------------------"

        comp_rm_header <- paste0(
            "Compounds removed for being constant:\n",
            "-------------------------------------\n")
        ms_comp_rm <- format(
            c("Mass-to-charge",
              "--------------",
              format_float(summ_info$cmp_rm$mtoz)))
        chg_comp_rm <- format(
            c("Charge",
              "------",
              format_float(summ_info$cmp_rm$chg)))
    }
    else {
        rm_notif <- rm_notif_bar <- ""
    }

    # Mass spectrometry and bioactivity data dimensions
    data_dim_vals <- format_int( unlist(summ_info$data_dim) )
    names(data_dim_vals) <- c("nroi", "ncomp", "nbio")

    # Column of mass spectrometry fraction names in region of interest
    ms_roi_nm  <- format(
        c("Mass spec.",
          "----------",
          summ_info$region_nm$ms))

    # Column of bioactivity fraction names in region of interest
    bio_rio_nm <- format(
        c("Bioactivity",
          "-----------",
          summ_info$region_nm$bio))

    # Parameter arguments to rankEN.  Note that we don't need nsmall as an
    # argument to the inner format call because since there is only one number
    # we do in fact want to ignore trailing 0's in this case.
    param_arg_vals <- c(
        format(summ_info$lambda, big.mark=","),
        ifelse(summ_info$pos_only, "yes", "no"),
        ifelse(is.null(summ_info$ncomp), "all", format(summ_info$ncomp, big.mark=",")))
    names(param_arg_vals) <- c("lambda", "pos_only", "ncomp")

    # Ordered compound rankings m/z values
    comp_rank_mtoz <- format(
        c("Mass-to-charge",
          "--------------",
          format_float( mtoz[1:print_len] )))

    # Ordered charge rankings charge values
    comp_rank_chg <- format(
        c("Charge",
          "------",
          format_int( charge[1:print_len] )))

    # Ordered charge rankings correlation values
    comp_rank_corr <- format(
        c("Correlation",
          "-----------",
          format(round(comp_cor[1:print_len], 4), nsmall=4)))

    # Compound rankings header
    if (length(mtoz) <= max_comp_print) {
        comp_rank_header <-
            paste0("Compounds in order of entrance (all compounds, earliest at top):\n",
                   "----------------------------------------------------------------\n")
    } else {
        header <-
            sprintf("Compounds in order of entrance (first %s compounds, earliest at top):\n",
                    formatC(max_comp_print, format="d", big.mark=","))
        comp_rank_header <- paste0(header, paste0(rep("-", nchar(header) - 1), collapse=""), "\n")
    }


    # Begin string construction ----------------------------

    # Dimensions of the data supplied to rankEN
    data_dim <- sprintf(
        paste0("Data dimensions%s:\n",
               "---------------%s-\n",
               "    region of interest:      %s\n",
               "    candidate compounds:     %s\n",
               "    bioactivity replicates:  %s\n",
               "\n"),
        rm_notif,
        rm_notif_bar,
        data_dim_vals["nroi"],
        data_dim_vals["ncomp"],
        data_dim_vals["nbio"])

    # Names of compounds removed for being constant in roi
    if (! is.null(summ_info$cmp_rm)) {
        compound_removal <- paste0(
            comp_rm_header,
            paste0("    ", ms_comp_rm, "    ", chg_comp_rm, "\n", collapse=""),
            "\n")
    }
    else {
        compound_removal <- ""
    }

    # Names of fractions included in region of interest
    region_of_interest <- paste0(
        "Fractions included in region of interest:\n",
        "-----------------------------------------\n",
        paste0("    ", ms_roi_nm, "    ", bio_rio_nm, "\n", collapse=""),
        "\n")

    # Parameter arguments supplied to rankEN
    param_args <- sprintf(
        paste0("Parameter arguments provided to rankEN:\n",
               "---------------------------------------\n",
               "    Quadratic penalty parameter:          %s\n",
               "    Consider only positive correlations:  %s\n",
               "    Max number of candidate compounds:    %s\n",
               "\n"),
        param_arg_vals["lambda"],
        param_arg_vals["pos_only"],
        param_arg_vals["ncomp"])

    # Compound rankings header
    comp_rank <- paste0(
        comp_rank_header,
        paste0("    ", comp_rank_mtoz, "    ", comp_rank_chg,
               "    ", comp_rank_corr, "\n", collapse=""),
        "\n")

    c(newl = "\n",
      ddim = data_dim,
      crem = compound_removal,
      regi = region_of_interest,
      args = param_args,
      rank = comp_rank)
}
