
#' Filter compounds from mass spectrometry data
#'
#' Filters mass spectrometry data using a set of criteria, described in
#' \code{Details}. Returns an object of classes \code{\link{msDat}} and
#' \code{filterMS}.
#'
#' @param msObj An object class \code{\link{msDat}}.  Note that this includes
#'     objects created by the functions \code{binMS} and \code{msDat}.
#'
#' @param region A vector either of mode character or mode numeric.  If numeric
#'     then the entries should provide the indices for the region of interest in
#'     the mass spectrometry data provided as the argument for \code{msObj}.  If
#'     character then the entries should uniquely specify the region of interest
#'     through partial string matching (see criterion 1, 4).
#'
#' @param border Either a character string \code{"all"}, or a character string
#'     \code{"none"}, or a length-1 or length-2 numeric value specifying the
#'     number of fractions to either side of the region of interest to comprise
#'     the bordering region.  If a single numeric value, then this is the number
#'     of fractions to each side of the region of interest; if it is two values,
#'     then the first value is the number of fractions to the left, and the
#'     second value is the number of fractions to the right.  If there are not
#'     enough fractions in either direction to completely span the number of
#'     specified fractions, then all of the available fractions to the side in
#'     question are considered to be part of the bordering region (see criterion
#'     2).
#'
#' @param bord_ratio A single nonnegative numeric value.  A value of 0 will not
#'     admit any compounds, while a value greater than 1 will admit all
#'     compounds (see criterion 2).
#'
#' @param min_inten A single numeric value.  A value less than the minimum mass
#'     spectrometry value in the data will admit all compounds (see criterion
#'     4).
#'
#' @param max_chg A single numeric value specifying the maximum charge which a
#'     compound may exhibit (see criterion 5)
#'
#' @details Attempts to filter out candidate compounds via subject-matter
#'     knowledge, with the goal of removing spurious noise from downstream
#'     models.  The criteria for the downstream inclusion of a candidate
#'     compound is listed below.
#'
#'     \enumerate{
#'
#'     \item The m/z intensity maximum must fall inside the range of the
#'         bioactivity region of interest
#'
#'     \item The ratio of the m/z intensity of a species in the areas bordering
#'         the region of interest and the species maximum intensity must be less
#'         than \code{bord_ratio}.  When there is no bordering area then it is
#'         taken to mean that all observations satisfy this criterion.
#'
#'     \item The immediately right adjacent fraction to its maximum intensity
#'         fraction for a species must have a non-zero abundance.  In the case
#'         of ties for the maximum, it is the fraction immediately to the right
#'         of the rightmost maximum fraction which cannot have zero abundance.
#'         When the fraction with maximum intensity is the rightmost fraction in
#'         the data for an observation, then it is taken to mean that the
#'         observation satisfies this criterion.
#'
#'     \item At least 1 fraction in the region of interest must have intensity
#'         greater than \code{min_inten}
#'
#'     \item Compound charge state must be less than or equal to \code{max_chg}
#'
#'     }
#'
#' @return Returns an object of class \code{filterMS} which inherits from
#'     \code{msDat}.  This object is a \code{list} with elements described
#'     below.  The class is equipped with a \code{print}, \code{summary}, and
#'     \code{extractMS} function.
#'
#'     \describe{
#'
#'     \item{\code{msDatObj}}{ An object of class \code{\link{msDat}} such that
#'         the encapsulated mass spectrometry data corresponds to each of the
#'         candidate compounds that satisfed each of the criteria.  If no
#'         criteria are satisfied then \code{NULL} is returned. }
#'
#'     \item{\code{cmp_by_crit}}{ A list containing \code{data.frame}s, one for
#'         each criterion. Each row (if any) in one of the
#'         sub-\code{data.frame}s contains the mass-to-charge and charge
#'         information for a candidate compound that satisfies the criterion
#'         represented by the \code{data.frame}; all of the compounds that
#'         satisfied the criterion are included in the data.  The
#'         \code{data.frame}s are named \code{c1}, ..., \code{c5}, etc
#'         corresponding to criterion 1, ..., criterion 5. }
#'
#'     \item{\code{summ_info}}{ A list containing information pertaining to the
#'          filtering process; for use by the summary function. }
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
#' # Filter out potential candidate compounds
#' filter_out <- filterMS(msObj = ms,
#'                        region = paste0("VO_", 17:25),
#'                        border = "all",
#'                        bord_ratio = 0.01,
#'                        min_inten = 1000,
#'                        max_chg = 7)
#'
#' # print, summary function
#' filter_out
#' summary(filter_out)
#'
#' # Extract filtered mass spectrometry data as a matrix or msDat object
#' filter_matr <- extractMS(msObj = filter_out, type = "matrix")
#' filter_msDat <- extractMS(msObj = filter_out, type = "matrix")
#'
#' @export


filterMS <- function(msObj, region, border="all", bord_ratio=0.05, min_inten=1000, max_chg=7L) {

    # Check validity of arguments
    filterMS_check_valid(msObj, region, border, bord_ratio, min_inten, max_chg)

    # Number of criterion
    nCrit <- 5L

    # Create pointers to mass spec variables for convenience
    msDatObj <- extractMS(msObj, "msDat")
    ms <- msDatObj$ms
    mtoz <- msDatObj$mtoz
    chg <- msDatObj$chg

    # Dimensions of mass spectrometry data
    ms_nr <- NROW(ms)
    ms_nc <- NCOL(ms)

    # Create region index variable
    regIdx <- extract_idx(ms, region, TRUE)

    # Create border index, i.e. the indices that surround the region of interest
    borIdx <- filterMS_border_idx(border, regIdx, ms_nc)

    # All indices in either region or border
    allIdx <- union(regIdx, borIdx)
    lastIdx <- max(allIdx)

    # maxIdx: the column index per row (i.e. per observation) of the fraction
    # containing the maximum intensity level.  The rightmost column index is
    # chosen in the case of ties.  Note that we only consider columns in either
    # the region of interest or the bordering region.
    maxIdx <- apply(ms[, allIdx, drop=FALSE], 1, function(x) {
        last_in_allIdx <- utils::tail(which(x == max(x)), 1L)
        allIdx[last_in_allIdx]
    })

    # Evaluate criteria predicates.  Note that all() returns true when there are
    # no values as happens for criteria 2 when borIdx is integer(0), which is
    # the desired behavior.
    critBool <- data.frame( array(dim=c(ms_nr, nCrit)) )
    row_seq <- seq_len(ms_nr)
    critBool[, 1L] <- maxIdx %in% regIdx
    critBool[, 2L] <- sapply(row_seq, function(i) all(ms[i, borIdx] < bord_ratio * ms[i, maxIdx[i]]))
    critBool[, 3L] <- sapply(row_seq, function(i) (ms[i, min(maxIdx[i] + 1L, lastIdx)] > 0))
    critBool[, 4L] <- sapply(row_seq, function(i) any(ms[i, regIdx] > min_inten))
    critBool[, 5L] <- (chg <= max_chg)

    # Create a vector of indices which satisfy every criterion
    keepIdx <- which( Reduce("&", critBool) )

    # Extract mass spec fraction names
    ms_nm <- colnames(ms)
    if (is.null(ms_nm)) {
        ms_nm <- as.character(seq_len(ms_nc))
    }

    # Create filtered mass spectrometry data
    if (length(keepIdx) > 0) {
        msDatObj <- msDat(ms[keepIdx, , drop=FALSE], mtoz[keepIdx], chg[keepIdx])
    }
    else {
        msDatObj <- NULL
        warning("There are no compounds that met all of the criteria\n", call.=FALSE)
    }

    # Create mass-to-charge and charge datasets for each criterion
    cmp_by_cr <- stats::setNames(vector("list", nCrit), paste0("c", seq_len(nCrit)))
    for (j in seq_len(nCrit)) {
        thisKeep <- critBool[, j]
        # note: data.frame can handle the case when thisKeep has length 0
        cmp_by_cr[[j]] <- data.frame(mtoz = mtoz[thisKeep],
                                     chg  = chg[thisKeep] )
    }

    # Construct return object
    outObj <- list( msDatObj  = msDatObj,
                   cmp_by_cr = cmp_by_cr,
                   summ_info = list(orig_dim   = c(ms_nr, ms_nc),
                                    reg_nm     = ms_nm[regIdx],
                                    bor_nm     = ms_nm[borIdx],
                                    border     = border,
                                    bord_ratio = bord_ratio,
                                    min_inten  = min_inten,
                                    max_chg    = max_chg ) )

    structure(outObj, class=c("filterMS", "msDat"))
}




#' Basic information for class \code{filterMS}
#'
#' Displays the number of candidate compounds left in the data after filtering
#'
#' @param x An object of class \code{\link{filterMS}}
#'
#' @param ... Arguments passed to dot-dot-dot are ignored
#'
#' @export


print.filterMS <- function(x, ...) {

    if (is.null(x$msDatObj)) {
        cat("An object of class \"filterMS\"; no observations ",
            "satisfied all of the inclusion criteria.\n", sep="")
    }
    else {
        cat("An object of class \"filterMS\" with ", format(NROW(x$msDatObj$ms), big.mark=","),
            " compounds and ", NCOL(x$msDatObj$ms), " fractions.\n", sep="")
    }
    cat("Use summary to see more details regarding the filtering process.\n",
        "Use extractMS to extract the filtered mass spectrometry data\n\n", sep="")
}




#' Overview of the filtering process
#'
#' Prints a description of the filtering process.  Displays arguments chosen for
#' the \code{filterMS} constructor, how many candidate compounds were chosen for
#' each criterion, and how many candidate compounds were chosen overall.
#'
#' @param object An object of class \code{\link{filterMS}}
#'
#' @param ... Arguments passed to dot-dot-dot are ignored
#'
#' @export


summary.filterMS <- function(object, ...) {
    cat(format(object), sep="")
}




# Return a vector of strings supplying the output for summary.filterMS

format.filterMS <- function(x, ...) {

    # Add pointers to summ_info variables for convenience
    orig_dim   <- x$summ_info$orig_dim
    reg_nm     <- x$summ_info$reg_nm
    bor_nm     <- x$summ_info$bor_nm
    border     <- x$summ_info$border
    bord_ratio <- x$summ_info$bord_ratio
    min_inten  <- x$summ_info$min_inten
    max_chg    <- x$summ_info$max_chg
    ncmp_by_cr <- sapply(x$cmp_by_cr, nrow)
    msDatObj   <- x$msDatObj

    # A bar (i.e. ----) to place underneath the region of interest header
    roi_bar <- paste0( rep("-", 53 + nchar(length(reg_nm))), collapse="" )

    # The region of interest names concatenated together
    roi_nm_cat <- paste0("    ", reg_nm, "\n", collapse="")

    # Bordering region specification
    if (identical(border, "all")) {
        border_spec <- "\"all\""
    } else if (identical(border, "none")) {
        border_spec <- "\"none\""
    } else if (identical(length(border), 1L)) {
        border_spec <- paste0("each having length ", border, ":")
    } else {
        border_spec <- paste0("having lengths ", border[1L], " and ", border[2L], ":")
    }

    # A bar (i.e. ----) to place underneath the bordering region header
    bor_bar <- paste0( rep("-", 40 + nchar(border_spec)), collapse="" )

    # The bordering region names concatenated together
    if (identical(length(bor_nm), 0L)) {
        bor_nm_cat <- "    * no fraction names to show *\n"
    } else if (length(bor_nm) > 10) {
        bor_nm_cat <- "    * fraction names omitted for brevity *\n"
    } else {
        bor_nm_cat <- paste0("    ", bor_nm, "\n", collapse="")
    }

    # Filtering criteria strings with uniform width
    filt_crit <- format_float( c(min_inten, max_chg, bord_ratio) )
    names(filt_crit) <- c("min_inten", "max_chg", "bord_ratio")

    # Prior dimension strings with uniform width
    orig_dim_str <- format_int(orig_dim)
    names(orig_dim_str) <- c("compounds", "fractions")

    # Individual filtering criterion remaining m/z levels
    nlevels_fmt <- format(orig_dim[1L], big.mark=",")
    ncrit_remain <- format_int(ncmp_by_cr)
    names(ncrit_remain) <- paste0("crit", 1:5)

    # A bar (i.e. ----) to place underneath the ind. crit. section
    ncrit_bar <- paste0( rep("-", 77 + nchar(nlevels_fmt)), collapse="" )

    # Final number of remaining levels after all filtering
    nfinal <- ifelse(is.null(msDatObj), 0, NROW(msDatObj))
    nfinal_str <- formatC(nfinal, format="d", big.mark=",")


    # Begin string construction ----------------------------

    # The region of interest column names
    region_of_interest <- sprintf(
        paste0("The region of interest was specified as (%d fractions):\n",
               "%s\n",
               "%s\n"),
        length(reg_nm),
        roi_bar,
        roi_nm_cat)

    # The bordering region column names
    bordering_region <- sprintf(
        paste0("The bordering regions were specified as %s\n",
               "%s\n",
               "%s",
               "\n"),
        border_spec,
        bor_bar,
        bor_nm_cat)

    # The filtering criteria specification
    filtering_criteria <- sprintf(
        paste0("The filtering criteria was specified as:\n",
               "----------------------------------------\n",
               "    minimum intensity:      %s\n",
               "    maximum charge:         %s\n",
               "    bordering region ratio: %s\n",
               "\n"),
        filt_crit["min_inten"],
        filt_crit["max_chg"],
        filt_crit["bord_ratio"])

    # Dimension of data prior to filtering
    prior_dimen <- sprintf(
        paste0("The mass spectrometry data prior to filtering had:\n",
               "--------------------------------------------------\n",
               "    %s compounds\n",
               "    %s fractions\n",
               "\n"),
        orig_dim_str["compounds"],
        orig_dim_str["fractions"])

    # Remaining compounds after individual criterion filtering
    filtering_remaining <- sprintf( paste0(
        "Individually, each criterion reduced the %s m/z levels to the following number:\n",
        "%s\n",
        "    criterion 1:  %s    (fraction with max. abundance is in region of interest)\n",
        "    criterion 2:  %s    (fractions in bordering region have < %s%% of max. abundance)\n",
        "    criterion 3:  %s    (nonzero abundance in right adjacent fraction to max.)\n",
        "    criterion 4:  %s    (at least 1 intensity > %s in region of interest)\n",
        "    criterion 5:  %s    (must have charge <= %d)\n",
        "\n"),
        nlevels_fmt,
        ncrit_bar,
        ncrit_remain["crit1"],
        ncrit_remain["crit2"],
        format(round(100 * bord_ratio), big.mark=","),
        ncrit_remain["crit3"],
        ncrit_remain["crit4"],
        format(min_inten, big.mark=","),
        ncrit_remain["crit5"],
        max_chg)

    # Number of compounds that satisfy all filtering criterion
    filtering_final_amt <- sprintf(
        paste0("The total number of candidate compounds was reduced to:\n",
               "-------------------------------------------------------\n",
               "    %s\n",
               "\n"),
        nfinal_str)

    c(newl = "\n",
      regn = region_of_interest,
      bord = bordering_region,
      fcri = filtering_criteria,
      prir = prior_dimen,
      frem = filtering_remaining,
      ffin = filtering_final_amt)
}
