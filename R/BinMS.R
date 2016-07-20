
#' Consolidate mass spectrometry observations
#'
#' Combines mass spectrometry observations that are believed to belong to the
#' same underlying compound into a single observation.  In concept, the data
#' produced by the mass spectrometer may produce multiple reads for a single
#' compound; thus, \code{binMS} attempts to recover these underlying compounds
#' through a binning procedure, described in more detail in \code{Details}.
#'
#' @inheritParams msDat
#'
#' @param mass The information for the mass need not be provided, as it can be
#'     derived using the mass-to-charge and charge information; in this case the
#'     parameter should be given its default, i.e. \code{NULL}.  If however the
#'     information for mass is already included in the dataset in hand, then
#'     providing it to the function will be slightly more efficient then
#'     re-performing the calculations.  The information for the \code{charge}
#'     parameter can be provided in the same manner as for the mass-to-charge
#'     values.
#'
#' @param time_peak_reten The information for the \code{time_peak_reten}
#'     parameter can be provided in the same manner as for the mass-to-charge
#'     and other information; this paramater specifies the time at which the
#'     peak retention level of the compound was achieved.
#'
#' @param time_range A length-2 numeric vector specifying the lower bound and
#'     upper bound (inclusive) of allowed peak retention time occurance for an
#'     observation to be included in the consolidation process.
#'
#' @param mass_range A length-2 numeric vector specifying the lower bound and
#'     upper bound (inclusive) of allowed mass for an observation to be included
#'     in the consolidation process.
#'
#' @param charge_range A length-2 numeric vector specifying the lower bound and
#'     upper bound (inclusive) of allowed electrical charge state for an
#'     observation to be included in the consolidation process.
#'
#' @param mtoz_diff A single numerical value such that any two observations with
#'     a larger absolute difference between their mass-to-charge values are
#'     considered to have originated from different underlying compounds.  Two
#'     observations with a smaller absolute difference between their
#'     mass-to-charge values could potentially be considered to originate from
#'     the same underlying compound, contingent on other criteria also being
#'     met.  Nonnegative values are allowed; such a value has the effect of not
#'     consolidating any groups, and consequently reduces the function to a
#'     filtering routine only.
#'
#' @param time_diff A single numerical value such that any two observations with
#'     a larger absolute difference between their peak elution times are
#'     considered to have originated from different underlying compounds.  Two
#'     observations with a smaller absolute difference between their peak
#'     elution times could potentially be considered to originate from the same
#'     underlying compound, contingent on other criteria also being met.
#'     Nonnegative values are allowed; such a value has the effect of not
#'     consolidating any groups, and consequently reduces the function to a
#'     filtering routine only.
#'
#' @details The algorithm described in what follows attempts to combines mass
#'     spectrometry observations that are believed to belong to the same
#'     underlying compound into a single observation for each compound.  There
#'     are two conceptually separate steps.
#'
#'     The first step is as follows.  All observations must satisfy each of the
#'     following criteria for inclusion in the binning process.
#'
#'     \enumerate{
#'
#'     \item Each observation must have its peak elution time occur during the
#'         interval specified by \code{time_range}
#'
#'     \item Each observation must have a mass that falls within the interval
#'         specified by \code{mass_range}
#'
#'     \item Each observation must have an electrical charge state that falls
#'         within the interval specified by \code{charge_range}
#'
#'     }
#'
#'     Once that a set of observations satisfying the above criteria is
#'     obtained, then a second step attempts to combine observations believed to
#'     belong to the same underlying compound.  The algorithm considers two
#'     observations that satisfy each of the following criteria to belong to the
#'     same compound.
#'
#'     \enumerate{
#'
#'     \item The absolute difference in Daltons of the mass-to-charge value
#'         between the two observations is less the the value specified by
#'         \code{mtoz_diff}
#'
#'      \item The absolute difference of the peak elution time between the two
#'          observations is less than the value specified by \code{time_pr_diff}
#'
#'     \item The electrical charge state must be the same for the two observations
#'
#'     }
#'
#'     Then the binning algorithm is defined as follows.  Consider an
#'     observation that satisfies the inclusion criteria; this observation is
#'     compaired pairwise with every other observation that satisfies the
#'     inclusion criteria.  If a pair of observations satisfies the criteria
#'     determining them to belong to the same underlying compound then the two
#'     observations are merged into a single observation.  The two previous
#'     compounds are removed from the working set, and the process starts over
#'     with the newly created observation.  The process repeats until no other
#'     observation in the working set meets the criteria determining it to
#'     belong to the same underlying compound as that of the current
#'     observation; at this point it is considered that all observations
#'     belonging to the compound have been found, and the process starts over
#'     with a new observation.
#'
#'     The merging process has not yet been defined; it is performed by
#'     averaging the mass-to-charge values and peak elution times, and summing
#'     the mass spectrometry intensities at each fraction.  Although
#'     observations are merged pairwise, when multiple observations are combined
#'     in a sequence of pairings, the averages are given equal weight for all of
#'     the observations.  In other words, if a pair of observations are merged,
#'     and then a third observation is merged with the new observation created
#'     by combining the original two, then the mass-to-charge value and peak
#'     elution time values of the new observation are obtained by summing the
#'     values for each of the three original observations and dividing by three.
#'     The merging process for more than three observations is conducted
#'     similarly.
#'
#'     Having described the binning algorithm, it is apparent that there are
#'     scenarios in which the order in which observations are merged affects the
#'     outcome of the algorithm.  Since it seems that a minumum requirement of
#'     any binning algorithm is that the algorithm is invariant to the ordering
#'     of the observations in the data, this algorithm abides by the following
#'     rules.  The observations in the data are sorted in increasing order by
#'     mass-to-charge value, peak elution time, and electical charge state,
#'     respectively.  Then when choosing an observation to compare to the rest
#'     of the set, we start with the observation at the top of the sort
#'     ordering, and compare it one-at-a-time to the other elements in the set
#'     according to the same ordering.  When a consolidated observation is
#'     complete in that no other observation left in the working set satisfies
#'     the merging criteria, then this consolidated observation can be removed
#'     from consideration for all future merges.

#' @return Returns an object of class \code{binMS} which inherits from
#'     \code{msDat}.  This object is a \code{list} with elements described
#'     below.  The class is equipped with a \code{print}, \code{summary}, and
#'     \code{extractMS} function.
#'
#'     \describe{
#'
#'     \item{\code{msDatObj}}{ An object of class \code{\link{msDat}} that
#'         encapsulates the mass spectrometry data for the consolidated data. }
#'
#'     \item{\code{summ_info}}{ A list containing information pertaining to the
#'         consolidation process; for use by the summary function. }
#'
#'     }
#'
#' @export


binMS <- function(mass_spec,
                  mtoz,
                  charge,
                  mass = NULL,
                  time_peak_reten,
                  ms_inten = NULL,
                  time_range,
                  mass_range,
                  charge_range,
                  mtoz_diff,
                  time_diff) {

    # Check for validity of the form of the arguments
    binMS_check_valid_input(mass_spec, mtoz, charge, mass, time_peak_reten, ms_inten, time_range,
                            mass_range, charge_range, mtoz_diff, time_diff)

    # Translate arguments to values
    dmtoz <- extract_var(mass_spec, mtoz)
    dcharge <- extract_var(mass_spec, charge)
    dtime_pr <- extract_var(mass_spec, time_peak_reten)
    # Check if we need to calculate mass from scratch
    if (is.null(mass)) {
        dmass <- charge * (dmtoz - 1.007825)
        # Can't pass NULL to extract_var so create a variable
        mass <- dmass
    } else {
        dmass <- extract_var(mass_spec, mass)
    }
    dms_inten <- extract_var(mass_spec, ms_inten, TRUE, mtoz, charge, mass, time_peak_reten)

    ## Step 1: construct sorted indices of rows in the data that satisfy the
    ## mass, time of peak retention, and charge criteria

    # Calculate inclusion criteria status
    time_pr_bool <- (time_range[1] <= dtime_pr) & (dtime_pr <= time_range[2])
    mass_bool <- (mass_range[1] <= dmass) & (dmass <= mass_range[2])
    charge_bool <- (charge_range[1] <= dcharge) & (dcharge <= charge_range[2])
    keepIdx <- which( Reduce("&", list(time_pr_bool, mass_bool, charge_bool)) )

    # Number of observations satisfying the inclusion criteria
    n_befbin <- length(keepIdx)

    ## Step 2: reduce data by row criteria and by selecting columns needed in
    ## future.  Transpose so as to conform to column major order.

    # case: some observations satisfied the inclusion criteria
    if (!identical(n_befbin, 0L)) {

        # Remove m/z levels that didn't meet inclusion criteria
        dmtoz <- dmtoz[keepIdx]
        dcharge <- dcharge[keepIdx]
        dtime_pr <- dtime_pr[keepIdx]
        dmass <- dmass[keepIdx]
        dms_inten <- dms_inten[keepIdx, , drop=FALSE]

        # Obtain indices of ordered data
        sortIdx <- order(dmtoz, dtime_pr, dcharge)

        # ms: transpose of the mass spec data
        ms <- t( dms_inten[sortIdx, , drop=FALSE] )
        # info: data used to identify and combine the mass-to-charge levels
        info_orig <- matrix(c(dmtoz, dcharge, dtime_pr, dmass), nrow=4, byrow=TRUE)
        info_orig <- info_orig[, sortIdx]

        ## Step 3: perform binning of observations that meet the similarity
        ## criteria based on mass-to-charge values, time of peak retention, and
        ## charge state

        # Allocate memory for binned data.  We keep the data in two matrices;
        # one for what we call the information which are the values that
        # determine whether or not to combine two m/z levels, and another for
        # the mass spectrometry abundances.  The reason for separating the two
        # is that the information data gets averaged within a bin while the ms
        # abundances are summed.
        info_bin <- matrix(nrow=4, ncol=n_befbin)
        ms_bin <- matrix(nrow=nrow(ms), ncol=n_befbin)
        row.names(ms_bin) <- colnames(dms_inten)

    } # end keepIdx not empty case

    # Give names to rows (i.e. the variables) in the binned data information
    # matrix so as to make the algorithm more readable.  The variables
    # correspond to the m/z value, charge, time of peak retention, and mass
    # respectively.
    rmtoz  <- 1L
    rchg   <- 2L
    rtime  <- 3L
    rmass  <- 4L

    binCtr   <- 1L  # Index for next binned data entry
    origCtr  <- 1L  # Index for next orig data to be considered for binned data
    innerCtr <- 1L  # Index for orig data within inner loop that looks for all
                    # occurences of data within the allowed m/z difference
    nbin     <- 1L  # Number of mass-to-charge levels in current bin

    # Create an out-of-bounds charge value used to signal that a mass-to-charge
    # level has already been included as part of a bin
    flagv <- charge_range[1] - 1

    # Each iteration compares an m/z level to see if we can start a new bin.  If
    # the current level of the iteration isn't part of a bin, then we start a
    # new bin and compare to larger m/z levels, adding levels to the bin when
    # criteria are met, and ending when m/z levels gets out of range.

    while (origCtr <= n_befbin) {

        # case: current row was already included as part of a bin; move on to
        # next m/z level in the pre-binned data (recall that when a row was
        # already included as part of a bin then the charge is set to
        # charge_range[1] - 1 as a signal)
        if (info_orig[rchg, origCtr] < charge_range[1]) {
            origCtr <- origCtr + 1L
            next
        }

        # If we've made it here then we've found a m/z level in the original
        # data not yet included in any previous bin.  Begin a new bin starting
        # with this level as its first entry.
        info_bin[, binCtr] <- info_orig[, origCtr]
        ms_bin[, binCtr] <- ms[, origCtr]
        # Set innerCtr to be the index of the first m/z level to consider for
        # adding to the new bin
        innerCtr <- origCtr + 1L
        # Initial number of m/z levels in the bin
        nbin <- 1L

        # Each iteration compares a m/z value to the current m/z level until the
        # values being compared against have a m/z outside of the allowed
        # range. If all of the criteria for the iteration are met, then the m/z
        # level is combined with the current bin.
        #
        # Note that when comparisons are made that we divide by nbin (i.e. the
        # number of compounds in the bin) since we are only summing observations
        # from the original data when placing it in the bin (and taking the mean
        # occurs after the (following inner) loop) ends.

        while ( (innerCtr <= n_befbin) &&
                (info_orig[rmtoz, innerCtr] - (info_bin[rmtoz, binCtr] / nbin) < mtoz_diff) ) {

            # case: criteria met.  Add m/z level to current bin.
            if ( (info_orig[rchg, innerCtr] == (info_bin[rchg, binCtr] / nbin))
                && ( abs(info_orig[rtime, innerCtr] - (info_bin[rtime, binCtr] / nbin))
                < time_diff ) ) {

                # add m/z level to current bin
                info_bin[, binCtr] <- info_bin[, binCtr] + info_orig[, innerCtr]
                ms_bin[, binCtr] <- ms_bin[, binCtr] + ms[, innerCtr]

                # signal that m/z level is already part of a bin
                info_orig[rchg, innerCtr] <- flagv
                nbin <- nbin + 1L

            }
            # else: criteria not me; noop

            # Update the counter indexing the next m/z level to consider
            # for adding to the current bin
            innerCtr <- innerCtr + 1L

        } # end compare current bin to next m/z level loop

        # Take the average of the binned m/z levels (note: don't need to do this
        # for ms_bin as these values are strictly summed)
        info_bin[, binCtr] <- info_bin[, binCtr] / nbin

        # We've found every m/z level in the original data that belongs in the
        # binned data.  Update the counter indexing the next m/z level from the
        # original data, and update the counter indexing the next available bin
        # to start placing consolidated data into.
        origCtr <- origCtr + 1L
        binCtr <- binCtr + 1L

    } # end iterate over original m/z levels loop

    # Construct msDat object from binned data
    nbinned <- binCtr - 1L
    if (!identical(nbinned, 0L)) {

        # It is possible that the binned data can be out of order.  This happens
        # when the minimum observation of one bin is smaller than that of
        # another, but the mean is larger.  binCtr indexes the next available
        # bin, so 1 less is the number of bins.
        resortIdx <- order(info_bin[rmtoz, 1:nbinned],
                           info_bin[rtime, 1:nbinned],
                           info_bin[rchg, 1:nbinned])

        # Construct msDat object
        msDatObj <- msDat(t(ms_bin[, resortIdx,  drop=FALSE]),
                          info_bin[rmtoz, resortIdx],
                          info_bin[rchg, resortIdx])
    }
    else {
        warning("No observations satisfied all of the inclusion criteria", call.=FALSE)
        msDatObj <- NULL
    }

    # Info for use by summary.binMS to describe binning process
    summ_info <- list(n_tot        = nrow(mass_spec),
                      n_time_pr    = sum(time_pr_bool),
                      n_mass       = sum(mass_bool),
                      n_charge     = sum(charge_bool),
                      n_tiMaCh     = length(keepIdx),
                      n_binned     = nbinned,
                      time_range   = time_range,
                      mass_range   = mass_range,
                      charge_range = charge_range,
                      mtoz_diff    = mtoz_diff,
                      time_diff    = time_diff)

    # Construct binMS object
    outObj <- list(msDatObj  = msDatObj,
                   summ_info = summ_info)

    structure(outObj, class=c("binMS", "msDat"))
}




#' Print routine for class \code{binMS}
#'
#' Prints the number of m/z levels and fractions of the resultant mass
#' spectrometry data
#'
#' @param x An object of class \code{\link{binMS}}
#'
#' @param ... Arguments passed to dot-dot-dot are ignored
#'
#' @export

print.binMS <- function(x, ...) {

    if (is.null(x$msDatObj)) {
        cat("An object of class \"binMS\"; no observations satisfied all of the inclusion criteria.\n")
    }
    else {
        cat("An object of class \"binMS\" with ", NROW(x$msDatObj$ms), " compounds and ",
            NCOL(x$msDatObj$ms), " fractions.\n", sep="")
    }
    cat("Use summary to see more details regarding the consolidation process.\n",
        "Use extractMS to extract the consolidated mass spectrometry data.\n\n", sep="")
}




#' Overview of the binning process
#'
#' Prints a text description of the binning process.  Displays arguments passed
#' to the \code{binMS} routine, how many m/z levels were chosen for each
#' criterion, how many candidate compounds were chosen overall, and how many
#' candidate compounds were obtained after consolidation.
#'
#' @param object An object of class \code{\link{binMS}}
#'
#' @param ... Arguments passed to dot-dot-dot are ignored
#'
#' @export


summary.binMS <- function(object, ...) {

    # Add pointers to summ_info variables for convenience
    n_tot        <- object$summ_info$n_tot
    n_time_pr    <- object$summ_info$n_time_pr
    n_mass       <- object$summ_info$n_mass
    n_charge     <- object$summ_info$n_charge
    n_tiMaCh     <- object$summ_info$n_tiMaCh
    n_binned     <- object$summ_info$n_binned
    time_range   <- object$summ_info$time_range
    mass_range   <- object$summ_info$mass_range
    charge_range <- object$summ_info$charge_range
    mtoz_diff    <- object$summ_info$mtoz_diff
    time_diff    <- object$summ_info$time_diff

    cat("\n",
        "The mass spectrometry data prior to binning had:\n",
        "------------------------------------------------\n",
        "    ", format(n_tot, big.mark=","), " m/z levels\n",
        "\n", sep="")

    cat("The inclusion criteria was specified as follows:\n",
        "------------------------------------------------\n",
        "    time of peak retention:  between ", time_range[1],   " and ", time_range[2], "\n",
        "    mass:                    between ", mass_range[1],   " and ", mass_range[2], "\n",
        "    charge:                  between ", charge_range[1], " and ", charge_range[2], "\n",
        "\n", sep="")

    nlev <- sapply(c("n_time_pr", "n_mass", "n_charge"), function(x) {
        formatC(get(x), format="d", big.mark=",")
    })
    len_nlev <- nchar(nlev)
    max_nlev <- max( len_nlev )
    combin <- formatC(n_tiMaCh, format="d", big.mark=",")
    cat("The number of remaining m/z levels after filtering by the inclusion criteria was:\n",
        "---------------------------------------------------------------------------------\n",
        "    time of peak retention:  ", rep(" ", max_nlev - len_nlev[1]), nlev[1], "\n",
        "    mass:                    ", rep(" ", max_nlev - len_nlev[2]), nlev[2], "\n",
        "    charge:                  ", rep(" ", max_nlev - len_nlev[3]), nlev[3], "\n",
        "    all combined:            ", rep(" ", max_nlev - nchar(combin)), combin, "\n",
        "\n", sep="")

    cat("m/z levels were consolidated when each of the following criteria were met:\n",
        "--------------------------------------------------------------------------\n",
        "    (i)   m/z levels were no more than ", mtoz_diff, " units apart\n",
        "    (ii)  the time peak retention occured no farther apart than ", time_diff, " units\n",
        "    (iii) the charge states were the same\n",
        "\n", sep="")

    cat("After consolidating the m/z levels, there were:\n",
        "-----------------------------------------------\n",
        "    ", formatC(n_binned, format="d", big.mark=","), " levels\n",
        "\n", sep="")
}
