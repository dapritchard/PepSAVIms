
#' Consolidate mass spectrometry observations
#'
#' Combines mass spectrometry observations that are believed to belong to the
#' same underlying compound into a single observation.  In concept, the data
#' produced by the mass spectrometer may produce multiple reads for a single
#' compound; thus, \code{binMS} attempts to recover these underlying compounds
#' through a binning procedure, described in more detail in \code{Details}.
#'
#' @param mass_spec Either a matrix or data frame providing the entirety of the
#'   mass spectrometry data.  Must contain columns with data for at least the
#'   following: ******
#'
#' @param mtoz
#'
#' @param charge
#' @param mass
#'
#' @param time_peak_reten
#' @param ms_inten
#' @param time_pr_range
#' @param mass_range
#' @param charge_range
#' @param mtoz_diff
#' @param time_pr_diff
#' @param ...
#' @return out
#'
#' @export

binMS <- function(mass_spec, mtoz, charge, mass, time_peak_reten, ms_inten, time_range,
                  mass_range, charge_range, mtoz_diff, time_diff) {

  # If necessary read data from file into a dataframe (or leave as is)
  # if( !is.matrix(mass_spec) && !is.data.frame(mass_spec) ) {
  #   raw_data <- read_csv(file=mass_spec, ...)
  # } else {
  #   raw_data <- mass_spec
  # }

  ## Step 1: construct sorted indices of rows in the data that satisfy the
  ## mass, time of peak retention, and charge criteria

  # Calculate row (i.e. mass-to-charge level) criteria status
  time_pr_bool <- ( (time_range[1] <= mass_spec[, time_peak_reten])
                    & (mass_spec[, time_peak_reten] <= time_range[2]) )

  mass_bool <- ( (mass_range[1] <= mass_spec[, mass])
                 & (mass_spec[, mass] <= mass_range[2]) )

  charge_bool <- ( (charge_range[1] <= mass_spec[, charge])
                   & (mass_spec[, charge] <= charge_range[2]) )

  keepIdx <- which( Reduce("&", list(time_pr_bool, mass_bool, charge_bool)) )

  # Obtain indices of ordered data after removing m/z levels that didn't meet
  # row criteria
  sortIdx <- order(mass_spec[keepIdx, mtoz],
                   mass_spec[keepIdx, time_peak_reten],
                   mass_spec[keepIdx, charge])

  ## Step 2: reduce data by row criteria and by selecting columns needed in
  ## future.  Transpose so as to conform to column major order.

  # ms: transpose of the mass spec data
  ms <- t( mass_spec[keepIdx, ms_inten, drop=FALSE] )
  ms <- ms[, sortIdx, drop=FALSE]
  # info: data used to identify and combine the mass-to-charge levels
  info_orig <- t( mass_spec[keepIdx, c(mtoz, charge, time_peak_reten, mass)] )
  info_orig <- info_orig[, sortIdx]


  ## Step 3: perform binning of observations that meet the similarity criteria
  ## based on mass-to-charge values, time of peak retention, and charge state

  # Allocate memory for binned data.  We keep the data in two matrices; one for
  # what we call the information which are the values that determine whether or
  # not to combine two m/z levels, and another for the mass spectrometry
  # abundances.  The reason for separating the two is that the information data
  # gets averaged within a bin while the ms abundances are summed.
  n_befbin <- length(keepIdx)
  info_bin <- matrix(nrow=4, ncol=n_befbin)
  ms_bin <- matrix(nrow=nrow(ms), ncol=n_befbin)

  # Give names to rows (i.e. the variables) in the binned data information
  # matrix so as to make the algorithm more readable.  The variables correspond
  # to the m/z value, charge, time of peak retention, and mass respectively.
  rmtoz  <- 1L
  rchg   <- 2L
  rtime  <- 3L
  rmass  <- 4L

  binCtr   <- 1L  # Index for next binned data entry
  origCtr  <- 1L  # Index for next orig data to be considered for binned data
  innerCtr <- 1L  # Index for orig data within inner loop that looks for all
                  # occurences of data within the allowed m/z difference
  nbin     <- 1L  # Number of mass-to-charge levels in current bin

  # Create an out-of-bounds charge value used to signal that a mass-to-charge level
  # has already been included as part of a bin
  flagv <- charge_range[1] - 1

  # Each iteration compares an m/z level to see if we can start a new bin.  If
  # the current level of the iteration isn't part of a bin, then we start a new
  # bin and compare to larger m/z levels, adding levels to the bin when criteria
  # are met, and ending when m/z levels gets out of range.
  
  while (origCtr <= n_befbin) {

    # case: current row was already included as part of a bin; move on to next
    # m/z level in the pre-binned data (recall that when a row was already
    # included as part of a bin then the charge is set to charge_range[1] - 1 as
    # a signal)
    if (info_orig[rchg, origCtr] < charge_range[1]) {
      origCtr <- origCtr + 1L
      next
    }

    # If we've made it here then we've found a m/z level in the original data
    # not yet included in any previous bin.  Begin a new bin starting with this
    # level as its first entry.
    info_bin[, binCtr] <- info_orig[, origCtr]
    ms_bin[, binCtr] <- ms[, origCtr]
    # Set innerCtr to be the index of the first m/z level to consider for adding
    # to the new bin
    innerCtr <- origCtr + 1L
    # Initial number of m/z levels in the bin
    nbin <- 1L

    # Each iteration compares a m/z value to the current m/z level until the
    # values being compared against have a m/z outside of the allowed range. If
    # all of the criteria for the iteration are met, then the m/z level is
    # combined with the current bin.
    #
    # Note that when comparisons are made that we divide by nbin (i.e. the number
    # of compounds in the bin) since we are only summing observations from the
    # original data when placing it in the bin (and taking the mean occurs after
    # the (following inner) loop) ends.

    while ( (innerCtr <= n_befbin)
            && (info_orig[rmtoz, innerCtr] - (info_bin[rmtoz, binCtr] / nbin) < mtoz_diff) ) {

      # case: criteria met.  Add m/z level to current bin.
      if ( (info_orig[rchg, innerCtr] == (info_bin[rchg, binCtr] / nbin))
           && (abs(info_orig[rtime, innerCtr] - (info_bin[rtime, binCtr] / nbin))
             < time_diff) ) {

        # add m/z level to current bin
        info_bin[, binCtr] <- info_bin[, binCtr] + info_orig[, innerCtr]
        ms_bin[, binCtr] <- ms_bin[, binCtr] + ms[, innerCtr]

        # signal that m/z level is already part of a bin
        info_orig[rchg, innerCtr] <- flagv
        nbin <- nbin + 1L
        
      }
      # else: criteria not me; noop
              
      # Update the counter indexing the next m/z level to consider for adding to
      # the current bin
      innerCtr <- innerCtr + 1L
              
    } # end compare current bin to next m/z level loop

    # Take the average of the binned m/z levels (note: don't need to do this
    # for ms_bin as these values are strictly summed)
    info_bin[, binCtr] <- info_bin[, binCtr] / nbin
    
    # We've found every m/z level in the original data that belongs in the
    # binned data.  Update the counter indexing the next m/z level from the
    # original data, and update the counter indexing the next available bin to
    # start placing consolidated data into.
    origCtr <- origCtr + 1L
    binCtr <- binCtr + 1L
    
  } # end iterate over original m/z levels loop

  # It is possible that the binned data can be out of order.  This happens when
  # the minimum observation of one bin is smaller than that of another, but the
  # mean is larger.  binCtr indexes the next available bin, so 1 less is the
  # number of bins.
  nbinned <- binCtr - 1L
  resortIdx <- order(info_bin[rmtoz, 1:nbinned],
                     info_bin[rtime, 1:nbinned],
                     info_bin[rchg, 1:nbinned])

  # Construct msDat object
  msObj <- msDat(t(ms_bin[, resortIdx,  drop=FALSE]),
                 info_bin[rmtoz, resortIdx],
                 info_bin[rchg, resortIdx])

  # Info for use by summary.binMS to describe binning process
  summ_info <- list(n_tot = nrow(mass_spec),
                    n_time_pr = sum(time_pr_bool),
                    n_mass = sum(mass_bool),
                    n_charge = sum(charge_bool),
                    n_tiMaCh = length(keepIdx),
                    n_binned = nbinned,
                    time_range = time_range,
                    mass_range = mass_range,
                    charge_range = charge_range,
                    mtoz_diff = mtoz_diff,
                    time_diff = time_diff)

  # Construct binMS object
  outObj <- list(msObj     = msObj,
                 summ_info = summ_info)
  structure(outObj, class="binMS")
}




#' Overview of the binning process
#'
#' Prints a text description of the binning process.  Displays arguments
#' chosen for the \code{binMS} constructor, how many candidate compounds were
#' chosen for each criterion, and how many candidate compounds were chosen
#' overall.
#'
#' @export


summary.binMS <- function(binObj) {


  cat("\n",
      "The mass spectrometry data prior to binning had:\n",
      "------------------------------------------------\n",
      "    ", format(orig_dim[1], width=5, big.mark=","), " m/z levels\n",
      "\n", sep="")

  cat("The inclusion criteria was specified as follows:\n",
      "------------------------------------------------\n",
      "    time of peak retention:  between ", time_range[1],   " and ", time_range[2], "\n",
      "    mass:                    between ", mass_range[1],   " and ", mass_range[2], "\n",
      "    charge:                  between ", charge_range[1], " and ", charge_range[2], "\n",
      "\n", sep="")

  nlev <- sapply(c("n_time_pr", "n_mass", "n_chg"), function(x) {
    formatC(get(x), format="d", big.mark=",")
  })
  len_nlev <- nchar(nlev)
  max_nlev <- max( len_nlev )
  combin <- formatC(n_tiMaCh, format="d", big.mark=",")
  cat("The number of remaining compounds after filtering by the inclusion criteria was:\n"
      "--------------------------------------------------------------------------------\n"
      "    time of peak retention:  ", rep(" ", max_nlev - len_nlev[1]), nlev[1], "\n",
      "    mass:                    ", rep(" ", max_nlev - len_nlev[2]), nlev[2], "\n",
      "    charge:                  ", rep(" ", max_nlev - len_nlev[3]), nlev[3], "\n",
      "    all combined:            ", rep(" ", max_nlev - nchar(combin)), combin, "\n",
      "\n", sep="")
  
  cat("m/z levels were consolidated when each of the following criteria were met:\n",
      "--------------------------------------------------------------------------\n",
      "    m/z levels no different than ", mtoz_diff, " units apart\n",
      "    the time peak retention occured no farther apart than ", time_diff, " units\n",
      "    the charge states were the same\n"
      "\n", sep="")

  cat("After consolidating the m/z levels, there were:\n",
      "-----------------------------------------------\n",
      "    ", formatC(n_tot, format="d", big.mark=","), "levels\n",
      "\n", sep="")  
}
