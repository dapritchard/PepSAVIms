
#' Consolidate mass spectrometry observations
#'
#' Combines mass spectrometry observations that are believed to belong to the
#' same underlying compound into a single observation.  In concept, the data
#' produced by the mass spectrometer may produce multiple reads for a single
#' compound; thus, \code{binMS} attempts to recover these underlying compounds
#' through a binning procedure, described in more detail in \code{Details}.
#'
#' @param mass_spec
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

binMS <- function(mass_spec, mtoz, charge, mass, time_peak_reten, ms_inten, time_pr_range,
                  mass_range, charge_range, mtoz_diff, time_pr_diff, ...) {

  # If necessary read data from file into a dataframe (or leave as is)
  if( !is.matrix(mass_spec) && !is.data.frame(mass_spec) ) {
    raw_data <- read_csv(file=mass_spec, ...)
  } else {
    raw_data <- mass_spec
  }

  ## Step 1: construct (sorted) indices of rows in the data that satisfy the
  ## mass, time of peak retention, and charge criteria

  # Calculate row (i.e. mass-to-charge level) criteria status
  time_pr_bool <- ( (time_pr_range[1] <= raw_data[, time_peak_reten])
                    & (raw_data[, time_peak_reten] <= time_pr_range[2]) )

  mass_bool <- ( (mass_range[1] <= raw_data[[mass]])
                 & (raw_data[[mass]] <= mass_range[2]) )

  charge_bool <- ( (charge_range[1] <= raw_data[[charge]])
                   & (raw_data[[charge]] <= charge_range[2]) )

  keepIdx <- which( Reduce("&", list(time_pr_bool, mass_bool, charge_bool)) )

  # Obtain indices of ordered data after removing m/z levels that didn't meet
  # row criteria
  sortIdx <- order(raw_data[keepIdx, mtoz],
                   raw_data[keepIdx, time_peak_reten],
                   raw_data[keepIdx, charge])

  ## Step 2: reduce data by row criteria and by selecting columns needed in
  ## future.  Transpose so as to conform to column major order.

  # ms: transpose of the mass spec data
  ms <- t( raw_data[keepIdx, ms_inten] )
  ms <- ms[, sortIdx]
  # info: data used to identify and combine the mass-to-charge levels
  info_orig <- t( raw_data[keepIdx, c(mtoz, charge, time_peak_reten, mass)] )
  info_orig <- info_orig[, sortIdx]


  ## Step 3: perform binning of observations that meet the similarity criteria
  ## based on mass-to-charge values, time of peak retention, and charge state

  # Allocate memory for binned data
  n_bef_comb <- length(keepIdx)
  info_bin <- matrix(nrow=4, ncol=n_bef_comb)
  ms_bin <- matrix(nrow=nrow(ms), ncol=n_bef_comb)

  # Row numbers for the information array (after transposing)
  rmtoz  <- 1L
  rchg   <- 2L
  rtime  <- 3L
  rmass  <- 4L

  binIdx     <- 1L  # Index for binned data
  origDatIdx <- 1L  # Index for original data
  cmpMtozIdx <- 1L  # Index for original data within inner loop that looks for
                    # all occurences of data within the allowed m/z difference
  ncmp       <- 1L  # Number of mass-to-charge levels in current bin

  # Create an illegal charge value used to signal that a mass-to-charge level
  # has already been included as part of a bin
  flagv <- charge_range[1] - 1

  # Each iteration compares an m/z level to see if we can start a new bin.  If
  # the current level of the iteration isn't part of a bin, then we start a new
  # bin and compare to larger m/z levels, adding levels to the bin when criteria
  # are met, and ending when m/z levels gets out of range.
  
  while (origDatIdx <= n_bef_comb) {

    # case: current row was already included as part of a bin; move on to next
    # m/z level in the pre-binned data (recall that when a row was already
    # included as part of a bin then the charge is set to (charge_range - 1)
    # as a signal)
    if (info_orig[rchg, origDatIdx] < charge_range[1]) {
      origDatIdx <- origDatIdx + 1L
      next
    }

    # Begin a new bin with current m/z level
    info_bin[, binIdx] <- info_orig[, origDatIdx]
    ms_bin[, binIdx] <- ms[, origDatIdx]
    cmpMtozIdx <- origDatIdx + 1L
    ncmp <- 1L

    # Each iteration compares a m/z value to the current m/z level until the
    # values being compared against have a m/z outside of the allowed range. If
    # all of the criteria for the iteration are met, then the m/z level is
    # combined with the current bin.
    #
    # Note that when comparisons are made that we divide by ncmp (i.e. the number
    # of compounds in the bin) since we are only summing observations from the
    # original data when placing it in the bin (and taking the mean occurs after
    # the (following inner) loop) ends.

    while ( (cmpMtozIdx <= n_bef_comb)
            && (info_orig[rmtoz, cmpMtozIdx] - (info_bin[rmtoz, binIdx] / ncmp)
              < mtoz_diff) ) {

      # case: criteria met.  Add m/z level to current bin.
      if ( (info_orig[rchg, cmpMtozIdx] == (info_bin[rchg, binIdx] / ncmp))
           && (abs(info_orig[rtime, cmpMtozIdx] - (info_bin[rtime, binIdx] / ncmp))
             < time_pr_diff) ) {

        # add m/z level to current bin
        info_bin[, binIdx] <- info_bin[, binIdx] + info_orig[, cmpMtozIdx]
        ms_bin[, binIdx] <- ms_bin[, binIdx] + ms[, cmpMtozIdx]

        # signal that m/z level is already part of a bin
        info_orig[rchg, cmpMtozIdx] <- flagv
        ncmp <- ncmp + 1L
      }
      # case: criteria not me
      else {
        # noop
      }

      cmpMtozIdx <- cmpMtozIdx + 1L
    } # end compare current bin to next m/z level loop

    info_bin[, binIdx] <- info_bin[, binIdx] / ncmp
    origDatIdx <- origDatIdx + 1L
    binIdx <- binIdx + 1L
  } # end start a new bin loop (when m/z level not already in a bin)

  
  list(info_bin[, seq_len(binIdx - 1L)], ms_bin[, seq_len(binIdx - 1L)])
}
