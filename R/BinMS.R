
# raw_data <- read.csv("/home/dpritch/Documents/Projects/Hicks Bioactivity/Data/MS_First1000.csv",
#                      skip=2, row.names=1)
#
#
#
#
# x <- read.csv("/home/dpritch/Documents/Projects/Hicks Bioactivity/Data/Sweet violet peptide ion data.csv",
#               skip=2, row.names=1, nrows=10)

# library(readr)
#
# file_path <- "/home/dpritch/Documents/Projects/Hicks Bioactivity/Data/Sweet violet peptide ion data.csv"
# mtoz <- "m/z"
# charge <- "Charge"
# mass <- "Mass"
# time_peak_reten <- "Retention time (min)"
# ms_inten <- 33:70
# time_pr_range <- c(14, 45)
# mass_range <- c(2000, 15000)
# charge_range <- c(2L, 10L)
# mtoz_diff <- 0.05
# time_pr_diff <- 1
#
# #raw_data <- read_csv(file_path, skip=2)
#
# x <- binMS(file_path, mtoz, charge, mass, time_peak_reten, ms_inten, time_pr_range,
#       mass_range, charge_range, mtoz_diff, time_pr_diff, skip=2)

#' Consolidate mass spectrometry observations
#'
#' Combines mass spectrometry observations that are believed to belong to the
#' same true observation into a single observation.  In concept, the data
#' produced by the mass spectrometer may contain multiple reads for a single
#' underlying compound.  Thus, this function attempts to recover these
#' underlying compounds through a binning procedure, described in more detail
#' in \code{Details}.
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
#' @author Pritchard

binMS <- function(mass_spec, mtoz, charge, mass, time_peak_reten, ms_inten, time_pr_range,
                  mass_range, charge_range, mtoz_diff, time_pr_diff, ...) {

  # If necessary read data from file into a dataframe (or leave as is)
  if( !is.matrix(mass_spec) && !is.data.frame(mass_spec) ) {
    raw_data <- read_csv(file=mass_spec, ...)
  } else {
    raw_data <- mass_spec
  }

  # Calculate row (mass-to-charge level) criteria status
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

  # Reduce data by row criteria and by selecting columns needed in future.  Transpose
  # so as to conform to column major order.

  # ms: transpose of the mass spec data
  ms <- t( raw_data[keepIdx, ms_inten] )
  ms <- ms[, sortIdx]
  # info: data used to identify and combine the mass-to-charge levels
  info <- t( raw_data[keepIdx, c(mtoz, charge, time_peak_reten, mass)] )
  info <- info[, sortIdx]


  ## Binning mass-to-charge values

  # Allocate storage for binned data
  n <- length(keepIdx)
  info_bin <- matrix(nrow=4, ncol=n)
  ms_bin <- matrix(nrow=nrow(ms), ncol=n)

  # Row numbers for the information array (after transposing)
  rmtoz   <- 1L
  rcharge <- 2L
  rtime   <- 3L
  rmass   <- 4L

  binIdx  <- 1L    # Track index for binned data
  currIdx <- 1L    # Track index for current mass-to-charge level in loop
  nextIdx <- 1L    # Track index for next mass-to-charge to compare current against
  ncmp    <- 1     # Number of mass-to-charge levels in current bin

  # Create an illegal charge value used to signal that a mass-to-charge level
  # has already been included as part of a bin
  flagv <- charge_range[1] - 1

  # Each iteration compares an m/z level to see if we can start a new bin.  If
  # the current level of the iteration isn't part of a bin, then we start a new
  # bin and compare to larger m/z levels, adding levels to the bin when criteria
  # are met, and ending when m/z levels gets out of range.
  while (currIdx <= n) {

    # case: current row was already included as part of a bin.  Move on to next
    # m/z level
    if (info[rcharge, currIdx] < charge_range[1]) {
      currIdx <- currIdx + 1L
      next
    }

    # Begin a new bin with current m/z level
    info_bin[, binIdx] <- info[, currIdx]
    ms_bin[, binIdx] <- ms[, currIdx]
    nextIdx <- currIdx + 1L
    ncmp <- 1

    # Each iteration compares a m/z value to the current m/z level until the
    # values being compared against have a m/z outside of the allowed range. If
    # all of the criteria for the iteration are met, then the m/z level is
    # combined with the current bin

    while ( (nextIdx <= n)
            && (info[rmtoz, nextIdx] - (info_bin[rmtoz, binIdx] / ncmp) < mtoz_diff) ) {

      # case: criteria met.  Add m/z level to current bin.
      if ( (info[rcharge, nextIdx] == info_bin[rcharge, binIdx])
           && (abs(info[rtime, nextIdx] - (info_bin[rtime, binIdx] / ncmp)) < time_pr_diff) ) {

        # add m/z level to current bin
        info_bin[, binIdx] <- info_bin[, binIdx] + info[, nextIdx]
        ms_bin[, binIdx] <- ms_bin[, binIdx] + ms[, nextIdx]

        # signal that m/z level is already part of a bin
        info[rcharge, nextIdx] <- flagv
        ncmp <- ncmp + 1
      }
      # case: criteria not me
      else {
        # noop
      }

      info_bin[, binIdx] <- info_bin[, binIdx] / ncmp
      nextIdx <- nextIdx + 1
    } # end compare current bin to next m/z level loop

    currIdx <- currIdx + 1
    binIdx <- binIdx + 1
  } # end start a new bin loop (when m/z level not already in a bin)

  list(info_bin[, seq_len(binIdx)], ms_bin[, seq_len(binIdx)])
}









# y <- as.matrix()
#
# path <- "~/Documents/Projects/Hicks Bioactivity/20160218_VO vs. E. faecium bioactivity.csv"
#
#
# binned_dat <- read.csv("/home/dpritch/Documents/Projects/Hicks Bioactivity/20160215_CLK_BAP_VO_binned_short.csv",
#                        row.names=1)
#
# msObj <- msDat(binned_dat, 1, 2)
#
# filterObj <- filterMS(msObj, paste0("_", 14:22), bord_rat=0.10, min_inten=-1, max_chg=7L)

