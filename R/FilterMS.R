
#' Filter compounds from mass spectrometry data
#'
#' Desc here ******************
#'
#' @param msObj An object of class \code{\link{msDat}}
#'
#' @param region asdfas
#'
#' @param bor_size
#'
#' @param min_inten
#'
#' @param max_chg

# 1. The m/z intensity maxima must fall inside the range of the bioactivity “area of interest”
# 2. The m/z intensity of species meeting the first criteria must be < 5 % of its respective maximum peak
#      intensity in the areas bordering said “area of interest”
# 3. There must be non-zero abundance in the right adjacent fraction to the maximum intensity fraction
# 4. The intensity must be > 1,000 in active window
# 5. All charge states > +7 are excluded

filterMS <- function(msObj, region, bor_size=2L, min_inten=1000, max_chg=7L) {

  # TODO: check type of region, bor_size

  ms <- msObj$ms
  mtoz <- msObj$mtoz
  chg <- msObj$chg

  ms_nr <- nrow(ms)
  ms_nc <- ncol(ms)
  bor_size <- as.integer(bor_size)

  # Create region index variable: regIdx
  regIdx <- reg_to_idx(msObj, NULL, region, "ms")

  # Create border index variable: borIdx
  reg_lo <- head(regIdx, 1)
  if (reg_lo > 1L) {
    bef_lo <- max(1L, reg_lo - bor_size)
    bef_hi <- reg_lo - 1L
    bef_seq <- seq(bef_lo, bef_hi)
  }
  else {
    bef_seq <- integer(0)
  }
  reg_hi <- tail(regIdx, 1)
  if (reg_hi < ms_nc) {
    aft_lo <- reg_hi + 1L
    aft_hi <- min(reg_hi + bor_size, ms_nc)
    aft_seq <- seq(aft_lo, aft_hi)
  }
  else {
    aft_seq <- integer(0)
  }
  borIdx <- c(bef_seq, aft_seq)

  # maxIdx: the column index per row (and hence fraction) of the maximum intensity level
  maxIdx <- apply(ms, 1, which.max)
  crit1_bool <- maxIdx %in% regIdx

  #crit1_bool <- logical(ms_nr)
  crit2_bool <- logical(ms_nr)
  crit3_bool <- logical(ms_nr)
  crit4_bool <- logical(ms_nr)
  crit5_bool <- logical(ms_nr)

  for (i in seq_len(ms_nr)) {

    #crit1_bool[i] <- which.max(ms[i, ]) %in% regIdx
    crit2_bool[i] <- (crit1_bool[i] && all(20 * ms[i, borIdx] < ms[i, maxIdx[i]]))
    crit3_bool[i] <- (ms[i, min(maxIdx[i] + 1, ms_nc)] > 0)
    crit4_bool[i] <- any(ms[i, ] > min_inten)
    crit5_bool[i] <- (chg[i] <= max_chg)
  }

  crit_nm <- paste0("crit", 1:5, "_bool")

  critIdx <- lapply(crit_nm, function(x) which(get(x)))
  names(critIdx) <- paste0("crit_", 1:5)


  #keepIdx <- which( Reduce("&", lapply(crit_nm, get)) )
  keepIdx <- which( Reduce("&", list(crit1_bool, crit2_bool, crit3_bool, crit4_bool, crit5_bool)) )
  ms_nm <- colnames(ms)
  if (is.null(ms_nm)) {
    ms_nm <- as.character(seq_len(ms_nc))
  }

  outObj <- list( msObj     = msDat(ms[keepIdx, ], mtoz[keepIdx], chg[keepIdx]),
                  orig_dim  = c(ms_nr, ms_nc),
                  critIdx   = critIdx,
                  reg_nm    = ms_nm[regIdx],
                  bor_nm    = ms_nm[borIdx],
                  bor_len   = bor_size,
                  min_inten = min_inten,
                  max_chg   = max_chg )

  structure(outObj, class="filterMS")
}




summary.filterMS <- function(filtObj) {

  cat("\n",
      "The mass spectrometry data prior to filtering had:\n",
      "--------------------------------------------------\n",
      "    ", format(filtObj$orig_dim[1], width=5, big.mark=","), " compounds\n",
      "    ", format(filtObj$orig_dim[2], width=6, big.mark=","), " fractions\n",
      "\n", sep="")

  cat("The region of interest was specified as (", length(filtObj$reg_nm), " fractions):\n",
      rep("-", 53 + nchar(length(filtObj$reg_nm))), "\n", sep="")
  for (nm in filtObj$reg_nm) {
    cat(nm, "\n", sep="")
  }
  cat("\n")

  cat("The bordering regions were specified as having length ",
      filtObj$bor_len, ", corresponding to:\n",
      rep("-", 73 + nchar(filtObj$bor_len)), "\n", sep="")
  for (nm in filtObj$bor_nm) {
    cat(nm, "\n", sep="")
  }
  cat("\n")

  cat("The minimum intensity was specified as: ",
      format(filtObj$min_inten, width=6, big.mark=","), "\n",
      "The maximum charge was specified as:    ", format(filtObj$max_chg, width=7), "\n",
      "\n", sep="")

}



