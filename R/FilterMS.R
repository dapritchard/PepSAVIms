
#' Filter compounds from mass spectrometry data
#'
#' Filters mass spectrometry data using a set of criteria, described in
#' \code{Details}. Returns an object of classes \code{\link{msDat}} and
#' \code{filterMS}.
#'
#' @param msObj An object either of class \code{binMS} or of class \code{\link{msDat}}
#'
#' @param region A vector either of mode character or mode numeric.  If numeric
#'   then the entries should provide the indices for the region of interest in
#'   the mass spectrometry data in the argument for \code{msObj}.  If character
#'   then the entries should uniquely specify the region of interest through
#'   partial string matching (see criterion 1, 4).
#'
#' @param border Either a character string \code{"all"}, or a character string
#'   \code{"none"}, or a length-1 or length-2 numeric value specifying the
#'   number of fractions to either side of the region of interest to comprise
#'   the bordering region.  If a single numeric value, then this is number of
#'   fractions to each side of the region of interest; if it is two values, the
#'   the first value is the number of fractions to the left, and the second
#'   value is the number of fractions to the right.  If there are not enough
#'   fractions in either direction to completely span the number of specified
#'   fractions, then all of the available fractions to the side in question are
#'   considered to be part of the bordering region (see criterion 2).
#'
#' @param bord_ratio A single nonnegative numeric value.  A value of 0 will not
#'   admit any compounds, while a value greater than 1 will admit
#'   all compounds (see criterion 2).
#'
#' @param min_inten A single numeric value.  A value less than the minimum mass
#'   spectrometry value in the data will admit all compounds (see criterion 4).
#'
#' @param max_chg A single numeric value specifying the maximum charge which a
#'   compound may exhibit (see criterion 5)
#'
#' @details Attempts to filter out candidate compounds via subject-matter
#'   knowledge, with the goal of removing spurious noise from downstream models.
#'   The criteria for the downstream inclusion of a candidate compound is listed
#'   below.
#'
#'   \enumerate{
#'
#'   \item The m/z intensity maximum must fall inside the range of the
#'   bioactivity region of interest
#'
#'   \item The ratio of the m/z intensity of a species in the areas bordering
#'   the region of interest and the species maximum intensity must be less than
#'   \code{bord_ratio}
#'
#'   \item The immediately right adjacent fraction to its maximum intensity
#'   fraction for a species must have a non-zero abundance.  In the case of ties
#'   for the maximum, it is the fraction immediately to the right of the
#'   rightmost maximum fraction which cannot have zero abundance.
#'
#'   \item At least 1 fraction in the region of interest must have intensity
#'   greater than \code{min_inten}
#'
#'   \item Compound charge state must be less than or equal to \code{max_chg}
#'
#'   }
#'
#' @return Returns an object of classes \code{filterMS} and \code{\link{msDat}}.
#'   This object is a \code{list} with elements described below.  The classes
#'   are equipped with a summary function.
#'
#'   \describe{
#'
#'   \item{\code{msObj}}{ An object of class \code{\link{msDat}} such that the
#'   encapsulated mass spectrometry data corresponds to each of the candidate
#'   compounds that satisfed each of the criteria.  If no criteria are satisfied
#'   then \code{NULL} is returned. }
#'
#'   \item{\code{cmp_by_crit}}{ A list containing \code{data.frame}s, one for
#'   each criterion. Each row (if any) in one of the sub-\code{data.frame}s
#'   contains the mass-to-charge and charge information for a candidate compound
#'   that satisfies the criterion represented by the \code{data.frame}; all of
#'   the compounds that satisfied the criterion are included in the data.  The
#'   \code{data.frame}s are named \code{c1}, ..., \code{c5}, etc corresponding
#'   to criterion 1, ..., criterion 5. }
#'
#'   \item{\code{summ_info}}{ A list containing information pertaining to the
#'   filtering process; for use by the summary function. }
#'
#'   }
#'
#' @export


filterMS <- function(msObj, region, border="all", bord_ratio=0.05, min_inten=1000, max_chg=7L) {

  # Check validity of arguments
  filterMS_check_valid(msObj, region, border, bord_ratio, min_inten, max_chg)

  # Number of criterion
  nCrit <- 5

  # Create pointers to mass spec variables for convenience
  msDatObj <- extract_msDat(msObj, "msDat")
  ms <- msDatObj$ms
  mtoz <- msDatObj$mtoz
  chg <- msDatObj$chg

  # Dimensions of mass spectrometry data
  ms_nr <- nrow(ms)
  ms_nc <- ncol(ms)

  # Create region index variable
  regIdx <- filterMS_getRegionIdx(region, ms)
  
  # Create border index, i.e. the indices that surround the region of interest
  borIdx <- filterMS_getBorderIdx(border, regIdx, ms_nc)

  # maxIdx: the column index per row (and hence fraction) of the maximum intensity level.
  # The rightmost column index is chosen in the case of ties so as
  maxIdx <- apply(ms, 1, function(x) tail(which(x == max(x)), 1))

  # Evaluate criteria predicates
  critBool <- data.frame( array(dim=c(ms_nr, nCrit)) )
  row_seq <- seq_len(ms_nr)
  critBool[, 1] <- maxIdx %in% regIdx
  critBool[, 2] <- sapply(row_seq, function(i) all(ms[i, borIdx] < bord_ratio * ms[i, maxIdx[i]]))
  critBool[, 3] <- sapply(row_seq, function(i) (ms[i, min(maxIdx[i] + 1, ms_nc)] > 0))
  critBool[, 4] <- sapply(row_seq, function(i) any(ms[i, regIdx] > min_inten))
  critBool[, 5] <- (chg <= max_chg)

  # Create a vector of indices which satisfy every criterion
  keepIdx <- which( Reduce("&", critBool) )
  ms_nm <- colnames(ms)
  if (is.null(ms_nm)) {
    ms_nm <- as.character(seq_len(ms_nc))
  }

  # Create filtered mass spectrometry data
  if (length(keepIdx) > 0) {
    msObj <- msDat(ms[keepIdx, , drop=FALSE], mtoz[keepIdx], chg[keepIdx])
  }
  else {
    msObj <- NULL
    warning("There are no compounds that met all of the criteria\n", call.=FALSE)
  }

  # Create mass-to-charge and charge datasets for each criterion
  cmp_by_cr <- setNames(vector("list", nCrit), paste0("c", seq_len(nCrit)))
  for (j in seq_len(nCrit)) {
    thisKeep <- critBool[, j]
    # note: data.frame can handle the case when thisKeep has length 0
    cmp_by_cr[[j]] <- data.frame( mtoz = mtoz[thisKeep],
                                  chg  = chg[thisKeep] )
  }

  # Construct return object
  outObj <- list( msObj    = msObj,
                  cmp_by_cr = cmp_by_cr,
                  summ_info = list( orig_dim   = c(ms_nr, ms_nc),
                                    reg_nm     = ms_nm[regIdx],
                                    bor_nm     = ms_nm[borIdx],
                                    border     = border,
                                    bord_ratio = bord_ratio,
                                    min_inten  = min_inten,
                                    max_chg    = max_chg ) )

  structure(outObj, class="filterMS")
}




#' Overview of the filtering process
#'
#' Prints a text description of the filtering process.  Displays arguments
#' chosen for the \code{filterMS} constructor, how many candidate compounds were
#' chosen for each criterion, and how many candidate compounds were chosen
#' overall.
#'
#' @export


summary.filterMS <- function(filtObj) {

  # Add variables to current environment for convenience: orig_dim, reg_nm,
  # bor_nm, border, bord_ratio, min_inten, max_chg
  list2env(filtObj$summ_info, envir=environment())

  cat("\n",
      "The mass spectrometry data prior to filtering had:\n",
      "--------------------------------------------------\n",
      "    ", format(orig_dim[1], width=5, big.mark=","), " compounds\n",
      "    ", format(orig_dim[2], width=6, big.mark=","), " fractions\n",
      "\n", sep="")

  cat("The region of interest was specified as (", length(reg_nm), " fractions):\n",
      rep("-", 53 + nchar(length(reg_nm))), "\n", sep="")
  for (nm in reg_nm) {
    cat(nm, "\n", sep="")
  }
  cat("\n")

  if ( identical(border, "all") ) {
    cat("- The bordering regions were specified as:  everthing not the region of interest\n")
  }
  else if ( (identical(border, "none")) || all(border == 0) ) {
    cat("- The bordering regions were specified as:  no bordering regions\n")
  }
  else if ( identical(length(border), 1L) ) {
    cat("The bordering regions were specified as each having length ",
        border, ", corresponding to:\n",
        rep("-", 78 + nchar(border)), "\n", sep="")
    for (nm in bor_nm) {
      cat(nm, "\n", sep="")
    }
    cat("\n")
  }
  else {
    cat("The bordering regions were specified as having lengths ",
        border[1], " and ", border[2], ", corresponding to:\n",
        rep("-", 79 + sum(nchar(border))), "\n", sep="")
    for (nm in bor_nm) {
      cat(nm, "\n", sep="")
    }
    cat("\n")
  }

  mi <- format(min_inten, big.mark=",")
  cat("- The minimum intensity was specified as:   ", mi, "\n",
      "- The maximum charge was specified as:",
      rep(" ", 5 + nchar(mi), sep=""),  max_chg, "\n",
      "- The bordering region ratio was:  ",
      rep(" ", 8 + nchar(mi), sep=""), format(bord_ratio, digits=2, nsmall=2), "\n",
      "\n", sep="")

  cat("Individually, each criterion reduced the ",
      q <- format(orig_dim[1], big.mark=","),
      " fractions to the following number:\n", sep="")
  cat(rep("-", 80 + length(q)), "\n", sep="")

  ncri <- sapply(filtObj$cmp_by_cr, function(x) format(nrow(x), big.mark=","))
  plen <- sapply(ncri, nchar)
  mlen <- max(plen)
  bordperc <- paste0(format(100 * bord_ratio, digits=0), "%")

  cat("Criterion 1:  ", rep(" ", mlen - plen[1]), ncri[1],
      "    (fraction with maximum abundance is in region of interest)\n",
      "Criterion 2:  ", rep(" ", mlen - plen[2]), ncri[2],
      "    (fractions in bordering region have < ", bordperc, " of maximum abundance)\n",
      "Criterion 3:  ", rep(" ", mlen - plen[3]), ncri[3],
      "    (nonzero abundance in right adjacent fraction to maximum)\n",
      "Criterion 4:  ", rep(" ", mlen - plen[4]), ncri[4],
      "    (at least 1 intensity > ", format(min_inten, big.mark=","), " in region of interest)\n",
      "Criterion 5:  ", rep(" ", mlen - plen[5]), ncri[5],
      "    (must have charge <= ", max_chg, ")\n",
      "\n", sep="")

  totcmp <- ifelse(is.null(filtObj$msObj), 0, format(nrow(filtObj$msObj$ms), big.mark=","))
  cat("The total number of candidate compounds was reduced to:\n",
      "-------------------------------------------------------\n",
      rep(" ", 14 + mlen - nchar(totcmp)), totcmp,      
      "\n\n", sep="")
}


