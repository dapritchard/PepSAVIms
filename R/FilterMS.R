
#' Filter compounds from mass spectrometry data
#'
#' Filters mass spectrometry data using a set of criteria, described in
#' \code{Details}. Returns an object of classes \code{\link{msDat}} and
#' \code{filterMS}.
#'
#' @param msObj An object of class \code{\link{msDat}}
#'
#' @param region A vector either of mode numeric or character.  If numeric then
#'   the entries should provide the indices for the region of interest in the
#'   mass spectrometry data in the argument for \code{msObj}.  If character then
#'   the entries should uniquely specify the region of interest through partial
#'   string matching.
#'
#' @param border Either a character string with the value \code{"all"}, or a
#'   single numeric value specifying the number of fractions to either side of
#'   the region of interest to comprise the bordering region
#'
#' @param bor_rat A single numeric value between 0 and 1, inclusive
#'
#' @param min_inten A single numeric value greater than or equal to 0
#'
#' @param max_chg A single numeric value specifying the maximum charge which a
#'   compound may exhibit
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
#'   \code{bor_rat}
#'
#'   \item The immediately right adjacent fraction to its maximum intensity
#'   fraction for a species must have a non-zero abundance
#'
#'   \item Each fraction in the region of interest must have intensity greater
#'   than \code{min_inten}
#'
#'   \item Compound charge state must be less than or equal to \code{max_chg}
#'
#'   }
#'
#' @return Returns an object of class of classes \code{link{msDat}} and
#'   \code{filterMS}.  This object is a \code{list} with elements described
#'   below.  The classes are equipped with a summary function.
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
#'   \code{data.frame}s are named \code{c1}, \code{c2}, etc corresponding to
#'   criterion 1, criterion 2, and so on. }
#'
#'   \item{\code{summ_info}}{ A list containing information pertaining to the
#'   filtering process; for use by the summary function. }
#'
#'   }
#'
#' @export


filterMS <- function(msObj, region, border="all", bor_rat=0.05, min_inten=1000, max_chg=7L) {

  # TODO: check validity of arguments

  # Number of criterion
  nCrit <- 5

  ms <- msObj$ms
  mtoz <- msObj$mtoz
  chg <- msObj$chg

  # Dimensions of mass spectrometry data
  ms_nr <- nrow(ms)
  ms_nc <- ncol(ms)

  # Create region index variable
  regIdx <- reg_to_idx(msObj, NULL, region, "ms")

  # Create border index, i.e. the indices that surround the region of interest
  borIdx <- getBorderIdx(border, regIdx,ms_nc)

  # maxIdx: the column index per row (and hence fraction) of the maximum intensity level
  maxIdx <- apply(ms, 1, which.max)

  # critBool: containers for critia, one row for each compound and one column for each
  critBool <- data.frame( array(dim=c(ms_nr, nCrit)) )
  row_seq <- seq_len(ms_nr)
  critBool[, 1] <- maxIdx %in% regIdx
  critBool[, 2] <- sapply(row_seq, function(i) all(ms[i, borIdx] < bor_rat * ms[i, maxIdx[i]]))
  critBool[, 3] <- sapply(row_seq, function(i) (ms[i, min(maxIdx[i] + 1, ms_nc)] > 0))
  critBool[, 4] <- sapply(row_seq, function(i) any(ms[i, ] > min_inten))
  critBool[, 5] <- sapply(row_seq, function(i) chg[i] <= max_chg)

  # Create a vector of indices which satisfy every criterion
  keepIdx <- which( Reduce("&", critBool) )
  ms_nm <- colnames(ms)
  if (is.null(ms_nm)) {
    ms_nm <- as.character(seq_len(ms_nc))
  }

  # Create filtered mass spectrometry data
  if (length(keepIdx) > 0) {
    msObj <- msDat(ms[keepIdx, ], mtoz[keepIdx], chg[keepIdx])
  }
  else {
    msObj <- NULL
    warning("There are no compounds that met all of the criteria\n", call.=FALSE)
  }

  # Create mass-to-charge and charge datasets for each criterion
  cmp_by_cr <- setNames(vector("list", nCrit), paste0("c", seq_len(nCrit)))
  for (j in seq_len(nCrit)) {
    thisKeep <- critBool[, j]
    # note: data frame can handle the case when thisKeep is empty
    cmp_by_cr[[j]] <- data.frame( mtoz = mtoz[thisKeep],
                                  chg  = chg[thisKeep] )
  }

  outObj <- list( msObj    = msObj,
                  cmp_by_cr = cmp_by_cr,
                  summ_info = list( orig_dim  = c(ms_nr, ms_nc),
                                    reg_nm    = ms_nm[regIdx],
                                    bor_nm    = ms_nm[borIdx],
                                    border    = border,
                                    min_inten = min_inten,
                                    max_chg   = max_chg ) )

  structure(outObj, class=c("filterMS", "msDat"))
}




#' Overview of the filtering process
#'
#' Prints a text description of the filtering process.  Displays arguments
#' chosen for the \code{filterMS} constructor and how many candidate compounds
#' were chosen for each criterion, as well as overall.
#'
#' @export


summary.filterMS <- function(filtObj) {

  # Add variables to current environment: orig_dim, reg_nm, bor_nm, border, min_inten, max_chg
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
  else {
    cat("The bordering regions were specified as each having length ",
        border, ", corresponding to:\n",
        rep("-", 78 + nchar(border)), "\n", sep="")
    for (nm in bor_nm) {
      cat(nm, "\n", sep="")
    }
    cat("\n")
  }

  cat("- The minimum intensity was specified as: ",
      format(min_inten, width=6, big.mark=","), "\n",
      "- The maximum charge was specified as:    ", format(max_chg, width=7), "\n",
      "\n", sep="")

  cat("Individually, each criterion reduced the ",
      q <- format(orig_dim[1], big.mark=","),
      " fractions to the following number:\n", sep="")
  cat(rep("-", 80 + length(q)), "\n", sep="")

  ncri <- sapply(filtObj$cmp_by_cr, function(x) format(nrow(x), big.mark=","))
  plen <- sapply(ncri, nchar)
  mlen <- max(plen)

  cat("Criterion 1:  ", rep(" ", mlen - plen[1]), ncri[1],
      "    (maximum in region of interest)\n",
      "Criterion 2:  ", rep(" ", mlen - plen[2]), ncri[2],
      "    (< 5% of maximum in bordering areas)\n",
      "Criterion 3:  ", rep(" ", mlen - plen[3]), ncri[3],
      "    (nonzero abundance in right adj. fraction to maximum)\n",
      "Criterion 4:  ", rep(" ", mlen - plen[4]), ncri[4],
      "    (must be an intensity > ", format(min_inten, big.mark=","), " in region of interest)\n",
      "Criterion 5:  ", rep(" ", mlen - plen[5]), ncri[5],
      "    (must have charge <= ", max_chg, ")\n",
      "\n", sep="")

  cat("The total number of candidate compounds was reduced to:\n",
      "-------------------------------------------------------\n",
      "    ", ifelse(is.null(filtObj$msObj), 0, format(nrow(filtObj$msObj$ms), big.mark=",")),
      "\n\n", sep="")

  return (NULL)
}




getBorderIdx <- function(border, regIdx, ms_nc) {

  if ( is.character(border) ) {
    if ( !identical(border, "all") ) {
      stop("If border is of type character then it must have value \"all\"")
    }
    borIdx <- setdiff(seq_len(ms_nc), regIdx)
  }
  # case: border is numeric
  else {
    bsize <- as.integer(border)
    borIdx <- getBorderIdx_numeric(bsize, regIdx, ms_nc)
  }

  return (borIdx)
}




getBorderIdx_numeric <- function(bsize, regIdx, ms_nc) {

  if ( !identical(length(bsize), 1L) ) {
    stop("If border is of mode numeric then it must have length 1")
  }
  else if (bsize < 1L) {
    stop("Value of border must be greater than or equal to 1")
  }

  # Create border index variable: borIdx
  reg_lo <- head(regIdx, 1)
  if (reg_lo > 1L) {
    bef_lo <- max(1L, reg_lo - bsize)
    bef_hi <- reg_lo - 1L
    bef_seq <- seq(bef_lo, bef_hi)
  }
  else {
    bef_seq <- integer(0)
  }
  reg_hi <- tail(regIdx, 1)
  if (reg_hi < ms_nc) {
    aft_lo <- reg_hi + 1L
    aft_hi <- min(reg_hi + bsize, ms_nc)
    aft_seq <- seq(aft_lo, aft_hi)
  }
  else {
    aft_seq <- integer(0)
  }

  borIdx <- c(bef_seq, aft_seq)
}



