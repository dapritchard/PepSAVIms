
#' Ranks compounds using the Elastic Net path
#'
#' Returns identifying information for the compounds in the order in which they
#' first enter the Elastic Net model
#'
#' @param msDat An object of class \code{\link{msDat}} containing the mass
#'   spectrometry data and identifying information
#'
#' @param bioact Either a numeric vector, a matrix, or a data frame providing
#'   bioactivity data.  If a numeric vector, then it is assumed that each entry
#'   corresponds to a particular fraction.  If the data is 2-dimensional, then
#'   it is assumed that each column corresponds to a particular fraction, and
#'   that each row corresponds to a particular bioactivity replicate.
#'
#' @param region Either \code{NULL}, a non-list vector, a matrix, or a list
#'   containing exactly two atomic vectors, providing information specifying
#'   which fractions are to be included in the Lasso model.  Note that a data
#'   frame satisfies these requirements.
#'
#'   If \code{NULL}, then it is assumed that all fractions included in the data
#'   are to be used in the model.  This requires that the number of fractions in
#'   the data for the parameter passed to \code{bioact} be the same as the
#'   number of fractions in the mass spectrometry data for the paramter passed
#'   to \code{msDat}.  It further assumes that the \code{k}-th column must refer
#'   to the same fraction for both the mass spectrometry data and the
#'   bioactivity data for every \code{k}.
#'
#'   If a non-list vector then it must be either a numeric or charcter vector,
#'   such that the vector specifies which columns (and hence which fractions) to
#'   include in the model.  If the vector is numeric, then the desired columns
#'   are specified by number, and if the vector is character, then the desired
#'   columns are specified by name (partial matching is allowed).  Note that
#'   this assumes that the corresponding columns in the mass spectrometry data
#'   and the bioactivity data refer to the same fractions.
#'
#'   If a matrix then it must be either numeric or character with exactly two
#'   columns - one is to be named \code{ms} and the other is to be named
#'   \code{bio}.  The \code{ms} column specifies the columns (and hence
#'   fractions) to include in the model from the mass spectrometry data, either
#'   as a vector of the column numbers or the column names.  The \code{bio}
#'   vector specifies the columns (and hence fractions) to include in the model
#'   from the bioactivity data, either as a vector of the column numbers or the
#'   column names.  It is assumed that two entries from a given row refer to the
#'   same fraction.
#'
#'   If a list, then it must be a list with two named non-list vectors of equal
#'   length - one is to be named \code{ms} and the other is to be named
#'   \code{bio} (note that a \code{n x 2} data frame satisfies this
#'   requirement).  The \code{ms} vector specifies the columns (and hence
#'   fractions) to include in the model from the mass spectrometry data, either
#'   as a vector of the column numbers or the column names.  The \code{bio}
#'   vector specifies the columns (and hence fractions) to include in the model
#'   from the bioactivity data, either as a vector of the column numbers or the
#'   column names.  It is assumed that the column from the mass spectrometry
#'   data specified by the \code{k}-th value in the \code{ms} vector corresponds
#'   to the same fraction as the column specified by the \code{k}-th value in
#'   the \code{bio} vector, for each \code{k}.
#'
#' @details Note that in the current incarnation of rankLasso, the solution to
#'   the Lasso path is the same when using either the average of bioactivity
#'   replicates or individual replicates.  A parameter, useAve, used to be
#'   offered - but since the result is the same either way, now it is just set
#'   \code{TRUE}.  The rest of the code is left unchanged, in the event another
#'   way is found to use individual replicates.
#'
#' @export

# TODO: document function output




rankEN <- function(msObj, bioact, region_ms=NULL, region_bio=NULL, lambda,
                   ncomp=NULL, pos_only=TRUE) {

  # Mung data into the right form ----------------------------------------------
  
  # Obtain msDat obj.  Error checking performed in extractMS.
  msDatObj <- extractMS(msObj, type="msDat")

  # Check that remaining arguments are of the right type
  rankEN_check_valid_input(bioact, region_ms, region_bio, lambda, ncomp)
  
  # Ensure (by coercion if necessary) that bioact is in matrix form
  bioMat <- rankEN_bioact_to_matrix(bioact)
  
  # Extract the region of interest for the ms data and the bioactivity data
  ms <- extract_var(region_ms, msDatObj$ms, TRUE)
  bio <- extract_var(region_bio, bioMat, TRUE)
  # Check for missing and that dimensions match
  rankEN_check_regr_args(ms, bio)

  # Convert ms to form with rows are an observation (i.e. fraction) and cols are
  # a variable (i.e. a compound)
  ms_t <- t(ms)
  
  # Obtain the mean of the bioactivity replicates and convert to a vector
  bio_vec <- colMeans(bio)
  

  # Fit model ------------------------------------------------------------------

  # Calculate the elastic net path
  enet_fit <- elasticnet::enet(ms_t, bio_vec, lambda)
  

  # Extract compound entrance results ------------------------------------------

  # Obtain indices for compounds as they first enter the model
  comp_idx <- rankEN_comp_entrance(enet_fit, ncmp)

  # Correlation of chosen compounds
  comp_cor <- rankEN_comp_cor(ms_t, bio_vec)

  # Filter compounds by correlation and by how many we wish to keep
  comp_idx_out <- rankEN_filter_compIdx(comp_idx, comp_cor, ncomp, pos_only)


  # Construct return object ----------------------------------------------------

  # Extract mass spec and bioactivity fraction names  
  ms_nm <- colnames(ms)
  if (is.null(ms_nm)) {
    ms_nm <- as.character(seq_len(ms_nc))
  }
  bio_nm <- colnames(bio)
  if (is.null(bio_nm)) {
    bio_nm <- as.character(seq_len(bio_nc))
  }  

  # Create info for the summary function
  summ_info <- list(
    orig_dim  = list(ms  = dim(msDatObj$ms),
                     bio = dim(bioMat)),
    reduc_dim = list(ms  = dim(ms),
                     bio = dim(bio)),
    region_nm = list(ms  = ms_nm,
                     bio = bio_nm),
    lambda    = lambda,
    ncomp     = ncomp,
    pos_only  = pos_only
  )

  # Construct output object
  outDat <- list( mtoz      = msDatObj$mtoz[comp_idx_out],
                  charge    = msDatObj$chg[comp_idx_out],
                  comp_cor  = comp_cor[comp_idx_out],
                  enet_fit  = enet_fit,
                  summ_info = summ_info )

  structure(outDat, class="rankLasso")
}




summary.rankLasso <- function(rl_out) {

  # Create a link for convenience
  data_desc <- rl_out$data_desc

  # Print mass spectrometry and bioactivity data dimensions
  cat(sep="",
      "\n",
      "Mass spectrometry data:\n",
      "-----------------------\n",
      "  ", format(data_desc$msDim[1], width=6, big.mark=","), " compounds\n",
      "  ", format(data_desc$msDim[2], width=6, big.mark=","), " fractions\n",
      "\n",
      "Bioactivity data:\n",
      "-----------------\n",
      "  ", format(data_desc$bioDim[1], width=6, big.mark=","), " replicates\n",
      "  ", format(data_desc$bioDim[2], width=6, big.mark=","), " fractions\n",
      "\n")


  # TODO: have to consider case where (data_desc$regionNm$ms == NULL) and similarly for bio

  blkCol <- rep("    ", length(data_desc$regionIdx$ms))
  regDat <- setNames( data.frame( data_desc$regionNm$ms,
                                  data_desc$regionIdx$ms,
                                  blkCol,
                                  data_desc$regionNm$bio,
                                  data_desc$regionIdx$bio ),
                      c("M.S. Names", "  M.S. Index", "", "Bio. Names", "  Bio. Index") )

  # Print a table with the names and indices for mass spectrometry and bioactivity data
  cat(sep="",
      "Fractions included in region of interest:\n",
      "-----------------------------------------\n")
  print(regDat, row.names=FALSE)
  cat("\n")

  # Create a table with the mass-to-charge, charge, and correlation values for
  # chosen compounds
  sigCmp <- setNames( data.frame(rl_out$mtoz,
                                 rl_out$charge,
                                 rl_out$cmp_cor),
                      c("Mass-to-charge", "  Charge", "  Correlation") )
  sig_dig <- nchar( length( data_desc$cmpIdx ) )

  # Print a table with the mass-to-charge, charge, and correlation values for
  # chosen compounds
  cat(sep="",
      "Compounds chosen by the model (", length(data_desc$cmpIdx), " total):\n",
      "-------------------------------",        rep("-", sig_dig), "--------\n")
  print(sigCmp, row.names=FALSE, digits=4)
  cat("\n")

}





