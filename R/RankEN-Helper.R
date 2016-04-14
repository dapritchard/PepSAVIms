
rankEN_bioact_to_matrix <- function(bioact) {

  # case: bioact a data.frame
  if (is.data.frame(bioact)) {
    bioact <- as.matrix(bioact)
  }
  # case: bioact a vector; turn into a 1 x n matrix
  else if (identical(attr(bioact, "dim"), NULL)) {
    attr(bioact, "dim") <- c(1, length(bioact))
  }

  # Ensure that number of rows > 0, number of cols > 0
  if ( !(nrow(bioact) >= 1L) ) {
    stop("bioact must have 1 or more rows", call.=FALSE)
  }
  else if ( !(ncol(bioact) >= 1L) ) {
    stop("bioact must have 1 or more columns", call.=FALSE)
  }

  return (bioact)
}




rankEN_check_regr_args <- function(ms, bio) {

  # Ensure no missing in region of interest
  if (anyNA(ms)) {
    stop("mass spec data in region of interest cannot contain any missing", call.=FALSE)
  }
  else if (anyNA(bioact)) {
    stop("bioact cannot contain any missing", call.=FALSE)
  }

  # Ensure that the dimensions match
  if ( !identical(ncol(ms), ncol(bio)) ) {
    stop("Number of fractions for mass spectrometry (", ncol(ms), ") ",
         "does not match the number of fractions for bioactivity (", ncol(bio), ")")
  }
}




# Obtain an integer vector indexing (some of the) compounds by the order in
# which they first enter the Lasso model, from first to last
#
# The enet slot 'actions' is a list of all the actions taken in the model path
# in terms of either adding or removing a predictor variable from model.  If the
# variable is added to the model, then the column index is inserted into the
# next position in the list.  If the variable is removed from the model, then -1
# times the column index of the variable is inserted into the next position in
# the list.

rankEN_comp_entrance <- function(enet_fit, ncmp) {
  
  actions <- unlist(enet_fit$actions)
  # enet (and lars) adds an entry to the end of the actions list signaling the
  # number of actions performed
  actions <- actions[-length(actions)]

  unique(actions[actions > 0])
}




rankEN_comp_cor <- function(ms_t, bio_vec) {
  apply(ms_t, 2, function(x) cor(x, bio_vec))
}




rankEN_filter_compIdx <- function(comp_idx, comp_cor, ncomp, pos_only) {

  # Sort correlations by the order in which the covariate entered the model
  cor_by_idx <- comp_cor[comp_idx]

  # case: NULL value for comp_idx means to keep all values
  if (is.null(ncomp)) {
    ncomp <- length(comp_idx)
  }

  # case: remove all nonpositive correlations
  if (pos_only) {
    comp_idx <- comp_idx[cor_by_idx > 0]
  }

  # case: keep only the first ncomp (possibly positive) compounds
  if (ncomp < length(comp_idx)) {
    comp_idx <- comp_idx[1:ncomp]
  }

  comp_idx  
}




# rankEN_data_desc <- function(ms, bioMat, regionIdx, cmpIdx) {

#   #   # Extract mass spec fraction names
#   # ms_nm <- colnames(ms)
#   # if (is.null(ms_nm)) {
#   #   ms_nm <- as.character(seq_len(ms_nc))
#   # }

#   msDim <- dim(ms)
#   bioDim <- dim(bioMat)

#   regionNm <- list( ms  = colnames(ms)[regionIdx$ms],
#                     bio = colnames(bioMat)[regionIdx$bio] )

#   data_desc <- list( msDim     = msDim,
#                      bioDim    = bioDim,
#                      regionNm  = regionNm,
#                      regionIdx = regionIdx,
#                      cmpIdx    = cmpIdx )
#   return (data_desc)
# }




rankEN_check_valid_input <- function(bioact, region_ms, region_bio, lambda, ncomp) {

  ## Check for missing arguments

  # note: other vars have defaults
  if (missing(bioact)) {
    stop("Must provide an argument for bioact", call.=FALSE)
  }
  else if (missing(lambda)) {
    stop("Must provide an argument for lambda", call.=FALSE)
  }    

  ## Check bioact

  if (!is.numeric(bioact) && !is.data.frame(bioact)) {
    stop("bioact must be either a data.frame or of mode numeric", call.=FALSE)
  }
  # Check that bioact is a vector or matrix
  else if (is.numeric(bioact)) {
    bio_dim <- attr(bioact, "dim")
    if (!identical(bio_dim, NULL) && !identical(length(bio_dim), 2L)) {
      stop("If of mode numeric, then bioact must be a vector or a matrix", call.=FALSE)
    }
  }
  # Don't check for NAs yet b/c we only care about missing in region of interest

  ## Check region_ms, region_bio

  for (var_nm in c("region_ms", "region_bio")) {
    x <- get(var_nm)
    if (!is.null(x) && !is.numeric(x) && !is.character(x)) {
      stop(var_nm, " must be either NULL or of mode numeric or character", call.=FALSE)
    }
    else if (anyNA(x)) {
      stop(var_nm, " cannot contain any missing", call.=FALSE)
    }
  }

  ## Check lambda

  if (!is.numeric(lambda)) {
    stop("lambda must be a numeric value", call.=FALSE)
  }
  else if ( !identical(length(lambda), 1L) ) {
    stop("lambda must have length of 1", call.=FALSE)
  }
  else if (is.na(lambda) || (lambda < 0) || (1 < lambda)) {
    stop("lambda must be between 0 and 1", call.=FALSE)
  }

  ## Check ncomp

  if (!is.null(ncomp)) {
    if (!is.numeric(ncomp)) {
      stop("ncomp must be either NULL or a numeric value", call.=FALSE)
    }
    else if ( !identical(length(ncomp), 1L) ) {
      stop("If non-NULL then ncomp must have length of 1", call.=FALSE)
    }
    else if (is.na(ncomp) || (ncomp < 1)) {
      stop("If non-NULL then ncomp must be >= 1", call.=FALSE)
    }
  }
}
