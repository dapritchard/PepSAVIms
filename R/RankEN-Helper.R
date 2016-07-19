
rankEN_vector_to_matrix <- function(vec) {

  # case: not a vector
  if ( !(is.vector(vec) && !is.list(vec)) ) {
    return (vec)
  }

  attr(vec, "dim") <- c(1, length(vec))
  attr(vec, "dimnames") <- list(NULL, attr(vec, "names"))
  attr(vec, "names") <- NULL
  vec
}



rankEN_check_regr_args <- function(ms, bio) {

  # Ensure that the dimensions match
  if ( !identical(ncol(ms), ncol(bio)) ) {
    stop("Number of fractions for mass spectrometry (", ncol(ms), ") ",
         "does not match the number of fractions for bioactivity (", ncol(bio), ")", call.=FALSE)
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

rankEN_comp_entrance <- function(enet_fit) {

  actions <- unlist(enet_fit$actions)
  # enet (and lars) adds an entry to the end of the actions list signaling the
  # number of actions performed
  actions <- actions[-length(actions)]

  unique(actions[actions > 0])
}




rankEN_comp_cor <- function(ms_t, bio_vec) {
  apply(ms_t, 2, function(x) stats::cor(x, bio_vec))
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




rankEN_check_valid_input <- function(msObj, bioact, region_ms, region_bio,
                                     lambda, pos_only, ncomp) {

  ## Check for missing arguments

  all_var_nm <- c("msObj", "bioact", "region_ms", "region_bio", "lambda", "pos_only", "ncomp")
  for (var_nm in all_var_nm) {
    if (!eval(substitute(hasArg(var_nm)))) {
      stop("Must provide an argument for ", var_nm, call.=FALSE)
    }
    # Check that an object exists for provided argument
    tryCatch(get(var_nm), error = function(err) {
      err <- as.character(err)
      obj_nm <- regmatches(err, gregexpr("(?<=\')(.*?)(?=\')", err, perl=TRUE))[[1L]]
      stop("object \'", obj_nm, "\' not found for ", var_nm, call.=FALSE)
    })
  }

  ## Check bioact

  if (!is.numeric(bioact) && !is.data.frame(bioact)) {
    stop("bioact must be either a data.frame or of mode numeric", call.=FALSE)
  }
  # Check that bioact is a vector or matrix.  Note that we don't check for NAs
  # or non-numeric yet b/c we only care about problems in the region of interest.
  else if (is.numeric(bioact)) {
    bio_dim <- attr(bioact, "dim")
    if (!identical(bio_dim, NULL) && !identical(length(bio_dim), 2L)) {
      stop("If of mode numeric, then bioact must be a vector or a matrix", call.=FALSE)
    }
    # We require matrices to have number of cols >= 2 for extract_var(), and
    # when bioact is a vector then it is turned into a (1 x p) matrix later
    else if (identical(bio_dim, NULL) && (length(bioact) < 2L)) {
      stop("If a numeric vector, then bioact must have length >= 2", call.=FALSE)
    }
  }
  # Note: don't check for NAs or non-numeric yet b/c we only care about problems
  # in the region of interest

  ## Check region_ms, region_bio

  for (var_nm in c("region_ms", "region_bio")) {
    x <- get(var_nm)
    if (!is.null(x) && !is.numeric(x) && !is.character(x)) {
      stop(var_nm, " must be either NULL or either of mode numeric or character", call.=FALSE)
    }
    else if (anyNA(x)) {
      stop(var_nm, " cannot contain any missing", call.=FALSE)
    }
  }

  ## Check lambda

  if (!is.numeric(lambda)) {
    stop("lambda must be a numeric value", call.=FALSE)
  }
  else if (anyNA(lambda)) {
    stop("lambda cannot contain any missing", call.=FALSE)
  }
  else if ( !identical(length(lambda), 1L) ) {
    stop("lambda must be an atomic vector of length 1", call.=FALSE)
  }
  else if (is.na(lambda) || (lambda < 0)) {
    stop("lambda cannot be smaller than 0", call.=FALSE)
  }

  ## Check pos_only

  if (!identical(pos_only, TRUE) && !identical(pos_only, FALSE)) {
    stop("pos_only must be either TRUE or FALSE", call.=FALSE)
  }

  ## Check ncomp

  if (!is.null(ncomp)) {
    if (!is.numeric(ncomp)) {
      stop("ncomp must be either NULL or a numeric value", call.=FALSE)
    }
    else if (anyNA(ncomp)) {
      stop("If non-NULL then ncomp cannot contain any missing", call.=FALSE)
    }
    else if ( !identical(length(ncomp), 1L) ) {
      stop("If non-NULL then ncomp must be an atomic vector of length 1", call.=FALSE)
    }
    else if (is.na(ncomp) || (ncomp < 1)) {
      stop("If non-NULL then ncomp must be >= 1", call.=FALSE)
    }
  }
}
