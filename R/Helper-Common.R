
# `````````````````````````````````````````````````````` #
#  Provide helper functions common to multiple routines  #
# ...................................................... #

extract_var <- function(callArg, dataObs, expect_matr_bool=FALSE) {
  
  nVars <- NCOL(dataObs)

  if (is.null(callArg)) {
    # TODO
    return (1)
  }
  # case: were already passed the variable
  else if (identical(length(callArg), nVars)) {
    return (callArg)
  }
  # else: the variable is included as part of a data.frame or matrix

  # Check that there are no duplicates
  if ( !identical(length(unique(callArg)), length(callArg)) ) {
    stop("region cannot have any duplicate values", call.=FALSE)
  }

  # case: numeric region provided
  if (is.numeric(callArg)) {
    return ( num_to_idx(callArg, dataObs, expect_matr_bool) )
  }
  # case: character region provided
  else {
    return ( char_to_idx(callArg, dataObs, expect_matr_bool) )
  }
}




num_to_idx <- function(callArg, dataObs, expect_matr_bool) {
  
  # Check valid input values
  if ((min(regVec) < 1) || (max(regVec) > NCOL(dataObs))) {
    stop("out of bounds region value provided", call.=FALSE)
  }

  if (expect_matr_bool) {
    return ( dataObs[, callArg, drop=FALSE] )
  }
  else {
    return ( dataObs[, callArg] )
  }
}




char_to_idx <- function(callArg, dataObs, expect_matr_bool) {

  var_nm <- colnames(dataObs)
  if (is.null(var_nm)) {
    stop("Variable specified by name but data columns not equipped with names")
  }

  varIdx <- sapply(callArg, function(nm) {
    
    # Number of matches for current element of regVec in fraction names
    matchIdx <- grep(nm, var_nm)
    nMatch <- length(matchIdx)

    # Check that current name has exactly one match
    if (identical(nMatch, 0L)) {
      stop(paste0("name provided not in data:  ", nm), call.=FALSE)
    }
    else if (!identical(nMatch, 1L)) {
      stop(paste0("name provided had multiple matches in data", nm), call.=FALSE)
    }

    return (matchIdx)
  })

  if (expect_matr_bool) {
    return ( dataObs[, callArg, drop=FALSE] )
  }
  else {
    return ( dataObs[, callArg] )
  }
}
