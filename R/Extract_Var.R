
extract_var <- function(callArg, dataObs, expect_matr_bool=FALSE, ...) {
  
  # Extract variable name from callArg; used for error messages
  arg_nm <- deparse( substitute(callArg) )

  ## Figure out in which way the data is presented to us

  if (is.null(callArg)) {
    argType <- "null"
  }
  else if (identical(length(callArg), NROW(dataObs))) {
    argType <- "data"
  }
  else if (is.numeric(callArg)) {
    argType <- "nume"
  }
  else if (is.character(callArg)) {
    argType <- "char"
  }
  else {
    stop(arg_nm, " is not of the right type")
  }

  ## Call the appropriate fcn and return the extracted data

  switch(argType,
         null = extract_null(dataObs, arg_nm, ...),
         data = extract_data(callArg, arg_nm),
         nume = extract_nume(callArg, dataObs, arg_nm, expect_matr_bool),
         char = extract_char(callArg, dataObs, arg_nm, expect_matr_bool),
         stop("Invalid switch type"))
}




extract_null <- function(dataObs, arg_nm, ...) {

  other_vars <- list(...)

  # Obtain a list with each element a length-0 or length-1 vector providing the
  # index in dataObs corresponding to the variable from ... being passed
  other_idx_list <- lapply(other_vars, function(callArg) {
    this_arg_nm <- deparse( substitute(callArg) )

    if (identical(length(callArg), NROW(dataObs))) {
      return (integer(0))
    }
    else if (is.numeric(callArg)) {
      extract_check_valid(callArg, this_arg_nm, TRUE)
      return ( extract_num_to_idx(callArg, dataObs, this_arg_nm) )
    }
    else if (is.character(callArg)) {
      extract_check_valid(callArg, this_arg_nm, TRUE)
      return ( extract_char_to_idx(callArg, dataObs, this_arg_nm) )
    }
    else {
      stop(this_arg_nm , " is not of the right type")
    }
  })
  # Obtain an integer vector
  other_idx <- unlist(other_idx_list)

  # Return dataObs after removing non-ms columns and converting to a matrix
  varIdx <- setdiff(seq_len( NCOL(dataObs) ), other_idx)
  extract_idx_to_data(varIdx, dataObs, arg_nm, TRUE)
}




extract_data <- function(callArg, arg_nm) {
  
  if (!is.numeric(callArg)) {
    stop("If ", var_nm, " is the same length of the data then it must be of ",
         "mode numeric", call.=FALSE)
  }

  return (callArg)
}




extract_nume <- function(callArg, dataObs, arg_nm, expect_matr_bool) {

  extract_check_valid(callArg, arg_nm, expect_matr_bool)

  varIdx <- extract_num_to_idx(callArg, dataObs, arg_nm)

  extract_idx_to_data(varIdx, dataObs, arg_nm, expect_matr_bool)
}




extract_char <- function(callArg, dataObs, arg_nm, expect_matr_bool) {

  extract_check_valid(callArg, arg_nm, expect_matr_bool)

  varIdx <- extract_char_to_idx(callArg, dataObs, arg_nm)

  extract_idx_to_data(varIdx, dataObs, arg_nm, expect_matr_bool)
}




extract_num_to_idx <- function(callArg, dataObs, arg_nm, expect_matr_bool) {
  
  # Check valid input values
  if ((min(callArg) < 1) || (max(callArg) > NCOL(dataObs))) {
    stop(arg_nm, ":  out of bounds region value provided", call.=FALSE)
  }

  return (callArg)
}




extract_char_to_idx <- function(callArg, dataObs, arg_nm, expect_matr_bool) {

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
      stop("name provided not in data:  ", nm, "from ", arg_nm, call.=FALSE)
    }
    else if (!identical(nMatch, 1L)) {
      stop("name provided had multiple matches in data:  ", nm, "from ", arg_nm, call.=FALSE)
    }

    return (matchIdx)
  })
}




extract_check_valid <- function(callArg, arg_nm, expect_matr_bool) {

  # Check that there are no duplicates
  if ( !identical(length(unique(callArg)), length(callArg)) ) {
    stop(arg_nm, " cannot have any duplicate values", call.=FALSE)
  }
  
  # Check that vector data has only 1 entry
  else if (!expect_matr_bool && !identical(length(callArg), 1L)) {
    stop(arg_nm, " must have length 1", call.=FALSE)
  }
}




extract_idx_to_data <- function(varIdx, dataObs, arg_nm, expect_matr_bool) {

  # Subset to desired data
  outDat <- dataObs[, varIdx, drop=!expect_matr_bool]

  ## Ensure that data is numeric

  # case: data is a data.frame; convert to a matrix
  if (is.data.frame(outDat)) {
    if ( !all( sapply(outDat, is.numeric) ) ) {
      stop("The data specified for ", arg_nm, "must be numeric", call.=FALSE)
    }
    outDat <- as.matrix(outDat)
  }
  # case: data is atomic; noop
  else {
    if ( !is.numeric(outDat) ) {
      stop("The data specified for ", arg_nm, "must be numeric", call.=FALSE)
    }
  }

  if (any( is.na(outDat) )) {
    stop("Data corresponding to ", var_nm, "contains missing")
  }
  return (outDat)
}




# extract_var <- function(callArg, dataObs, expect_matr_bool=FALSE) {
#   arg_nm <- deparse( substitute(callArg) )
  
#   # case: callArg itself is the variable
#   if (identical(length(callArg), NROW(dataObs))) {
#     if (!is.numeric(callArg)) {
#       stop("If ", var_nm, " is the same length of the data then it must be of ",
#            "mode numeric", call.=FALSE)
#     }
#     return (callArg)
#   }

#   # case: variable specified by either column names or column indices
#   else if (is.numeric(callArg) || is.character(callArg)) {
#     return ( extract_var_char_num(callArg, dataObs, arg_nm, expect_matr_bool) )
#   }

#   # case: not the variable itself, or column names or column indices; should
#   # not have called this function
#   else {
#     stop("One of the preconditions not met before calling this function; "
#          arg_nm, " is not of the right type")
#   }
# }

  


# extract_var_char_num <- function(callArg, dataObs, arg_nm, expect_matr_bool=FALSE) {

#   ## Error checking: validity of callArg arguments
  
#   # Check that there are no duplicates
#   if ( !identical(length(unique(callArg)), length(callArg)) ) {
#     stop(arg_nm, " cannot have any duplicate values", call.=FALSE)
#   }
#   # Check that vector data has only 1 entry
#   else if (!expect_matr_bool && !identical(length(callArg), 1L)) {
#     stop(arg_nm, " must have length 1", call.=FALSE)
#   }

#   ## Obtain variable indices using callArg

#   # case: numeric region provided
#   if (is.numeric(callArg)) {
#     varIdx <- num_to_idx(callArg, dataObs, arg_nm, expect_matr_bool)
#   }
#   # case: character region provided
#   else Idxif (is.character(callArg)) {
#     var <- char_to_idx(callArg, dataObs, arg_nm, expect_matr_bool)
#   }
#   # case: neither numeric or character; should not have called this function
#   else {
#     stop("One of the preconditions not met before calling this function; "
#          arg_nm, " is not of the right type")
#   }

#   ## Extract data from dataObs using variable indices

#   outDat <- dataObs[, callArg, drop=expect_matr_bool]

#   ## Error checking: ensure that data is numeric (hence not a data.frame either)

#   # case: data is a data.frame; convert to a matrix
#   if (is.data.frame(outDat)) {
#     if ( !all( sapply(outDat, is.numeric) ) ) {
#       stop("The data specified for ", arg_nm, "must be numeric", call.=FALSE)
#     }
#     outDat <- as.matrix(outDat)
#   }
#   # case: data is atomic; noop
#   else {
#     if ( !is.numeric(outDat) ) {
#       stop("The data specified for ", arg_nm, "must be numeric", call.=FALSE)
#     }
#   }

#   return (outDat)
# }




# extract_ms <- function(ms_inten, dataObs, expect_matr_bool=FALSE, ...) {

#   if ( !is.null(ms_inten) ) {
#     return ( extract_var(ms_inten, dataObs, TRUE) )
#   }

#   # else: ms_inten is null.  Have to find which variables in dataObs (if any)
#   # belong to non-ms data (i.e. other variables such as mtoz etc), and then
#   # return dataObs minus these other vars.
  
#   other_vars <- list(...)
#   remIdx <- 1
  
    
# }
