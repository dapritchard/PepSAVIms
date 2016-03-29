
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

    if (identical(length(callArg), NROW(dataObs))) {
      return (integer(0))
    }
    else if (is.numeric(callArg)) {
      extract_check_valid(callArg, dataObs, "other var", TRUE)
      return ( extract_num_to_idx(callArg, dataObs, "other var") )
    }
    else if (is.character(callArg)) {
      extract_check_valid(callArg, dataObs, "other var", TRUE)
      return ( extract_char_to_idx(callArg, dataObs, "other var") )
    }
    else {
      stop("other var" , " is not of the right type")
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

  extract_check_valid(callArg, dataObs, arg_nm, expect_matr_bool)

  varIdx <- extract_num_to_idx(callArg, dataObs, arg_nm)

  extract_idx_to_data(varIdx, dataObs, arg_nm, expect_matr_bool)
}




extract_char <- function(callArg, dataObs, arg_nm, expect_matr_bool) {

  extract_check_valid(callArg, dataObs, arg_nm, expect_matr_bool)
  
  varIdx <- extract_char_to_idx(callArg, dataObs, arg_nm)
  
  extract_idx_to_data(varIdx, dataObs, arg_nm, expect_matr_bool)
}




extract_num_to_idx <- function(callArg, dataObs, arg_nm, expect_matr_bool) {

  # Check valid input values
  if ((min(callArg) < 1) || (max(callArg) > NCOL(dataObs))) {
    stop(arg_nm, ":  out of bounds variable value provided", call.=FALSE)
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
    matchIdx <- grep(nm, var_nm, fixed=TRUE)
    nMatch <- length(matchIdx)

    # Check that current name has exactly one match
    if (identical(nMatch, 0L)) {
      stop("name provided not in data:  ", nm, " from ", arg_nm, call.=FALSE)
    }
    else if (!identical(nMatch, 1L)) {
      stop("name provided had multiple matches in data:  ", nm, "from ", arg_nm, call.=FALSE)
    }

    return (matchIdx)
  })
}




extract_check_valid <- function(callArg, dataObs, arg_nm, expect_matr_bool) {

  # Check that data is of the right type
  if (!is.matrix(dataObs) && !is.data.frame(dataObs)) {
    stop("data must be a matrix or data.frame")
  }
  
  # Check that there are no duplicates
  else if ( !identical(length(unique(callArg)), length(callArg)) ) {
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

  # Ensure that data is numeric
  extract_checkifnumeric(outDat, arg_nm)

  # case: row names are lost in coercion from data.frame to atomic
  if (is.data.frame(dataObs)) {
    if (!expect_matr_bool) {
      names(outDat) <- row.names(dataObs)
    }
    else {
      outDat <- as.matrix(outDat)
      row.names(outDat) <- row.names(dataObs)
    }
  }

  # Data may not have any missingness
  if (any( is.na(outDat) )) {
    stop("Data corresponding to ", arg_nm, " contains missing")
  }
  
  return (outDat)
}




extract_checkifnumeric <- function(outDat, arg_nm) {

  # case: data is a data.frame; convert to a matrix
  if (is.data.frame(outDat)) {
    if ( !all( sapply(outDat, is.numeric) ) ) {
      stop("The data specified for ", arg_nm, " must be numeric", call.=FALSE)
    }
  }
  # case: data is atomic; noop
  else {
    if ( !is.numeric(outDat) ) {
      stop("The data specified for ", arg_nm, " must be numeric", call.=FALSE)
    }
  }
}
