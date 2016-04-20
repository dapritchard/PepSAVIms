
# Extract from dataObs as specified by callArg ---------------------------------

# extract_var is the interface that user-level routines will likely call when
# they want the data in dataObs that are specified by callArg

extract_var <- function(callArg, dataObs, expect_matr_bool=FALSE, ...) {
  
  # Extract variable name from callArg and dataObs; used for error messages
  arg_nm <- deparse( substitute(callArg) )
  dat_nm <- deparse( substitute(dataObs) )

  ## Figure out in which way the data is presented to us

  if (is.null(callArg)) {
    argType <- "null"
  }
  else if (identical(length(callArg), NROW(dataObs))) {
    argType <- "data"
  }
  else if (is.numeric(callArg)) {
    argType <- "num"
  }
  else if (is.character(callArg)) {
    argType <- "char"
  }
  else {
    stop(arg_nm, " is not of the right type")
  }

  ## Call the appropriate  subroutine and return the extracted data

  switch(argType,
         null = extract_null(dataObs, arg_nm, ...),
         data = extract_data(callArg, arg_nm),
         num  = extract_num(callArg, dataObs, expect_matr_bool, arg_nm, dat_nm),
         char = extract_char(callArg, dataObs, expect_matr_bool, arg_nm, dat_nm),
         stop("Invalid switch type"))
}




extract_null <- function(dataObs, arg_nm, ...) {

  # other_vars <- list(...)

  # # Obtain a list with each element a length-0 or length-1 vector providing the
  # # index in dataObs corresponding to the variable from ... being passed
  # other_idx_list <- lapply(other_vars, function(callArg) {

  #   if (identical(length(callArg), NROW(dataObs))) {
  #     return (integer(0))
  #   }
  #   else if (is.numeric(callArg)) {
  #     extract_check_valid(callArg, dataObs, "other var", TRUE)
  #     return ( extract_num_to_idx(callArg, dataObs, "other var") )
  #   }
  #   else if (is.character(callArg)) {
  #     extract_check_valid(callArg, dataObs, "other var", TRUE)
  #     return ( extract_char_to_idx(callArg, dataObs, "other var") )
  #   }
  #   else {
  #     stop("other var is not of the right type", call.=FALSE)
  #   }
  # })
  # # Obtain an integer vector
  # other_idx <- unlist(other_idx_list)

  # # Obtain variable indices after removing indices corresponding to other variables
  # varIdx <- setdiff(seq_len( NCOL(dataObs) ), other_idx)
  varIdx <- extract_null_to_idx(dataObs, arg_nm, ...)

  # Check that we haven't used all the data on other variables
  if (identical(length(varIdx), 0L)) {
    stop("There cannot be 0 columns left for ", arg_nm,
         " after removing data for other variables", call.=FALSE)
  }

  # Return dataObs after removing non-ms columns and converting to a matrix
  extract_idx_to_data(varIdx, dataObs, TRUE, arg_nm, dat_nm)
}




extract_data <- function(callArg, arg_nm) {
  
  if (!is.numeric(callArg)) {
    stop("If ", arg_nm, " is the same length of the data then it must be of ",
         "mode numeric", call.=FALSE)
  }
  else if (anyNA(callArg)) {
    stop(arg_nm, " cannot have any missing")
  }

  return (callArg)
}




extract_num <- function(callArg, dataObs, expect_matr_bool, arg_nm, dat_nm) {

  extract_check_valid(callArg, dataObs, expect_matr_bool, arg_nm, dat_nm)

  varIdx <- extract_num_to_idx(callArg, dataObs, arg_nm)

  extract_idx_to_data(varIdx, dataObs, expect_matr_bool, arg_nm, dat_nm)
}




extract_char <- function(callArg, dataObs, expect_matr_bool, arg_nm, dat_nm) {

  extract_check_valid(callArg, dataObs, expect_matr_bool, arg_nm, dat_nm)

  varIdx <- extract_char_to_idx(callArg, dataObs, arg_nm)
  
  extract_idx_to_data(varIdx, dataObs, expect_matr_bool, arg_nm, dat_nm)
}




extract_idx_to_data <- function(varIdx, dataObs, expect_matr_bool, arg_nm, dat_nm) {

  # Subset to desired data
  outDat <- dataObs[, varIdx, drop=!expect_matr_bool]

  # Ensure that data is numeric and w/o missing
  extract_check_if_numer_nomiss(outDat, arg_nm, dat_nm)

  # case: dataObs is a data.frame; convert to matrix or vector.  However, row
  # names are lost in coercion from data.frame to atomic so we have to recover
  # manually.
  if (is.data.frame(dataObs)) {
    # case: return object is to be a vector; the elements can be accessed named
    # using names()
    if (!expect_matr_bool) {
      names(outDat) <- row.names(dataObs)
    }
    # case: outDat a matrix; use row.names to set dimnames
    else {
      outDat <- as.matrix(outDat)
      row.names(outDat) <- row.names(dataObs)
    }
  }

  # # Ensure that data does not have any missingness
  # if (anyNA(outDat)) {
  #   stop("Data corresponding to ", arg_nm, " contains missing", call.=FALSE)
  # }
  
  return (outDat)
}




# Extracting indices section ---------------------------------------------------

# extract_idx is the interface that internal routines will likely call when they
# need an index set of the columns in dataObs that are specified by callArg

extract_idx <- function(callArg, dataObs, expect_matr_bool=FALSE, ...) {

  # Extract variable name from callArg; used for error messages
  arg_nm <- deparse( substitute(callArg) )
  
  extract_check_valid(callArg, dataObs, expect_matr_bool, arg_nm, dat_nm)

  ## Figure out in which way the data is presented to us
  
  if (is.null(callArg)) {
    argType <- "null"
  }
  else if (is.numeric(callArg)) {
    argType <- "num"
  }
  else if (is.character(callArg)) {
    argType <- "char"
  }
  else {
    stop(arg_nm, " is not of the right type")
  }

  ## Call the appropriate subroutine and return the extracted data
  
  switch(argType,
         null = extract_null_to_idx(callArg, dataObs, arg_nm, ...),
         num  = extract_num_to_idx(callArg, dataObs, expect_matr_bool, arg_nm, dat_nm),
         char = extract_char_to_idx(callArg, dataObs, expect_matr_bool, arg_nm, dat_nm),
         stop("Invalid switch type"))
}




extract_null_to_idx <- function(dataObs, arg_nm, ...) {

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
      stop("other var is not of the right type", call.=FALSE)
    }
  })
  # Obtain an integer vector
  other_idx <- unlist(other_idx_list)

  # Return variable indices after removing indices corresponding to other vars
  setdiff(1:NCOL(dataObs), other_idx)
}




extract_num_to_idx <- function(callArg, dataObs, arg_nm) {

  # Check valid input values
  for (k in callArg) {
    if ((k < 1) || (k > NCOL(dataObs))) {
      stop("out of bounds value ", k, " provided for ", arg_nm, call.=FALSE)
    }
  }

  return (callArg)
}




extract_char_to_idx <- function(callArg, dataObs, arg_nm) {

  var_nm <- colnames(dataObs)
  if (is.null(var_nm)) {
    stop("Variable specified by name but data columns not equipped with names", call.=FALSE)
  }

  varIdx <- sapply(callArg, function(nm) {
    
    # Number of matches for current element of regVec in fraction names
    matchIdx <- grep(nm, var_nm, fixed=TRUE)
    nMatch <- length(matchIdx)

    # Check that current name has exactly one match
    if (identical(nMatch, 0L)) {
      stop("name provided not in data - ", nm, " element in ", arg_nm, call.=FALSE)
    }
    else if (!identical(nMatch, 1L)) {
      stop("name provided had multiple matches in data - ", nm, " element in ",
           arg_nm, call.=FALSE)
    }

    return (matchIdx)
  })
}




# Argument validity checking ---------------------------------------------------

extract_check_valid <- function(callArg, dataObs, expect_matr_bool, arg_nm, dat_nm) {

  # Check that data is of the right type
  if (!is.matrix(dataObs) && !is.data.frame(dataObs)) {
    stop(dat_nm, " must be a matrix or data.frame", call.=FALSE)
  }
  
  # Check that there are no duplicates
  else if ( !identical(length(unique(callArg)), length(callArg)) ) {
    stop(arg_nm, " cannot have any duplicate values", call.=FALSE)
  }
  
  # Check that input corresponding to vector data has exactly 1 entry (assume we
  # already know it is not n-length)
  else if (!expect_matr_bool && !identical(length(callArg), 1L)) {
    stop(arg_nm, " must have length 1 or length equal to the number of observations", call.=FALSE)
  }

  # Check that input corresponding to matrix data has at least 1 entry
  else if (expect_matr_bool && !(length(callArg) >= 1L)) {
    stop("If non-NULL, then ", arg_nm, " must have length >= 1", call.=FALSE)
  }
}




extract_check_if_numer_nomiss <- function(outDat, arg_nm, dat_nm) {

  # case: data is a data.frame; ensure each column is numeric w/o NAs
  if (is.data.frame(outDat)) {
    if ( !all( sapply(outDat, is.numeric) ) ) {
      stop(dat_nm, " cannot have any non-numeric in the region provided by ", arg_nm, call.=FALSE)
    }
    else if ( any( sapply(outDat, anyNA) ) ) {
      stop(dat_nm, " cannot have any missing in the region provided by ", arg_nm, call.=FALSE)
    }
  }
  # case: data is atomic
  else {
    if ( !is.numeric(outDat) ) {
      stop(dat_nm, " cannot have any non-numeric in the region provided by ", arg_nm, call.=FALSE)
    }
    else if (anyNA(outDat)) {
      stop(dat_nm, " cannot have any missing in the region provided by ", arg_nm, call.=FALSE)
    }
  }
}
