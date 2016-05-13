
# Extraction functionality overview --------------------------------------------


# There are two functions in this file which are designed to be directly invoked
# by calling routines - the remaining functions are subroutines for these main
# routines.  The routines designed for direct access by calling routines are the
# following:
#
#   extract_var - returns a subset of the data
#
#   extract_idx - returns indices corresponding to the data which is to be
#                 subsetted
#
# Note that the main routines are the ones that perform the error
# checking, and it is assumed that the arguments are valid when calling the
# subroutines.
#
# The routines in this section have a desired type of a matrix or data.frame for
# data_obs, a numeric or character vector for var_specify, and logical for
# expect_matr




# Extracting data section ------------------------------------------------------


# extract the data in data_obs that is specified by var_specify.
#
# PRE: when we have a 1-row data_obs, and an index argument to var_specify then
# this will erroneously treat the index as data.  Thus a potential calling
# function should not use this routine if this is a possibility.

extract_var <- function(data_obs, var_specify, expect_matr=FALSE, ...) {
  
  # Extract variable name from var_specify and data_obs; used for error messages
  spec_nm <- deparse( substitute(var_specify) )
  dat_nm <- deparse( substitute(data_obs) )

  # Check arguments for type and size validity
  extract_check_valid(data_obs, var_specify, expect_matr, dat_nm, spec_nm)

  ## Figure out in which way the data is presented to us

  if (is.null(var_specify)) {
    if (identical(expect_matr, FALSE)) {
      stop("is var_specify is NULL, then expect_matr must be TRUE")
    }
    argstr <- "null"
  }
  # Note: the second condition in the following case is to handle the case where
  # the length of the region is the same as the length of the data (e.g. if
  # there are 3 replicates in bioactivity data and region is length 3).  This
  # condition is usually able to disambiguate the cases b/c we only expect
  # var_specify to be data (as opposed to column names or indices) in the case
  # of vector arguments (see precondition regarding 1-row data).
  else if ((identical(length(var_specify), NROW(data_obs)))
           && identical(expect_matr, FALSE)) {
    argstr <- "data"
  }
  else if (is.numeric(var_specify)) {
    argstr <- "num"
  }
  else if (is.character(var_specify)) {
    argstr <- "char"
  }
  else {
    stop(spec_nm, " is not of the right type")
  }

  ## Call the appropriate  subroutine and return the extracted data

  switch(argstr,
         null = extract_null(data_obs, dat_nm, spec_nm, ...),
         data = extract_data(var_specify, spec_nm),
         num  = extract_num(data_obs, var_specify, expect_matr, dat_nm, spec_nm),
         char = extract_char(data_obs, var_specify, expect_matr, dat_nm, spec_nm),
         stop("Invalid switch type"))
}




# Return data_obs after removing any columns (i.e. variables in data_obs) that
# belong to other variables as specified in dot-dot-dot

extract_null <- function(data_obs, dat_nm, spec_nm, ...) {

  # Obtain variable indices after removing indices corresponding to other
  # variables (as passed on through variadic arg)
  var_idx <- extract_null_to_idx(data_obs, dat_nm, ...)

  # Check that we haven't used all the data on other variables.  Require 2
  # columns for consistency since we demand this for the following call to
  # extract_idx_to_data.
  if (length(var_idx) < 2L) {
    stop("There must be at least 2 columns left for ", dat_nm,
         " after removing data for other variables", call.=FALSE)
  }

  # Return data_obs after removing non-ms columns and converting to a matrix
  extract_idx_to_data(data_obs, var_idx, TRUE, dat_nm, spec_nm)
}




# var_specify is itself a variable so just return it (i.e. is the identity
# function)

extract_data <- function(var_specify, spec_nm) {
  
  if (!is.numeric(var_specify)) {
    stop("If ", spec_nm, " is the same length of the data then it must be of ",
         "mode numeric", call.=FALSE)
  }
  else if (anyNA(var_specify)) {
    stop(spec_nm, " cannot have any missing")
  }

  return (var_specify)
}




# Return the column(s) in data_obs as specified by a vector of indices provided
# to var_specify

extract_num <- function(data_obs, var_specify, expect_matr, dat_nm, spec_nm) {

  var_idx <- extract_num_to_idx(data_obs, var_specify, dat_nm, spec_nm)

  extract_idx_to_data(data_obs, var_idx, expect_matr, dat_nm, spec_nm)
}




# Return the column(s) in data_obs as specified by a vector of column names
# provided to var_specify

extract_char <- function(data_obs, var_specify, expect_matr, dat_nm, spec_nm) {

  var_idx <- extract_char_to_idx(data_obs, var_specify, dat_nm, spec_nm)
  
  extract_idx_to_data(data_obs, var_idx, expect_matr, dat_nm, spec_nm)
}




# Returns a *numeric object* providing the column(s) in data_obs as specified by
# a vector of indices provided to var_idx

extract_idx_to_data <- function(data_obs, var_idx, expect_matr, dat_nm, spec_nm) {

  # Subset to desired data
  out_dat <- data_obs[, var_idx, drop=!expect_matr]

  # Ensure that data is numeric and w/o missing
  extract_check_if_numer_nomiss(out_dat, dat_nm, spec_nm)

  # case: data_obs is a data.frame; convert to matrix or vector.  However, row
  # names are lost in coercion from data.frame to atomic so we have to recover
  # manually.
  if (is.data.frame(data_obs)) {
    # case: return object is to be a vector; the elements can be accessed named
    # using names()
    if (!expect_matr) {
      names(out_dat) <- row.names(data_obs)
    }
    # case: out_dat a matrix; use row.names to set dimnames
    else {
      out_dat <- as.matrix(out_dat)
      row.names(out_dat) <- row.names(data_obs)
    }
  }
  
  return (out_dat)
}




# Extracting indices section ---------------------------------------------------


# Returns a vector of indices corresponding to columns (i.e. variables) in
# data_obs as specified by var_specify
#
# PRE: assumes that var_specify is not a data vector (i.e. is of type NULL,
# numeric, or character).  This is in contrast to extract_var() where a data
# vector can be passed as the argument to var_specify.

extract_idx <- function(data_obs, var_specify, expect_matr=FALSE, ...) {

  # Extract variable name from var_specify; used for error messages
  dat_nm <- deparse( substitute(data_obs) )
  spec_nm <- deparse( substitute(var_specify) )

  # Check arguments for type and size validity
  extract_check_valid(data_obs, var_specify, expect_matr, dat_nm, spec_nm)

  # Figure out in which way the data is presented to us
  
  if (is.null(var_specify)) {
    if (!identical(expect_matr, TRUE)) {
      stop("is var_specify is NULL, then expect_matr must be TRUE", call.=FALSE)
    }
    argstr <- "null"
  }
  else if (is.numeric(var_specify)) {
    argstr <- "num"
  }
  else if (is.character(var_specify)) {
    argstr <- "char"
  }
  else {
    # Should have restricted to the above cases by extract_check_valid
    stop("Should not be able to reach this part of the program")
  }

  # Call the appropriate subroutine and return the extracted data
  switch(argstr,
         null = extract_null_to_idx(data_obs, spec_nm, ...),
         num  = extract_num_to_idx(data_obs, var_specify, dat_nm, spec_nm),
         char = extract_char_to_idx(data_obs, var_specify, dat_nm, spec_nm),
         stop("Should not be able to reach this part of the program"))
}




# Returns a vector of indices corresponding to columns in data_obs after
# removing any columns (i.e. variables in data_obs) that belong to other
# variables as specified in dot-dot-dot

extract_null_to_idx <- function(data_obs, dat_nm, ...) {

  other_vars <- list(...)

  # Obtain a list with each element a length-0 or length-1 vector providing the
  # index in data_obs corresponding to the variable from dot-dot-dot being
  # passed.  When no args are passed to dot-dot-dot then an empty list is
  # returned.
  other_idx_list <- lapply(other_vars, function(rem_var) {

    if (identical(length(rem_var), NROW(data_obs))) {
      return (integer(0))
    }
    else if (is.numeric(rem_var)) {
      extract_check_valid(data_obs, rem_var, TRUE, dat_nm, "remove_var")
      return ( extract_num_to_idx(data_obs, rem_var, dat_nm, "remove_var") )
    }
    else if (is.character(rem_var)) {
      extract_check_valid(data_obs, rem_var, TRUE, dat_nm, "remove_var")
      return ( extract_char_to_idx(data_obs, rem_var, dat_nm, "remove_var") )
    }
    else {
      stop("remove_var must either be of mode numeric or mode character", call.=FALSE)
    }
  })
  # Obtain an integer vector; returns NULL if other_idx_list is empty
  other_idx <- unlist(other_idx_list)

  # Return variable indices after removing indices corresponding to other vars.
  # If other_idx is NULL then setdiff returns the first argument.
  setdiff(1:NCOL(data_obs), other_idx)
}




# Return a vector of indices corresponding to column(s) in data_obs as specified
# by a vector of indices provided to var_specify; in other words this is the
# identity function

extract_num_to_idx <- function(data_obs, var_specify, dat_nm, spec_nm) {

  # Check valid input values
  for (k in var_specify) {
    if ((k < 1) || (k > NCOL(data_obs))) {
      stop("out of bounds index ", k, " provided for ", spec_nm,
           " relative to ", dat_nm, call.=FALSE)
    }
  }  
  # Check that there are no duplicates in var_specify (note: integer(0) fails
  # the conditional as desired)
  if ( !identical(length(unique(var_specify)), length(var_specify)) ) {
    stop(spec_nm, " cannot have any duplicate values", call.=FALSE)
  }

  # As a side effect this strips away vector names
  as.integer(var_specify)
}




# Return a vector of indices corresponding to the column(s) in data_obs as
# specified by a vector of column names provided to var_specify

extract_char_to_idx <- function(data_obs, var_specify, dat_nm, spec_nm) {

  var_nm <- colnames(data_obs)
  if (is.null(var_nm)) {
    stop("Variable specified by name but data columns not equipped with names", call.=FALSE)
  }
  
  # Check that there are no duplicates in var_specify (note: integer(0) fails
  # the conditional as desired)
  if ( !identical(length(unique(var_specify)), length(var_specify)) ) {
    stop(spec_nm, " cannot have any duplicate values", call.=FALSE)
  }

  out_idx <- sapply(var_specify, function(nm) {
    
    # Number of matches for current element of var_specify in fraction names
    matchIdx <- grep(nm, var_nm, fixed=TRUE)
    nMatch <- length(matchIdx)

    # Check that current name has exactly one match
    if (identical(nMatch, 0L)) {
      stop("column names in ", dat_nm, " do not contain ", nm, " element in ", spec_nm, call.=FALSE)
    }
    else if (!identical(nMatch, 1L)) {
      stop("name provided had multiple matches in data - ", nm, " element in ",
           spec_nm, call.=FALSE)
    }

    return (matchIdx)
  })

  # Return indices after stripping variable names
  setNames(out_idx, NULL)
}




# Argument validity checking ---------------------------------------------------


# Check arguments for type and size validity

extract_check_valid <- function(data_obs, var_specify, expect_matr, dat_nm, spec_nm) {

  # Check that data_obs is of the right type
  if (!is.matrix(data_obs) && !is.data.frame(data_obs)) {
    stop(dat_nm, " must be a matrix or data.frame", call.=FALSE)
  }
  else if (identical(NROW(data_obs), 0L)) {
    stop(dat_nm, " must have number of rows no less than 1", call.=FALSE)
  }
  # Can't have 1 column b/c then it is ambiguous whether a length 1 numeric
  # vector is data or is a column index
  else if (NCOL(data_obs) < 2L) {
    stop(dat_nm, " must have number of columns no less than 2", call.=FALSE)
  }

  # Check that var_specify is of the right form
  if (!is.null(var_specify) && !is.numeric(var_specify) && !is.character(var_specify)) {
    stop(spec_nm, " must either be NULL or either of mode numeric or mode character", call.=FALSE)
  }
  
  if (!is.null(var_specify)) {
    spec_len <- length(var_specify)
    
    # Check that there are no missing in var_specify
    if (anyNA(var_specify)) {
      stop(spec_nm, " cannot contain any missing", call.=FALSE)
    }
    
    # Check that input corresponding to matrix data has at least 1 entry
    if ( identical(length(var_specify), 0L) ) {
      stop("If non-NULL, then ", spec_nm, " must have length no less than 1", call.=FALSE)
    }
    
    # Check that input corresponding to vector data has either exactly 1 entry or
    # exactly the number of entries in the data (i.e. is providing a variable)
    if (!expect_matr && !(identical(spec_len, 1L) || identical(spec_len, nrow(data_obs)))) {
      stop(spec_nm, " must have length 1 or length equal to the number of observations", call.=FALSE)
    }
  }

  # Check that expect_matr is of the right form
  if (!identical(expect_matr, TRUE) && !identical(expect_matr, FALSE)) {
    stop("expect_matr must be either TRUE or FALSE", call.=FALSE)
  }
}




# Ensure that out_dat has columns that are each numeric and contain no missing
# (worded like this because out_dat may either be atomic or a data.frame)

extract_check_if_numer_nomiss <- function(out_dat, dat_nm, spec_nm) {

  # case: data is a data.frame; ensure each column is numeric w/o NAs
  if (is.data.frame(out_dat)) {
    if ( !all( sapply(out_dat, is.numeric) ) ) {
      stop(dat_nm, " cannot have any non-numeric in the region provided by ", spec_nm, call.=FALSE)
    }
    else if ( any( sapply(out_dat, anyNA) ) ) {
      stop(dat_nm, " cannot have any missing in the region provided by ", spec_nm, call.=FALSE)
    }
  }
  # case: data is atomic
  else {
    if ( !is.numeric(out_dat) ) {
      stop(dat_nm, " cannot have any non-numeric in the region provided by ", spec_nm, call.=FALSE)
    }
    else if (anyNA(out_dat)) {
      stop(dat_nm, " cannot have any missing in the region provided by ", spec_nm, call.=FALSE)
    }
  }
}
