
# Note: since extract_idx and extract_var are not user-level functions, we
# assume that we have already checked for argument missingness / existence in
# the calling function (and hence do not check for these issues in any of the
# routines in the Extract_Var.R file).


obs <- matrix(1:50, nrow=10, dimnames=list(1:10, paste0("var", 1:5)))

# v2 <- obs[, 2]

var_idx <- c(2L, 4L)
other_var_num <- setdiff(1:ncol(obs), var_idx)

var_nm <- paste0("r", c(2, 4))
var_nm_zero_len <- character(0L)
var_nm_nomatch <- c("var2", "varxyz")
var_nm_dups <- c("var2", "var2")
# other_char <- paste0("r", other_num)
# col_idx_NA <- replace(col_idx, 1, NA)
var_idx_too_small <- replace(var_idx, 1L, -1L)
var_idx_too_big <- replace(var_idx, 1L, 100L)
var_idx_zero_len <- integer(0L)
var_idx_dups <- c(2L, 2L)
var_idx_NA <- replace(var_idx, 1L, NA_integer_)
# extract_char_no_match <- c("r2", "z4")
# col_idx_zero_len <- integer(0)
# extract_char_zero_len <- character(0)
# col_idx_dup <- rep(2L, 2)
# extract_char_dup <- rep("r2", 2)



# obs_NA_out <- obs_NA_in <- obs
# obs_NA_out[2, 5] <- NA
# obs_NA_in[2, 4] <- NA

# obs_short <- obs[1:3, ]
obs_zero_row <- obs[integer(0), ]
obs_zero_row_df <- data.frame(obs_zero_row)
obs_zero_col <- obs[, integer(0)]
obs_one_col <- obs[, 1L, drop=FALSE]


# obs_df <- obs_df_char_out <- obs_df_char_in <- data.frame(obs)
# obs_df_char_out$var1 <- rep("a", 10)
# obs_df_char_in$var2 <- rep("b", 10)



#
#  
# 

context("extract_idx method")

test_that("extract_idx: valid input with various forms", {

  # Tests are described in terms of their input to var_specify, with a few
  # exceptions at the end

  # numeric length 1: integer
  expect_identical(extract_idx(obs, 1L),
                   1L)
  
  # numeric length 1: integer with names attribute
  expect_identical(extract_idx(obs, setNames(1L, "xyz")),
                   1L)
  
  # numeric length 1: double
  expect_identical(extract_idx(obs, 1),
                   1L)
  
  # numeric length 1: double not approximately an integer
  expect_identical(extract_idx(obs, 1.3),
                   1L)
  
  # numeric length > 1: integer
  expect_identical(extract_idx(obs, var_idx, TRUE),
                   c(2L, 4L))
  
  # numeric length > 1: integer (data_obs a data.frame)
  expect_identical(extract_idx(obs_df, var_idx, TRUE),
                   c(2L, 4L))

  # numeric length 1: specify as a matrix
  expect_identical(extract_idx(obs, 1L, TRUE),
                   1L)

  # character length 1: name
  expect_identical(extract_idx(obs, "var1"),
                   1L)

  # character length 1: specify as a matrix
  expect_identical(extract_idx(obs, "var2", TRUE),
                   2L)

  # character length 1: vector with names attribute
  expect_identical(extract_idx(obs, setNames("var1", "xyz")),
                   1L)

  # character length > 1: names
  expect_identical(extract_idx(obs, var_nm, TRUE),
                   c(2L, 4L))

  # character length > 1: names (data_obs a data.frame)
  expect_identical(extract_idx(obs_df, var_nm, TRUE),
                   c(2L, 4L))

  # NULL: no other vars to dot-dot-dot
  expect_identical(extract_idx(obs, NULL, TRUE), 1:5)

  # NULL: single numeric vector to dot-dot-dot
  expect_identical(extract_idx(obs, NULL, TRUE, other_var_num),
                   c(2L, 4L))

  # NULL: multiple numeric vectors to dot-dot-dot
  expect_identical(extract_idx(obs, NULL, TRUE, 1, 2.5, pi, 4L),
                   5L)

  # NULL: single character vector to dot-dot-dot
  expect_identical(extract_idx(obs, NULL, TRUE, other_char),
                   c(2L, 4L))

  # NULL: a numeric vector and a character vector to dot-dot-dot
  expect_identical(extract_idx(obs, NULL, TRUE, 2, other_char),
                   4L)

  # NULL: with vars passed to dot-dot-dot (data_obs a data.frame)
  expect_identical(extract_idx(obs_df, NULL, TRUE, 2, other_char),
                   4L)
})


test_that("extract_idx: input of the wrong type", {

  wrong_type <- Sys.Date()

  # data_obs wrong type
  expect_error(extract_idx(wrong_type, col_idx),
               "wrong_type must be a matrix or data.frame" )

  # var_specify wrong type
  expect_error(extract_idx(obs, wrong_type),
               "wrong_type must either be NULL or either of mode numeric or mode character" )

  # expect_matr wrong type
  expect_error(extract_idx(obs, var_idx, wrong_type),
               "expect_matr must be either TRUE or FALSE" )

  # dot-dot-dot arg wrong type (first arg)
  expect_error(extract_idx(obs, NULL, TRUE, wrong_type),
               "remove_var must either be of mode numeric or mode character" )

  # dot-dot-dot arg wrong type (second arg)
  expect_error(extract_idx(obs, NULL, TRUE, other_var_num, wrong_type),
               "remove_var must either be of mode numeric or mode character" )
})


test_that("extract_idx: right type but illegal values", {

  # data_obs: 0 row matrix
  expect_error(extract_idx(obs_zero_row, col_idx, TRUE),
               "obs_zero_row must have number of rows no less than 1")

  # data_obs: 0 row data.frame
  expect_error(extract_idx(obs_zero_row_df, col_idx, TRUE),
               "obs_zero_row must have number of rows no less than 1")

  # note: 1-row matrices allowed, in contrast to 1-column matrices

  # data_obs: 0 col matrix
  expect_error(extract_idx(obs_zero_col, col_idx, TRUE),
               "obs_zero_col must have number of columns no less than 2")

  # data_obs: 1 col matrix
  expect_error(extract_idx(obs_one_col, col_idx, TRUE),
               "obs_one_col must have number of columns no less than 2")         

  # note: we don't check data_obs for valid values (i.e. no NAs or non-numeric
  # entries) with extract_idx because we never subset the data.  In other words,
  # if the data is needed then the checking will presumably occur when it is
  # extracted, via some other function.

  # var_specify: zero-length (numeric)
  expect_error(extract_idx(obs, var_idx_zero_len),
               "If non-NULL, then var_idx_zero_len must have length no less than 1")

  # var_specify: zero-length (character)
  expect_error(extract_idx(obs, var_nm_zero_len),
               "If non-NULL, then var_nm_zero_len must have length no less than 1")

  # var_specify: index too small
  expect_error( extract_idx(obs, var_idx_too_small, TRUE),
               "out of bounds index -1 provided for var_idx_too_small relative to obs")

  # var_specify: index too big
  expect_error( extract_idx(obs, var_idx_too_big, TRUE),
               "out of bounds index 100 provided for var_idx_too_big relative to obs")

  # var_specify: name doesn't match any in data_obs column names
  expect_error(extract_idx(obs, var_nm_nomatch, TRUE),
               "column names in obs do not contain varxyz element in var_nm_nomatch")

  # var_specify: duplicates (numeric)
  expect_error(extract_idx(obs, var_idx_dups, TRUE),
               "var_idx_dups cannot have any duplicate values")
  
  # var_specify: duplicates (character)
  expect_error(extract_idx(obs, var_nm_dups, TRUE),
               "var_nm_dups cannot have any duplicate values")

  # var_specify: wrong length wrt expect_matr
  expect_error(extract_idx(obs, var_idx),
               "var_idx must have length 1 or length equal to the number of observations")

  # var_specify NA's
  expect_error(extract_idx(obs, var_idx_NA),
               "var_idx_NA cannot contain any missing")
})



















## extract_var

# num
expect_identical( extract_var(obs, 1L), obs[, 1] )
expect_identical( extract_var(obs, 1L, TRUE), obs[, 1, drop=FALSE] )
expect_identical( extract_var(obs, col_idx, TRUE), obs[, col_idx] )
expect_identical( extract_var(obs_NA_out, col_idx, TRUE), obs[, col_idx] )


# char
expect_identical( extract_var(obs, "var2"), obs[, 2] )
expect_identical( extract_var(obs, "var3", TRUE), obs[, 3, drop=FALSE] )
expect_identical( extract_var(obs, extract_char, TRUE), obs[, col_idx] )
expect_identical( extract_var(obs_NA_out, extract_char, TRUE), obs[, col_idx] )


# NULL
expect_identical( extract_var(obs, NULL, TRUE), obs )
expect_identical( extract_var(obs, NULL, TRUE, other_num), obs[, col_idx] )
expect_identical( extract_var(obs, NULL, TRUE, 1, 2.5, pi, 4L), obs[, 5, drop=FALSE] )
expect_identical( extract_var(obs, NULL, TRUE, other_char), obs[, col_idx] )
expect_identical( extract_var(obs, NULL, TRUE, 2, other_char), obs[, 4, drop=FALSE])

# with data.frame
expect_identical( extract_var(obs_df, col_idx, TRUE), obs[, col_idx] )
expect_identical( extract_var(obs_df, extract_char, TRUE), obs[, col_idx] )
expect_identical( extract_var(obs_df, NULL, TRUE, other_char), obs[, col_idx] )
expect_identical( extract_var(obs_df_char_out, col_idx, TRUE), obs[, col_idx] )
expect_identical( extract_var(obs_df_char_out, extract_char, TRUE), obs[, col_idx] )
expect_identical( extract_var(obs_df_char_out, NULL, TRUE, other_char), obs[, col_idx] )

# data
expect_identical( extract_var(obs, v2), v2 )
expect_identical( extract_var(obs_df_char_out, v2), v2 )
# If handed data that is same length as the data matrix then return the data.
# Compare this to the next case.
expect_identical( extract_var(obs_short, 2:4), 2:4 )
# If handed a region specifier that is same length as the data matrix then
# return the data matrix subset to the region.  Compare this to the previous
# case.
expect_identical( extract_var(data.frame(obs_short), 2:4, TRUE), obs_short[, 2:4] )

# TODO: 1-row data



# Illegal vars --------------------------------------------------










  






# data_obs NA in specified region
extract_var(obs_NA_in, col_idx, TRUE)
extract_var(obs_df_char_in, col_idx, TRUE)

# var_specify NA in specified region

# var_specify invalid indices
expect_error( extract_var(obs, col_idx_too_small, TRUE),
              "out of bounds index -1 provided for col_idx_too_small relative to obs" )
expect_error( extract_var(obs, col_idx_too_big, TRUE),
              "out of bounds index 100 provided for col_idx_too_big relative to obs" )

# var_specify column names don't match

# var_specify duplicates in vec

# var_specify data vector is of wrong mode or has missing


# Wrong choices of expect_matr for wrong situation



# # Differences
# expect_identical( extract_var(obs_short, c(2, 3, 4)), 2:4 )



# # f <- function(arg) {
# #   arg_nm <- deparse( substitute(arg) )
# #   print(exists("arg"))
# # }















