
obs <- matrix(1:50, nrow=10, dimnames=list(1:10, paste0("var", 1:5)))

v2 <- obs[, 2]

extract_num <- c(2L, 4L)
extract_char <- paste0("r", c(2, 4))
other_num <- setdiff(1:ncol(obs), extract_num)
other_char <- paste0("r", other_num)
extract_num_NA <- replace(extract_num, 1, NA)
extract_num_too_big <- replace(extract_num, 1, 100L)
extract_num_too_small <- replace(extract_num, 1, -1L)
extract_char_no_match <- c("r2", "z4")
extract_num_zero_len <- integer(0)
extract_char_zero_len <- character(0)
extract_num_dup <- rep(2L, 2)
extract_char_dup <- rep("r2", 2)

obs_NA_out <- obs_NA_in <- obs
obs_NA_out[2, 5] <- NA
obs_NA_in[2, 4] <- NA

obs_short <- obs[1:3, ]
obs_zero_row <- obs[integer(0), ]
obs_zero_row_df <- data.frame(obs_zero_dim)
obs_one_col <- obs[, 1L]


obs_df <- obs_df_char_out <- obs_df_char_in <- data.frame(obs)
obs_df_char_out$var1 <- rep("a", 10)
obs_df_char_in$var2 <- rep("b", 10)



## extract_idx

# num
expect_identical( extract_idx(obs, 1L), 1L )
expect_identical( extract_idx(obs, setNames(1L, "xyz")), 1L )
expect_identical( extract_idx(obs, 1), 1L )
expect_identical( extract_idx(obs, 1.3), 1L )
expect_identical( extract_idx(obs, extract_num, TRUE), c(2L, 4L) )

# char
expect_identical( extract_idx(obs, "var1"), 1L )
expect_identical( extract_idx(obs, "var2", TRUE), 2L )
expect_identical( extract_idx(obs, setNames("var1", "xyz")), 1L )
expect_identical( extract_idx(obs, extract_char, TRUE), c(2L, 4L) )

# NULL
expect_identical( extract_idx(obs, NULL, TRUE), 1:5 )
expect_identical( extract_idx(obs, NULL, TRUE, other_num), c(2L, 4L) )
expect_identical( extract_idx(obs, NULL, TRUE, 1, 2.5, pi, 4L), 5L )
expect_identical( extract_idx(obs, NULL, TRUE, other_char), c(2L, 4L) )
expect_identical( extract_idx(obs, NULL, TRUE, 2, other_char), 4L )
                 
# expect_matr as TRUE for 1-vector
expect_identical( extract_idx(obs, 1L, TRUE), 1L )
# with data.frame
expect_identical( extract_idx(obs_df, extract_num, TRUE), c(2L, 4L) )
expect_identical( extract_idx(obs_df, extract_char, TRUE), c(2L, 4L) )
expect_identical( extract_idx(obs_df, NULL, TRUE, 2, other_char), 4L )



## extract_var

# num
expect_identical( extract_var(obs, 1L), obs[, 1] )
expect_identical( extract_var(obs, 1L, TRUE), obs[, 1, drop=FALSE] )
expect_identical( extract_var(obs, extract_num, TRUE), obs[, extract_num] )
expect_identical( extract_var(obs_NA_out, extract_num, TRUE), obs[, extract_num] )


# char
expect_identical( extract_var(obs, "var2"), obs[, 2] )
expect_identical( extract_var(obs, "var3", TRUE), obs[, 3, drop=FALSE] )
expect_identical( extract_var(obs, extract_char, TRUE), obs[, extract_num] )
expect_identical( extract_var(obs_NA_out, extract_char, TRUE), obs[, extract_num] )


# NULL
expect_identical( extract_var(obs, NULL, TRUE), obs )
expect_identical( extract_var(obs, NULL, TRUE, other_num), obs[, extract_num] )
expect_identical( extract_var(obs, NULL, TRUE, 1, 2.5, pi, 4L), obs[, 5, drop=FALSE] )
expect_identical( extract_var(obs, NULL, TRUE, other_char), obs[, extract_num] )
expect_identical( extract_var(obs, NULL, TRUE, 2, other_char), obs[, 4, drop=FALSE])

# with data.frame
expect_identical( extract_var(obs_df, extract_num, TRUE), obs[, extract_num] )
expect_identical( extract_var(obs_df, extract_char, TRUE), obs[, extract_num] )
expect_identical( extract_var(obs_df, NULL, TRUE, other_char), obs[, extract_num] )
expect_identical( extract_var(obs_df_char_out, extract_num, TRUE), obs[, extract_num] )
expect_identical( extract_var(obs_df_char_out, extract_char, TRUE), obs[, extract_num] )
expect_identical( extract_var(obs_df_char_out, NULL, TRUE, other_char), obs[, extract_num] )

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

# Note: since extract_idx and extract_var are not user-level functions, we
# assume that we have already checked for argument missingness / existence in
# the calling function (and hence do not check for these issues in any of the
# routines in the Extract_Var.R file).

wrong_type <- Sys.Date()

# wrong types for args
expect_error( extract_idx(wrong_type, extract_num),
              "wrong_type must be a matrix or data.frame" )
expect_error( extract_idx(obs, wrong_type),
              "wrong_type must either be NULL or either of mode numeric or mode character" )
expect_error( extract_idx(obs, extract_num, wrong_type),
              "expect_matr must be either TRUE or FALSE" )
expect_error( extract_idx(obs, NULL, TRUE, wrong_type),
              "remove_var must either be of mode numeric or mode character" )
expect_error( extract_idx(obs, NULL, TRUE, other_num, wrong_type),
              "remove_var must either be of mode numeric or mode character" )

# data_obs has 0 or 1 dimension
extract_idx(obs_zero_row, extract_num, TRUE)
extract_idx(obs_zero_row_df, extract_num, TRUE)
extract_idx(obs_one_col, extract_num, TRUE)

# note that we don't check data_obs for valid values with extract_idx because we never subset the data


# var_specify invalid indices
expect_error( extract_idx(obs, extract_num_too_small, TRUE),
              "out of bounds index -1 provided for extract_num_too_small relative to obs" )
expect_error( extract_idx(obs, extract_num_too_big, TRUE),
              "out of bounds index 100 provided for extract_num_too_big relative to obs" )

# var_specify column names don't match
extract_idx(obs, extract_char_no_match, TRUE)

# var_specify duplicates
extract_idx(obs, extract_num_dup, TRUE)
extract_idx(obs, extract_char_dup, TRUE)

# var_specify zero-length
extract_idx(obs, extract_num_zero_len)
extract_idx(obs, extract_char_zero_len)

# var_specify wrong length wrt expect_matr
extract_idx(obs, extract_num)

# var_specify NA's
extract_idx(obs, extract_num_NA)




  






# data_obs NA in specified region
extract_var(obs_NA_in, extract_num, TRUE)
extract_var(obs_df_char_in, extract_num, TRUE)

# # var_specify NA in specified region

# # var_specify invalid indices
# expect_error( extract_var(obs, extract_num_too_small, TRUE),
#               "out of bounds index -1 provided for extract_num_too_small relative to obs" )
# expect_error( extract_var(obs, extract_num_too_big, TRUE),
#               "out of bounds index 100 provided for extract_num_too_big relative to obs" )

# # var_specify column names don't match

# # var_specify duplicates in vec

# # var_specify data vector is of wrong mode or has missing


# # Wrong choices of expect_matr for wrong situation



# # Differences
# expect_identical( extract_var(obs_short, c(2, 3, 4)), 2:4 )



# # f <- function(arg) {
# #   arg_nm <- deparse( substitute(arg) )
# #   print(exists("arg"))
# # }















