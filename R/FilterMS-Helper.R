
# filterMS_getRegionIdx <- function(region, ms) {
  
#   extract_check_valid(region, ms, "region", TRUE)
  
#   if (is.character(region)) {
#     regIdx <- extract_char_to_idx(region, ms, "region", TRUE)
#   } else if (is.numeric(region)) {
#     regIdx <- extract_num_to_idx(region, ms, "region", TRUE)
#   } else {
#     stop("Shouldn't reach here!  Please send a bug report")
#   }
  
#   return (regIdx)
# }




filterMS_border_idx <- function(border, regIdx, ms_nc) {

  if ( is.character(border) ) {
    if ( identical(border, "all") ) {
      borIdx <- setdiff(seq_len(ms_nc), regIdx)
    }
    else if ( identical(border, "none") ) {
      borIdx <- integer(0)
    }
    else {
      stop("Shouldn't reach here!  Please send a bug report")
    }
  }
  # case: border is numeric
  else {
    bsize <- as.integer(border)
    borIdx <- filterMS_border_idx_num(bsize, regIdx, ms_nc)
  }

  return (borIdx)
}




filterMS_border_idx_num <- function(bsize, regIdx, ms_nc) {

  if ( !(identical(length(bsize), 1L) || identical(length(bsize), 2L)) ) {
    stop("If border is of mode numeric then it must have length 1 or 2", call.=FALSE)
  }
  else if ( identical(length(bsize), 1L) ) {
    bsize <- rep(bsize, 2)
  }

  if (any(bsize < 0L)) {
    stop("If border is of mode numeric, then the values must be nonnegative", call.=FALSE)
  }

  # Create border index variable: borIdx
  reg_lo <- head(regIdx, 1)
  if ((reg_lo > 1L) && (bsize[1] >= 1L)) {
    bef_lo <- max(1L, reg_lo - bsize[1])
    bef_hi <- reg_lo - 1L
    bef_seq <- seq(bef_lo, bef_hi)
  }
  else {
    bef_seq <- integer(0)
  }
  reg_hi <- tail(regIdx, 1)
  if ((reg_hi < ms_nc) && (bsize[2] >= 1L)) {
    aft_lo <- reg_hi + 1L
    aft_hi <- min(reg_hi + bsize[2], ms_nc)
    aft_seq <- seq(aft_lo, aft_hi)
  }
  else {
    aft_seq <- integer(0)
  }

  borIdx <- c(bef_seq, aft_seq)
}




filterMS_check_valid <- function(msObj, region, border, bord_ratio, min_inten, max_chg) {

  ## Check msObj

  if (missing(msObj)) {
    stop("Must provide an argument for msObj", call.=FALSE)
  }
  else if (!("msDat" %in% class(msObj))) {
    stop("msObj must be of class \"msDat\"", call.=FALSE)
  }

  ## Check region

  if (missing(region)) {
    stop("Must provide an argument for region", call.=FALSE)
  }
  else if ( !(is.character(region) || is.numeric(region)) ) {
    stop("region must be either of mode character or numeric", call.=FALSE)
  }
  else if (length(region) == 0L) {
    stop("region must have length > 0", call.=FALSE)
  }

  ## Check border

  if ( !(is.character(border) || is.numeric(border)) ) {
    stop("border must be an atomic vector with either ",
         "mode character or mode numeric", call.=FALSE)
  }
  else if ( is.character(border) && !(identical(border, "all") || identical(border, "none")) ) {
    stop("If border is of type character, then it must ",
         "have value \"all\" or \"none\"", call.=FALSE)
  }
  else if ( is.numeric(border) ) {
    if ( !(identical(length(border), 1L) || identical(length(border), 2L)) ) {
      stop("border must have length 1 or 2", call.=FALSE)
    }
    else if ( any(as.integer(border) < 0L) ) {
      stop("The value of border must be greater than or equal to 0", call.=FALSE)
    }
  }

  ## Check bord_ratio

  if ( !is.numeric(bord_ratio) ) {
    stop("bord_ratio must be of mode numeric", call.=FALSE)
  }
  else if ( !identical(length(bord_ratio), 1L) ) {
    stop("bord_ratio must be of length 1", call.=FALSE)
  }
  else if (bord_ratio < 0) {
    stop("bord_ratio must be nonnegative", call.=FALSE)
  }

  ## Check min_inten

  if ( !is.numeric(min_inten) ) {
    stop("min_inten must be of mode numeric", call.=FALSE)
  }
  else if ( !identical(length(min_inten), 1L) ) {
    stop("min_inten must be of length 1", call.=FALSE)
  }


  ## Check max_chg

  if ( !is.numeric(max_chg) ) {
    stop("max_chg must be of mode numeric", call.=FALSE)
  }
  else if ( !identical(length(max_chg), 1L) ) {
    stop("max_chg must be of length 1", call.=FALSE)
  }

}



