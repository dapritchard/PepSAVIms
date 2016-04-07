
# ````````````````````````````````````` #
#  Create data and 'true' filterMS obj  #
# ..................................... #

# Create simple mass spec dataset ----------------------------------------------

# Create mass_spec data with 8 fractions
region <- paste0("frac", 3:4)
border <- "all"
bord_ratio <- 0.05
min_inten <- 1000
max_chg <- 7L

good_ms <- matrix( c(   25,   50, 1500, 1750,   50,   25,
                        85,   65, 2000, 2000,   20,    0,
                        40,   45, 1001, 1099,    5,    0,
                       100,    0, 2555, 2200,    0,   90 ), nrow=4, byrow=TRUE)

bad_1 <- bad_2 <- bad_3 <- bad_4 <- bad_5 <- good_ms

# Maximum outside of region
bad_1[, 1] <- 3000

# Bordering region has values with ratio greater than bord_ratio in comparison to
# the maximum value
bad_2[, 2] <- 500

# The fraction next to the maximum value must be nonzero
bad_3 <- apply(good_ms, 1, function(x) {
  maxIdx <- tail(which(x == max(x)), 1)
  replace(x, maxIdx + 1, 0)
})
bad_3 <- t(bad_3)

# Each fraction in region must have value greater than min_inten
bad_4[, 3:4] <- 999

# Create ms data.  The last block of data is to correspond to a charge value > 7.
ms <- rbind(good_ms,
            bad_1,
            bad_2,
            bad_3,
            bad_4,
            good_ms)
colnames(ms) <- paste0("frac", 1:6)

# Create mass-to-charge and charge values
mtoz_vals <- seq_len( nrow(ms) )
chg_vals <- rep(3:8, each=nrow(good_ms))


# Create 'true' filterMS object ------------------------------------------------

# The following are created by visually inspecting the data

keep_cmp_idx <- list(
  c1 = setdiff(1:24, 5:8),
  c2 = setdiff(1:24, c(5:12, c(17, 18, 20))),
  c3 = setdiff(1:24, c(8, 13:16, 20)),
  c4 = setdiff(1:24, 17:20),
  c5 = setdiff(1:24, 21:24)
)

true_cmp_by_cr <- lapply(keep_cmp_idx, function(j)
  data.frame(mtoz=mtoz_vals[j], chg=chg_vals[j])
)

true_msObj <- msDat(ms[1:4, ], mtoz_vals[1:4], chg_vals[1:4])

true_summ_info <- list(
  orig_dim   = c(24L, 6L),
  reg_nm     = region,
  bor_nm     = paste0("frac", c(1:2, 5:6)),
  border     = border,
  bord_ratio = bord_ratio,
  min_inten  = min_inten,
  max_chg    = max_chg
)

true_filterMS <- list( msObj     = true_msObj,
                       cmp_by_cr = true_cmp_by_cr,
                       summ_info = true_summ_info )
class(true_filterMS) <- c("filterMS", "msDat")


# Create some variants on the first data ---------------------------------------

# These two datasets give the same results as the original via different input.
# Specifying borders >= 2 is equivalent to specifying "all".  We do have to
# change one of the items in the summ_info.

true_filterMS_num1 <- true_filterMS
true_filterMS_num1$summ_info$border <- 2

true_filterMS_num2 <- true_filterMS
true_filterMS_num2$summ_info$border <- c(6, 7)

# Create a dataset without any successful compounds by setting the charge > 7
# for every compound

summ_info_chg_sm <- true_summ_info
summ_info_chg_sm$max_chg <- 0L

cmp_by_empty5 <- true_cmp_by_cr
cmp_by_empty5$c5 <- data.frame( mtoz = mtoz_vals[integer(0)],
                                chg  = chg_vals[integer(0)] )

true_filterMS_noCmp <- list( msObj   = NULL,
                             cmp_by_cr = cmp_by_empty5,
                             summ_info = summ_info_chg_sm )
class(true_filterMS_noCmp) <- c("filterMS", "msDat")


# Create filterMS object from function -----------------------------------------

msObj <- msDat(ms, mtoz_vals, chg_vals)
fobj <- filterMS(msObj, region, border, bord_ratio, min_inten, max_chg)

fobj_noCmp <- suppressWarnings( filterMS(msObj, region, border, bord_ratio, min_inten, max_chg=0L) )




# ``````````````````` #
#    Begin testing    #
# ................... #

context("filterMS method")

# Test for simple dataset ------------------------------------------------------

test_that("filterMS: msDat obj", {
  expect_identical( fobj$msObj$ms,   true_msObj$ms )
  expect_identical( fobj$msObj$mtoz, true_msObj$mtoz )
  expect_identical( fobj$msObj$chg,  true_msObj$chg )
  expect_identical( fobj$msObj,      true_msObj )
})

test_that("filterMS: cmp_by_cr", {
  expect_identical( fobj$cmp_by_cr$c1, true_cmp_by_cr$c1 )
  expect_identical( fobj$cmp_by_cr$c2, true_cmp_by_cr$c2 )
  expect_identical( fobj$cmp_by_cr$c3, true_cmp_by_cr$c3 )
  expect_identical( fobj$cmp_by_cr$c4, true_cmp_by_cr$c4 )
  expect_identical( fobj$cmp_by_cr$c5, true_cmp_by_cr$c5 )
  expect_identical( fobj$cmp_by_cr,    true_cmp_by_cr )
})

test_that("filterMS: summ_info", {
  expect_identical( fobj$summ_info$orig_dim,  true_summ_info$orig_dim )
  expect_identical( fobj$summ_info$reg_nm,    true_summ_info$reg_nm )
  expect_identical( fobj$summ_info$bor_nm,    true_summ_info$bor_nm )
  expect_identical( fobj$summ_info$border,    true_summ_info$border )
  expect_identical( fobj$summ_info$min_inten, true_summ_info$min_inten )
  expect_identical( fobj$summ_info$max_chg,   true_summ_info$max_chg )
})

# Includes some of the variants created earlier

test_that("filterMS: overall", {
  expect_identical( fobj,                              true_filterMS )
  expect_identical( filterMS(msObj, c("ac3", "ac4")),  true_filterMS )
  expect_identical( filterMS(msObj, c(3, 4)),          true_filterMS )
  expect_identical( filterMS(msObj, c(3, 4), 2),       true_filterMS_num1 )
  expect_identical( filterMS(msObj, c(3, 4), c(6, 7)), true_filterMS_num2 )
  expect_identical( fobj_noCmp,                        true_filterMS_noCmp )
})


# Test region and borders assignment -------------------------------------------

ms2 <- matrix(1:36, ncol=12)
ms2[13:18] <- 1500
colnames(ms2) <- paste0("ms_chr_", 1:12)
msObj2 <- msDat(ms2, 1:3, rep(3L, 3))

fobj_region_char <- filterMS(msObj2, paste0("chr_", 5:6))
fobj_region_nume <- filterMS(msObj2, 5:6)
fobj_border_all <- filterMS(msObj2, 5:6, "all")
fobj_border_none <- filterMS(msObj2, 5:6, "none")
fobj_border_num1 <- filterMS(msObj2, 5:6, 3)
fobj_border_num2 <- filterMS(msObj2, 5:6, c(3, 4))
fobj_border_0 <- filterMS(msObj2, 5:6, 0)

test_that("filterMS: region and borders", {
  expect_identical( fobj_region_char$summ_info$reg_nm, paste0("ms_chr_", 5:6) )
  expect_identical( fobj_region_nume$summ_info$reg_nm, paste0("ms_chr_", 5:6) )
  expect_identical( fobj_border_all$summ_info$bor_nm, paste0("ms_chr_", c(1:4, 7:12)) )
  expect_identical( fobj_border_all$summ_info$border, "all" )
  expect_identical( fobj_border_none$summ_info$bor_nm, character(0) )
  expect_identical( fobj_border_none$summ_info$border, "none" )
  expect_identical( fobj_border_num1$summ_info$bor_nm, paste0("ms_chr_", c(2:4, 7:9)) )
  expect_identical( fobj_border_num1$summ_info$border, 3 )
  expect_identical( fobj_border_num2$summ_info$bor_nm, paste0("ms_chr_", c(2:4, 7:10)) )
  expect_identical( fobj_border_num2$summ_info$border, c(3, 4) )
  expect_identical( fobj_border_0$summ_info$bor_nm, character(0) )
  expect_identical( fobj_border_0$summ_info$border, 0 )
})


# Test for invalid input -------------------------------------------------------

test_that("filterMS: missing input", {
  expect_error( filterMS(region=region),
                "Must provide an argument for msObj" )
  expect_error( filterMS(msObj),
                "Must provide an argument for region" )
  expect_error( filterMS(22, region),
                "msObj must be of class \"msDat\"" )
})

test_that("filterMS: invalid region", {
  expect_error( filterMS(msObj, list()),
                "region must be either of mode character or numeric" )
  expect_error( filterMS(msObj, integer(0)),
                "region must have length > 0" )
  expect_error( filterMS(msObj, character(0)),
                "region must have length > 0" )
  expect_error( filterMS(msObj, c(1, 1)),
                "region cannot have any duplicate values" )
  expect_error( filterMS(msObj, c("a", "a")),
                "region cannot have any duplicate values" )
  expect_error( filterMS(msObj, -99),
                "out of bounds value provided for region" )
  expect_error( filterMS(msObj, 1e10),
                "out of bounds value provided for region" )
  expect_error( filterMS(msObj, c("frac1", "frac_not_in_ms")),
                "name provided not in data - frac_not_in_ms element in region" )
  expect_error( filterMS(msObj, "frac"),
                "name provided had multiple matches in data - frac element in region" )
})

test_that("filterMS: invalid border", {
  expect_error( filterMS(msObj, 3:4, list()),
                "border must be an atomic vector with either mode character or mode numeric" )
  expect_error( filterMS(msObj, 3:4, "some"),
                "If border is of type character, then it must have value \"all\" or \"none\"" )
  expect_error( filterMS(msObj, 3:4, c("all", "none")),
                "If border is of type character, then it must have value \"all\" or \"none\"" )
  expect_error( filterMS(msObj, 3:4, 1:3),
                "border must have length 1 or 2" )
  expect_error( filterMS(msObj, 3:4, -3),
                "The value of border must be greater than or equal to 0" )
  expect_error( filterMS(msObj, 3:4, c(3, -9)),
                "The value of border must be greater than or equal to 0" )
})

test_that("filterMS: invalid bord_ratio", {
  expect_error( filterMS(msObj, 3:4, bord_ratio=list()), "bord_ratio must be of mode numeric" )
  expect_error( filterMS(msObj, 3:4, bord_ratio="none"), "bord_ratio must be of mode numeric" )
  expect_error( filterMS(msObj, 3:4, bord_ratio=-1), "bord_ratio must be nonnegative" )
})

test_that("filterMS: invalid min_inten", {
  expect_error( filterMS(msObj, 3:4, min_inten=list()), "min_inten must be of mode numeric" )
  expect_error( filterMS(msObj, 3:4, min_inten="no_inten"), "min_inten must be of mode numeric" )
  expect_error( filterMS(msObj, 3:4, min_inten=numeric(0)), "min_inten must be of length 1" )
  expect_error( filterMS(msObj, 3:4, min_inten=1:2), "min_inten must be of length 1" )
})

test_that("filterMS: max_chg", {
  expect_error( filterMS(msObj, 3:4, max_chg=list()), "max_chg must be of mode numeric" )
  expect_error( filterMS(msObj, 3:4, max_chg="no_inten"), "max_chg must be of mode numeric" )
  expect_error( filterMS(msObj, 3:4, max_chg=integer(0)), "max_chg must be of length 1" )
  expect_error( filterMS(msObj, 3:4, max_chg=1:2), "max_chg must be of length 1" )
})



