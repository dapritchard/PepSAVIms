
fullDat <- matrix(rep(1:9, each=10),
                  nrow=10,
                  dimnames = list(as.character(1:10),
                                  c("mtoz",
                                    "charge",
                                    "mass",
                                    "time_pr",
                                    paste0("ms_", 1:5))))

fulldf <- data.frame(fullDat)
fulldf_let <- data.frame(fulldf, let=letters[1:10])


makecol <- function(x) {
  attr(x, "names") <- as.character(1:length(x))
  x
}

ms <- matrix(rep(5:9, each=10),
             nrow=10,
             dimnames = list(as.character(1:10),
                             paste0("ms_", 1:5)))

full_no_mass <- fullDat[, setdiff(colnames(fullDat), "mass")]
fulldf_no_mass <- fulldf[, setdiff(colnames(fulldf), "mass")]




# ``````````````````` #
#    Begin testing    #
# ................... #

context("filterMS method")

test_that("extract_var: var by name", {
  expect_identical( extract_var("mtoz", fullDat), makecol(rep(1L, 10)) )
  expect_identical( extract_var("charge", fullDat), makecol(rep(2L, 10)) )
  expect_identical( extract_var("mass", fullDat), makecol(rep(3L, 10)) )
  expect_identical( extract_var("time_pr", fullDat), makecol(rep(4L, 10)) )
})

test_that("extract_var: var by index", {
  expect_identical( extract_var(1, fullDat), makecol(rep(1L, 10)) )
  expect_identical( extract_var(2, fullDat), makecol(rep(2L, 10)) )
  expect_identical( extract_var(3, fullDat), makecol(rep(3L, 10)) )
  expect_identical( extract_var(4, fullDat), makecol(rep(4L, 10)) )
})

test_that("extract_var: var by data", {
  expect_identical( extract_var(makecol(rep(1L, 10)), fullDat), makecol(rep(1L, 10)) )
})

test_that("extract_var: ms", {
  expect_identical( extract_var(NULL, fullDat, TRUE, 1, 2, 3, 4), ms )
  expect_identical( extract_var(NULL, fullDat, TRUE, "mtoz", 2, "mass", 4), ms )
  expect_identical( extract_var(NULL, full_no_mass, TRUE, "mtoz", 2,
                                makecol(rep(3L, 10)), "time_pr"), ms )
  expect_identical( extract_var(NULL, ms, TRUE, makecol(rep(1L, 10)), makecol(rep(2L, 10)),
                                makecol(rep(3L, 10)), makecol(rep(3L, 10))), ms )
  expect_identical( extract_var(NULL, fullDat[, 1:5], TRUE, 1, 2, 3, 4), ms[, 1, drop=FALSE] )
  expect_identical( extract_var(5:9, fullDat, TRUE, 1, 2, 3, 4), ms )
  expect_identical( extract_var(paste0("ms_", 1:5), fullDat, TRUE, 1, 2, 3, 4), ms )
})

test_that("extract_var: data.frame", {
  expect_identical( extract_var("mtoz", fulldf_let), makecol(rep(1L, 10)) )
  expect_identical( extract_var(2, fulldf_let), makecol(rep(2L, 10)) )
  expect_identical( extract_var(3, fulldf_let), makecol(rep(3L, 10)) )
  expect_identical( extract_var(NULL, fulldf, TRUE, 1, 2, 3, 4), ms )
  expect_identical( extract_var(NULL, fulldf, TRUE, "mtoz", 2, "mass", 4), ms )
  expect_identical( extract_var(NULL, fulldf_no_mass, TRUE, "mtoz", 2,
                                makecol(rep(3L, 10)), "time_pr"), ms )
  expect_identical( extract_var(NULL, ms, TRUE, makecol(rep(1L, 10)), makecol(rep(2L, 10)),
                                makecol(rep(3L, 10)), makecol(rep(3L, 10))), ms )
  expect_identical( extract_var(NULL, fulldf[, 1:5], TRUE, 1, 2, 3, 4), ms[, 1, drop=FALSE] )
  expect_identical( extract_var(5:9, fulldf_let, TRUE, 1, 2, 3, 4), ms )
  expect_identical( extract_var(paste0("ms_", 1:5), fulldf_let, TRUE, 1, 2, 3, 4), ms )
})

test_that("extract_var: illegal args", {
  expect_error( extract_var(list(), fullDat), "is not of the right type" )
  expect_error( extract_var(c(1, 2), fullDat), "must have length 1" )
  expect_error( extract_var(c("a", "b"), fullDat), "must have length 1" )
  expect_error( extract_var("a", list()), "data must be a matrix or data.frame" )
  expect_error( extract_var(100, fullDat), "out of bounds variable value provided")
  expect_error( extract_var(c(1, 2), fullDat), "must have length 1")
  expect_error( extract_var(c("ms_1", "ms_2"), fullDat), "must have length 1")
  expect_error( extract_var(c(1, 1), fullDat, TRUE), "cannot have any duplicate values")
  expect_error( extract_var("nomatch", fullDat), "name provided not in data" )
  expect_error( extract_var(c("mtoz", "mtoz"), fullDat, TRUE), "cannot have any duplicate values")
  expect_error( extract_var("ms", fullDat), "name provided had multiple matches in data" )
  expect_error( extract_var("let", fulldf_let), "must be numeric" )
  fullDat_na <- fullDat
  fullDat_na[5, 1] <- NA
  expect_error( extract_var(1, fullDat_na), "contains missing" )
  expect_error( extract_var(NULL, fulldf_let, TRUE, "nomatch", 2, 3, 4),
                "name provided not in data" )
  expect_error( extract_var(NULL, fulldf_let, TRUE, 100, 2, 3, 4),
                "out of bounds variable value provided" )
  expect_error( extract_var(NULL, fulldf_let, TRUE, list(), 2, 3, 4),
                "var is not of the right type" )
})  








