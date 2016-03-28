
fullDat <- matrix(rep(1:9, each=10),
                  nrow=10,
                  dimnames = list(as.character(1:10),
                                  c("mtoz",
                                    "charge",
                                    "mass",
                                    "time_pr",
                                    paste0("ms_", 1:5))))


makecol <- function(x) {
  attr(x, "names") <- as.character(1:length(x))
  x
}

ms <- matrix(rep(5:9, each=10),
             nrow=10,
             dimnames = list(as.character(1:10),
                             paste0("ms_", 1:5)))
full_no_mass <- fullDat[, setdiff(colnames(fullDat), "mass")]




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
                                makecol(rep(3L, 10)), "time_pr"),
                    ms )
  expect_identical( extract_var(NULL, ms, TRUE, makecol(rep(1L, 10)), makecol(rep(2L, 10)),
                                makecol(rep(3L, 10)), makecol(rep(3L, 10))),
                    ms )
})

# TODO: test for illegal arguments / bad data etc.







