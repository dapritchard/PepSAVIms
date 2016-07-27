
# ``````````````````` #
#  Load testing data  #
# ................... #

# See construct-data/cd-binMS.R for the file used to sample / create the data
load("../data/data-binMS-summary.RData")



# ````````````````` #
#  Perform testing  #
# ................. #

context("summary function for binMS")

test_that("binMS summary: compare outputs from binMS.format", {

    # Test output for data taken from binMS testing
    expect_identical(out_v1, target_v1)

    # Change some of the numbers from binMS object used for 1st test
    expect_identical(out_v2, target_v2)
})
