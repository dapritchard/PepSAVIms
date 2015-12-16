
context("msDat constructor")


trueMatOut <- structure(list( ms   = matrix(11:25, nrow=5),
                              mtoz = 1:5,
                              chg  = 6:10 ), class="msDat")

trueDfOut <- structure(list( ms   = setNames(data.frame(matrix(11:25, nrow=5)), c("X3", "X4", "X5")),
                              mtoz = 1:5,
                              chg  = 6:10 ), class="msDat")

msMat <- matrix(1:25, nrow=5)
msDf <- data.frame(msMat)
mtoz <- 1:5
chg <- 6:10

expect_identical( msDat(msMat, 1, 2), trueMatOut )
expect_identical( msDat(msMat[, -1], mtoz, 1), trueMatOut )
expect_identical( msDat(msMat[, -2], 1, chg), trueMatOut )
expect_identical( msDat(msMat[, -c(1,2)], mtoz, chg), trueMatOut )

expect_identical( msDat(msDf, 1, 2), trueDfOut )
expect_identical( msDat(msDf[, -1], mtoz, 1), trueDfOut )
expect_identical( msDat(msDf[, -2], 1, chg), trueDfOut )
expect_identical( msDat(msDf[, -c(1,2)], mtoz, chg), trueDfOut )
expect_identical( msDat(msDf, "X1", 2), trueDfOut )
expect_identical( msDat(msDf, 1, "X2"), trueDfOut )
expect_identical( msDat(msDf, "X1", "X2"), trueDfOut )


