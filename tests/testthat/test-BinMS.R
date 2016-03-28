
# ``````````` #
#  Load data  #
# ........... #

# Contains only the object testMS.  See Construct-Data-BinMS.R for the script
# used to create the data.

load("tests/testthat/Data-BinMS.RData")




# ```````````````````````````````````````` #
#  Calculate the binned data for test set  #
# ........................................ #

# Consider the following settings for consolidation:
#   peak retention time:  14 - 45
#   mass:  2000 - 15000
#   charge: 2 - 10
#   mass-to-charge difference:  0.05 Da
#   peak retention time difference:  1 minute
#   
# Under these settings, the data was constructed so that there are 8 bins, each
# consisting of a consecutive block of 5 observations.  I.e. the first bin
# consists of observations 1,...,5, the second bin consists of observations
# 6,...,10, etc. for the first 40 observations.  Observations 41,...,55 do not
# pass the inclusion criteria.

# Indices in data for each bin
idx <- lapply(1:8, function(j) seq(5 * (j - 1) + 1, 5 * j))
# Consolidate into bins
binDat <- t( sapply(idx, function(j) colMeans(testMS[j, ])) )
# Reorder the binned data
binDat <- binDat[order(binDat[, "mtoz"], binDat[, "time"], binDat[, "chg"]), ]

# Construct msDat object
msObj <- msDat(binDat[, "ms", drop=FALSE], binDat[, "mtoz"], binDat[, "chg"])
# Construct summary function information.  See Construct-Data-BinMS.R for where
# these values come from.
summ_info <- list(n_tot = 55,
                    n_time_pr = 50,
                    n_mass = 50,
                    n_charge = 50,
                    n_tiMaCh = 40,
                    n_binned = nbinned,
                    time_range = time_range,
                    mass_range = mass_range,
                    charge_range = charge_range,
                    mtoz_diff = mtoz_diff,
                    time_diff = time_diff)

out <- binMS(testMS, "mtoz", "chg", "mass", "time", "ms", c(14, 45), c(2e3, 15e3), c(2, 10),
             0.05, 1)
