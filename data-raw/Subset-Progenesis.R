
# Arg specifications determined by a visual inspection of the csv file.  The
# reason for the nrows argument is b/c the last row of data is incomplete.

progenesis <- read.csv("data-raw/20150427_CLK_BAP_VO_peptide_ion_data.csv",
                       skip=2, row.names=1, nrows=30799, check.names=FALSE)


# Keep m/z, Retention time (min), Mass, and Charge columns, corresponding to the
# 1, 2, 4, and 5-th columns respectively.

# Want to subset to the "intensity, condition 1" fractions: correspond to
# columns 77-111

# Fraction 32 has versions Sample0 and Sample1.  We don't want to use Sample0,
# which corresponds to the 98-th column.

mass_spec <- progenesis[, c(1:2, 4:5, 77:97, 99:111)]
names(mass_spec)[ grepl("Sample1", names(mass_spec)) ] <- "20150207_CLK_BAP_VO_32"


# Saves mass_spec in the file data/mass_spec.rda
devtools::use_data(mass_spec)
