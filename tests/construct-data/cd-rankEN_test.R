
# `````````````````````````````````` #
#  Load a saved dataset for testing  #
# .................................. #

# Load saved simulated data.  See object sim_args in the RData file for the
# arguments used to generate the data.
load("../data/data-rankEN_sim.RData")

# msDatObj: a (200 x 50) mass spec data object
msDatObj <- testDat$msDat

# bioact: a (4 x 50) bioactivity object
bioact <- testDat$bioact

# region_ix: the indices corresponding to the region of interest used in
# simulation
region_idx <- sim_args$regIdx








# region_idx_NA <- replace(region_idx, 1, NA)



# biodf_NA_out_region <- biodf_NA_in_region <-
#   biodf_char_out_region <- biodf_char_in_region <- data.frame(bioact)
# biodf_NA_out_region[1, 1] <- NA
# biodf_char_out_region$char_col <- rep("a", nrow(bioact))
# biodf_NA_in_region[1, region_idx[1]] <- NA
# biodf_char_in_region[, region_idx[1]] <- rep("a", nrow(bioact))

# # Create a single numeric value with missing
# na_num <- numeric(1)
# na_num[1] <- NA

# na_log <- logical(1)
# na_log[1] <- NA

# # Create data restriced to region of interest
# msDat_region_only <- msDatObj[, region_idx]
# bioact_region_only <- bioact[, region_idx]

# colnamesMS(testDat$msDat) <- paste0("mass_spec", 1:50)
# colnames(testDat$bioact) <- paste0("bioact", 1:50)



# ``````````````````````` #
#  Fit elastic net model  #
# ....................... #

# Explanetory data
ms_regr <- t( msDatObj$ms[, region_idx] )

# Outcome values
bio_regr <- colMeans(bioact[, region_idx])

# Arbitrary choice of lambda (quadratic penalty parameter)
lambda <- 0.1

# Fit model
enet_fit <- elasticnet::enet(x=ms_regr, y=bio_regr, lambda=lambda)




# ```````````````````````````````` #
#  Obtain compound entrance lists  #
# ................................ #

# Obtain action index.  This is a vector of the compound indices as they enter /
# leave the model; a positive value means it entered, and a negative value means
# it exited.
actions <- unlist(enet_fit$actions)
# Last action number is spurious (it gives the number of total fractions)
actions <- actions[-length(actions)]

# Obtain correlation values of proposed compounds
comp_cor <- apply(ms_regr, 2, function(x) cor(x, bio_regr))

# Indices in the order that the compounds first entered the model
enter_idx <- unique( actions[actions > 0] )
# Indices restricted to compounds positively correlated with bioactivity
pos_idx <- enter_idx[comp_cor[enter_idx] > 0]
# First 10 indices to enter the model
allcomp_10_idx <- enter_idx[ 1:min(length(pos_idx), 10L) ]
# First 10 indices to enter the model
pos_10_idx <- pos_idx[ 1:min(length(pos_idx), 10L) ]

# Indices in the order that the compounds first entered the model for various
# specifications of pos_only and ncomp
cmpidx <- list(all    = enter_idx,
               pos    = pos_idx,
               all_10 = allcomp_10_idx,
               pos_10 = pos_10_idx)




# `````````````````````````````` #
#  Create 'true' rankEN objects  #
# .............................. #

# Computes four 'true' objects; each is a particular combination of pos_only and
# ncomp specifications

# Dimensions of the data used for analysis
data_dim  = list(reg  = nrow(ms_regr),
                 comp = ncol(ms_regr),
                 repl = nrow(bioact))

# Column names for region of interest
region_nm <- list(ms  = colnames(msDatObj$ms)[region_idx],
                  bio = colnames(bioact)[region_idx])

# Create summ_info objects for each combination of pos_only and ncomp
generate_summ_info <- function(pos_only, ncomp) {
  list(
    data_dim  = data_dim,
    region_nm = region_nm,
    lambda    = 0.1,
    pos_only  = pos_only,
    ncomp     = ncomp
  )
}
summ_info <- list(
  all    = generate_summ_info(FALSE, NULL),
  pos    = generate_summ_info(TRUE,  NULL),
  all_10 = generate_summ_info(FALSE, 10),
  pos_10 = generate_summ_info(TRUE,  10)
)

# Construct rankEN objects for each combination of pos_only and ncomp
true_rankEN <- mapply(cmpidx, summ_info, SIMPLIFY=FALSE, FUN=function(idx, info) {
  print(idx)
  print(info$pos_only)
  list(mtoz      = msDatObj$mtoz[idx],
       charge    = msDatObj$chg[idx],
       comp_cor  = comp_cor[idx],
       enet_fit  = enet_fit,
       summ_info = info)
})
# Provide names for list elements
names(true_rankEN) <- names(cmpidx)
# Equip objects with rankEN class attribute
true_rankEN <- lapply(true_rankEN, structure, class="rankEN")

# true_default <- list(mtoz      = msDatObj$mtoz[pos_idx],
#                      charge    = msDatObj$chg[pos_idx],
#                      comp_cor  = comp_cor[pos_idx],
#                      enet_fit  = enet_default,
#                      summ_info = summ_info_default)
# class(true_default) <- "rankEN"

# true_biovec_default <- true_default
# true_biovec_default$summ_info$data_dim$repl <- 1L


# true_keep_10 <- list(mtoz      = msDatObj$mtoz[pos_10_idx],
#                      charge    = msDatObj$chg[pos_10_idx],
#                      comp_cor  = comp_cor[pos_10_idx],
#                      enet_fit  = enet_default,
#                      summ_info = summ_info_pos_10)
# true_allcomp <- list(mtoz      = msDatObj$mtoz[enter_idx],
#                      charge    = msDatObj$chg[enter_idx],
#                      comp_cor  = comp_cor[enter_idx],
#                      enet_fit  = enet_default,
#                      summ_info = summ_info_allcomp)
# true_keep_10_allcomp <- list(mtoz      = msDatObj$mtoz[allcomp_10_idx],
#                              charge    = msDatObj$chg[allcomp_10_idx],
#                              comp_cor  = comp_cor[allcomp_10_idx],
#                              enet_fit  = enet_default,
#                              summ_info = summ_info_allcomp_10)
# class(true_keep_10) <- class(true_allcomp) <- class(true_keep_10_allcomp) <- "rankEN"



# Create rankLasso object using function ---------------------------------------

# A few variations of pos_only and ncomp; matrix for bioact and indices for
# region specifiers
rankEN_default <- rankEN(msDatObj, bioact, region_idx, region_idx, lambda)
rankEN_keep_10 <- rankEN(msDatObj, bioact, region_idx, region_idx, lambda, , 10L)
rankEN_allcomp <- rankEN(msDatObj, bioact, region_idx, region_idx, lambda, FALSE)
rankEN_keep_10_allcomp <- rankEN(msDatObj, bioact, region_idx, region_idx, lambda, FALSE, 10L)

# bioact as a data.frame or vector
rankEN_biodf_default <- rankEN(msDatObj, data.frame(bioact), region_idx, region_idx, lambda)
rankEN_biovec_default <- rankEN(msDatObj, colMeans(bioact), region_idx, region_idx, lambda)

# region specifiers as character vectors
rankEN_ms_reg_char_default <-
  rankEN(msDatObj, bioact, paste0("ms", 21:30), region_idx, lambda)
rankEN_bio_reg_char_default <-
  rankEN(msDatObj, bioact, region_idx, paste0("bio", 21:30), lambda)
rankEN_biovec_reg_char_default <-
  rankEN(msDatObj, colMeans(bioact), region_idx, paste0("bio", 21:30), lambda)

# region specifiers as NULL
rankEN_ms_region_null_allcomp <-
  rankEN(msDat_region_only, bioact, NULL, region_idx, lambda, FALSE)
rankEN_bio_region_null_keep10 <-
  rankEN(msDatObj, bioact_region_only, region_idx, NULL, lambda, , 10L)
rankEN_biodf_region_null_default <-
  rankEN(msDatObj, data.frame(bioact_region_only), region_idx, NULL, lambda)
rankEN_biovec_region_null_default <-
  rankEN(msDatObj, colMeans(bioact_region_only), region_idx, NULL, lambda)

# bioact as data.frame / vector + character region specifier
rankEN_biodf_region_char_default <-
  rankEN(msDatObj, data.frame(bioact), region_idx, paste0("bio", 21:30), lambda)
rankEN_biovec_region_char_default <-
  rankEN(msDatObj, colMeans(bioact), region_idx, paste0("bio", 21:30), lambda)

# bioact with missing or character outside region of interest
rankEN_NA_outside_default <-
  rankEN(msDatObj, biodf_NA_out_region, region_idx, region_idx, lambda)
rankEN_char_outside_default <-
  rankEN(msDatObj, biodf_char_out_region, region_idx, region_idx, lambda)

# region as a matrix
# rankEN_ms_region_matr_default <-
#   rankEN(msDatObj, bioact, region_idx_matr, region_idx, lambda)
# rankEN_biovec_region_matr_default <-
#   rankEN(msDatObj, colMeans(bioact), region_idx, region_idx_matr, lambda)
