
# ms_mat <- matrix(1:15, nrow=5L, dimnames=list(1:5, paste0("ms", 1:3)))

# msDatObj <- msDat(ms_mat, 1:5*100, rep(7, 5))
# filterObj <- structure(list(msDatObj     = msDatObj,
#                             comp_by_crit = NULL,
#                             summ_info    = NULL),
#                        class = c("filterMS", "msDat"))


# extractMS(filterObj)

# dimnames(filterObj)
# colnames(filterObj)
# row.names(filterObj)

# extractMS(filterObj[1:2, ])

# dim(filterObj)
# nrow(filterObj)
# ncol(filterObj)


# colnames(filterObj)[3] <- "xyz"


# modified <- filterObj
# modified[2, ] <- 67:69
# extractMS(modified)


