
#' Ranks compounds using the Lasso path
#'
#' Returns a vector containing the compounds in the order in which they first
#' enter the Lasso model
#'
#' @param ms The mass spectrometry data
#'
#' @param act The bioactivity data


rankLasso <- function(ms, act, msIden, region, useAve=TRUE, byName=TRUE) {
  # Transpose of mass spec data; converts to a matrix
  msDat <- t( ms[, -msIden] )
  # Transpose of bioactivity data; converts to a matrix
  actDat <- t(act)
  # Number of bioactivity replicates
  numRep <- ncol(actDat)


  # Shorten names to fraction numbers
  substrRight <- function(x, n) {
    substr(x, nchar(x) - n + 1, nchar(x))
  }
  row.names(msDat) <- substrRight(row.names(msDat), 2)
  row.names(actDat) <- substrRight(row.names(actDat), 2)

  # Indexes the desired region
  regionIdx <- which( row.names(msDat) %in% as.character(region) )

  if (useAve) {
    explan <- msDat[regionIdx, ]
    actAve <- rowMeans(actDat)
    respon <- actAve[regionIdx]
  }
#   else {
#     expanIdx <- regionIdx[rep(1:length(region), each=numRep)]
#     explan <- ms[expanIdx, ]
#     respon <- c( t(activity[regionIdx, ]) )
#   }

  fit <- lars(x=as.matrix(explan), y=as.matrix(respon))
  actions <- unlist(fit$actions)
  ctr <- 1
  actionLen <- length(actions)
  entryVec <- integer()

  for (k in 1:length(actions)) {
    currTerm <- actions[k]

    theCond <- ((currTerm > 0) & !(currTerm %in% entryVec))

    if (theCond) {
      entryVec[ctr] <- currTerm
      ctr <- ctr + 1
      if (ctr > numRank)
        break
    }
  }

  return ( as.integer( colnames(explan)[entryVec] ) )
}
