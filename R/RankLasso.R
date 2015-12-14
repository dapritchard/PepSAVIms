
#' Ranks compounds using the Lasso path
#'
#' Returns a vector containing the compounds in the order in which they first
#' enter the Lasso model
#'
#' @param msDat The \emph{mass spectrometry} data
#'
#' @param actDat The bioactivity data

rankLasso <- function(msDat, activity, actAve, region, useAve, useFirst, numRank) {
  numRep <- ncol(activity)

  regionIdx <- which( row.names(msDat) %in% as.character(region) )

  if (useAve) {
    explan <- msDat[regionIdx, ]
    respon <- actAve[regionIdx]
  }
  else {
    expanIdx <- regionIdx[rep(1:length(region), each=numRep)]
    explan <- msDat[expanIdx, ]
    respon <- c( t(activity[regionIdx, ]) )
  }

  fit <- lars(x=as.matrix(explan), y=as.matrix(respon))
  actions <- unlist(fit$actions)
  ctr <- 1
  actionLen <- length(actions)
  entryVec <- integer()

  for (k in 1:length(actions)) {
    currTerm <- actions[k]

    if (useFirst)
      theCond <- ((currTerm > 0) & !(currTerm %in% entryVec))
    else
      theCond <- ((currTerm > 0) & !(-currTerm %in% actions[k:actionLen]))

    if (theCond) {
      entryVec[ctr] <- currTerm
      ctr <- ctr + 1
      if (ctr > numRank)
        break
    }
  }

  return ( as.integer( colnames(explan)[entryVec] ) )
}
