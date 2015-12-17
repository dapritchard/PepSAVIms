
#' Ranks compounds using the Lasso path
#'
#' Returns a vector containing information for the compounds in the order in
#' which they first enter the Lasso model
#'
#' @param msDat An object of class \code{\link{msDat}} containing the mass
#'   spectrometry data and identifying information
#'
#' @param bioact Either a non-list vector, a matrix, or a data frame providing
#'   bioactivity data.  If a non-list vector, then it is assumed that each entry
#'   corresponds to a particular fraction.  If the data is 2-dimensional, then
#'   it is assumed that each column corresponds to a particular fraction, and
#'   that each row corresponds to a particular bioactivity replicate.
#'
#' @param region Either \code{NULL}, a non-list vector, or a list containing
#'   exactly two atomic vectors, providing information specifying which
#'   fractions are to be included in the Lasso model.
#'
#'   If \code{NULL}, then it is assumed that all fractions included in the data
#'   are to be used in the model.  This requires that the number of fractions in
#'   the data for the parameter passed to \code{bioact} be the same as the
#'   number of fractions in the mass spectrometry data for the paramter passed
#'   to \code{msDat}.  It further assumes that the \code{k}-th column must refer
#'   to the same fraction for both the mass spectrometry data and the
#'   bioactivity data for every \code{k}.
#'
#'   If a non-list vector then it must be either a numeric or charcter vector,
#'   such that the vector specifies which columns (and hence which fractions) to
#'   include in the model.  If the vector is numeric, then the desired columns
#'   are specified by number, and if the vector is character, then the desired
#'   columns are specified by name.  Note that this assumes that the
#'   corresponding columns in the mass spectrometry data and the bioactivity
#'   data refer to the same fractions.  In other words if the \code{k}-th column
#'   is selected, then the \code{k}-th column must refer to the same fraction
#'   for both the mass spectrometry data and the bioactivity data.
#'
#'   If a list, then it must be a list with two named non-list vectors of equal
#'   length - one is to be named \code{ms} and the other is to be named
#'   \code{bio}.  The \code{ms} vector specifies the columns (and hence
#'   fractions) to include in the model from the mass spectrometry data, either
#'   as a vector of the column numbers or the column names.  The \code{bio}
#'   vector specifies the columns (and hence fractions) to include in the model
#'   from the bioactivity data, either as a vector of the column numbers or the
#'   column names.  It is assumed that the column from the mass spectrometry
#'   data specified by the \code{k}-th value in the \code{ms} vector corresponds
#'   to the same fraction as the column specified by the \code{k}-th value in
#'   the \code{bio} vector, for each \code{k}.
#'
#' @param useAve A logical value specifying whether or not to average replicate
#'   bioactivity observations.  Ignored if only one bioactivity observation is
#'   provided.


rankLasso <- function(msDat, bioact, region=NULL, useAve=TRUE) {

  # check if msDat of class msDat

  # check that bioact is the right form

  regionIdx <- getRegionIdx()

  # Transpose of mass spec data
  ms <- t( msDat$ms )
  # Transpose of bioactivity data;
  act <- getAct()


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
