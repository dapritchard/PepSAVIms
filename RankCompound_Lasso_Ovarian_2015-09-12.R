getwd()
setwd("/Users/Christine/Documents/Hicks Lab/2015 Data")
library(lars)

# Read data
# 'msData': mass spec data
msData <- read.csv(file="2015-09-12 Compounds Ovarian.csv", row.names=1)
# 'activity": bioactivity data
activity <- read.csv(file="2015-09-12 Bioactivity Ovarian.csv", skip=1, nrows=3, row.names=1)
actAve <- apply(X=activity, MARGIN=2, FUN=mean)

# Columns containing intensity data start at 3 (cols 1-2 are identifiers)
dataCols <- seq(3, ncol(msData))
# Identifier cols for mass spec data
msIden <- msData[, 1:2]
# Transpose puts mass spec data into fraction x compound form
msData <- t( msData[, dataCols] )
# Transpose puts activity data into fraction x replicate form
activity <- t( activity )

# Shorten names to fraction numbers
substrRight <- function(x, n) {
  substr(x, nchar(x)-n+1, nchar(x))
}
row.names(msData) <- substrRight(row.names(msData), 2)
row.names(activity) <- substrRight(row.names(activity), 2)



# Get compounds in 'potential effect size order' -------------------------------

getCmpds <- function(msData, activity, actAve, region, useAve, useFirst, numRank) {
  numRep <- ncol(activity)
  
  regionIdx <- which( row.names(msData) %in% as.character(region) )
  
  if (useAve) {
    explan <- msData[regionIdx, ]
    respon <- actAve[regionIdx]
  }
  else {
    expanIdx <- regionIdx[rep(1:length(region), each=numRep)]
    explan <- msData[expanIdx, ]
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




# Create a list of ordered compound effects ------------------------------------

region <- rep(list(11:16, 17:24), each=2)
useAveBool <- rep(c(TRUE, FALSE), times=2)
ordCmpList <- vector("list", 4)


for (k in 1:length(ordCmpList)) {
  
  ordCmpList[[k]] <- getCmpds(msData, activity, actAve, region[[k]], 
                              useAveBool[k], TRUE, numRank=20)
}




# Create table of ordered compound effects -------------------------------------

tabLen <- max(sapply(X=ordCmpList, FUN=length))
ordCmpTab <- data.frame(matrix(nrow=tabLen, ncol=length(ordCmpList)))
for (j in 1:length(ordCmpList)) {
  thisEl <- ordCmpList[[j]]
  ordCmpTab[, j] <- c(thisEl, rep(NA, tabLen - length(thisEl)))
}

key <- matrix(c(rep(c(TRUE, FALSE), each=2), rep(c(TRUE, FALSE), times=2)), nrow=2, byrow=TRUE)
dimnames(key) <- list(c("Region 11-16", "Use Ave."), paste("Col", 1:4))

write.csv(x=ordCmpTab, file="/Users/Christine/Documents/Hicks Lab/2015 Data/Ordered_Compounds_Ovarian.csv", na="")






