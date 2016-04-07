
extractMS <- function(msObj, type="matrix") {

  # Check args are of the right type
  if (missing(msObj) {
    stop("Must provide an argument for msObj")
  }
  else if (!is.character(type)) {
    stop("type must have mode character")
  }

  # class() returns a character vector regardless of input
  class_nm <- class(msObj)
  msDatObj <- switch(class_nm,
                     binMS    = msObj$msObj,
                     filterMS = msObj$msObj,
                     msDat    = msObj,
                     stop("msObj must be an object of class ",
                          "\"binMS\", \"filterMS\", or \"msDat\""))

  # Return data in the desired form
  switch(type,
         matrix = with(msDatObj, cbind(mtoz, chg, ms)),
         msDat  = msDatObj,
         stop("type must have a value of \"matrix\" or \"msDat\""))
}




namesMS <- function(msObj) {
  # Error-checking performed in extractMS
  msDatObj <- extractMS(msObj, "msDat")
  # returns NULL when msDatObj is NULL
  colnames(msDatObj$ms)
}



