
# Pretty-prints ints so that they have uniform width and are right-justified
#
# PRE: assumes that vals is of type numeric.  If floating point then the values
# are coerced to integer internally by formatC.

format_int <- function(vals) {
    format_width( formatC(vals, format="d", big.mark=",") )
}




# Pretty-prints floats so that they have uniform width and line up at the
# decimal point
#
# PRE: assumes that vals is of type numeric (call to trunc will throw an error
# otherwise)

format_float <- function(vals) {

    # Integer part of the numbers
    intpart <- formatC(trunc(vals), format="d", big.mark=",")
    # The above command loses the negative sign for a value such as -0.4 since
    # it truncates to 0
    for (i in seq_along(vals)) {
        if ((-1 < vals[i]) && (vals[i] < 0)) {
            intpart[i] <- "-0"
        }
    }
    intpart <- format_width(intpart)

    # Decimal part of the numbers
    decpart <- sapply(strsplit(as.character(vals), "\\."), function(x) {
        if (identical(length(x), 2L)) {
            return ( paste0(".", x[2L], collapse="") )
        } else {
            return ( "" )
        }
    })
    decpart <- format_width(decpart, FALSE)

    paste0(intpart, decpart)
}




# Pads the elements of strvec with blanks so that they are of uniform number of
# characters.  align_right == TRUE right justifies the elements, i.e. places the
# padding to the left, and vice-versa for align_right == FALSE.
#
# PRE: assumes strvec is of type character and align_right is of type logical
# with length 1

format_width <- function(strvec, align_right=TRUE) {

    # Num chars per element of strvec
    strvec_nchar <- nchar(strvec)
    # Num chars to be added as padding to achieve uniform width
    padlen <- max(strvec_nchar) - strvec_nchar
    # Each element is a string consisting of the number of blanks in padlen
    pad <- sapply(padlen, function(n) paste0(rep(" ", n), collapse=""))

    # Return strvec with padding added to the appropriate side
    if (align_right) {
        return( paste0(pad, strvec) )
    }
    else {
        return( paste0(strvec, pad) )
    }
}
