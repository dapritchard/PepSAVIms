
msDat_check_valid_input <- function(mass_spec, mtoz, charge, ms_inten) {

    ## Check for missing arguments

    all_var_nm <- c("mass_spec", "mtoz", "charge", "ms_inten")
    for (var_nm in all_var_nm) {
        if (!eval(substitute(hasArg(var_nm)))) {
            stop("Must provide an argument for ", var_nm, call.=FALSE)
        }
        # Check that an object exists for provided argument
        tryCatch(get(var_nm), error = function(err) {
            err <- as.character(err)
            obj_nm <- regmatches(err, gregexpr("(?<=\')(.*?)(?=\')", err, perl=TRUE))[[1L]]
            stop("object \'", obj_nm, "\' not found for ", var_nm, call.=FALSE)
        })
    }

    ## Check mass_spec

    if (!is.matrix(mass_spec) && !is.data.frame(mass_spec)) {
        stop("mass_spec must be either a matrix or data.frame", call.=FALSE)
    }
    # The decision to not allow 0 rows is a design decision, and the reason that
    # we can't have 1 row is because that if a number is supplied then we don't
    # know if it is the value of the data or an index being supplied
    else if (NROW(mass_spec) <= 1L) {
        stop("mass_spec must have 2 or more rows", call.=FALSE)
    }
    # Don't allow 0 columns in data as a design decision
    else if (identical(NCOL(mass_spec), 0L)) {
        stop("mass_spec cannot have 0 columns", call.=FALSE)
    }

    ## Check mtoz, charge

    for (var_nm in c("mtoz", "charge")) {
        x <- get(var_nm)
        if (!is.numeric(x) && !is.character(x)) {
            stop(var_nm, " must be either of mode numeric or character", call.=FALSE)
        }
    }

    ## Check ms_inten

    if (!is.null(ms_inten) && !is.numeric(ms_inten) && !is.character(ms_inten)) {
        stop("ms_inten must be either NULL or of mode numeric or character", call.=FALSE)
    }
}
