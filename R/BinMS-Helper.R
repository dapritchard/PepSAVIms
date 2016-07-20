
binMS_check_valid_input <- function(mass_spec, mtoz, charge, mass, time_peak_reten, ms_inten,
                                    time_range, mass_range, charge_range, mtoz_diff, time_diff) {

    ## Check for missing arguments

    all_var_nm <- c("mass_spec", "mtoz", "charge", "mass", "time_peak_reten",
                    "ms_inten", "time_range", "mass_range", "charge_range",
                    "mtoz_diff", "time_diff")
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

    ## Check mtoz, charge, time_peak_reten

    for (var_nm in c("mtoz", "charge", "time_peak_reten")) {
        x <- get(var_nm)
        if (!is.numeric(x) && !is.character(x)) {
            stop(var_nm, " must be either of mode numeric or character", call.=FALSE)
        }
    }

    ## Check mass, ms_inten

    for (var_nm in c("mass", "ms_inten")) {
        x <- get(var_nm)
        if (!is.null(x) && !is.numeric(x) && !is.character(x)) {
            stop(var_nm, " must be either NULL or of mode numeric or character", call.=FALSE)
        }
    }

    ## Check time_range, mass_range, charge_range

    for (var_nm in c("time_range", "mass_range", "charge_range")) {
        x <- get(var_nm)
        if (!is.numeric(x)) {
            stop(var_nm, " must be of mode numeric", call.=FALSE)
        }
        else if (!identical(length(x), 2L)) {
            stop(var_nm, " must have a length of 2", call.=FALSE)
        }
        else if (anyNA(x)) {
            stop(var_nm, " cannot have any missing", call.=FALSE)
        }
        else if (isTRUE(x[1] >= x[2])) {
            stop("The values of ", var_nm, " must be in increasing order", call.=FALSE)
        }
    }

    ## Check mtoz_diff, time_diff

    for (var_nm in c("mtoz_diff", "time_diff")) {
        x <- get(var_nm)
        if (!is.numeric(x)) {
            stop(var_nm, " must be of mode numeric", call.=FALSE)
        }
        else if (!identical(length(x), 1L)) {
            stop(var_nm, " must have a length of 1", call.=FALSE)
        }
        else if (anyNA(x)) {
            stop(var_nm, " cannot have any missing", call.=FALSE)
        }
    }
}
