#' extract and save retention time deviation models from AdductSpec class 
#' object
#'
#' @param object an 'AdductSpec' class object or full path to a 
#' .RData file of the 'AdductSpec' object
#' @param outputDir character full path to a directory to save the
#' .RData file
#' (defaults to the current working directory if unsupplied).
#' @return save a .RData file containing the rt deviation models 
#' and returns to the
#' workspace.
#' @usage rtDevModelSave(object = NULL, outputDir = NULL)
rtDevModelSave <- function(object = NULL, outputDir = NULL) {
    # error handling
    if (is.character(object)) {
        message("loading AdductSpec .RData file...Please wait.\n")
        flush.console()
        objectName <- load(object, envir = environment())
        # if different object then assign object name eval parse
        object <- eval(parse(text = objectName))
    }
    if (!is(object,'AdductSpec')) {
        stop("object is not an \"AdductSpec\" class object")
    }
    if (is.null(outputDir)) {
        outputDir <- paste0(getwd(), "/")
    }
    # extract rtDevModels
    rtDevModels <- object@rtDevModels
    names(rtDevModels) <- basename(object@file.paths)
    rm(object)
    # save as
    save(rtDevModels, file = paste0(outputDir, ifelse(grepl("/$", outputDir),
    "", "/"), "rtDevModels.RData"))
}  # end function