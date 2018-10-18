#' adductQuantif class 
#' The adductQuantif class contains a peak integral matrix,
#' peak ranges and region of integration, the isotopic distribution identified
#' for each integrated peak and the target table of peaks integrated. 
#' 
#' @include 00adductSpec-class.R
#' @slot peakQuantTable a matrix containing the peak 
#' integration results and 
#' consisting of a row for each peak
#' identified in each sample (e.g 200 samples
#' and 50 targets 200 * 50 = 10,000 rows) 
#' @slot peakIdData list of peak IDs
#' @slot predIsoDist list of predicted Iso distances
#' @slot targTable dataframe target table
#' @slot file.paths character path to file
#' @slot Parameters dataframe of specified parameters
#' @return peak integral matrix,
#' peak ranges and region of integration, the isotopic distribution identified
#' for each integrated peak and the target table of peaks integrated
#' and their corresponding MS1 scan isotopic patterns
#' @name adductQuantif-class
#' @rdname adductQuantif-class
#' @aliases show,c,adductQuantif,adductQuantif-class 
#' @exportClass adductQuantif
#' @author JL Hayes \email{jlhayes1982@gmail.com}
setClass("adductQuantif",
        representation(peakQuantTable = "matrix",
                        peakIdData = 'list',
                        predIsoDist = 'list',
                        targTable = 'data.frame',
                        file.paths = "character",
                        Parameters = "data.frame"))
#' @rdname adductQuantif-class
#' @aliases c,adductQuantif-method  
x <- NULL                  
setMethod("c", signature(x= "adductQuantif"), function(x, ...) {
                            elements = list(x, ...)
                            # error handling check if all adductSpec object
                            if (any(vapply(elements, function(ele) 
                                is(ele,'adductQuantif'),
                                FUN.VALUE=logical(1)) == FALSE)) {
                                    stop("all elements must be an 
                                    adductQuantif class object")
                                }
                                emptyAdductQuantif <- new("adductQuantif")     
                                for (i in seq_len(length(elements))) {
                                    if 
                                    (ncol(emptyAdductQuantif@peakQuantTable) 
                                    == 0) {
                                        emptyAdductQuantif@peakQuantTable <- 
                                        elements[[i]]@peakQuantTable
                                    } else {
                                        emptyAdductQuantif@peakQuantTable <- 
                    rbind(emptyAdductQuantif@peakQuantTable, 
                                    elements[[i]]@peakQuantTable)
                                    }
                                    emptyAdductQuantif@peakIdData <- 
                                    c(emptyAdductQuantif@peakIdData,
                                elements[[i]]@peakIdData)
                                    emptyAdductQuantif@predIsoDist <- 
                                    c(emptyAdductQuantif@predIsoDist,
                                elements[[i]]@predIsoDist)
                                    # keep unique
                                    emptyAdductQuantif@predIsoDist <- 
                                    emptyAdductQuantif@predIsoDist[
                        duplicated(names(emptyAdductQuantif@predIsoDist)) 
                                    ==FALSE]
                                    # target table
                                    emptyAdductQuantif@targTable <- 
                                    rbind(emptyAdductQuantif@targTable,
                                    elements[[i]]@targTable)
                                    # remove duplicates
                                    uniEntries <-
                                    apply(emptyAdductQuantif@targTable[
                                    , seq_len(3)],
                                    1, paste, collapse = "")
                                    emptyAdductQuantif@targTable <- 
                                    emptyAdductQuantif@targTable[
                                    duplicated(uniEntries) == 
                                    FALSE, , drop = FALSE]
                                    emptyAdductQuantif@file.paths <- 
                                    c(emptyAdductQuantif@file.paths, 
                                elements[[i]]@file.paths)        
                                }
                                return(emptyAdductQuantif)
                            })  # end function 