#' AdductSpec class 
#'
#' The AdductSpec class contains dynamic noise filtered composite MS/MS spectra
#' and their corresponding MS1 scan isotopic patterns. Produced by
#' AdductSpecGen() from mzXML files.
#'
#' @slot adductMS2spec list of adduct MS2 spectras
#' @slot groupMS2spec list of group MS2 spectras
#' @slot metaData dataframe of metadata from mzXML
#' @slot aaResSeqs matrix of amino acid sequences
#' @slot specPepMatches list of spectra peptide matches
#' @slot specPepCompSpec list of comp spectra peptide matches
#' @slot sumAdductType dataframe of adduct types
#' @slot Peptides dataframe of peptides under study
#' @slot rtDevModels list of rtDevModels
#' @slot targetTable dataframe target table
#' @slot file.paths character of file path
#' @slot Parameters dataframe of parameters
#' @return dynamic noise filtered composite MS/MS spectra
#' and their corresponding MS1 scan isotopic patterns
#' @section Methods:
#' \describe{
#' \item{c}{\code{signature(object = "AdductSpec")}: Concatenates the 
#' spectra information.}
#' }
#' @name AdductSpec-class
#' @rdname AdductSpec-class
#' @aliases show,c,AdductSpec-class
#' @exportClass AdductSpec
#' @exportMethod "c"
#' @author JL Hayes \email{jlhayes1982@gmail.com}
setClass("AdductSpec",
        representation(adductMS2spec = "list",
                        groupMS2spec = "list",
                        metaData = "data.frame",
                        aaResSeqs = 'matrix',
                        specPepMatches = "list",
                        specPepCompSpec = 'list',
                        sumAdductType = "data.frame",
                        Peptides = "data.frame",
                        rtDevModels = "list",
                        targetTable = "data.frame",
                        file.paths = "list",
                        Parameters = "data.frame"))
#' @rdname AdductSpec-class
#' @aliases c,AdductSpec-method
# set method concatenate
x <- NULL  
setMethod("c", signature(x = "AdductSpec"), function(x, ...){
    elements = list(x, ...)
    # error handling check if all AdductSpec object
    if(any(vapply(elements, function(ele) is(ele,'AdductSpec'),
    FUN.VALUE=logical(1)) == FALSE)){
    stop('all elements must be an AdductSpec class object')
    }
    emptyAdductSpec <- new('AdductSpec')
    # bind together results
    # do not include any group info or other information
    for (i in seq_len(length(elements))){
    emptyAdductSpec@adductMS2spec <- c(emptyAdductSpec@adductMS2spec,  
    elements[[i]]@adductMS2spec)
    metaDataTmp <- elements[[i]]@metaData
    #metaDataTmp$aboveMinPeaks <- NULL
    metaDataTmp$msPrecursor_group <- NULL
    metaDataTmp$interMSMSrtGroups <- NULL
    metaDataTmp$predRtLoess <- NULL
    metaDataTmp$MS2groupFreq <- NULL
    metaDataTmp$MS2groupFreqAbove <- NULL
    emptyAdductSpec@metaData <- rbind(emptyAdductSpec@metaData, metaDataTmp)
    emptyAdductSpec@file.paths <- c(emptyAdductSpec@file.paths, 
    elements[[i]]@file.paths)
    }
    message('Grouping, retention time correction and composite spectra
    identification must be repeated in the concatenated "AdductSpec" class 
    object...\n')
    flush.console()
    return(emptyAdductSpec)
}) # end function

setGeneric("Specfile.paths", function(object) standardGeneric("Specfile.paths"))
setMethod("Specfile.paths", 
    signature="AdductSpec",
    definition = function(object){
         return(object@file.paths)
})
     
setGeneric("Specfile.paths<-", function(object,
        value) standardGeneric("Specfile.paths<-"))    
setMethod("Specfile.paths<-","AdductSpec",function(object, value) {
    object@file.paths <- value
    if (validObject(object))
    return(object)   
})

setGeneric("adductMS2spec", function(object) standardGeneric("adductMS2spec"))
setMethod("adductMS2spec", 
    signature="AdductSpec",
    definition = function(object){
        return(object@file.paths)
    })
setGeneric("adductMS2spec<-", function(object,
        value) standardGeneric("adductMS2spec<-"))    
setMethod("adductMS2spec<-","AdductSpec",function(object, value) {
    object@file.paths <- value
    if (validObject(object))
    return(object)   
})

setGeneric("metaData", function(object) standardGeneric("metaData"))
setMethod("metaData", 
    signature="AdductSpec",
    definition = function(object){
        return(object@file.paths)
    })
setGeneric("metaData<-", function(object,
        value) standardGeneric("metaData<-"))    
setMethod("metaData<-","AdductSpec",function(object, value) {
    object@file.paths <- value
    if (validObject(object))
    return(object)   
})

setGeneric("Parameters", function(object) standardGeneric("Parameters"))
setMethod("Parameters", 
    signature="AdductSpec",
    definition = function(object){
        return(object@file.paths)
    })
setGeneric("Parameters<-", function(object,
        value) standardGeneric("Parameters<-"))    
setMethod("Parameters<-","AdductSpec",function(object, value) {
    object@file.paths <- value
    if (validObject(object))
    return(object)   
})

setGeneric("groupMS2spec", function(object) standardGeneric("groupMS2spec"))
setMethod("groupMS2spec", 
    signature="AdductSpec",
    definition = function(object){
        return(object@file.paths)
    })
setGeneric("groupMS2spec<-", function(object,
        value) standardGeneric("groupMS2spec<-"))    
setMethod("groupMS2spec<-","AdductSpec",function(object, value) {
    object@file.paths <- value
    if (validObject(object))
    return(object)   
})

