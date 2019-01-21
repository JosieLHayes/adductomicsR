#'output peak table from AdductQuantif object
#'
#'@param object a 'AdductQuantif' class object
#'@param outputDir character full path to a directory to output the peak to
#'default is the current working directory
#'@return a peaktable with number of rows equal to the number of adducts
#'quantified and 14 peak group information columns plus a number of columns
#'equal to the number of samples quantified. 
#'The peak table is saved as a csv file in the output directory
#'named: adductQuantif_peakList_'todays date'.csv. 
#'The peak table is also returned to
#'the R session and can be assigned to an object.
#'@usage outputPeakTable(object = NULL, outputDir = NULL)
#'@examples
#'\dontrun{
#'outputPeakTable(object=get(load(paste0(system.file("extdata", 
#'package= "adductData"), '/adductQuantResults.Rdata'))))}
#'@export
outputPeakTable <- function(object = NULL, outputDir = NULL) {
    if (is.null(object)) {
        stop("object is missing with no default")
    }
    if (!is(object,"AdductQuantif")) {
        stop("object is not an \"AdductQuantif\" class object")
    }
    if (is.null(outputDir)) {
        outputDir <- getwd()
    }
    outputDir <- paste0(outputDir, "/")
    # output table
    peakQuantTable <- object@peakQuantTable
    fileNamesTmp <- unique(peakQuantTable[, "file"])
    peakQuantDf <- data.frame(matrix(object@peakQuantTable[, "peakArea"],
    byrow = FALSE, 
    ncol = length(fileNamesTmp)), stringsAsFactors = FALSE)
    colnames(peakQuantDf) <- fileNamesTmp
    # weighted average of column names
    reqPeakColumns <- c("massAcc", "rtDev", "peakWidth", "peakArea",
    "height", "mzmed", 
    "mzmax", "mzmin", "corrRtMed", "corrRtMin", "corrRtMax", "rtmed",
    "rtmin", "rtmax")
    resMatrixTmp <- matrix(0, nrow = nrow(peakQuantDf), 
    ncol = length(reqPeakColumns))
    colnames(resMatrixTmp) <- reqPeakColumns
    for (i in seq_len(length(reqPeakColumns))) {
        # change to weighted mean with by
        nonZeroIndx <- as.numeric(peakQuantTable[, reqPeakColumns[i]]) != 0
        wMeanColTmp <- tapply(as.numeric(peakQuantTable[,
        reqPeakColumns[i]])[nonZeroIndx], 
        as.factor(peakQuantTable[, "featureName"])[nonZeroIndx], mean)
        tmpidx <- match(peakQuantTable[seq_len(nrow(resMatrixTmp)), 
        "featureName"],
        names(wMeanColTmp))
        resMatrixTmp[, i] <- as.numeric(wMeanColTmp)[tmpidx]
    }
    peakQuantDf <- data.frame(cbind(object@targTable, 
    resMatrixTmp, peakQuantDf), 
    stringsAsFactors = FALSE)
    write.csv(peakQuantDf, paste0(outputDir, "adductQuantif_peakList_", 
    Sys.Date(), 
    ".csv"), row.names = FALSE)
    return(peakQuantDf)
}  # end function