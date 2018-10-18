#' integrate a peak from a peak table with peak start and peak end 
#' retention times
#' @param peakTable a table of at least 5 columns:
#' \enumerate{
#' \item mass-to-charge.
#' \item intensity
#' \item adjusted retention time
#' \item raw retention time
#' \item scan numbers
#' }
#' @param peakStart retention time for peak start (in seconds).
#' @param peakEnd retention time for peak end (in seconds).
#' @param expMass expected mass-to-charge of target.
#' @param expRt expected retention time of target (in seconds).
#' @usage peakIntegrate(peakTable = NULL, peakStart = NULL, 
#' peakEnd = NULL, expMass = NULL, 
#' expRt = NULL)
#' @return list with peak and peak table
peakIntegrate <- function(peakTable = NULL, peakStart = NULL, 
peakEnd = NULL, expMass = NULL, 
expRt = NULL) {
    if (all(c("peaks", "troughs") %in% colnames(peakTable))) {
        peakTable <- peakTable[{
            peakTable$troughs %in% "intPeakTrough"
        } == FALSE, ]
        peakTable$peaks <- NULL
        peakTable$troughs <- NULL
    }
    peakStart <- which.min(abs(peakTable[, 3] - peakStart))
    peakEnd <- which.min(abs(peakTable[, 3] - peakEnd))
    peakTable$troughs <- ""
    peakTable$peaks <- ""
    peakApIndx <- {
        peakStart:peakEnd
    }[which.max(peakTable[peakStart:peakEnd, 2])]
    peakTable$peaks[peakApIndx] <- "intPeak"

    psDf <- data.frame(matrix(0, nrow = 1, ncol = ncol(peakTable)), 
stringsAsFactors = FALSE)
    colnames(psDf) <- colnames(peakTable)
    psDf$troughs <- "intPeakTrough"
    peDf <- psDf
    psDf$AdjRetentionTime <- peakTable[peakStart, 3] - 0.01
    psDf$retentionTime <- peakTable[peakStart, 4] - 0.01
    peDf$AdjRetentionTime <- peakTable[peakEnd, 3] + 0.01
    peDf$retentionTime <- peakTable[peakEnd, 4] + 0.01
    peakTable <- rbind(peakTable[seq_len({
        peakStart - 1
        }), ], psDf, peakTable[peakStart:peakEnd, ], peDf, peakTable[{
        peakEnd + 1
        }:nrow(peakTable), ])
    # sort by adj rt
    peakTable <- peakTable[order(peakTable[, 3]), ]
    # peak Quantification peak to zero
    tmpIndx <- which(peakTable$troughs == "intPeakTrough")
    tmpIndxV <- tmpIndx[1]:tmpIndx[2]
    # which min peak start intensity or peak end intensity
    peakTime <- peakTable[tmpIndxV, 3]
    peakIntensity <- peakTable[tmpIndxV, 2]
    n <- length(tmpIndxV)
    x <- vector(mode = "numeric", length = n)
    for (k in seq_len(n)) {
        h <- (k%%n) + 1
        x[k] <- peakTime[k] * peakIntensity[h] - peakTime[h] * peakIntensity[k]
    }

    resV <- vector("numeric", 16)
    names(resV) <- c("massAcc", "rtDev", "peakWidth", "scanRange", 
    "nScans", "peakArea", 
    "height", "mzmed", "mzmax", "mzmin", "corrRtMed", "corrRtMax",
    "corrRtMin", 
    "rtmed", "rtmax", "rtmin")
    resV["peakWidth"] <- round({
        peakTable[tmpIndx[2], 3] - peakTable[tmpIndx[1], 3]
        }/60, 2)
        scRange <- peakTable[tmpIndxV, 5]
        zeroIdx <- scRange != 0
        scRange <- scRange[zeroIdx]
        resV["scanRange"] <- paste0(range(scRange), collapse = "-")
        resV["nScans"] <- length(scRange)
        resV["peakArea"] <- round(abs(sum(x)/2), 2)
        peakApIndx <- peakTable$peaks == "intPeak"
        resV["massAcc"] <- round({
            {
                peakTable[peakApIndx, 1] - expMass
            }/expMass
        } * 1e+06, 4)
            resV["rtDev"] <- round({
                peakTable[peakApIndx, 3]/60
            } - expRt/60, 2)
                resV["height"] <- round(peakTable[peakApIndx, 2], 2)
                resV["mzmed"] <- round(as.numeric(peakTable[peakApIndx, 1]), 4)
                resV["mzmax"] <- round(max(as.numeric(peakTable[tmpIndxV[
                zeroIdx], 1])), 4)
                resV["mzmin"] <- round(min(as.numeric(peakTable[tmpIndxV[
                zeroIdx], 1])), 4)
                resV["corrRtMed"] <- round(as.numeric(peakTable[
                peakApIndx, 3])/60, 2)
                resV["corrRtMax"] <- round(max(as.numeric(peakTable[
                tmpIndxV[zeroIdx],
                3])/60), 2)
                resV["corrRtMin"] <- round(min(as.numeric(peakTable[
                tmpIndxV[zeroIdx],
                3])/60), 2)
                resV["rtmed"] <- round(as.numeric(peakTable[peakApIndx, 4])/
                60, 2)
                resV["rtmax"] <- round(max(as.numeric(peakTable[tmpIndxV[
                zeroIdx], 
                4])/60), 2)
                resV["rtmin"] <- round(min(as.numeric(peakTable[tmpIndxV[
                zeroIdx], 
                4])/60), 2)
                return(list(intRes = resV, peakTable = peakTable))
            }  # end function