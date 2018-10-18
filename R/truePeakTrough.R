#' true peak and trough detection
#' @param peaksTmp character vector with indices of detected peaks 
#' from findPeaks
#' @param troughsTmp character vector with indices of detected 
#' troughs from findPeaks
#' @param peakRangeTmp matrix of the peak range data with at least 
#' 3 columns (1. mass-to-charge, 2. intensity, 3. retention time)
#' @param minRes numeric minimum percentage left/right resolution
#' @param lrRes logical if true both the left and right troughs 
#' must be above the minRes else
#' the peak will be discounted. (default = FALSE i.e. if only the 
#' left or right trough is less than minRes then the peak will be retained)
#' @usage truePeakTrough(peaksTmp = NULL, troughsTmp = NULL, peakRangeTmp = 
#' NULL, minRes = 50, lrRes = FALSE)
#' @return a named numeric of both the peaks and troughs fitting the criteria.
truePeakTrough <- function(peaksTmp = NULL, troughsTmp = NULL, 
peakRangeTmp = NULL, 
minRes = 50, lrRes = FALSE) {
    minRes <- 100/50
    # identify intra-peak troughs and remove
    truePeaks <- unlist(lapply(seq_len(length(peaksTmp)), 
    function(peakInd) {
        # left trough
        nearLTrough <- which(troughsTmp - peaksTmp[peakInd] < 0)
        nearLTrough <- nearLTrough[which.min(abs(troughsTmp[nearLTrough] -
        peaksTmp[peakInd]))]
        lTroughTmp <- peakRangeTmp[peaksTmp[peakInd], 2]/
        peakRangeTmp[troughsTmp[nearLTrough],  2]
        # right trough
        nearRTrough <- which(troughsTmp - peaksTmp[peakInd] > 0)
        nearRTrough <- nearRTrough[which.min(abs(troughsTmp[nearRTrough] -
        peaksTmp[peakInd]))]
        rTroughTmp <- peakRangeTmp[peaksTmp[peakInd], 2]/
        peakRangeTmp[troughsTmp[nearRTrough],2]
        if (lrRes == FALSE) {
            truePeak <- ifelse(lTroughTmp <= minRes & rTroughTmp <=
        minRes, 0, peakInd)
        } else {
            truePeak <- ifelse(lTroughTmp <= minRes | rTroughTmp <=
        minRes, 0, peakInd)
        }
        lTroughTmp <- ifelse(lTroughTmp >= minRes, nearLTrough, 0)
        rTroughTmp <- ifelse(rTroughTmp >= minRes, nearRTrough, 0)
        resTmp <- c(lTroughTmp, truePeak, rTroughTmp)
        names(resTmp) <- c("trough", "peak", "trough")
        return(resTmp)
    }))

        truePeaks <- truePeaks[truePeaks != 0]
        return(truePeaks)
    }  # end function