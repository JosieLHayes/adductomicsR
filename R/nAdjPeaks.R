#' remove lower intensity adjacent peaks
#' @param peaksTmp character vector with indices of
#' detected peaks from findPeaks
#' @param troughsTmp character vector with indices
#' of detected troughs from findPeaks
#' @param peakRangeTmp matrix of the peak range data
#' with at least 3 columns (1. mass-to-charge, 2. intensity, 3. retention time)
#' @usage nAdjPeaks(peaksTmp = NULL, troughsTmp = NULL, peakRangeTmp = NULL)
#' @return peaksTmp but with lower intensity adjacent
#' peaks between the same troughs removed
nAdjPeaks <- function(peaksTmp = NULL,
                      troughsTmp = NULL,
                      peakRangeTmp = NULL) {
    names(peaksTmp) <- rep("peaks", length(peaksTmp))
    names(troughsTmp) <- rep("troughs", length(troughsTmp))
    pTCombTmp <- sort(c(peaksTmp, troughsTmp))
    rlencodTmp <- rle(names(pTCombTmp))
    # multiple peaks
    multPeaksTmp <- rlencodTmp$lengths > 1 & rlencodTmp$values ==
        "peaks"
    if (any(multPeaksTmp)) {
        indxTmp <- cumsum(rlencodTmp$lengths)[multPeaksTmp]
        remPeaks <- numeric()
        for (k in seq_along(indxTmp)) {
            p1 <- (indxTmp[k] - 1)
            p2 <- (p1 + rlencodTmp$lengths[multPeaksTmp][k]) - 1
            mPindxTmp <- seq(p1, p2, 1)
            mPtmp <- pTCombTmp[mPindxTmp]
            remPeaks <-
                c(remPeaks, mPtmp[which.min(peakRangeTmp[mPtmp, 2])])
        }
        # remove any peaks necessary
        peaksTmp <- peaksTmp[-which(peaksTmp %in% remPeaks)]
    }
    return(peaksTmp)
}  # end function