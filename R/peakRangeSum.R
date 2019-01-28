#' raw eic signal intensity and mass summation and spike removal.
#'
#' @param peakRange matrix consisting of 5 columns:
#' \enumerate{
#' \item mass-to-charge values
#' \item intensity
#' \item retention time (in seconds)
#' \item scan number
#' }
#' @param spikeScans numeric number of scans <= a spike.
#' Any peaks <= this value will be removed (default = 2).= FALSE
#' @param rtDevModel loess model to correct retention times.
#' @param gaussAlpha numeric alpha value for \code{\link{smth.gaussian}}
#' of smoother package. If supplied
#' gaussian smoothing will be performed (suggested value = 16).
#' @param maxEmptyRt numeric maximum size of empty retention time
#' beyond which missing values
#' will be zero-filled
#'
#' @return matrix with masses and intensities summed by
#' retention time and retention time
#' correction based on the loess model supplied, the matrix has spikes
#' removed (consecutive non-zero intensity values <= spikeScans in length),
#' empty time segments are zero filled (> 3 seconds), optionally gaussian
#' smoothed using the link{\code{smth.gaussian}} function of the smoother
#' package
#' and is also subset based on the minimum and maximum retention time windows
#' supplied (rtWin).
#' The returned matrix consists of 5 columns:
#' \enumerate{
#' \item average mass-to-charge values by unique retention time
#' in supplied peakRange table
#' \item maximum intensity values by unique retention time in supplied
#' peakRange table
#' \item loess model corrected retention times
#' \item original retention time values
#' \item scan number by unique retention time in supplied peakRange table
#' }
#' @usage peakRangeSum(peakRange = NULL, spikeScans = 2, rtDevModel = NULL,
#' gaussAlpha = NULL,
#' maxEmptyRt = 7)
peakRangeSum <-
    function(peakRange = NULL,
            spikeScans = 2,
            rtDevModel =
                NULL,
            gaussAlpha = NULL,
            maxEmptyRt = 7) {
        intTmp <- tapply(peakRange[, 2, drop = FALSE], peakRange[, 3, drop =
                                                                FALSE], max)
        massTmp <- tapply(peakRange[, 1, drop = FALSE], peakRange[, 3, drop = 
                                                                FALSE], mean)
        scansTmp <- unique(peakRange[, 4, drop = FALSE])
        peakRange <-
            cbind(massTmp, intTmp, as.numeric(names(intTmp)), scansTmp)
        # regions great than n seconds then add zeros
        emptyRts <- diff(peakRange[, 3, drop = FALSE])
        emptyRts <- c(0.5, emptyRts)
        names(emptyRts) <- peakRange[, 3, drop = FALSE]
        emptyRtsIndx <- emptyRts >= maxEmptyRt
        if (any(emptyRtsIndx)) {
            byTmp <- min(emptyRts)
            fromTmp <-
                as.numeric(names(emptyRts))[which(emptyRtsIndx) - 1] + byTmp
            toTmp <-
                as.numeric(names(emptyRts))[which(emptyRtsIndx)] - byTmp
            emptyRts <- unlist(mapply(
                seq,
                from = fromTmp,
                to = toTmp,
                MoreArgs = list(by = byTmp),
                SIMPLIFY = FALSE
            ))
            peakRange <-
                rbind(peakRange, cbind(
                    rep(0, length(emptyRts)),
                    rep(0, length(emptyRts)),
                    emptyRts,
                    rep(0, length(emptyRts))
                ))
            peakRange <-
                peakRange[order(peakRange[, 3]), , drop = FALSE]
        }
        # replace only one scans (spikes) with zero
        spikes <- rle(peakRange[, 2] == 0)
        spikesIndx <- cumsum(spikes$lengths)[spikes$lengths <=
                                                spikeScans & spikes$values 
                                            == FALSE]
        spikesIndx <- do.call(c, lapply(spikesIndx,
                                        function(spike)
                                            seq(spike - (spikeScans -
                                                            1), spike, 1)))
        peakRange <- cbind(peakRange[, 1, drop = FALSE],
                            ifelse(seq_len(nrow(peakRange)) %in% spikesIndx,
                                    0, peakRange[, 2, drop = FALSE]),
                            peakRange[, 3:4, drop = FALSE])
        
        # adjust retention times based on retention time deviation loess model
        #  if
        # available
        if (!is.null(rtDevModel)) {
            rtDevsTmp <- predict(rtDevModel, newdata = peakRange[, 3] / 60)
            peakRange <- cbind(peakRange[, seq_len(2), drop = FALSE],
                                peakRange[, 3, drop = FALSE] -
                                    (rtDevsTmp * 60),
                                peakRange[, 3, drop = FALSE],
                                peakRange[, 4, drop = FALSE])
        } else {
            peakRange <- cbind(peakRange[, seq_len(3), drop = FALSE],
                                peakRange[, 3, drop = FALSE],
                                peakRange[, 4, drop = FALSE])
        }
        peakRange <- peakRange[order(peakRange[, 3]), , drop = FALSE]
        # if desired gaussian smooth
        if (nrow(peakRange) > 1) {
            if (!is.null(gaussAlpha)) {
                peakRange[, 2] <- smoother::smth.gaussian(peakRange[, 2],
                                                            alpha = gaussAlpha,
                                                            tails = TRUE)
            }
        }
        return(peakRange)
    }  # end function