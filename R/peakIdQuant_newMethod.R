#' Adduct Peak quant
#'
#' @param mzTmp expected mass to charge of target
#' @param rtTmp expected retention time (in minutes) of target
#' @param peakRangeRtSub matrix MS1 scans covering entire chromatographic range
#' within which to identify peaks of interest. Contains the following
#' three columns
#' column 1 = mass, column 2 = intensity, column 3 = retention time,
#' column 4 = scan number.
#' @param rtDevModel loess retention time deviation model for the file.
#' @param isoPat named numeric containing the expected mass differences
#' between isotopes for the peptide
#' of interest.
#' @param isoPatPred matrix output from the \code{\link{IsotopicDistribution}}
#' function with additional 'id' column.
#' @param minSimScore numeric minimum dot product score for consideration
#' (must be between 0-1, default = 0.96).
#' @param maxPpm numeric ppm value for EIC extraction and integration.
#' @param gaussAlpha numeric alpha value for \code{\link{smth.gaussian}}
#' of smoother package. If supplied gaussian smoothing will be performed
#' (suggested value = 16).
#' @param spikeScans numeric number of scans that constitute a spike.
#' @param minPeakHeight numeric minimum peak height, default 5000
#' @param maxRtDrift numeric maximum retention time drift, default 20 secs
#' @param showPlots boolean for whether plots should be produced
#' @param isoWindow numeric isowindow size, default 10
#' @param maxGapMs1Scan maximum MS1 scan gap, default 5
#' @param intMaxPeak boolean integrate maximum peak
#' @description peak must be at least 50% resolved from overlapping peaks.
#' i.e. the peaks trough
#' must be at least 50% of the peak apex intensity for the peak to be
#' considered sufficiently resolved.
#' @usage peakIdQuant_newMethod(mzTmp = NULL, rtTmp = NULL,
#' peakRangeRtSub = NULL, rtDevModel = NULL, isoPat = NULL,
#' isoPatPred = NULL, minSimScore = 0.96, maxPpm = 4,
#' gaussAlpha = 16, spikeScans = 2, minPeakHeight = 5000,
#' maxRtDrift = 20, showPlots = FALSE,
#' isoWindow = 10, maxGapMs1Scan = 5, intMaxPeak = FALSE)
#' @return list
peakIdQuant_newMethod <- function(mzTmp = NULL,
                                  rtTmp = NULL,
                                  peakRangeRtSub = NULL,
                                  rtDevModel = NULL,
                                  isoPat = NULL,
                                  isoPatPred = NULL,
                                  minSimScore = 0.96,
                                  maxPpm = 4,
                                  gaussAlpha = 16,
                                  spikeScans = 2,
                                  minPeakHeight = 5000,
                                  maxRtDrift = 20,
                                  showPlots = FALSE,
                                  isoWindow = 10,
                                  maxGapMs1Scan = 5,
                                  intMaxPeak = FALSE) {
    # empty results list
    resList <- list()
    minMz <- mzTmp - (mzTmp / 1e+06) * maxPpm
    maxMz <- mzTmp + (mzTmp / 1e+06) * maxPpm
    # 1. monoisotopic mass peak range detect
    peakRangeMIM <- peakRangeRtSub[peakRangeRtSub[, 1] <
                    maxMz && peakRangeRtSub[,1] > minMz, , drop = FALSE]
    if (nrow(peakRangeMIM) > 0)
    {
        # sum signal rt
        peakRangeMIM <- peakRangeSum(
            peakRangeMIM,
            spikeScans,
            rtDevModel,
            gaussAlpha = gaussAlpha,
            maxEmptyRt = maxGapMs1Scan
        )
        # id troughs and peaks
        troughsMIM <- findPeaks(-peakRangeMIM[, 2], m = 1)
        peaksMIM <- findPeaks(peakRangeMIM[, 2], m = 1)
        # 1. subset by min peak height
        peaksMIM <-
            peaksMIM[peakRangeMIM[peaksMIM, 2] >= minPeakHeight]
        # 2. subset by max Rt drift
        peaksMIM <- peaksMIM[abs(peakRangeMIM[peaksMIM, 3] -
                                     rtTmp) <= maxRtDrift]
        if (length(peaksMIM) > 0)
        {
            #if min peak index larger than min valley then add to
            # troughs
            if (min(troughsMIM) > min(peaksMIM)) {
                troughsMIM <- c(1, troughsMIM)
            }
            if (max(troughsMIM) < max(peaksMIM)) {
                troughsMIM <- c(troughsMIM, nrow(peakRangeMIM))
            }
            # 1st round # identify true well resolved peaks
            truePeaks <-
                truePeakTrough(peaksMIM, troughsMIM,
                               peakRangeMIM)
            # subset peaksTmp and troughsMIM
            trueTroughs <- unique(truePeaks[grep("trough",
                                                 names(truePeaks))])
            truePeaks <- unique(truePeaks[grep("peak",
                                               names(truePeaks))])
            
            if (length(truePeaks) > 0 &&
                length(trueTroughs) > 0)
            {
                peaksMIM <- peaksMIM[truePeaks]
                troughsMIM <- troughsMIM[trueTroughs]
                
                if (min(troughsMIM) > min(peaksMIM)) {
                    troughsMIM <- c(1, troughsMIM)
                }
                if (max(troughsMIM) < max(peaksMIM)) {
                    troughsMIM <- c(troughsMIM, nrow(peakRangeMIM))
                }
                # if more than one peak between
                #two troughs then take the most intense
                peaksMIM <-
                    nAdjPeaks(peaksMIM, troughsMIM,
                              peakRangeMIM)
                # 2nd round # identify true well resolved peaks
                truePeaks <-
                    truePeakTrough(peaksMIM, troughsMIM,
                                   peakRangeMIM,
                                   lrRes = TRUE)
                # subset peaksTmp and troughsMIM
                trueTroughs <-
                    unique(truePeaks[grep("trough",
                                          names(truePeaks))])
                truePeaks <-
                    unique(truePeaks[grep("peak",
                                          names(truePeaks))])
                if (length(truePeaks) > 0 &&
                    length(trueTroughs) > 0)
                {
                    peaksMIM <- peaksMIM[truePeaks]
                    troughsMIM <-
                        troughsMIM[trueTroughs]
                    if (showPlots == TRUE) {
                        par(mfrow = c(1, 1))
                        plot(
                            peakRangeMIM[, 3],
                            peakRangeMIM[, 2],
                            type = "l",
                            col = 1,
                            ylab = "intensity",
                            xlab = "RT (secs)",
                            ylim = c(
                                0,
                                max(peakRangeMIM[peaksMIM, 2]) *
                                    (100 / isoPatPred$percent[1] +
                                         0.5)
                            )
                        )
                        points(peakRangeMIM[peaksMIM, 3],
                               peakRangeMIM[peaksMIM,
                                            2], col = "red")
                        points(peakRangeMIM[troughsMIM, 3],
                               peakRangeMIM[troughsMIM,
                                            2], col = "blue")
                        abline(
                            v = c(
                                peakRangeMIM[peaksMIM, 3] -
                                    isoWindow,
                                peakRangeMIM[peaksMIM,
                                             3] + isoWindow
                            ),
                            lty = (seq_len(
                                length(peaksMIM)
                            )) + 1,
                            lwd = 0.1
                        )
                        abline(
                            v = peakRangeMIM[peaksMIM, 3],
                            lty = (seq_len(
                                length(peaksMIM)
                            )) +
                                1,
                            lwd = 0.1
                        )
                        abline(
                            v = c(rtTmp - maxRtDrift,
                                  rtTmp + maxRtDrift),
                            lty = 2,
                            col = "red"
                        )
                    }
                    # if min peak index larger than min valley
                    # then add to troughs
                    if (min(troughsMIM) > min(peaksMIM)) {
                        troughsMIM <- c(1, troughsMIM)
                    }
                    if (max(troughsMIM) < max(peaksMIM)) {
                        troughsMIM <- c(troughsMIM,
                                        nrow(peakRangeMIM))
                    }
                    peaksTableTmp <- cbind("iso0",
                                           peakRangeMIM[peaksMIM
                                                        , , drop = FALSE])
                    colnames(peaksTableTmp) <-
                        c("isoPat",
                          "mass",
                          "intensity",
                          "predRtLoess",
                          "RT",
                          "seqNum")
                    # dist iso table
                    distIsoTable <-
                        matrix(nrow = 0,
                               ncol = ncol(peaksTableTmp))
                    colnames(distIsoTable) <-
                        c("isoPat",
                          "mass",
                          "intensity",
                          "predRtLoess",
                          "RT",
                          "seqNum")
                    # 2. id isotope peaks (first 5)
                    for (isoIt in 2:6) {
                        peakRangeTmp <- peakRangeRtSub[peakRangeRtSub[, 1] <
                                       (maxMz + isoPat[isoIt]) &
                                       peakRangeRtSub[, 1] >
                                       (minMz + isoPat[isoIt]), , drop =
                                       FALSE]
                        # sum signal rt
                        if (nrow(peakRangeTmp) > 0) {
                            peakRangeTmp <-
                                peakRangeSum(peakRangeTmp,
                                             spikeScans,
                                             rtDevModel)
                            peaksTmp <-
                                findPeaks(peakRangeTmp[, 2], m = 1)
                            
                            if (!is.null(peaksTmp)) {
                                if (showPlots == TRUE) {
                                    points(
                                        as.numeric(peakRangeTmp[, 3]),
                                        as.numeric(peakRangeTmp[,
                                                                2]),
                                        type = "l",
                                        col =
                                            isoIt
                                    )
                                    points(peakRangeTmp[peaksTmp, 3],
                                           peakRangeTmp[peaksTmp,
                                                        2],
                                           col = "red")
                                }
                                peakRangeTmp <-
                                    cbind(names(isoPat)[isoIt],
                                          peakRangeTmp)
                                colnames(peakRangeTmp) <-
                                    c(
                                        "isoPat",
                                        "mass",
                                        "intensity",
                                        "predRtLoess",
                                        "RT",
                                        "seqNum"
                                    )
                                # rbind dist iso info
                                distIsoTable <-
                                    rbind(distIsoTable,
                                          peakRangeTmp)
                                peakRangeTmp <-
                                    peakRangeTmp[peaksTmp, ,
                                                 drop = FALSE]
                                peaksTableTmp <-
                                    rbind(peaksTableTmp,
                                          peakRangeTmp)
                            }
                        }
                    }
                    # 3. hierarchically cluster isotope peaks
                    if (nrow(peaksTableTmp) > 1) {
                        hr <- fastcluster::hclust.vector(
                            as.numeric(peaksTableTmp[,
                                                     4]),
                            method = "median",
                            members = NULL
                        )
                        # cut tree according to absolute rt
                        #difference of less than isoWindow
                        #seconds
                        rtGroups <-
                            cutree(hr, h = isoWindow)
                        nIsoDetected <-
                            table(rtGroups[duplicated(paste0(
                                peaksTableTmp[,1], "_", rtGroups)) == FALSE])
                        nIsoDetected <-
                            nIsoDetected[names(nIsoDetected)
                                %in% rtGroups[seq_len(length(peaksMIM))]]
                        # at least 3 peaks detected
                        nIsoDetected <-
                            nIsoDetected[nIsoDetected >= 3]
                        if (length(nIsoDetected) > 0)  {
                            # subset peaks table
                            indxTmp <-
                                rtGroups %in%
                                as.numeric(names(nIsoDetected))
                            rtGroups <-
                                rtGroups[indxTmp]
                            peaksTableTmp <-
                                peaksTableTmp[indxTmp, ,
                                              drop = FALSE]
                            # 4. calculate dot product
                            # similarity
                            #with predicted iso distribution
                            # empty
                            # matrix for results
                            dotProdRes <-
                                matrix(0,
                                       ncol = 9,
                                       nrow = length(unique(rtGroups)))
                            dotProdRes[, 1] <-
                                unique(rtGroups)
                            # loop through diff rt groups and
                            # calc dot prod
                            for (rtGr in seq_len(nrow(dotProdRes))) {
                                indxTmp <- which(rtGroups %in%
                                                     unique(rtGroups)[rtGr])
                                # split in to a list of each
                                # isotope
                                splitF <- 
                                    as.factor(peaksTableTmp[, 1])[indxTmp]
                                levels(splitF) <-
                                    paste0("iso", 0:5)
                                indxTmp <-
                                    split(indxTmp,
                                          splitF)
                                # replace integer(0) with zero
                                indxTmp[vapply(indxTmp,
                                               length,
                                               FUN.VALUE = 
                                                   numeric(1)) == 0] <- 0
                                # all combination of different
                                # MIM and isotope
                                #peaks
                                allComb <-
                                    do.call(expand.grid,
                                            indxTmp)
                                # calculate dot prods for all
                                # combinations
                                dotProds <-
                                    t(apply(allComb, 1,
                                            function(indxComb)
                                {
                                    nonZeroTmp <- which(indxComb != 0)
                                    isoPattPercent <-
                                        rep(0, 6)
                                    # subset intensity values
                                    isoPattInt <-
                                        as.numeric(peaksTableTmp[as.numeric(
                                            indxComb[nonZeroTmp]), 3])
                                    # relative intensities
                                    isoPattPercent[nonZeroTmp] <-
                                        (isoPattInt / max(isoPattInt)) * 100
                                    dotProd <-
                                        as.vector((
                                            isoPattPercent %*%
                                                isoPatPred$percent[seq_len(6)]
                                        ) /
                                            (sqrt(
                                                sum(isoPattPercent ^ 2)
                                            ) *
                                                sqrt(
                                                    sum(isoPatPred$percent[
                                                        seq_len(6)] ^ 2)
                                                ))
                                        )
                                    # if rel int within 10%
                                    # expected first 3 isotopes
                                    #then accept
                                    corrPattTmp <-
                                        abs(isoPatPred$percent[seq_len(6)] -
                                                isoPattPercent)[seq_len(3)]
                                    corrPattTmp <-
                                        (
                                            all(
                                                corrPattTmp < (isoPattPercent *
                                                                0.1)[
                                                                seq_len(3)]
                                            ) | all(
                                                sign(diff(
                                                    isoPattPercent
                                                ))[seq_len(2)] ==
                                                    isoPatPred$sign[2:3]
                                            )
                                        ) * 1
                                    dotProd <-
                                        c(corrPattTmp,
                                          dotProd)
                                }))
                                # concatenate results add
                                #to table
                                dotProdRes[rtGr, 2:9] <-
                                    c(
                                        ifelse(
                                            any(dotProds[,
                                                         1] == 1),
                                            1,
                                            max(dotProds[, 2])
                                        ),
                                        max(dotProds[, 2]),
                                        as.numeric(allComb[which.max(
                                            dotProds[,2]), , drop = FALSE])
                                    )
                            }  # end rtGr loop
                            
                            # 5. identify best peak match
                            #and integrate
                            #subset
                            #by min similarity score
                            indxTmp <- which(dotProdRes[, 2]
                                             >= minSimScore)
                            dotProdRes <-
                                dotProdRes[indxTmp, ,
                                           drop = FALSE]
                            if (nrow(dotProdRes) > 0) {
                                # subset by highest
                                #intensity (hk
                                #peptide) or
                                #lowest retention time
                                #difference
                                # with expected
                                if (intMaxPeak == TRUE) {
                                    indxTmp <-
                                        which.max(as.numeric(peaksTableTmp[
                                            indxTmp, 3]))
                                } else {
                                    indxTmp <- which.min(abs(
                                        as.numeric(peaksTableTmp[indxTmp,
                                                                 4]) - rtTmp
                                    ))
                                }
                                dotProdRes <-
                                    dotProdRes[indxTmp, ,
                                               drop = FALSE]
                                # subset peaksTable by
                                #retention time
                                #group
                                indxTmp <-
                                    dotProdRes[, 4:9]
                                indxTmp <-
                                    indxTmp[indxTmp != 0]
                                peaksTableTmp <-
                                    peaksTableTmp[indxTmp, ,
                                                  drop = FALSE]
                                # add column names peak
                                #range mim
                                colnames(peakRangeMIM) <-
                                    c(
                                        "mass",
                                        "intensity",
                                        "AdjRetentionTime",
                                        "retentionTime",
                                        "seqNum"
                                    )
                                # convert to data.frame
                                peakRangeMIM <-
                                    as.data.frame(peakRangeMIM,
                                                  stringsAsFactors = FALSE)
                                # id peak apex
                                peakApex <-
                                    peaksMIM[dotProdRes[, 1]]
                                peakStart <-
                                    peakRangeMIM[max(troughsMIM[troughsMIM <
                                                                    peakApex]),
                                                 "AdjRetentionTime"]
                                peakEnd <-
                                    peakRangeMIM[min(troughsMIM[troughsMIM >
                                                                    peakApex]),
                                                 "AdjRetentionTime"]
                                # integrate peak
                                peakIntRes <-
                                    peakIntegrate(peakRangeMIM,
                                                  peakStart,
                                                  peakEnd,
                                                  mzTmp,
                                                  rtTmp)
                                peakRangeMIM <-
                                    peakIntRes$peakTable
                                # results columns
                                peakApIdx <-
                                    peakRangeMIM$peaks ==
                                    "intPeak"
                                resultsTmp <-
                                    c(
                                        dotProd =
                                            as.numeric(dotProdRes[,
                                                                  3]),
                                        isoDet = paste0(
                                            c("MIM",
                                              peaksTableTmp[2:nrow(
                                                peaksTableTmp),1]),
                                            collapse = "; "
                                        ),
                                        avOrApex = "apexSpec",
                                        peakIntRes$intRes
                                    )
                                # add id for peak in
                                #peaksTableTmp
                                peaksTableTmp[, "isoPat"] <-
                                    paste0(peaksTableTmp[,
                                                         "isoPat"], "_peak")
                                # plot m/z profile
                                peaksTableTmp <-
                                    rbind(peaksTableTmp,
                                          distIsoTable)
                                colnames(peaksTableTmp) <-
                                    c(
                                        "isoPat",
                                        "mass",
                                        "intensity",
                                        "AdjRetentionTime",
                                        "retentionTime",
                                        "seqNum"
                                    )
                                resList <-
                                    list(
                                        distRange =
                                            peaksTableTmp,
                                        peakRangeMIM = peakRangeMIM,
                                        results = resultsTmp
                                    )
                            }  # if any above min dot
                            #product score
                        }  # if nIsotopes detected equal
                        #zero
                    }  # if nrow peaks greater than 1
                }  # if length true peaks > 0 and length
                #true troughs >
                #0 2nd round
            }  # if length true peaks > 0 and length true
            #troughs > 0
            #1st round
        }  # if length peaks MIM > 0
    }  # if nrow peakMIM > 0
    if (length(resList) == 0) {
        row.names(peakRangeMIM) <- NULL
        peakRangeMIM <-
            data.frame(peakRangeMIM,
                       stringsAsFactors = FALSE)
        colnames(peakRangeMIM) <-
            c("mass",
              "intensity",
              "AdjRetentionTime",
              "retentionTime",
              "seqNum")
        peakRangeMIM$peaks <- ""
        peakRangeMIM$troughs <- ""
        peakRangeMIM <-
            peakRangeMIM[, c(
                "mass",
                "intensity",
                "AdjRetentionTime",
                "retentionTime",
                "seqNum",
                "peaks",
                "troughs"
            )]
        resList <-
            list(peakRangeMIM = peakRangeMIM)
    }
    return(resList)
}  # end function