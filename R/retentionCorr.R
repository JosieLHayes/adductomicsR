#' loess-based retention time deviation correction
#'
#' @param adductSpectra AdductSpec object
#' @param smoothingSpan numeric. fixed smoothing span, argument to loess.
#' If argument is not supplied then optimal smoothing span is 
#' calculated for each file seperately.
#' @param nMissing numeric. maximum number of missing files for a 
#' MS/MS scan group to be
#' utilized in the loess retention time deviation model. 
#' Roughly 15 percent missing values is a good starting point 
#' (e.g. nMissing=10 for 68 samples).
#' @param nExtra numeric maximum number of extra scans above
#' the total number of
#' files for a MS/MS scan group to be utilized in the 
#' loess retention time deviation model.
#' If a MS/MS scan group consists of many scans far 
#' in excess of the number of files
#' then potentially MS/MS scans from large tailing peaks or
#' isobars may be erroneously
#' grouped together and used to adjust retention time incorrectly.
#' @param folds numeric. number of cross validation steps to 
#' perform in identifying
#' optimal smoothing span parameter (see: bisoreg package
#' for more details)
#' @param outputFileDir character full path to a directory 
#' to save the output images
#' @usage retentionCorr(adductSpectra = NULL, 
#' smoothingSpan = NULL, nMissing = 1, 
#' nExtra = 1, folds = 7, outputFileDir = NULL)
#' @return LOESS RT models as adductSpectra AdductSpec object
retentionCorr <- function(adductSpectra = NULL, 
smoothingSpan = NULL, nMissing = 1, 
nExtra = 1, folds = 7, outputFileDir = NULL) {

    # error handling
    if (is.null(adductSpectra)) {
        stop("argument adductSpectra is missing with no default.")
    } else if (!is(adductSpectra,'AdductSpec')) {
        stop("argument adductSpectra is not an AdductSpec class object.")
    }

    metaDataTmp <- metaData(adductSpectra)
    # single point rt drift
    if (!is.null(metaDataTmp$intStdRtDrift)) {
        metaDataTmp$retentionTime <- as.numeric(metaDataTmp$retentionTime) + 
        (as.numeric(metaDataTmp$intStdRtDrift) * -1)
    } else {
        metaDataTmp$retentionTime <- as.numeric(metaDataTmp$retentionTime)
    }

    nFiles <- length(Specfile.paths(adductSpectra))
    nFilesPerGroup <- tapply(metaDataTmp$mzXMLFile, 
    as.factor(metaDataTmp$interMSMSrtGroups), 
    function(MSMSgroup) {
        length(unique(MSMSgroup))
    })[-1]
    wellBehaved <- names(nFilesPerGroup)[which(nFilesPerGroup >=
    (nFiles - nMissing))]
    # check n extra
    nExtraScans <- table(metaDataTmp$interMSMSrtGroups)
    nExtraScans <- nExtraScans[nExtraScans < {
        nFiles + nExtra
    }]
    wellBehaved <- wellBehaved[wellBehaved %in% names(nExtraScans)]
    if (length(wellBehaved) == 0) {
        stop("Not enough well behaved MS/MS peak groups for retention 
        time correction: consider reducing the value of the nMissing 
        parameter and/or increasing the nExtra parameter.")
    }
    wellBehavedMeta <- metaDataTmp[metaDataTmp$interMSMSrtGroups %in%
    wellBehaved, , drop = FALSE]
    medianRts <- tapply(wellBehavedMeta$retentionTime, 
    wellBehavedMeta$interMSMSrtGroups, 
    median)/60
    minMaxRt <- c(min(metaDataTmp$retentionTime)/60, 
    max(metaDataTmp$retentionTime)/60)
    rtSeqTmp <- seq(minMaxRt[1], minMaxRt[2], 0.1)
    # plot deviation from RT loess
    rtDevModels(adductSpectra) <- vector("list", nFiles)
    minMaxRtDf <- matrix(0, ncol = 2, nrow = nFiles)
    rtDevAllTmp <- vector("list", nFiles)
    if (!is.null(smoothingSpan)) {
        message("calculating LOESS fit 
        (fixed smoothing span: ", smoothingSpan, ") retention drift (n=", 
        length(wellBehaved), " well-behaved retention time groups)...\n")
    } else {
        message("calculating LOESS fit (", folds, "-fold CV) 
        retention drift (n=", 
        length(wellBehaved), " well-behaved retention time groups)...\n")
    }
    flush.console()
    metaDataTmp$predRtLoess <- 0
    pb <- txtProgressBar(min = 0, max = nFiles, style = 3)
    for (i in seq_len(nFiles)) {
        setTxtProgressBar(pb, i)
        fileNameTmp <- basename(Specfile.paths(adductSpectra))[i]
        fileIndx <- metaDataTmp$mzXMLFile %in% fileNameTmp
        fileMetaTmp <- metaDataTmp[fileIndx, , drop = FALSE]
        # mean rt each MS/MS rt group
        meanRtAll <- tapply(fileMetaTmp$retentionTime, 
        as.factor(fileMetaTmp$interMSMSrtGroups), mean)
        indxWellBeTmp <- match(wellBehaved, names(meanRtAll))
        indxWellBeTmp <- indxWellBeTmp[complete.cases(indxWellBeTmp)]
        rtsFileTmp <- meanRtAll[indxWellBeTmp]/60

        # deviation of rts from median
        deviationMed <- rtsFileTmp - medianRts[complete.cases(
        match(names(medianRts),names(rtsFileTmp)))]
        # add min and max gradient rt drifts to prevent wild predictions at
        #start and finish
        deviationMed <- c(deviationMed, deviationMed[which.min(rtsFileTmp)],
        deviationMed[which.max(rtsFileTmp)])
        rtsFileTmp <- c(rtsFileTmp, minMaxRt)
        # loess model retentionTime and retentionTime deviation optimal loess
        # predict(adductSpectra@rtDevModels[[i]], newdata=rtSeqTmp)
        if (!is.null(smoothingSpan)) {
            rtDevModels(adductSpectra)[[i]] <- loess(deviationMed ~ rtsFileTmp,
            span = smoothingSpan, 
        surface = "direct")
        } else {
            rtDevModels(adductSpectra)[[i]] <- loessWrapperMod(rtsFileTmp, 
        deviationMed, folds = folds)
        }
        # predict retention time for plotting
        metaDataTmp$predRtLoess[fileIndx] <- {
            {
                metaDataTmp$retentionTime[fileIndx]/60
            } - predict(rtDevModels(adductSpectra)[[i]], newdata = 
        metaDataTmp$retentionTime[fileIndx]/60)
        } * 60
        # min/max deviation
        minMaxRtDf[i, ] <- c(min(deviationMed), max(deviationMed))
        # deviation values
        lTmp <- length(rtsFileTmp)
        rtDevAllTmp[[i]] <- cbind(c(rtsFileTmp[seq_len((lTmp - 2))] - 
        deviationMed[seq_len((lTmp - 
        2))], rtsFileTmp[(lTmp - 1):lTmp]), deviationMed)
    }
    if (!is.null(outputFileDir)) {
        png(paste0(outputFileDir, "/rtDevPlot.png"))
    }
    plot(rtSeqTmp, rep(0, length(rtSeqTmp)), xlim = c(minMaxRt[1],minMaxRt[2]), 
    ylim = c(min(minMaxRtDf[, 1]), max(minMaxRtDf[, 2])), 
    xlab = "retentionTime (min)", 
    ylab = "deviation median retentionTime (min)",
    main = paste0("retentionTime deviation (n=", 
    length(wellBehaved), " groups, max. n=", 
    nMissing, " missing files)"), 
    type = "l")
    for (i in seq_len(nFiles)) {
        lines(x = rtSeqTmp, 
        y = predict(rtDevModels(adductSpectra)[[i]], newdata = rtSeqTmp), 
    col = i + 1)
        rtDevTmp <- rtDevAllTmp[[i]]
        points(rtDevTmp[-c(nrow(rtDevTmp) - 1, nrow(rtDevTmp)), ], 
        col = i + 1, pch = 19)    
    }
    abline(h = rep(0, length(rtSeqTmp)), col = "blue")
    if (!is.null(outputFileDir)) {
        dev.off()
    }
    # plot adjusted
    wellBehavedMeta <- metaDataTmp[metaDataTmp$interMSMSrtGroups %in% 
    wellBehaved, , drop = FALSE]
    medianRts <- tapply(wellBehavedMeta$predRtLoess, 
    wellBehavedMeta$interMSMSrtGroups, median)/60
    minMaxRt <- c(min(metaDataTmp$predRtLoess)/60, 
    max(metaDataTmp$predRtLoess)/60)
    rtSeqTmp <- seq(minMaxRt[1], minMaxRt[2], 0.1)
    # deviation from loess adjusted median values
    message("calculating deviation from loess-adjusted median values.\n")
    flush.console()
    for (i in seq_len(nFiles)) {
        setTxtProgressBar(pb, i)
        fileNameTmp <- basename(Specfile.paths(adductSpectra))[i]
        fileIndx <- metaDataTmp$mzXMLFile %in% fileNameTmp
        fileMetaTmp <- metaDataTmp[fileIndx, , drop = FALSE]
        # mean rt each MS/MS rt group
        meanRtAll <- tapply(fileMetaTmp$predRtLoess, 
        as.factor(fileMetaTmp$interMSMSrtGroups), mean)
        
        indxWellBeTmp <- match(wellBehaved, names(meanRtAll))
        indxWellBeTmp <- indxWellBeTmp[complete.cases(indxWellBeTmp)]
        
        rtsFileTmp <- meanRtAll[indxWellBeTmp]/60
        
        # deviation of rts from median
        deviationMed <- rtsFileTmp - medianRts[complete.cases(
        match(names(medianRts), names(rtsFileTmp)))]
        deviationMed <- c(deviationMed, deviationMed[which.min(rtsFileTmp)],
        deviationMed[which.max(rtsFileTmp)])
        rtsFileTmp <- c(rtsFileTmp, minMaxRt)
        # min/max deviation
        minMaxRtDf[i, ] <- c(min(deviationMed), max(deviationMed))
        # deviation values
        lTmp <- length(rtsFileTmp)
        rtDevAllTmp[[i]] <- cbind(c(rtsFileTmp[seq_len((lTmp - 2))] - 
        deviationMed[seq_len((lTmp - 
        2))], rtsFileTmp[(lTmp - 1):lTmp]), deviationMed)
    }
    if (!is.null(outputFileDir)) {
        png(paste0(outputFileDir, "/adjRtPlot.png"))
    }
    plot(rtSeqTmp, rep(0, length(rtSeqTmp)), xlim = c(minMaxRt[1],minMaxRt[2]), 
    ylim = c(min(minMaxRtDf[, 1]), max(minMaxRtDf[, 2])), 
    xlab = "retentionTime (min)", 
    ylab = "deviation median retentionTime (min)", 
    main = paste0("loess-adjusted retentionTime deviation (n=", 
    length(wellBehaved), " groups, max. n=", nMissing,
    " missing files)"), 
    type = "l")
    for (i in seq_len(nFiles)) {
        rtDevTmp <- rtDevAllTmp[[i]]
        points(rtDevTmp[-c(nrow(rtDevTmp) - 1, nrow(rtDevTmp)), ], 
    col = i + 1, pch = 19)
    }
    abline(h = rep(0, length(rtSeqTmp)), col = "blue")
    if (!is.null(outputFileDir)) {
        dev.off()
    }
    return(adductSpectra)
}