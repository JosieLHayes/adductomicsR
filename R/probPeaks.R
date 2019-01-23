#' potentially problematic peak identification
#' @param object an 'AdductQuantif' class object
#' @param nTimesMad numeric number of median absolute
#' deviations to identify potential
#' problem peaks.
#' @param metrics character string column names of metrics with
#' which to identify
#' potential problem peaks or a list with individual nTimesMad arguments
#' and with
#' list element names corresponding to column names of metrics.
#' @param ... further arguments to \code{\link{mad}}
#' @usage probPeaks(object = NULL, nTimesMad = 3,
#' metrics = c("nMadDotProdDistN",
#' "nMadSkewness", "nMadKurtosis", "nMadRtGroupDev",
#' "nMadPeakArea", "duplicates"))
#' @return 'AdductQuantif' class object
probPeaks <- function(object = NULL,
                      nTimesMad = 3,
                      metrics = c(
                          "nMadDotProdDistN",
                          "nMadSkewness",
                          "nMadKurtosis",
                          "nMadRtGroupDev",
                          "nMadPeakArea",
                          "duplicates"
                      )) {
    # error handling
    if (is.null(object)) {
        stop("argument object is missing with no default.")
    } else if (!is(object, 'AdductQuantif')) {
        stop("argument object is not AdductQuantif class object.")
    }
    allMetricColNames <- c(
        "nMadDotProdDistN",
        "nMadSkewness",
        "nMadKurtosis",
        "nMadRtGroupDev",
        "nMadPeakArea",
        "duplicates"
    )
    if (!is.list(metrics)) {
        if (!all(metrics %in% allMetricColNames)) {
            stop(
                "metrics argument name(s) must be one
                or a combination of:\n1. nMadDotProdDistN\n2. nMadSkewness\n3.
                nMadKurtosis\n4. nMadRtGroupDev\n5. nMadPeakArea\n"
            )
        }
        } else {
            if (!all(names(metrics) %in% allMetricColNames)) {
                stop(
                    "metrics list name(s) must be one or a combination of:\n1.
                    nMadDotProdDistN\n2.
                    nMadSkewness\n3. nMadKurtosis\n4. nMadRtGroupDev\n5.
                    nMadPeakArea\n6. duplicates\n"
                )
            }
            nTimesMad <- as.numeric(unlist(metrics))
            metrics <- names(metrics)
            }
    
    indxTmp <- as.numeric(peakQuantTable(object)[, "peakArea"]) != 0
    resultsTmp <- matrix(0, ncol = 11, nrow = length(indxTmp))
    colnames(resultsTmp) <- c(
        "dotProdDistN",
        "skewness",
        "kurtosis",
        "nMadRtGroupDev",
        "nMadPeakArea",
        "nMadDotProdDistN",
        "nMadSkewness",
        "nMadKurtosis",
        "nMadDuplicates",
        "possOutPeak",
        "peakOutMetrics"
    )
    
    quantTableTmp <- peakQuantTable(object)
    # remove previous results if necessary
    quantTableTmp <-
        quantTableTmp[, setdiff(colnames(quantTableTmp),
                                colnames(resultsTmp))]
    # median absolute deviation of the rt deviation
    madRtDevTmp <-
        tapply(abs(as.numeric(quantTableTmp[indxTmp, "rtDev"])),
               quantTableTmp[indxTmp,
                             "featureName"], function(feat)
                                 mad(feat))
    matchIndx <-
        match(quantTableTmp[, "featureName"], names(madRtDevTmp))
    # calculate median peak area for the group
    peakAreaMed <-
        tapply(as.numeric(quantTableTmp[indxTmp, "peakArea"]),
               quantTableTmp[indxTmp,
                             "featureName"], median)
    peakAreaMed <- peakAreaMed[matchIndx]
    # deviation of peak area from the group median
    peakAreaDev <-
        as.numeric(quantTableTmp[, "peakArea"]) - peakAreaMed
    madPeakAreaTmp <-
        tapply(as.numeric(quantTableTmp[indxTmp, "peakArea"]),
               quantTableTmp[indxTmp,
                             "featureName"], function(feat)
                                 mad(feat))
    madRtDevTmp <- madRtDevTmp[matchIndx]
    madPeakAreaTmp <- madPeakAreaTmp[matchIndx]
    
    for (i in which(indxTmp)) {
        # setTxtProgressBar(pb, i)
        peakRangeTmp <- peakIdData(object)[[i]]$peakRange
        peakIntIndx <-
            which(peakRangeTmp$troughs == "intPeakTrough")
        if (length(peakIntIndx) == 0) {
            next
        }
        peakIntTmp <-
            peakRangeTmp[min(peakIntIndx):max(peakIntIndx), ,
                         drop = FALSE]
        # calculate guassian
        apexTmp <- which(peakIntTmp$peaks == "intPeak")
        seqTmp <- seq(-4, 4, length = nrow(peakIntTmp))
        normGaussian <- dnorm(seqTmp, sd = 2)
        # calc inner product autoscaled normal and integrated peak
        normGaussian <- normGaussian / max(normGaussian) * 100
        scaledPeakIntTmp <-
            peakIntTmp$intensity / max(peakIntTmp$intensity) * 100
        # inner prod
        resultsTmp[i, "dotProdDistN"] <-
            as.vector((scaledPeakIntTmp %*%
                           normGaussian) /
                          (sqrt(sum(
                              scaledPeakIntTmp ^ 2
                          )) *
                              sqrt(sum(
                                  normGaussian ^ 2
                              ))))
        # skewness
        n <- length(scaledPeakIntTmp)
        resultsTmp[i, "skewness"] <- (sum((
            scaledPeakIntTmp -
                mean(scaledPeakIntTmp)
        ) ^ 3) / n) /
            (sum((
                scaledPeakIntTmp -
                    mean(scaledPeakIntTmp)
            ) ^ 2) / n) ^ (3 / 2)
        # kurtosis
        resultsTmp[i, "kurtosis"] <- n * sum((scaledPeakIntTmp -
                                                  mean(scaledPeakIntTmp)) ^ 4) /
            (sum((scaledPeakIntTmp - mean(scaledPeakIntTmp)
            ) ^ 2) ^ 2)
        # number of median absolute deviations from the Rt deviation
        resultsTmp[i, "nMadRtGroupDev"] <-
            abs(as.numeric(quantTableTmp[i,
                                         "rtDev"])) / madRtDevTmp[i]
        # number of median absolute deviations from the peak area
        resultsTmp[i, "nMadPeakArea"] <-
            abs(peakAreaDev[i]) / madPeakAreaTmp[i]
    }
    # madness skewness
    madSkewness <- mad(resultsTmp[indxTmp, "skewness"])
    devMedianSkewness <- ifelse(resultsTmp[, "skewness"] == 0, 0,
                                abs(abs(resultsTmp[,
                                                   "skewness"]) - median(
                                                    resultsTmp[
                                                    indxTmp, "skewness"])))
    resultsTmp[, "nMadSkewness"] <- devMedianSkewness / madSkewness
    # mad kurtosis
    madKurtosis <- mad(resultsTmp[indxTmp, "kurtosis"])
    devMedianKurtosis <- ifelse(resultsTmp[, "kurtosis"] == 0, 0,
                                abs(resultsTmp[,
                                               "kurtosis"] - median(
                                                    resultsTmp[
                                                    indxTmp, "kurtosis"])))
    resultsTmp[, "nMadKurtosis"] <- devMedianKurtosis / madKurtosis
    # mad dotprod
    madDotProd <- mad(resultsTmp[indxTmp, "dotProdDistN"])
    devMedianDotProd <-
        ifelse(resultsTmp[, "dotProdDistN"] == 0, 0,
               abs(resultsTmp[,
                              "dotProdDistN"] - median(resultsTmp[
                                indxTmp, "dotProdDistN"])))
    resultsTmp[, "nMadDotProdDistN"] <- devMedianDotProd / madDotProd
    
    # duplicate ratios
    nPeakGrs <- length(unique(quantTableTmp[, "featureName"]))
    seqTmp <- seq_len(nrow(quantTableTmp))
    nSamps <- length(seqTmp) / nPeakGrs
    if (nSamps %% 2 != 0) {
        stop(
            "Cannot compute duplicate injection peak area ratios
            as there are an odd number of samples (n=",
            nSamps,
            ").\n"
            )
    }
    splSeq <- rep(c(rep(1, nPeakGrs), rep(2, nPeakGrs)), nSamps)
    inj1 <- seqTmp[splSeq == 1]
    inj2 <- seqTmp[splSeq == 2]
    
    peakRatios <- as.numeric(quantTableTmp[inj1, "peakArea"]) /
        as.numeric(quantTableTmp[inj2,
                                 "peakArea"])
    # possible problem peaks logical
    possOutPeaksTmp <-
        apply(resultsTmp[, allMetricColNames, drop = FALSE],
              2, function(colTmp)
                  colTmp >=
                  nTimesMad)
    resultsTmp[, "possOutPeak"] <-
        (rowSums(possOutPeaksTmp[, metrics,
                                 drop = FALSE]) >=
             1) * 1
    message(
        "Summary nPeaks > ",
        nTimesMad,
        " * median absolute deviation (mad):\n\n",
        paste0(metrics, " ",
               colSums(possOutPeaksTmp[, metrics, drop = FALSE]), "\n")
    )
    
    
    par(mfrow = c(3, 2))
    for (j in seq_len(5)) {
        colNameTmp <- colnames(resultsTmp)[j]
        yTmp <- resultsTmp[indxTmp, colNameTmp]
        orderTmp <- order(yTmp)
        possProbIndx <- possOutPeaksTmp[indxTmp, j] + 1
        plot(
            yTmp[orderTmp],
            ylab = colNameTmp,
            col = possProbIndx[orderTmp],
            main = paste0(
                ifelse(colnames(possOutPeaksTmp)[j] %in%
                           metrics, "***", ""),
                colNameTmp,
                " (median = ",
                round(median(yTmp), 2),
                ")"
            ),
            col.main = ifelse(
                colnames(possOutPeaksTmp)[j] %in% metrics,
                "blue",
                "black"
            )
        )
    }
    plot(
        c(0, 1),
        c(0, 1),
        ann = FALSE,
        bty = "n",
        type = "n",
        xaxt = "n",
        yaxt = "n"
    )
    text(
        x = 0.5,
        y = 0.5,
        paste("*** selected metrics used
              to identify possible problem peaks"),
        cex = 1.3,
        col = "blue"
        )
    text(
        x = 0.5,
        y = 0.3,
        paste("peaks > ", nTimesMad,
              " * MAD of each metric"),
        cex = 1.3,
        col = "red"
    )
    message(
        "based on a maximum number of median
        absolute deviation(s) of ",
        nTimesMad,
        ":\n\n",
        sum(resultsTmp[, "possOutPeak"]),
        " peaks ",
        "(",
        round(sum(resultsTmp[,
                             "possOutPeak"]) / sum(indxTmp) * 100, 1),
        "%)",
        " of a total of ",
        sum(indxTmp),
        " peaks quantified were potentially non-gaussian,
        long tailed, or deviant from the peak group retention
        time or peak area.\n"
    )
    
    nMadColNames <-
        colnames(resultsTmp)[grep("nMad", colnames(resultsTmp))]
    resultsTmp[, "peakOutMetrics"] <-
        apply(resultsTmp[, nMadColNames], 1,
              function(rowTmp) {
                  paste0(nMadColNames[rowTmp >= nTimesMad], collapse = "; ")
              })
    peakQuantTable(object) <-
        cbind(resultsTmp[, c("possOutPeak",
                             "peakOutMetrics")], quantTableTmp,
              resultsTmp[, c(
                  "dotProdDistN",
                  "skewness",
                  "kurtosis",
                  "nMadRtGroupDev",
                  "nMadPeakArea",
                  "nMadDotProdDistN",
                  "nMadSkewness",
                  "nMadKurtosis"
              )])
    return(object)
        }  # end function