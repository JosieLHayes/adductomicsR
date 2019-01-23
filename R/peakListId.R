#' peak list Identification
#' @param adductSpectra AdductSpec object
#' param peakList numeric vector of peak masses
#' param exPeakMass numeric internal standard peak mass
#' @param frag.delta integer delta mass accuracy difference.
#' @param minSpecEx numeric the minimum percentage of the total ion
#' current explained
#' by the internal standard fragments (default = 40). Sometime spectra are not
#' identified due to this cutoff being set too high. If unexpected datapoints
#' have been interpolated then reduce this value.
#' @param maxRtDrift numeric the maximum retention time drift (in seconds)
#' to identify MS/MS spectrum scans (default = 360).
#' param outputPlotDir character string of output directory
#' (e.g. internal standard IAA-T3 peak list = peakList= c(290.21, 403.30,
#' 516.38, 587.42, 849.40, 884.92, 958.46, 993.97, 1050.52, 1107.06, 1209.73,
#' 1337.79, 1465.85))
#' @param peakList numeric vector of peak masses
#' @param exPeakMass numeric mass of explained peak
#' @param minPeaksId numeric minimum number of peaks IDed
#' @param maxPpmDev numeric ppm deviation
#' @param allScans boolean include all scans
#' @param closestMassByFile boolean closest mass in files
#' @param outputPlotDir character string for output plot directory
#' @usage peakListId(adductSpectra = NULL, peakList = c(290.21, 403.3,
#' 516.38, 587.42, 849.4, 884.92, 958.46, 993.97, 1050.52, 1107.06,
#' 1209.73, 1337.79,
#' 1465.85), exPeakMass = 834.7769, frag.delta = 1, minPeaksId = 7,
#' minSpecEx = 50, maxRtDrift = 360, maxPpmDev = 200, allScans = TRUE,
#' closestMassByFile = TRUE, outputPlotDir = NULL)
#' @return dataframe peak list
peakListId <-
    function(adductSpectra = NULL,
             peakList = c(
                 290.21,
                 403.3,
                 516.38,
                 587.42,
                 849.4,
                 884.92,
                 958.46,
                 993.97,
                 1050.52,
                 1107.06,
                 1209.73,
                 1337.79,
                 1465.85
             ),
             exPeakMass = 834.7769,
             frag.delta = 1,
             minPeaksId = 7,
             minSpecEx = 50,
             maxRtDrift = 360,
             maxPpmDev = 200,
             allScans = TRUE,
             closestMassByFile = TRUE,
             outputPlotDir = NULL) {
        if (is.null(adductSpectra)) {
            stop("argument adductSpectra is missing with no default.")
        } else if (!is(adductSpectra, 'AdductSpec')) {
            stop("argument adductSpectra is not an AdductSpec class object.")
        }
        stopifnot(is.numeric(peakList))
        # add to parameters
        Parameters(adductSpectra)[, 'idLevel'] <-
            ifelse(allScans == FALSE,
                   "compSpec", "rawSpec")
        # if the id is on the raw scans rather than composite spectra
        if (Parameters(adductSpectra)[, 'idLevel'] == "compSpec") {
            spectraTmp <- groupMS2spec(adductSpectra)
        } else {
            spectraTmp <-
                unlist(adductMS2spec(adductSpectra), recursive = FALSE)
        }
        # mass and rts of interest
        idxTmp <- match(names(spectraTmp),
                        paste0(
                            metaData(adductSpectra)[, 'mzXMLFile'],
                            ".MS2spectra.",
                            metaData(adductSpectra)[, 'seqNum']
                        ))
        precursorMz <- metaData(adductSpectra)[, 'precursorMZ'][idxTmp]
        # all scans within maxPpm Deviation (default = 200 ppm)
        ppmIdx <-
            abs((precursorMz - exPeakMass) / exPeakMass * 1e+06) < maxPpmDev
        spectraTmp <- spectraTmp[ppmIdx]
        # extract vectors of mass and rt
        massV <- unlist(lapply(spectraTmp, function(subL)
            subL[, 1]))
        namesTmp <- unlist(mapply(
            rep,
            names(spectraTmp),
            each = vapply(spectraTmp, nrow, FUN.VALUE = numeric(1))
        ))
        names(massV) <- namesTmp
        intV <- unlist(lapply(spectraTmp, function(subL)
            subL[, 2]))
        names(intV) <- namesTmp
        peakNos <-
            unlist(lapply(spectraTmp, function(x)
                seq_len(nrow(x))))
        names(peakNos) <- namesTmp
        massOrder <- order(massV)
        intOrder <- order(intV, decreasing = TRUE)
        # sort by intensity
        massV <- massV[intOrder]
        intV <- intV[intOrder]
        namesTmp <- namesTmp[intOrder]
        peakNos <- peakNos[intOrder]
        # add order information to massV vector
        names(massV) <-
            paste0(names(massV), ";", seq_len(length(massV)))
        peakList <- sort(peakList)
        names(peakList) <- paste("peak", seq_len(length(peakList)))
        message("identifying target ion list...\n")
        
        pb <- txtProgressBar(max = length(peakList), style = 3)
        peakListMatches <- as.numeric()
        for (j in seq_len(length(peakList))) {
            setTxtProgressBar(pb, j)
            # which peakList fragments match
            indxTmp <- which(abs(massV - peakList[j]) < frag.delta)
            names(indxTmp) <- rep(names(peakList)[j], length(indxTmp))
            peakListMatches <- c(peakListMatches, indxTmp)
        }
        peakListMatches <-
            cbind(namesTmp[peakListMatches],
                  intV[peakListMatches],
                  massV[peakListMatches],
                  names(peakListMatches),
                  peakListMatches)
        sumTic <- tapply(intV, namesTmp, sum)
        # add sum of comp spectrum TICs
        peakListMatches <- cbind(peakListMatches,
                                 sumTic[match(peakListMatches[, 1],
                                              names(sumTic))])
        # add missing comp spectra if necessary
        resColNamesTmp <- c(
            "mzXMLFile",
            "seqNum",
            "name",
            "peakId",
            "peakNo",
            "precursorMz",
            "retentionTime",
            "Freq",
            "totPercentSumInt"
        )
        specPepMatchesTmp <- as.data.frame(matrix(
            "",
            ncol = length(resColNamesTmp),
            nrow = nrow(peakListMatches)
        ), stringsAsFactors = FALSE)
        colnames(specPepMatchesTmp) <- resColNamesTmp
        indxTmp <- match(peakListMatches[, 1],
                         paste0(
                             metaData(adductSpectra)[, 'mzXMLFile'],
                             ".MS2spectra.",
                             metaData(adductSpectra)[, 'seqNum']
                         ))
        specPepMatchesTmp$precursorMz <-
            metaData(adductSpectra)[, 'precursorMZ'][indxTmp]
        specPepMatchesTmp$retentionTime <-
            metaData(adductSpectra)[, 'retentionTime'][indxTmp]
        specPepMatchesTmp$seqNum <-
            metaData(adductSpectra)[, 'seqNum'][indxTmp]
        specPepMatchesTmp$mzXMLFile <- 
            metaData(adductSpectra)[, 'mzXMLFile'][indxTmp]
        specPepMatchesTmp$peakNo <-
            peakNos[as.numeric(peakListMatches[, 5])]
        peakPercentSumInt <- (as.numeric(peakListMatches[, 2]) /
                                  as.numeric(peakListMatches[,
                                                             6])) * 100
        totPercSum <-
            tapply(peakPercentSumInt, peakListMatches[, 1], sum)
        specPepMatchesTmp$totPercentSumInt <-
            as.numeric(totPercSum[match(peakListMatches[,
                    1], names(totPercSum))])
        specPepMatchesTmp$peakId <-
            names(peakList)[match(peakListMatches[, 4],
                                  names(peakList))]
        FreqTmp <-
            tapply(specPepMatchesTmp$peakId, peakListMatches[, 1],
                   function(x)
                       length(unique(x)))
        specPepMatchesTmp$Freq <-
            as.numeric(FreqTmp[match(peakListMatches[, 1],
                                     names(FreqTmp))])
        specPepMatchesTmp$name <- peakListMatches[, 1]
        
        # subset by min peaks id'ed
        specPepMatchesTmp <- specPepMatchesTmp[as.numeric(
            specPepMatchesTmp$Freq) >= minPeaksId, , drop = FALSE]
        specPepMatchesTmp <- specPepMatchesTmp[as.numeric(
            specPepMatchesTmp$totPercentSumInt) >=minSpecEx, , drop = FALSE]
        # remove extremes from median Rt
        specPepMatchesTmp <- specPepMatchesTmp[abs(
            specPepMatchesTmp$retentionTime -median(
                specPepMatchesTmp$retentionTime)) < maxRtDrift, , drop = FALSE]
        
        # if more than 1 scan id'ed per file go for closest in mass
        if (closestMassByFile == TRUE) {
            precursorMzs <- specPepMatchesTmp$precursorMz
            names(precursorMzs) <- row.names(specPepMatchesTmp)
            minAbsMassDiff <-
                tapply(precursorMzs, specPepMatchesTmp$mzXMLFile,
                       function(x)
                           names(x)[which.min(abs(x -
                                                      exPeakMass))])
            closestMassScans <-
                specPepMatchesTmp$name[row.names(specPepMatchesTmp)
                                       %in% minAbsMassDiff]
            specPepMatchesTmp <-
                specPepMatchesTmp[specPepMatchesTmp$name %in%
                                      closestMassScans,
                                  , drop = FALSE]
        } else {
            rts <- specPepMatchesTmp$retentionTime
            names(rts) <- row.names(specPepMatchesTmp)
            medRt <- median(rts)
            minAbsRtDiff <- tapply(rts, specPepMatchesTmp$mzXMLFile,
                                   function(x)
                                       names(x)[which.min(abs(x -medRt))])
            closestRtScans <-
                specPepMatchesTmp$name[row.names(specPepMatchesTmp)
                                       %in% minAbsRtDiff]
            specPepMatchesTmp <-
                specPepMatchesTmp[specPepMatchesTmp$name %in%
                                      closestRtScans, , drop = FALSE]
        }
        # if necessary output plots
        if (!is.null(outputPlotDir)) {
            fileNames <- unique(specPepMatchesTmp$name)
            message(
                "saving ",
                length(fileNames),
                " plots in the output directory:\n",
                outputPlotDir,
                "\n\n"
            )
            
            pb <- txtProgressBar(max = length(fileNames), style = 3)
            for (j in seq_len(length(fileNames))) {
                setTxtProgressBar(pb, j)
                indxTmp <-
                    which(specPepMatchesTmp$name %in% fileNames[j])
                png(paste0(outputPlotDir, fileNames[j], ".png"))
                plot(
                    spectraTmp[[fileNames[j]]],
                    xlab = "m/z",
                    ylab = "intensity",
                    xlim = c(0,
                             1600),
                    type = "h",
                    main = paste0(
                        "preMz: ",
                        round(specPepMatchesTmp$precursorMz[indxTmp[1]],
                              4),
                        " RT: ",
                        round(specPepMatchesTmp$retentionTime[indxTmp[1]] / 60,
                              2),
                        " sumIntEx: ",
                        round(specPepMatchesTmp$totPercentSumInt[indxTmp[1]],
                              2),
                        " ppmDiff: ",
                        round((specPepMatchesTmp$precursorMz[indxTmp[1]] -
                                   exPeakMass) / exPeakMass * 1e+06,
                              2
                        )
                    )
                )
                points(
                    spectraTmp[[fileNames[j]]][specPepMatchesTmp$peakNo[
                        indxTmp], , drop = FALSE],
                    type = "h",
                    col = "red",
                    lwd = 2
                )
                dev.off()
            }
            message("...DONE\n")
            
        }
        
        return(specPepMatchesTmp)
    }  # end function