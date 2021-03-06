#' Adduct quantification for adductomicsR
#'@import foreach
#'@import ade4
#'@import pastecs
#'@import DT
#'@import fpc
#'@import stats
#'@import utils
#'@import graphics
#'@import grDevices
#'@import methods
#'@import datasets
#'@import rvest
#'@import reshape2
#'@param nCores number of cores to use for analysis.
#'If NULL thenumber of cores detected will be used.
#'@param targTable is the fullpath to the target table.
#'See inst/extdata/examplePeptideTargetTable.csv for an example.
#'@param intStdRtDrift the maximum drift for the internal standard in seconds.
#'Default = NULL and therefore no RT correction is applied to the internal
#'standard.
#'@param rtDevModels is the full path to the rtDevModels.RData file
#'from rtDevModels(). default is NULL and therefore has no RT correction.
#'@param filePaths required list of mzXML files for analysis.
#'If all files are in the
#'same directory these can be accessed using
#'list.files('J:\\parentdirectory\\directoryContainingfiles',
#'pattern='.mzXML', all.files=FALSE, full.names=TRUE).
#'@param quantObject character string for filepath to an
#'AdductQuantif object to be integrated.
#'@param maxPpm numeric for the maximum parts per million to be used.
#'@param minSimScore a numeric between 0
# and 1 for peak quality stringency.
#'@param spikeScans a numeric for the number of scans that a spike
#'must be seen in for it to be integrated. Default is 2.
#'@param indivAdduct numeric vector of AdductQuantif targets to re-integrate
#'@param minPeakHeight numeric to determine the minimum height for a
#'peak to be integrated. Default
#'is set low at 100.
#'@param maxRtDrift numeric for the maximum retention time
#'drift to be considered. Default is 20.
#'@param maxRtWindow numeric in
#'seconds for the retention time window (total window will be 2 times
#'this value)
#'@param isoWindow numeric for the pepide isotope window in seconds,
#'default is 80
#'@param hkPeptide is capitalized string for the housekeeping peptide.
#'The default is 'LVNEVTEFAK' from human serum albumin.
#'@param gaussAlpha numeric for the gaussian smoothing parameter to
#'smooth the peaks. Default is 16.
#'Output is an adductQuantf object saved to the working directory
#'@description reads mzXML files from a directory, corrects RT according
#'to RT correction model and quantifies peaks.
#'@return adductQuant object
#'@examples
#'\dontrun{
#'adductQuant(nCores=4, targTable=paste0(system.file("extdata",
#'package = "adductomicsR"),'/exampletargTable2.csv'), intStdRtDrift=30,
#'rtDevModels=system.file("extdata", "rtDevModels.RData", package =
#'"adductData"),
#'filePaths=list.files(system.file("extdata", package =
#'"adductData"),pattern=".mzXML", all.files=FALSE,full.names=TRUE),
#'quantObject=NULL,indivAdduct=NULL,maxPpm=5,minSimScore=0.8,spikeScans=1,
#'minPeakHeight=100,maxRtDrift=20,maxRtWindow=240,isoWindow=80,
#'hkPeptide='LVNEVTEFAK', gaussAlpha=16)
#'}
#'@usage adductQuant(nCores = NULL, targTable = NULL,
#'intStdRtDrift = NULL, rtDevModels = NULL,
#'filePaths = NULL, quantObject = NULL, indivAdduct = NULL, maxPpm = 4,
#'minSimScore = 0.8, spikeScans = 2, minPeakHeight = 100, maxRtDrift = 20,
#'maxRtWindow = 120, isoWindow = 80,
#'hkPeptide = "LVNEVTEFAK", gaussAlpha = 16)
#'@export
adductQuant <- function(nCores = NULL,
                        targTable = NULL,
                        intStdRtDrift = NULL,
                        rtDevModels = NULL,
                        filePaths = NULL,
                        quantObject = NULL,
                        indivAdduct = NULL,
                        maxPpm = 4,
                        minSimScore = 0.8,
                        spikeScans = 2,
                        minPeakHeight = 100,
                        maxRtDrift = 20,
                        maxRtWindow = 120,
                        isoWindow = 80,
                        hkPeptide = "LVNEVTEFAK",
                        gaussAlpha = 16) {
    if (is.null(rtDevModels)) {
        stop("Please provide an rtDevModels file...\n")
    }
    if (is.character(rtDevModels)) {
        rtDevModelsDir <- dirname(rtDevModels)
        message("loading rtDevModels .RData file...Please wait.\n")
        objectName <- load(rtDevModels, envir = environment())
        # if different object then assign object name eval parse
        rtDevModels <- eval(parse(text = objectName))
    }
    # error handling
    if (is.null(quantObject)) {
        if (is.null(targTable)) {
            stop("Please provide a target table file ")
        }
        if (is.character(targTable)) {
            # read in table from character string if necessary
            targTable <- as.data.frame(data.table::fread(targTable,
                            header = TRUE), stringsAsFactors = FALSE)
        }
        if (!is.data.frame(targTable)) {
            stop("targTable is not a data.frame")
        }
        if (!all(c("mass", "RT", "peptide", "chargeState") %in%
                colnames(targTable))) {
            stop(
                "The Quantification target table must have the
                column names:
                \"mass\", \"RT\", \"peptide\", \"chargeState\"
                (the retention time must appear in minutes
                and the peptide(s) as single letter sequences)"
            )
        }
        if (is.null(filePaths)) {
            stop("a character vector of file paths to mzXML
                or mzML files must be supplied.")
        }
        # peptide isotopic distribution prediction
        indxTmp <- !duplicated(targTable$peptide)
        pepSeqs <- as.character(targTable$peptide[indxTmp])
        eleForm <-
            lapply(pepSeqs, OrgMassSpecR::ConvertPeptide,
                   IAA = FALSE)
        eleForm <-
            lapply(seq_along(eleForm), function(ele) {
                eleForm[[ele]]$H <- eleForm[[ele]]$H +
                    targTable$chargeState[indxTmp][ele]
                return(eleForm[[ele]])
            })
        predIsoDist <-
            mapply(IsotopicDistributionMod,
                    eleForm,
                    targTable$chargeState[indxTmp],
                    SIMPLIFY = FALSE)
        predIsoDist <- lapply(predIsoDist,
                            function(x)
                                data.frame(
                                    x,
                                    id = c("MIM",
                                            paste0("iso", seq_len((
                                                nrow(x) - 1
                                            )))),
                                    sign = c(0, sign(diff(x$percent)))
                                ))
        names(predIsoDist) <- pepSeqs
        # min/ max rt diff
        minMaxRt <- c((min(as.numeric(
            targTable$RT
        )) * 60) -
            maxRtWindow,
        (max(as.numeric(
            targTable$RT
        )) *
            60) + maxRtWindow)
        # empty list for results
        results <-
            matrix(0,
                ncol = 24,
                nrow = length(filePaths)
                * nrow(targTable))
        colnames(results) <-
            c(
                "file",
                "featureName",
                "expMass",
                "expRt",
                "peptide",
                "dotProd",
                "isoDet",
                "avOrApex",
                "massAcc",
                "rtDev",
                "peakWidth",
                "scanRange",
                "nScans",
                "peakArea",
                "height",
                "mzmed",
                "mzmax",
                "mzmin",
                "corrRtMed",
                "corrRtMax",
                "corrRtMin",
                "rtmed",
                "rtmin",
                "rtmax"
            )
        # add file and feature names to table
        results[, "file"] <- rep(basename(filePaths),
                                each = nrow(targTable))
        results[, "featureName"] <- rep(paste0(
            "M",
            round(targTable$mass, 3),
            "_RT",
            round(targTable$RT, 2)
        ),
        length(filePaths))
        results[, "expMass"] <-
            rep(targTable$mass, length(filePaths))
        results[, "expRt"] <-
            rep(targTable$RT, length(filePaths))
        results[, "peptide"] <- rep(targTable$peptide,
                                    length(filePaths))
        resultsList <- vector("list", nrow(results))
        # integration sequence
        intSeq <- seq_len(nrow(targTable))
        } else {
            # error handling
            if (!is(quantObject, 'AdductQuantif')) {
                stop("argument object is not an
                    AdductQuantif class object")
            } else if (is.null(indivAdduct)) {
                stop("a numeric vector of AdductQuantif
                    targets to re-integrate must be provided")
            }
            targTable <- targTable(quantObject)
            intSeq <- indivAdduct
            minMaxRt <- c((min(as.numeric(
                targTable$RT
            )) * 60) -
                maxRtWindow, (max(as.numeric(
                    targTable$RT
                )) * 60) +
                maxRtWindow)
            # extract previous results
            predIsoDist <- predIsoDist(quantObject)
            results <- peakQuantTable(quantObject)
            # empty the results rows necessary
            keepCols <-
                grepl("file|featureName|expMass|expRt|peptide",
                    colnames(results))
            keepCols <- which(keepCols == FALSE)
            emptyResultRows <- rep(((seq_len(
                length(file.paths(quantObject))
            ) - 1)
            * nrow(targTable)),
            length(indivAdduct)) + rep(indivAdduct,
                                    each = length(file.paths(quantObject)))
            results[emptyResultRows, keepCols] <- 0
            resultsList <- peakIdData(quantObject)
            filePaths <- file.paths(quantObject)
            }
    for (i in seq_along(filePaths)) {
        filePathTmp <- filePaths[i]
        # mass drift correction
        massDriftFile <- gsub("\\.mzXML", ".massDrift.csv",
                            filePathTmp)
        if (file.exists(massDriftFile)) {
            massDriftTable <- read.csv(massDriftFile,
                                    header = TRUE,
                                    stringsAsFactors = FALSE)
        }
        MSfile <- mzR::openMSfile(filePathTmp)
        metaSubTmp <- mzR::header(MSfile)
        # max gap ms1 scans
        maxGapMs1Scan <- max(diff(metaSubTmp$retentionTime[
            metaSubTmp$msLevel == 1])) + 1
        # single point rt drift
        if (!is.null(intStdRtDrift)) {
            metaSubTmp$retentionTime <- as.numeric(metaSubTmp$retentionTime)
            + (intStdRtDrift * -1)
        }
        # scan numbers
        indxScans <-
            which(
                metaSubTmp$retentionTime < minMaxRt[2] &
                    metaSubTmp$retentionTime >
                    minMaxRt[1] & metaSubTmp$msLevel == 1
            )
        peakRangeAll <-
            do.call(rbind, mzR::peaks(MSfile, indxScans))
        # add retention times
        peakRangeAll <- cbind(
            peakRangeAll,
            rep(metaSubTmp$retentionTime[indxScans],
                metaSubTmp$peaksCount[indxScans]),
            rep(metaSubTmp$seqNum[indxScans],
                metaSubTmp$peaksCount[indxScans])
        )
        # retention time deviation model
        if (!is.null(rtDevModels)) {
            rtDevModel <- rtDevModels[[i]]
        } else {
            rtDevModel <- NULL
        }
        nCores <- ifelse(is.null(nCores), 1, nCores)
        if (nCores > 1) {
            pmt <- proc.time()
            message(paste0(
                "Starting SNOW cluster with ",
                nCores,
                " local sockets...\n"
            ))
            message("identifying and quantifying ",
                    nrow(targTable),
                    " features...\n")
            cl <- parallel::makeCluster(nCores)
            doSNOW::registerDoSNOW(cl)
            # foreach and dopar from foreach package
            resultsTmp <- foreach(j =
                                    seq_len(nrow(targTable))) %dopar% {
                  mzTmp <- as.numeric(targTable[j, 1])
                  rtTmp <- as.numeric(targTable[j, 2]) * 60
                  isoPat <- 1.0032 / targTable$chargeState[j]
                  isoPatPred <- predIsoDist[[as.character(targTable$peptide[
                    j])]]
                  isoPat <-
                      seq(0, isoPat * (nrow(isoPatPred) - 1),
                        isoPat)
                  names(isoPat) <- isoPatPred$id
                  # retention time deviation model
                  if (!is.null(rtDevModel)) {
                      rtDevTmp <- predict(rtDevModel, newdata = rtTmp / 60)
                      predRtTmp <- rtTmp + (rtDevTmp * 60)
                  } else {
                      predRtTmp <- rtTmp
                  }
                  # different parameters for hk peptide
                  hkIndx <-
                      targTable[j, "peptide"] %in% hkPeptide
                  if (hkIndx) {
                      peakRangeRtSub <- peakRangeAll[which(peakRangeAll[
                        , 3] < {
                        predRtTmp + 240
                    } & peakRangeAll[, 3] > {
                        predRtTmp - 240
                    }), , drop = FALSE]
                  } else {
                    peakRangeRtSub <- peakRangeAll[which(
                        peakRangeAll[, 3] <
                            (predRtTmp +
                                maxRtWindow) & peakRangeAll[, 3] >
                            (predRtTmp - maxRtWindow)
                    ),
                    , drop = FALSE]
                }
                # mass drift correction
                if (file.exists(massDriftFile)) {
                    peakRangeRtSub[, 1] <- peakRangeRtSub[, 1] - {
                        {
                            peakRangeRtSub[, 1] / 1e+06
                        } * massDriftTable[peakRangeRtSub[, 4],
                                            "ppmDrift"]
                    }
                }
                resultTMP <-
                    peakIdQuant_newMethod(
                        mzTmp = mzTmp,
                        rtTmp = rtTmp,
                        peakRangeRtSub = peakRangeRtSub,
                        rtDevModel = rtDevModel,
                        isoPat = isoPat,
                        isoPatPred = isoPatPred,
                        minSimScore = minSimScore,
                        maxPpm = maxPpm,
                        gaussAlpha = gaussAlpha,
                        spikeScans = spikeScans,
                        minPeakHeight = minPeakHeight,
                        maxRtDrift = ifelse(hkIndx, 60, maxRtDrift),
                        isoWindow = isoWindow,
                        maxGapMs1Scan = maxGapMs1Scan,
                        intMaxPeak = hkIndx
                    )
            }
            # stop SNOW cluster
            parallel::stopCluster(cl)
            proc.time() - pmt
            for (res in seq_along(resultsTmp)) {
                resultRow <- ((i - 1) * nrow(targTable)) + res
                if (length(resultsTmp[[res]]) == 3) {
                    results[resultRow, 6:ncol(results)] <-
                        resultsTmp[[res]]$results
                    resultsTmp[[res]]$results <- NULL
                }
                resultsList[[resultRow]] <-
                    resultsTmp[[res]]
            }
            # single threaded
        } else {
            pmt <- proc.time()
            for (j in intSeq) {
                resultRow <- ((i - 1) * nrow(targTable)) + j
                mzTmp <- as.numeric(targTable[j, 1])
                rtTmp <- as.numeric(targTable[j, 2]) * 60
                isoPat <- 1.0032 / targTable$chargeState[j]
                isoPatPred <- predIsoDist[[as.character(targTable$peptide[j])]]
                isoPat <-
                    seq(0, isoPat * (nrow(isoPatPred) - 1),
                        isoPat)
                names(isoPat) <- isoPatPred$id
                # retention time deviation model
                if (!is.null(rtDevModel)) {
                    rtDevTmp <- predict(rtDevModel, newdata = rtTmp / 60)
                    predRtTmp <- rtTmp + (rtDevTmp * 60)
                } else {
                    predRtTmp <- rtTmp
                }
                # different parameters for hk peptide
                hkIndx <-
                    targTable[j, "peptide"] %in% hkPeptide
                if (hkIndx) {
                    peakRangeRtSub <- peakRangeAll[which(peakRangeAll[, 3] < {
                        predRtTmp + 240
                    } & peakRangeAll[, 3] > {
                        predRtTmp - 240
                    }), , drop = FALSE]
                } else {
                    peakRangeRtSub <- peakRangeAll[which(
                        peakRangeAll[, 3] <
                            (predRtTmp + maxRtWindow) &
                            peakRangeAll[, 3] >
                            (predRtTmp - maxRtWindow)
                    ),
                    , drop = FALSE]
                }
                # mass drift correction
                if (file.exists(massDriftFile)) {
                    peakRangeRtSub[, 1] <- peakRangeRtSub[, 1] - {
                        {
                            peakRangeRtSub[, 1] / 1e+06
                        } * massDriftTable[peakRangeRtSub[, 4],
                                            "ppmDrift"]
                    }
                }
                resultsList[[resultRow]] <-
                    peakIdQuant_newMethod(
                        mzTmp
                        = mzTmp,
                        rtTmp = rtTmp,
                        peakRangeRtSub = peakRangeRtSub,
                        rtDevModel = rtDevModel,
                        isoPat = isoPat,
                        isoPatPred = isoPatPred,
                        minSimScore = minSimScore,
                        maxPpm = maxPpm,
                        spikeScans = spikeScans,
                        minPeakHeight = minPeakHeight,
                        maxRtDrift = ifelse(hkIndx, 60, maxRtDrift),
                        isoWindow = isoWindow,
                        maxGapMs1Scan = maxGapMs1Scan,
                        intMaxPeak = hkIndx
                    )
                ########################### PLOTS #############
                colnames(peakRangeRtSub) = c("MIM", "int",
                                            "adjRT", "seq Number")
                print(head(peakRangeRtSub))
            }  # end feature loop
            proc.time() - pmt
        }  # if parallel
    }  # end file loop
    # create new adductQuant object
    if (!is.null(quantObject)) {
        object <- quantObject
    } else {
        object <- new("AdductQuantif")
        predIsoDist(object) <- predIsoDist
        targTable(object) <- targTable
        file.paths(object) <- filePaths
    }
    peakQuantTable(object) <- results
    peakIdData(object) <- resultsList
    save(object, file = "adductQuantResults.Rdata")
    return(object)
    }  # end Function
setMethod("show", "AdductQuantif", function(object) {
    if (length(file.paths(object)) > 0) {
        message(
            "A \"AdductQuantif\" class object derived from ",
            length(file.paths(object)),
            " MS files \n\n"
        )
        message(
            "Consisting of:\n",
            sum(peakQuantTable(object)[,
                                        "peakArea"] != "0"),
            " quantified peaks\n",
            "and ",
            sum(peakQuantTable(object)[, "peakArea"] ==
                    "0"),
            " missing peaks (i.e. zero peak area)\n",
            "Derived from a target list of: ",
            nrow(targTable(object)),
            " targets.\n\n"
        )
        if (any(grepl("^possOutPeak$",
                    colnames(peakQuantTable(object))))) {
            message("Potentially outlying peaks:")
            message(
                sum(peakQuantTable(object)[, "possOutPeak"] == 1),
                paste0("(", round((
                    sum(peakQuantTable(object)[,
                                            "possOutPeak"] == "1") /
                        sum(peakQuantTable(object)[, "peakArea"] !="0") * 100
                ), 1), "%)", collapse = ""),
                "of a total of",
                sum(peakQuantTable(object)[,
                                        "peakArea"] != "0"),
                "peaks quantified were potentially non-gaussian,
                long tailed,
                or deviant from the peak group
                retention time or peak area.\n\n"
            )
        }
        memsize <- object.size(object)
        message("Memory usage:", signif(memsize / 2 ^ 20, 3), "MB")
        } else {
            message("A new empty\"AdductQuantif\" class object")
        }
    })
setMethod("c", signature(x = "AdductQuantif"), function(x, ...) {
    elements = list(...)
    # error handling check if all adductSpec object
    if (any(vapply(elements, function(ele)
        is(ele, 'AdductQuantif'),
        FUN.VALUE = logical(1)) == FALSE)) {
        stop("all elements must be an AdductQuantif class object")
    }
    emptyAdductQuantif <- new("AdductQuantif")
    for (i in seq_along(elements)) {
        if (ncol(peakQuantTable(emptyAdductQuantif))==0) {
            peakQuantTable(emptyAdductQuantif) <-
                peakQuantTable(elements[[i]])
        } else {
            peakQuantTable(emptyAdductQuantif) <-
                rbind(peakQuantTable(emptyAdductQuantif),
                    peakQuantTable(elements[[i]]))
        }
        peakIdData(emptyAdductQuantif) <-
            c(peakIdData(emptyAdductQuantif),
            peakIdData(elements[[i]]))
        predIsoDist(emptyAdductQuantif) <-
            c(predIsoDist(emptyAdductQuantif),
            predIsoDist(elements[[i]]))
        # keep unique
        predIsoDist(emptyAdductQuantif) <-
            predIsoDist(emptyAdductQuantif)[!duplicated(names(predIsoDist(
                emptyAdductQuantif)))]
        # target table
        targTable(emptyAdductQuantif) <-
            rbind(targTable(emptyAdductQuantif),
                targTable(elements[[i]]))
        # remove duplicates
        uniEntries <-
            apply(targTable(emptyAdductQuantif)[, seq_len(3)],
                1, paste, collapse = "")
        targTable(emptyAdductQuantif) <-
            targTable(emptyAdductQuantif)[duplicated(uniEntries) ==
                                            FALSE, , drop = FALSE]
        file.paths(emptyAdductQuantif) <-
            c(file.paths(emptyAdductQuantif),
            file.paths(elements[[i]]))
    }
    return(emptyAdductQuantif)
    })  # end function