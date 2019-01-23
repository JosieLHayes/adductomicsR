#' group MS/MS precursor masses
#'
#' @description hierarchically cluster ms/ms precursor scans within and across
#' samples, according to a m/z and retention time error.
#' @return a list identical to adductSpectra containing an additional list
#' element:
#' @param adductSpectra AdductSpec object
#' @param nCores numeric the number of cores to use for parallel computation. 
#' The default is to use all available cores detected using the function 
#' parallel::detectCores()
#' @param maxRtDrift numeric for the maximum rentention time
#' drift to be considered. Default is 20.
#' @param ms1mzError numeric maximum MS1 mass:charge error
#' @param ms2mzError numeric maximum MS2 mass:charge error
#' @param dotProdClust logical remove previous dot prod clustering results
#' @param minDotProd numeric. Minimum mean dot product spectral similarity 
#' score to keep a spectrum within an MS/MS group (default = 0.8).
#' @param disMetric metric to use for distance in clustering
#' @param fclustMethod method to use for the fclust function
#' @param compSpecGen logical for whether composite spectra generation
#' is necessary
#' @usage ms2Group(adductSpectra = NULL, nCores = NULL, maxRtDrift = NULL, 
#' ms1mzError = 0.1, ms2mzError = 1, dotProdClust = TRUE, minDotProd = 0.8,
#' fclustMethod = "median", disMetric = "euclidean", compSpecGen = TRUE,
#' adjPrecursorMZ = TRUE)
#' @param adjPrecursorMZ logical for precursor mass:charge adjustment
#'
ms2Group <- function(adductSpectra = NULL, nCores = parallel::detectCores(), 
maxRtDrift = NULL, ms1mzError = 0.1, 
ms2mzError = 1, dotProdClust = TRUE, minDotProd = 0.8, 
fclustMethod = "median", 
disMetric = "euclidean", compSpecGen = TRUE, adjPrecursorMZ = TRUE) {
    # error handling
    if (is.null(adductSpectra)) {
        stop("argument adductSpectra is missing with no default")
    } else if (!is(adductSpectra,'AdductSpec')) {
    
        stop("argument adductSpectra is not an AdductSpec class object.")
    }
    # check if the parameters are the same and stop if re-grouping unneccessary
    if (ncol(Parameters(adductSpectra)) > 0) {
        if (all(c("maxRtDrift", "ms1mzError", "dotProdClust") %in% 
            colnames(Parameters(adductSpectra)))) {
                if (!is.null(maxRtDrift)) {
                    if (Parameters(adductSpectra)$maxRtDrift == 
                        maxRtDrift & 
                    Parameters(adductSpectra)$ms1mzError == 
                        ms1mzError &
                         Parameters(adductSpectra)[,'dotProdClust'] ==
                        FALSE)
                        {
                            stop("Grouping and composite spectrum generation is 
                            unneccessary as maxRtDrift and ms1mzError 
                            arguments are un-changed from 
                            the previous round of grouping...\n")
                        }
                    }
                }
            }
            # empty peptide identification results slots
            specPepMatches(adductSpectra) <- vector("list", 0)
            sumAdductType(adductSpectra) <- data.frame()
            Peptides(adductSpectra) <- data.frame()
            metaDataTmp <- metaData(adductSpectra)
            metaDataTmp$orderTmp <- seq_len(nrow(metaDataTmp))
            indxTmp <- metaDataTmp$aboveMinPeaks == 1
            metaDataTmp <- metaDataTmp[indxTmp, ]
            # all spectra 1 list
            allSpectra <- unlist(adductMS2spec(adductSpectra), recursive =
            FALSE)
            # check if mz error has changed if not then do not re-do 
            # hierarchical clustering of mass
            hclustMass <- TRUE
            if (!is.null(Parameters(adductSpectra)$ms1mzError)) {
                hclustMass <- ifelse(Parameters(adductSpectra)$ms1mzError == 
                ms1mzError, FALSE, 
            TRUE)
            }
            if (hclustMass == TRUE) {
                message("Hierarchically clustering ", 
                prettyNum(nrow(metaDataTmp),
                big.mark = ","), 
                " MS/MS precursors and cutting dendrogram according 
                to a mass error of ", 
            ms1mzError, "...\n")
                flush.console()
                # if necessary adjust the mass drift
                if ({
                    "adjPrecursorMZ" %in% colnames(metaDataTmp)
                } == FALSE) {
                    adjPrecursorMZ <- FALSE
                }
                if (adjPrecursorMZ == TRUE) {
                    # add ppm mass drift value
                    precursorMZs <- metaDataTmp$precursorMZ
                } else {
                    precursorMZs <- metaDataTmp$precursorMZ
                }

                hr <- fastcluster::hclust.vector(precursorMZs, 
                metric = disMetric, 
                method = fclustMethod)
                # m/z error
                metaDataTmp$msPrecursor_group <- 0
                metaDataTmp$msPrecursor_group <- cutree(hr, h = ms1mzError)
            }
            # order by msPrecursor group
            orderMsPreC <- order(metaDataTmp$msPrecursor_group)
            metaDataTmp <- metaDataTmp[orderMsPreC, ]
            allSpectra <- allSpectra[orderMsPreC]
            # adjust RT if necessary adjust with intStd
            if (!is.null(metaDataTmp$intStdRtDrift)) {
                retentionTime <- as.numeric(metaDataTmp$retentionTime) + 
                (as.numeric(metaDataTmp$intStdRtDrift) * 
            -1)
            } else {
                retentionTime <- as.numeric(metaDataTmp$retentionTime)
            }
            if (length(rtDevModels(adductSpectra)) > 0) {
                message("Adjusting retention time based on loess model...")
                flush.console()
                metaDataTmp$predRtLoess <- 0
                for (x in seq_len(length(Specfile.paths(adductSpectra)))) {
                    indxFileTmp <- which(metaDataTmp$mzXMLFile %in% 
                    basename(Specfile.paths(adductSpectra))[x])
                    metaDataTmp$predRtLoess[indxFileTmp] <- 
                    retentionTime[indxFileTmp]-
                    (predict(rtDevModels(adductSpectra)[[x]], newdata = 
                    as.numeric(retentionTime[indxFileTmp])/60) * 
                60)
                }
                retentionTime <- as.numeric(metaDataTmp$predRtLoess)
            }
    
            # if unneccessary do not re do retention time clustering
            hclustRt <- TRUE
            if (!is.null(Parameters(adductSpectra)$maxRtDrift)) {
                if (!is.null(maxRtDrift)) {
                    hclustRt <- ifelse(Parameters(adductSpectra
                        )[,'maxRtDrift'] == maxRtDrift, FALSE, TRUE)
                } else {
                    hclustRt <- FALSE
                }
            }
            if (hclustRt == TRUE) {
                # matching ms/ms scans across samples
                message("clustering ms/ms scans across samples based on a
                maximum RT drift of ", maxRtDrift, ".\n")
                flush.console()
                interMSMSrtGroups <- tapply(retentionTime, 
                metaDataTmp$msPrecursor_group, 
                function(x) {
                    if (length(x) > 1) {
                        cutree(fastcluster::hclust.vector(x, metric = 
                        disMetric,
                        method = fclustMethod), 
                    h = maxRtDrift)
                    } else {
                        1
                    }
                })
                    # interMSMS groups
                    metaDataTmp$interMSMSrtGroups <- 
                    paste0(metaDataTmp$msPrecursor_group,
                    "_", unlist(interMSMSrtGroups))
                    message(length(unique(metaDataTmp$interMSMSrtGroups)), 
                    " MS/MS peak groups identified at maxRtDrift of ", 
                maxRtDrift, " seconds")
                    flush.console()
                }
                if (dotProdClust == TRUE) 
                    {
                        # spectral similarity all ms/ms groups # order by 
                        #adduct MS2 file
                        #names
                        spectrumIds <- paste0(metaDataTmp$mzXMLFile,
                        ".MS2spectra",
                        ".", metaDataTmp$seqNum)
                        allSpectra <- allSpectra[spectrumIds]
            
                        # remove previous dot prod clustering results
                        metaDataTmp$interMSMSrtGroups <- gsub("\\..+", "",
                    metaDataTmp$interMSMSrtGroups)
                        # intra-MSMS group spectral similarity (dot product) 
                        # clustering
                        message("intra-MS/MS group (n=", 
                        length(unique(metaDataTmp$interMSMSrtGroups)), 
                        " RT + m/z groups) spectral similarity
                        (dot product) clustering...\n")
                        flush.console()
                        specByGroup <- split(allSpectra, 
                        metaDataTmp$interMSMSrtGroups)
                        spectraList <- NULL
                        groupName <- NULL
                        dotProdSimClust <- function(spectraList = NULL, 
                        groupName = NULL, 
                        minDotProdThresh = minDotProd, 
                        binSizeMS2 = 1) {
                            # if length 1 return nothing
                            clusterId <- groupName
                            if (length(spectraList) > 1) {
                                specNamesVecTmp <- rep(names(spectraList), 
                                vapply(spectraList, 
                                nrow, FUN.VALUE=numeric(1)))
                                allSpecTmp <- do.call(rbind, spectraList)
                                maxMass <- floor(max(allSpecTmp[, 1])) + 10
                                # padded integer labels
                                message("Calculating dot product matrix ", 
                                length(unique(specNamesVecTmp)), 
                                " spectra\n")
                                flush.console()
                                # padded integer labels
                                labelsTmp <- paste0("(", seq(binSizeMS2, 
                                (maxMass - binSizeMS2), 
                                binSizeMS2), ",", seq((2 * binSizeMS2), 
                                maxMass, binSizeMS2), 
                                "]")
                                massBinsIndivTmp <- cut(allSpecTmp[, 1], 
                                breaks = seq(binSizeMS2, 
                                maxMass, binSizeMS2), labels = labelsTmp)
                                # empty bins
                                indivSpecVec <- tapply(allSpecTmp[, 2], 
                                paste0(specNamesVecTmp, 
                                massBinsIndivTmp), sum)
                                # identify any absent bins
                                allBinNames <- 
                                paste0(rep(unique(specNamesVecTmp), 
                                each = length(labelsTmp)), 
                                rep(labelsTmp, length(unique(
                                specNamesVecTmp))))
                                # add absent bins as zeros
                                allBinsTmp <- rep(0, length(allBinNames))
                                names(allBinsTmp) <- allBinNames
                                # ensure indivSpecVec is in right order
                                allBinsTmp[match(names(indivSpecVec), 
                                allBinNames)] <-
                                indivSpecVec
                                indivSpecMat <- matrix(allBinsTmp, byrow = 
                                FALSE, 
                                nrow = length(labelsTmp))
                                # mean all pairwise dotproducts
                                dotProdMat <- crossprod(indivSpecMat)
                                sqrtMatrixTmp <- matrix(sqrt(colSums(
                                indivSpecMat^2)), 
                                nrow = nrow(dotProdMat), 
                                ncol = ncol(dotProdMat), byrow = TRUE)
                                dotProdsTmp <- dotProdMat/(sqrtMatrixTmp * 
                                diag(sqrtMatrixTmp))
                                hr <- fastcluster::hclust(as.dist(1 - 
                                dotProdsTmp), 
                                method = "single")
                                clusterId <- cutree(hr, h = {
                                    1 - minDotProdThresh
                                })
                                    clusterId <- paste0(groupName, ".", 
                                clusterId)
                                }
                                return(clusterId)
                            }  # end dot product hierarchical clustering
                            if (!is.null(nCores)) {
                                # multi-threaded
                                pmt <- proc.time()
                                message(paste0("Starting SNOW cluster with ",
                                nCores,
                                " local sockets...\n"))
                                flush.console()
                                cl <- parallel::makeCluster(nCores)
                                doSNOW::registerDoSNOW(cl)
                                dotProdResults <- foreach(spectraList = 
                                specByGroup, 
                                groupName = names(specByGroup), 
                                .combine="c",
                                .packages=c("fastcluster"))%dopar%{
                                dotProdSimClust(spectraList, groupName, 
                                binSizeMS2 = 
                                ms2mzError)
                                }
                                # stop SNOW cluster
                                parallel::stopCluster(cl)
                                proc.time() - pmt
                                # sort dot prod results to reinsert in correct
                                # order
                                indxSpecTmp <- match(do.call(c, lapply(
                                specByGroup, names)), 
                                spectrumIds)
                                dotProdResults[indxSpecTmp] <- dotProdResults
                            } else {
                                dotProdResults <- vector("character", 
                                nrow(metaDataTmp))
                                for (i in seq_len(length(specByGroup))) {
                                    message(i, " of ", length(specByGroup), 
                                    " m/z + RT groups.\n")
                                    flush.console()
                                    dotProdResTmp <- dotProdSimClust(
                                    spectraList = 
                                    specByGroup[[i]], 
                                    groupName = names(specByGroup)[i], 
                                    binSizeMS2 = ms2mzError)
                                    indxSpecTmp <- match(names(
                                    specByGroup[[i]]), spectrumIds)
                                    dotProdResults[indxSpecTmp] <- 
                                    dotProdResTmp 
                                } 
                            }
                            # intra-MSMS group spectral similarity (dot
                            # product) clustering
                            metaDataTmp$interMSMSrtGroups <- dotProdResults
                            message(prettyNum(length(unique(
                            metaDataTmp$interMSMSrtGroups))-1, 
                            big.mark = ","), " inter-MS/MS groups 
                            following spectral
                            similarity (dot product) clustering...\n")
                            flush.console()
                        }  
                        if (compSpecGen == TRUE) 
                            {
                                # order by adduct MS2 file names
                                spectrumIds <- paste0(metaDataTmp$mzXMLFile,
                                ".MS2spectra",
                                ".", metaDataTmp$seqNum)
                                allSpectra <- allSpectra[spectrumIds]
                                compSpecTable <- table(
                                metaDataTmp$interMSMSrtGroups)
                                message("Generating ", sum(compSpecTable > 1), 
                                " composite spectra...\n")
                                flush.console()
                                # if nCores not null then parallel comp
                                if (!is.null(nCores)) {
                                    # all single spectra groups first
                                    interMSMScompSpec <- vector("list", 
                                    length(unique(
                                    metaDataTmp$interMSMSrtGroups)))
                                    names(interMSMScompSpec)<-unique(
                                    metaDataTmp$interMSMSrtGroups)
                                    singleGroupIndx <- match(names(
                                    compSpecTable)[
                                    compSpecTable == 1], 
                                    metaDataTmp$interMSMSrtGroups)
                                    listIndx<-match(
                                    metaDataTmp$interMSMSrtGroups[
                                    singleGroupIndx], 
                                    names(interMSMScompSpec))
                                    interMSMScompSpec[listIndx]<-allSpectra[
                                    paste0(metaDataTmp$mzXMLFile[
                                    singleGroupIndx], 
                                    ".MS2spectra", ".", metaDataTmp$seqNum[
                                    singleGroupIndx])]
                                    multSampNames <- names(compSpecTable)[
                                    compSpecTable > 1]
                                    # mult sample index
                                    multSampIndx <- match(names(
                                    interMSMScompSpec), multSampNames)
                                    message(paste0("Starting SNOW cluster 
                                    with ", nCores, " local sockets...\n"))
                                    flush.console()
                                    cl <- parallel::makeCluster(nCores)
                                    doSNOW::registerDoSNOW(cl)
                                    compSpecTmp <- foreach(j = multSampNames,
                                    .packages = c("fastcluster", 
                                    "adductomics")) %dopar% {
                                        subSpecNamesTmp <- spectrumIds[
                                        metaDataTmp$interMSMSrtGroups
                                        %in% j]
                                        spec.df <- do.call(rbind, allSpectra[
                                        subSpecNamesTmp])
                                        return(signalGrouping(spec.df, 
                                        ms2mzError, minPeaks = 0))
                                    }
                                    # stop SNOW cluster
                                    parallel::stopCluster(cl)
                                    interMSMScompSpec[which(!is.na(
                                    multSampIndx))]<-compSpecTmp[
                                    multSampIndx[!is.na(multSampIndx)]]
                                } else {
                                    pmt <- proc.time()
                                    interMSMScompSpec <- vector("list",
                                    length(unique(
                                    metaDataTmp$interMSMSrtGroups)))
                                    names(interMSMScompSpec)<-unique(
                                    metaDataTmp$interMSMSrtGroups)
                                    singleGroupIndx <- match(names(
                                    compSpecTable)[
                                    compSpecTable == 1],
                                    metaDataTmp$interMSMSrtGroups)
                                    listIndx <- match(
                                    metaDataTmp$interMSMSrtGroups[
                                    singleGroupIndx],
                                    names(interMSMScompSpec))
                                    interMSMScompSpec[listIndx]<-allSpectra[
                                    paste0(metaDataTmp$mzXMLFile[
                                    singleGroupIndx], 
                                    ".MS2spectra", ".", metaDataTmp$seqNum[
                                    singleGroupIndx])]
                                    compSpecIndx <- names(compSpecTable)[
                                    compSpecTable > 1]
                                    pb<-txtProgressBar(min = 0, 
                                    max = length(compSpecIndx),style = 3)
                                    for (j in seq_len(length(compSpecIndx))) {
                                        setTxtProgressBar(pb, j)
                                        intGroup <- compSpecIndx[j]
                                        listIndxTmp <- which(names(
                                        interMSMScompSpec) == intGroup)
                                        subSpecNamesTmp <- spectrumIds[
                                        metaDataTmp$interMSMSrtGroups
                                        %in% intGroup]
                                        spec.df <- do.call(rbind, allSpectra[
                                        subSpecNamesTmp])
                                        interMSMScompSpec[[listIndxTmp]] <- 
                                        signalGrouping(spec.df,
                                        ms2mzError, 
                                    minPeaks = 0)
                                    }
                                    proc.time() - pmt
                                }  # nCores
                                # add inter sample spectra to adductspectra
                                groupMS2spec(adductSpectra
                                    ) <- lapply(
                                interMSMScompSpec,
                                function(subL) {
                                    colnames(subL) <- c("mass", "intensity")
                                    return(data.frame(subL, stringsAsFactors = 
                                FALSE))
                                })
                                }  
                                metaDataTmp <- metaDataTmp[order(
                                metaDataTmp$orderTmp), ]
                                metaDataTmp$orderTmp <- NULL
                                # add mod metaData back to AdductSpec object 
                                # add additional columns if necc
                                missCol <- setdiff(colnames(metaDataTmp), 
                                colnames(metaData(
                                    adductSpectra)))
                                metaData(adductSpectra
                                    )[, missCol] <- ""
                                nameIndx <- match(colnames(metaDataTmp), 
                            colnames(metaData(adductSpectra)))
                            metaData(adductSpectra) <- 
                                metaData(adductSpectra)[, nameIndx]
                                metaData(adductSpectra)[indxTmp, ] <- 
                                metaDataTmp
                                # add parameters to parameters slot
                                argsTmp <- c("maxRtDrift", "ms1mzError",
                                "ms2mzError", "dotProdClust", 
                                "minDotProd", 
                                "fclustMethod", "disMetric")
                                paramTmp <- data.frame(ifelse(!is.null(
                                maxRtDrift), maxRtDrift, "NULL"),
                                ms1mzError,  ms2mzError, dotProdClust, 
                                minDotProd, fclustMethod, 
                                disMetric, stringsAsFactors = FALSE)
                                colnames(paramTmp) <- argsTmp
                                if (ncol(Parameters(adductSpectra)) == 0) {
                                    Parameters(adductSpectra) <- paramTmp
                                } else {
                                    Parameters(adductSpectra)[, argsTmp] <- 
                                    paramTmp
                                }
                                return(adductSpectra)
                            }  # end function