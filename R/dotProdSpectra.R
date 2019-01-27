#' dot product calculation
#' @description hierarchical clustering 
#' (complete method see \code{\link[stats:hclust]{hclust}}).
#' Dissimilarity metric based on 1-dot product spectral similarity. 
#' Retention time and
#' mass groups are therefore further subdivided based on spectral similarity. 
#' If outlying mass spectra have been erroneously grouped then these will 
#' be reclassified.
#' The dopProdSimClust function calculates the dot product and returns a 
#'  cluster id. This function is used later in the script on all spectra
#' in order to group the spectra (removing peak tails etc)
#' @param adductSpectra adductSpec object
#' @param minDotProdSpec numeric minimum dot product score
#' @param nCores numeric the number of cores to use for parallel computation. 
#' The default is to use all available cores detected using 
#' the function parallel::detectCores()
#' @param maxGroups numeric maximum number of groups to include from 
#' the dendrogram.
#' @usage dotProdSpectra(adductSpectra = NULL, nCores = NULL, 
#' minDotProdSpec = 0.8, maxGroups = 10)
#' @return adductSpectra adductSpec object
dotProdSpectra <- function(adductSpectra = NULL, nCores = NULL, 
minDotProdSpec = 0.8, maxGroups = 10) {
    # error handling
    if (is.null(adductSpectra)) {
        stop("argument adductSpectra is missing with no default")
    } else if (!is(adductSpectra,'adductSpec')) {
        stop("argument adductSpectra is not an adductSpec class object.")
    }
    if (is.null(metaData(adductSpectra)[,'interMSMSrtGroups'])) {
        stop("spectral grouping using the ?ms2Group function has 
    not yet been performed...")
    }
    # all constituent spectra
    indivSpecs <- unlist(adductMS2spec(adductSpectra), recursive = FALSE)

    uniqueMSMSgroups <- unique(metaData(adductSpectra)[,'interMSMSrtGroups'])
    uniqueMSMSgroups <- uniqueMSMSgroups[uniqueMSMSgroups != ""]
    freqMSMSgroup <- table(metaData(adductSpectra)[,'interMSMSrtGroups'])
    freqMSMSgroup <- freqMSMSgroup[freqMSMSgroup != ""]
    freqMSMSgroup <- names(freqMSMSgroup)[freqMSMSgroup > 1]

    uniqueMSMSgroups <- uniqueMSMSgroups[uniqueMSMSgroups %in% freqMSMSgroup]

    metaDataTmp <- metaData(adductSpectra)[
    metaData(adductSpectra)[,'interMSMSrtGroups']
    %in% uniqueMSMSgroups, , drop = FALSE]
    spectrumIds <- paste0(metaDataTmp$mzXMLFile, ".MS2spectra", ".", 
    metaDataTmp$seqNum)
    indivSpecs <- indivSpecs[spectrumIds]
    indivSpecs <- split(indivSpecs, metaDataTmp$interMSMSrtGroups)
    dotProdSimClust <- function(spectraList = NULL, compSpecName = NULL, 
    maxGroupsDend = maxGroups, 
    minDotProdThresh = minDotProdSpec) {
        specNamesVecTmp <- rep(names(spectraList), vapply(spectraList, nrow), 
        FUN.VALUE=numeric(1))
        allSpecTmp <- do.call(rbind, spectraList)
        labelsTmp <- paste0(sprintf("(%04d", seq_len(1999)), ",", 
        sprintf("%04d", 2:2000), 
        "]")
        massBinsIndivTmp <- cut(allSpecTmp[, 1], breaks = seq(1, 2000, 1), 
        labels = labelsTmp)
        # empty bins
        emptyBins <- cut(seq(2, 2000, 1), breaks = seq(1, 2000, 1), 
        labels = labelsTmp)
        indivSpecVec <- tapply(allSpecTmp[, 2], paste0(specNamesVecTmp, 
        massBinsIndivTmp), 
        sum)
        # identify any absent bins
        allBinNames <- paste0(rep(unique(specNamesVecTmp), 
        each = length(emptyBins)), 
        rep(emptyBins, length(unique(specNamesVecTmp))))
        # add absent bins as zeros
        allBinsTmp <- rep(0, length(allBinNames))
        names(allBinsTmp) <- allBinNames
        allBinsTmp[which(allBinNames %in% names(indivSpecVec))] <- indivSpecVec
        indivSpecMat <- matrix(allBinsTmp, byrow = FALSE, nrow =
        length(emptyBins))
        # mean all pairwise dotproducts dotProdMat 
        dotProdMat <- crossprod(indivSpecMat)
        sqrtMatrixTmp <- matrix(sqrt(colSums(indivSpecMat^2)), 
        nrow = nrow(dotProdMat), 
        ncol = ncol(dotProdMat), byrow = TRUE)
        subDotProds <- dotProdMat/(sqrtMatrixTmp * diag(sqrtMatrixTmp))
        hr <- fastcluster::hclust(as.dist(1 - subDotProds), 
        method = "complete")
        clusterId <- cutree(hr, h = 1 - minDotProdThresh)
        if (length(unique(clusterId)) > maxGroupsDend) {
            clusterId <- cutree(hr, k = maxGroupsDend)
        }
        # add the within group means of dot products to cluster ids
        diag(subDotProds) <- NA
        namesClustId <- vapply(seq_along(clusterId), function(x) {
            subGroupIndxTmp <- clusterId == clusterId[x]
            if (sum(subGroupIndxTmp) > 1) {
                return(mean(subDotProds[x, ], na.rm = TRUE))
            } else {
                return(1)
            }
        }, FUN.VALUE=numeric(1))
        clusterId <- paste0(compSpecName, ".", clusterId)
        names(clusterId) <- namesClustId
        # }
    return(clusterId)
    }  # end function
    nCores <- ifelse(is.null(nCores), 1, nCores)
    if (nCores > 1) {
        # multi-threaded
        pmt <- proc.time()
        message(paste0("Starting SNOW cluster with ", nCores, 
        " local sockets...\n"))
        
        cl <- parallel::makeCluster(nCores)
        doSNOW::registerDoSNOW(cl)
        dotProdResults <- foreach(i = seq_along(indivSpecs),
        .combine = "c",
        .packages = c("fastcluster")) %dopar% 
        {
            dotProdSimClust(spectraList = indivSpecs[[i]], 
            compSpecName = names(indivSpecs)[i])
        }
        # stop SNOW cluster
        parallel::stopCluster(cl)
        proc.time() - pmt
        # sort dot prod results to reinsert in correct order
        indxSpecTmp <- match(do.call(c, lapply(indivSpecs, names)), 
        spectrumIds)
        dotProdResults[indxSpecTmp] <- dotProdResults
    } else {
        dotProdResults <- vector("character", length(spectrumIds))
        pmt <- proc.time()

        pb <- txtProgressBar(max = length(indivSpecs), style = 3)
        for (i in seq_along(indivSpecs)) {
            setTxtProgressBar(pb, i)
            dotProdResTmp <- dotProdSimClust(spectraList = indivSpecs[[i]], 
            compSpecName = names(indivSpecs)[i])
            indxSpecTmp <- match(names(indivSpecs[[i]]), spectrumIds)
            dotProdResults[indxSpecTmp] <- dotProdResTmp
            names(dotProdResults)[indxSpecTmp] <- names(dotProdResTmp)        
        }  
    }
    proc.time() - pmt
    prevMetaData <- metaData(adductSpectra)
    prevSpectrumIds <- paste0(prevMetaData$mzXMLFile, ".MS2spectra", ".", 
    prevMetaData$seqNum)
    tmpIndx <- match(prevSpectrumIds, spectrumIds)
    prevMetaData$interMSMSrtGroups[which(!is.na(tmpIndx))]<-
    dotProdResults[tmpIndx[!is.na(tmpIndx)]]
    # add mean dot prod rt group
    prevMetaData$meanDotProd_mzRtGroup <- 0
    prevMetaData$meanDotProd_mzRtGroup[which(!is.na(tmpIndx))]<-
    as.numeric(names(dotProdResults)[tmpIndx[!is.na(tmpIndx)]])
    metaData(adductSpectra) <- prevMetaData
    return(adductSpectra)
}  # end function