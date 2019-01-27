#'
#' spectral similarity based adducted peptide
#' identification for adductomicsR
#'
#' @param MS2Dir character a full path to a directory containing
#' either .mzXML or .mzML data
#' @param nCores numeric the number of cores to use for parallel computation.
#' The default
#' is to use all available cores detected using the function
#' parallel::detectCores()
#' @param rtDevModels a list object or a full path to an RData file
#' containing the retention time deviation models for the dataset.
#' @param topIons numeric the number of most intense ions to consider for the
#' basepeak to fragment mass difference calculation (default = 100).
#' Larger values will slightly increase computation time, however when the
#' modified/variable ions happen to be low abundance this value should be set
#' high to ensure these fragment ions are considered.
#' @param topIntIt numeric the number of most intense peaks to calculate the
#' peak to peak mass differences from (default = 5 i.e.
#' the base peak and the next
#' 4 most intense ions greater than 10 daltons in mass from one
#' another will be considered
#' the multiple iterations increase computation time but in the case
#' that the peptide spectrum
#' is contaminated/chimeric or the variable ions are of lower intensity
#' this parameter should be
#' increased).
#' @param minDotProd numeric minimum dot product similarity score
#' (cosine) between the
#' model spectra's variable ions and the corresponding intensities
#' of the basepeak to fragment ion
#' mass differences identified in the experimental spectrum scans
#' (default = 0.8). Low values
#' will greatly increase the potential for false positive peptide annotations.
#' @param precCh integer charge state of precursors (default = 3).
#' @param minSNR numeric the minimum signal to noise
#' ratio for a fragment ion to be considered. The noise level for
#' each fixed or variable
#' ion is calculated by taking the median of the bottom half of ion intensities
#' within the locality of the fragment ion. The locality is defined
#' as within +/-
#' 100 Daltons of the fragment ion.
#' @param minRt numeric the minimum retention time (in minutes)
#' within which to
#' identify peptide spectra (default=20).
#' @param maxRt numeric the maximum retention time (in minutes)
#' within which to
#' identify peptide spectra (default=45).
#' @param minFixed numeric the minimum number of fixed fragment ions
#' that must have
#' been identified in a spectrum for it to be considered.
#' @param minMz numeric the minimum mass-to-charge ratio of a precursor ion.
#' @param maxMz numeric the maximum mass-to-charge ration of a precursor ion.
#' @param minIdScore numeric the minimum identification
#' score this is an average
#' score of all of the 7 scoring metrics (default=0.4).
#' @param modelSpec character full path to a model spectrum file (.csv).
#' Alternatively
#' built in model tables (in the extdata directory) can be used by just
#' supplying the one letter amino acid code
#' for the peptide (currently available are: "ALVLIAFAQYLQQCPFEDHVK" and
#' "RHPYFYAPELLFFAK"). If supplying a custom table it must
#' consist of the following mandatory columns ("mass", "intensity",
#' "ionType" and "fixed or variable").
#' \enumerate{
#' \item mass - m/z of fragment ions.
#' \item intensity - intensity of fragment ions can be either
#' relative or absolute intensity
#' \item ionType - the identity of the B and Y fragments can optionally
#' added here (e.g. [b6]2+, [y2]1+)
#' or if not known such as for mixed disulfates
#' this column can also contain empty fields.
#' \item fixed or variable - this column contains whether a fragment ion
#' should be considered either 'fixed', 'variable' (i.e. modified) or if it is
#' an empty field it will not be considered.
#' }
#' As default the following model spectra are included in the external data
#' directory of the adductomics package:
#' \enumerate{
#' \item 'modelSpectrum_ALVLIAFAQYLQQCPFEDHVK.csv'
#' \item 'modelSpectrum_RHPYFYAPELLFFAK.csv'
#' }
#' @param groupMzabs numeric after hierarchical clustering of the spectra the
#' dendrogram will be cut at this height (in Da) generating the mass groups.
#' @param groupRtDev numeric after hierarchical clustering of the spectra the
#' dendrogram will be cut at this height (in minutes) generating the
#' retention time groups.
#' @param possFormMzabs numeric the maximum absolute mass difference for
#' matching adduct mass to possible
#' formulae.
#' @param minMeanSpecSim numeric minimum mean dot product
#'similarity score (cosine) between
#' the spectra of a group identified by hierarchical clustering.
#' This parameter
#' is set to prevent erroneous clustering of dissimilar spectra
#' (default = 0.7).
#' @param idPossForm integer if = 1 then the average adduct masses
#' of each spectrum
#' group will be matched against an internal database of possible
#' formula to generate
#' hypotheses. The default 0 mean this will not take place as the
#' computation is
#' potentially time consuming.
#' @param outputPlotDir character (default = NULL) where to save 
#'  plots. Default option of NULL does not save plots.
#' @examples
#' \dontrun{
#' specSimPepId(MS2Dir=system.file("extdata",package
#' ="adductData"),nCores=4,rtDevModels=paste0(system.file("extdata",
#' package ="adductData"),'/rtDevModels.RData'))
#' }
#' @usage specSimPepId(MS2Dir=NULL,nCores=parallel::detectCores(),
#' rtDevModels=NULL, topIons=100, topIntIt=5,minDotProd=0.8, precCh=3,
#' minSNR=3,minRt=20, maxRt=35, minIdScore=0.4,minFixed=3, minMz=750,
#' maxMz=1000,modelSpec=c('ALVLIAFAQYLQQCPFEDHVK','RHPYFYAPELLFFAK'),
#' groupMzabs=0.005, groupRtDev=0.5, possFormMzabs=0.01,
#' minMeanSpecSim=0.7,idPossForm=0, outputPlotDir = NULL)
#' @return dataframe of putative adducts
#' @export
specSimPepId <-
    function(MS2Dir = NULL,
             nCores = NULL,
             rtDevModels = NULL,
             topIons = 100,
             topIntIt = 5,
             minDotProd = 0.8,
             precCh = 3,
             minSNR = 3,
             minRt = 20,
             maxRt = 35,
             minIdScore = 0.4,
             minFixed = 3,
             minMz = 750,
             maxMz = 1000,
             modelSpec = c('ALVLIAFAQYLQQCPFEDHVK',
                           'RHPYFYAPELLFFAK'),
             groupMzabs = 0.005,
             groupRtDev = 0.5,
             possFormMzabs = 0.01,
             minMeanSpecSim = 0.7,
             idPossForm = 0,
             outputPlotDir = NULL) {
        binSizeMS2 = 3
        noiseBin = 200
        minVarPeaks = 0
        if (length(modelSpec) == 2) {
            modelSpec <- modelSpec[1]
            builtInModSpec <- system.file(
                "extdata",
                c(
                    'modelSpectrum_ALVLIAFAQYLQQCPFEDHVK.csv',
                    'modelSpectrum_RHPYFYAPELLFFAK.csv'
                ),
                package = "adductomicsR"
            )
            buildInIndx <-
                grepl(paste0('modelSpectrum_', modelSpec, '\\.csv'),
                      basename(builtInModSpec))
            if (any(buildInIndx)) {
                modelSpec <- builtInModSpec[buildInIndx]
            }
        }
        nCores <- ifelse(is.null(nCores), 1, nCores)
        
        if (is.null(MS2Dir)) {
            stop("Please provide an .mzXML data directory")
        }
        
        if (is.null(rtDevModels)) {
            stop("Please provide an rtDevModels object")
        }
        
        if (is.character(rtDevModels)) {
            rtDevModelsDir <- dirname(rtDevModels)
            message('loading rtDevModels .RData file...Please wait.\n')
            objectName <- load(rtDevModels, envir = environment())
            rtDevModels <- eval(parse(text = objectName))
        }
        pepTmp <- gsub('.+_|\\.csv$', '', basename(modelSpec))
        # calculate the peptide mass
        pepForm <- OrgMassSpecR::ConvertPeptide(pepTmp, IAA = FALSE)
        pepMass <- OrgMassSpecR::MonoisotopicMass(pepForm)
        message(
            'Monoisotopic mass of peptide (',
            pepTmp,
            ') = ',
            prettyNum(pepMass, big.mark = ','),
            '\n'
        )
        # save a parameters file
        todaysDate <- gsub('-', '', Sys.Date())
        parameters <-
            cbind(
                binSizeMS2,
                topIons,
                minDotProd,
                minRt,
                maxRt,
                minVarPeaks,
                minFixed,
                peptide = pepTmp,
                modelSpectrum = basename(modelSpec),
                groupMzabs = groupMzabs,
                groupRtDev = groupRtDev,
                possFormMzabs = possFormMzabs,
                minMeanSpecSim = minMeanSpecSim
            )
        # read model spectrum
        modelSpec <- as.data.frame(read.csv(
            modelSpec,
            header = TRUE,
            stringsAsFactors = FALSE
        ))
        # png
        if (!is.null(outputPlotDir)) {
        png(
            paste0(outputPlotDir, '/', 'modelSpec.png'),
            width = 2200,
            height = 2000,
            res = 275
        )
        par(bg = 'black', fg = 'white')
        plot(
            modelSpec[, seq_len(2)],
            type = 'h',
            main = paste0('Model Spectrum ', pepTmp),
            lwd = 2,
            col = 'white',
            col.axis = 'white',
            col.lab = 'white',
            col.main = 'white',
            col.sub = 'white'
        )
        points(modelSpec[modelSpec$fixed.or.variable == 'fixed', seq_len(2)],
               type = 'h',
               col = 'sandybrown',
               lwd = 2)
        points(modelSpec[modelSpec$fixed.or.variable == 'variable',seq_len(2)],
               type = 'h',
               col = 'yellowgreen',
               lwd = 2)
        points(modelSpec[modelSpec$fixed.or.variable == '', seq_len(2)],
               type = 'h',
               col = 'purple',
               lwd = 2)
        text(
            modelSpec[modelSpec$fixed.or.variable != '', seq_len(2)],
            labels = modelSpec$ionType[modelSpec$fixed.or.variable != ''],
            pos = 3,
            col = 'white'
        )
        legend(
            'topright',
            c('fixed', 'variable', 'not considered'),
            lwd = c(2, 2, 2),
            col = c('sandybrown', 'yellowgreen', 'purple')
        )
        dev.off()
        }
        fixedIons <-
            modelSpec$mass[modelSpec$fixed.or.variable == 'fixed']
        fixedIonNames <-
            modelSpec$ionType[modelSpec$fixed.or.variable == 'fixed']
        modelSpec <-
            modelSpec[modelSpec$fixed.or.variable == 'variable',]
        modelSpec <- modelSpec[order(modelSpec[, 1]),]
        
        # diff between base peak and other ions
        diffSpec <-
            cbind(modelSpec[, 1] - modelSpec[which.max(modelSpec[, 2]), 1],
                  modelSpec[, 2])
        maxMass <-
            floor(max(diffSpec[, 1], na.rm = TRUE)) + {
                binSizeMS2 * 2
            }
        minMass <-
            floor(min(diffSpec[, 1], na.rm = TRUE)) - {
                binSizeMS2 * 2
            }
        # bin them according to expected mass accuracy
        labelsTmp <-
            paste0('(',
                   seq(minMass, (maxMass - binSizeMS2), binSizeMS2),
                   ',',
                   seq((minMass + binSizeMS2), maxMass, binSizeMS2),
                   ']')
        massBinsDiffSpec <-
            cut(diffSpec[, 1],
                breaks = seq(minMass, maxMass,
                             binSizeMS2),
                labels = labelsTmp)
        # empty bins
        massBinsDiffSpec <- tapply(diffSpec[, 2], massBinsDiffSpec, max)
        massBinsDiffSpec[is.na(massBinsDiffSpec)] <- 0
        modelSpecIndx <- massBinsDiffSpec != 0
        # massBinsDiffSpec <- massBinsDiffSpec[modelSpecIndx]
        ms2Files <-
            list.files(MS2Dir, pattern = '\\.mzXML$', full.names = TRUE)
        # reorder ms2Files based on rtDevModel names
        tmpIdx <- match(names(rtDevModels), basename(ms2Files))
        ms2Files <- ms2Files[tmpIdx]
        # save parameters
        parameters <- cbind(parameters, nMS2Files = length(ms2Files))
        write.csv(
            parameters,
            paste0(
                dirname(ms2Files[1]),
                '/specSimPepId_parameters_',
                todaysDate,
                '.csv'
            ),
            row.names = FALSE
        )
        massDriftFiles <-
            gsub('\\.mzXML$|\\.mzML$', '.massDrift.csv', ms2Files)
    
        if (nCores < 2) {
            pb <- txtProgressBar(max = length(ms2Files), style = 3)
            allResults <- data.frame(stringsAsFactors = FALSE)
            for(i in seq_len(length(ms2Files))) {
                setTxtProgressBar(pb, i)
                ms2File <- mzR::openMSfile(ms2Files[i])
                metaData <- mzR::header(ms2File)
                if (file.exists(massDriftFiles[i])) {
                    massDriftTmp <- read.csv(massDriftFiles[i],
                                             header = TRUE,
                                             stringsAsFactors = FALSE)
                    metaData <- cbind(metaData,
                                    massDriftTmp[, c(
                                        'ppmDrift', 'adjPrecursorMZ')])
                } else {
                    metaData$ppmDrift <- NA
                    metaData$adjPrecursorMZ <- metaData$precursorMZ
                }
                
                ms2Indx <- which(metaData$msLevel == 2 &
                                     metaData$precursorCharge %in% precCh)
                # rt range
                metaData$retentionTime <- metaData$retentionTime / 60
                ms2Indx <- ms2Indx[metaData$retentionTime[ms2Indx] >=
                                    minRt & metaData$retentionTime[
                                        ms2Indx] <= maxRt]
                # if fixed ions detected
                fixedDetectIons <- lapply(mzR::peaks(ms2File, ms2Indx),
                          function(x) {
                              detIons <- lapply(fixedIons,
                                                function(y) {
                                  # calculate local noise level
                                  localArTmp <-
                                      x[x[, 1] > {
                                          y - {
                                              noiseBin / 2
                                          }
                                      } & x[, 1] <
                                      {
                                          y + {
                                              noiseBin / 2
                                          }
                                      }, 2]
                                  noiseLevTmp <- median(
                                      localArTmp)
                                  snrTmp <- x[, 2] / noiseLevTmp
                                  inMass <-
                                      which({
                                          abs(x[, 1] - y) <= binSizeMS2
                                      } &
                                      {
                                          snrTmp >= minSNR
                                      })
                                  if (length(inMass) > 0) {
                                    inMass <- inMass[which.max(
                                    snrTmp[inMass])]
                                    names(inMass) <- round(snrTmp[inMass], 2)
                                  }
                                  return(inMass)
                              })
                          })#{{x[, 2]/max(x[, 2])} * 100} >= minRelFixed))})
                names(fixedDetectIons) <- ms2Indx
                fixedDetectIndx <- unlist(lapply(fixedDetectIons,
                                                 function(x)
                                             sum(
                                                 vapply(x, function(y)
                                                     length(y) > 0,
                                                     FUN.VALUE = logical(1))
                                             ) >= minFixed))
                fixedDetectIndx <- ms2Indx[fixedDetectIndx]
                # calc sim of top 5 base peaks to example spec
                dpMat <- do.call(cbind, lapply(mzR::peaks(ms2File,
                fixedDetectIndx), function(x) {
                              mostIntIdx <- order(x[, 2], decreasing = TRUE)[
                                  seq_len(topIons)]
                              x <- x[mostIntIdx, , drop = FALSE]
                              x <- x[order(x[, 1]),]
                              topIntIdx <- order(x[, 2], decreasing = TRUE)
                              top5 <- x[topIntIdx,]
                              # more than 10 daltons from base peak
                              topIntIdx <-
                                  topIntIdx[c(1, which(abs(top5[1, 1] -
                                      top5[, 1]) > 10)[seq_len({
                                                    topIntIt - 1
                                                    })])]
                              # generate scaled spec
                              diffSpecTop5 <- vapply(topIntIdx, function(y) {
                                  diffTmp <- x[, 1] - x[y, 1]
                                  binTmp <-
                                      cut(
                                          diffTmp,
                                          breaks = seq(minMass, maxMass,
                                                       binSizeMS2),
                                          labels = labelsTmp
                                      )
                                  sumTmp <- tapply(x[, 2], binTmp, max)
                                  massTmp <- tapply(x[, 1], binTmp, mean)
                                  resTmp <- c(massTmp, sumTmp)
                              }, FUN.VALUE = numeric(c(nlevels(
                                  cut(
                                      rnorm(1, length(x)),
                                      breaks = seq(minMass, maxMass,binSizeMS2)
                                  )
                              ) *
                                  2)))
                                                  }))
                colnames(dpMat) <-
                    paste0(rep(fixedDetectIndx,
                               each = topIntIt),
                           '_',
                           rep(seq_len(topIntIt), length(fixedDetectIndx)))
                dpMat[is.na(dpMat)] <- 0
                massesMat <-
                    dpMat[seq_len(length(massBinsDiffSpec)),]
                dpMat <-
                    dpMat[-{
                        seq_len(length(massBinsDiffSpec))
                    },]
                #subset to only include ions in model spec
                massesMat <- massesMat[modelSpecIndx,]
                dpMat <- dpMat[modelSpecIndx,]
                dpMat <- cbind(modelSpec = massBinsDiffSpec[modelSpecIndx],
                            dpMat)
                dotProdMat <- crossprod(dpMat)
                sqrtMatrixTmp <-
                    matrix(
                        sqrt(colSums(dpMat ^ 2)),
                        nrow = nrow(dotProdMat),
                        ncol = ncol(dotProdMat),
                        byrow = TRUE
                    )
                dotProdMat <-
                    dotProdMat / (sqrtMatrixTmp * diag(sqrtMatrixTmp))
                dpIndx <- dotProdMat[, 1] >= minDotProd
                propIded <-
                {
                    apply(dpMat, 2, function(x)
                        sum(x != 0))
                } >= minVarPeaks
                dpIndx <- which(dpIndx & propIded)
                if (length(dpIndx) > 1) {
                    spec2Plot <- colnames(dpMat)[dpIndx[-1]]
                    scanIdx <-
                        as.numeric(gsub('_.+', '', spec2Plot))
                    maxDp <- tapply(spec2Plot, scanIdx,
                                    function(x)
                                        x[which.max(dotProdMat[1, x])])
                    scanIdx <- as.numeric(names(maxDp))
                    resTable <- metaData[scanIdx,]
                    # keep specific columns
                    resTable <- resTable[, c(
                        'precursorMZ',
                        'adjPrecursorMZ',
                        'ppmDrift',
                        'retentionTime',
                        'precursorScanNum',
                        'seqNum',
                        'totIonCurrent',
                        'precursorCharge'
                    )]
                    colnames(resTable) <-
                        c(
                            'MIM',
                            'adjMIM',
                            'ppmDrift',
                            'RT',
                            'precursorScanNum',
                            'MS2ScanNum',
                            'TIC',
                            'precursorCharge'
                        )
                    resTable$peptide <- pepTmp
                    resTable$MS2FileName <-
                        basename(ms2Files[i])
                    resTable$plotURL <- ''
                    # resTable$plotHyperLink <- ''
                    resTable$meanSNRFixed <- 0
                    resTable$SNRFixed <- 0
                    resTable$meanSNRVar <- 0
                    resTable$SNRVar <- 0
                    resTable$nFixedIons <- 0
                    resTable$nVarIons <- 0
                    resTable$fixedDetected <- ''
                    resTable$varDetected <- ''
                    resTable$dotProdScore <-
                        dotProdMat[1, maxDp]
                    resTable$adjRT <-
                        predict(rtDevModels[[i]], newdata = resTable$RT)
                    resTable$adjRT <-
                        resTable$RT - resTable$adjRT
                    adductMasses <- mapply(
                        OrgMassSpecR::MonoisotopicMass,
                        charge = unique(resTable$precursorCharge),
                        MoreArgs = list(formula = pepForm)
                    )
                    names(adductMasses) <- unique(resTable$precursorCharge)
                    resTable$adductMass <- 
                    {
                        {   resTable$MIM -
                            { 
                                1.00782504 / resTable$precursorCharge
                            }
                        } -
                        adductMasses[as.character(
                            resTable$precursorCharge)]
                    } * 
                    resTable$precursorCharge
                    resTable$massUnModPep <- adductMasses[as.character(
                        resTable$precursorCharge)]
                    # subset for adduct mass
                    massIdx <- resTable$adjMIM >= minMz &
                        resTable$adjMIM <= maxMz
                    resTable <-
                        resTable[massIdx, , 
                        drop = FALSE]
                    maxDp <- maxDp[massIdx]
                    scanIdx <- scanIdx[massIdx]
                    resTable$fixedIonIdx <- ''
                    resTable$varIonIdx <- ''
                    resTable$fixedIonIntTmp <- 0
                    resTable$varIonIntTmp <- 0
                    resTable$fixedIonMzTmp <- 0
                    resTable$varIonMzTmp <- 0
                    resTable$relIntBPScore <- 0
                    resTable$propEx <- 0
                    # multiply by charge state
                    resTable <-
                        resTable[, c(
                            'MIM',
                            'adjMIM',
                            'ppmDrift',
                            'RT',
                            'adjRT',
                            'massUnModPep',
                            'adductMass',
                            'precursorCharge',
                            'MS2FileName',
                            'plotURL',
                            'peptide',
                            'precursorScanNum',
                            'MS2ScanNum',
                            'TIC',
                            'dotProdScore',
                            'meanSNRFixed',
                            'SNRFixed',
                            'meanSNRVar',
                            'SNRVar',
                            'nFixedIons',
                            'nVarIons',
                            'fixedDetected',
                            'varDetected',
                            'fixedIonIdx',
                            'varIonIdx',
                            'fixedIonIntTmp',
                            'varIonIntTmp',
                            'fixedIonMzTmp',
                            'varIonMzTmp',
                            'relIntBPScore',
                            'propEx'
                        )]
                    outDir <-
                        paste0(gsub('\\.mzXML$', '', ms2Files[i]),
                               '_adductID')
                    suppressWarnings(dir.create(outDir))
                    outDir <- paste0(outDir, '/', pepTmp)
                    suppressWarnings(dir.create(outDir))
                    for (j in seq_len(length(maxDp))) {
                        specTmp <- mzR::peaks(ms2File, scanIdx[j])
                        iTmp <- dpMat[, maxDp[j]]
                        peaksDet <- iTmp != 0
                        iTmp <- iTmp[peaksDet]
                        fragTypes <-
                            modelSpec$ionType[peaksDet]
                        varPeakIdx <-
                            which(specTmp[, 2] %in% iTmp)
                        # calculate SNR y series
                        snrVarPeaks <-
                            unlist(lapply(varPeakIdx,
                                          function(x) {
                                              # calculate local noise level
                                              localArTmp <-
                                                  specTmp[specTmp[, 1] >
                                                  {
                                                      specTmp[x, 1] -
                                                      {
                                                          noiseBin / 2
                                                      }
                                                  } & specTmp[, 1] < {
                                                      specTmp[x, 1] + {
                                                          noiseBin / 2
                                                      }
                                                  }, 2]
                                              noiseLevTmp <-
                                                  median(localArTmp)
                                              snrTmp <-
                                                  specTmp[x, 2] / noiseLevTmp
                                          }))
                        varSNRIdx <-
                            snrVarPeaks >= minSNR
                        varPeakIdx <-
                            varPeakIdx[varSNRIdx]
                        snrVarPeaks <-
                            snrVarPeaks[varSNRIdx]
                        
                        if (length(varPeakIdx) >= minVarPeaks) {
                            resTable$varDetected[j] <- paste0(fragTypes[
                                varSNRIdx],collapse = '; ')
                            resTable$nVarIons[j] <-
                                sum(varSNRIdx)
                            resTable$varIonIdx[j] <-
                                paste0(varPeakIdx,
                                       collapse = '; ')
                            resTable$varIonIntTmp[j] <-
                                paste0(specTmp[varPeakIdx, 2],
                                       collapse = '; ')
                            resTable$varIonMzTmp[j] <-
                                paste0(specTmp[varPeakIdx, 1],
                                       collapse = '; ')
                            resTable$meanSNRVar[j] <-
                                round(mean(snrVarPeaks), 2)
                            resTable$SNRVar[j] <-
                                paste0(round(snrVarPeaks, 2),
                                       collapse = '; ')
                            relIntBPScore <-
                                max(specTmp[varPeakIdx, 2],
                                     na.rm=TRUE) / max(specTmp[, 2]
                                         , na.rm=TRUE)
                            suppressWarnings(dir.create(paste0(outputPlotDir, 
                                '/scanplots/')))
                            plotName <- paste0( 
                                '/spectrumGroups/',
                                basename(ms2Files[i]),
                                ' scan ',
                                resTable$MS2ScanNum[j],
                                ' M',
                                round(resTable$MIM[j], 4),
                                '_RT',
                                round(resTable$RT[j], 2),
                                ' dp ',
                                round(dotProdMat[1,
                                                 maxDp[j]], 2),
                                ' varPeakDet ',
                                sum(peaksDet)
                            )
                            fixIonTmp <-
                                fixedDetectIons[[as.character(scanIdx[j])]]
                            names(fixIonTmp) <- paste0(
                                fixedIonNames, '_', names(fixIonTmp))
                            if (is.list(fixIonTmp)) {
                                fixIonTmp <- unlist(fixIonTmp)
                            }
                            snrFixTmp <- as.numeric(gsub(
                                '.+_\\.|.+_', '', names(fixIonTmp)))
                            resTable$plotURL[j] <-
                                paste0(outputPlotDir,
                                       '/', plotName, '.png')
                            #png
                            if (!is.null(outputPlotDir)) {
                            png(
                                resTable$plotURL[j],
                                width = 2200,
                                height = 2000,
                                res = 275
                            )
                            specTmp <-
                                specTmp[order(specTmp[, 1]),]
                            par(bg = 'white', fg = 'black')
                            suppressMessages(
                                plot(
                                    specTmp,
                                    type = 'h',
                                    xlab = 'm/z',
                                    ylab = 'intensity',
                                    main = plotName,
                                    col = 'darkgrey'
                                )
                            )
                            }
                            fixTmp <-
                                specTmp[fixIonTmp, ,
                                        drop = FALSE]
                            rownames(fixTmp) <-
                                gsub('\\+.+', '+',
                                     names(fixIonTmp))
                            #label fixed ions
                            if (!is.null(outputPlotDir)) {
                            points(
                                fixTmp,
                                col = 'blue',
                                type = 'h',
                                lwd = 2
                            )
                            text(
                                fixTmp,
                                labels =
                                    row.names(fixTmp),
                                col = 'blue',
                                pos = 3,
                                srt = 90
                            )
                            }
                            resTable$fixedDetected[j] <-
                                paste0(names(fixTmp), collapse = '; ')
                            resTable$nFixedIons[j] <-
                                length(fixTmp)
                            resTable$fixedIonIdx[j] <-
                                paste0(fixIonTmp, collapse = '; ')
                            resTable$fixedIonIntTmp[j] <-
                                paste0(specTmp[fixIonTmp, 2], collapse = '; ')
                            snrFixTmp <- as.numeric(gsub('.+_\\.', '',
                                                         names(fixIonTmp)))
                            resTable$fixedIonMzTmp[j] <-
                                paste0(specTmp[fixIonTmp, 1], collapse = '; ')
                            resTable$meanSNRFixed[j] <-
                                round(mean(snrFixTmp), 2)
                            resTable$SNRFixed[j] <-
                                paste0(round(snrFixTmp, 2), collapse = '; ')
                            resTable$relIntBPScore[j] <-
                                max(c(
                                    relIntBPScore,
                                    max(specTmp[fixIonTmp, 2], na.rm =
                                            TRUE) /
                                        max(specTmp[, 2], na.rm = TRUE)
                                ))
                            sumIntEx <-
                                sum(specTmp[unique(c(
                                    varPeakIdx, fixIonTmp)), 2])
                            resTable$propEx[j] <-
                                sumIntEx /
                                sum(specTmp[, 2])
                            #label variable ions
                            varTmp <- specTmp[varPeakIdx, , drop =
                                                  FALSE]
                            rownames(varTmp) <- gsub('\\+.+', '+',
                                                     (fragTypes[varSNRIdx]))
                            if (!is.null(outputPlotDir)) {
                            points(
                                varTmp,
                                col = 'red',
                                type = 'h',
                                lwd = 2
                            )
                            text(
                                varTmp,
                                labels =
                                    fragTypes,
                                col = 'red',
                                srt = 90,
                                pos = 3
                            )
                            legend(
                                'topright',
                                c('fixed', 'variable'),
                                lwd = c(2, 2),
                                col = c('blue', 'red')
                            )
                            dev.off()
                        
                        }
                        }
                    }
            
            
        
                    # filter by minimum total spectrum signal to
                    #noise
                    write.csv(
                        resTable,
                        paste0(
                            outDir,
                            '/adductID_',
                            pepTmp,
                            '_',
                            gsub('\\.mzXML$', '',
                                 basename(ms2Files[i])),
                            '.csv'
                        ),
                        row.names = FALSE
                    )

                    resTable <-
                        resTable[resTable$meanSNRVar != 0,, drop = FALSE]
                    
                
                    }
            
                        
                    allResults <-
                        rbind(allResults, resTable)
            }        
            
                    
        # MULTITHREADED
        } else { 
            message(paste0("Starting SNOW cluster with ", nCores,
                           " local sockets..."))
            cl <- parallel::makeCluster(nCores, outfile = '')
            doSNOW::registerDoSNOW(cl)
            
            message("Identifying ",
                    pepTmp,
                    " spectra in ",
                    length(ms2Files),
                    " MS2 files.\n")
            
            progress <- function(n)
                message(paste0(
                    n,
                    ' of ',
                    length(ms2Files),
                    ' complete (',
                    basename(ms2Files)[n],
                    ').\n'
                ))
            opts <- list(progress = progress)
            # foreach and dopar from foreach package
            allResults <-
                foreach(
                    i = seq_len(length(ms2Files)),
                    .packages = c('mzR', 'data.table'),
                    .options.snow = opts,
                    .verbose = FALSE
                ) %dopar% {
                    ms2File <- mzR::openMSfile(ms2Files[i])
                    metaData <- mzR::header(ms2File)
                    if (file.exists(massDriftFiles[i])) {
                        massDriftTmp <- as.data.frame(
                            data.table::fread(
                                massDriftFiles[i],
                                header = TRUE,
                                stringsAsFactors = FALSE
                            )
                        )
                        metaData <- cbind(metaData,
                                        massDriftTmp[, c(
                                            'ppmDrift', 'adjPrecursorMZ')])
                    } else {
                        metaData$ppmDrift <- NA
                        metaData$adjPrecursorMZ <-
                            metaData$precursorMZ
                    }
                    ms2Indx <- which(metaData$msLevel == 2 &
                                        metaData$precursorCharge %in% precCh)
                    # rt range
                    metaData$retentionTime <-
                        metaData$retentionTime / 60
                    ms2Indx <-
                        ms2Indx[metaData$retentionTime[ms2Indx]
                                >= minRt &
                                    metaData$retentionTime[ms2Indx] <= maxRt]
                    # if fixed ions detected
                    fixedDetectIons <-
                        lapply(mzR::peaks(ms2File,
                                          ms2Indx),
                       function(x) {
                           detIons <- lapply(fixedIons, function(y) {
                               # calculate local noise level
                               localArTmp <- x[x[, 1] >
                               {
                                   y - {
                                       noiseBin / 2
                                   }
                               } & x[, 1] <
                               {
                                   y + {
                                       noiseBin / 2
                                   }
                               }, 2]
                               noiseLevTmp <- median(localArTmp)
                               snrTmp <- x[, 2] / noiseLevTmp
                               inMass <- which({
                                   abs(x[, 1] - y) <=
                                       binSizeMS2
                               } & {
                                   snrTmp >= minSNR
                               })
                               if (length(inMass) > 0) {
                                   inMass <- inMass[which.max(snrTmp[inMass])]
                                   names(inMass) <-
                                       round(snrTmp[inMass], 2)
                               }
                               return(inMass)
                           })
                       })
                    names(fixedDetectIons) <- ms2Indx
                    fixedDetectIndx <- unlist(lapply(fixedDetectIons,
                                     function(x)
                                         sum(
                                             vapply(x, function(y)
                                                 length(y)
                                                 > 0, FUN.VALUE = logical(1))
                                         ) >= minFixed))
                    fixedDetectIndx <-
                        ms2Indx[fixedDetectIndx]
                    # calc sim of top 5 base peaks to example spec
                    dpMat <-
                        do.call(cbind, lapply(mzR::peaks(ms2File,
             fixedDetectIndx), function(x) {
                 mostIntIdx <- order(x[, 2], decreasing =
                                         TRUE)[seq_len(topIons)]
                 x <- x[mostIntIdx, , drop = FALSE]
                 x <- x[order(x[, 1]),]
                 topIntIdx <-
                     order(x[, 2], decreasing =
                               TRUE)
                 top5 <- x[topIntIdx,]
                 # more than 10 daltons from base peak
                 topIntIdx <-
                     topIntIdx[c(1, which(abs(top5[1, 1] - top5[, 1]) >
                                              10)[seq_len({
                                                  topIntIt - 1
                                              })])]
                 diffSpecTop5 <-
                     vapply(topIntIdx,
                            function(y) {
                                diffTmp <- x[, 1] - x[y, 1]
                                binTmp <-
                                    cut(
                                        diffTmp,
                                        breaks = seq(minMass,
                                                     maxMass, binSizeMS2),
                                        labels = labelsTmp
                                    )
                                sumTmp <-
                                    tapply(x[, 2], binTmp, max)
                                massTmp <-
                                    tapply(x[, 1], binTmp, mean)
                                resTmp <- c(massTmp, sumTmp)
                            }, FUN.VALUE = numeric(c(nlevels(
                                cut(
                                    rnorm(1, length(x)),
                                    breaks = seq(minMass, maxMass, binSizeMS2)
                                )
                            ) *
                                2)))
                                                         }))
                    colnames(dpMat) <- paste0(rep(fixedDetectIndx,
                                                each = topIntIt),
                                            '_',rep(seq_len(topIntIt), length(
                                                fixedDetectIndx)))
                    dpMat[is.na(dpMat)] <- 0
                    massesMat <-
                        dpMat[seq_len(length(massBinsDiffSpec)),]
                    dpMat <- dpMat[-{
                        seq_len(length(massBinsDiffSpec))
                    },]
                    #subset to only include ions in model spec
                    massesMat <-
                        massesMat[modelSpecIndx,]
                    dpMat <- dpMat[modelSpecIndx,]
                    dpMat <-
                        cbind(modelSpec = massBinsDiffSpec[
                            modelSpecIndx], dpMat)
                    
                    dotProdMat <- crossprod(dpMat)
                    sqrtMatrixTmp <-
                        matrix(
                            sqrt(colSums(dpMat ^ 2)),
                            nrow = nrow(dotProdMat),
                            ncol = ncol(dotProdMat),
                            byrow = TRUE
                        )
                    
                    dotProdMat <-
                        dotProdMat / (sqrtMatrixTmp *
                                          diag(sqrtMatrixTmp))
                    dpIndx <-
                        dotProdMat[, 1] >= minDotProd
                    propIded <-
                    {
                        apply(dpMat, 2, function(x)
                            sum(x != 0))
                    } >= minVarPeaks
                    dpIndx <-
                        which(dpIndx & propIded)
                    if (length(dpIndx) > 1) {
                        spec2Plot <- colnames(dpMat)[dpIndx[-1]]
                        scanIdx <-
                            as.numeric(gsub('_.+', '',
                                            spec2Plot))
                        maxDp <-
                            tapply(spec2Plot, scanIdx,
                                   function(x)
                                       x[which.max(dotProdMat[1,
                                                              x])])
                        scanIdx <-
                            as.numeric(names(maxDp))
                        resTable <-
                            metaData[scanIdx,]
                        # keep specific columns
                        resTable <-
                            resTable[, c(
                                'precursorMZ',
                                'adjPrecursorMZ',
                                'ppmDrift',
                                'retentionTime',
                                'precursorScanNum',
                                'seqNum',
                                'totIonCurrent',
                                'precursorCharge'
                            )]
                        colnames(resTable) <-
                            c(
                                'MIM',
                                'adjMIM',
                                'ppmDrift',
                                'RT',
                                'precursorScanNum',
                                'MS2ScanNum',
                                'TIC',
                                'precursorCharge'
                            )
                        resTable$peptide <- pepTmp
                        resTable$MS2FileName <-
                            basename(ms2Files[i])
                        resTable$plotURL <- ''
                        #resTable$plotHyperLink <- ''
                        resTable$meanSNRFixed <- 0
                        resTable$SNRFixed <- 0
                        resTable$meanSNRVar <- 0
                        resTable$SNRVar <- 0
                        resTable$nFixedIons <- 0
                        resTable$nVarIons <- 0
                        resTable$fixedDetected <- ''
                        resTable$varDetected <- ''
                        resTable$dotProdScore <-
                            dotProdMat[1, maxDp]
                        resTable$adjRT <- predict(rtDevModels[[i]],
                                                  newdata = resTable$RT)
                        resTable$adjRT <-
                            resTable$RT -
                            resTable$adjRT
                        adductMasses <-
                            mapply(
                                OrgMassSpecR::MonoisotopicMass,
                                charge = unique(resTable$precursorCharge),
                                MoreArgs = list(formula = pepForm)
                            )
                        names(adductMasses) <-
                            unique(resTable$precursorCharge)
                        resTable$adductMass <-
                        {
                            {
                                resTable$MIM -
                                {
                                    1.00782504 / resTable$precursorCharge
                                }
                            } -
                                adductMasses[as.character(
                                    resTable$precursorCharge)]
                        } *
                            resTable$precursorCharge
                        resTable$massUnModPep <-
                            adductMasses[as.character(resTable$precursorCharge)]
                        # subset for adduct mass
                        massIdx <-
                            resTable$adjMIM >=
                            minMz &
                            resTable$adjMIM <= maxMz
                        resTable <-
                            resTable[massIdx, ,
                                     drop = FALSE]
                        maxDp <- maxDp[massIdx]
                        scanIdx <- scanIdx[massIdx]
                        resTable$fixedIonIdx <- ''
                        resTable$varIonIdx <- ''
                        resTable$fixedIonIntTmp <- 0
                        resTable$varIonIntTmp <- 0
                        resTable$fixedIonMzTmp <- 0
                        resTable$varIonMzTmp <- 0
                        resTable$relIntBPScore <- 0
                        resTable$propEx <- 0
                        # multiply by charge state
                        resTable <-
                            resTable[, c(
                                'MIM',
                                'adjMIM',
                                'ppmDrift',
                                'RT',
                                'adjRT',
                                'massUnModPep',
                                'adductMass',
                                'precursorCharge',
                                'MS2FileName',
                                'plotURL',
                                'peptide',
                                'precursorScanNum',
                                'MS2ScanNum',
                                'TIC',
                                'dotProdScore',
                                'meanSNRFixed',
                                'SNRFixed',
                                'meanSNRVar',
                                'SNRVar',
                                'nFixedIons',
                                'nVarIons',
                                'fixedDetected',
                                'varDetected',
                                'fixedIonIdx',
                                'varIonIdx',
                                'fixedIonIntTmp',
                                'varIonIntTmp',
                                'fixedIonMzTmp',
                                'varIonMzTmp',
                                'relIntBPScore',
                                'propEx'
                            )]
                        # create result directory
                        outDir <-
                            paste0(gsub('\\.mzXML$', '',
                                        ms2Files[i]), '_adductID')
                        #outDir <- gsub('/', '\\\\', outDir)
                        suppressWarnings(dir.create(outDir))
                        outDir <-
                            paste0(outDir, '/', pepTmp)
                        suppressWarnings(dir.create(outDir))
                        for (j in seq_len(length(maxDp))) {
                            specTmp <- mzR::peaks(ms2File, scanIdx[j])
                            iTmp <- dpMat[, maxDp[j]]
                            peaksDet <- iTmp != 0
                            iTmp <- iTmp[peaksDet]
                            fragTypes <- modelSpec$ionType[peaksDet]
                            varPeakIdx <- which(specTmp[, 2] %in% iTmp)
                            # calculate SNR y series
                            snrVarPeaks <- unlist(lapply(varPeakIdx,
                                     function(x) {
                                         # calculate local noise level
                                         localArTmp <- specTmp[specTmp[, 1] >
                                         {
                                             specTmp[x, 1] - {
                                                 noiseBin / 2
                                             }
                                         } &
                                             specTmp[, 1] <
                                             {
                                                 specTmp[x, 1] +
                                                 {
                                                     noiseBin / 2
                                                 }
                                             }, 2]
                                         noiseLevTmp <- median(localArTmp)
                                         snrTmp <- specTmp[x, 2] /
                                             noiseLevTmp
                                     }))
                            varSNRIdx <- snrVarPeaks >=
                                minSNR
                            varPeakIdx <- varPeakIdx[varSNRIdx]
                            snrVarPeaks <- snrVarPeaks[varSNRIdx]
                            if (length(varPeakIdx) >
                                minVarPeaks) {
                                resTable$varDetected[j] <-
                                    paste0(fragTypes[varSNRIdx], collapse =
                                               '; ')
                                resTable$nVarIons[j] <-
                                    sum(varSNRIdx)
                                resTable$varIonIdx[j] <-
                                    paste0(varPeakIdx,
                                           collapse = '; ')
                                resTable$varIonIntTmp[j] <-
                                    paste0(specTmp[varPeakIdx, 2], collapse
                                           = '; ')
                                resTable$varIonMzTmp[j] <-
                                    paste0(specTmp[varPeakIdx, 1], collapse =
                                               '; ')
                                resTable$meanSNRVar[j] <-
                                    round(mean(snrVarPeaks), 2)
                                resTable$SNRVar[j] <- paste0(
                                    round(snrVarPeaks, 2),collapse = '; ')
                                relIntBPScore <- max(specTmp[varPeakIdx, 2],
                                                     na.rm = TRUE) /
                                    max(specTmp[, 2],
                                        na.rm = TRUE)
                                        
                                        
                            suppressWarnings(dir.create(paste0(outputPlotDir, 
                                            '/scanplots/')))
                                plotName <- paste0( 
                                    '/spectrumGroups/',
                                    basename(ms2Files[i]),
                                    ' scan ',
                                    resTable$MS2ScanNum[j],
                                    ' M',
                                    round(resTable$MIM[j],
                                          4),
                                    '_RT',
                                    round(resTable$RT[j], 2),
                                    ' dp ',
                                    round(dotProdMat[1,
                                                     maxDp[j]], 2),
                                    ' varPeakDet ',
                                    sum(peaksDet)
                                )
                                fixIonTmp <-
                                    fixedDetectIons[[as.character(scanIdx[j])]]
                                names(fixIonTmp) <-
                                    paste0(fixedIonNames,
                                           '_', names(fixIonTmp))
                                if (is.list(fixIonTmp)) {
                                    fixIonTmp <- unlist(fixIonTmp)
                                }
                                snrFixTmp <- as.numeric(gsub('.+_\\.|.+_',
                                            '', names(fixIonTmp)))
                                resTable$plotURL[j] <-
                                    paste0(outputPlotDir, '/',
                                           plotName, '.png')
                                # png
                                if (!is.null(outputPlotDir)) {
                                png(
                                    resTable$plotURL[j],
                                    width = 2200,
                                    height = 2000,
                                    res = 275
                                )
                                specTmp <- specTmp[order(specTmp[, 1]),]
                                par(bg = 'white', fg = 'black')
                                suppressMessages(
                                    plot(
                                        specTmp,
                                        type = 'h',
                                        xlab = 'm/z',
                                        ylab = 'intensity',
                                        main = plotName,
                                        col = 'darkgrey'
                                    )
                                )
                            }
                                fixTmp <- specTmp[fixIonTmp, , drop =
                                                      FALSE]
                                rownames(fixTmp) <- gsub('\\+.+', '+',
                                                         names(fixIonTmp))
                                if (!is.null(outputPlotDir)) {                         
                                #label fixed ions
                                points(
                                    fixTmp,
                                    col = 'blue',
                                    type = 'h',
                                    lwd = 2
                                )
                                text(
                                    fixTmp,
                                    labels =
                                        row.names(fixTmp),
                                    col = 'blue',
                                    pos = 3,
                                    srt = 90
                                )
                            }
                                resTable$fixedDetected[j] <-
                                    paste0(names(fixTmp),
                                           collapse = '; ')
                                resTable$nFixedIons[j] <-
                                    length(fixTmp)
                                resTable$fixedIonIdx[j] <-
                                    paste0(fixIonTmp,
                                           collapse = '; ')
                                resTable$fixedIonIntTmp[j] <-
                                    paste0(specTmp[fixIonTmp,
                                                   2], collapse = '; ')
                                resTable$fixedIonMzTmp[j] <-
                                    paste0(specTmp[fixIonTmp, 1],
                                           collapse = '; ')
                                snrFixTmp <- as.numeric(gsub('.+_\\.', '',
                                                             names(fixIonTmp)))
                                resTable$meanSNRFixed[j] <-
                                    round(mean(snrFixTmp), 2)
                                resTable$SNRFixed[j] <-
                                    paste0(round(snrFixTmp, 2),
                                           collapse = '; ')
                                resTable$relIntBPScore[j] <-
                                    max(c(
                                        relIntBPScore,
                                        max(specTmp[fixIonTmp, 2],
                                            na.rm = TRUE) /
                                            max(specTmp[, 2],
                                                na.rm = TRUE)
                                    ))
                                sumIntEx <- sum(specTmp[unique(c(varPeakIdx,
                                            fixIonTmp)), 2])
                                resTable$propEx[j] <-
                                    sumIntEx /
                                    sum(specTmp[, 2])
                                #label variable ions
                                varTmp <- specTmp[varPeakIdx, , drop =
                                                      FALSE]
                                rownames(varTmp) <- gsub('\\+.+', '+',
                                                        (fragTypes[varSNRIdx]))
                                if (!is.null(outputPlotDir)) {                
                                points(
                                    varTmp,
                                    col = 'red',
                                    type = 'h',
                                    lwd = 2
                                )
                                text(
                                    varTmp,
                                    labels =
                                        fragTypes,
                                    col = 'red',
                                    srt = 90,
                                    pos = 3
                                )
                                legend(
                                    'topright',
                                    c('fixed', 'variable'),
                                    lwd = c(2, 2),
                                    col = c('blue', 'red')
                                )
                                dev.off()
                            }
                            }
                        }
                        # save results table
                        write.csv(
                            resTable,
                            paste0(
                                outDir,
                                '/adductID_',
                                pepTmp,
                                '_',
                                gsub('\\.mzXML$', '',
                                     basename(ms2Files[i])),
                                '.csv'
                            ),
                            row.names = FALSE
                        )
                        # filter by minimum total spectrum
                        #signal to
                        #noise
                        resTable <- resTable[resTable$meanSNRVar != 0, ,
                                             drop = FALSE]
                        return(resTable)
                    }
                }
            # stop SNOW cluster
            parallel::stopCluster(cl)
            allResults <- do.call(rbind, allResults)
        }

        # Spectrum grouping and formula generation
        message(
            paste0(
                'Grouping spectra and
                identifying possible atomic
                formulae from adduct mass.\n'
            )
            )
        message(
            paste0(
                'Hierarchically clustering
                spectra and grouping based
                on a mass error of ',
                groupMzabs,
                ' and a retention time deviation of ',
                groupRtDev,
                '.\n'
            )
        )
        if (is.null(allResults$adjMIM)) {
            hr <- fastcluster::hclust.vector(allResults$MIM, 
                                             metric = 'euclidean',
                                             method = 'median')
        } else {
            hr <- fastcluster::hclust.vector(allResults$adjMIM,
                                             metric = 'euclidean',
                                             method = 'median')
        }
        allResults$msPrecursor_group <- 0
        hval <- groupMzabs
        nclust <- (nrow(allResults) - 1):1 # Number of clusters at each step
        # Find the number of clusters that first
        #exceeds hval
        k <- max(nclust[which(diff(hr$height < hval) ==
                                  -1) + 1], na.rm=TRUE)
        allResults$msPrecursor_group <- cutree(hr,
                                               k = k)
        
        # sort by mass group
        allResults <- allResults[order(allResults$msPrecursor_group),]
        # rt deviation
        interMSMSrtGroups <- tapply(allResults$adjRT,
                                    allResults$msPrecursor_group, function(x) {
                                        if (length(x) > 1) {
                                            cutree(
                                                fastcluster::hclust.vector(
                                                    x, metric = 'euclidean',
                                                    method = 'median'),
                                                    h = groupRtDev
                                            )
                                        } else {
                                            1
                                        }
                                    })
        # interMSMS groups
        allResults$MSMSGroups <- paste0(allResults$msPrecursor_group,
                                        "_",
                                        unlist(interMSMSrtGroups))
        allResults$propExScore <-
            allResults$propEx /
            max(allResults$propEx, na.rm = TRUE)
        # 3. number of fixed and var
        #identified 1= most 0=lowest
        allResults$nFixVarScore <-allResults$nFixedIons + allResults$nVarIons
        allResults$nFixVarScore <- allResults$nFixVarScore /
            max(allResults$nFixVarScore, na.rm = TRUE)
        # 4. number of consecutive fixed and
        #var ions
        consecFixed <- unlist(lapply(strsplit(allResults$fixedDetected, '; '),
                                     function(x) {
                                         ionNos <- sort(as.numeric(
                                            gsub('\\[[A-z]|\\].+', '', x)))
                                            nConsec <- sum(diff(ionNos) == 1)
                                     }))
        consecVar <- unlist(lapply(strsplit(allResults$varDetected,
                                            '; '),
                                   function(x) {
                                       ionNos <- sort(as.numeric(
                                        gsub('\\[[A-z]|\\].+', '', x)))
                                        nConsec <- sum(diff(ionNos) == 1)
                                   }))
        allResults$consecIonsScore <-
            consecFixed + consecVar
        allResults$consecIonsScore <-
            allResults$consecIonsScore /
            max(allResults$consecIonsScore,
                na.rm = TRUE)
        # 5. snr fixed and var ions
        allResults$meanSNRScore <-
            allResults$meanSNRFixed +
            allResults$meanSNRVar
        allResults$meanSNRScore <-
            allResults$meanSNRScore /
            max(allResults$meanSNRScore,
                na.rm = TRUE)
        # 6. dot product with model spec (just as is)?
        # allResults$dotProdScore
        # calculate mean score
        allResults$idScore <-
            rowMeans(allResults[, grep('Score$', colnames(allResults))],
                     na.rm = TRUE)
        # check the dot prod
        #similarity for each group
        dpIndxGroup <- lapply(unique(allResults$MSMSGroups),
          function(x) {
              idxTmp <-
                  allResults$MSMSGroups %in% x
              if (sum(idxTmp) == 1) {
                  dpIndx <- which(idxTmp)[1]
              } else if (sum(idxTmp) == 2) {
                  dpIndx <- which(idxTmp)[which.max(allResults$idScore[idxTmp])]
              } else {
                  fixMassTmp <-
                      suppressMessages(strsplit(
                          allResults$fixedIonMzTmp[idxTmp],'; '))
                  names(fixMassTmp) <- paste0(seq_len(length(fixMassTmp)), '_')
                  fixMassTmp <- unlist(fixMassTmp)
                  
                  fixIntTmp <-
                      suppressMessages(strsplit(
                        allResults$fixedIonIntTmp[idxTmp], '; '))
                  names(fixIntTmp) <- paste0(seq_len(length(fixIntTmp)), '_')
                  fixIntTmp <- unlist(fixIntTmp)
                  
                  fixSpec <- cbind(as.numeric(
                    fixMassTmp), as.numeric(fixIntTmp))
                  colnames(fixSpec) <- c('mass', 'int')
                  row.names(fixSpec) <- gsub('_.+', '', names(fixMassTmp))
                  
                  varMassTmp <- strsplit(allResults$varIonMzTmp[idxTmp], '; ')
                  names(varMassTmp) <- paste0(seq_len(length(varMassTmp)), '_')
                  varMassTmp <- unlist(varMassTmp)
                  
                  varIntTmp <- strsplit(allResults$varIonIntTmp[idxTmp], '; ')
                  names(varIntTmp) <- paste0(seq_len(length(varIntTmp)), '_')
                  varIntTmp <- unlist(varIntTmp)
                  
                  varSpec <- cbind(as.numeric(varMassTmp), 
                                as.numeric(varIntTmp))
                  colnames(varSpec) <- c('mass', 'int')
                  row.names(varSpec) <- gsub('_.+', '', names(varMassTmp))
                  
                  allSpecTmp <- rbind(fixSpec, varSpec)
                  suppressWarnings(allSpecTmp <- allSpecTmp[order(as.numeric(
                    row.names(allSpecTmp))),])
                  
                  dotProdMat <- suppressMessages(dotProdMatrix(
                    allSpecTmp, row.names(allSpecTmp), binSizeMS2 =
                        binSizeMS2))
                  dpIndx <- which(idxTmp)[rowMeans(
                      dotProdMat) >= minMeanSpecSim]
              }
              return(dpIndx)
          })
        dpIndxGroup <- unlist(dpIndxGroup)
        
        allResults <- allResults[dpIndxGroup, , drop = FALSE]
        
        freqTmp <- table(allResults$MSMSGroups)
        allResults$MSMSgroup_freq <- freqTmp[allResults$MSMSGroups]
        
        # calculate frequency score and recalculate average score
        # 7. frequency detected
        
        allResults <- allResults[!(is.na(allResults$MSMSgroup_freq)), ]
        
        allResults$freqDetScore <- allResults$MSMSgroup_freq /
            max(as.numeric(unlist(allResults$MSMSgroup_freq)), na.rm = TRUE)
        # calculate mean score
        allResults$idScore <- NULL
        allResults$idScore <- rowMeans(
            allResults[, grep('Score$',colnames(allResults))], na.rm = TRUE)
        # minimum id score
        minScIdx <- allResults$idScore >= minIdScore
        allResults <- allResults[minScIdx, , drop = FALSE]
        
        freqTmp <- table(allResults$MSMSGroups)
        allResults$MSMSgroup_freq <- freqTmp[allResults$MSMSGroups]
        
        
        message(length(unique(allResults$MSMSGroups)),
                ' MS/MS peak groups were identified.\n')
        
        # average MIM, rt and adduct mass
        if (is.null(allResults$adjMIM)) {
            avMass <- tapply(allResults$MIM, allResults$MSMSGroups, mean)
            allResults$MSMSgroup_averageMIM <- avMass[allResults$MSMSGroups]
            # ppm range
            ppmRange <-
                tapply(allResults$MIM, allResults$MSMSGroups, function(x) {
                    ppmTmp <- range(x)
                    round({
                        {
                            {
                                ppmTmp[2] - ppmTmp[1]
                            } / ppmTmp[2]
                        } * 1E06
                    }, 2)
                })
            allResults$ppmRange_minMaxGroup <- ppmRange[allResults$MSMSGroups]
            # ppm diff mean
            allResults$ppmDiff_groupMean <- round({
                {
                    allResults$MIM -
                        allResults$MSMSgroup_averageMIM
                } / allResults$MIM
            } * 1E06, 2)
        } else {
            avMass <- tapply(allResults$adjMIM, allResults$MSMSGroups, mean)
            allResults$MSMSgroup_averageMIM <- avMass[allResults$MSMSGroups]
            # ppm range
            ppmRange <-
                tapply(allResults$adjMIM, allResults$MSMSGroups, function(x) {
                    ppmTmp <- range(x)
                    round({
                        {
                            {
                                ppmTmp[2] - ppmTmp[1]
                            } / ppmTmp[2]
                        } * 1E06
                    }, 2)
                })
            allResults$ppmRange_minMaxGroup <-
                ppmRange[allResults$MSMSGroups]
            # ppm diff mean
            allResults$ppmDiff_groupMean <- round({
                {
                    allResults$adjMIM -
                        allResults$MSMSgroup_averageMIM
                } / allResults$adjMIM
            } * 1E06, 2)
        }
        avRt <- tapply(allResults$adjRT, allResults$MSMSGroups, mean)
        allResults$MSMSgroup_averageAdjRT <- avRt[allResults$MSMSGroups]
        rtRange <-
            tapply(allResults$adjRT, allResults$MSMSGroups, function(x) {
                rtRTmp <- range(x)
                return(paste0(round(rtRTmp[2] - rtRTmp[1], 2), collapse = '-'))
            })
        allResults$rtRange_minMaxGroup <- rtRange[allResults$MSMSGroups]
        avAddMass <-
            tapply(allResults$adductMass, allResults$MSMSGroups, mean)
        allResults$MSMSgroup_averageAdductMass <-
            avAddMass[allResults$MSMSGroups]
        # generate plots
        if (!is.null(outputPlotDir)) {
        message(paste0('Generating ', length(unique(
            allResults$MSMSGroups
        )),
        ' plots of spectrum groups'))
        
        groupDir <-
            paste0(outputPlotDir, '/spectrumGroups_', pepTmp, '/')
        suppressWarnings(dir.create(groupDir))
        allGroups <- unique(allResults$MSMSGroups)
        gPlotNames <- paste0(
            groupDir,
            '/spectrumGroup_',
            allGroups,
            '_M',
            round(avMass[allGroups], 3),
            'RT',
            round(avRt[allGroups], 2),
            '.png'
        )
        names(gPlotNames) <- allGroups
        allResults$MSMSgroup_plot <- gPlotNames[allResults$MSMSGroups]
        
        pb <- txtProgressBar(max = length(allGroups), style = 3)
        for (p in seq_len(length(allGroups))) {
            setTxtProgressBar(pb, p)
            png(
                gPlotNames[p],
                width = 2200,
                height = 2000,
                res = 275
            )
            par(mfrow = c(2, 1))
            tmpIndx <- allResults$MSMSGroups == allGroups[p]
            xlimTmp <-
                range(c(allResults$RT[tmpIndx], allResults$adjRT[tmpIndx]))
            suppressMessages(
                plot(
                    allResults$RT[tmpIndx],
                    allResults$MIM[tmpIndx],
                    ylab = 'm/z',
                    xlab = 'RT (mins)',
                    main = paste0(
                        'raw Rt (window = ',
                        paste0(round(diff(
                            range(allResults$RT[tmpIndx])
                        ), 2),
                        collapse = '-'),
                        ' mins)'
                    ),
                    xlim = xlimTmp
                )
            )
            suppressMessages(
                plot(
                    allResults$adjRT[tmpIndx],
                    allResults$MIM[tmpIndx],
                    ylab = 'm/z',
                    xlab = 'adj.RT (mins)',
                    main = paste0(
                        'corrected Rt (window = ',
                        paste0(round(diff(
                            range(allResults$adjRT[tmpIndx])
                        ), 2),
                        collapse = '-'),
                        ' mins)'
                    ),
                    xlim = xlimTmp
                )
            )
            dev.off()
        close(pb)
        }
        #
        # # plot all
        png(
            paste0(dirname(gPlotNames), '/allGroups.png'),
            width = 2200,
            height = 2000,
            res = 275
        )
        par(mfrow = c(2, 1))
        colsTmp <- as.numeric(as.factor(allResults$MSMSGroups))
        xlimTmp <- range(c(allResults$RT[!is.na(allResults$RT)],
                           allResults$adjRT[!is.na(allResults$adjRT)]))
        suppressMessages(
            plot(
                allResults$RT,
                allResults$MIM,
                col = colsTmp,
                ylab = 'm/z',
                xlab = 'RT (mins)',
                main = 'raw Rt',
                xlim = xlimTmp
            )
        )
        suppressMessages(
            plot(
                allResults$adjRT,
                allResults$MIM,
                col = colsTmp,
                ylab = 'm/z',
                xlab = 'adj.RT (mins)',
                main = 'corrected Rt',
                xlim = xlimTmp
            )
        )
        dev.off()
    }
        # write final results
        allResults$MSMSGroups <- NULL
        allResults$MSMSGroupsName <- paste0(
            'M',
            round(allResults$MSMSgroup_averageMIM, 4),
            '_RT',
            round(allResults$MSMSgroup_averageAdjRT, 2),
            '_n',
            allResults$MSMSgroup_freq
        )
        allResults$fixedIonMzTmp <- NULL
        allResults$fixedIonIntTmp <- NULL
        allResults$varIonMzTmp <- NULL
        allResults$varIonIntTmp <- NULL
        
        write.csv(
            allResults,
            paste0(
                dirname(ms2Files[1]),
                '/allResults_',
                pepTmp,
                '_',
                todaysDate,
                '.csv'
            ),
            row.names = FALSE
        )
        
    } # end function