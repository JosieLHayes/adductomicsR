#' Constructor of adductSpec object deconvolute spectra 
#' MS2 and MS1 levels
#' @param mzXmlDir character a full path to a directory containing 
#' either .mzXML or .mzML data
#' @param runOrder character a full path to a csv file specifying the 
#' runorder for each of the files
#' the first column must contain the precise file name and the second 
#' column an integer representing the precise run order.
#' @param nCores numeric the number of cores to use for parallel computation. 
#' The default
#' is to use all available cores detected using the function 
#' parallel::detectCores()
#' @param intStdMass numeric vector of the mass of the internal standard. 
#' Default is the mass of 
#  the T3 peptide of HSA.
#' @param intStdPeakList numeric vector of masses for the internal 
#' standard peaks
#' @param TICfilter numeric minimimum total ion current of an MS/MS scan.
#' Any MS/MS scan below this value will be filtered out (default=0).
#' @param DNF dynamic noise filter minimum signal to noise threshold 
#' (default = 2), calculated as the ratio between the linear model predicted 
#' intensity value and the actual intensity.
#' @param minInt numeric minimum intensity value
#' @param minPeaks minimum number of signal peaks following dynamic 
#' noise filtration (default = 5).
#' @param intStd_MaxMedRtDrift numeric the maximum retention time drift 
#' window (in seconds) to identify internal standard MS/MS spectrum scans
#' (default = 600).
#' @param intStd_MaxPpmDev numeric the maximum mass accuracy window (in ppm).
#' to identify internal standard MS/MS spectrum scans (default = 200 ppm).
#' @param minSpecEx numeric the minimum percentage of the total ion 
#' current explained
#' by the internal standard fragments (default = 40). Sometime spectra are not
#' identified due to this cutoff being set too high. If unexpected datapoints
#' have been interpolated then reduce this value.
#' @param outputPlotDir character string for the output directory for plots, 
#' default is working directory.
#' @description reads mzXML files from a directory extracts metadata info,
#' groups ion signals with \code{\link{signalGrouping}}, filters noise 
#' dynamically
#' \code{\link{dynamicNoiseFilter}} and identifies precursor ion charge state,
#' by isotopic pattern.
#' @examples 
#' adductSpecGen(mzXmlDir=system.file("extdata", package = 
#' "adductData"), runOrder=paste0(system.file("extdata", package = 
#' "adductomicsR"),'/runOrder.csv'), nCores=4, 
#' intStdMass=834.77692,intStdPeakList=c(290.21, 403.30, 516.38, 
#' 587.42,849.40, 884.92, 958.46, 993.97,1050.52, 1107.06, 1209.73, 
#' 1337.79,1465.85), TICfilter=10000, DNF=2, minInt=300, 
#' minPeaks=5,intStd_MaxMedRtDrift=360, 
#' intStd_MaxPpmDev=200,minSpecEx=40,outputPlotDir=NULL)
#' @usage adductSpecGen(mzXmlDir=NULL, runOrder=NULL, nCores=NULL,
#' intStdMass=834.77692,intStdPeakList=c(290.21, 403.30, 516.38,
#' 587.42,849.40, 884.92, 958.46, 993.97,1050.52, 1107.06, 1209.73,
#' 1337.79,1465.85),TICfilter=10000, DNF=2, minInt=300,
#' minPeaks=5,intStd_MaxMedRtDrift=360, intStd_MaxPpmDev=200,minSpecEx=40,
#' outputPlotDir=NULL)
#' @return adductSpec object
#' @export
adductSpecGen <- function(mzXmlDir=NULL, runOrder=NULL, nCores=NULL,
        intStdMass=834.77692,intStdPeakList=c(290.21, 403.30, 516.38,
        587.42,849.40, 884.92, 958.46, 993.97,1050.52, 1107.06, 1209.73,
        1337.79,1465.85),TICfilter=10000, DNF=2, minInt=300,
        minPeaks=5,intStd_MaxMedRtDrift=360, intStd_MaxPpmDev=200,minSpecEx=40,
        outputPlotDir=NULL){
# set global options
options(stringsAsFactors = FALSE)
# if is.null mzXmlDir select mzXML file containing directory
if(is.null(mzXmlDir)){
message("Select your .mzXML data directory")
flush.console()
mzXmlDir <- tcltk::tk_choose.dir(default = "",
caption = "Select your .mzXML data directory")
}
# if is.null runOrder select .csv file
if(is.null(runOrder)){
message("Select your .csv file containing the precise run order")
flush.console()
runOrder <- tcltk::tclvalue(tcltk::tkgetOpenFile(title="select your run   
order csv file"))
}  
if(is.character(runOrder)){
if(!grepl('\\.csv$', runOrder)){
stop('runOrder table file must be a .csv')
}
runOrder <- as.data.frame(data.table::fread(runOrder, header=TRUE,       
    stringsAsFactors = FALSE))
}
if(!is.data.frame(runOrder)){
stop("runOrder table must be supplied as either a full path to a csv file 
or as a data.frame")
}
if(!is.character(runOrder[, 1])){
stop('first column of runorder table must contain file name (is           
character)')
}
if(!is.integer(runOrder[, 2])){
stop('second column of runOrder table must contain a run order integer (is 
integer)')
}
    # identify all mzXML files in raw-data directory
    MS2files <- list.files(path=mzXmlDir, pattern = "*.mzXML$",
    full.names=TRUE,
    recursive=FALSE)
    message(paste0(length(MS2files),
    " MS (.mzXML) files were detected within the directory..."))
    flush.console()
    # sort ms2 files by run order info
    idxTmp <- match(gsub('\\.mzXML$', '', runOrder[, 1]), gsub('\\.mzXML$',
    '',      basename(MS2files)))
    if(any(is.na(idxTmp))){
    warning('files in runOrder table not detected in mzXMLDir:\n',            
    paste0(runOrder[is.na(idxTmp), 1], collapse='\n'), '\n', immediate. = TRUE)
    runOrder <- runOrder[!is.na(idxTmp), , drop=FALSE]
    idxTmp <- idxTmp[!is.na(idxTmp)]
    }
    MS2files <- MS2files[idxTmp]
    MS2files <- MS2files[order(runOrder[, 2])]
    # directory character lengths
    dirCharLengths <- vapply(dirname(MS2files), nchar, FUN.VALUE=numeric(1))
    # error handling level directory structure
    if(length(unique(dirCharLengths)) != 1){
    stop("directories containing .mzXML files are not level...e.g.\n",
        names(dirCharLengths)[which.min(dirCharLengths)], "\n",
        names(dirCharLengths)[which.max(dirCharLengths)], "\n")
    } 
    if(!is.null(nCores)){
        # start time to record computation time
        pmt <- proc.time()
        message(paste0("Starting SNOW cluster with ", nCores, " local          
        sockets..."))
        flush.console()
        cl <- parallel::makeCluster(nCores, outfile='')
        doSNOW::registerDoSNOW(cl)
        message("Performing Savitsky-Golay baseline subtraction, dynamic noise 
        filtration and adductSpec object construction...\n")
        flush.console()
        progress <- function(n) cat(paste0(n, ' of ', length(MS2files),
        ' complete (', basename(MS2files)[n],
        ').\n'))
        opts <- list(progress=progress)
        massDriftFiles <- gsub('\\.mzXML$|\\.mzML$', '.massDrift.csv',
        MS2files)
        # foreach and dopar from foreach package
        Results <- foreach(i=seq_len(length(MS2files)),  .packages=c('mzR',    
        "RcppEigen", 'pracma'),.options.snow=opts) %dopar% {
        Res.tmp <- spectraCreate(MS2file=MS2files[i],
        TICfilter=TICfilter, DNF=DNF, minInt=minInt,
        minPeaks=minPeaks)
        # if present add mass drift files
        if(file.exists(massDriftFiles[i])){
        massDriftTmp <- read.csv(massDriftFiles[i], header = TRUE,             
        stringsAsFactors = FALSE)
        ms1Level <- Res.tmp[[1]]$msLevel == 1
        ms1MassDrift <- massDriftTmp[ms1Level, , drop=FALSE]
        ms1MassDrift$ppmDrift[ms1MassDrift[, 3] == 0] <- NA
        Res.tmp[[1]]$ppmDrift <- NA
        Res.tmp[[1]]$ppmDrift[ms1Level] <-                                     
        zoo::na.spline(ms1MassDrift$ppmDrift)
        rleLevel <- rle(Res.tmp[[1]]$msLevel == 2)
        ms2Level <- which(ms1Level == FALSE)
        ms1Level <- which(ms1Level)
        precMassesIdx <- ms1Level[{ms1Level + 1} %in% ms2Level]
        Res.tmp[[1]]$adjPrecursorMZ <- 0
        ppmVals <- mapply(rep, Res.tmp[[1]]$ppmDrift[precMassesIdx],           
        each=rleLevel$lengths[rleLevel$values == TRUE], SIMPLIFY = FALSE)
        ppmVals <- do.call(c, ppmVals)
        Res.tmp[[1]]$adjPrecursorMZ[ms2Level] <-                               
        Res.tmp[[1]]$precursorMZ[ms2Level] -                                   
        {{Res.tmp[[1]]$precursorMZ[ms2Level]/1E06} * ppmVals}
        Res.tmp[[1]]$ppmDrift[ms2Level] <- ppmVals
        massDriftTmp$precursorMZ <- Res.tmp[[1]]$precursorMZ
        massDriftTmp$adjPrecursorMZ <- Res.tmp[[1]]$adjPrecursorMZ
        massDriftTmp$ppmDrift <- Res.tmp[[1]]$ppmDrift
        write.csv(massDriftTmp, massDriftFiles[i], row.names = FALSE)
        }
        return(Res.tmp)
        }
        # stop SNOW cluster
        parallel::stopCluster(cl)
        proc.time() - pmt
    } else {
    # create list to store results
        Results <- vector("list", length(MS2files))
        massDriftFiles <- gsub('\\.mzXML$|\\.mzML$', '.massDrift.csv',         
        MS2files)
        for(i in seq_len(length(MS2files))){
        Res.tmp <- spectraCreate(MS2file=MS2files[i],
                                TICfilter=TICfilter, DNF=DNF, minInt=minInt,
                                minPeaks=minPeaks)
        # if present add mass drift files
        if(file.exists(massDriftFiles[i])){
        massDriftTmp <- read.csv(massDriftFiles[i], header = TRUE,             
            stringsAsFactors = FALSE)
        ms1Level <- Res.tmp[[1]]$msLevel == 1
        ms1MassDrift <- massDriftTmp[ms1Level, , drop=FALSE]
        ms1MassDrift$ppmDrift[ms1MassDrift[, 3] == 0] <- NA
        Res.tmp[[1]]$ppmDrift <- NA
        Res.tmp[[1]]$ppmDrift[ms1Level] <-                                     
        zoo::na.spline(ms1MassDrift$ppmDrift)
        rleLevel <- rle(Res.tmp[[1]]$msLevel == 2)
        ms2Level <- which(ms1Level == FALSE)
        ms1Level <- which(ms1Level)
        precMassesIdx <- ms1Level[{ms1Level + 1} %in% ms2Level]
        Res.tmp[[1]]$adjPrecursorMZ <- 0
        ppmVals <- mapply(rep, Res.tmp[[1]]$ppmDrift[precMassesIdx],           
            each=rleLevel$lengths[rleLevel$values == TRUE], SIMPLIFY = FALSE)
        ppmVals <- do.call(c, ppmVals)
        Res.tmp[[1]]$adjPrecursorMZ[ms2Level] <-                               
        Res.tmp[[1]]$precursorMZ[ms2Level] -                                   
        {{Res.tmp[[1]]$precursorMZ[ms2Level]/1E06} * ppmVals}
        Res.tmp[[1]]$ppmDrift[ms2Level] <- ppmVals
        massDriftTmp$precursorMZ <- Res.tmp[[1]]$precursorMZ
        massDriftTmp$adjPrecursorMZ <- Res.tmp[[1]]$adjPrecursorMZ
        massDriftTmp$ppmDrift <- Res.tmp[[1]]$ppmDrift
        write.csv(massDriftTmp, massDriftFiles[i], row.names = FALSE)
    }
    Results[[i]] <- Res.tmp
    }
}
    # add MS2 file names to results list
    names(Results) <- basename(MS2files)
    Results <- unlist(Results, recursive=FALSE)
    adductSpecTmp <- new("adductSpec")
    # add in MS2 spectra
    adductSpecTmp@adductMS2spec <- Results[grep('MS2spectra', names(Results))]
    # add MS1 isotopes
    #adductSpecTmp@adductMS1iso <- Results[grep('MS1isotopes', names(Results))]
    # rbind meta data and add in file names
    metaDataTmp <- Results[grep('metaData', names(Results))]
    fileNamesTmp <- rep(gsub('\\.metaData$', '', names(metaDataTmp)),
        vapply(metaDataTmp, nrow, numeric(1)))
    metaDataTmp <- do.call(rbind, metaDataTmp)
    metaDataTmp <- data.frame(mzXMLFile=fileNamesTmp, metaDataTmp, 
    stringsAsFactors=FALSE)
    adductSpecTmp@metaData <- metaDataTmp
    adductSpecTmp@file.paths <- MS2files
    adductSpecTmp@Parameters <- data.frame(nCores, ifelse(is.null(intStdMass), 
    NA, intStdMass), TICfilter, DNF, minInt, minPeaks, intStd_MaxMedRtDrift,  
    intStd_MaxPpmDev, stringsAsFactors=FALSE)
    # if necessary id internal standard
    if(!is.null(intStdMass)){
    message('Identifying internal standard MS/MS spectra within a m/z window  
    of ', intStd_MaxPpmDev, ' ppm and a retention time window of ',           
    round(intStd_MaxMedRtDrift/60, 2), ' minutes...\n')
    flush.console()
    intStdScans <- peakListId(adductSpecTmp, peakList=intStdPeakList,         
    minPeaksId=7,
    minSpecEx=minSpecEx, outputPlotDir=NULL,
    exPeakMass=intStdMass, maxRtDrift=intStd_MaxMedRtDrift,
    maxPpmDev = intStd_MaxPpmDev, closestMassByFile=FALSE)
    # remove duplicates
    intStdScans <- intStdScans[duplicated(intStdScans$mzXMLFile) == FALSE, ,  
    drop=FALSE]
    # take average ppmDev of all less than 5 ppm
    # adjust theoretical then calculate deviation from observed intstd mass.
    ppmDiffFile <- ((intStdScans$precursorMz - intStdMass)/ intStdMass) * 1E06
    obsStdMass <- mean(intStdScans$precursorMz[abs(ppmDiffFile) <= 5])
    ppmDiffFile <- ((intStdScans$precursorMz - obsStdMass)/ obsStdMass) * 1E06
    medRtDiffFile <- intStdScans$retentionTime -                    
    median(intStdScans$retentionTime)
    intStdScans$IntStd_ppmDev <- ppmDiffFile
    intStdScans$IntStd_medRtDev <- medRtDiffFile / 60
    # match names with file path names
    runOrderTmp <- match(intStdScans$mzXMLFile, 
    basename(adductSpecTmp@file.paths))
    # sort table
    intStdScans <- intStdScans[order(runOrderTmp), , drop=FALSE]
    nFilesTmp <- seq_len(length(adductSpecTmp@file.paths))
    missingFiles <- setdiff(nFilesTmp, runOrderTmp)
    # if neccessary cubic spline interpolate missing files with zoo::na.spline
    if(length(missingFiles) > 0){
    message('cubic spline interpolating ', length(missingFiles), ' missing 
    values...\n')
    flush.console()
    # add NA for missing files
    ppmDiffFile <- c(ppmDiffFile, rep(NA, length(missingFiles)))
    medRtDiffFile <- c(medRtDiffFile, rep(NA, length(missingFiles)))
    runOrderTmp <- c(runOrderTmp, missingFiles)
    # sort by runorder
    ppmDiffFile <- ppmDiffFile[order(runOrderTmp)]
    medRtDiffFile <- medRtDiffFile[order(runOrderTmp)]
    runOrderTmp <- runOrderTmp[order(runOrderTmp)]
    # interpolate ppmdiff and rtdiff
    ppmDiffFile <- zoo::na.spline(ppmDiffFile)
    medRtDiffFile <- zoo::na.spline(medRtDiffFile)
    }
    # add to meta data
    indxTmp <- match(adductSpecTmp@metaData$mzXMLFile,
    basename(adductSpecTmp@file.paths))
    adductSpecTmp@metaData$intStdRtDrift <- medRtDiffFile[indxTmp]
    adductSpecTmp@metaData$intStdPpmDrift <- ppmDiffFile[indxTmp]
    # save plot
    png(paste0(ifelse(is.null(outputPlotDir), paste0(getwd(), '/'),
    outputPlotDir),'internalStandard_plots.png'), width = 961, height = 587)
    par(mfrow=c(2, 2))
    plot(runOrderTmp, medRtDiffFile / 60, xlab='run order',
        ylab='deviation from median Rt (mins)')
    # highlight inputed values
    if(length(missingFiles) > 0){
    points(missingFiles, medRtDiffFile[missingFiles] / 60, col='red', pch=19)
    }
    plot(runOrderTmp, ppmDiffFile, xlab='run order', ylab='ppm difference int. 
    std.')
    # highlight inputed values
    if(length(missingFiles) > 0){
    points(missingFiles, ppmDiffFile[missingFiles], col='red', pch=19)
    }
    plot(ppmDiffFile, medRtDiffFile / 60, ylab='deviation from median Rt (
        mins)', xlab='ppm difference int. std.', main=paste('medRT dev 
        range:', round(min(medRtDiffFile / 60), 2), '-', 
        round(max(medRtDiffFile/ 60), 2), 'ppm dev range:', 
        round(min(ppmDiffFile), 2), '-', round(max(ppmDiffFile), 2)))
    if(length(missingFiles) > 0){
    points(ppmDiffFile[missingFiles], medRtDiffFile[missingFiles] / 60, 
        col='red', pch=19)
    }
    plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', 
    yaxt = 'n')
    text(x = 0.5, y = 0.5, paste0("red data-points: cubic spline interpolated 
    files (n=",length(missingFiles), ") where"),
    cex = 1.3, col = "red")
    text(x = 0.5, y = 0.3, paste("internal standard MS2 spectrum were not 
    identified"),cex = 1.3, col = "red")
    dev.off()
    # write table
    intStdScans$peakId <- NULL
    intStdScans$peakNo <- NULL
    intStdScans$plotFile <- paste0(ifelse(is.null(outputPlotDir), 
    paste0(getwd(), '/'), outputPlotDir),intStdScans$name, '.png')
    intStdScans$retentionTime <- intStdScans$retentionTime / 60
    write.csv(intStdScans, paste0(ifelse(is.null(outputPlotDir), 
    paste0(getwd(), '/'), outputPlotDir),'internalStd_scansTable.csv'), 
    row.names = FALSE)
    message('...DONE')
    flush.console()
    }
    return(adductSpecTmp)
} # end function
    setMethod("show", "adductSpec", function(object) {
    if(length(object@file.paths) > 0){
    cat("A \"adductSpec\" class object derived from", 
    length(object@file.paths),"MS2 files \n\n")
    cat("Consisting of:\n", sum(object@metaData$msLevel == 1), "raw MS1 
    scans\n", sum(object@metaData$msLevel == 2), "raw MS2 scans of which",
    sum(object@metaData$aboveMinPeaks == TRUE), "scans were above the noise 
    following filtration.\n\n")
    cat("m/z range:", round(object@metaData$lowMZ[1], digits=3), "-",
    round(object@metaData$highMZ[1], digits=3), "\n\n")
    if(length(object@groupMS2spec) > 0){
    cat('\nInter-file grouping:\n')
    cat(length(object@groupMS2spec), 'inter-file composite spectra groups 
    were identified.\n')
    }
    cat("\nAdductome identification:\n")
    # cat()
    memsize <- object.size(object)
    cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
    } else {
    cat("A new empty\"adductSpec\" class object")
    }
})
    # set method concatenate
    setMethod("c", signature(x = "adductSpec"), function(x, ...){
    elements = list(x, ...)
    # error handling check if all adductSpec object
    if(any(vapply(elements, function(ele) is(ele,'adductSpec'),
    FUN.VALUE=logical(1)) == FALSE)){
    stop('all elements must be an adductSpec class object')
    } 
    emptyAdductSpec <- new('adductSpec')
    # bind together results
    # do not include any group info or other information
    for (i in seq_len(length(elements))){
    emptyAdductSpec@adductMS2spec <- c(emptyAdductSpec@adductMS2spec,        
    elements[[i]]@adductMS2spec)
    metaDataTmp <- elements[[i]]@metaData
    #metaDataTmp$aboveMinPeaks <- NULL
    metaDataTmp$msPrecursor_group <- NULL
    metaDataTmp$interMSMSrtGroups <- NULL
    metaDataTmp$predRtLoess <- NULL
    metaDataTmp$MS2groupFreq <- NULL
    metaDataTmp$MS2groupFreqAbove <- NULL
    emptyAdductSpec@metaData <- rbind(emptyAdductSpec@metaData, metaDataTmp)
    emptyAdductSpec@file.paths <- c(emptyAdductSpec@file.paths,          
    elements[[i]]@file.paths)
    }
    message('Grouping, retention time correction and composite spectra 
    identification must be repeated in the concatenated "adductSpec" class 
    object...\n')
    flush.console()
    return(emptyAdductSpec)
}) # end function
