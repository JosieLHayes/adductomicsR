#' MS/MS spectrum grouping and retention time deviation modelling for 
#' adductomicsR
#' @param MS2Dir character a full path to a directory containing either
#' .mzXML or .mzML data
#' @param runOrder character a full path to a csv file specifying the 
#' runorder for each of the files
#' the first column must contain the precise file name and the second 
#' column an integer representing the precise run order.
#' @param nCores numeric the number of cores to use for parallel computation.
#' The default is to use all available cores detected using the function 
#' parallel::detectCores()
#' @param TICfilter numeric minimimum total ion current of an MS/MS scan. 
#' Any MS/MS scan below this value will be filtered out (default=0).
#' @param intStdPeakList character a comma seperated list of expected 
#' fragment ions
#' for the internal standard spectrum (no white space).
#' @param intStdMass numeric expected mass-to-charge ratio of 
#' internal standard
#' precursor (default = 834.77692).
#' @param intStd_MaxMedRtDrift numeric the maximum retention time 
#' drift window (in seconds)
#' to identify internal standard MS/MS spectrum scans (default = 600).
#' @param intStd_MaxPpmDev numeric the maximum mass accuracy window (in ppm)
#' to identify internal standard MS/MS spectrum scans (default = 200 ppm).
#' @param minSpecEx numeric the minimum percentage of the total ion current
#' explained by the internal standard fragments (default = 40). 
#' Sometimes spectra are not
#' identified due to this cutoff being set too high. If unexpected datapoints
#' have been interpolated then reduce this value.
#' @param minDotProd numeric. Minimum mean dot product spectral similarity 
#' score to keep a spectrum within an MS/MS group (default = 0.8).
#' @param percMissing numeric. percentage of missing files 
#' for a MS/MS scan group to be
#' utilized in the loess retention time deviation model. 
#' Roughly 15 percent missing values (default = 15\%) is a good starting point 
#' (e.g. nMissing=10 for 68 samples).
#' @param percExtra numeric percentage of extra scans above the total number of
#' files for a MS/MS scan group to be utilized in the 
#' loess retention time deviation model.
#' If a MS/MS scan group consists of many scans far in 
#' excess of the number of files
#' then potentially MS/MS scans from large tailing peaks or 
#' isobars may be erroneously
#' grouped together and used to adjust retention time incorrectly 
#' (default = 100\% i.e. the peak group can only have one scan per file, 
#' this value can be increased if two or more consecutive scans for 
#' example can be considered).
#' @param smoothingSpan numeric. fixed smoothing span, 
#' argument to \code{\link{loess}}.
#' If argument is not supplied then optimal smoothing span
#' is calculated for each file seperately using 7-fold CV.
#' @param saveRtDev integer (default = 1) should just the retention time
#' deviation model be saved (TRUE = 1) or the adductSpec class object 
#'(FALSE = 0) as .RData workspace files.
#' @examples
#' rtDevModelling(MS2Dir=hubCache(temp),nCores=4,runOrder=paste0(
#' system.file("extdata",package="adductomicsR"),
#' '/runOrder.csv'))
#' @return LOESS RT models as adductSpectra adductSpec object
#' @import adductData
#' @usage rtDevModelling(MS2Dir = NULL, runOrder = NULL, 
#' nCores = parallel::detectCores(), TICfilter = 0, 
#' intStdPeakList=c(290.21, 403.30, 516.38, 587.42,849.40, 884.92, 958.46,
#' 993.97,1050.52, 1107.06, 1209.73, 1337.79,1465.85), 
#' intStdMass = 834.77692, intStd_MaxMedRtDrift = 600, intStd_MaxPpmDev = 200,
#' minSpecEx = 40, 
#' minDotProd = 0.8, percMissing = 15, percExtra = 100, smoothingSpan = 0.8,
#' saveRtDev = 1)
#' @export
rtDevModelling <- function(MS2Dir = NULL, runOrder = NULL, 
nCores = parallel::detectCores(), 
TICfilter = 0, 
intStdPeakList=c(290.21, 403.30, 516.38, 587.42,
849.40, 884.92, 958.46, 993.97,
1050.52, 1107.06, 1209.73, 1337.79,
1465.85), 
intStdMass = 834.77692, intStd_MaxMedRtDrift = 600, intStd_MaxPpmDev = 200,
minSpecEx = 40, 
minDotProd = 0.8, percMissing = 15, percExtra = 100, smoothingSpan = 0.8,
saveRtDev = 1) {

    #intStdPeakList <- as.numeric(strsplit(intStdPeakList, ",")[[1]])
    if (is.null(MS2Dir)) {
        stop("Please provide an .mzXML data directory")   
    }
    fSlashIndx <- grepl("/$", MS2Dir)
    MS2Dir <- ifelse(fSlashIndx == FALSE, paste0(MS2Dir, "/"), MS2Dir)
    if (is.null(runOrder)) {
        stop("Please provide an run order file")  
    }

    adductSpectra <- adductSpecGen(mzXmlDir = MS2Dir, nCores = nCores, 
    runOrder = runOrder, DNF = 0, minInt = 0, TICfilter = TICfilter, 
    intStdPeakList = intStdPeakList, 
    intStd_MaxMedRtDrift = intStd_MaxMedRtDrift, 
    intStd_MaxPpmDev = intStd_MaxPpmDev, 
    outputPlotDir = MS2Dir, intStdMass = intStdMass, 
    minSpecEx = minSpecEx)
    maxRtDrift_intStd <- max(abs(adductSpectra@metaData$intStdRtDrift)) * 1.5
    #adductSpectra@metaData$intStdPpmDrift <- NULL
    #adductSpectra@metaData$intStdRtDrift <- NULL
    adductSpectra <- ms2Group(adductSpectra, nCores, 
    ms1mzError = 0.01, ms2mzError = 1, 
    maxRtDrift = maxRtDrift_intStd, dotProdClust = TRUE, 
    compSpecGen = FALSE, 
    minDotProd = minDotProd)
    nExtra <- round(length(adductSpectra@file.paths) * {
        percExtra/100
        }, 0)
        nMissing <- round(length(adductSpectra@file.paths) * {
            percMissing/100
        }, 0)
            # # retention time correction
            adductSpectra <- retentionCorr(adductSpectra, nMissing = nMissing, 
            nExtra = nExtra, 
        smoothingSpan = smoothingSpan, outputFileDir = MS2Dir)
            if (saveRtDev == 1) {
                cat(paste0("saving retention-time deviation 
                model .RData file:\n",
                MS2Dir, 
                "rtDevModels.RData\n\nThis file is neccessary for the adduct 
            identification and adduct quantification workflows."))
                rtDevModelSave(adductSpectra, outputDir = MS2Dir)
            } else {
                cat(paste0("saving adductSpec class object as an
                .RData file:\n", 
                MS2Dir, 
                "adductSpectra.RData\n\n"))
                save(adductSpectra, file = paste0(MS2Dir,
                "adductSpectra.RData"))
            }
        }  # end function