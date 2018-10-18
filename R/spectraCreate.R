#' Deconvolute both MS2 and MS1 levels scans adductomics
#'
#' @param MS2file character vector of mzXML file locations
#' @param TICfilter numeric minimimum total ion current of an MS/MS scan. Any
#' MS/MS scan below this value will be filtered out (default=0).
#' @param DNF dynamic noise filter minimum signal to noise threshold 
#' (default = 2), calculated as the ratio between the linear model predicted 
#' intensity value and the actual intensity.
#' @param minInt numeric minimum intensity value
#' @param minPeaks minimum number of signal peaks following dynamic 
#' noise filtration (default = 5).
#' @examples spectraCreate(MS2file = paste0(system.file("extdata", package = 
#' "adductData"),'/ORB35017.mzXML'), TICfilter = 10000, DNF= 2, minInt =
#' 100, minPeaks = 5)
#' @usage spectraCreate(MS2file = NULL, TICfilter = 10000, DNF = 2, minInt = 
#' 100, minPeaks = 5)
#' @return list of MS2 spectra
#'@export
spectraCreate <- function(MS2file = NULL, TICfilter = 10000, 
DNF = 2, minInt = 100, 
minPeaks = 5) {
    # MS2 file name
    MS2fileName <- basename(MS2file)
    message(paste0("Reading ", MS2fileName, "...\n"))
    flush.console()
    # read MS2 file
    MS2file <- mzR::openMSfile(MS2file)
    message("...DONE")
    flush.console()
    message("extracting metaData from MS2 file...\n")
    flush.console()
    # extract info from MS/MS file
    fileInfoDf <- mzR::header(MS2file)
    message("...DONE\n")
    flush.console()
    # remove all MS/MS scans where the TIC is less than 
    #the minimum TIC threshold set
    # by the user
    message(paste0("Of a total of ", 
    length(which(fileInfoDf$msLevel == 2)), " MS2 spectra...\n"))
    flush.console()
    # index ms2 scan and above TIC filter
    fileInfoDf$MS2TICfilt <- (fileInfoDf$msLevel == 2 &
    fileInfoDf$totIonCurrent > 
    TICfilter) * 1
    MS2TICfilt.indx <- which(fileInfoDf$MS2TICfilt == 1)
    message(paste0(length(MS2TICfilt.indx),
    " MS2 spectra were above the TIC filter of ", 
    TICfilter))
    flush.console()
    if (length(MS2TICfilt.indx) == 0) {
        # cond if no scan above the TIC filter
        warning(paste0("No MS2 levels scans above TIC filter of ",
        TICfilter, " in ", 
        MS2fileName, ",  reduce the TIC filter parameter or check that
        the file has been converted 
        to the mzXML format correctly."), 
        immediate. = TRUE)
        flush.console()
        message("...moving to next MS2 file")
        flush.console()
    } else {

        # applying loess variable baseline noise removal
        message("applying savitsky-golay baseline subtraction...\n")
        flush.console()

        savGolRemSpectra <- lapply(MS2TICfilt.indx, function(scan) {
            specDfTmp <- mzR::peaks(MS2file, scan)
            savGolTmp <- pracma::savgol(specDfTmp[, 2], 51)
            specDfTmp[, 2] <- specDfTmp[, 2] - savGolTmp
            specDfTmp <- specDfTmp[specDfTmp[, 2] > 0, , drop = FALSE]
        })

            if (DNF > 0) {
                message("Applying dynamic noise filter...\n")
                flush.console()
            }
            NoiseFiltSpectra <- lapply(savGolRemSpectra, function(spec) {
                dynamicNoiseFilter(spec, DNF, minPeaks, minInt)
            })
                # id all above min/ noise
                aboveMinPeaksIndx <- vapply(NoiseFiltSpectra,
                function(spec) spec$aboveMinPeaks, FUN.VALUE=logical(1))

                if (all(aboveMinPeaksIndx == FALSE)) {
                    warning(" No composite spectra were found 
                    above the DNF filter of ", 
                    DNF, " or any peaks above the minimum 
                    intensity (minInt) ", minInt, 
                    " or minimum peaks setting of ", minPeaks, " in ", 
                    MS2fileName,
                    " consider reducing one of more of these 
                    parameter settings and/or", 
                    " check that the file has been converted 
                    to the mzXML format correctly.", 
                    immediate. = TRUE)
                    flush.console()
                    message("...moving to next MS2 file")
                    flush.console()
                } else {
                    # unlist and rbind metaData and add to file.info dataframe
                    NoiseFiltSpectra <- unlist(NoiseFiltSpectra, recursive=
                    FALSE)
                    metaData <- do.call("rbind", 
                    NoiseFiltSpectra[which(names(NoiseFiltSpectra) == 
                "metaData")])
                    metaData <- data.frame(metaData, 
                    aboveMinPeaks = aboveMinPeaksIndx * 1)
                    fileInfoDf[, colnames(metaData)] <- 0
                    fileInfoDf[MS2TICfilt.indx, colnames(metaData)] <- 
                    metaData

                    message("...done")
                    flush.console()

                    # spectra
                    specNames <- fileInfoDf$acquisitionNum[
                    fileInfoDf$aboveMinPeaks 
                    == 1]
                    MS2spectra <- NoiseFiltSpectra[which(names(
                    NoiseFiltSpectra) ==
                    "Above.noise")][aboveMinPeaksIndx]
                    names(MS2spectra) <- specNames
                    return(list(metaData = fileInfoDf, MS2spectra = 
                    MS2spectra))
                    close(MS2file)
                }  # cond if no spectra are above noise filter
            }  # cond if no scan above the TIC filter
        }  # end func