#' filter samples with low QC and features with large missing values
#' Removes adducts that have not been integrated with many missing values
#' and provides QC on samples
#' @param adductTable character a full path to the peaktable with number of 
#' rows equal to the number of adducts from outputPeakTable()
#' which starts with adductQuantif_peakList_
#' @param percMissing numeric percentage threshold to remove adducts with 
#' missing values. Default
#' is 51. It is recommended to use just over the number of samples in the 
#' smallest
#' group of your study. 51 is used as default for a 50:50 case control study
#' @param HKPmass numeric mass for the housekeeping peptide. Must be the same  
#' asthat in the adduct table. max 2 decimal places.
#' default= 575.3 for the LVNEVTEFAK peptide
#' @param quantPeptideMass numeric mass for the peptide for which adducts are 
#' being quantified, 
#  must be the same as that in the adduct table, max 2 decimal places. 
#' Default is
#' 811.7 for the ALVLIAFAQYLQQCPFEDHVK peptide
#' @param remHKPzero logical if TRUE removes all samples where the 
#' housekeeping peptide is 0. default= FALSE
#' @param remQuantPepzero logical if TRUE removes all samples where the 
#' peptide under quantification is 0. default= FALSE
#' @param remHKPlow logical if TRUE removes all samples where the 
#' housekeeping peptide has an area less than 100000. 
#' default= TRUE.This is recommended because this peak should be large. 
#' If the HKP
#' has been mis-identified quantification of all adducts will be affected.
#' @param outputDir character path to results directory
#' output is a csv file with only adducts and samples that passed filter. 
#' Remaining adducts can be quantified manually 
#' however it is recommended to rescale the quantification results and include
#' the quantification method as a covariate in downstream analysis.
#' @examples filterAdductTable(adductTable=paste0(system.file("extdata",
#' package="adductomicsR"),'/example_adductQuantif_peakList.csv'), percMissing
#' =51,HKPmass = "575.3", quantPeptideMass = "811.7",
#' remHKPzero=FALSE,remQuantPepzero = FALSE, remHKPlow = FALSE, outputDir =
#' NULL)
#' @usage filterAdductTable(adductTable = NULL, percMissing = 51, HKPmass = 
#' "575.3", quantPeptideMass = "811.7", remHKPzero = FALSE, remQuantPepzero 
#' =FALSE, remHKPlow = FALSE, outputDir = NULL)
#'@return csv file
#'@export
filterAdductTable <- function(adductTable = NULL, percMissing = 51, 
HKPmass = "575.3", 
quantPeptideMass = "811.7", remHKPzero = FALSE, 
remQuantPepzero = FALSE, remHKPlow = FALSE, 
outputDir = NULL) {  
    if (is.null(adductTable)) {
        stop("adduct file is missing with no default")
    }
    if (is.null(outputDir)) {
        outputDir <- getwd()
    }
    outputDir <- paste0(outputDir, "/")
    data <- read.csv(adductTable, header = TRUE)

    # remove adducts that have not been quantified correctly
    to_rem <- which(is.na(data$massAcc))
    if (length(to_rem) > 0) {
        data = data[-which(is.na(data$massAcc)), ]
    }
    rownames(data) <- paste0("M", round(data$mass, 2), "_RT", round(data$RT,2))

    # only peak areas with unique adduct ID
    data2 <- data[, grep("mzXML", colnames(data))]

    # remove adducts where >51% have missing values
    perc.missing <- apply(data2, 1, function(x) 
    (length(which(x == 0))/ncol(data2)) * 100)
    to_rem <- which(perc.missing >= percMissing)
    if (length(to_rem) > 0) {
        data2 = data2[-to_rem, ]
    }

    # identify HKP peptide
    hkp <- suppressWarnings(as.numeric(data[grep(HKPmass, 
    rownames(data2)), ]))

    # remove samples where housekeeping not detected
    if (remHKPzero == TRUE) {
        to_rem <- which(hkp == 0)
        if (length(to_rem) > 0) {
            data2 = data2[, -to_rem]
        }
    }

    # remove samples where the housekeeping peptide is very low (but not zero)
    if (remHKPlow == TRUE) {
        hkp <- as.numeric(data2[grep(HKPmass, rownames(data2)), ])
        to_rem <- which(hkp > 0 & hkp < 1e+05)
        if (length(to_rem) > 0) {
            data2 = data2[, -to_rem]
        }
    }

    # remove samples where hkp did not pass QC and also the T3 was not detected
    if (remQuantPepzero == TRUE) {
        T3 <- as.numeric(data2[grep(quantPeptideMass, rownames(data2)), ])
        to_rem <- which(T3 == 0)
        if (length(to_rem) > 0) {
            data2 = data2[, -to_rem]
        }
    }

    write.csv(data2, paste0(outputDir, "adductQuantif_peakList_", 
    Sys.Date(), "_filtered.csv"), 
    row.names = TRUE)

}  # end function