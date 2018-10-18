#' modified \code{\link[OrgMassSpecR]{Digest}} function 
#' (from OrgMassSpecR package)
#' @description allows maxCharge to be set to calculate precursor m/z
#' @param sequence a character string representing the amino acid sequence.
#' @param enzyme is the enzyme to perform in silico digestion with
#' @param missed the maximum number of missed cleavages. 
#' Must be an integer of 0 (default)
#' or greater. An error will result if the specified number of 
#' missed cleavages is
#' greater than the maximum possible number of missed cleavages.
#' @param maxCharge numeric max charge charge for predicted precursor m/z 
#' @param IAA logical. TRUE specifies iodoacetylated cysteine and FALSE 
#'specifies unmodified
#' cysteine. Used only in determining the elemental formula, 
#' not the three letter codes.
#' @param custom list of custom masses
#' @param N15 logical indicating if the nitrogen-15 isotope should be used
#' in place of the default
#' nitrogen-14 isotope.
#' calculation 
#' @examples digestMod('MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIA',
#' enzyme = "trypsin", missed = 0, maxCharge = 8,IAA = TRUE, N15 = FALSE, 
#' custom = list()) 
#' @details see \code{\link{Digest}} for details of further function arguments.
#' @return dataframe
#' @usage digestMod(sequence, enzyme = "trypsin", missed = 0, 
#' maxCharge = 8,IAA = TRUE, N15 = FALSE, custom = list())
#' @export
digestMod <- function(sequence, enzyme = "trypsin", missed = 0, maxCharge = 8,
IAA = TRUE, N15 = FALSE, custom = list()) {
    dsAAnos <- strsplit(gsub("[A-Z]|ds", "", sequence), "\\[|\\]")[[1]]
    dsAAnos <- dsAAnos[dsAAnos != ""]
    dsRes <- as.numeric(unique(strsplit(paste0(dsAAnos, 
collapse = "-"), "-")[[1]]))
    sequence <- gsub("\\[ds|[0-9]|\\]|-", "", sequence)
    seq_vector <- strsplit(sequence, split = "")[[1]]
    end_position <- length(seq_vector)
    if (enzyme == "trypsin") {
        if (seq_vector[end_position] == "K" | seq_vector[end_position] == 
            "R") {
            seq_vector[end_position] <- "!"
            seq_string <- paste(seq_vector, collapse = "")
        } else seq_string <- sequence
        seq_string <- gsub("KP", "!P", seq_string)
        seq_string <- gsub("RP", "!P", seq_string)
        seq_vector <- strsplit(seq_string, split = "")[[1]]
        stop <- grep("K|R", seq_vector)
        start <- stop + 1
    }
    if (enzyme == "trypsin.strict") {
        if (seq_vector[end_position] == "K" | seq_vector[end_position] ==
            "R") {
            seq_vector[end_position] <- "!"
            seq_string <- paste(seq_vector, collapse = "")
        } else seq_string <- sequence
        seq_vector <- strsplit(seq_string, split = "")[[1]]
        stop <- grep("K|R", seq_vector)
        start <- stop + 1
    }
    if (enzyme == "pepsin") {
        if (seq_vector[end_position]=="F" | seq_vector[end_position]=="L" |
            seq_vector[end_position] == 
            "W" | seq_vector[end_position]=="Y" | seq_vector[end_position]==
            "A" |
            seq_vector[end_position] == "E" | seq_vector[end_position] ==
            "Q") {
                seq_vector[end_position] <- "!"
            }
            stop <- grep("F|L|W|Y|A|E|Q", seq_vector)
            start <- stop + 1
        }
        if (enzyme != "trypsin" & enzyme != "trypsin.strict" & 
            enzyme != "pepsin") 
            stop("undefined enzyme, defined enzymes are trypsin, 
        trypsin.strict, and pepsin")
            if (length(stop) == 0) 
                warning("sequence does not contain cleavage sites")
                if (missed > length(stop)) 
                    stop("number of specified missed cleavages is greater 
                than the maximum possible")
                    cleave <- function(sequence, start, stop, misses) {
                        peptide <- substring(sequence, start, stop)
                        mc <- rep(misses, times = length(peptide))
                        result <- data.frame(peptide, start, stop, mc, 
                    stringsAsFactors = FALSE)
                        return(result)
                    }
                    start <- c(1, start)
                    stop <- c(stop, end_position)
                    results <- cleave(sequence, start, stop, 0)
                    if (missed > 0) {
                        for (i in seq_len(missed)) {
                            start_tmp <- start[seq_len((length(start) - i))]
                            stop_tmp <- stop[(1 + i):length(stop)]
                            peptide <- cleave(sequence, start_tmp, stop_tmp, i)
                            results <- rbind(results, peptide)
                        }
                    }
                    # add in whole protein sequence
                    wholeProt <- data.frame(peptide = sequence, start = 1, 
                stop = end_position, mc = "whole sequence")
                    results <- rbind(results, wholeProt)
                    C <- 12
                    H <- 1.0078250321
                    O <- 15.9949146221
                    S <- 31.97207069
                    N <- ifelse(N15 == TRUE, 15.0001088984, 14.0030740052)
                    proton <- 1.007276466
                    residueMass <- function(residue) {
                        if (residue == "A") 
                            mass = C * 3 + H * 5 + N + O
                        if (residue == "R") 
                            mass = C * 6 + H * 12 + N * 4 + O
                        if (residue == "N") 
                            mass = C * 4 + H * 6 + N * 2 + O * 2
                        if (residue == "D") 
                            mass = C * 4 + H * 5 + N + O * 3
                        if (residue == "E") 
                            mass = C * 5 + H * 7 + N + O * 3
                        if (residue == "Q") 
                            mass = C * 5 + H * 8 + N * 2 + O * 2
                        if (residue == "G") 
                            mass = C * 2 + H * 3 + N + O
                        if (residue == "H") 
                            mass = C * 6 + H * 7 + N * 3 + O
                        if (residue == "I") 
                            mass = C * 6 + H * 11 + N + O
                        if (residue == "L") 
                            mass = C * 6 + H * 11 + N + O
                        if (residue == "K") 
                            mass = C * 6 + H * 12 + N * 2 + O
                        if (residue == "M") 
                            mass = C * 5 + H * 9 + N + O + S
                        if (residue == "F") 
                            mass = C * 9 + H * 9 + N + O
                        if (residue == "P") 
                            mass = C * 5 + H * 7 + N + O
                        if (residue == "S") 
                            mass = C * 3 + H * 5 + N + O * 2
                        if (residue == "T") 
                            mass = C * 4 + H * 7 + N + O * 2
                        if (residue == "W") 
                            mass = C * 11 + H * 10 + N * 2 + O
                        if (residue == "Y") 
                            mass = C * 9 + H * 9 + N + O * 2
                        if (residue == "V") 
                            mass = C * 5 + H * 9 + N + O
                        if (residue == "C" & IAA == FALSE) 
                            mass = C * 3 + H * 5 + N + O + S
                        if (residue == "C" & IAA == TRUE) 
                            mass <- ifelse(N15 == FALSE, C * 5 + H * 8 + N *
                            2 + O * 2 + S, C * 5 + H * 8 + N + 14.0030740052 +
                            O * 2 + S)
                        if (length(custom) != 0) 
                            for (i in seq_len(
                            length(custom$code))) if (residue==custom$code[i]) 
                            mass = custom$mass[i]
                            return(mass)
                        }
                        mz <- vector("list", length = nrow(results))
                        for (i in seq_len(nrow(results))) {
                            peptide_vector <- strsplit(results$peptide[i],
                            split = "")[[1]]
                            peptide_mass <- sum(vapply(peptide_vector, 
                            residueMass, 
                            FUN.VALUE=numeric(1)))
                            mz[[i]] <- round((peptide_mass + H * 2 + O + 
                            (c(seq_len(maxCharge)) * 
                            proton))/c(seq_len(maxCharge)), 
                        digits = 3)
                        }
                        mz <- as.data.frame(do.call("rbind", mz))
                        names(mz) <- paste0("mz", seq_len(maxCharge))
                        results <- cbind(results, mz)
                        return(results)
                    }
