#' identify peaks
#' @description identifies peaks in a vector of intensities.
#' @param x numeric vector of intensities.
#' @param m number of peaks to identify
#' @usage findPeaks(x, m = 3)
#' @return string of peaks
findPeaks <- function(x, m = 3) {
    diffts <- diff(x, na.pad = FALSE)
    shape <- diff(sign(diffts))
    pks <- vapply(which(shape < 0), FUN = function(i) {
        z <- i - m + 1
        z <- ifelse(z > 0, z, 1)
        w <- i + m + 1
        w <- ifelse(w < length(x), w, length(x))
        if (all(x[c(z:i, (i + 2):w)] <= x[i + 1])) 
            return(i + 1) else return(numeric(0))
        }, FUN.VALUE=numeric(1))
        pks <- unlist(pks)
        pks
    }  # end function