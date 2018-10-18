#' modified function from package OrgMassSpecR
#' @param formula list of character strings representing elemental formula
#' @param charge numeric for charge of the element
#' @return dataframe of a spectrum
#' @usage IsotopicDistributionMod(formula = list(), charge = 1)
#' @examples IsotopicDistributionMod(formula=list("CH3CH2OH","H2O"),charge = 1)
#' @export
IsotopicDistributionMod <- function(formula = list(), charge = 1) {
    if (charge == 0) 
        stop("a charge of zero is not allowed")
        inputFormula <- list(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0, Br = 0, 
        Cl = 0, 
        F = 0, Si = 0)
        inputFormula[names(formula)] <- formula
        simulation <- function(inputFormula) {
            massCarbon <- sum(sample(c(12, 13.0033548378), size = 
            inputFormula$C, 
            replace = TRUE, 
            prob = c(0.9893, 0.0107)))
            massHydrogen <- sum(sample(c(1.0078250321, 2.014101778), 
            size = inputFormula$H, 
            replace = TRUE, prob = c(0.999885, 0.000115)))
            massNitrogen <- sum(sample(c(14.0030740052, 15.0001088984), 
            size = inputFormula$N, 
            replace = TRUE, prob = c(0.99632, 0.00368)))
            massOxygen <- sum(sample(c(15.9949146221, 16.9991315, 17.9991604), 
            size = inputFormula$O, 
            replace = TRUE, prob = c(0.99757, 0.00038, 0.00205)))
            massSulfer <- sum(sample(c(31.97207069, 32.9714585, 
            33.96786683, 35.96708088), 
            size = inputFormula$S, replace = TRUE, prob = c(0.9493, 0.0076,
            0.0429, 2e-04)))
            massPhosphorus <- inputFormula$P * 30.97376151
            massBromine <- sum(sample(c(78.9183376, 80.916291), 
            size = inputFormula$Br, 
            replace = TRUE, prob = c(0.5069, 0.4931)))
            massChlorine <- sum(sample(c(34.96885271, 36.9659026), 
            size = inputFormula$Cl, 
            replace = TRUE, prob = c(0.7578, 0.2422)))
            massFluorine <- inputFormula$F * 18.9984032
            massSilicon <- sum(sample(c(27.9769265327, 28.97649472, 
            29.97377022), 
            size = inputFormula$Si, 
            replace = TRUE, prob = c(0.922297, 0.046832, 0.030872)))
            massMolecule <- sum(massCarbon, massHydrogen, massNitrogen, 
            massOxygen, massSulfer, 
            massPhosphorus, massBromine, massChlorine, 
            massFluorine, massSilicon)
            mz <- massMolecule/abs(charge)
            return(mz)
        }
        sim <- replicate(10000, expr = simulation(inputFormula))
        b <- seq(from = min(sim) - (1/(2 * abs(charge))), to = max(sim) + 1,
        by = 1/abs(charge))
        bins <- cut(sim, breaks = b)
        mz <- round(tapply(sim, bins, mean), digits = 5)
        intensity <- as.vector(table(bins))
        spectrum <- data.frame(mz, intensity)
        spectrum <- spectrum[spectrum$intensity != 0, ]
        spectrum$percent <- with(spectrum, round(intensity/max(intensity) * 
        100,
        digits = 2))
        row.names(spectrum) <- seq_len((nrow(spectrum)))
        return(spectrum)
    }