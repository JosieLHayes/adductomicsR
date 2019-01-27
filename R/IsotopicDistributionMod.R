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
    inputFormula <-
        list(
            C = 0,
            H = 0,
            N = 0,
            O = 0,
            S = 0,
            P = 0,
            Br = 0,
            Cl = 0,
            F = 0,
            Si = 0
        )
    inputFormula[names(formula)] <- formula
    simulation <- function(inputFormula) {
        data = read.table(paste0(system.file("extdata", 
        package = "adductomicsR"), '/ElementMasses.txt'), h=TRUE)
        row.names(data)= data[,1]
        data = data[,-1]
        massCarbon <- sum(sample(
            as.numeric(data['C', 1:2]),
            size = inputFormula$C,
            replace = TRUE,
            prob = as.numeric(data['C', 5:6])
        ))
        massHydrogen <-
            sum(sample(
                as.numeric(data['H', 1:2]),
                size = inputFormula$H,
                replace = TRUE,
                prob = as.numeric(data['H', 5:6])
            ))
        massNitrogen <-
            sum(sample(
                as.numeric(data['N', 1:2]),
                size = inputFormula$N,
                replace = TRUE,
                prob = as.numeric(data['N', 5:6])
            ))
        massOxygen <-
            sum(sample(
                as.numeric(data['O', 1:3]),
                size = inputFormula$O,
                replace = TRUE,
                prob = as.numeric(data['O', 5:7])
            ))
        massSulfer <- sum(sample(
            as.numeric(data['S', 1:4]),
            size = inputFormula$S,
            replace = TRUE,
            prob = as.numeric(data['S', 5:8])
        ))
        massPhosphorus <- inputFormula$P * as.numeric(data['P', 1])
        massBromine <- sum(sample(
            as.numeric(data['Br', 1:2]),
            size = inputFormula$Br,
            replace = TRUE,
            prob = as.numeric(data['Br', 5:6])
        ))
        massChlorine <- sum(sample(
            as.numeric(data['Cl', 1:2]),
            size = inputFormula$Cl,
            replace = TRUE,
            prob = as.numeric(data['Cl', 5:6])
        ))
        massFluorine <- inputFormula$F * as.numeric(data['F', 1])
        massSilicon <- sum(sample(
            as.numeric(data['Si', 1:3]),
            size = inputFormula$Si,
            replace = TRUE,
            prob = as.numeric(data['Si', 5:7])
        ))
        massMolecule <-
            sum(
                massCarbon,
                massHydrogen,
                massNitrogen,
                massOxygen,
                massSulfer,
                massPhosphorus,
                massBromine,
                massChlorine,
                massFluorine,
                massSilicon
            )
        mz <- massMolecule / abs(charge)
        return(mz)
    }
    sim <- replicate(10000, expr = simulation(inputFormula))
    b <-
        seq(
            from = min(sim) - (1 / (2 * abs(charge))),
            to = max(sim) + 1,
            by = 1 / abs(charge)
        )
    bins <- cut(sim, breaks = b)
    mz <- round(tapply(sim, bins, mean), digits = 5)
    intensity <- as.vector(table(bins))
    spectrum <- data.frame(mz, intensity)
    spectrum <- spectrum[spectrum$intensity != 0,]
    spectrum$percent <-
        with(spectrum, round(intensity / max(intensity) *
                                 100,
                             digits = 2))
    row.names(spectrum) <- seq_len((nrow(spectrum)))
    return(spectrum)
}