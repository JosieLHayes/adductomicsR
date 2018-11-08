library(adductomicsR)
# MS2Dir = system.file("extdata", package = "adductData")
# specSimPepId(MS2Dir = system.file("extdata", package = "adductData"),
# nCores=4, rtDevModels = paste0(system.file("extdata",
# package = "adductData"),'/rtDevModels.RData'))
data = read.csv(paste0(system.file("extdata",package = "adductomicsR"),'/allResults_ALVLIAFAQYLQQCPFEDHVK_example.csv'))
testthat::test_that('specSimPep output is correct', {
    expect_that(data, is_a('data.frame'))
    # expect_equal(nlevels(factor(data$MS2FileName)), length(list.files(MS2Dir,
#     pattern="mzXML")))
})
