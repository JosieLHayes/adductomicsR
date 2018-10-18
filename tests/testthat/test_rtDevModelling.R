library(adductomicsR)
library(adductData)
rtDevModelling(MS2Dir = system.file("extdata", package = "adductData"),
nCores=4, runOrder = paste0(system.file("extdata", package = "adductomicsR"),
 '/runOrder.csv'))
load(paste0(system.file("extdata", package = "adductData"),"/rtDevModels.RData"))
MS2Dir = system.file("extdata", package = "adductData")
test_that('RT deviation models are correct', {
    expect_that(rtDevModels, is_a('list'))
    expect_equal(length(rtDevModels), length(list.files(MS2Dir,
    pattern="mzXML")))
})

 
    
    
