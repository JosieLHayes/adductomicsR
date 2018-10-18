adductomicsR
===========
A package created Josie Hayes.

Processes peptide adducts from LC-MS data. The software interrogates tandem mass spectra to perform retention time drift corrections, untargeted putative adduct detection with MS2 plot outputs for quality control, and peak area quantification relative to a chemically neutral peptide.  
 The **adductomicsR** package performs peak-picking, retention time alignment, grouping, peak table output pre-processing, and peak area quantification. 

Installation
===============

The package requires an R version of at least 3.4. It is also recommended to install mzR manually prior to installing the package if it is not already installed on your system.
```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("mzR")
```

Note that if you receive an error "fatal error: netcdf.h" in installing mzR, homebrew has not been able to install netcdf. To install mzR download it from  https://github.com/Unidata/netcdf-c/releases/v4.4.1.1 (using configure, make, make install) and then install mzR again.


The latest development version and all other package dependencies can be installed with one-line of code directly from GitHub using the devtools package. First ensure devtools is installed, instructions can be found here: https://github.com/hadley/devtools
```{r}
library(devtools)
devtools::install_github('JosieLHayes/adductomicsR', dependencies=c("Depends", "Imports", "Suggests"))
```

Vignette
========
The vignette can be viewed here https://github.com/JosieLHayes/adductomicsR/blob/master/vignettes/adductomicsRWorkflow.Rmd with 2 example mzXML files acquired on a LTQ Orbitrap XL HRMS coupled with a Dionex Ultimate® 3000 nanoflow LC system via a Flex Ion nano-electrospray-ionization source and converted to mzXML using MSConvert(http://proteowizard.sourceforge.net/)).

The *adductomicsR* package has thus far only been tested with a LTQ Orbitrap XL HRMS on computers running Windows, OSX and Linux operating systems but depending on interest could be readily extended to other instrument manufacturers.

Details 
=======
The R package utilizes [xcms](https://bioconductor.org/packages/release/bioc/html/xcms.html), [CAMERA](https://bioconductor.org/packages/release/bioc/html/CAMERA.html), [MetMSLine](https://github.com/WMBEdmands/MetMSLine), and many other packages to attempt to implement LC-MS adduct identification and quantification.

A major impetus for development of this package was to provide an open-source pipeline to identify protein adducts on a peptide of interest. Our laboratory has extensive experience in identification and quantification of putative adducts to the Cys34 of human serum albumin (https://www.ncbi.nlm.nih.gov/pubmed/27684351, https://www.ncbi.nlm.nih.gov/pubmed/27936627, https://www.ncbi.nlm.nih.gov/pubmed/29350914, https://www.ncbi.nlm.nih.gov/pubmed/29538615). These analyses used Xcalibur https://www.thermofisher.com/order/catalog/product/OPTON-30487, a proprietry software from Thermo Fisher Scientific, to acquire MS1 and MS2 spectra.  


******
The *adductomicsR* workflow consists of a retention time correction step (optional), `rtDevModeling`, a adduct identification step `specSimPepId`, and a putative adduct quantification step `adductQuant`. A target table can be created for `adductQuant` from the results of `specSimPepId` using `generateTargTable` and the `adductQuant` result object can be processed and filtered using the `outputPeaktable` and `filterAdductTable` respectively.

******

**rtDevModeling** - Performs MS/MS spectrum grouping and loess retention time deviation modeling. Requires as imput a directory path where the mzXML files are and a path to a run order file. Examples mzXML files are available at https://berkeley.box.com/s/fnhttc87v4mn1x50nvckpt99999y7uhl and an example run order file for these examples is in inst/extdata. Information on the internal standard (for Cys34 we use isotopic T3 adducted with iodoacetamide) must be provided here - a list (no white space) of expected fragment ions
for the internal standard spectrum and the expected mass-to-charge ratio of the internal standard precursor (default = 834.77692, for Cys34)  In addition the internal standard retention time drift window (in seconds) can be specified by the user (default 200-600 ppm). 

This function produces a plot for the internal standard RT, ppm difference and deviation from the median across the run order to highlight retention time drift. This plot from a previous dataset https://www.ncbi.nlm.nih.gov/pubmed/27936627 shows that towards the end of the run the retention time drops.

![example_RT](https://github.com/JosieLHayes/adductomicsR/blob/master/inst/extdata/internalStandard_plots.png)

A plot of the adjusted retention time for each retention time (seconds) of this study shows that retention time deviated specifically at certain times. This may indicate caution should be taken when results are reported at these retention times, and may be due to washes and instrument related artifacts that occur during the run.

![example_RTdev](https://github.com/JosieLHayes/adductomicsR/blob/master/inst/extdata/adjRtPlot.png)

******

**specSimPepId** performs spectral similarity based adducted peptide identification. It takes as input the `rtDevModeling` object and a directory path where the mzXML files are. A retention time window within which to identify spectra can be specified using minRT and maxRT (default 20-45 minutes). Similarly a mass-to-charge window can also be specified using minMz and maxMz (defaults 750-1000). A model spectrum file for the peptide under study must be provided to perform spectral similarity to. Built in model tables (in the extdata directory) can be used by specifying the path to the table (currently available are: "ALVLIAFAQYLQQCPFEDHVK" and "RHPYFYAPELLFFAK"). If supplying a custom table it must consist of the following mandatory columns ("mass", "intensity", "ionType" and "fixed or variable"). This function also performs grouping using hierarchical clustering of the spectra. The mass-to-charge ratio and RT threshold for cutting the tree can be specified using groupMzabs and groupRtDev respectively. 

This function produces an MS2 plot for each adduct in each scan. This is saved in the mzXML directory in a separate directory for each sample ending in _adductID. These should be used to visually inspect 2-3 plots for each adduct group identified to remove false positives. A plot of the model spectrum provided is also saved in the mzXML directory for comparison. An example plot for adduct A40 from dataset https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5555296/ is shown below.

![example_MS2](https://github.com/JosieLHayes/adductomicsR/blob/master/inst/extdata/ORB35018.mzXML%20scan%201616%20M870.4355_RT26.39%20dp%200.97%20varPeakDet%208.png)

In addition a plot of the mass-to-charge vs the RT and adjusted RT is produced by this function. Each group, assigned using the grouping thresholds the user provided) is colored differently. These plots are provided within the mzXML directory in a directory labeled spectrumGroups_[peptide]. The plot of all groups for the dataset https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5555296/ is shown below. This shows that some groups, such as those at m/z 850, should be merged into one group as they represent tails of the same peak.

![example_groups](https://github.com/JosieLHayes/adductomicsR/blob/master/inst/extdata/allGroups.png)

**generateTargTable** can be used to generate a target table from these results for the quantification step. It is recommended that the MS2 plots and spectrum grouping plots are used to remove false positives and merge groups that are tails of the same peak prior to quantification.


******

**adductQuant** quantifies putative adducts by peak area. The putative adducts must be provided in the form of a target table which can be manually generated or produced from **generateTargTable**. Two example target tables are provided in inst/extdata. The **rtDevModeling** object should also be provided. The maximum parts per million to be used for peak integration is specified by the user (default 4), increasing this will merge peaks and lower resolution. The number of scans that a spike must be seen in for it to be integrated as a peak can also be specified with spikeScans (default 2). The maximum retention time drift default is 20 seconds and can be altered by the user, and the maximum retention time window to search in is set at 120 seconds. A string for the amino acid sequence of a chemically neutral peptide ('housekeeping peptide') of the protein under study must also be provided. The default is LVNEVTEFAK for Cys34. It is recommended to also include this in the target table (automatically done using the generateTargTable function) so that peak area ratios relative to the housekeeping peptide can be calculated. The result is an adductQuant object. This can be converted to a peak table using **outputPeakTable** and filtered using **filterAdductTable**.


******


License
=============
The *adductomicsR* package is licensed under Artistic License 2.0

 
