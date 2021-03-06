GLAUG2014
=========

**GLAUG2014** contains files used during the Great Lakes Acoustic Users Group 2014 Workshop on Trawl Performance.
I am currently working on recreating **artiFISHal** as an R package, [here](https://github.com/JVAdams/GLAUG2014).
But the files in **GLAUG2014** should still continue to function as they did during the workshop.

**artiFISHal** is a pelagic fish community simulator developed using the [R programming language](http://www.r-project.org/). 
It can be used to create artificial lakes and populate them with known numbers of fish (identified in species-size groups) to mimic pelagic fish communities. 
**artiFISHal** can then be used to sample these fish with virtual acoustic and midwater trawl surveys. 

[Yule et al. (2013)](http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2013-0072#.U1KYxPldXTQ) used **artiFISHal** 
to evaluate several different approaches for estimating the biomass of pelagic fish species in the Great Lakes.

This information is also available on the U. S. Geological Survey - Great Lakes Science Center [website](http://www.glsc.usgs.gov/artifishal). 
As major changes are made to the code, a new version number will be assigned to **artiFISHal**, and it will be posted to that site.

## Instructions 

#### To use **artiFISHal**, follow the steps outlined below.

1. Download three files from this repository:
 * `Inputs.xls`
 * `SimulateFish.r`
 * `SampleFish.r`

1. Add a folder to the same directory as `Inputs.xls` for storing the output files

1. Edit the three worksheets in the Excel workbook, `Inputs.xls`.
 * `SimLake` - edit the inputs for the shape of the artificial lake you wish to create
 * `SimFish` - edit the inputs for the fish population you wish to create
 * `Sample` - edit the inputs for the acoustic and midwater trawl survey

1. Edit the inputs in the R script `SimulateFish.r`
 * `dir` - directory for reading input data (location of `Inputs.xls`)
 * `subdir` - directory for saving output data
 * `sel.lk` - the lake you wish to create (`LC` in `Inputs.xls` tab `SimFish`)
 * `cap.no.fish` - set the maximum number of fish
 * `save.plots` - specify whether to save plots to a pdf or show them on the screen

1. Download the latest version of [R](http://www.r-project.org/)
 * Install these packages:
  * `XLConnect` - *note that in order for* `XLConnect` *to work in R for Windows 64-bit, you must have Java installed for Windows 64-bit*
  * `class`
  * `MASS`

1. Run the R script `SimulateFish.r`, which will create four output files (one comma delimited data file, one pdf document, and two Rdata files)
 * `Truth-lake#.csv` - number and biomass of true population of fish in lake
 * `Diagnostic-lake#.pdf` - diagnostic plots of fish population
 * `inputs-lake#.Rdata` - collection of R objects
 * `fish-lake#.Rdata` - true population of fish in lake, one row per fish

1. Edit the inputs in the R script `SampleFish.r`
 * `dir` - directory for reading input data (location of `Inputs.xls`)
 * `subdir` - directory for saving output data
 * `sel.run` - the survey run you wish to execute (`run.id` in `Inputs.xls` tab `Sample`)
 * `save.plots` - specify whether to save plots to a pdf or show them on the screen
 * `min.no.fish.per.haul` - set the minimum number of fish required for a valid trawl sample

1. Run the R script `SampleFish.r`, which will create five output files (four comma delimited data files and one pdf document)
 * `ACTargets-lake#-run#.csv` - individual acoustic targets
 * `ACSummaryI-lake#-run#.csv` - summary of acoustic targets by interval (I)
 * `ACSummaryIL-lake#-run#.csv` - summary of acoustic targets by interval (I) and layer (L)
 * `MTRCatch-lake#-run#.csv` - individual fish captured in midwater trawl
 * `Survey-lake#-run#.pdf` - diagnostic plots of acoustics and midwater trawl sampling

## Additional files

#### Several additional files were created for the Great Lakes Acoustic Users Group 2014 Workshop on Trawl Performance
* `Availability.r` - interactive graph of availability relative to water depth for choosing coefficients
* `Selectivity.r` - interactive graph of trawl selectivity relative to fish length for choosing coefficients
* `Catchability Inputs.xlsx` - example spreadsheet for inputting availability and selectivity coefficients
* `CatchComb.r` - apply catchability to the midwater trawl catch, combine with acoustic info (using regression trees), and summarize
* `NearestTrawlByDepthStrata.r` **(NOT USED)** - assign the id of the nearest midwater trawl to each acoustic interval/layer