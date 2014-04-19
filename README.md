artiFISHal
==========

artiFISHal is a pelagic fish community simulator developed using the R programming language. It can be used to create artificial lakes and populate them with known numbers of fish (identified in species-size groups) to mimic pelagic fish communities. artiFISHal can then be used to sample these fish with virtual acoustic and midwater trawl surveys. Yule et al. (2013) used artiFISHal to evaluate several different approaches for estimating the biomass of pelagic fish species in the Great Lakes.

# Instructions 

## To use artiFISHal, follow the steps outlined below.

1. Download three files from this repository:
 * Inputs.xls
  * SimulateFish.r
 * SampleFish.r

1. Edit the three worksheets in the Excel workbook, Inputs.xls.
..* SimLake -- edit the inputs for the shape of the artificial lake you wish to create
..* SimFish -- edit the inputs for the fish population you wish to create
..* Sample -- edit the inputs for the acoustic and midwater trawl survey

1. Edit the inputs in the R script SimulateFish.r
..* dir -- directory for reading input data (location of Inputs.xls)
..* subdir -- directory for saving output data
..* sel.lk -- the lake you wish to create (LC in Inputs.xls tab SimFish)
..* cap.no.fish -- set the maximum number of fish
..* save.plots -- specify whether to save plots to a pdf or show them on the screen

1. Run the R script SimulateFish.r, which will create four output files, one comma delimited data file, one pdf document, and two Rdata files
..* Truth-lake#.csv -- number and biomass of true population of fish in lake
..* Diagnostic-lake#.pdf -- diagnostic plots of fish population
..* inputs-lake#.Rdata -- collection of R objects
..* fish-lake#.Rdata -- true population of fish in lake, one row per fish

1. Edit the inputs in the R script SampleFish.r
..* dir -- directory for reading input data (location of Inputs.xls)
..* subdir -- directory for saving output data
..* sel.run -- the survey run you wish to execute (run.id in Inputs.xls tab Sample)
..* save.plots -- specify whether to save plots to a pdf or show them on the screen
..* min.no.fish.per.haul -- set the minimum number of fish required for a valid trawl sample

1. Run the R script SampleFish.r, which will create five output files, four comma delimited data files and one pdf document
..* ACTargets-lake#-run#.csv -- individual acoustic targets
..* ACSummaryI-lake#-run#.csv -- summary of acoustic targets by interval (I)
..* ACSummaryIL-lake#-run#.csv -- summary of acoustic targets by interval (I) and layer (L)
..* MTRCatch-lake#-run#.csv -- individual fish captured in midwater trawl
..* Survey-lake#-run#.pdf -- diagnostic plots of acoustics and midwater trawl sampling
