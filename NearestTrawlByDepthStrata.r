# NearestTrawlByDepthStrata.r

##########################################################################
#############################  input values  #############################
##########################################################################

# directory for reading input data (Inputs.xls)
dir <- "C:/Temp"

# directory for reading previously saved Rdata files (inputs-lake# and fish-lake#) 
# and saving new output data (csv files ACTargets-lake#-run#, ACSummaryIL-lake#-run#, ACSummaryI-lake#-run#, 
# MTRCatch-lake#-run#, and pdf file Survey-lake#-run#)
subdir <- "C:/Temp/Output"

# input depth cutoff to use for each lake
# for example, if you have three lakes, you might specify depth cutoffs like this
depthcutoff <- c(30, 20, 20)

##########################################################################
############################  data crunching  ############################
##########################################################################

# load the MASS, XLConnect, and class packages
if(!require(MASS)) stop("R package MASS is required")
if(!require(XLConnect)) stop("R package XLConnect is required")
if(!require(class)) stop("R package class is required")

# read in inputs file
wb <- loadWorkbook(paste0(dir, "/Inputs.xls"))

# get unique collection of lakes and runs
sampinput <- readWorksheet(wb, sheet="Sample", startRow=17)
if(sum(duplicated(sampinput$run.id))>0) stop(paste("Each row of the Sample tab in Inputs.xls must have a unique run.id.  No duplicates!"))

allfiles <- list.files(subdir)

for(k in 1:dim(sampinput)[1]) {
	lakeruncsvname <- paste0("-lake", sampinput$LC[k], "-run", sampinput$run.id[k], ".csv")

	ACsmryILname <- paste0("ACSummaryIL", lakeruncsvname)
	MTRname <- paste0("MTRCatch", lakeruncsvname)
	if(ACsmryILname %in% allfiles & MTRname %in% allfiles) {

		ACsmryIL <- read.csv(paste(subdir, ACsmryILname, sep="/"))
		MTR <- read.csv(paste(subdir, MTRname, sep="/"))

		# define number of sampling events
		nS <- sampinput$nS[k]

		# define depth strata
		zstrat <- depthcutoff[sampinput$LC[k]]

		# summarize the midwater trawl location (like sampinfo in earlier programs)
		suM <- sort(unique(MTR$MTRid))
		sampinf <- MTR[match(suM, MTR$MTRid), c("Event", "ACid", "MTRid", "MTReast", "ACMTRnorth", "MTRfdep")]

		# find the nearest midwater trawl to each acoustic interval/layer by x and y within z strata
		# calculations must be done separately for each sampling event and depth stratum
		ACsmryIL$nearest.MTRid <- NA
		for(i in 1:nS) {
			selAC <- ACsmryIL$Event==i & !is.na(ACsmryIL$fdep) & (ACsmryIL$fdep + 5) <= zstrat
			selsi <- sampinf$Event==i & sampinf$MTRfdep <= zstrat
			# if there are no MTRs in that depth zone, use all MTRs for that event
			if(sum(selsi) < 0.5) selsi <- sampinf$Event==i
			# if there are no MTRs in that event, no MTR can be assigned to that AC
			if(sum(selsi) > 0.5) {
				ACsmryIL$nearest.MTRid[selAC] <- as.numeric(as.character(knn1(sampinf[selsi, c("MTReast", "ACMTRnorth")], 
					ACsmryIL[selAC, c("east", "north")], sampinf$MTRid[selsi])))
				}

			selAC <- ACsmryIL$Event==i & !is.na(ACsmryIL$fdep) & (ACsmryIL$fdep + 5) > zstrat
			selsi <- sampinf$Event==i & sampinf$MTRfdep > zstrat
			# if there are no MTRs in that depth zone, use all MTRs for that event
			if(sum(selsi) < 0.5) selsi <- sampinf$Event==i
			# if there are no MTRs in that event, no MTR can be assigned to that AC
			if(sum(selsi) > 0.5) {
				ACsmryIL$nearest.MTRid[selAC] <- as.numeric(as.character(knn1(sampinf[selsi, c("MTReast", "ACMTRnorth")], 
					ACsmryIL[selAC, c("east", "north")], sampinf$MTRid[selsi])))
				}
			}

		# add size of catch from nearest trawl to ACsmryIL
		ACsmryIL$near.mtr.catch <- recode(ACsmryIL$nearest.MTRid, sort(unique(MTR$MTRid)), table(MTR$MTRid))

		write.csv(ACsmryIL, paste(subdir, paste("ACSummaryILwNearestTrawl", lakeruncsvname, sep=""), sep="/"))
		}
}
