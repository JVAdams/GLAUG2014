# CatchComb.r - 
# apply catchability (availability and selectivity) to the midwater trawl catch
# combine acoustic and (adjusted) midwater trawl data 
#	using classification trees with spatial component (e.g., MTReast, ACMTRnorth, MTRd2sh, MTRbdep)

##########################################################################
#############################  input values  #############################
##########################################################################

# directory for reading input data (Inputs.xls)
dir <- "C:/Temp"

# directory for reading previously saved Rdata files (inputs-lake# and fish-lake#) 
# and saving new output data (csv files ACTargets-lake#-run#, ACSummaryIL-lake#-run#, ACSummaryI-lake#-run#, 
# MTRCatch-lake#-run#, and pdf file Survey-lake#-run#)
subdir <- "C:/Temp/Output"

# select just one catchability scenario ???
sel.sce <- 1

##########################################################################
############################  data crunching  ############################
##########################################################################

# load the rpart package
if(!require(rpart)) stop("R package rpart is required")

# read in inputs file
wb <- loadWorkbook(paste0(dir, "/Inputs.xls"))

# generate survey design according to the guidelines laid out in the tab referenced below
sampinput <- readWorksheet(wb, sheet="Sample", startRow=17)

# cycle through all runs (all lakes)
ur <- unique(sampinput$run.id)

for(i in seq(ur)) {
	sel.run <- ur[i]
	just <- sampinput[sampinput$run.id==sel.run, ]
	sel.lk <- just$LC

	ac <- read.csv(paste0(subdir, "/ACSummaryIL-lake", sel.lk, "-run", sel.run, ".csv"), as.is=T)
	mtr <- read.csv(paste0(subdir, "/MTRCatch-lake", sel.lk, "-run", sel.run, ".csv"), as.is=T)

	#######################
	# apply catchability (selectivity AND availability) functions to "perfect" MTR catch
	# Catchability Inputs.xls

	dbllogit <- function(len, L501, SR1, L502, SR2) {
		# double logistic regression
		# allowing for parabolic-shaped probabilities
		# function provided by Kresimir Williams, NOAA-AFSC, 22 April 2014
		((1 + exp(2*log(3)*(L501 - len)/SR1))^-1) * ((1 + exp(2*log(3)*(L502 - len)/-SR2))^-1)
		}
	#######################

	ue <- sort(unique(ac$Event))
	ug <- sort(unique(mtr$sp))

	results <- expand.grid(lake=sel.lk, group=ug, event=ue, dens=NA, bio=NA)

	for(k in ue) {

		# fit a classification tree to the MTR data for Event k
		treek <- rpart(as.factor(sp) ~ MTReast + ACMTRnorth + MTRd2sh + MTRbdep + MTRfdep + MTRd2bot, 
			data=mtr[mtr$Event==k, ], control=list(cp=0.05, minsplit=10, minbucket=5))

		# create a new data frame with just the AC data from Event k
		ack <- ac[ac$Event==k, ]

		# make the variable names the same as in the MTR data
		ack.mnames <- ack
		names(ack.mnames)[match(c("east", "north", "d2sh", "botdep", "fdep", "d2bot"), names(ack))] <- 
			c("MTReast", "ACMTRnorth", "MTRd2sh", "MTRbdep", "MTRfdep", "MTRd2bot")

		# mean density for each species
		# use the fitted tree to predict species composition of each AC interval/layer
		pred.props <- predict(treek, newdata=ack.mnames)
		mean.dens <- apply(apply(ack$nperha*pred.props, 2, tapply, list(ack$ACid, ack$interval), sum), 2, mean)

		# mean biomass for each species
		# calculate the mean weight of each group at each node of the fitted tree
		mwt <- tapply(mtr$wt[mtr$Event==k], list(treek$where, mtr$sp[mtr$Event==k]), mean)
		mwt[is.na(mwt)] <- 0
		# character variable of proportions corresponding to each node
		node <- apply(format(round(predict(treek), 6)), 1, paste, collapse="-")[match(sort(unique(treek$where)), treek$where)]
		props <- predict(treek)[match(sort(unique(treek$where)), treek$where), ]
		pred.node <- match(apply(format(round(pred.props, 6)), 1, paste, collapse="-"), node)
		pred.mwtprop <- (mwt*props)[pred.node, ]
		mean.bio <- apply(apply(ack$nperha*pred.mwtprop, 2, tapply, list(ack$ACid, ack$interval), sum), 2, mean)

		for(g in seq(ug)) {
			sel <- results$event==k & results$group==ug[g]
			results$dens[sel] <- mean.dens[ug[g]]
			results$bio[sel] <- mean.bio[ug[g]]/1000
			}
		}

	results$dens[is.na(results$dens)] <- 0
	results$bio[is.na(results$bio)] <- 0

	write.csv(results, paste0(subdir, "/ResultsSummary-lake", sel.lk, "-run", sel.run, "-scenario", sel.sce, ".csv"), row.names=FALSE)
	}
