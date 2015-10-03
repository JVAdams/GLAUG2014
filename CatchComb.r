# CatchComb.r - 
# apply catchability (availability and selectivity) to the midwater trawl catch
# combine acoustic and (adjusted) midwater trawl data 
#  using classification trees with spatial component (e.g., MTReast, ACMTRnorth, MTRd2sh, MTRbdep)

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

# read in catchability inputs files
wb <- loadWorkbook(paste0(dir, "/Catchability Inputs.xlsx"))
cat <- readWorksheet(wb, sheet=getSheets(wb)[1], startRow=16)


# catchability (selectivity AND availability) functions
dbllogit <- function(len, L501, SR1, L502, SR2) {
  # double logistic regression allowing for parabolic-shaped probabilities
  # function provided by Kresimir Williams, NOAA-AFSC, 22 April 2014
  ((1 + exp(2*log(3)*(L501 - len)/SR1))^-1) * ((1 + exp(2*log(3)*(L502 - len)/-SR2))^-1)
  }
dbllogit2 <- function(fdep, D501, SD1, d2bot, D502, SD2) {
  # double logistic regression allowing for parabolic-shaped probabilities
  # function provided by Kresimir Williams, NOAA-AFSC, 22 April 2014
  ((1 + exp(2*log(3)*(D501 - fdep)/SD1))^-1) * ((1 + exp(2*log(3)*(D502 - d2bot)/SD2))^-1)
  }
if(FALSE) {
  windows()
  par(mfrow=c(2, 2))
  x <- 1:300
  plot(x, dbllogit(x, -Inf, 10, 220, 10), ylab="y", las=1, main="L501 = -Inf")
  plot(x, dbllogit(x, 70, 10, Inf, 10), ylab="y", las=1, main="L502 = Inf")
  x <- 1:100
  plot(x, dbllogit2(x, -Inf, 5, 100-x, 10, 5), ylab="y", las=1, main="D501 = -Inf")
  plot(x, dbllogit2(x, 10, 5, 100-x, -Inf, 5), ylab="y", las=1, main="D502 = -Inf")
  }

# cycle through all runs (all lakes)
ur <- unique(sampinput$run.id)

for(i in seq(ur)) {
  sel.run <- ur[i]
  just <- sampinput[sampinput$run.id==sel.run, ]
  sel.lk <- just$LC

  ac <- read.csv(paste0(subdir, "/ACSummaryIL-lake", sel.lk, "-run", sel.run, ".csv"), as.is=T)
  mtr <- read.csv(paste0(subdir, "/MTRCatch-lake", sel.lk, "-run", sel.run, ".csv"), as.is=T)
  cat.sub <- cat[cat$LC==sel.lk, ]
  mtr.cat <- merge(mtr, cat.sub, by.x="sp", by.y="G", all.x=TRUE)
  mtr.cat$p.avail <- with(mtr.cat, dbllogit2(fdep=f.fdep, D501=D501, SD1=Dslope1, d2bot=f.d2bot, D502=D502, SD2=Dslope2))
  mtr.cat$p.selec <- with(mtr.cat, dbllogit(len=len, L501=L501, SR1=Lslope1, L502=L502, SR2=Lslope2))
  mtr.cat$p.catch <- mtr.cat$p.avail*mtr.cat$p.selec

  # apply catchability (selectivity AND availability) functions to "perfect" MTR catch
  mtr.cat$keep <- sapply(mtr.cat$p.catch, function(p) sample(0:1, size=1, replace=TRUE, prob=c(1-p, p)))

  ue <- sort(unique(ac$Event))
  ug <- sort(unique(mtr.cat$sp))

  results <- expand.grid(lake=sel.lk, group=ug, event=ue, dens=NA, bio=NA)

  for(k in ue) {

    # create a new data frame with just the AC data from Event k
    ack <- ac[ac$Event==k, ]

    # make the variable names the same as in the MTR data
    ack.mnames <- ack
    names(ack.mnames)[match(c("east", "north", "d2sh", "botdep", "fdep", "d2bot"), names(ack))] <- 
      c("MTReast", "ACMTRnorth", "MTRd2sh", "MTRbdep", "MTRfdep", "MTRd2bot")

    # subset the midwater trawl data
    mtr.sub <- mtr.cat[mtr.cat$Event==k & mtr.cat$keep==1, ]

    # only do these calculations if there were "keep" fish in the midwater trawl
    # without "keep" fish, no species-specific density and biomass can be estimated
    if(dim(mtr.sub)[1] > 0) {
      # fit a classification tree to the MTR data for Event k
      treek <- rpart(as.factor(sp) ~ MTReast + ACMTRnorth + MTRd2sh + MTRbdep + MTRfdep + MTRd2bot, 
        data=mtr.sub, control=list(cp=0.05, minsplit=10, minbucket=5))

      # mean density for each species
      # use the fitted tree to predict species composition of each AC interval/layer
      pred.props <- predict(treek, newdata=ack.mnames)
      mean.dens <- apply(apply(ack$nperha*pred.props, 2, tapply, list(ack$ACid, ack$interval), sum), 2, mean)

      # mean biomass for each species
      # calculate the mean weight of each group at each node of the fitted tree
      mwt <- tapply(mtr.sub$wt, list(treek$where, mtr.sub$sp), mean)
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

      results$dens[is.na(results$dens)] <- 0
      results$bio[is.na(results$bio)] <- 0

      }
    }

  write.csv(results, paste0(subdir, "/ResultsSummary-lake", sel.lk, "-run", sel.run, "-scenario", sel.sce, ".csv"), row.names=FALSE)
  }
