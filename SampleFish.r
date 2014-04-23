################################################################################################
###                                        artiFISHal                                        ###
### U.S. Geological Survey (USGS) Computer Program "artiFISHal" version 2014-04              ###
### Written by Jean V. Adams, USGS - Great Lakes Science Center, Ann Arbor, Michigan, USA    ###
### Use of this program was first described in a publication by Yule et al. (2013)           ###
###      in the Canadian Journal of Fisheries and Aquatic Sciences                           ###
### Written in programming language R (R Core Team, 2013, www.R-project.org),                ###
###      version 3.1.0 (2014-04-10)                                                          ###
### Run on a PC with Intel(R) Core(TM) I7-4600m CPU @ 2.90 GHz processor, 16.0 GB RAM,       ###
###      and Microsoft Windows 7 Enterprise operating system 2009 Service Pack 1             ###
### Source code is available from Jean V. Adams, jvadams@usgs.gov                            ###
################################################################################################

################################################################################################
###                                        artiFISHal                                        ###
### Disclaimer:                                                                              ###
### Although this program has been used by the USGS, no warranty, expressed or implied, is   ###
### made by the USGS or the United States Government as to the accuracy and functioning of   ###
### the program and related program material nor shall the fact of distribution constitute   ###
### any such warranty, and no responsibility is assumed by the USGS in connection therewith. ###
################################################################################################



################################################################################################
###                                        artiFISHal                                        ###
### SampleFish.r                                                                             ###
### Sample (survey) simulated population of fish (created with SimulateFish.r) with down-    ###
### looking acoustics (cross-lake, west-to-east) and midwater trawls (west-to-east)          ###
################################################################################################



##########################################################################
#############################  input values  #############################
##########################################################################

# directory for reading input data (Inputs.xls)
dir <- "C:/Temp"

# directory for reading previously saved Rdata files (inputs-lake# and fish-lake#) 
# and saving new output data (csv files ACTargets-lake#-run#, ACSummaryIL-lake#-run#, ACSummaryI-lake#-run#, 
# MTRCatch-lake#-run#, and pdf file Survey-lake#-run#)
subdir <- "C:/Temp/Output"

# select just one RUN number for surveying the artificial population
sel.run <- 1

# TRUE = save diagnostic plots to a pdf, FALSE = show them on the screen
save.plots <- TRUE

# set the minimum number of fish per midwater trawl (others will be tossed out)
min.no.fish.per.haul <- 2



##########################################################################
############################  data crunching  ############################
##########################################################################

# load the MASS and XLConnect packages
if(!require(MASS)) stop("R package MASS is required")
if(!require(XLConnect)) stop("R package XLConnect is required")

# read in inputs file
wb <- loadWorkbook(paste0(dir, "/Inputs.xls"))

# generate survey design according to the guidelines laid out in the tab referenced below
sampinput <- readWorksheet(wb, sheet="Sample", startRow=17)
if(sum(duplicated(sampinput$run.id))>0) stop(paste("Each row of the Sample tab in Inputs.xls must have a unique run.id.  No duplicates!"))

# rename input values
runrow <- sampinput$run.id==sel.run
# determine the lake number from the run number
sel.lk <- sampinput$LC[runrow]

# inputs from simulation program
load(file=paste(subdir, "/inputs-lake", sel.lk, ".Rdata", sep=""))

# clear all graphics windows
graphics.off()
if(save.plots) pdf(paste(subdir, "/Survey-lake", sel.lk, "-run", sel.run, ".pdf", sep=""), width=9, height=6.5, 
	title="Survey", paper="USr")

# function to recode values
recode <- function(x, old, new, must.match=TRUE) {
        partial <- match(x, old)
        if(must.match) new[partial] else ifelse(!is.na(partial), new[partial], x)
        }

# function to identify first occurence of different value
first <- function(x) {
    l <- length(x)
    c(1, 1-(x[-1]==x[-l]))
    }


attach(sampinput[runrow, ])

# define number of sampling events
nS.orig <- nS
nS <- max(nS.orig, 2)



###  SAMPLING GEAR & PROCESSING  ###

# increase the number of initial trawls, so we'll still have enough when we get rid of those that 
# hit bottom and those with fewer than 2 fish
nT2 <- ceiling(6*nT/nA)*nA


ACid <- 1:(nA*nS)
ACspace <- diff(northr)/nA

sampinfo <- as.data.frame(matrix(NA, nrow=nT2*nS, ncol=10, 
	dimnames=list(NULL, c("Event", "ACid", "ACnorth", "MTRgrp", "MTRid", "MTReast", "MTRbdep", "MTRfdep", "MTRd2sh", "MTRd2bot"))))
sampinfo$Event <- rep(1:nS, rep(nT2, nS))
sampinfo$ACid <- rep(1:(nA*nS), rep(nT2/nA, nA*nS))
sampinfo$MTRgrp <- rep(1:(nT2/nA), nA*nS)
sampinfo$MTRid <- seq(along=sampinfo$MTRgrp)

set.seed(sampinput$seed[sampinput$run.id==sel.run])
for(k in 1:nS) {
	sel <- sampinfo$Event==k
	ACstart <- runif(1, northr[1], northr[2])
	ACnorth <- sort(unique(c(seq(ACstart, northr[1], -ACspace), seq(ACstart, northr[2], ACspace))))
	sampinfo$ACnorth[sel] <- rep(ACnorth, rep(nT2/nA, nA))
	sampinfo$MTReast[sel] <- runif(nT2, eastr[1]+tow.l/2, eastr[2]-tow.l/2)

	# calculate midwater trawl measures for the MIDPOINT of the trawl
	sampinfo$MTRbdep[sel] <- zfromx(x=sampinfo$MTReast[sel], maxz=maxbotdep, eastr=eastr, ints=ints, slopes=slopes)

	# random selection of midwater trawl depth
	sampinfo$MTRfdep[sel] <- runif(nT2, 0+(trawl.h/2), sampinfo$MTRbdep[sel]-(trawl.h/2))

	sampinfo$MTRd2sh[sel] <- dfromx(x=sampinfo$MTReast[sel], d2shr.we=d2shr.we, eastr=eastr, mid.d=mid.d)
	sampinfo$MTRd2bot[sel] <- sampinfo$MTRbdep[sel] - sampinfo$MTRfdep[sel]

	# make sure acoustic transect cones don't overlap
	mindist.allowed <- botdepr[2]*tan(pi*cone.deg/2/360)
	mindist.observed <- min(diff(sort(ACnorth)))
	if(mindist.observed < mindist.allowed) stop("Acoustic transects are too close together.")
	}

# make sure that all of the acoustic transects are kept ... even if they have no midwater trawls associated with them
ACsampinfo <- sampinfo[first(sampinfo$ACid)==1, 1:3]

# make sure that midwater trawl doesn't hit bottom
MTRbdep1 <- zfromx(x=sampinfo$MTReast - (tow.l/2), maxz=maxbotdep, eastr=eastr, ints=ints, slopes=slopes)
MTRbdep2 <- zfromx(x=sampinfo$MTReast + (tow.l/2), maxz=maxbotdep, eastr=eastr, ints=ints, slopes=slopes)
hits.bottom <- (MTRbdep1 - (sampinfo$MTRfdep + trawl.h/2)) < 0 | (MTRbdep2 - (sampinfo$MTRfdep + trawl.h/2)) < 0

# make sure that midwater trawl doesn't extend too far east/west/north/south
extends.out <- sampinfo$ACnorth - trawl.w/2 < northr[1] |
	sampinfo$ACnorth + trawl.w/2 > northr[2] |
	sampinfo$MTReast - tow.l/2 < eastr[1] |
	sampinfo$MTReast + tow.l/2 > eastr[2]

sampinfo <- sampinfo[!hits.bottom & !extends.out, ]

rm(ACid, ACspace, sel, ACstart, ACnorth, mindist.allowed, mindist.observed, MTRbdep1, MTRbdep2, hits.bottom, extends.out)



attach(ACsampinfo)



# acoustic transects

# convert half of the down-looking angle from degrees to radians
half.cone.rad <- 2*pi*(cone.deg/2)/360

# select only those targets within the volume of space sampled by the acoustic transect (triangular prism)
sua <- sort(unique(ACid))
AC <- data.frame(matrix(NA, nrow=0, ncol=13, dimnames=list(NULL, c("Event", "ACid", "ACnorth", 
	"sp", "f.east", "f.north", "f.d2sh", "f.botdep", "f.fdep", "f.d2bot", "len", "wt", "ts"))))
# acoustic slice cone cushion, based on angle of transducer and maximum depth
cushion <- botdepr[2]*tan(half.cone.rad)

# simulated population (fish)
load(file=paste(subdir, "/fish-lake", sel.lk, ".Rdata", sep=""))

# create a single row of missing values that looks just like the "fish" data frame
# this will be used to maintain a row of information on an acoustic transect, even if it captures no targets
nofish <- fish[1, ]
ncol <- dim(fish)[2]
nofish[1, 1:ncol] <- rep(NA, ncol)
rm(ncol)

# acoustic transects
for(j in sua) {
	# make sure that the acoustic slice does not extend further north or south than our sample space
	if(ACnorth[j] > (northr[2] - cushion) | ACnorth[j] < (northr[1] + cushion)) 
		stop(paste("ACid = ", ACid[j], ", ACnorth = ", ACnorth[j], 
			": The cone of the acoustic slice extends farther north or south than the boundary of our simulated lake."))
	sel <- (abs(fish$f.north - ACnorth[j])) < (fish$f.fdep*tan(half.cone.rad))

	if(sum(sel)>0) {
		temp <- data.frame(ACsampinfo[rep(j, sum(sel)), c("Event", "ACid", "ACnorth")], fish[sel, ])
		AC <- rbind(AC, temp)
		} else {
		temp <- data.frame(ACsampinfo[rep(j, 1), c("Event", "ACid", "ACnorth")], nofish)
		AC <- rbind(AC, temp)
		}
	}

rm(cushion)

detach(ACsampinfo)
attach(sampinfo)



# midwater trawl tows
# select only those targets within the volume of space sampled by the midwater trawl (rectangular prism)
sut <- sort(unique(MTRid))
MTRbig <- data.frame(matrix(NA, nrow=0, ncol=20, dimnames=list(NULL, 
	c("Event", "ACid", "ACnorth", "MTRgrp", "MTRid", "MTReast", "MTRbdep", "MTRfdep", "MTRd2sh", "MTRd2bot", 
	"sp", "f.east", "f.north", "f.d2sh", "f.botdep", "f.fdep", "f.d2bot", "len", "wt", "ts"))))
for(m in sut) {
	j <- match(m, MTRid)
	sel <- fish$f.north >= (ACnorth[j] - trawl.w/2) & fish$f.north <= (ACnorth[j] + trawl.w/2) &
		fish$f.east >= (MTReast[j] - tow.l/2) & fish$f.east <= (MTReast[j] + tow.l/2) &
		fish$f.fdep <= (MTRfdep[j] + trawl.h/2) & fish$f.fdep >= (MTRfdep[j] - trawl.h/2)
	if(sum(sel) >= min.no.fish.per.haul) {
		temp <- cbind(sampinfo[rep(j, sum(sel)), ], fish[sel, ])
		MTRbig <- rbind(MTRbig, temp)
		}
	}



### select only needed number of trawls (remember, we sampled many more than what we needed!) ##

# create an index of trawls
mtr.indx <- MTRbig[first(MTRbig$MTRid)==1, c("Event", "ACid", "MTRgrp", "MTRid")]
mtr.indx$mtr.grp <- unlist(lapply(table(mtr.indx$ACid), sample))
# create a random variable for sorting
y <- runif(length(mtr.indx$mtr.grp))
# sort the index of trawls by event abd group number
mtr.indx.sort <- mtr.indx[order(mtr.indx$Event, mtr.indx$mtr.grp, y), ]

# count up the number of trawls in each event
f <- first(mtr.indx.sort$Event)
x <- rep(NA, length(f))
for(ix in 1:length(f)) {
	x[ix] <- if(f[ix]==1) 1 else x[ix-1]+1
	}
mtr.indx.sort$mtr.count <- x
# keep the first nT trawls in each event
mtr.indx.keep <- mtr.indx.sort[mtr.indx.sort$mtr.count <= nT, ]

MTR <- MTRbig[MTRbig$MTRid %in% mtr.indx.keep$MTRid, ]
old.n <- c("ACnorth")
new.n <- c("ACMTRnorth")
names(MTR)[match(old.n, names(MTR))] <- new.n

sampinfo.sub <- sampinfo[sampinfo$MTRid %in% mtr.indx.keep$MTRid, ]

rm(mtr.indx, y, mtr.indx.sort, f, x, mtr.indx.keep, old.n, new.n)



detach(sampinfo)



# categorize targets interval (along transect) and layer (in water column), both in m

int.breaks <- seq(eastr[1], eastr[2]+interval.-1, interval.)
int.mids <- int.breaks[-1] - diff(int.breaks)/2
AC$interval <- cut(AC$f.east, include.lowest=TRUE, breaks=int.breaks, labels=FALSE)
AC$east <- int.mids[AC$interval]

lay.breaks <- seq(0, botdepr[2]+layer.-1, layer.)
lay.mids <- lay.breaks[-1] - diff(lay.breaks)/2
AC$layer <- cut(AC$f.fdep, include.lowest=TRUE, breaks=lay.breaks, labels=FALSE)
AC$fdep <- lay.mids[AC$layer]



# create a matrix of all possible interval-by-layer combinations for each ACid
# this will be used to ensure that interval-by-layers with no fish are included in later summaries

# each AC transect is the same length, so each will have all of the possible intervals
all.ints <- seq(from=length(int.mids))

# each AC transect goes over the same depth profile, so we can determine the max depth for each interval
# then see whether the max layer included

easts <- seq(eastr[1], eastr[2], length=1000)
depth.contour <- zfromx(x=easts, maxz=maxbotdep, eastr=eastr, ints=ints, slopes=slopes)
depth.cont.int <- as.numeric(cut(easts, include.lowest=TRUE, breaks=int.breaks))
all.maxes <- tapply(depth.contour, depth.cont.int, max)

max.lays <- as.numeric(cut(all.maxes, include.lowest=TRUE, breaks=lay.breaks))

# full matrix of all ACids, all intervals, and all layers
full.mat <- expand.grid(layer=1:max(max.lays), interval=all.ints, ACid=ACsampinfo$ACid)
full.mat <- merge(full.mat, ACsampinfo[, c("ACid", "Event", "ACnorth")])

sub.mat <- merge(data.frame(interval=all.ints, max.layer=max.lays), full.mat, all=TRUE)
sub.mat <- sub.mat[sub.mat$layer <= (sub.mat$max.layer + 0.5), c("Event", "ACid", "ACnorth", "interval", "layer")]
sub.mat$east <- int.mids[sub.mat$interval]
sub.mat$fdep <- lay.mids[sub.mat$layer]
old.n <- c("ACnorth")
new.n <- c("north")
names(sub.mat)[match(old.n, names(sub.mat))] <- new.n

rm(easts, depth.contour, depth.cont.int, all.maxes, max.lays, full.mat) #all.ints



# add "range weight" to the AC data (Yule 2000)
# a weighting variable to account for different volumes sampled as a function of range (dist. from ducer in m) and ducer half angle
AC$rng.wt <- 1/(2*AC$f.fdep*tan(half.cone.rad))

# summarize acoustic data by interval and layer
AC$indx <- paste(AC$ACid, AC$interval, AC$layer, sep="-")
sum.rw <- tapply(AC$rng.wt, AC$indx, sum)
ord <- match(names(sum.rw), AC$indx)
ACsmryIL <- AC[ord, c("Event", "ACid", "interval", "layer", "ACnorth", "east", "fdep")]
ACsmryIL$sum.rw <- sum.rw
rm(sum.rw, ord)
ACsmryIL <- ACsmryIL[order(ACsmryIL$ACid, ACsmryIL$interval, ACsmryIL$layer), ]

old.n <- c("ACnorth")
new.n <- c("north")
names(ACsmryIL)[match(old.n, names(ACsmryIL))] <- new.n
ACsmryIL$nperha <- 10000 * ACsmryIL$sum.rw/interval.
rm(old.n, new.n)

ACsmryIL <- merge(ACsmryIL, sub.mat, all=TRUE)
ACsmryIL$nperha[is.na(ACsmryIL$nperha)] <- 0
ACsmryIL$sum.rw[is.na(ACsmryIL$sum.rw)] <- 0

ACsmryIL$botdep <- zfromx(x=ACsmryIL$east, maxz=maxbotdep, eastr=eastr, ints=ints, slopes=slopes)
ACsmryIL$d2bot <- ACsmryIL$botdep - ACsmryIL$fdep
ACsmryIL$d2sh <- dfromx(x=ACsmryIL$east, d2shr.we=d2shr.we, eastr=eastr, mid.d=mid.d)

# summarize acoustic data by interval only
ACsmryI <- aggregate(ACsmryIL[, c("sum.rw", "nperha")], 
	ACsmryIL[, c("Event", "ACid", "interval", "north", "east", "botdep", "d2sh")], sum)
ACsmryI <- ACsmryI[order(ACsmryI$Event, ACsmryI$ACid, ACsmryI$interval), ]



ac.out.vars <- c("Event", "ACid", "interval", "layer", "f.east", "f.north", "f.d2sh", "f.botdep", "f.fdep", "f.d2bot", "ts", "rng.wt")
acsumil.out.vars <- c("Event", "ACid", "interval", "layer", "east", "north", "d2sh", "botdep", "fdep", "d2bot", "nperha")
acsumi.out.vars <- c("Event", "ACid", "interval", "east", "north", "d2sh", "botdep", "nperha")
mtr.out.vars <- c("Event", "ACid", "MTRid", "MTReast", "ACMTRnorth", "MTRd2sh", "MTRbdep", "MTRfdep", "MTRd2bot", "sp", "len", "wt")

#### only save *1* sampling event, if that is what was originally specified ###
if(nS.orig==1) {
	write.csv(AC[AC$Event==1, ac.out.vars], 
		paste(subdir, paste("ACTargets-lake", sel.lk, "-run", sel.run, ".csv", sep=""), sep="/"), row.names=FALSE)
	write.csv(ACsmryIL[ACsmryIL$Event==1, acsumil.out.vars], 
		paste(subdir, paste("ACSummaryIL-lake", sel.lk, "-run", sel.run, ".csv", sep=""), sep="/"), row.names=FALSE)
	write.csv(ACsmryI[ACsmryI$Event==1, acsumi.out.vars], 
		paste(subdir, paste("ACSummaryI-lake", sel.lk, "-run", sel.run, ".csv", sep=""), sep="/"), row.names=FALSE)
	write.csv(MTR[MTR$Event==1, mtr.out.vars], 
		paste(subdir, paste("MTRCatch-lake", sel.lk, "-run", sel.run, ".csv", sep=""), sep="/"), row.names=FALSE)
	} else {
	write.csv(AC[, ac.out.vars], 
		paste(subdir, paste("ACTargets-lake", sel.lk, "-run", sel.run, ".csv", sep=""), sep="/"), row.names=FALSE)
	write.csv(ACsmryIL[, acsumil.out.vars], 
		paste(subdir, paste("ACSummaryIL-lake", sel.lk, "-run", sel.run, ".csv", sep=""), sep="/"), row.names=FALSE)
	write.csv(ACsmryI[, acsumi.out.vars], 
		paste(subdir, paste("ACSummaryI-lake", sel.lk, "-run", sel.run, ".csv", sep=""), sep="/"), row.names=FALSE)
	write.csv(MTR[, mtr.out.vars], 
		paste(subdir, paste("MTRCatch-lake", sel.lk, "-run", sel.run, ".csv", sep=""), sep="/"), row.names=FALSE)
	}
rm(ac.out.vars, acsumil.out.vars, acsumi.out.vars, mtr.out.vars)



attach(sampinfo.sub)


sel <- AC$Event==1
ts.scaled <- (AC$ts[sel] - tsr[1])/diff(tsr)
sua <- sort(unique(AC$ACid[sel]))
sut <- sort(unique(MTRid[Event==1]))


# top view of acoustic transects and midwater trawls - to scale
if(!save.plots) windows(w=9, h=6.5)
par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(5.1, 4.1, 4.1, 2.1))
eqscplot(1, 1, type="n", xlim=eastr/1000, ylim=northr/1000, axes=FALSE, xlab="Easting  (km)", ylab="Northing  (km)", 
	main=paste(lkinfo$Lake[sel.lk], "- Top View - drawn to scale"))
axis(1)
axis(2)
arrows(rep(eastr[1], length(sua))/1000, AC$ACnorth[match(sua, AC$ACid[sel])]/1000, 
	rep(eastr[2], length(sua))/1000, AC$ACnorth[match(sua, AC$ACid[sel])]/1000, length=0, lwd=2, col="blue") 
polygon(eastr[c(1, 2, 2, 1)]/1000, northr[c(1, 1, 2, 2)]/1000)
if(length(sut)>0) {
	for(m in seq(along=sut)) {
		sel2 <- MTRid==sut[m]
		polygon((MTReast[sel2]+c(-1, 1, 1, -1)*tow.l/2)/1000, (ACnorth[sel2]+c(1, 1, -1, -1)*trawl.w/2)/1000, border="red", lwd=3)
		}
	rm(sel2)
	}


# plot of each acoustic transect with outline of midwater trawl tows
if(!save.plots) windows(w=9, h=6.5)
par(mfrow=n2mfrow(length(sua)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
catch.tots <- table(MTR$MTRid)
for(j in seq(along=sua)) {
	plot(AC$f.east/1000, -AC$f.fdep, type="n", xlab="", ylab="")
	sel3 <- AC$ACid==sua[j]
	points(AC$f.east[sel3]/1000, -AC$f.fdep[sel3], cex=2*ts.scaled, col="blue")
	points(AC$f.east[sel3]/1000, -AC$f.botdep[sel3], pch=16, cex=0.5)
	sut2 <- sort(unique(MTRid[ACid==sua[j]]))
	mtext(paste("id=", sua[j], "\nn=", sum(sel3), sep=""), side=1, line=-1.5, adj=0.98, font=2, cex=par("cex"))
	if(length(sut2)>0) {
		for(m in seq(along=sut2)) {
			sel2 <- MTRid==sut2[m]
			polygon((MTReast[sel2]+c(-1, 1, 1, -1)*tow.l/2)/1000, -MTRfdep[sel2]+c(1, 1, -1, -1)*trawl.h/2, border="red", col="white")
			text(MTReast[sel2]/1000, -MTRfdep[sel2], catch.tots[sel2], font=2, col="red")
			}
		rm(sel2)
		}
	}
mtext("Easting  (km)", side=1, outer=TRUE)
mtext("Fishing depth  (m)", side=2, outer=TRUE)
mtext(paste(lkinfo$Lake[sel.lk], "- Acoustic Transects and Midwater Trawls"), side=3, outer=TRUE, font=2)


# plot of each midwater trawl catch
if(length(sut)>0) {
	barz <- tapply(!is.na(MTR$len), list(MTR$sp, 50*floor(MTR$len/50), MTR$MTRid), sum)
	barz[is.na(barz)] <- 0
	barz <- barz[, , dimnames(barz)[[3]] %in% sut, drop=FALSE]
	barz <- barz[apply(barz, 1, sum)>0, , , drop=FALSE]
	sus <- dimnames(barz)[[1]]
	if(!save.plots) windows(w=9, h=6.5)
	par(mfrow=n2mfrow(length(sut)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
	for(m in seq(along=sut)) {
		barplot(barz[, , m], ylim=c(0, max(apply(barz, 2:3, sum))), col=1:50, 
			names.arg=paste(dimnames(barz[, , 1])[[2]], "+", sep=""))
		mtext(paste("id=", sut[m], "\nn=", sum(barz[, , m]), sep=""), side=3, line=-2.5, adj=0.98, font=2, cex=par("cex"))
		box()
		}
	mtext("Length  (mm)", side=1, outer=TRUE)
	mtext("Frequency", side=2, outer=TRUE)
	mtext(paste(lkinfo$Lake[sel.lk], "- Midwater Trawl Catch"), side=3, outer=TRUE, font=2)
	legend("topleft", sus, fill=seq(sus))
	rm(barz, sus)
	}

detach(sampinfo.sub)

detach(sampinput[runrow, ])

if(save.plots) graphics.off()



# mean number of targets per transect
mean(table(AC$ACid))

# mean number of fish per trawl
mean(table(MTR$MTRid))

# remove all objects from current working directory
rm(AC, ACsampinfo, ACsmryI, ACsmryIL, botdepr, catch.tots, d2shr.we, dfromx, eastr, first, fish, half.cone.rad, 
	ints, ix, j, k, lkinfo, m, maxbotdep, mid.d, min.no.fish.per.haul, MTR, MTRbig, nofish, northr, nS, 
	nS.orig, nT2, recode, sampinfo, sampinfo.sub, sampinput, save.plots, sel, sel.lk, sel.run, sel3, slopes, sua, dir, subdir, sub.mat, 
	sut, sut2, temp, ts.scaled, tsr, xfromz, zfromx, wb, runrow)
