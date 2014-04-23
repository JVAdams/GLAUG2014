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
### SimulateFish.r                                                                           ###
### Create a simulated population of fish in an artificial lake                              ###
################################################################################################



##########################################################################
#############################  input values  #############################
##########################################################################

# directory for reading input data (Inputs.xls)
dir <- "C:/Temp"

# directory for saving output data 
# (.Rdata file inputs-lake# and fish-lake#, pdf file Diagnostics-lake#, and csv file Truth-lake#)
subdir <- "C:/Temp/Output"

# select just one lake number (LC) for creation of artificial population
sel.lk <- 1

# cap the total number of fish at 5 million for a memory of ~ 2 GB (2047 MB)
# this limit can be increased if you have more memory available in R
# you can check the memory available by submitting this command: memory.limit()
cap.no.fish <- 5e6

# TRUE = save diagnostic plots to a pdf, FALSE = show them on the screen
save.plots <- TRUE



##########################################################################
############################  data crunching  ############################
##########################################################################

# load packages
if(!require(MASS)) stop("R package MASS is required")
if(!require(XLConnect)) stop("R package XLConnect is required")

# clear all graphics windows
graphics.off()
if(save.plots) pdf(paste(subdir, "/Diagnostics-lake", sel.lk, ".pdf", sep=""), width=9, height=6.5, title="Diagnostics", paper="USr")

# read in inputs file
wb <- loadWorkbook(paste0(dir, "/Inputs.xls"))

# generate lakes according to the guidelines laid out in the tab referenced below
lkinfo <- readWorksheet(wb, sheet="SimLake", startRow=17)
if(any(lkinfo$LC != lkinfo$LC)) stop("Lake codes (LC) must be sequentially numbered: 1, 2, 3, etc.")

# generate fish (targets) according to the guidelines laid out in the tab referenced below
spinfo <- readWorksheet(wb, sheet="SimFish", startRow=22)

# determine distance to shore given easting
dfromx <- function(x, d2shr.we, eastr, mid.d) {
	d <- rep(NA, length(x))
	d[!is.na(x) & x <= mid.d] <-  d2shr.we[1] + x[!is.na(x) & x <= mid.d]
	d[!is.na(x) & x >  mid.d] <- eastr[2] + d2shr.we[2] - x[!is.na(x) & x > mid.d]
	d
	}

# determine bottom depth given easting
zfromx <- function(x, maxz, eastr, ints, slopes) {
	z <- rep(NA, length(x))
	z[!is.na(x) & x <= eastr[2]/3] <- slopes[1]*x[!is.na(x) & x <= eastr[2]/3] + ints[1]
	z[!is.na(x) & x >  eastr[2]/3] <- slopes[2]*x[!is.na(x) & x >  eastr[2]/3] + ints[2]
	ifelse(z > maxz, maxz, z)
	}

# determine easting given bottom depth
# incorporate random assignment to west or east shore if not specified
xfromz <- function(z, maxz, ints, slopes, shore="random") {
	x <- rep(NA, length(z))
	side <- if(shore[1]=="random") sample(0:1, length(z), replace=TRUE) else shore
	x[!is.na(z) & side==0] <- (z[!is.na(z) & side==0] - ints[1]) / slopes[1]
	x[!is.na(z) & side==1] <- (z[!is.na(z) & side==1] - ints[2]) / slopes[2]
	x[z >= maxz] <- runif(sum(z >= maxz), (maxz - ints[1]) / slopes[1], (maxz - ints[2]) / slopes[2])
	x
	}


# select data from one lake for creation of artificial population
spinfo <- spinfo[spinfo$LC==sel.lk, ]

# total number of fish
TotNFish <- min(c(lkinfo$TotNFish[sel.lk], cap.no.fish))
spinfo$Nfish <- floor(spinfo$Nprop * TotNFish / sum(spinfo$Nprop))
spinfo <- spinfo[spinfo$Nfish > 0, ]

# minimum bottom depth, maximum bottom depth, and "vertex" bottom depth in m (z direction)
botdepr <- c(lkinfo$BotDepMin[sel.lk], lkinfo$BotDepMax[sel.lk], lkinfo$BotDepVertex[sel.lk])

# assign maximum bottom depth to an object
maxbotdep <- botdepr[2]

# easting (m), range, (x direction)
eastr <- c(0, 1000*lkinfo$LkWidth[sel.lk])

# northing (m), range, (y direction)
northr <- c(0, 1000*lkinfo$LkLength[sel.lk])

# slopes and intercepts of west and east shores
slopes <- c( (botdepr[3]-botdepr[1])/(eastr[2]/3), (botdepr[1]-botdepr[3])/(2*eastr[2]/3) )
ints <- c(botdepr[1], -eastr[2]*slopes[2]+botdepr[1])

# distance from west and east shore excluded from lake in m (x direction)
d2shr.we <- c(0+ints[1]/slopes[1], -ints[2]/slopes[2]-eastr[2])

# distance to shore (m), range, (x direction)
mid.d <- ((eastr[1]-d2shr.we[1]) + (eastr[2]+d2shr.we[2]) ) / 2
d2shr <- c(min(d2shr.we), d2shr.we[1]+mid.d)

# target strength (db), range
tsr <- c(lkinfo$TSMin[sel.lk], lkinfo$TSMax[sel.lk])

# save selected objects for use in sampling programs
save(lkinfo, dfromx, zfromx, xfromz, botdepr, maxbotdep, eastr, northr, tsr, slopes, ints, d2shr.we, mid.d,
	file=paste(subdir, "/inputs-lake", sel.lk, ".Rdata", sep=""))



###  FISH  (targets)  ###

# check that exactly ONE of these two variables (WD, D2B) were set to NA
look <- spinfo[, c("WD", "D2B")]
not1na <- apply(is.na(look), 1, sum) != 1
if(sum(not1na)>0) stop(paste("Rows ", paste(seq_along(not1na)[not1na], collapse=", "), 
	".  Either a mean fishing depth or a mean distance to bottom MUST be specified, but NOT BOTH!", sep=""))
rm(look, not1na)

nrowz <- dim(spinfo)[1]
start.i <- (c(0, cumsum(spinfo$Nfish))+1)[1:nrowz]
end.i <- cumsum(spinfo$Nfish)
 
totfish <- sum(spinfo$Nfish)

attach(spinfo)

fish <- data.frame(sp=rep(G, Nfish), f.east=NA, f.north=NA, f.d2sh=NA, f.botdep=NA, f.fdep=NA, f.d2bot=NA, len=NA, wt=NA, ts=NA)
set.seed(lkinfo$Seed[sel.lk])
for(i in seq(nrowz)) {

	# easting available? no, then yes
	if(is.na(E[i])) {
		f.east <- runif(Nfish[i], eastr[1], eastr[2])
		f.d2sh <- dfromx(x=f.east, d2shr.we=d2shr.we, eastr=eastr, mid.d=mid.d)
		f.botdep <- zfromx(x=f.east, maxz=maxbotdep, eastr=eastr, ints=ints, slopes=slopes)
		} else {
		f.east <- rnorm(Nfish[i], 1000*E[i], 1000*EE[i]*E[i])
		f.d2sh <- dfromx(x=f.east, d2shr.we=d2shr.we, eastr=eastr, mid.d=mid.d)
		f.botdep <- zfromx(x=f.east, maxz=maxbotdep, eastr=eastr, ints=ints, slopes=slopes)
		}

	# northing available? no, then yes
	if(is.na(N[i])) {
		f.north <- runif(Nfish[i], northr[1], northr[2])
		} else {
		f.north <- rnorm(Nfish[i], 1000*N[i], 1000*NE[i]*N[i])
		}

	# distance to bottom available? no, then yes (if not, fishing depth is)
	if(is.na(D2B[i])) {
		f.fdep <- rnorm(Nfish[i], WD[i], WDE[i]*WD[i])
		f.d2bot <- f.botdep - f.fdep
		} else {
		f.d2bot <- rnorm(Nfish[i], D2B[i], D2BE[i]*D2B[i])
		f.fdep <- f.botdep - f.d2bot
		}

	# generate lengths from gamma distribution
	shape <- 1/ZE[i]^2
	scale <- Z[i]/shape
	len <- rgamma(Nfish[i], shape=shape, scale=scale)
	# predict weight from length-weight regression coefficients and add error
	wt. <- LWC1[i]*len^LWC2[i]
	wt <- wt. + rnorm(Nfish[i], 0, LWCE[i]*wt.)
	# predict target strength from length-ts regression coefficients and add error
	ts. <- TSC1[i] + TSC2[i]*log10(len/10)
	ts <- ts. + rnorm(Nfish[i], 0, TSCE[i]*abs(ts.))

	fish[start.i[i]:end.i[i], -1] <- cbind(f.east, f.north, f.d2sh, f.botdep, f.fdep, f.d2bot, len, wt, ts)
	}


rm(i, start.i, end.i, f.east, f.north, f.d2sh, f.botdep, f.fdep, f.d2bot, shape, scale, len, wt., wt, ts., ts)

bad.d2sh <- fish$f.d2sh < d2shr[1] | fish$f.d2sh > d2shr[2] | is.na(fish$f.d2sh)
bad.east <- fish$f.east < eastr[1] | fish$f.east > eastr[2]

bad.north <- fish$f.north < northr[1] | fish$f.north > northr[2]
bad.botdep <- fish$f.botdep < 0 | fish$f.botdep > maxbotdep
bad.fdep <- fish$f.fdep < 0 | fish$f.fdep > fish$f.botdep
bad.lenwt <- fish$len < 0 | fish$wt < 0
bad.ts <- fish$ts < tsr[1] | fish$ts > tsr[2]
n <- dim(fish)[1]

# if you are losing a lot of fish, you can use this to try and determine why
if(FALSE) {
sum(bad.d2sh)/n
sum(bad.east)/n
sum(bad.north)/n
sum(bad.botdep)/n
sum(bad.fdep)/n
sum(bad.lenwt)/n
sum(bad.ts)/n
}

# get rid of rows that are beyond the bounds we set or that have other problems
bad <- bad.d2sh | bad.east | bad.north | bad.botdep | bad.fdep | bad.lenwt | bad.ts
fish <- fish[!bad, ]

rm(bad.d2sh, bad.east, bad.north, bad.botdep, bad.fdep, bad.lenwt, bad.ts, bad)

save(fish, file=paste(subdir, "/fish-lake", sel.lk, ".Rdata", sep=""))

detach(spinfo)



###  diagnostic plots  ###



# a random selection of 1,000 fish (total)
n <- dim(fish)[1]
pick <- if(n<1001) 1:n else sample(1:n, 1000)
attach(fish[pick, ])
sus <- sort(unique(sp))

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(5.1, 4.1, 4.1, 2.1))
plot(f.east/1000, -f.fdep, type="n", xlim=eastr/1000, ylim=c(-maxbotdep, 0), 
	xlab="Easting  (km)", ylab="Water depth  (m)", main=paste(lkinfo$Lake[sel.lk], "- Side View"))
lines(c(0, xfromz(z=rep(maxbotdep-0.01, 2), maxz=maxbotdep, ints=ints, slopes=slopes, shore=0:1), eastr[c(2, 2, 1, 1)])/1000, 
	-c(botdepr[1], rep(maxbotdep, 2), botdepr[1], 0, 0, botdepr[1]))

for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	text(f.east[sel]/1000, -f.fdep[sel], sp[sel], col=i)
	}

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(5.1, 4.1, 4.1, 2.1))
plot(f.east/1000, f.north/1000, type="n", xlim=eastr/1000, ylim=northr/1000, 
	xlab="Easting  (km)", ylab="Northing  (km)", main=paste(lkinfo$Lake[sel.lk], "- Top View"))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	text(f.east[sel]/1000, f.north[sel]/1000, sp[sel], col=i)
	}

detach(fish[pick, ])



# a random selection of 250 fish FROM EACH SPECIES

rows.sp <- split(seq(along=fish$sp), fish$sp)
pick <- unlist(lapply(rows.sp, function(x) sample(x, min(250, length(x)))))
attach(fish[pick, ])
sus <- sort(unique(sp))

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	plot(f.east/1000, -f.fdep, type="n", xlim=eastr/1000, ylim=c(-maxbotdep, 0), xlab="", ylab="")
	lines(c(0, xfromz(z=rep(maxbotdep-0.01, 2), maxz=maxbotdep, ints=ints, slopes=slopes, shore=0:1), eastr[c(2, 2, 1, 1)])/1000, 
		-c(botdepr[1], rep(maxbotdep, 2), botdepr[1], 0, 0, botdepr[1]))
	text(f.east[sel]/1000, -f.fdep[sel], sp[sel], col=i)
	}
mtext("Easting  (km)", side=1, outer=TRUE)
mtext("Water depth  (m)", side=2, outer=TRUE)
mtext(paste(lkinfo$Lake[sel.lk], "- Side View"), side=3, outer=TRUE, font=2)

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	plot(f.east/1000, f.north/1000, type="n", xlim=eastr/1000, ylim=northr/1000, xlab="", ylab="")
	text(f.east[sel]/1000, f.north[sel]/1000, sp[sel], col=i)
	}
mtext("Easting  (km)", side=1, outer=TRUE)
mtext("Northing  (km)", side=2, outer=TRUE)
mtext(paste(lkinfo$Lake[sel.lk], "- Top View"), side=3, outer=TRUE, font=2)

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(5.1, 4.1, 4.1, 2.1))
plot(len, -f.fdep, type="n", xlab="Fish length  (mm)", ylab="Water depth  (m)", 
	main=paste(lkinfo$Lake[sel.lk], "- Size at Depth"))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	text(len[sel], -f.fdep[sel], sp[sel], col=i)
	}

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(5.1, 4.1, 4.1, 2.1))
plot(ts, -f.fdep, type="n", xlab="Target strength  (dB)", ylab="Water depth  (m)", 
	main=paste(lkinfo$Lake[sel.lk], "- Size at Depth"))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	text(ts[sel], -f.fdep[sel], sp[sel], col=i)
	}

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	plot(len[sel], wt[sel], xlab="", ylab="")
	mtext(sus[i], side=3, adj=0.1, line=-2, font=2)
	}
mtext("Length  (mm)", side=1, outer=TRUE)
mtext("Weight  (mm)", side=2, outer=TRUE)
mtext(paste(lkinfo$Lake[sel.lk], "- Length-Weight Relation"), side=3, outer=TRUE, font=2)

detach(fish[pick, ])



# histograms of all fish

attach(fish)

sus <- sort(unique(sp))

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	hist(len[sel], nclass=25, col="gray", xlab="", ylab="", main="")
	box()
	mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
	}
mtext("Length  (mm)", side=1, outer=TRUE)
mtext("Frequency", side=2, outer=TRUE)
mtext(paste(lkinfo$Lake[sel.lk], "- Length Distribution"), side=3, outer=TRUE, font=2)

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	hist(wt[sel], nclass=25, col="gray", xlab="", ylab="", main="")
	box()
	mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
	}
mtext("Weight  (g)", side=1, outer=TRUE)
mtext("Frequency", side=2, outer=TRUE)
mtext(paste(lkinfo$Lake[sel.lk], "- Weight Distribution"), side=3, outer=TRUE, font=2)

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	hist(ts[sel], nclass=25, col="gray", xlab="", ylab="", main="")
	box()
	mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
	}
mtext("Target strength  (dB)", side=1, outer=TRUE)
mtext("Frequency", side=2, outer=TRUE)
mtext(paste(lkinfo$Lake[sel.lk], "- TS Distribution"), side=3, outer=TRUE, font=2)

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	hist(f.east[sel]/1000, nclass=25, col="gray", xlim=eastr/1000, xlab="", ylab="", main="")
	box()
	mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
	}
mtext("Easting  (km)", side=1, outer=TRUE)
mtext("Frequency", side=2, outer=TRUE)
mtext(paste(lkinfo$Lake[sel.lk], "- Easting Distribution"), side=3, outer=TRUE, font=2)

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	hist(f.north[sel]/1000, nclass=25, col="gray", xlim=northr/1000, xlab="", ylab="", main="")
	box()
	mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
	}
mtext("Northing  (km)", side=1, outer=TRUE)
mtext("Frequency", side=2, outer=TRUE)
mtext(paste(lkinfo$Lake[sel.lk], "- Northing Distribution"), side=3, outer=TRUE, font=2)

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	hist(f.d2sh[sel], nclass=25, col="gray", xlim=d2shr, xlab="", ylab="", main="")
	box()
	mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
	}
mtext("Distance to Shore  (m)", side=1, outer=TRUE)
mtext("Frequency", side=2, outer=TRUE)
mtext(paste(lkinfo$Lake[sel.lk], "- Distance to Shore Distribution"), side=3, outer=TRUE, font=2)

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	hist(f.fdep[sel], nclass=25, col="gray", xlim=c(0, maxbotdep), xlab="", ylab="", main="")
	box()
	mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
	}
mtext("Fishing Depth  (m)", side=1, outer=TRUE)
mtext("Frequency", side=2, outer=TRUE)
mtext(paste(lkinfo$Lake[sel.lk], "- Fishing Depth Distribution"), side=3, outer=TRUE, font=2)

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	hist(f.d2bot[sel], nclass=25, col="gray", xlim=c(0, maxbotdep), xlab="", ylab="", main="")
	box()
	mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
	}
mtext("Distance to Bottom  (m)", side=1, outer=TRUE)
mtext("Frequency", side=2, outer=TRUE)
mtext(paste(lkinfo$Lake[sel.lk], "- Distance to Bottom Distribution"), side=3, outer=TRUE, font=2)

if(!save.plots) windows(w=9, h=6.5)
par(mfrow=n2mfrow(length(sus)), oma=c(2, 2, 2, 0), mar=c(3, 3, 1, 1))
for(i in seq(along=sus)) {
	sel <- sp==sus[i]
	hist(f.botdep[sel], nclass=25, col="gray", xlim=c(0, maxbotdep), xlab="", ylab="", main="")
	box()
	mtext(sus[i], side=3, adj=0.9, line=-2, font=2)
	}
mtext("Bottom Depth  (m)", side=1, outer=TRUE)
mtext("Frequency", side=2, outer=TRUE)
mtext(paste(lkinfo$Lake[sel.lk], "- Bottom Depth Distribution"), side=3, outer=TRUE, font=2)

detach(fish)



if(save.plots) graphics.off()



# total number and weight of each species in population
truth <- cbind(n=table(fish$sp), kg=tapply(fish$wt, fish$sp, sum)/1000)
truth <- as.data.frame(rbind(truth, Total=apply(truth, 2, sum, na.rm=TRUE)))

truth$dens <- truth$n/lkinfo$LkArea[lkinfo$LC==sel.lk]
truth$bio <- truth$kg/lkinfo$LkArea[lkinfo$LC==sel.lk]

# save the true numbers and weights of each species in the population to a csv file
write.csv(truth, file=paste(subdir, "/Truth-lake", sel.lk, ".csv", sep=""))

# print out of summary numbers
truth

# remove all objects from current working directory
rm(botdepr, cap.no.fish, d2shr, d2shr.we, dfromx, dir, eastr, fish, i, ints, lkinfo, maxbotdep, mid.d, n, northr, nrowz, 
	pick, rows.sp, save.plots, sel, sel.lk, slopes, spinfo, subdir, sus, totfish, TotNFish, truth, tsr, xfromz, zfromx, wb)
