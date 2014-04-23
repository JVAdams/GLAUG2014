# Use a slider to explore availability using double logistic function

# load the rpanel package
if(!require(rpanel)) stop("R package rpanel is required")

dbllogit2 <- function(fdep, D501, SD1, d2bot, D502, SD2) {
	# double logistic regression
	# allowing for parabolic-shaped probabilities
	# function provided by Kresimir Williams, NOAA-AFSC, 22 April 2014
	((1 + exp(2*log(3)*(D501 - fdep)/SD1))^-1) * ((1 + exp(2*log(3)*(D502 - d2bot)/SD2))^-1)
	}

# probability graphing function
double.draw2 <- function(panel) {
	d2bot <- 100 - panel$x
	y <- dbllogit2(fdep=panel$x, D501=panel$D501, SD1=panel$SD1, d2bot=d2bot, D502=panel$D502, SD2=panel$SD2)
	plot(panel$x, y, las=1, xlim=c(0, 100), ylim=0:1, axes=FALSE, type="n", xlab="Distance from bottom  (m)", ylab="Availability")
	mtext("Distance from surface  (m)", side=3, line=3)
	axis(1, at=100-pretty(d2bot), labels=pretty(d2bot))
	axis(2)
	axis(3)
	box()
	abline(v=c(panel$D501, 100-panel$D502), col="gray", lwd=2)
	abline(h=c(0, 0.5, 1), col="gray", lwd=2)
	lines(spline(panel$x, y, 1000), lwd=3)
	panel
	}

# plot it, with a slider to adjust coeficients of the double logistic function
windows(w=7, h=5)
par(mar=c(4, 4, 4, 1))
panel <- rp.control(x=seq(0, 100, 0.1), D501=5.1, SD1=10.1, D502=5.1, SD2=20.1)
rp.slider(panel, D501, 0.1, 40, resolution=0.1, showvalue=T, action=double.draw2, title="D501")
rp.slider(panel, SD1, 0.1, 200, resolution=0.1, showvalue=T, action=double.draw2, title="Dslope1")
rp.slider(panel, D502, 0.1, 40, resolution=0.1, showvalue=T, action=double.draw2, title="D502")
rp.slider(panel, SD2, 0.1, 200, resolution=0.1, showvalue=T, action=double.draw2, title="Dslope2")
