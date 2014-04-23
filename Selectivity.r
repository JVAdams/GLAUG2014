# Use a slider to explore selectivity using double logistic function

# load the rpanel package
if(!require(rpanel)) stop("R package rpanel is required")

dbllogit <- function(len, L501, SR1, L502, SR2) {
	# double logistic regression
	# allowing for parabolic-shaped probabilities
	# function provided by Kresimir Williams, NOAA-AFSC, 22 April 2014
	((1 + exp(2*log(3)*(L501 - len)/SR1))^-1) * ((1 + exp(2*log(3)*(L502 - len)/-SR2))^-1)
	}

# probability graphing function
double.draw <- function(panel) {
	y <- dbllogit(len=panel$x, L501=panel$L501, SR1=panel$SR1, L502=panel$L502, SR2=panel$SR2)
	plot(panel$x, y, las=1, xlim=range(panel$x[y>0.001]), ylim=0:1, type="n", xlab="Length  (mm)", ylab="Selectivity")
	abline(v=c(panel$L501, panel$L502), col="gray", lwd=2)
	abline(h=c(0, 0.5, 1), col="gray", lwd=2)
	lines(spline(panel$x, y, 1000), lwd=3)
	panel
	}

# plot it, with a slider to adjust coefficients of the double logistic function
windows(w=7, h=5)
par(mar=c(4, 4, 1, 1))
panel <- rp.control(x=1:2000, L501=10.1, SR1=10.1, L502=200.1, SR2=20.1)
rp.slider(panel, L501, 0.1, 1000, resolution=0.1, showvalue=T, action=double.draw, title="L501")
rp.slider(panel, SR1, 0.1, 200, resolution=0.1, showvalue=T, action=double.draw, title="Lslope1")
rp.slider(panel, L502, 0.1, 1000, resolution=0.1, showvalue=T, action=double.draw, title="L502")
rp.slider(panel, SR2, 0.1, 200, resolution=0.1, showvalue=T, action=double.draw, title="Lslope2")
