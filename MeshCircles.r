# determining the size of fish that can pass through a given mesh
# modification of MATLAB code from Kresimir Williams, NOAA-AFSC, 23 April 2014

# no inputs required

##########################################################################
############################  data crunching  ############################
##########################################################################

# load the plotrix, MASS, and rpanel packages
if(!require(plotrix)) stop("R package plotrix is required")
if(!require(MASS)) stop("R package MASS is required")
if(!require(rpanel)) stop("R package rpanel is required")

bigcircle <- function(barmesh, pctopen) {
	# find "max" body depth of fish that can fit through the specified mesh in a left-right orientation
	# barmesh = the length of one side of a square mesh in inches
	# pctopen = the ratio of the height of a partially open mesh to its width (should be between 0 and 1)
	# half width of partially open mesh
	W <- barmesh / sqrt(1 + pctopen^2)
	# half height of partially open mesh
	H <- W * pctopen
	# half angle of mesh peak (should be > 90 degrees, or pi/2)
	theta <- atan(1/pctopen)
	# radius of biggest circle that can fit inside right (or upper) half of partially open mesh
	r <- H * tan(theta/2)
	# fish body depth in inches
	bodydepth <- 4*r
	list(barmesh=barmesh, pctopen=pctopen, W=W, H=H, r=r, bodydepth=bodydepth)
	}

addfish <- function(radius, leftright=TRUE) {
	# draw a "fish" on the mesh ... consisting of two tangent circles, with an overlaying square
	draw.ellipse(x=radius*leftright, y=radius*!leftright, a=radius, b=radius, angle=0, col="gray", border=NA)
	draw.ellipse(x=-radius*leftright, y=-radius*!leftright, a=radius, b=radius, angle=0, col="gray", border=NA)
	polygon(c(-radius, radius, radius, -radius), c(-radius, -radius, radius, radius), col="gray", border=NA)
	}

drawmesh <- function(panel) {
	# draw a single mesh with the biggest fish that can fit through, both horizontally and vertically
	horiz <- bigcircle(panel$barmesh, panel$pctopen)
	vert <- bigcircle(panel$barmesh, 1/panel$pctopen)
	both <- c(horiz, r2=vert$r, bodydepth2=vert$bodydepth)
	attach(both)
	flength <- round(bodydepth*panel$ldratio*25.4)
	flength2 <- round(bodydepth2*panel$ldratio*25.4)
	eqscplot(0, 0, type="n", xlim=c(-W, W), ylim=c(-H, H), las=1, xlab="", ylab="", 
		main=paste0('Horizontal\nbody depth = ', round(bodydepth, 2), '", length = ', round(flength), ' mm'))
	addfish(r)
	polygon(c(0, W, 0, -W), c(-H, 0, H, 0), lwd=2, density=0)
	eqscplot(0, 0, type="n", xlim=c(-W, W), ylim=c(-H, H), las=1, xlab="", ylab="", 
		main=paste0('Vertical\nbody depth = ', round(bodydepth2, 2), '", length = ', round(flength2), ' mm'))
	addfish(r2, FALSE)
	polygon(c(0, W, 0, -W), c(-H, 0, H, 0), lwd=2, density=0)
	mtext(paste0(100*pctopen, '% open, bar mesh = ', barmesh, '", len/bdy dep = ', panel$ldratio), side=3, outer=TRUE, cex=1.2, line=0.5)
	mtext("Inches", side=1, outer=TRUE)
	mtext("Inches", side=2, outer=TRUE)
	detach(both)
	panel
	}

# paste0(names(both), "=both$", names(both), collapse=", ")

windows(h=9, w=6.5)
par(mfrow=c(2, 1), mar=c(3, 3, 3, 1), oma=c(1, 1, 2, 0))
panel <- rp.control(barmesh=1.5, pctopen=0.3, ldratio=5)
rp.slider(panel, pctopen, 0, 1, resolution=0.01, showvalue=T, action=drawmesh, 
	title="Relative opening of mesh (height/width)                                                                   ")
rp.slider(panel, barmesh, 0.25, 10, resolution=0.01, showvalue=T, action=drawmesh, title="Bar mesh (inches)")
rp.slider(panel, ldratio, 1, 20, resolution=0.1, showvalue=T, action=drawmesh, title="Fish length to body depth ratio")
