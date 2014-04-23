# determining the size of fish that can pass through a given mesh
# modification of MATLAB code from Kresimir Williams, NOAA-AFSC, 23 April 2014

##########################################################################
#############################  input values  #############################
##########################################################################

# the extent to which the mesh is stretched while fishing, expressed as height as a percentage of width
pctopen <- 0.35

# mesh in inches
mesh <- 2.5

# fish length to body depth ratio
ldratio <- 8

# fish body depth
bodydepth <- 1.25

##########################################################################
############################  data crunching  ############################
##########################################################################

# load the plotrix package
if(!require(plotrix)) stop("R package plotrix is required")

polar2cart <- function(r, theta, degrees=FALSE) {
	# http://stackoverflow.com/questions/16351178/r-converting-cartesian-coordinates-to-polar-coordinates-and-then-calculating-d
	# convert degrees to radians (dividing by 360/2*pi, or multiplying by pi/180)
	if(degrees) {
		rads <- theta*pi/180
		} else {
		rads <- theta
		}
	# convert to cartesian coordinates
	x <- r * sin(rads)
	y <- r * cos(rads)
	c(x=x, y=y)
	}

length <- bodydepth*ldratio*25.4

xy <- polar2cart(mesh, asin(mesh/mesh*pctopen))
x <- xy["x"]
y <- xy["y"]

windows()
par(mar=c(3, 3, 3, 1))
eqscplot(c(0, y, 0, -y, 0), c(-x, 0, x, 0, -x), type="l", las=1, xlab="", ylab="", 
	main=paste0('pctopen = ', 100*pctopen, '%, mesh = ', mesh, 
	'",\nldratio = ', ldratio, ', bodydepth = ', bodydepth, '", length = ', length, ' mm'), lwd=2)
draw.ellipse(x=0, y=0, a=bodydepth, b=bodydepth/2, angle=0, col="gray", lwd=2)
