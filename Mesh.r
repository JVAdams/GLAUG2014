# determining the size of fish that can pass through a given mesh
# modification of MATLAB code from Kresimir Williams, NOAA-AFSC, 23 April 2014

##########################################################################
#############################  input values  #############################
##########################################################################

# the extent to which the mesh is stretched while fishing, expressed as height / width
pctopen <- 0.2

# bar mesh in inches (stretch mesh = bar mesh * 2)
barmesh <- 1.5

##########################################################################
############################  data crunching  ############################
##########################################################################

# load the plotrix, MASS, and rpanel packages
if(!require(plotrix)) stop("R package plotrix is required")
if(!require(MASS)) stop("R package MASS is required")
if(!require(rpanel)) stop("R package rpanel is required")

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

# fish length to body depth ratio
ldratio <- 4.8

# fish body depth in inches
bodydepth <- 1

strmesh <- 2 * barmesh

# mesh and ellipse graphing function
drawmesh <- function(panel) {
  flength <- round(panel$bodydepth*panel$ldratio*25.4)
  L <- sqrt(4*panel$barmesh^2 / (panel$pctopen^2 + 1))
  H <- L * panel$pctopen
  eqscplot(c(0, L, 0, -L, 0)/2, c(-H, 0, H, 0, -H)/2, type="l", las=1, xlab="Inches", ylab="Inches", 
    main=paste0('mesh opening = ', 100*panel$pctopen, '%, bar mesh = ', panel$barmesh, 
    '",\nfish len/dpth ratio = ', panel$ldratio, ', body dpth = ', panel$bodydepth, '", len = ', flength, ' mm'), lwd=2)
  draw.ellipse(x=0, y=0, a=panel$bodydepth/2, b=panel$bodydepth/4, angle=0, col="gray", lwd=2)
  panel
  }

# plot it, with a slider to adjust coefficients of the double logistic function
windows()
par(mar=c(4, 4, 3, 1))
panel <- rp.control(bodydepth=bodydepth, ldratio=ldratio, barmesh=barmesh, pctopen=pctopen)
rp.slider(panel, bodydepth, 0.01, 5, resolution=0.01, showvalue=T, action=drawmesh, title="Body depth (inches)")
rp.slider(panel, ldratio, 1, 20, resolution=0.1, showvalue=T, action=drawmesh, title="Fish length to depth ratio")
