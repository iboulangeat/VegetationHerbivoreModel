colo <-function(alpha = 255){
  co = c(T = rgb(254, 127, 0, alpha, maxColorValue = 255), V= rgb(12, 12, 12, alpha, maxColorValue = 255) , R= rgb(32, 201, 69, alpha, maxColorValue = 255), H = rgb(0, 0, 0, alpha, maxColorValue = 255), B = rgb(38, 140, 248, alpha, maxColorValue = 255))
  return(co)
}
##### graph  template  #######
plot.template <- function(xlim, ylim, xticks, yticks, xlab, ylab, col.grid = "white", col.bg = "gray", plot.axis = TRUE, ...)
{
  ## Empty plot
  par(family = "serif", col.axis = 1)
  plot(xlim, ylim, type = "n", axes = F, ann = F, ...)
  ## Background
  rect(xlim[1], ylim[1], xlim[2], ylim[2], col = col.bg, border = par()$col.axis)
  for (i in seq(xlim[1], xlim[2], length.out = xticks[2])) segments(i,ylim[1], i, ylim[2], col = col.grid , lty = 3)
  for (i in seq(xlim[1], xlim[2], length.out = xticks[1])) segments(i,ylim[1], i, ylim[2],col = col.grid )
  for (i in seq(ylim[1], ylim[2], length.out = yticks[2])) segments(xlim[1], i, xlim[2],i,col = col.grid , lty = 3)
  for (i in seq(ylim[1], ylim[2], length.out = yticks[1])) segments(xlim[1], i, xlim[2],i,col = col.grid )
  if(plot.axis)
  {
    axis(1, pos = 0, at = seq(xlim[1], xlim[2], length.out = xticks[1]), col = par()$col.axis, labels = format(seq(xlim[1], xlim[2], length.out = xticks[1])), font = 2)
    axis(1, pos = 0, at = seq(xlim[1], xlim[2], length.out = xticks[2]), labels = F, lwd = 0, tck = -0.01, lwd.ticks = 1, col.ticks = par()$col.axis)
    axis(side = 2, pos = xlim[1], at = seq(ylim[1], ylim[2], length.out = yticks[1]), col = par()$col.axis, labels = format(seq(ylim[1], ylim[2], length.out = yticks[1])), font = 2, las = 2)
    axis(side = 2, pos = xlim[1], at = seq(ylim[1], ylim[2], length.out = yticks[2]), labels = F, lwd = 0, tck = -0.01, lwd.ticks = 1, col.ticks = par()$col.axis)
  }

  mtext(xlab, 1, line = 1.5, font = 2, col = par()$col.axis, cex = par()$cex)
  mtext(ylab, 2, line = 1.5, font = 2, col = par()$col.axis, cex = par()$cex)
}

plot.eq <- function(newEq.all, newEq.veg, O2B, B2T, ymin, ymax)
{
  par(mfrow = c(3,1), cex = 1.2, mar = c(3,4,3,1))

  plot.template(xlim = c(0,50), ylim = c(0,1), xticks = c(6,11), yticks = c(6,11), xlab = "increasing temperature", ylab = "Vegetation proportion", col.bg = "white", col.grid = NA, plot.axis = F)
  lines(newEq.veg[,"T"]~gradient, col = colo()["T"], type = "l", lwd = 1.5)
  lines(newEq.veg[,"B"]~gradient, col = colo()["B"], type = "l", lwd = 1.5)
  lines(newEq.veg[,"S"]~gradient, col = colo()["S"], type = "l", lwd = 1.5)
  lines(newEq.veg[,"G"]~gradient, col = colo()["G"], type = "l", lwd = 1.5)
  axis(1, pos=0, at = seq(0,50,le=6), labels = seq(-4,6,2), cex.axis=.9)
  axis(2, pos=0, at = seq(0,1,.2), labels = seq(0,1,.2), cex.axis=.9, las = 2)
  mtext(at=10, side = 3, line = 0, text = "boreal",adj=0,las=1, cex=1)
  mtext(at=35, side = 3, line = 0, text = "temperate",adj=0,las=1, cex=1)
  title("(a)")

  plot.template(xlim = c(0,50), ylim = c(0,1), xticks = c(6,11), yticks = c(6,11), xlab = "increasing temperature", ylab = "Vegetation proportion", col.bg = "white", col.grid = NA, plot.axis = F)
  lines(newEq.all[,"T"]~gradient, col = colo()["T"], type = "l", lwd = 1.5)
  lines(newEq.all[,"B"]~gradient, col = colo()["B"], type = "l", lwd = 1.5)
  lines(newEq.all[,"S"]~gradient, col = colo()["S"], type = "l", lwd = 1.5)
  lines(newEq.all[,"G"]~gradient, col = colo()["G"], type = "l", lwd = 1.5)
  axis(1, pos=0, at = seq(0,50,le=6), labels = seq(-4,6,2), cex.axis=.9)
  axis(2, pos=0, at = seq(0,1,.2), labels = seq(0,1,.2), cex.axis=.9, las =2)

  rect(0,0,O2B-2,1,col=colo(70)["G"], border = NA)
  rect(O2B-1,0,B2T-2,1,col=colo(70)["B"], border = NA)
  rect(B2T-1,0,50,1,col=colo(70)["T"], border = NA)
  mtext(at=0, side = 3, line = 0, text = "open",adj=0,las=1, cex=1)
  mtext(at=13, side = 3, line = 0, text = "boreal",adj=0,las=1, cex=1)
  mtext(at=37, side = 3, line = 0, text = "temperate",adj=0,las=1, cex=1)
  points(O2B-1, 1, pch=19)
  points(B2T-1, 1, pch=19)
  title("(b)")

  par(cex=1.2)
  plot.template(xlim = c(0,50), ylim = c(ymin,ymax), xticks = c(6,11), yticks = c(6,11), xlab = "increasing temperature", ylab = "herbivore biomass (tons per km2)", col.bg = NA, col.grid = NA,plot.axis=FALSE)
  axis(1, pos=0, at = seq(0,50,le=6), labels = seq(-4,6,2), cex.axis=.9)
  axis(2, pos = 0, at = seq(0,2,0.5), labels = seq(0,2,0.5))
  lines(newEq.all[,"H"]/1000~gradient, col = 1, type = "l", lwd = 1.5)
  lines(newEq.veg[,"H"]/1000~gradient, col = 1, lwd = 1.5, lty = 2)
  points(O2B-1, 2, pch=19)
  points(B2T-1, 2, pch=19)
  title("(c)")#

}

#=====================
# plot trajectory to equilibrium after solving equilibrium
#=====================
plotEq <- function(eq)
{
  trace.mat = eq$mat
  plot.template(xlim = c(0, nrow(trace.mat)), ylim = c(0,1), xticks = c(6,11), yticks = c(6,11), xlab = "time", ylab = "vegetation proportions", col.bg = "white", col.grid = NA)
  points(trace.mat[,1], col=colo()["T"], cex=.1)
  points(trace.mat[,2], col=colo()["S"], cex=.1)
  points(trace.mat[,3], col=colo()["B"], cex=.1)
  par(new = TRUE)
  plot( c(0, nrow(trace.mat)), c(0,max(trace.mat[,4])), type = "n", axes = F, ann = F)
  lines(trace.mat[,4], col=colo()["H"], lwd=.8)
}

## title plot
ltitle=function(x,backcolor="#336600",forecolor="#ccff99",cex=2,ypos=0.4, xpos=0.5, srt=0)
{
  plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE)
  polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col=backcolor,border=NA)
  text(x=xpos,y=ypos,pos=4,cex=cex,labels=x,col=forecolor, srt=srt, font=2)
}
