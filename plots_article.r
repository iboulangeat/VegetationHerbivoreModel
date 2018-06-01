source("plot_fct.r")
##--- calibration
plot_impacts <- function()
{
  par(mfrow = c(1,2))
  plot(impact~pressure, type = "l", xlab = "Browser impact on S", ylab = "impact multiplicative factor", xlim = c(0,1), ylim = c(0,1))
  abline(v = h, lty=2)
  abline(h=0.5, lty = 3)
  text(h+0.05, 1,"h", cex=.8)
  plot(impact~c(0,H,Inf), type = "l", xlab = "Browser density", ylab = "impact multiplicative factor", ylim = c(0,1), xaxt = "n")
  axis(1, at = seq(0.5,15,l=6)*37500, labels = seq(0.5,15,l=6))
}

plot_fitnesses <- function()
{
  par(mar = c(5, 3, 3, 1), mfrow = c(1,2), cex = 1.2)
  plot(fit[["S"]]~gradient, type = "n", col = colo()["S"], ylim = c(0,.3),yaxt = "n", xlab = "annual average temperature", xaxt="n")
  axis(1, at = seq(0,50,le=6), labels = temp.range)
  title("(a)", font = 2)
  lines(fit[["B"]]~gradient, col = colo()["B"])
  lines(fit[["T"]]~gradient, col = colo()["T"])
  legend(x = min(gradient), y = .3, bg = NA, box.col=NA,  col = colo()[c("B", "T")], legend = names(colo()[c("B", "T")]), lwd = 1.5, ncol = 5, xjust = 0, yjust = 1, cex = 1)
  mtext("fitness", side =2, line = 1)

  par(mar = c(5, 5, 3, 1))
  plot(c0.var~gradient, type = "l", col = 1, ylim = c(0,1), ylab = "parameter values", xaxt="n", xlab = "annual average temperature")
  title("(b)", font = 2)
  axis(1, at = seq(0,50,le=6), labels = temp.range)
  lines(k.var~gradient, col = 1, lty=2)
  legend(x = min(gradient), y = 1, bg = NA, box.col=NA,  col = 1,  lty = 1:2,legend = c("c","k"), lwd = 1.5, ncol = 1, xjust = 0, yjust = 1, cex = 1.1)
}
##--- sensitivity analysis
plot_sa.a <- function()
{
  par(cex = 2)
  varImpPlot(mod1, main = "",
	labels = c("c","a",expression(d[T]),expression(d[B]),expression(d[T]/d[B]), "k"))
  title("(a)", font = 2)
}

plot_sa.b <- function()
{
  par(cex = 2)
  varImpPlot(mod2, main = "",
	labels = c("a","k",expression(d[T]/d[B]),expression(d[T]),expression(d[B]),"c"))
  title("(b)", font = 2)
}
##--- supp gradients
plot_suppGradients <- function()
{
  colos = c(1,"darkgreen","orange","blue", "red")
  layout(matrix(1:8, nrow = 4, byrow = F), width = c(1,1), height=c(1,1,1,.3))

  par(cex=1.2, mar=c(1,4,3,1))
  plot.template(xlim = c(0,50), ylim = c(0,1), xticks = c(6,11), yticks = c(6,11), xlab = "", ylab = "Proportion of mature tree stands", col.bg = NA, col.grid = NA,plot.axis=FALSE)
  axis(1, pos=0, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos = 0, at = seq(0,1,0.2), labels = seq(0,1,0.2))
  forest = analysis.energy$newEq.all[,"T"] + analysis.energy$newEq.all[,"B"]
  lines(forest~gradient, col = colos[3], type = "l", lwd = 1.5)
  forest = analysis.herbivory$newEq.all[,"T"] + analysis.herbivory$newEq.all[,"B"]
  lines(forest~gradient, col = colos[4], type = "l", lwd = 1.5)
  forest = analysis.season$newEq.all[,"T"] + analysis.season$newEq.all[,"B"]
  lines(forest~gradient, col = colos[5], type = "l", lwd = 1.5)
  forest = analysis.base$newEq.all[,"T"] + analysis.base$newEq.all[,"B"]
  lines(forest~gradient, col = colos[1], type = "l", lwd = 1.5)
  mtext("(a)", 3, font = 2, cex=1.5)

  par(cex=1.2, mar=c(1,4,3,1))
  plot.template(xlim = c(0,50), ylim = c(0,0.01), xlab = "", ylab = "Asymptotic resilience", xticks = c(6,11), yticks = c(6,11), col.bg = NA, col.grid = NA, plot.axis=FALSE)
  axis(1, pos=0, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos = 0, at = seq(0,0.01,0.002), labels = seq(0,0.01,0.002))
  lines(-analysis.energy$lambdamax.all~gradient, col = colos[3], type = "l", lwd = 1.5)
  lines(-analysis.herbivory$lambdamax.all~gradient, col = colos[4], type = "l", lwd = 1.5)
  lines(-analysis.season$lambdamax.all~gradient, col = colos[5], type = "l", lwd = 1.5)
  lines(-analysis.base$lambdamax.all~gradient, col = colos[1], type = "l", lwd = 1.5)
  mtext("(b)", 3, font = 2, cex=1.5)

  par(cex=1.2, mar=c(1,4,3,1))
  plot.template(xlim = c(0,50), ylim = c(-0.3,0), xlab = "", ylab = "Initial resilience", xticks = c(6,11), yticks = c(6,11), col.bg = NA, col.grid = NA, plot.axis=FALSE)
  axis(1, pos=-0.3, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos = 0, at = seq(-0.3,0,0.05), labels = seq(-0.3,0,0.05))
  lines(-analysis.energy$reac.all~gradient, col = colos[3], type = "l", lwd = 1.5)
  lines(-analysis.herbivory$reac.all~gradient, col = colos[4], type = "l", lwd = 1.5)
  lines(-analysis.season$reac.all~gradient, col = colos[5], type = "l", lwd = 1.5)
  lines(-analysis.base$reac.all~gradient, col = colos[1], type = "l", lwd = 1.5)
  mtext("(c)", 3, font = 2, cex=1.5)

  par(cex=1.2, mar=c(1,4,1,1))
  plot(1:10,1:10, ylab="", xlab="", xaxt="n", yaxt="n", type="n", bty="n")
  mtext("Temperature", 1, line = -1.5, font = 2, col = 1, cex = 1.2)

  par(cex=1.2, mar=c(1,4,3,1))
  plot.template(xlim = c(0,50), ylim = c(0,0.15), xlab = "", ylab = "Exposure", xticks = c(6,11), yticks = c(6,11), col.bg = NA, col.grid = NA, plot.axis=FALSE)
  axis(1, pos=0, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos = 0, at = seq(0,0.15,0.05), labels = seq(0,0.15,0.05))
  lines(c(NA,analysis.energy$deltaN.all)~gradient, col = colos[3], type = "l", lwd = 1.5)
  lines(c(NA,analysis.herbivory$deltaN.all)~gradient, col = colos[4], type = "l", lwd = 1.5)
  lines(c(NA,analysis.season$deltaN.all)~gradient, col = colos[5], type = "l", lwd = 1.5)
  lines(c(NA,analysis.base$deltaN.all)~gradient, col = colos[1], type = "l", lwd = 1.5)
  mtext("(d)", 3, font = 2, cex=1.5)

  par(cex=1.2, mar=c(1,4,3,1))
  plot.template(xlim = c(0,50), ylim = c(0,8000), xlab = "", ylab = "Sensitivity", xticks = c(6,11), yticks = c(6,11), col.bg = NA, col.grid = NA, plot.axis=FALSE)
  axis(1, pos=0, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos = 0, at = seq(0,8000,2000), labels = seq(0,8000,2000))
  lines(c(NA,analysis.energy$deltaT.all)~gradient, col = colos[3], type = "l", lwd = 1.5)
  lines(c(NA,analysis.herbivory$deltaT.all)~gradient, col = colos[4], type = "l", lwd = 1.5)
  lines(c(NA,analysis.season$deltaT.all)~gradient, col = colos[5], type = "l", lwd = 1.5)
  lines(c(NA,analysis.base$deltaT.all)~gradient, col = colos[1], type = "l", lwd = 1.5)
  mtext("(e)", 3, font = 2, cex=1.5)

  par(cex=1.2, mar=c(1,4,3,1))
  plot.template(xlim = c(0,50), ylim = c(0,35), xlab = "", ylab = "Cumulative changes in state", xticks = c(6,11), yticks = c(6,11), col.bg = NA, col.grid = NA, plot.axis=FALSE)
  axis(1, pos=0, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos = 0, at = seq(0,35,5), labels = seq(0,35,5))
  lines(c(NA,analysis.energy$integraleV.all)~gradient, col = colos[3], type = "l", lwd = 1.5)
  lines(c(NA,analysis.herbivory$integraleV.all)~gradient, col = colos[4], type = "l", lwd = 1.5)
  lines(c(NA,analysis.season$integraleV.all)~gradient, col = colos[5], type = "l", lwd = 1.5)
  lines(c(NA,analysis.base$integraleV.all)~gradient, col = colos[1], type = "l", lwd = 1.5)
  mtext("(f)", 3, font = 2, cex=1.5)

  par(cex=1.2, mar=c(1,4,1,1))
  plot(1:10,1:10, ylab="", xlab="", xaxt="n", yaxt="n", type="n", bty="n")
  mtext("Temperature \n(final temperature)", 1, line = -1, font = 2, col = 1, cex = 1.2)
}
##--- Fisichelli
plot_fisichelli <- function()
{
  par(mfrow=c(1,2))
  colo = c("blue", "orange")
  plot(yy,fitT, type="l", col= colo[1], ylab="theoretical fitness value",
                        xlab = "mean annual temperature")
  lines(yy, fitB, col=colo[2], type="l")
  legend("right", col=colo, lwd=1, legend=c("T", "B"), bty="n")
  plot(yy,seq(6,10.5,le=length(fitT)), type="l", col= colo[1], ylab="approx. heigth growth",
                        xlab = "mean annual temperature", ylim = c(6, 10))
  lines(yy, seq(9,6,le=length(fitT)), col=colo[2])
  legend("right", col=colo, lwd=1,
                legend=c("Temperate", "Boreal"), bty="n")
  title("Fisichelli")
}
##--- main
plot_equilibrium <- function(ylim.H=3)
{
  par(mfrow = c(3,1), cex = 1.2, mar = c(3,4,3,1))
  plot.template(xlim = c(0,50), ylim = c(0,1), xticks = c(6,11), yticks = c(6,11), xlab = "Temperature", ylab = "Vegetation proportion", col.bg = "white", col.grid = NA, plot.axis = F)
  lines(newEq.veg[,"T"]~gradient, col = colo()["T"], type = "l", lwd = 1.5)
  lines(newEq.veg[,"B"]~gradient, col = colo()["B"], type = "l", lwd = 1.5)
  lines(newEq.veg[,"R"]~gradient, col = colo()["R"], type = "l", lwd = 1.5)
  lines(newEq.veg[,"V"]~gradient, col = colo()["V"], type = "l", lwd = 1.5)
  axis(1, pos=0, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos=0, at = seq(0,1,.2), labels = seq(0,1,.2), cex.axis=.9, las = 2)
  points(B2T.veg-1, 1, pch=1)
  title("(a)")
  plot.template(xlim = c(0,50), ylim = c(0,1), xticks = c(6,11), yticks = c(6,11), xlab = "Temperature", ylab = "Vegetation proportion", col.bg = "white", col.grid = NA, plot.axis = F)
  lines(newEq.all[,"T"]~gradient, col = colo()["T"], type = "l", lwd = 1.5)
  lines(newEq.all[,"B"]~gradient, col = colo()["B"], type = "l", lwd = 1.5)
  lines(newEq.all[,"R"]~gradient, col = colo()["R"], type = "l", lwd = 1.5)
  lines(newEq.all[,"V"]~gradient, col = colo()["V"], type = "l", lwd = 1.5)
  axis(1, pos=0, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos=0, at = seq(0,1,.2), labels = seq(0,1,.2), cex.axis=.9, las =2)
  points(Hmax-1, 1, pch=8)
  points(O2B-1, 1, pch=19)
  points(B2T-1, 1, pch=19)
  title("(b)")
  par(cex=1.2)
  plot.template(xlim = c(0,50), ylim = c(0,ylim.H), xticks = c(6,11), yticks = c(6,11), xlab = "Temperature", ylab = "browser population biomass (tons per km2)", col.bg = NA, col.grid = NA,plot.axis=FALSE)
  axis(1, pos=0, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos = 0, at = seq(0,ylim.H,0.5), labels = seq(0,ylim.H,0.5))
  lines(newEq.all[,"H"]/1000~gradient, col = 1, type = "l", lwd = 1.5)
  lines(newEq.veg[,"H"]/1000~gradient, col = 1, lwd = 1.5, lty = 2)
  points(Hmax-1, ylim.H, pch=8)
  points(O2B-1, ylim.H, pch=19)
  points(B2T-1, ylim.H, pch=19)
  title("(c)")
}
##---
plot_transient <- function()
{
  layout(matrix(1:8, nrow = 4, byrow = F), width = c(1,1), height=c(1,1,1,.3))

  par(cex=1.4, mar=c(1,4,3,1))
  plot.template(xlim = c(0,50), ylim = c(0,1), xticks = c(6,11), yticks = c(6,11), xlab = "", ylab = "Proportion of mature tree stands", col.bg = NA, col.grid = NA,plot.axis=FALSE)
  axis(1, pos=0, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos = 0, at = seq(0,1,0.2), labels = seq(0,1,0.2))
  lines(newEq.all.For~gradient, col = 1, type = "l", lwd = 1.5)
  lines(newEq.veg.For~gradient, col = 1, lwd = 1.5, lty = 2)
  points(Hmax-1, newEq.all.For[which.max(newEq.all[,"H"])], pch=8)
  points(O2B-1, newEq.all.For[which(dominance=="B") [1]], pch=19)
  points(B2T-1, newEq.all.For[which(dominance=="T") [1]], pch=19)
  points(B2T.veg-1, newEq.veg.For[which(dominance.veg=="T") [1]], pch=1)
  mtext("(a)", 3, font = 2, cex=1.5)

  par(cex=1.4, mar=c(1,4,3,1))
  plot.template(xlim = c(0,50), ylim = c(0,0.01), xlab = "", ylab = "Asymptotic resilience", xticks = c(6,11), yticks = c(6,11), col.bg = NA, col.grid = NA, plot.axis=FALSE)
  axis(1, pos=0, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos = 0, at = seq(0,0.01,0.002), labels = seq(0,0.01,0.002))
  lines(-lambdamax.all~gradient, col = 1, type = "l", lwd = 1.5)
  lines(-lambdamax.veg~gradient, col = 1, lwd = 1.5, lty = 2)
  points(Hmax-1, -lambdamax.all[which.max(newEq.all[,"H"])], pch=8)
  points(O2B-1, -lambdamax.all[which(dominance=="B") [1]], pch=19)
  points(B2T-1, -lambdamax.all[which(dominance=="T") [1]], pch=19)
  points(B2T.veg-1, -lambdamax.veg[which(dominance.veg=="T") [1]], pch=1)
  mtext("(b)", 3, font = 2, cex=1.5)

  par(cex=1.4, mar=c(1,4,3,1))
  plot.template(xlim = c(0,50), ylim = c(-0.3,0), xlab = "", ylab = "Initial resilience", xticks = c(6,11), yticks = c(6,11), col.bg = NA, col.grid = NA, plot.axis=FALSE)
  axis(1, pos=-0.3, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos = 0, at = seq(-0.3,0,0.05), labels = seq(-0.3,0,0.05))
  lines(-reac.all~gradient, col = 1, type = "l", lwd = 1.5)
  lines(-reac.veg~gradient, col = 1, lwd = 1.5, lty = 2)
  points(Hmax-1, -reac.all[which.max(newEq.all[,"H"])], pch=8)
  points(O2B-1, -reac.all[which(dominance=="B") [1]], pch=19)
  points(B2T-1, -reac.all[which(dominance=="T") [1]], pch=19)
  points(B2T.veg-1, -reac.veg[which(dominance.veg=="T") [1]], pch=1)
  mtext("(c)", 3, font = 2, cex=1.5)

  par(cex=1.4, mar=c(1,4,1,1))
  plot(1:10,1:10, ylab="", xlab="", xaxt="n", yaxt="n", type="n", bty="n")
  mtext("Temperature", 1, line = -1.5, font = 2, col = 1, cex = 1.4)

  par(cex=1.4, mar=c(1,4,3,1))
  plot.template(xlim = c(0,50), ylim = c(0,0.15), xlab = "", ylab = "Exposure", xticks = c(6,11), yticks = c(6,11), col.bg = NA, col.grid = NA, plot.axis=FALSE)
  axis(1, pos=0, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos = 0, at = seq(0,0.15,0.05), labels = seq(0,0.15,0.05))
  lines(c(NA,deltaN.all)~gradient, col = 1, type = "l", lwd = 1.5)
  lines(c(NA, deltaN.veg)~gradient, col = 1, lwd = 1.5, lty = 2)
  points(Hmax-1, deltaN.all[which.max(newEq.all[,"H"])-1], pch=8)
  points(O2B-1, deltaN.all[which(dominance=="B") [1] -1], pch=19)
  points(B2T-1, deltaN.all[which(dominance=="T") [1] -1], pch=19)
  points(B2T.veg-1, deltaN.veg[which(dominance.veg=="T") [1] -1], pch=1)
  mtext("(d)", 3, font = 2, cex=1.5)

  par(cex=1.4, mar=c(1,4,3,1))
  plot.template(xlim = c(0,50), ylim = c(0,8000), xlab = "", ylab = "Sensitivity", xticks = c(6,11), yticks = c(6,11), col.bg = NA, col.grid = NA, plot.axis=FALSE)
  axis(1, pos=0, at = seq(0,50,le=6), labels =temp.range, cex.axis=.9)
  axis(2, pos = 0, at = seq(0,8000,2000), labels = seq(0,8000,2000))
  lines(c(NA,deltaT.all)~gradient, col = 1, type = "l", lwd = 1.5)
  lines(c(NA, deltaT.veg)~gradient, col = 1, lwd = 1.5, lty = 2)
  points(Hmax-1, deltaT.all[which.max(newEq.all[,"H"])-1], pch=8)
  points(O2B-1, deltaT.all[which(dominance=="B") [1] -1], pch=19)
  points(B2T-1, deltaT.all[which(dominance=="T") [1] -1], pch=19)
  points(B2T.veg-1, deltaT.veg[which(dominance.veg=="T") [1] -1], pch=1)
  mtext("(e)", 3, font = 2, cex=1.5)

  par(cex=1.4, mar=c(1,4,3,1))
  plot.template(xlim = c(0,50), ylim = c(0,35), xlab = "", ylab = "Cumulative state changes", xticks = c(6,11), yticks = c(6,11), col.bg = NA, col.grid = NA, plot.axis=FALSE)
  axis(1, pos=0, at = seq(0,50,le=6), labels = temp.range, cex.axis=.9)
  axis(2, pos = 0, at = seq(0,35,5), labels = seq(0,35,5))
  lines(c(NA,integraleV.all)~gradient, col = 1, type = "l", lwd = 1.5)
  lines(c(NA, integraleV.veg)~gradient, col = 1, lwd = 1.5, lty = 2)
  points(Hmax-1, integraleV.all[which.max(newEq.all[,"H"])-1], pch=8)
  points(O2B-1, integraleV.all[which(dominance=="B") [1] -1], pch=19)
  points(B2T-1, integraleV.all[which(dominance=="T") [1] -1], pch=19)
  points(B2T.veg-1, integraleV.veg[which(dominance.veg=="T") [1] -1], pch=1)
  mtext("(f)", 3, font = 2, cex=1.5)

  par(cex=1.4, mar=c(1,4,1,1))
  plot(1:10,1:10, ylab="", xlab="", xaxt="n", yaxt="n", type="n", bty="n")
  mtext("Temperature \n(final temperature)", 1, line = -1, font = 2, col = 1, cex = 1.4)
}
##---
plot_correlations <- function()
{
  layout(matrix(c(1,1,1,1,2,3,3,3,4,4,5,5,5,5,5), nrow = 3, byrow = TRUE), width = c(1,1,1,1,1), height=c(1,1, lcm(2)))
  par(cex.lab = 0.9, mar = c(4,4,4,2))

  # panel1
  xx= matrix(c(cRSn, cLSn, cDnSn, cDtSn), ncol = 2, byrow = T)
  par(lwd = 1)
  barplot(t(abs(xx)), beside = TRUE, ylim = c(0,1), ylab = "Absolute spearman correlation", main = "Correlation with cumulative state changes", space = c(0.1, 0.9))
  par(lwd = 2)
  barplot(t(abs(xx)), beside = TRUE, col = 1, density = c(0, 5, 5,5,0,0,0,0), border = NA, add = T, space = c(0.1, 0.9))
  axis(1, at = c(2, 5, 8, 11), labels = c("initial\n resilience", "asymptotic\n resilience", "exposure", "sensitivity"), tick =FALSE, cex.axis = 1.2)
  text(c(1.5,2.5), y=0.35, "n.s")
  box(which = "plot", col = "grey", lwd = 2)

  # panel2
  xx= matrix(c(cLR), ncol = 2, byrow = T)
  par(lwd = 1)
  barplot(t(abs(xx)), beside = TRUE, ylim = c(0,1), ylab = "Absolute spearman correlation", main = "Initial and\n asymptotic resilience", bty = "o", space = c(0.1, 0.9))
  par(lwd = 2)
  barplot(t(abs(xx)), beside = TRUE, col = 1, density = 5, border = NA, add = T, space = c(0.1, 0.9))
  box(which = "plot", col = "grey", lwd = 2)

  # panel3
  xx= matrix(c(cRDt, cLDt, cDnDt), ncol = 2, byrow = T)
  par(lwd = 1)
  barplot(t(abs(xx)), beside = TRUE, ylim = c(0,1), ylab = "Absolute spearman correlation", main = paste("Correlation with sensitivity (return time)"), space = c(0.1, 0.9))
  par(lwd = 2)
  barplot(t(abs(xx)), beside = TRUE, col = 1, density = c(0,0,5,5,0,0), border = NA, add = T, space = c(0.1, 0.9))
  axis(1, at = c(2, 5, 8), labels = c("initial\n resilience", "asymptotic\n resilience", "exposure"), tick =FALSE, cex.axis = 1.2)
  text(c(2.5), y=0.2, "n.s")
  box(which = "plot", col = "grey", lwd = 2)

  # panel4
  xx= matrix(c(cDnR, cDnL), ncol = 2, byrow = T)
  par(lwd = 1)
  barplot(t(abs(xx)), beside = TRUE, ylim = c(0,1), ylab = "Absolute spearman correlation", main = "Correlation with exposure", space = c(0.1, 0.9))
  par(lwd = 2)
  barplot(t(abs(xx)), beside = TRUE, col = 1, density = c(0,5, 5,5), border = NA, add = T, space = c(0.1, 0.9))
  axis(1, at = c(2, 5), labels = c("initial\n resilience", "asymptotic\n resilience"), tick =FALSE, cex.axis = 1.2)
  text(c(1.5,2.5), y=0.1, "n.s")
  box(which = "plot", col = "grey", lwd = 2)

  #panel 5
  par(mar = c(0.5, 0.1, 0.1, 0.1), bg="white")
  plot(1:10,1:10, type = "n", bty="n", xaxt = "n", yaxt="n")
  legend("top", legend = c("no herbivores", "with herbivores", "positive correlation", "negative correlation"), fill = c("grey20", "grey80", "white", 1), density =c(NA,NA,NA,10) , ncol=4, bty="n", pt.cex = 3, cex =1.3)
}
##---
plot_trajectories <- function()
{
  par(mfrow = c(1,3))
  par(cex=1, mar=c(4,4,3,1),lwd = 1)
  #range(simu2.veg[,"R"])
  miny = 0.093
  maxy = 0.096
  step = 0.001
  plot.template(xlim = c(0,500), ylim = c(miny,maxy), xlab = "Time", ylab = "Regeneration state proportion", xticks = c(6,11), yticks = c(6,11), col.bg = NA, col.grid = NA, plot.axis=FALSE)
  axis(1, pos = miny, at = seq(0,500,100), labels = seq(0,500,100))
  axis(2, pos = 0, at = seq(miny,maxy,step), labels = seq(miny,maxy,step))
  lines(simu2.veg[,"R"], col = 1)
  title("Without browser")
  abline(h=newEq.veg[35,"R"], lty = 2)
  abline(h=newEq.veg[36,"R"], lty = 2)

  #range(simu2[,"R"])
  miny = 0.203
  maxy = 0.208
  step = 0.002
  plot.template(xlim = c(0,500), ylim = c(miny,maxy), xlab = "Time", ylab = "Regeneration state proportion", xticks = c(6,11), yticks = c(6,11), col.bg = NA, col.grid = NA, plot.axis=FALSE)
  axis(1, pos = miny, at = seq(0,500,100), labels = seq(0,500,100))
  axis(2, pos = 0, at = seq(miny,maxy,step), labels = seq(miny,maxy,step))
  lines(simu2[,"R"], col = 1)
  title("With browser")
  abline(h=newEq.all[35,"R"], lty = 2)
  abline(h=newEq.all[36,"R"], lty = 2)

  par(mar=c(4,4,1,1))
  #range(simu2[,"H"])
  miny = 2565
  maxy = 2585
  step = 5
  plot.template(xlim = c(0,500), ylim = c(miny,maxy), xlab = "Time", ylab = "Herbivore biomass", xticks = c(6,11), yticks = c(6,11), col.bg = NA, col.grid = NA, plot.axis=FALSE)
  axis(1, pos = miny, at = seq(0,500,100), labels = seq(0,500,100))
  axis(2, pos = 0, at = seq(miny,maxy,step), labels = seq(miny,maxy,step))
  lines(simu2[,"H"], col = 1)
  abline(h=newEq.all[35,"H"], lty = 2)
  abline(h=newEq.all[36,"H"], lty = 2)
}
##--- pR
plot_correlations_pR <- function()
{
  layout(matrix(c(1,1,1,1,2,3,3,3,4,4,5,5,5,5,5), nrow = 3, byrow = TRUE), width = c(1,1,1,1,1), height=c(1,1, lcm(2)))
  par(cex = 1, mar = c(4,4,4,2))

  # panel1
  xx= matrix(c(cRSn, cLSn, cDnSn, cDtSn), ncol = 2, byrow = T)
  par(lwd = 1)
  barplot(t(abs(xx)), beside = TRUE, ylim = c(0,1), ylab = "Absolute spearman correlation", main = "Correlation with cumulative state changes", space = c(0.1, 0.9))
  par(lwd = 2)
  barplot(t(abs(xx)), beside = TRUE, col = 1, density = c(0, 5, 5,5,0,0,0,0), border = NA, add = T, space = c(0.1, 0.9))
  axis(1, at = c(2, 5, 8, 11), labels = c("initial\n resilience", "asymptotic\n resilience", "exposure", "sensitivity"), tick =FALSE, cex.axis = 1.1)
  text(c(1.5), y=0.2, "n.s")
  box(which = "plot", col = "grey", lwd = 2)

  # panel2
  xx= matrix(c(cLR), ncol = 2, byrow = T)
  par(lwd = 1)
  barplot(t(abs(xx)), beside = TRUE, ylim = c(0,1), ylab = "Absolute spearman correlation", main = "Initial and\n asymptotic resilience", bty = "o", space = c(0.1, 0.9))
  par(lwd = 2)
  barplot(t(abs(xx)), beside = TRUE, col = 1, density = c(5,0), border = NA, add = T, space = c(0.1, 0.9))
  box(which = "plot", col = "grey", lwd = 2)

  # panel3
  xx= matrix(c(cRDt, cLDt, cDnDt), ncol = 2, byrow = T)
  par(lwd = 1)
  barplot(t(abs(xx)), beside = TRUE, ylim = c(0,1), ylab = "Absolute spearman correlation", main = paste("Correlation with sensitivity (return time)"), space = c(0.1, 0.9))
  par(lwd = 2)
  barplot(t(abs(xx)), beside = TRUE, col = 1, density = c(0,5,5,5,0,0), border = NA, add = T, space = c(0.1, 0.9))
  axis(1, at = c(2, 5, 8), labels = c("initial\n resilience", "asymptotic\n resilience", "exposure"), tick =FALSE, cex.axis = 1.1)
  text(c(8.5), y=0.3, "n.s")
  box(which = "plot", col = "grey", lwd = 2)

  # panel4
  xx= matrix(c(cDnR, cDnL), ncol = 2, byrow = T)
  par(lwd = 1)
  barplot(t(abs(xx)), beside = TRUE, ylim = c(0,1), ylab = "Absolute spearman correlation", main = "Correlation with exposure", space = c(0.1, 0.9))
  par(lwd = 2)
  barplot(t(abs(xx)), beside = TRUE, col = 1, density = c(0,5, 5,5), border = NA, add = T, space = c(0.1, 0.9))
  axis(1, at = c(2, 5), labels = c("initial\n resilience", "asymptotic\n resilience"), tick =FALSE, cex.axis = 1.1)
  text(c(1.5,2.5,5.5), y=0.2, "n.s")
  box(which = "plot", col = "grey", lwd = 2)

  #panel 5
  par(mar = c(0.5, 0.1, 0.1, 0.1), bg="white")
  plot(1:10,1:10, type = "n", bty="n", xaxt = "n", yaxt="n")
  legend("top", legend = c("no herbivores", "with herbivores", "positive correlation", "negative correlation"), fill = c("grey20", "grey80", "white", 1), density =c(NA,NA,NA,10) , ncol=4, bty="n", pt.cex = 3, cex =1.3)

}
