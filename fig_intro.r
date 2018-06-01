
xx1 = seq(1,12,le=100)
xx2 = seq(2,10,le=100)

yy1 = dchisq(xx1, df=5)
yy2 = dchisq(xx2, df=5)

plot(xx1,yy1, type = "l", xlim = c(1,10))
lines(xx2,yy2)


pdf("../graphs/fig_intro.pdf", width=12, height=5)

layout(matrix(c(1,2,1,3), byrow = T, nrow= 2), width = c(3,2))

xlim = range(xx2)
ylim = range(yy2)

par(mar = c(4,5,2,2), cex=1.3)

plot(xx2,yy2, type = "l", xlim = c(2,12), bty = "n", xaxt="n", yaxt="n", ylab = "" ,xlab = "", ylim = c(0,1.2*max(yy2)))
polygon(c(xx2, rev(xx2)), c(rep(min(yy2),length(xx2)), rev(yy2)), density = NA, col = "grey")
lines(xx2,yy2)

arrows(xlim[1], ylim[1], xlim[1], 1.2*ylim[2], le = .1)
arrows(xlim[1], ylim[1],1.1* xlim[2], ylim[1], le = .1)

axis(2, at = median(yy2), labels= "State",line=-1, lwd=0, lwd.ticks=0, las=3, cex.axis=1)
axis(1, at = median(xx2), labels= "Time",line=-2.5, lwd=0, lwd.ticks=0, cex.axis = 1)
#delta N
axis(2, at = c(min(yy2), yy2[1]), labels= c(expression(N[1]),expression(N[0])),line=-1, lwd=0, lwd.ticks=1, las=2, cex.axis=1)

arrows(1, ylim[1], 1, yy2[1], le = .1, code=3, xpd=NA)
mtext(text = expression(Delta[N]), side = 2, line = 1.9, at = 1.3*yy2[1]/2, las = 2, cex=1.7)

#delta T
axis(1, at = range(xx2), labels= c(expression(t[0]),expression(t[1])),line=-2.5, lwd=0, lwd.ticks=1, cex.axis = 1)
arrows(xlim[1], -0.01, xlim[2], -0.01, le = .1, code=3, xpd=NA)
mtext(text = expression(Delta[t]), side = 1, line = 0.5, at = 6, cex=1.7)


#REACTIVITY
text(x = quantile(xx2, 0.05), y = max(yy2)*1.1, labels = expression(-R[0]), cex=1.1)
arrows(xlim[1], yy2[1] ,quantile(xx2, 0.15), yy2[1]*1.25, le = .1, lwd=1.5)

#lambda
text(x = 0.9*xx2[which.min(yy2)], y = 2*min(yy2), labels = expression(R[infinity]), cex = 1.1)
arrows(quantile(xx2, 0.75), quantile(yy2, 0.2) ,xlim[2], min(yy2), le = .1, lwd=1.5)

#integrale
text(x = 1.2*xx2[which.max(yy2)], y = max(yy2)/2, labels = expression(integral(N(t) * dt)), bg = "white", cex = 1.1)

#dev.off()

###

#par(mfrow = c(2,1))
par(mar = c(1,1,2,2), cex=1.2)

xlim = range(xx2)

yy3 = exp(1-xx2)
ylim = range(yy3)

plot(xx2,yy3, type = "l", xlim = c(2,12), bty = "n", xaxt="n", yaxt="n", ylab = "" ,xlab = "", ylim = c(0,1.2*max(yy3)))
polygon(c(xx2, rev(xx2)), c(rep(min(yy3),length(xx2)), rev(yy3)), density = NA, col = "grey")
lines(xx2,yy3)
arrows(xlim[1], ylim[1], xlim[1], 1.2*ylim[2], le = .1)
arrows(xlim[1], ylim[1],1.1* xlim[2], ylim[1], le = .1)
title("case A", font = 2)
axis(2, at = mean(ylim), labels= "State",line=-1, lwd=0, lwd.ticks=0, las=3, cex.axis=1)
axis(1, at = median(xx2), labels= "Time",line=-1, lwd=0, lwd.ticks=0, cex.axis = 1)

ylim = range(yy1)
xlim = range(xx1)
plot(xx1,yy1, type = "l", xlim = xlim, bty = "n", xaxt="n", yaxt="n", ylab = "" ,xlab = "", ylim = c(0,1.2*max(yy1)))
polygon(c(xx1, rev(xx1)), c(rep(min(yy1),length(xx1)), rev(yy1)), density = NA, col = "grey")
lines(xx1,yy1)
arrows(xlim[1], ylim[1], xlim[1], 1.2*ylim[2], le = .1)
arrows(xlim[1], ylim[1],1.1* xlim[2], ylim[1], le = .1)
title("case B", font = 2)
axis(2, at = mean(ylim), labels= "State",line=-1, lwd=0, lwd.ticks=0, las=3, cex.axis=1)
axis(1, at = median(xx1), labels= "Time",line=-1, lwd=0, lwd.ticks=0, cex.axis = 1)

dev.off()
