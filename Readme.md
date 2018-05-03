---
title: "Main analysis"
author: "Isabelle Boulangeat"
date: "02/05/2018"
output: 
  html_document:
      keep_md: yes
      theme: cosmo
      highlight: tango
      number_sections: true
      toc: true
---

<!-- library(rmarkdown) -->
<!-- library(knitr) -->
<!-- knit("MainAnalysis.Rmd", "Readme.md") -->

# Set up

## Load parameters and model

Also load plot functions

```r
source("params.r")
source("model_fct.r")
source("plot_fct.r")
source("analysis_fct.r")
```

## Libraries


```r
library(rootSolve)
library(cluster)
```

## Test model


```r
T0 = c(T=0.3, S=0.3, B = 0.3, H=5*375)
out = solveEq(func = model, init =T0, parms = init.params, maxsteps = 10000)
(eq = out$eq)
```

```
##            T            S            B            H 
##    0.3250468    0.2510789    0.3977062 2617.3140121
```
#  Equilibrium along environmental gradient

## Vegetation alone


```r
calcEq.veg = data.frame(t(eqveg.fct.Vect(1:grad.div, init.params, par.name, par.clim)))
newEq.veg = data.frame(matrix(unlist(calcEq.veg[[1]]), ncol = 4, byrow = T))
colnames(newEq.veg) = c("T", "S", "B", "H")
newEq.veg[,"G"] = 1 - newEq.veg[,"S"] - newEq.veg[,"B"] - newEq.veg[,"T"]
lambdamax.veg = unlist(calcEq.veg[[2]])
reac.veg = unlist(calcEq.veg[[3]])
newEq.veg.h0 = newEq.veg
newEq.veg.h0[,"H"] = 0
newEq.veg.For = newEq.veg[,"T"] + newEq.veg[,"B"]
```

## Vegetation and herbivores


```r
calcEq.all = data.frame(t(eqall.fct.Vect(1:grad.div, init.params, par.name, par.clim, model)))
newEq.all = data.frame(matrix(unlist(calcEq.all[[1]]), ncol = 4, byrow = T))
colnames(newEq.all) = c("T", "S", "B", "H")
newEq.all[,"G"] = 1 - newEq.all[,"S"] - newEq.all[,"B"] - newEq.all[,"T"]
oscillations = data.frame(matrix(unlist(calcEq.all[[2]]), nrow = grad.div, byrow = TRUE))
colnames(oscillations) = c("T", "S", "B", "H")
lambdamax.all = unlist(calcEq.all[[3]])
reac.all = unlist(calcEq.all[[4]])
newEq.all.For = newEq.all[,"T"] + newEq.all[,"B"]
```

##  dominance zones


```r
Smax = which.max(newEq.all[,"S"])
Bmax = which.max(newEq.all[,"B"])
Hmax = which.max(newEq.all[,"H"])
open = newEq.all[,"G"]
dominance = c("T","B", "O")[apply(cbind(newEq.all[,c("T","B")],open), 1, which.max)]
open.veg = newEq.veg[,"G"]
dominance.veg = c("T","B", "O")[apply(cbind(newEq.veg[,c("T","B")],open.veg), 1, which.max)]
O2B = which(dominance=="B") [1]
B2T = which(dominance=="T") [1]
O2B.veg = which(dominance.veg=="B") [1]
B2T.veg = which(dominance.veg=="T") [1]
```

## Equilibrium along the environmental gradient: figure



![plot of chunk equilibrium](figure/equilibrium-1.png)

# Transition metrics

##  deltaN


```r
deltaN.all = deltaN.veg =  rep(NA, 50)
require(cluster)
dd2<-daisy(newEq.all[,-4], metric="euclidean")
dd.veg2 <- daisy(newEq.veg[,-4], metric="euclidean")
for (g in 1:50)
{
deltaN.all[g] = as.matrix(dd2)[g,g+1]
deltaN.veg[g] = as.matrix(dd.veg2)[g,g+1]
}
```


## Return to equilibrium after climate change


```r
simu.time = 10000 # max time of simu before reach eq

calcTime.all = calcTime.Vect(1:50, newEq.all, simu.time=10000, woH=FALSE,init.params=init.params, model)
deltaT.all = unlist(calcTime.all[1,])
integraleV.all = unlist(calcTime.all[2,]) 

calcTime.veg = calcTime.Vect(1:50, newEq.veg, simu.time=10000, woH=TRUE, init.params=init.params, model)
deltaT.veg = unlist(calcTime.veg[1,])
integraleV.veg = unlist(calcTime.veg[2,]) 
```

## All metrics of transient dynamics along the gradient: figure

![plot of chunk transient](figure/transient-1.png)

# Correlations between metrics

## Correlations and tests


```r
(cLDt = c(cor(-lambdamax.veg[-1], deltaT.veg, method="spearman"), cor(-lambdamax.all[-1], deltaT.all, method="spearman")) )
```

```
## [1] -1.0000000 -0.9789676
```

```r
cor.test(-lambdamax.veg[-1], deltaT.veg, method="spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  -lambdamax.veg[-1] and deltaT.veg
## S = 41650, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
## rho 
##  -1
```

```r
cor.test(-lambdamax.all[-1], deltaT.all,  method="spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  -lambdamax.all[-1] and deltaT.all
## S = 41212, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## -0.9789676
```

```r
(cLR = c(cor(-lambdamax.veg[-1], -reac.veg[-1], method="spearman"), cor(-lambdamax.all[-1], -reac.all[-1], method="spearman")) )
```

```
## [1] -0.5303721 -0.2250660
```

```r
cor.test(-lambdamax.veg[-1], -reac.veg[-1], method="spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  -lambdamax.veg[-1] and -reac.veg[-1]
## S = 31870, p-value = 9.527e-05
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## -0.5303721
```

```r
cor.test(-lambdamax.all[-1], -reac.all[-1], method="spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  -lambdamax.all[-1] and -reac.all[-1]
## S = 25512, p-value = 0.116
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## -0.225066
```

```r
(cDnSn = c(cor(deltaN.veg, integraleV.veg, method="spearman"), cor(deltaN.all, integraleV.all, method="spearman")) )
```

```
## [1] 0.8097479 0.7170708
```

```r
cor.test(deltaN.veg, integraleV.veg, method="spearman") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  deltaN.veg and integraleV.veg
## S = 3962, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.8097479
```

```r
cor.test(deltaN.all, integraleV.all, method="spearman") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  deltaN.all and integraleV.all
## S = 5892, p-value = 2.538e-08
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.7170708
```

```r
(cDtSn = c(cor(deltaT.veg, integraleV.veg, method="spearman"), cor(deltaT.all, integraleV.all, method="spearman")) )
```

```
## [1] 0.7510684 0.6601200
```

```r
cor.test(deltaT.veg, integraleV.veg, method="spearman") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  deltaT.veg and integraleV.veg
## S = 5184, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.7510684
```

```r
cor.test(deltaT.all, integraleV.all, method="spearman") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  deltaT.all and integraleV.all
## S = 7078, p-value = 4.428e-07
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##     rho 
## 0.66012
```

```r
(cRDt = c(cor(-reac.veg[-1], deltaT.veg, method="spearman"), cor(-reac.all[-1], deltaT.all, method="spearman")) )
```

```
## [1] 0.53037215 0.09618247
```

```r
cor.test(-reac.veg[-1], deltaT.veg, method="spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  -reac.veg[-1] and deltaT.veg
## S = 9780, p-value = 9.527e-05
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.5303721
```

```r
cor.test(-reac.all[-1], deltaT.all, method="spearman") # NS
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  -reac.all[-1] and deltaT.all
## S = 18822, p-value = 0.5053
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## 0.09618247
```

```r
(cLSn = c(cor(-lambdamax.veg[-1], integraleV.veg, method="spearman"), cor(-lambdamax.all[-1], integraleV.all, method="spearman")) )
```

```
## [1] -0.7510684 -0.5855942
```

```r
cor.test(-lambdamax.veg[-1], integraleV.veg, method="spearman") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  -lambdamax.veg[-1] and integraleV.veg
## S = 36466, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## -0.7510684
```

```r
cor.test(-lambdamax.all[-1], integraleV.all, method="spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  -lambdamax.all[-1] and integraleV.all
## S = 33020, p-value = 1.185e-05
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## -0.5855942
```

```r
(cRSn = c(cor(-reac.veg[-1], integraleV.veg, method="spearman"), cor(-reac.all[-1], integraleV.all, method="spearman")) )
```

```
## [1]  0.1469868 -0.2303481
```

```r
cor.test(-reac.veg[-1], integraleV.veg, method="spearman") # NS
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  -reac.veg[-1] and integraleV.veg
## S = 17764, p-value = 0.3074
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.1469868
```

```r
cor.test(-reac.all[-1], integraleV.all, method="spearman") # NS
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  -reac.all[-1] and integraleV.all
## S = 25622, p-value = 0.1075
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## -0.2303481
```

```r
(cDnDt = c(cor(deltaN.veg, deltaT.veg, method="spearman"), cor(deltaN.all, deltaT.all, method="spearman")) )
```

```
## [1] 0.7050660 0.3544298
```

```r
cor.test(deltaN.veg, deltaT.veg, method="spearman") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  deltaN.veg and deltaT.veg
## S = 6142, p-value = 5.087e-08
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##      rho 
## 0.705066
```

```r
cor.test(deltaN.all, deltaT.all, method="spearman") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  deltaN.all and deltaT.all
## S = 13444, p-value = 0.01195
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.3544298
```

```r
(cDnL = c(cor(deltaN.veg, -lambdamax.veg[-1], method="spearman"), cor(deltaN.all, -lambdamax.all[-1], method="spearman")) )
```

```
## [1] -0.705066 -0.332533
```

```r
cor.test(deltaN.veg, -lambdamax.veg[-1], method="spearman") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  deltaN.veg and -lambdamax.veg[-1]
## S = 35508, p-value = 5.087e-08
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## -0.705066
```

```r
cor.test(deltaN.all, -lambdamax.all[-1], method="spearman") 
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  deltaN.all and -lambdamax.all[-1]
## S = 27750, p-value = 0.01872
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## -0.332533
```

```r
(cDnR = c(cor(deltaN.veg, -reac.veg[-1], method="spearman"), cor(deltaN.all, -reac.all[-1], method="spearman")) )
```

```
## [1]  0.01954382 -0.01829532
```

```r
cor.test(deltaN.veg, -reac.veg[-1], method="spearman") # NS
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  deltaN.veg and -reac.veg[-1]
## S = 20418, p-value = 0.8927
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## 0.01954382
```

```r
cor.test(deltaN.all, -reac.all[-1], method="spearman") # NS
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  deltaN.all and -reac.all[-1]
## S = 21206, p-value = 0.8995
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##         rho 
## -0.01829532
```

## Correlations : figure

![plot of chunk correlations](figure/correlations-1.png)


# Simulations of trajectories

## Simulation runs

```r
simu1 = simulation(newEq.all, grad.start=15, grad.pars=15+1, 10000, woH=FALSE, init.params =init.params, model)
simu1.veg = simulation(newEq.veg, grad.start=15, grad.pars=15+1, 10000, woH=TRUE, init.params =init.params, model)
simu2 = simulation(newEq.all, grad.start=35, grad.pars=35+1, 10000, woH=FALSE, init.params =init.params, model)
simu2.veg = simulation(newEq.veg, grad.start=35, grad.pars=35+1, 10000, woH=TRUE, init.params =init.params, model)
```

## Trajectories' examples figure

![plot of chunk trajectories](figure/trajectories-1.png)

# Variation for a mix-feeder animal

## Test of the model


```r
T0 = c(T=0.3, S=0.3, B = 0.3, H=5*375)
require(rootSolve)
out = solveEq(func = modelG, init =T0, parms = init.params, maxsteps = 10000)
(eq = out$eq)
```

```
##            T            S            B            H 
##    0.3084180    0.2637219    0.4110798 2990.6327832
```

## Equilibrium along the gradient


```r
calcEq.all = data.frame(t(eqall.fct.Vect(1:grad.div, init.params, par.name, par.clim, modelG)))
newEq.all = data.frame(matrix(unlist(calcEq.all[[1]]), ncol = 4, byrow = T))
colnames(newEq.all) = c("T", "S", "B", "H")
newEq.all[,"G"] = 1 - newEq.all[,"S"] - newEq.all[,"B"] - newEq.all[,"T"]
oscillations = data.frame(matrix(unlist(calcEq.all[[2]]), nrow = grad.div, byrow = TRUE))
colnames(oscillations) = c("T", "S", "B", "H")
lambdamax.all = unlist(calcEq.all[[3]])
reac.all = unlist(calcEq.all[[4]])

Smax = which.max(newEq.all[,"S"])
Bmax = which.max(newEq.all[,"B"])
Hmax = which.max(newEq.all[,"H"])
open = newEq.all[,"G"]
dominance = c("T","B", "O")[apply(cbind(newEq.all[,c("T","B")],open), 1, which.max)]
open.veg = newEq.veg[,"G"]
dominance.veg = c("T","B", "O")[apply(cbind(newEq.veg[,c("T","B")],open.veg), 1, which.max)]
O2B = which(dominance=="B") [1]
B2T = which(dominance=="T") [1]
O2B.veg = which(dominance.veg=="B") [1]
B2T.veg = which(dominance.veg=="T") [1]
```



![plot of chunk equiMixFeed](figure/equiMixFeed-1.png)

## Transient metrics



```r
newEq.all.For = newEq.all[,"T"] + newEq.all[,"B"]

deltaN.all =rep(NA, 50)
require(cluster)
dd2<-daisy(newEq.all[,-4], metric="euclidean")
for (g in 1:50)
{
  deltaN.all[g] = as.matrix(dd2)[g,g+1]
}
simu.time = 10000 # max time of simu before reach eq
calcTime.all = calcTime.Vect(1:50, newEq.all, simu.time=10000, woH=FALSE,init.params=init.params, model=modelG)
deltaT.all = unlist(calcTime.all[1,])
integraleV.all = unlist(calcTime.all[2,]) 
```

![plot of chunk transientMisFeed](figure/transientMisFeed-1.png)


