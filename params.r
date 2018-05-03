#------------------------------------------
### parms
#------------------------------------------

init.params =  c(
a0 = 0.05, 
dT = 0.004, 
dB = 0.01,
c0 = 0.7,
k = 0.4,
es = 0.35, 
ew = 0.25, 
mp = .30,
taus = 10/358, 
tauw = 10/358, 
gmax = 1, 
gseason = 124, 
uS = 600*100 , 
uT =  200*100, 
l = .5, 
lB = .2,
uB =  200*100, 
nu = 100, 
ptresh = 0.05 ,
omega = 0.8, 
omega2 = 0.5, 
r = 10, 
uG = 800*100,
pGraz = 1
)


#------------------------------------------
# parameter variation -- climate change
#------------------------------------------
grad.div = 51
gradient = seq(0, grad.div-1, length.out = grad.div)

xx = seq(-2,8,len = grad.div)
fit= list()
fit[["T"]] = dnorm(xx, 5.2,2.1)
fit[["B"]] = dnorm(xx, 3,2.2)
fit[["S"]] = apply(cbind(fit[["B"]],fit[["T"]]), 1, function(x){(1*x[2]+1*x[1])})
 
# ---- tree competitivity climate change
temp.range = seq(-2,8,2)
c0.var = .8*fit[["S"]]/max(fit[["S"]])
k.var = fit[["T"]] / (fit[["B"]]+fit[["T"]])

c0max = which.max(c0.var)
par.name = c("c0", "k")
par.clim = data.frame(c0 = c0.var, k = k.var)
