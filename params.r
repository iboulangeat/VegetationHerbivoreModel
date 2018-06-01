#------------------------------------------
### parms
#------------------------------------------

init.params =  c(
aT = 0.05,
aB = 0.05,
dT = 0.004,
dB = 0.01,
c0 = 0.7,
k = 0.4,
es = 0.35,
ew = 0.25,
m = .30,
taus = 124* 10/358, # growing season =124 days
tauw = (365-124)*10/358,
gamma = 1,
nus = 100*124/365, # nu =100, growing season =124 days
nuw = 100*(1-124/365),
uR = 600*100 ,
uT =  200*100,
uB =  200*100,
f = .85,
fB = 0, # proportion B browsed in summer, not used here
p = 0.8, # preference T /B
pR = 0.5, # preference R / B or T
r = 10,
hT = 0.05 ,
hB = 0.05,
hR = 0.05,
uG = 500*100,
pG = 1
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
