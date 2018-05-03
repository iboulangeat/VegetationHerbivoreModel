#----------------------------------------------------------------------------
# SOLVE EQ FUNCTIONS
#----------------------------------------------------------------------------

solveEq <- function(func , init , parms , maxsteps, veget=F )
{
  
  nochange = 0
  
  trace.mat = matrix(NA, ncol  = length(init), nrow = maxsteps+1 )
  trace.mat[1,] = init
  state = init
  for (i in 1:maxsteps)
  {
    di = func(1, state, parms)
    state = state + di[[1]]
    if(state[4] <0) state[4]=0
    
    trace.mat[i+1,] = state
    if(!veget) {
      if(sum(abs(trace.mat[i,]-trace.mat[i-1,]))<1e-10) nochange =nochange+1
    }else{
      if(sum(abs(trace.mat[i,-4]-trace.mat[i-1,-4]))<1e-10) nochange =nochange+1
    }
    
    if(nochange>=10) break;
    
  }
  trace.mat = trace.mat[1:i,]
  
  return(list(eq = state, mat = trace.mat))
}

####------------------------------------------

solveEqH <- function(func , init , parms , maxsteps )
{
  
  nochange = 0
  
  trace.mat = matrix(NA, ncol  = length(init), nrow = maxsteps+1 )
  trace.mat[1,] = init
  state = c(T=as.numeric(init))
  for (i in 1:maxsteps)
  {
    di = func(1, state, parms)
    state = state + di[[1]]
    if(state <0) state=0
    trace.mat[i+1,] = state
    if(sum(abs(trace.mat[i,]-trace.mat[i-1,]))<1e-10) nochange =nochange+1
    
    if(nochange>=10) break;
  }
  trace.mat = trace.mat[1:i,]
  
  return(list(eq = state, mat = trace.mat))
}

### ---- 
#  CASE WHERE H=0 (only vegetation dynamics)
### ---


equi_veget <- function(alpha, c, k, deltaB, deltaT)
{
  Tstar = alpha*deltaB*k*(c*deltaB*k - c*deltaT*k + c*deltaT - deltaB*deltaT)/(c*(deltaB*k - deltaT*k + deltaT)*(alpha*deltaB*k - alpha*deltaT*k + alpha*deltaT + deltaB*deltaT))
  Bstar = -alpha*deltaT*(k - 1)*(c*deltaB*k - c*deltaT*k + c*deltaT - deltaB*deltaT)/(c*(deltaB*k - deltaT*k + deltaT)*(alpha*deltaB*k - alpha*deltaT*k + alpha*deltaT + deltaB*deltaT))
  Sstar = deltaB*deltaT*(c*deltaB*k - c*deltaT*k + c*deltaT - deltaB*deltaT)/(c*(deltaB*k - deltaT*k + deltaT)*(alpha*deltaB*k - alpha*deltaT*k + alpha*deltaT + deltaB*deltaT))
  eq = c(T = Tstar, S=Sstar, B=Bstar)
  return(eq)
}


#----------------------------------------------------------------------------
# RUN EQ, VEGET ONLY NO FEEDBACKS
#----------------------------------------------------------------------------

eqveg.fct <- function(g, init.params, par.name, par.clim)
{
  fparams = init.params
  fparams[par.name] = unlist(par.clim[g,])
  
  equi.veg = equi_veget(alpha=fparams["a0"], c = fparams["c0"], k = fparams["k"], deltaT = fparams["dT"], deltaB = fparams["dB"])
  names(equi.veg) = c("T", "S", "B")
  if(min(equi.veg >0))
  {
    require(rootSolve)
    out = runsteady(func = modelH, y =c(H=5*375), parms = c(fparams,equi.veg ), maxsteps = 10000)
    eq = out$y
  } else eq = 0
  
  newEq = c(equi.veg,eq)
  
  # equilibrium measures
  require(rootSolve)
  jacob = jacobian.full(y = c(equi.veg,H=0), func= model, parms=fparams)
  
  lambdamax = max(Re(eigen(jacob[-4,-4])$values))
  
  #reactivite -> potentiel à amplifier les fluctuations env.
  M = (jacob[-4,-4] + t(jacob[-4,-4]))/2
  reac = max(eigen(M)$values) 
  
  return(list(newEq, lambdamax, reac))
}

eqveg.fct.Vect = Vectorize(eqveg.fct, "g")

#----------------------------------------------------------------------------
# RUN EQ, VEGET AND HERBIVORES
#----------------------------------------------------------------------------

eqall.fct <- function(g, init.params, par.name, par.clim, model)
{
  fparams = init.params
  fparams[par.name] = par.clim[g,]
  
  T0 = c(T=0.3, S=0.3, B = 0.3, H=5*375)
  require(rootSolve)
  out = solveEq(func = model, init =T0, parms = fparams, maxsteps = 10000)
  out = solveEq(func = model, init =out$eq, parms = fparams, maxsteps = 1000)
  eq = apply(out$mat, 2, mean)
  osc = apply(out$mat, 2, sd)
  names(eq) = c("T", "S", "B", "H")
  
  # equilibrium measures
  require(rootSolve)
  jacob = jacobian.full(y = eq, func= model, parms=fparams)
  lambdamax = max(Re(eigen(jacob[-4,-4])$values))
  
  #reactivity
  M = (jacob[-4,-4] + t(jacob[-4,-4]))/2
  reac = max(eigen(M)$values) 
  
  return(list(eq, osc, lambdamax, reac))
}

eqall.fct.Vect = Vectorize(eqall.fct, "g")

#----------------------------------------------
# integrale et deltaT -- calculs
#----------------------------------------------
simulation <- function(newEq, grad.start, grad.pars, simu.time, woH = FALSE, init.params, model)
{
  simu.init = as.numeric(newEq[grad.start,1:4])
  names(simu.init) = colnames(newEq)[1:4]
  simu.params =init.params
  simu.params[par.name] = as.numeric(par.clim[grad.pars,])
  
  if(woH) simu.init[4] = 0
  
  simu = solveEq(model, simu.init, simu.params, maxsteps = simu.time, veget=TRUE)
  mat = data.frame(simu$mat)
  colnames(mat) = c("T", "S", "B", "H")
  mat$G = 1-mat$T - mat$S - mat$B
  return(mat)
}


calcTime <- function(g, newEq, simu.time, woH, init.params, model)
{
  simu = simulation(newEq, grad.start=g, grad.pars=g+1, simu.time, woH, init.params =init.params, model = model)
  simu = simu[-c((nrow(simu)-10):nrow(simu)),]
  deltaT = nrow(simu)
  # N.B. because the simu breaks when there is 10 times steps without total change in state >10e-10
  
  integrale.V = sum(abs(sweep(simu[,-4], 2, unlist(newEq[g+1,-4]) , "-")))
  return(list(deltaT, integrale.V))
  
}# END CALCTIME

calcTime.Vect = Vectorize(calcTime, "g")

#----------------------------------------------
# APPLY ANALYSIS FUNCTION FOR SA
#----------------------------------------------

applyAnalysis <- function(init.params, grad.div, gradient,par.clim, par.name, model)
{
  
  #------------------------------------------
  # run equilibrium and calc jacobian measures
  #------------------------------------------
  
  #vegetation alone
  #------------------
  
  calcEq.veg = data.frame(t(eqveg.fct.Vect(1:grad.div, init.params, par.name, par.clim)))
  newEq.veg = data.frame(matrix(unlist(calcEq.veg[[1]]), ncol = 4, byrow = T))
  colnames(newEq.veg) = c("T", "S", "B", "H")
  newEq.veg[,"G"] = 1 - newEq.veg[,"S"] - newEq.veg[,"B"] - newEq.veg[,"T"]
  lambdamax.veg = unlist(calcEq.veg[[2]])
  reac.veg = unlist(calcEq.veg[[3]])
  
  #vegetation and herbivores
  #---------------------------
  
  calcEq.all = data.frame(t(eqall.fct.Vect(1:grad.div, init.params, par.name, par.clim, model)))
  newEq.all = data.frame(matrix(unlist(calcEq.all[[1]]), ncol = 4, byrow = T))
  colnames(newEq.all) = c("T", "S", "B", "H")
  newEq.all[,"G"] = 1 - newEq.all[,"S"] - newEq.all[,"B"] - newEq.all[,"T"]
  oscillations = data.frame(matrix(unlist(calcEq.all[[2]]), nrow = grad.div, byrow = TRUE))
  colnames(oscillations) = c("T", "S", "B", "H")
  lambdamax.all = unlist(calcEq.all[[3]])
  reac.all = unlist(calcEq.all[[4]])
  
  # dominance
  #--------------------------------
  open = apply(newEq.all[,c("G","S")],1,sum)
  dominance = c("T","B", "O")[apply(cbind(newEq.all[,c("T","B")],open), 1, which.max)]
  dominance.veg = c("T","B", "O")[apply(cbind(newEq.veg[,c("T","B")],open), 1, which.max)]
  O2B = which(dominance=="B") [1]
  B2T = which(dominance=="T") [1]
  O2B.veg = which(dominance.veg=="B") [1]
  B2T.veg = which(dominance.veg=="T") [1]
  
  #------------------------------------------
  
  #--------------------------------
  #  deltaN
  #--------------------------------
  
  deltaN.all = deltaN.veg =  rep(NA, 50)
  require(cluster)
  dd2<-daisy(newEq.all[,-4], metric="euclidean")
  dd.veg2 <- daisy(newEq.veg[,-4], metric="euclidean")
  
  for (g in 1:50)
  {
    deltaN.all[g] = as.matrix(dd2)[g,g+1]
    deltaN.veg[g] = as.matrix(dd.veg2)[g,g+1]
  }
  
  #----------------------------------------------
  # integrale and deltaT 
  #----------------------------------------------
  
  simu.time = 10000 # max time of simu before reach eq
  
  calcTime.all = calcTime.Vect(1:50, newEq.all, simu.time=10000, woH=FALSE, init.params=init.params, model)
  deltaT.all = unlist(calcTime.all[1,])
  integraleV.all = unlist(calcTime.all[2,]) 
  
  calcTime.veg = calcTime.Vect(1:50, newEq.veg, simu.time=10000, woH=TRUE, init.params=init.params, model)
  deltaT.veg = unlist(calcTime.veg[1,])
  integrale.veg = unlist(calcTime.veg[2,])
  
  stats = list(newEq.veg=newEq.veg, newEq.all=newEq.all, osc=oscillations, 
               O2B=O2B, B2T=B2T, lambdamax.all=lambdamax.all , reac.all=reac.all , 
               deltaN.all=deltaN.all[,2] , deltaT.all=deltaT.all, 
               integraleV.all =integraleV.all )
  
  return(stats)
  
} # END FCT


