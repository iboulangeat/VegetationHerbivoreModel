
#______________________________________________________________________
# model
#------------------------------------------
model <- function(ti, states, parms)
{
  with(as.list(c(states, parms)), {

    V=1-T-R-B

    # preferences seedlings temperate/boreal
    epsilon = 0.0001
    omega = (T+epsilon) / (B+T+2*epsilon)
    q = p*omega / (p*omega + (1-omega)*(1-p))

    # preferences seedlings/mature
    phis_0 = (uR*R)*f / ((uR*R+uT*T)*f +uB*B*fB)
    phiw_0 = (uR*R*(1-f)) / ((uR*R + uT*T)*(1-f) + uB*B*(1-fB))
    phis = pR*phis_0 / (pR*phis_0+(1-pR)*(1-phis_0))
    phiw = pR*phiw_0 / (pR*phiw_0+(1-pR)*(1-phiw_0))

    if(H>1e-5) # population cannot recover if no reproduction is possible
    {
      Fs = (uR*R + uT*T)*f + uB*B*fB
      Fw = (uR*R + uT*T)*(1-f) + uB*B*(1-fB)

      # intakes rates
      irate_s = (taus*Fs/H) / (nus + Fs/H)
      if(irate_s>(Fs/H)) irate_s = Fs/H # maximum consumption

      irate_w = (tauw*Fw/H) / (nuw + Fw/H)
      if(irate_w>(Fw/H)) irate_w = Fw/H

      # summer gain, reproduction limitation
      Gs = H* min(gamma, es* irate_s )

      # winter gain
      Gw = H* ew* irate_w

      # impacts
      if(R>0)
      {
        # R intake
        Us = phis * Gs/es
        Uw = phiw * Gw/ew

        PR = (Us+Uw)/(uR*R)
        if(k==0) PT=0 else PT = PR*q/omega
        if(k==1) PB=0 else PB = PR*(1-q)/(1-omega)


        aHT = aT/(1+exp(r*(PT-hT)))
        aHB = aB/(1+exp(r*(PB-hB)))
        cH = c0/(1+exp(r*(PR-hR)))

      }else {
        aHT = aT
        aHB = aB
        cH = c0
      }

    }else {
      aHT = aT
      aHB = aB
      cH = c0
      Gs = Gw = 0
    }

    # vegetation model
    dT= aHT*R*k - dT*T
    dB= aHB*R*(1-k) -dB*B
    dR= (T+B)*cH*V - (aHT*k + aHB*(1-k))*R

    #herbivore model

    dH =  Gs + Gw - m * H


    return(list(c(dT, dR, dB, dH)))
  })
}

#_________________________________________________________________
# model vegetation with fixed herbivore
#------------------------------------------
modelH1 <- function(ti, states, parms)
{
  with(as.list(c(states, parms)), {

    V=1-T-R-B

    # preferences seedlings temperate/boreal
    epsilon = 0.0001
    omega = (T+epsilon) / (B+T+2*epsilon)
    q = p*omega / (p*omega + (1-omega)*(1-p))

    # preferences seedlings/mature
    phis_0 = (uR*R)*f / ((uR*R+uT*T)*f +uB*B*fB)
    phiw_0 = (uR*R*(1-f)) / ((uR*R + uT*T)*(1-f) + uB*B*(1-fB))
    phis = pR*phis_0 / (pR*phis_0+(1-pR)*(1-phis_0))
    phiw = pR*phiw_0 / (pR*phiw_0+(1-pR)*(1-phiw_0))

    if(H>1e-5) # population cannot recover if no reproduction is possible
    {
      Fs = (uR*R + uT*T)*f + uB*B*fB
      Fw = (uR*R + uT*T)*(1-f) + uB*B*(1-fB)

      # intakes rates
      irate_s = (taus*Fs/H) / (nus + Fs/H)
      if(irate_s>(Fs/H)) irate_s = Fs/H

      irate_w = (tauw*Fw/H) / (nuw + Fw/H)
      if(irate_w>(Fw/H)) irate_w = Fw/H

      # summer gain, reproduction limitation
      Gs = H* min(gamma, es* irate_s )

      # winter gain
      Gw = H* ew* irate_w

      # impacts
      if(R>0)
      {
        # R intake
        Us = phis * Gs/es
        Uw = phiw * Gw/ew

        PR = (Us+Uw)/(uR*R)
        if(k==0) PT=0 else PT = PR*q/omega
        if(k==1) PB=0 else PB = PR*(1-q)/(1-omega)


        aHT = aT/(1+exp(r*(PT-hT)))
        aHB = aB/(1+exp(r*(PB-hB)))
        cH = c0/(1+exp(r*(PR-hR)))

      }else {
        aHT = aT
        aHB = aB
        cH = c0
      }

    }else {
      aHT = aT
      aHB = aB
      cH = c0
      Gs = Gw = 0
    }

    # vegetation model
    dT= aHT*R*k - dT*T
    dB= aHB*R*(1-k) -dB*B
    dR= (T+B)*cH*V - (aHT*k + aHB*(1-k))*R

    return(list(c(dT, dR, dB)))
  })
}


#___________________________________________________________________
# model herbivore with fixed vegetation
#------------------------------------------

modelH <- function(ti, states, parms)
{
  with(as.list(c(states, parms)), {

    V=1-T-R-B

    if(H>1e-5)
    {

      Fs = (uR*R + uT*T)*f + uB*B*fB
      Fw = (uR*R + uT*T)*(1-f) + uB*B*(1-fB)

      # intakes rates
      irate_s = (taus*Fs/H) / (nus + Fs/H)
      if(irate_s>(Fs/H)) irate_s = Fs/H

      irate_w = (tauw*Fw/H) / (nuw + Fw/H)
      if(irate_w>(Fw/H)) irate_w = Fw/H

      # summer gain, reproduction limitation
      Gs = H* min(gamma, es* irate_s )

      # winter gain
      Gw = H* ew* irate_w


    }else {
      Gs = Gw = 0
    }

    #herbivore model

    dH =  Gs + Gw - m * H

    return(list(c(dH)))
  })
}

#______________________________________________________________________
# model vegettion only
#------------------------------------------
modelV <- function(ti, states, parms)
{
  with(as.list(c(states, parms)), {

    V=1-T-R-B

    # vegetation model
    dT= aT*R*k - dT*T
    dB= aB*R*(1-k) -dB*B
    dR= (T+B)*c0*V - (aT*k + aB*(1-k))*R

    return(list(c(dT, dR, dB)))
  })
}

#------------------------------------------------------------------
# model G : add grazing possibility
#--------------------------------------------------------------------
modelG <- function(ti, states, parms)
{
  with(as.list(c(states, parms)), {

    V=1-T-R-B

    # preferences seedlings temperate/boreal
    epsilon = 0.0001
    omega = (T+epsilon) / (B+T+2*epsilon)
    q = p*omega / (p*omega + (1-omega)*(1-p))

    # preferences seedlings/mature
    phis_0 = (uR*R)*f / ((uR*R+uT*T)*f +uB*B*fB)
    phiw_0 = (uR*R*(1-f)) / ((uR*R + uT*T)*(1-f) + uB*B*(1-fB))
    phis = pR*phis_0 / (pR*phis_0+(1-pR)*(1-phis_0))
    phiw = pR*phiw_0 / (pR*phiw_0+(1-pR)*(1-phiw_0))

    if(H>1e-5)
    {
      Fs = (uR*R + uT*T)*f + uB*B*fB
      FsG = min(uG*V, pG*Fs)
      Fs = Fs + FsG
      Fw = (uR*R + uT*T)*(1-f) + uB*B*(1-fB)

      # intakes rates
      irate_s = (taus*(Fs+FsG)/H) / (nus + (Fs+FsG)/H)
      if(irate_s>((Fs+FsG)/H)) irate_s = (Fs+FsG)/H

      irate_w = (tauw*Fw/H) / (nuw + Fw/H)
      if(irate_w>(Fw/H)) irate_w = Fw/H

      # summer gain, reproduction limitation
      Gs = H* min(gamma, es* irate_s )

      # winter gain
      Gw = H* ew* irate_w

      # impacts
      if(R>0)
      {
        # R intake
        Us = phis * (Gs/es) * (1 - FsG/Fs)
        Uw = phiw * Gw/ew

        PR = (Us+Uw)/(uR*R)
        if(k==0) PT=0 else PT = PR*q/omega
        if(k==1) PB=0 else PB = PR*(1-q)/(1-omega)


        aHT = aT/(1+exp(r*(PT-hT)))
        aHB = aB/(1+exp(r*(PB-hB)))
        cH = c0/(1+exp(r*(PR-hR)))

      }else {
        aHT = aT
        aHB = aB
        cH = c0
      }

    }else {
      aHT = aT
      aHB = aB
      cH = c0
      Gs = Gw = 0
    }

    # vegetation model
    dT= aHT*R*k - dT*T
    dB= aHB*R*(1-k) -dB*B
    dR= (T+B)*cH*V - (aHT*k + aHB*(1-k))*R

    #herbivore model

    dH =  Gs + Gw - m * H


    return(list(c(dT, dR, dB, dH)))
  })
}
