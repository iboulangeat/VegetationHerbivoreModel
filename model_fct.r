
#______________________________________________________________________
# model
#------------------------------------------
model <- function(ti, states, parms)
{
  with(as.list(c(states, parms)), {
    
    G=1-T-S-B 
    
    # preferences seedlings
    epsilon = 0.0001
    wt = (T+epsilon) / (B+T+2*epsilon)
    wb = (B+epsilon) / (B+T+2*epsilon)
    qt = omega*wt / (omega*wt + (1-omega)*wb)
    qb = (1-omega)*wb / (omega*wt + (1-omega)*wb)
    # preferences mature/seedlings
    summer.ws = (uS*S)*l / ((uS*S+uT*T)*l +uB*B*lB)
    winter.ws = (uS*S*(1-l)) / ((uS*S + uT*T)*(1-l) + uB*B*(1-lB))
    phis = omega2*summer.ws / (omega2*summer.ws+(1-omega2)*(1-summer.ws))
    phiw = omega2*winter.ws / (omega2*winter.ws+(1-omega2)*(1-winter.ws))
    
    if(H>1e-5) 
    {
      Fs = (uS*S + uT*T)*l + uB*B*lB
      Fw = (uS*S + uT*T)*(1-l) + uB*B*(1-lB)
      
      # intakes rates
      nus = nu*gseason/365
      Is = (taus*gseason*Fs/H) / (nus + Fs/H)
      
      if(Is>(Fs/H)) Is = Fs/H
      
      #reproduction limitation
      gains_s = min(gmax, es* Is )
      
      nuw = nu*(1-gseason/365)
      Iw = (tauw*(365-gseason)*Fw/H) / (nuw + Fw/H)
      
      if(Iw>(Fw/H)) Iw = Fw/H
      
      gains_w = ew* Iw 
      gains = gains_s * H + gains_w * H
      
      Us = Is *H
      Uw = Iw *H
      # impacts
      if(S>0)
      {
        # S intake
        UsS = Us * phis 
        UwS = Uw * phiw 
        
        if(k==0) PT=0 else PT = (UsS+UwS)*qt/(uS*S*wt) 
        if(k==1) PB=0 else PB = (UsS+UwS)*qb/(uS*S*wb)
        PS = (UsS+UwS)/(uS*S)
        
        aHT = a0/(1+exp(r*(PT-ptresh)))
        aHB = a0/(1+exp(r*(PB-ptresh)))
        cH = c0/(1+exp(r*(PS-ptresh))) 
        
      }else {
        aHT = a0
        aHB = a0
        cH = c0
      }
      
    }else {
      aHT = a0
      aHB = a0
      cH = c0
      Is = Iw = 0
      gains = 0
      
    }        
    
    # vegetation model
    dT= aHT*S*k - dT*T
    dB= aHB*S*(1-k) -dB*B 
    dS= (T+B)*cH*G - (aHT*k + aHB*(1-k))*S
    
    #herbivore model  
    
    dH =  gains - mp * H
    
    
    return(list(c(dT, dS, dB, dH)))
  })
}

#_________________________________________________________________
# model vegetation with fixed herbivore 
#------------------------------------------
modelH1 <- function(ti, states, parms)
{
  with(as.list(c(states, parms)), {
    
    G=1-T-S-B 
    
    # preferences seedlings
    epsilon = 0.0001
    wt = (T+epsilon) / (B+T+2*epsilon)
    wb = (B+epsilon) / (B+T+2*epsilon)
    qt = omega*wt / (omega*wt + (1-omega)*wb)
    qb = (1-omega)*wb / (omega*wt + (1-omega)*wb)
    # preferences mature/seedlings
    summer.ws = (uS*S)*l / ((uS*S+uT*T)*l +uB*B*lB)
    winter.ws = (uS*S*(1-l)) / ((uS*S + uT*T)*(1-l) + uB*B*(1-lB))
    phis = omega2*summer.ws / (omega2*summer.ws+(1-omega2)*(1-summer.ws))
    phiw = omega2*winter.ws / (omega2*winter.ws+(1-omega2)*(1-winter.ws))
    
    if(H>1e-5) 
    {
      Fs = (uS*S + uT*T)*l + uB*B*lB
      Fw = (uS*S + uT*T)*(1-l) + uB*B*(1-lB)
      
      # intakes rates
      nus = nu*gseason/365
      Is = (taus*gseason*Fs/H) / (nus + Fs/H)
      
      if(Is>(Fs/H)) Is = Fs/H
      
      #reproduction limitation
      gains_s = min(gmax, es* Is )
      
      
      nuw = nu*(1-gseason/365)
      Iw = (tauw*(365-gseason)*Fw/H) / (nuw + Fw/H)
      
      if(Iw>(Fw/H)) Iw = Fw/H
      
      gains_w = ew* Iw 
      gains = gains_s * H + gains_w * H
      
      Us = Is *H
      Uw = Iw *H
      # impacts
      if(S>0)
      {
        # S intake
        UsS = Us * phis 
        UwS = Uw * phiw 
        
        if(k==0) PT=0 else PT = (UsS+UwS)*qt/(uS*S*wt) 
        if(k==1) PB=0 else PB = (UsS+UwS)*qb/(uS*S*wb)
        PS = (UsS+UwS)/(uS*S)
        
        aHT = a0/(1+exp(r*(PT-ptresh)))
        aHB = a0/(1+exp(r*(PB-ptresh)))
        cH = c0/(1+exp(r*(PS-ptresh))) 
        
      }else {
        aHT = a0
        aHB = a0
        cH = c0
      }
      
    }else {
      aHT = a0
      aHB = a0
      cH = c0
      Is = Iw = 0
    }          
    
    # vegetation model
    dT= aHT*S*k - dT*T
    dB= aHB*S*(1-k) -dB*B 
    dS= (T+B)*cH*G - (aHT*k + aHB*(1-k))*S 
    
    return(list(c(dT, dS, dB)))
  })
}


#___________________________________________________________________
# model herbivore with fixed vegetation 
#------------------------------------------

modelH <- function(ti, states, parms)
{
  with(as.list(c(states, parms)), {
    
    G=1-T-S-B 
    
    if(H>1e-5) 
    {
      
      Fs = (uS*S + uT*T)*l + uB*B*lB
      Fw = (uS*S + uT*T)*(1-l) + uB*B*(1-lB)
      
      # intakes rates
      nus = nu*gseason/365
      Is = (taus*gseason*Fs/H) / (nus + Fs/H)
      
      if(Is>(Fs/H)) Is = Fs/H
      
      #reproduction limitation
      gains_s = min(gmax, es* Is )
      
      nuw = nu*(1-gseason/365)
      Iw = (tauw*(365-gseason)*Fw/H) / (nuw + Fw/H)
      
      
      if(Iw>(Fw/H)) Iw = Fw/H
      
      
      gains_w = ew* Iw 
      gains = gains_s * H + gains_w * H
      
      
    }else {
      gains = 0
    }        
    
    
    #herbivore model  
    
    dH =  gains - mp * H
    
    
    
    return(list(c(dH)))
  })
}

#______________________________________________________________________
# model vegettion only 
#------------------------------------------
modelV <- function(ti, states, parms)
{
  with(as.list(c(states, parms)), {
    
    G=1-T-S-B 
    
    # vegetation model
    dT= a0*S*k - dT*T
    dB= a0*S*(1-k) -dB*B 
    dS= (T+B)*c0*G - (a0*k + a0*(1-k))*S
    
    return(list(c(dT, dS, dB)))
  })
}

#------------------------------------------------------------------
# model G : add grazing possibility
#--------------------------------------------------------------------
modelG <- function(ti, states, parms)
{
  with(as.list(c(states, parms)), {
    
    G=1-T-S-B 
    
    # preferences seedlings T vs B
    epsilon = 0.0001
    wt = (T+epsilon) / (B+T+2*epsilon)
    wb = (B+epsilon) / (B+T+2*epsilon)
    qt = omega*wt / (omega*wt + (1-omega)*wb)
    qb = (1-omega)*wb / (omega*wt + (1-omega)*wb)
    # preferences mature/seedlings
    summer.ws = (uS*S)*l / ((uS*S+uT*T)*l +uB*B*lB)
    winter.ws = (uS*S*(1-l)) / ((uS*S + uT*T)*(1-l) + uB*B*(1-lB))
    phis = omega2*summer.ws / (omega2*summer.ws+(1-omega2)*(1-summer.ws))
    phiw = omega2*winter.ws / (omega2*winter.ws+(1-omega2)*(1-winter.ws))
    
    
    if(H>1e-5) 
    {
      FsG = min(uG*G, pGraz*((uS*S + uT*T)*l + uB*B*lB))
      Fs = (uS*S + uT*T)*l + uB*B*lB + FsG
      Fw = (uS*S + uT*T)*(1-l) + uB*B*(1-lB)
      
      # intakes rates
      nus = nu*gseason/365
      Is = (taus*gseason*Fs/H) / (nus + Fs/H)
      IsG = (taus*gseason*FsG/H) / (nus + FsG/H)
      
      
      if(Is>(Fs/H)) Is = Fs/H
      if(IsG>(FsG/H)) IsG = FsG/H
      
      #reproduction limitation
      gains_s = min(gmax, es* (Is+IsG) )
      
      nuw = nu*(1-gseason/365)
      Iw = (tauw*(365-gseason)*Fw/H) / (nuw + Fw/H)
      
      if(Iw>(Fw/H)) Iw = Fw/H
      
      gains_w = ew* Iw 
      gains = gains_s * H + gains_w * H
      
      Us = Is *H
      UsG = IsG *H
      Uw = Iw *H
      # impacts
      if(S>0)
      {
        # S intake
        UsS = Us * phis 
        UwS = Uw * phiw
        # G intake : UsG
        
        if(k==0) PT=0 else PT = (UsS+UwS)*qt/(uS*S*wt) 
        if(k==1) PB=0 else PB = (UsS+UwS)*qb/(uS*S*wb)
        PS = (UsS+UwS)/(uS*S) - UsG/(uG*G)
        
        aHT = a0/(1+exp(r*(PT-ptresh)))
        aHB = a0/(1+exp(r*(PB-ptresh)))
        cH = c0/(1+exp(r*(PS-ptresh))) 
        
      }else {
        aHT = a0
        aHB = a0
        cH = c0
      }
      
    }else {
      aHT = a0
      aHB = a0
      cH = c0
      Is = Iw = 0
      gains = 0
      
    }        
    
    # vegetation model
    dT= aHT*S*k - dT*T
    dB= aHB*S*(1-k) -dB*B 
    dS= (T+B)*cH*G - (aHT*k + aHB*(1-k))*S
    
    #herbivore model  
    
    dH =  gains - mp * H
    
    
    return(list(c(dT, dS, dB, dH)))
  })
}
