require(data.table)
require(stats4)
require(nlme)
require(mgcv)
require(gaussquad)

#============#
# PARAMETERS #
#============#
b0<-.5 # fraction of H type homes
nq<-20 # #of quadrature abscissa 

#========#
# DRIVER #
#========#
cottonmouth<-function(metro.code,metro.name) {
  
  para.names<-c(paste("lam.q",1:4,sep=""),paste("lam.b",1:2,sep=""), # buyer arrival rate
                paste("eta.q",1:4,sep=""),paste("eta.b",1:2,sep=""), # value premium
                paste("phi.q",1:4,sep=""),paste("phi.b",1:2,sep=""), # signal-to-noise ratio
                "sig")
  
  #=======================#
  # LOAD AND PREPARE DATA #
  #=======================#
  load(paste("CLEAN",metro.code,".RData",sep=""))
  
  #========================================#
  # STEP #2: MAXIMUM LIKELIHOOD ESTIMATION #
  #========================================#
  
  # demean log(Age) and log(SqFt)
  SX<-SX[,Age:=log(Age)-mean(log(Age))]
  SX<-SX[,SqFt:=log(SqFt)-mean(log(SqFt))]
  attach(SX)
  
  #==================#
  # FUNCTION: BUNDLE #
  #==================#
  bundle<-function(para,dasm) {
    
    # set parameter names
    names(para)<-para.names 
    
    # the bundle
    B<-list()               
    
    # important variables
    imp.var<-function(name,hedo,quarter,year) {
      ivar<-0
      if (quarter) {
        quarter.dummies<-para[paste(paste(name,"q",sep="."),1:4,sep="")] # quarter dummies
        ivar<-quarter.dummies[ListQuarter]
      }
      if (year) {
        year.dummies<-para[paste(paste(name,"y",sep="."),2005:2015,sep="")] # year dummies
        ivar<-ivar+year.dummies[ListYear-2004]
      }
      if (hedo) {
        ivar<-ivar+
          para[paste(name,"b1",sep=".")]*Age+   # age 
          para[paste(name,"b2",sep=".")]*SqFt   # size
      }
      ivar
    }
    
    B$lam<-imp.var("lam",TRUE,TRUE,FALSE) # buyer arrival rate
    B$eta<-imp.var("eta",TRUE,TRUE,FALSE) # value premium
    B$phi<-imp.var("phi",TRUE,TRUE,FALSE) # signal-to-noise ratio
    
    # important variables
    B$sig  <- para["sig"]               # error term
    B$lamH <- B$lam*(1-pnorm(-B$phi/2)) # rate at which H homes sell
    B$lamL <- B$lam*(1-pnorm(+B$phi/2)) # rate at which L homes sell
    
    # sale-to-list
    B$p<-function(t) (b0*B$eta*exp(-B$lamH*t)+(1-b0)*exp(-B$lamL*t))/
      ((b0*B$eta+1-b0)*(b0*exp(-B$lamH*t)+(1-b0)*exp(-B$lamL*t)))
    
    # days-on-market
    B$q<-function(t,l) dexp(t,rate=l)
    
    # return
    B
  }
  
  #======================#
  # FUNCTION: LIKELIHOOD #
  #======================#
  L<-function(para) {            
    B<-bundle(para)                            # build the bundle 
    f<-dnorm(SL,mean=B$p(DM),sd=B$sig)         # sale-to-list
    g<-b0*B$q(DM,B$lamH)+(1-b0)*B$q(DM,B$lamL) # days-on-market
    -sum(log(f*g))                             # minus log likelihood 
  }
  
  #===============#
  # INITIAL GUESS #
  #===============# 
  x0<-rep(0,length(para.names))                    # initial guess 
  names(x0)<-para.names                            # set parameter names
  x0[paste("lam.q",1:4,sep="")]<-rep(2/mean(DM),4) # buyer arrival rate (2/mean(DM))
  x0[paste("eta.q",1:4,sep="")]<-rep(1.05,4)       # home value spread (1.05)
  x0[paste("phi.q",1:4,sep="")]<-rep(1.00,4)       # signal-to-noise ratio (1.)
  x0["sig"]<-.1                                    # spherical error (.1)
  
  #===============================#
  # MAXIMUM LIKELIHOOD ESTIMATION #
  #===============================#
  print("MAXIMUM LIKELIHOOD ESTIMATION")
  fml<-optim(x0,L,hessian=TRUE,control=list(trace=1,maxit=1.e6))
  Sig<-solve(fml$hessian)                      # asymptotic variance
  FML<-cbind(fml$par,sqrt(diag(Sig)))          # point estimates
  print(FML)
  
  #===============================#
  # FUNCTION: MEAN DAYS-ON-MARKET #
  #===============================#
  mean.DM<-function(para,type,task) {
    B<-bundle(para) # compute the bundle
    dm.cntrfact<-1/B$lam
    dm.estimate<-1/((type=="H")*B$lamH+(type=="L")*B$lamL)
    switch(task,
           estimate=mean(dm.estimate),
           cntrfact=mean(dm.cntrfact),
           prctdiff=mean(dm.estimate/dm.cntrfact)-1)
  }
  
  #===========================#
  # FUNCTION: MEAN SALE PRICE #
  #===========================#
  Q<-as.matrix(laguerre.quadrature.rules(nq)[[nq]])
  mean.SP<-function(para,type,task) {
    B<-bundle(para)                      # compute the bundle
    vL<-LP/(b0*B$eta+(1-b0))             # value of L home
    vH<-B$eta*vL                         # value of H home
    lamT<-(type=="H")*B$lamH+(type=="L")*B$lamL
    sp.cntrfact<-ifelse(type=="H",mean(vH),mean(vL)) # counterfactual sale price  
    p<-function(t,j) (b0*B$eta[j]*exp(-B$lamH[j]*t)+(1-b0)*exp(-B$lamL[j]*t))/
      ((b0*B$eta[j]+1-b0)*(b0*exp(-B$lamH[j]*t)+(1-b0)*exp(-B$lamL[j]*t)))
    sp.estimate<-sapply(1:nrow(SX),function(j)Q[,2]%*%p(Q[,1]/lamT[j],j))
    sp.estimate<-mean(LP*sp.estimate)   
    switch(task,
           estimate=mean(sp.estimate),
           cntrfact=mean(sp.cntrfact),
           prctdiff=mean(sp.estimate/sp.cntrfact)-1)
  }
  
  #=================#
  # COUNTERFACTUALS #
  #=================#
  print("COUNTERFACTUALS")
  se<-function(g) sqrt(g%*%Sig%*%g/nrow(SX))
  DMCF<-c() # days-on-market table
  SPCF<-c() # sale price table
  for (type in c("H","L")) {
    DMCF<-rbind(DMCF,c(mean(DM),sd(DM)/sqrt(nrow(SX))))
    SPCF<-rbind(SPCF,c(mean(SP),sd(SP)/sqrt(nrow(SX))))
    for (task in c("estimate","cntrfact","prctdiff")) {
      print(paste(task,type))
      DMCF<-rbind(DMCF,c(mean.DM(fml$par,type,task),se(fdHess(fml$par,mean.DM,type,task)$gradient)))
      SPCF<-rbind(SPCF,c(mean.SP(fml$par,type,task),se(fdHess(fml$par,mean.SP,type,task)$gradient)))
    }
  }
  
  out<-list(fm.para = FML,                  # maximum likelihood estimates 
            fm.dmcf = DMCF,                 # days-on-market counterfactuals
            fm.spcf = SPCF,                 # sale price counterfactuals
            sx=SX[,.(LP,SP,SL,DM,Age,SqFt,Beds,Baths,Condo)])
  
  #=================#
  # MODEL FIT PLOTS #
  #=================#
  print("MODEL FIT PLOTS")
  B<-bundle(fml$par) # solution bundle
  B$tau<-90
  
  # sale-to-list
  B$p.das<-function(t) ifelse(t<=B$tau,
                             (b0*B$eta*exp(-B$lamH*t)+(1-b0)*exp(-B$lamL*t))/
                               ((b0*B$eta+1-b0)*(b0*exp(-B$lamH*t)+(1-b0)*exp(-B$lamL*t))),
                             1/(b0*B$eta+1-b0))
  
  # days-on-market
  B$q.das<-function(t,l) ifelse(t<=B$tau,
                               dexp(t,rate=l),
                               pexp(B$tau,rate=l,lower.tail=FALSE)*dexp(t-B$tau,rate=B$lam))  
  
  n<-nrow(SX)
  SX[,SL.fit:=loess(SL~DM,enp.target=19)$fitted] # sale
  SX<-SX[,.(.N,b0*mean(B$q(DM,B$lamH))+(1-b0)*mean(B$q(DM,B$lamL)),b0*mean(B$q.das(DM,B$lamH))+(1-b0)*mean(B$q.das(DM,B$lamL)),
            mean(SL.fit),mean(B$p(DM)),mean(B$p.das(DM))),by=DM]
  SX[,N:=loess(N~DM)$fitted]
  write.csv(SX,file=paste(metro.name,"FitPlot.csv",sep=""),row.names=FALSE) # write to file
  
  detach(SX)
  return(out)
}

metro.code<-readline(prompt="Enter the Metro Code:")
metro.name<-readline(prompt="Enter the Metro Name:")
X<-cottonmouth(metro.code,metro.name) 
save(X,file=paste(metro.name,"SolX.RData",sep="")) 
