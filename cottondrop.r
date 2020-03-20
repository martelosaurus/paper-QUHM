require(data.table)
require(lfe)

#=============#
# COTTON PATH #
#=============#
cotton.path<-function(path.type) {
  machine<-ifelse(substr(getwd(),10,17)=="joma1199","joma1199","martelosaurus")
  paste("c:/Users/",machine,"/Dropbox/ProjectLemon/Cottonmouth/",sep="")
}

#=========================#
# SUMMARY STATISTICS LOAD #
#=========================#
sum.stat.load<-function() {
  metro.names<-c("MLDENVER","CAROLINA","REALCOMP")
  M<-list()
  for (metro in metro.names) {
    load(paste("CLEAN",metro,".RData",sep=""))
    M<-c(M,list(SX))
  }
  M
}

#=============================#
# TABLE I: SUMMARY STATISTICS #
#=============================#
sum.stats<-function(X) { # DT is a list of data.tables
  metro<-function(x,type) { # loop over metros
    x<-x[,.(LP,SP,SL,DM,Age,SqFt)]
    x[,SL:=100*SL]
    if (type=="tiles") {
      x<-lapply(x,function(t) quantile(t,probs=c(.25,.50,.75))) 
      x<-lapply(x,format,digits=0,nsmall=0,big.mark=",",scientific=FALSE) # format numbers
    } else {
      x<-lapply(x,mean)    
      x<-lapply(x,format,digits=3,nsmall=1,big.mark=",",scientific=FALSE) # format numbers
    }
    # x is a list of rows
    x<-lapply(x,paste,collapse=" & ")
    

  }
  for (type in c("means","tiles")) {
    # variables 
    m1<-sapply(X,metro,type) 
    m1<-cbind(c("List Price (\\$k)","Sale Price (\\$k)","Sale-to-List (\\%)","Days-on-Market","Age (years)","Square Feet"),m1)
    m1<-paste(apply(m1,1,paste,collapse=" & "),collapse="\\\\ \n")
    # observations
    nobs<-paste("\\multicolumn{",ifelse(type=="tiles",3,1),"}{c}{",format(sapply(X,nrow),big.mark=","),"}",sep="")
    m2<-paste(c("\\#Listings",nobs),collapse=" & ")
    cat(m1,m2,sep="\\\\ \\hline \n",file=paste(cotton.path(""),"sumstats",type,".tex",sep=""))
  }
}

#===============================#
# FUNCITON: FORMAT COEFFICIENTS #
#===============================#
# returns a formatted vector of the form 
# (p.e.^stars,(s.e.),p.e.^stars,(s.e.),...,p.e.^stars,(s.e.))
format.coef<-function(m) { # m is a matrix [point estimates, standard errors]
  x<-cbind(m[,1:2],2*pt(-abs(m[,1]/m[,2]),Inf)) # p-values
  x<-round(x,digits=3)
  comp.stars<-function(t)paste(rep("*",findInterval(-t,-c(1,.05,.025,.005))-1),collapse="")
  coef.stars<-unlist(lapply(x[,3],comp.stars))
  x<-format(x,nsmall=3,trim=TRUE)
  x[,1]<-paste(x[,1],"$^{",coef.stars,"}$",sep="") # add stars to coefficients
  x[,2]<-paste("(",x[,2],")",sep="") # add parentheses to standard errors
  x<-as.vector(t(x[,1:2])) # reshape
}

#===============================#
# TABLE II: HEDONIC REGRESSIONS #
#===============================#
hedo.regs<-function(FM) { # FM is a list of felm objects
  
  # functions: observations and r-squared
  obs<-function(fm) sum(summary(fm)$df[1:2])
  ar2<-function(fm) summary(fm)$adj.r.squared
  
  # function: format.coef wrapper
  coef.wrapper<-function(fm) format.coef(summary(fm)$coefficients)
  
  # coefficient names
  coef.names<-list(c("Days-on-Market","",
                     "$\\log(\\text{Age})$","",
                     "$\\log(\\text{Square Feet})$",""))
  
  # coefficients for paper
  M1<-lapply(FM,coef.wrapper)                                    # coefficients and standard errors
  M1<-do.call(cbind,c(coef.names,M1))                            # add the coefficient names
  M1<-paste(apply(M1,1,paste,collapse=" & "),collapse="\\\\ \n") # collapse rows, then columns
  
  # coefficients for deck
  M1.deck<-lapply(FM[c(2,4,6)],coef.wrapper)                                     # coefficients and standard errors
  M1.deck<-do.call(cbind,c(list(c("Days-on-Market","")),M1.deck))                      # add the coefficient names
  M1.deck<-paste(apply(M1.deck[1:2,],1,paste,collapse=" & "),collapse="\\\\ \n") # collapse rows, then columns

  # dummies, fixed effects, r-squared, and number of observations
  M2<-c(paste(c("Quarter",rep(c("$\\times$",""),length(FM)/2)),collapse=" & "),
        paste(c("ZIP$\\times$Year$\\times$Quarter",rep(c("","$\\times$"),length(FM)/2)),collapse=" & "),
        paste(c("Adjusted R$^{2}$",paste(format(100*sapply(FM,ar2),digits=3),"\\%",sep="")),collapse=" & "),
        paste(c("\\#Listings",format(sapply(FM,obs),big.mark=",")),collapse=" & "))
  
  M2.deck<-c(paste(c("Age and Size Controls",rep("$\\times$",length(FM)/2)),collapse=" & "),
             paste(c("ZIP$\\times$Year$\\times$Quarter",rep("$\\times$",length(FM)/2)),collapse=" & "),
             paste(c("Adjusted R$^{2}$",paste(format(100*sapply(FM[c(2,4,6)],ar2),digits=3),"\\%",sep="")),collapse=" & "),
             paste(c("\\#Listings",format(sapply(FM[c(2,4,6)],obs),big.mark=",")),collapse=" & "))
  
  # write paper table to file
  cat(M1,paste(M2,collapse="\\\\ \n"),sep="\\\\ \\hline\n",file=paste(cotton.path(""),"hedo.tex",sep=""))  
  
  # write deck table to file
  cat(M1.deck,paste(M2.deck,collapse="\\\\ \n"),sep="\\\\ \\hline\n",file=paste(cotton.path(""),"hedodeck.tex",sep=""))  
  
}

#==========================#
# TABLE III: MLE ESTIMATES #
#==========================#
mle.para<-function(FM) { # FM is a list of matrices [pe, se] 

  # reshape the formatted matrices
  M<-do.call(cbind,lapply(FM,format.coef))    
  
  # load age and size means
  X<-sum.stat.load() # X is a list of data.tables 
  me<-do.call(cbind,lapply(X,function(x)x[,c(mean(Age),mean(SqFt))]))
  me[1,]<-10/me[1,]
  me[2,]<-400/me[2,]
  
  # build mle table
  mle.tab<-function(type,dest) {
    
    # indices and names
    if (type=="controls") {
      coef.names<-c("$\\log(\\text{Age})$","","$\\log(\\text{SqFt})$","")
      I<-list(lam=9:12,eta=21:24,phi=33:36) # age and size (for eta and phi)
    } else if (type=="dummies") {
      coef.names<-c("1st Quarter","","2nd Quarter","","3rd Quarter","","4th Quarter","")
      I<-list(lam=1:8,eta=13:20,phi=25:32)  # quarter dummies (for eta, phi, and lam)
    } else if (type=="sigma") {
      coef.names<-c("Constant","")
      I<-list(sig=37:38)
    }
    
    # build table
    for (k in 1:length(I)) { # loop over variables
      
      m1<-M[I[[k]],]                                                 # peel-off coefficients
      m1<-cbind(coef.names,m1)                                       # bind coefficient names
      m1<-paste(apply(m1,1,paste,collapse=" & "),collapse="\\\\ \n") # collapse rows, then columns
      
      m2<-c() # quarter dummies label
      m3<-c() # magnitudes label 
      m4<-c() # magnitudes
      
      # controls table in deck
      if (type=="controls"&&dest=="deck") {
        
        # labels
        m2<-paste(paste(c(paste("\\hyperlink{",names(I)[k],"dummiesdeck}{Quarter Dummies}",sep=""),
                          rep("$\\times$",3)),collapse=" & "),collapse="\\\\ \n")
        m3<-"& \\multicolumn{3}{c}{\\textbf{Magnitudes}}"
        
        # magnitudes
        m4<-lapply(FM,function(fm)fm[paste(rep(names(I)[k],each=2),c("b1","b2"),sep="."),1]) # coefficients
        m4<-lapply(1:3,function(j)format(round(m4[[j]]*me[,j],digits=3),nsmall=3,trim=TRUE))
        qd<-lapply(FM,function(fm)fm[paste(rep(names(I)[k],each=4),c("q1","q2","q3","q4"),sep="."),1])
        qd<-lapply(1:3,function(j)format(round(mean(qd[[j]][X[[j]][,ListQuarter]]),digits=3),nsmall=3,trim=TRUE))
        print(qd)
        m4<-lapply(1:3,function(j)c(qd[[j]],m4[[j]]))
        m4<-do.call(cbind,c(list(c("Mean Level","+10 Years","+400 SqFt")),m4))
        m4<-paste(apply(m4,1,paste,collapse=" & "),collapse="\\\\ \n")
      } 
      
      # write to file
      cat(m1,m2,m3,m4,sep="\\\\ \\hline\n",file=paste(cotton.path(""),names(I)[k],type,dest,".tex",sep="")) 
    }
  }
  
  for (j in c("controls","dummies","sigma")) {
    for (k in c("paper","deck")) {
      mle.tab(j,k)
    }
  }
}

#==================================#
# TABLES IV and V: COUNTERFACTUALS #
#==================================#
counter.fact<-function(M,file.name) { # M: a list of matrices [point estimate, standard error]
  
  # coefficient names
  coef.names<-c("Mean","","Estimated","","Counterfactual","","Percent Difference","")
  
  M1<-do.call(cbind,lapply(M,format.coef))    # reshape the formatted matrices
  M1<-list(M1[1:8,],M1[9:16,])
  for (k in 1:2) {
    m<-M1[[k]]
    m<-cbind(coef.names,m)
    m<-paste(apply(m,1,paste,collapse=" & "),collapse="\\\\ \n") # collapse rows, then columns
    cat(m,sep="\\\\ \\hline\n",file=paste(cotton.path(""),file.name,k,".tex",sep="")) # write to file
  }
}

#==================================#
# DRIVER, FIGURE, AND TABLES CALLS #
#==================================#

# SUMMARY STATISTICS 
sum.stats(sum.stat.load()) 

# HEDONIC REGRESSIONS
metro.names<-c("Denver","Charlotte","Detroit")
M<-list()
for (metro in metro.names) {
  load(paste(metro,"Hedo.RData",sep=""))
  M<-c(M,H)
}
hedo.regs(M)

# MAXIMUM LIKELIHOOD ESTIMATION
metro.names<-c("Denver","Charlotte","Detroit")
FM<-list()
for (metro in metro.names) {
  load(paste(metro,"Sol.RData",sep=""))
  FM<-c(FM,list(X))
}

# DISPLAY MEAN RESULTS
print(matrix(apply(sapply(FM,function(fm)fm$fm.dmcf[,1]),1,mean),nrow=4))
print(matrix(apply(sapply(FM,function(fm)fm$fm.spcf[,1]),1,mean),nrow=4))

# BUYER ARRIVAL RATE PLOT
write.csv(do.call(cbind,lapply(FM,function(fm)fm$fm.para[1:4,1])),file="BuyerArrivalPlot.csv",row.names=FALSE)

# DAYS-ON-MARKET COUNTERFACTUALS
write.csv(do.call(cbind,lapply(FM,function(fm)fm$fm.dmcf[c(4,8),1])),file="dmcfplot.csv",row.names=FALSE)

# SALE PRICE COUNTERFACTUALS
write.csv(do.call(cbind,lapply(FM,function(fm)fm$fm.spcf[c(4,8),1])),file="spcfplot.csv",row.names=FALSE)

# TABLES
mle.para(lapply(FM,function(x)x$fm.para))            # maximum likelihood estimates
counter.fact(lapply(FM,function(x)x$fm.dmcf),"dmcf") # days-on-market counterfactuals
counter.fact(lapply(FM,function(x)x$fm.spcf),"spcf") # sale price counterfactuals
