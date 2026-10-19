#Total Clusters
nCluster=3

#Library
library(rjags)
library(R2jags)
library(jagsUI)
library(MASS)
library(base)
library(stats)
library(dplyr)
library(data.table)
library(survival)
library(splines)
library(lme4)


#-------------------------------------------------------------------------------------------#
#File Path
fp="C:/Users/cssc/OneDrive/Desktop/Experiment JM/EXP_JAGS"

#Input path
path=paste(fp,"/Cluster_JM/Exp_work/DATA",sep="")
#Output path
path_out=paste(fp,"/JAGS_Output/Cluster_JM",sep="")
dir.create(path_out)


#Complete Data
data_org=readRDS(paste(path,"/data_0.rds",sep=""))

#Making every id to the same starting point.
min_w=data_org%>%
  group_by(Child.Id)%>%
  summarise(min_week=min(Week)*277)
data_org=merge(x=data_org, y=min_w,by=c("Child.Id"),all.x=T)

data_org$Week=data_org$Week*277-data_org$min_week+1
data_org$Last_follow_up_date=data_org$Last_follow_up_date*277-data_org$min_week+1
data_org$Time_to_event=data_org$Time_to_event*277-data_org$min_week+1


#Removing columns
data_org=data_org[,!colnames(data_org)%in%c("Height","Weight","BSA","min_week")]

data_org$id_new=as.numeric(data_org$id_new)
data_org$time_id=as.numeric(data_org$time_id)
data_org=data_org[
  with(data_org, order(id_new, time_id)),
]

#Removing all NA's  LYMPH= WBC-ANC
data_org=subset(data_org, !is.na(data_org$WBC-data_org$ANC))
#Keeping id's with > 5 time points  NA and updating the id names
id_set=unique(data_org$Child.Id)
dat=NULL
for(i in 1:length(id_set))
{
  d=subset(data_org, data_org$Child.Id==id_set[i])
  dat=rbind(dat, c(id_set[i],nrow(d), round(length(d$WBC[is.na(d$WBC)])*100/nrow(d),2)))
  i=i+1
}
dat=data.frame(dat)
rownames(dat)=1:nrow(dat)
dat$X2=as.numeric(dat$X2)
dat$X3=as.numeric(dat$X3)
dat1=subset(dat, dat$X3<=105)
dat1=subset(dat1, dat1$X2>=5)
data_org1=subset(data_org, data_org$Child.Id%in%unique(dat1$X1))
unique(data_org1$id_new)

dat1=data_org1[,c("Child.Id","id_new")]
dat1=dat1[!duplicated(dat1), ]
dat1$id_new=sapply(dat1$Child.Id,function(x)as.numeric(strsplit(x,"_")[[1]][2]))
dat1=dat1[with(dat1, order(id_new)), ]
dat1$id_new=1:nrow(dat1)
data_org1=merge(x=data_org1[,!colnames(data_org1)%in%c("id_new")], y=dat1, by=c("Child.Id"), all.x=T)
data_org1=data_org1[,c(1,26,2:25)]
data_org=data_org1

#Scaling the time (max week is now 271 weeks instead of previous 277)
max_week=max(max(data_org$Week),max(data_org$Last_follow_up_date),
             max(data_org$Time_to_event))
data_org$Week=data_org$Week/max_week
data_org$Time_to_event=data_org$Time_to_event/max_week
data_org$Last_follow_up_date=data_org$Last_follow_up_date/max_week

#----------------------------------------------------------------------------------------------#
#                                 Covariate transformation
#----------------------------------------------------------------------------------------------#

#Medicine (transformed --> untransformed )
#-------------------------------------------#
data_org$X6_mp=data_org$X6_mp*275
data_org$MTX=data_org$MTX*12.5

#Normalizing medicine
#---------------------#
data_org$X6_mp= (data_org$X6_mp)/median(data_org$X6_mp)
data_org$MTX= (data_org$MTX)/median(data_org$MTX)


#Normalizing fixed covariates
#-----------------------------#
data_org$WBC_high=log(exp(data_org$WBC_high)*21250)
data_org$WBC_high=data_org$WBC_high-median(unique(data_org$WBC_high))


data_org$Age=data_org$Age*4.280008 + 6.589453
data_org$Age=(data_org$Age-mean(unique(data_org$Age)))/sd(unique(data_org$Age))

#----------------------------------------------------------------------------------------------#
#                                 Response transformation
#----------------------------------------------------------------------------------------------#


#Center response y to 0 (untransform --> center to 0)
#----------------------------------------------------#
data_org$WBC=data_org$WBC +log(3200)
data_org$ANC=data_org$ANC+log(1564)

data_org$LYMPH=exp(data_org$WBC)-exp(data_org$ANC)
data_org$LYMPH=log(data_org$LYMPH)
#data_org$WBC=0
#data_org$WBC=data_org$WBC_new

data_org=data_org[,!colnames(data_org)%in%c("WBC_new")]
data_org$Platelet=data_org$Platelet+log(192000)

#Centering to 0
#---------------#
data_org$WBC=data_org$WBC-median(data_org$WBC[!is.na(data_org$WBC)])
data_org$ANC=data_org$ANC-median(data_org$ANC[!is.na(data_org$ANC)])
data_org$Platelet=data_org$Platelet-median(data_org$Platelet[!is.na(data_org$Platelet)])
data_org$LYMPH=data_org$LYMPH-median(data_org$LYMPH[!is.na(data_org$LYMPH)])

#Exclude variables Day 8 , Day 35 risk stratification and Status
#---------------#
data_org=data_org[,!colnames(data_org)%in%c("Day8","Day35","Status")]

#Getting intercept
data_org$Intercept=1

#saveRDS(data_org,paste(path,"/Transformed_DATA.rds",sep=""))


#----------------------------------------------------------------------------------------------#

#number of ids 
#-----------------------------------------------------------------------------------------------------
nId <- length(unique(data_org$id_new))

#Response
#-----------------------------------------------------------------------------------------------------#
y <- as.matrix(data_org[,c("id_new","time_id","WBC","ANC","Platelet","LYMPH")])

#------------------------------------------------------------------------------------------------------#

#Fixed covariates
X_fix=data_org[,c("id_new","Intercept","Age","Gender","Lineage","WBC_high","NCI","Bulky","CNS","Risk","Risk_present","Morph_CR1")]
X_fix=X_fix[!duplicated(X_fix$id_new),]
X_fix=X_fix[with(X_fix, order(id_new)), ]
rownames(X_fix)=NULL


#------------------------------------------------------------------------------------------------------#

#X origin
X <- data_org[,c("id_new","time_id","Intercept","Week","X6_mp","MTX")]
X$Week2<-X$Week^2
X<-X[,c("id_new","time_id","Intercept","Week","Week2","X6_mp","MTX")]
X<-merge(x=X, y=X_fix[,!colnames(X_fix)%in%c("Intercept")], by=c("id_new"), all.x=T)
X=X[with(X, order(id_new,time_id)), ]
rownames(X)<-NULL
X<-as.matrix(X)


#medicine covariate
medi=function(id, t)
{
  d=subset(data_org, data_org$id_new==id)
  d=d[with(d, order( time_id)),]
  z=NULL
  if(t <=max(d$Week) & t>=min(d$Week))
  {
    for(j in 1:(nrow(d)-1))
    {
      if(t>=d[j,c("Week")] & t<d[(j+1),c("Week")])
      {
        z = d[j,c("X6_mp","MTX")]
      }
    }
    if(t==max(d$Week))
    {
      z=d[nrow(d),c("X6_mp","MTX")]
    }
  }
  else{
    z=c(0,0)
  }
  return(as.numeric(z))
}

#----------------------------------------------------------------------------------------------#

#G-K points (Gauss-Konrod points)
Gk_points=as.matrix(cbind(c(-0.949107912342758524526189684047851, -0.741531185599394439863864773280788, -0.405845151377397166906606412076961, 0,
                            0.405845151377397166906606412076961, 0.741531185599394439863864773280788, 0.949107912342758524526189684047851, -0.991455371120812639206854697526329,
                            -0.864864423359769072789712788640926, -0.586087235467691130294144838258730, -0.207784955007898467600689403773245, 0.207784955007898467600689403773245,
                            0.586087235467691130294144838258730, 0.864864423359769072789712788640926, 0.991455371120812639206854697526329),
                          c(0.063092092629978553290700663189204, 0.140653259715525918745189590510238, 0.190350578064785409913256402421014,
                            0.209482141084727828012999174891714, 0.190350578064785409913256402421014, 0.140653259715525918745189590510238, 0.063092092629978553290700663189204,
                            0.022935322010529224963732008058970, 0.104790010322250183839876322541518, 0.169004726639267902826583426598550, 0.204432940075298892414161999234649,
                            0.204432940075298892414161999234649, 0.169004726639267902826583426598550, 0.104790010322250183839876322541518, 0.022935322010529224963732008058970)))

colnames(Gk_points)=c("cuts","weights")
Gk_points=data.frame(Gk_points)
Gk_points= Gk_points[with(Gk_points, order(cuts)),]
rownames(Gk_points)=1:nrow(Gk_points)

#G-K cuts
#sk
sk <- Gk_points$cuts


#G-K weights
#wk
wk <- Gk_points$weights
wk_r <- rep(Gk_points$weights, nId)


#-------------------------------------------------------------------------------------------------------#

#For survival part 1
#-------------------
XT=data_org[,c("id_new","Intercept","Time_to_event")]
XT=XT[!duplicated(XT$id_new),]
XT=XT[with(XT, order(id_new)), ]
rownames(XT)<-NULL
XT$X6_mp=NA
XT$MTX=NA

for( i in 1:nrow(XT))
{
  v=as.numeric(medi(as.numeric(XT[i,"id_new"]),as.numeric(XT[i,"Time_to_event"])))
  XT[i,c("X6_mp","MTX")]=v
  i=i+1
}
XT$time_id=0
XT$Time_to_event2=XT$Time_to_event^2
XT<-XT[,c("id_new","time_id","Intercept","Time_to_event","Time_to_event2","X6_mp","MTX")]
XT<-merge(x=XT, y=X_fix[,!colnames(X_fix)%in%c("Intercept")], by=c("id_new"), all.x=T)
XT=XT[with(XT, order(id_new,time_id)), ]
rownames(XT)<-NULL
XT<-as.matrix(XT)



#For survival part 2
#-------------------
XS=data_org%>%
  group_by(id_new)%>%
  summarise(Intercept=1,min_week=min(Week), Time_to_event=Time_to_event[1])
XS=XS[c(rep(1:nId,each=15)),]
XS$time_id=rep(1:15, nId)
XS=XS[with(XS, order(id_new)), ]
rownames(XS)<-NULL
XS$X6_mp=NA
XS$MTX=NA


XS$Week=(((Gk_points[XS$time_id,"cuts"]+1)*0.5 )*(XS$Time_to_event-XS$min_week))+XS$min_week
XS=XS[,c("id_new","time_id","Intercept","Week","X6_mp","MTX")]
for( i in 1:nrow(XS))
{
  v=as.numeric(medi(as.numeric(XS[i,"id_new"]),as.numeric(XS[i,"Week"])))
  XS$X6_mp[i]=v[1]
  XS$MTX[i]=v[2]
  
  i=i+1
}
XS$Week2=XS$Week^2
XS<-XS[,c("id_new","time_id","Intercept","Week","Week2","X6_mp","MTX")]
XS<-merge(x=XS, y=X_fix[,!colnames(X_fix)%in%c("Intercept")], by=c("id_new"), all.x=T)
XS=XS[with(XS, order(id_new,time_id)), ]
rownames(XS)<-NULL
XS<-as.matrix(XS)

#----------------------------------------------------------------------------------------------#
#Survival responses
y_s=data_org[,c("id_new","Relapsed","Time_to_event")]
y_s=y_s[!duplicated(y_s$id_new),]
y_s=y_s[with(y_s, order(id_new)), ]
rownames(y_s)=NULL





#------------------------------------------------------------------------------------------------------#

#Offset for longitudinal part of the model
total_inputs=data_org%>%
  group_by(id_new)%>%
  summarise(n_i=n())
offset=c(0,cumsum(total_inputs$n_i))+1

#Offset for Survival part
offset_s=(0:(nId))*15+1


#------------------------------------------------------------------------------------------------------#
#Constant
C=10^7
#------------------------------------------------------------------------------------------------------#
#y polynomial
y_polynomial<-as.matrix(data_org[,c("id_new","time_id","Intercept","Week")])
#The Baseline Hazard
survObject <- coxph(Surv( time =Time_to_event ,as.numeric(Relapsed)) ~id_new , data = y_s , x=TRUE , model=TRUE)
SurvInf <- survObject$y
typeSurvInf <-attr(SurvInf, "type")

#for right censored
if(typeSurvInf=="right")
{
  Time <- SurvInf[,"time"]
  Time[Time < 1e-04] <-1e-04
  nT <- length(Time)
  event <-SurvInf[,"status"]
  LongFormat <- FALSE
}


#P
P <-0.5*(Time-y_polynomial[offset[c(-length(offset))],"Week"])


st <- outer(P , (sk+1))+ matrix(rep(y_polynomial[offset[c(-length(offset))],"Week"],each=length(sk)),byrow=T,ncol=length(sk))
tpower <- function(x,t,p) ((x-t)^p)*as.numeric(x>t)
bbase <- function(x,xl,xr,ndx,deg)
{
  dx <- (xr-xl)/ndx
  kn=tail(head(seq(xl, xr, length.out = (ndx + 1)), -1), -1)
  kn <- kn[kn < max(Time)]
  knots<- sort(c(rep(range(Time, st), (deg + 1)), kn)) #Extended knots
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n) , diff = deg +1)/gamma((deg + 1)*dx^deg)
  B <- (-1)^(deg+1) *P%*%t(D)
  B
}

#P splines of degree = deg with (ndx +1) equispaced knots on [xl,xr]
BS_T <-bbase(x=Time ,xl=quantile(Time, probs =c(0.05),names=F) , xr=quantile(Time, probs =c(0.95),names=F) , ndx=14 , deg=3)

BS_S <- bbase(x=c(t(st)) ,xl=quantile(Time, probs =c(0.05),names=F) , xr=quantile(Time, probs =c(0.95),names=F) , ndx=14 , deg=3)


DD <-diag(ncol(BS_T))
M.gam.h0 <-crossprod(diff(DD,differences =2))+1e-06*DD

#------------------------------------------------------------------------------------------------------#
zeros=rep(0,nId)
ze=rep(0,6)



#------------------------------------------------------------------------------------------------------#
#                                         The composite cluster model                                  #
#------------------------------------------------------------------------------------------------------#



sink(paste(path_out,"JM_comp_Cl.model.txt", sep="/"))
cat("
model {

for(i in 1:nId) 
{
 
 #Random Effect
 #-------------
 for(cl in 1:nCluster)
  {
   u[i,1:6,cl] ~ dmnorm(zeros[1:6],D_inv_L[1:6,1:6,cl]) #Cluster wise random intercept and slope for all 3 responses
  }

 #Longitudinal Process
 #---------------------
for( j in offset[i]:(offset[i+1]-1))
 {
  for(cl in 1:nCluster)
  {
    #Mean Longitudinal
    mu_LYM[j,cl]= u[i,1,cl] + u[i,4,cl]*X[j,4] + inprod(beta_LYM[1:(ncX-2),cl], X[j,3:ncX]) 
    mu_ANC[j,cl]= u[i,2,cl] + u[i,5,cl]*X[j,4] + inprod(beta_ANC[1:(ncX-2),cl], X[j,3:ncX])
    mu_PLT[j,cl]= u[i,3,cl] + u[i,6,cl]*X[j,4] + inprod(beta_PLT[1:(ncX-2),cl], X[j,3:ncX])
    
  #Log density cluster wise of LYM, ANC, PLT
  ld[j,cl] <- (logdensity.norm(y[j,6], mu_LYM[j,cl], prec_LYM[cl]) +
               logdensity.norm(y[j,4], mu_ANC[j,cl], prec_ANC[cl]) +
               logdensity.norm(y[j,5], mu_PLT[j,cl], prec_PLT[cl]) )
  }
  } #loop of j
  
   for(cl in 1:nCluster)
  {
  
  #Survival part 1
  #---------------
   mu_TLYM[i,cl] =u[i,1,cl] + u[i,4,cl]*XT[i,4]  + inprod(beta_LYM[1:(ncXT-2),cl], XT[i,3:ncXT])
   mu_TANC[i,cl] =u[i,2,cl] + u[i,5,cl]*XT[i,4]  + inprod(beta_ANC[1:(ncXT-2),cl], XT[i,3:ncXT])
   mu_TPLT[i,cl] =u[i,3,cl] + u[i,6,cl]*XT[i,4]  + inprod(beta_PLT[1:(ncXT-2),cl], XT[i,3:ncXT])

   log.hT[i,cl]= (inprod(gam.h0[1:17,cl],BS_T[i,1:17]) #Baseline log hazard
                  + psi_LYM[cl]*mu_TLYM[i,cl]
                  + psi_ANC[cl]*mu_TANC[i,cl] 
                  + psi_PLT[cl]*mu_TPLT[i,cl]
                  + inprod(theta[1:(ncF-2),cl],X_fix[i,3:ncF])
                  ) 

 #Survival part 2
 #---------------
 for(l in offset_s[i]:(offset_s[i+1]-1) )
  {
    mu_SLYM[l,cl] =u[i,1,cl]  + u[i,4,cl]*XS[l,4] + inprod(beta_LYM[1:(ncXS-2),cl], XS[l,3:ncXS])
    mu_SANC[l,cl] =u[i,2,cl]  + u[i,5,cl]*XS[l,4] + inprod(beta_ANC[1:(ncXS-2),cl], XS[l,3:ncXS])
    mu_SPLT[l,cl] =u[i,3,cl]  + u[i,6,cl]*XS[l,4] + inprod(beta_PLT[1:(ncXS-2),cl], XS[l,3:ncXS])
                                             
    SurvLong[l,cl] = wk_r[l]* ( exp( (inprod(gam.h0[1:17,cl],BS_S[l,1:17]) #Baseline log hazard
                                  + psi_LYM[cl]*mu_SLYM[l,cl] 
                                  + psi_ANC[cl]*mu_SANC[l,cl] 
                                  + psi_PLT[cl]*mu_SPLT[l,cl]) ) )
  }
   #Baseline covariates in Survival part 2
   #--------------------------------------
   surv_base[i,cl] = exp(inprod(theta[1:(ncF-2),cl],X_fix[i,3:ncF]))
   
   
    #The likelihood clusterwise
    #--------------------------
    lik[i,cl] <-  exp((logdensity.mnorm(u[i,1:6,cl], zeros[1:6], D_inv_L[1:6,1:6,cl]) #Random Effect of LYM, ANC, PLT
                       +sum(ld[offset[i]:(offset[i+1]-1),cl])   #Longitudinal LYM, ANC, PLT
                       +(y_s[i,2]*log.hT[i,cl])   #Hazard
                       -( P[i]*surv_base[i,cl]*sum(SurvLong[offset_s[i]:(offset_s[i+1]-1),cl]) ) )) #Survival function

   #Cluster wise likelihood times their cluster probs
   #-------------------------------------------------
   pi_st[i,cl] = lik[i,cl]*pi[cl]
   
  }
 
  #Zeros trick in JAGS for JM
  #---------------------------------
  phi_bern[i] = inprod(lik[i,1:nCluster],pi[1:nCluster])/Const
  ones[i] ~ dbern(phi_bern[i])
  
  #Get cluster membership
  #-----------------------
  pi_star[i,1:nCluster] =pi_st[i,1:nCluster]/sum(pi_st[i,1:nCluster]) # Posterior Cluster probabilities
  clus[i] ~ dcat(pi_star[i,1:nCluster]) #Cluser membership based on pi_star

}#loop of i


  
#Prior
#-------------------------------------------------------------------------


#pi with dirichlet prior
#-----------------------
pi[1:nCluster] ~ ddirch(alpha[1:nCluster])
 

  for(c in 1:nCluster)
  {  

  #Error Variance
  #--------------
  prec_LYM[c] ~ dgamma(0.001,0.001)
  sig_LYM[c] = 1/prec_LYM[c]
  
  prec_ANC[c] ~ dgamma(0.001,0.001)
  sig_ANC[c] = 1/prec_ANC[c]
  
  prec_PLT[c] ~ dgamma(0.001,0.001)
  sig_PLT[c] = 1/prec_PLT[c]
  

  #Variance for Random Effect of LYM, ANC, PLT 
  #-------------------------------------------
  tau[1:6,1,c]~dmnorm(zeros[1:6],(Omega[1:6,1:6]))
  tau[1:6,2,c]~dmnorm(zeros[1:6],(Omega[1:6,1:6]))
  tau[1:6,3,c]~dmnorm(zeros[1:6],(Omega[1:6,1:6]))
  tau[1:6,4,c]~dmnorm(zeros[1:6],(Omega[1:6,1:6]))
  tau[1:6,5,c]~dmnorm(zeros[1:6],(Omega[1:6,1:6]))
  tau[1:6,6,c]~dmnorm(zeros[1:6],(Omega[1:6,1:6]))
  tau[1:6,7,c]~dmnorm(zeros[1:6],(Omega[1:6,1:6]))

  
  D_inv[1:6,1:6,c]<-tau[1:6,1:7,c] %*% t(tau[1:6,1:7,c]) + Omega[1:6,1:6]*0.001

  
  #Beta
  #-----
  for(h in 1:(ncX-2))
  {
  beta_LYM[h,c]~dnorm(0,0.001)
  beta_ANC[h,c]~dnorm(0,0.001)
  beta_PLT[h,c]~dnorm(0,0.001)
  }
  
  #Association parameter
  #----------------------
  psi_LYM[c] ~ dnorm(0,0.001)
  psi_ANC[c] ~ dnorm(0,0.001)
  psi_PLT[c] ~ dnorm(0,0.001)
 
  #Theta
  #-----
  for(h in 1:(ncF-2))
  {
  theta[h,c]~dnorm(0,0.001)
  }
  
  #Prior for penalized B-spline basis coefficients
  #------------------------------------------------
  gam.h0[1:17,c] ~  dmnorm(priorMean.gam.h0[],tau.gam.h0[c]*M.gam.h0[,])
  
  #Smoothing parameter
  #-------------------
  tau.gam.h0[c] ~ dgamma(a.gam.h0,b.gam.h0)
  
  }#End of cluster parameters
  
  
}
",fill = TRUE)
sink()

#Data list for Lagged JM
dat_list_JM =list('nId'=nId,'nCluster'=nCluster, 'ones'=rep(1,nId), 'zeros'=rep(0,6) ,'Const'=10^9, 
                  'y'=y,'y_s'=y_s, 'offset'=offset,'offset_s'=offset_s , 'BS_S'=BS_S , 'BS_T'=BS_T ,
                  'X'=X,'ncX'=ncol(X),'XS'=XS,'ncXS'=ncol(XS),'XT'=XT,'ncXT'=ncol(XT),'X_fix'=X_fix, 'ncF'=ncol(X_fix), 
                  'P'=P , 'wk_r'=wk_r , 
                  'alpha'=c(1,1,1,1,1), 
                  'Omega_AP'=diag(rep(1,4)) ,'Omega'=diag(rep(1,2)) ,
                  'priorMean.gam.h0'=seq(-8.5,0,len =ncol(BS_S)) , 'M.gam.h0'=M.gam.h0 ,
                  'a.gam.h0'=1 , 'b.gam.h0'=0.05)


Params_save=c('u', #Random effects
              'beta_LYM','beta_ANC','beta_PLT',  #The betas
              'sig_LYM', 'sig_ANC' ,'sig_PLT', #Error variances
              'D_inv', #Random effect covariance matrix inverse
              'psi_LYM','psi_ANC','psi_PLT','theta','gam.h0','tau.gam.h0', #Elements of Hazard
              'pi','pi_star','clus', #Cluster properties
              'lik' #Cluster wise likelihood
)

Sys.time()
#Run Lagged JM
jagsfit_JM <- jags(data=dat_list_JM, 
                   inits=NULL, 
                   parameters.to.save=Params_save,
                   model.file=paste(path_out,'JM_comp_Cl.model.txt', sep="/"),parallel = TRUE,
                   n.chains=1, n.iter=20000, n.burnin=5000,
                   n.thin=1,n.adapt = 5000,
                   DIC=TRUE)


saveRDS(jagsfit_JM,paste(path_out,paste("/JM_comp_Cl_init0_",nCluster,".rds",sep=""),sep=""))