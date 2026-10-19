#Total Clusters
nCluster=2

#Degree of polynomial
deg=3
poly_removecov=c("Week","Week2","Week3")[-c(1:deg)]
poly_removecovT=c("Time_to_event","Time_to_event2","Time_to_event3")[-c(1:deg)]


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
fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"

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
X$Week3<-X$Week^3
X<-X[,c("id_new","time_id","Intercept","Week","Week2","Week3","X6_mp","MTX")]
X<-merge(x=X, y=X_fix[,!colnames(X_fix)%in%c("Intercept")], by=c("id_new"), all.x=T)
X=X[with(X, order(id_new,time_id)), ]
X<-X[,!colnames(X)%in%c(poly_removecov)]
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
XT$Time_to_event3=XT$Time_to_event^3
XT<-XT[,c("id_new","time_id","Intercept","Time_to_event","Time_to_event2","Time_to_event3","X6_mp","MTX")]
XT<-merge(x=XT, y=X_fix[,!colnames(X_fix)%in%c("Intercept")], by=c("id_new"), all.x=T)
XT=XT[with(XT, order(id_new,time_id)), ]
XT<-XT[,!colnames(XT)%in%c(poly_removecovT)]
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
XS$Week3=XS$Week^3
XS<-XS[,c("id_new","time_id","Intercept","Week","Week2","Week3","X6_mp","MTX")]
XS<-merge(x=XS, y=X_fix[,!colnames(X_fix)%in%c("Intercept")], by=c("id_new"), all.x=T)
XS=XS[with(XS, order(id_new,time_id)), ]
XS<-XS[,!colnames(XS)%in%c(poly_removecov)]
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




#---------------------------------------------------------------------------------------------------------#
# JM Cluster
#---------------------------------------------------------------------------------------------------------#
sink(paste(path_out,"JM_comp_Cl_Lymph.model.txt", sep="/"))
cat("
model {

for(i in 1:nId) 
{
 
 #Random Effect
 #-------------
  u_AP[i,1:4] ~ dmnorm(zeros[1:4],D_inv_AP[1:4,1:4]) #Global random intercept and slope of ANC , PLT
  
 for(cl in 1:nCluster)
  {
   u_L[i,1:2,cl] ~ dmnorm(zeros[5:6],D_inv_L[1:2,1:2,cl]) #Cluster wise random intercept and slope  of LYM
  }
  
 
 #Longitudinal Process
 #---------------------
for( j in offset[i]:(offset[i+1]-1))
 {
  #Mean ANC 
  mu_ANC[j]= u_AP[i,1] + u_AP[i,3]*X[j,4] + inprod(beta_ANC[1:(ncX-2)], X[j,3:ncX])
  y[j,4] ~ dnorm(mu_ANC[j],prec_ANC) #Missing imputation of ANC
  
  #Mean PLT
  mu_PLT[j]= u_AP[i,2] + u_AP[i,4]*X[j,4]  + inprod(beta_PLT[1:(ncX-2)], X[j,3:ncX])
  y[j,5] ~ dnorm(mu_PLT[j],prec_PLT) #Missing imputation of PLT
 
  #Log density of ANC, PLT
  ##ld_AP[j] <- logdensity.norm(y[j,4], mu_ANC[j], prec_ANC) + logdensity.norm(y[j,5], mu_PLT[j], prec_PLT)
  
  #Mean LYM
  for(cl in 1:nCluster)
  {
    mu_LYM[j,cl]= u_L[i,1,cl] + u_L[i,2,cl]*X[j,4] + inprod(beta_LYM[1:(ncX-2),cl], X[j,3:ncX]) 
  
  #Log density cluster wise of LYM
  ld_LYM[j,cl] <- logdensity.norm(y[j,6], mu_LYM[j,cl], prec_LYM[cl])  
  }
  
  } #loop of j
  
  #Complete longitudinal effect of ANC, PLT
  #------------------------------------------
   ##loglik_AP[i] = (logdensity.mnorm(u_AP[i,1:4], zeros[1:4], D_inv_AP[1:4,1:4]) #Random Effect of ANC , PLT
                 ##  +sum(ld_AP[offset[i]:(offset[i+1]-1)])   #Longitudinal of ANC, PLT
                  ## )
                    
   #Global Effects of ANC and PLT on Hazard and Survival
   
   #mu_ANC(t) and mu_PLT(t) in hazard (Survival part 1) 
   #----------------------------------------------------
   mu_TANC[i] =u_AP[i,1] + u_AP[i,3]*XT[i,4]  + inprod(beta_ANC[1:(ncXT-2)], XT[i,3:ncXT])
   mu_TPLT[i] =u_AP[i,2] + u_AP[i,4]*XT[i,4] + inprod(beta_PLT[1:(ncXT-2)], XT[i,3:ncXT])
   
   #mu_ANC(t) and mu_PLT(t) in Survival probability (Survival part 2) 
   #-----------------------------------------------------------------
   for(l in offset_s[i]:(offset_s[i+1]-1) )
  {
    mu_SANC[l] =u_AP[i,1] + u_AP[i,3]*XS[l,4] + inprod(beta_ANC[1:(ncXS-2)], XS[l,3:ncXS])
    mu_SPLT[l] =u_AP[i,2] + u_AP[i,4]*XS[l,4] + inprod(beta_PLT[1:(ncXS-2)], XS[l,3:ncXS])
  }
   
   for(cl in 1:nCluster)
  {
  
  #Survival part 1
  #---------------
 
   mu_TLYM[i,cl] =u_L[i,1,cl] + u_L[i,2,cl]*XT[i,4]  + inprod(beta_LYM[1:(ncXT-2),cl], XT[i,3:ncXT]) 
                               
   log.hT[i,cl]= (inprod(gam.h0[1:17,cl],BS_T[i,1:17]) #Baseline log hazard
                  + psi_LYM[cl]*mu_TLYM[i,cl]
                  + psi_ANC[cl]*mu_TANC[i] + psi_PLT[cl]*mu_TPLT[i]
                  + inprod(theta[1:(ncF-2),cl],X_fix[i,3:ncF])
                  ) 
 
 
 #Survival part 2
 #---------------
 
 for(l in offset_s[i]:(offset_s[i+1]-1) )
  {
   
    mu_SLYM[l,cl] =u_L[i,1,cl]  + u_L[i,2,cl]*XS[l,4] + inprod(beta_LYM[1:(ncXS-2),cl], XS[l,3:ncXS]) 
                                             
    SurvLong[l,cl] = wk_r[l]* ( exp( (inprod(gam.h0[1:17,cl],BS_S[l,1:17]) #Baseline log hazard
                                  + psi_LYM[cl]*mu_SLYM[l,cl] + psi_ANC[cl]*mu_SANC[l] + psi_PLT[cl]*mu_SPLT[l]) ) )
  }
   #Baseline covariates in Survival part 2
   #--------------------------------------
   surv_base[i,cl] = exp(inprod(theta[1:(ncF-2),cl],X_fix[i,3:ncF]))
   
   
    #The likelihood clusterwise
    lik[i,cl] <-  exp((logdensity.mnorm(u_L[i,1:2,cl], zeros[5:6], D_inv_L[1:2,1:2,cl]) #Random Effect of LYM
                       +sum(ld_LYM[offset[i]:(offset[i+1]-1),cl])   #Longitudinal LYM
                      ## +loglik_AP[i]   #log likelihood of ANC , PLT
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


#Global priors for ANC, PLT
#---------------------------
tau_AP[1:4,1]~dmnorm(zeros[1:4],Omega_AP)
tau_AP[1:4,2]~dmnorm(zeros[1:4],Omega_AP)
tau_AP[1:4,3]~dmnorm(zeros[1:4],Omega_AP)
tau_AP[1:4,4]~dmnorm(zeros[1:4],Omega_AP)
tau_AP[1:4,5]~dmnorm(zeros[1:4],Omega_AP)

D_inv_AP[1:4,1:4]<-tau_AP[1:4,1:5] %*% t(tau_AP[1:4,1:5]) + Omega_AP[1:4,1:4]*0.001
#D_AP[1:4,1:4] = inverse(D_inv_AP[1:4,1:4])

#Betas of ANC and PLT
#--------------------

for(h in 1:(ncX-2))
  {
  beta_ANC[h]~dnorm(0,0.001)
  beta_PLT[h]~dnorm(0,0.001)
  }


#Error variance of ANC and PLT
#-----------------------------
prec_ANC ~ dgamma(0.001,0.001)
sig_ANC = 1/prec_ANC
prec_PLT ~ dgamma(0.001,0.001)
sig_PLT = 1/prec_PLT


#pi with dirichlet prior
#-----------------------
pi[1:nCluster] ~ ddirch(alpha[1:nCluster])
 

  for(c in 1:nCluster)
  {  

  #Error Variance
  #--------------
  prec_LYM[c] ~ dgamma(0.001,0.001)
  sig_LYM[c] = 1/prec_LYM[c]
  

  #Variance for Random Effect of LYM 
  #---------------------------------
  tau_L[1:2,1,c]~dmnorm(zeros[1:2],(Omega[1:2,1:2]))
  tau_L[1:2,2,c]~dmnorm(zeros[1:2],(Omega[1:2,1:2]))
  tau_L[1:2,3,c]~dmnorm(zeros[1:2],(Omega[1:2,1:2]))

  
  D_inv_L[1:2,1:2,c]<-tau_L[1:2,1:3,c] %*% t(tau_L[1:2,1:3,c]) + Omega[1:2,1:2]*0.001
 # D_L[1:2,1:2,c]<-inverse(D_inv_L[1:2,1:2,c])
  
  
  #Beta
  #-----
  for(h in 1:(ncX-2))
  {
  beta_LYM[h,c]~dnorm(0,0.001)
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
                  'alpha'=rep(1,nCluster), 
                  'Omega_AP'=diag(rep(1,4)) ,'Omega'=diag(rep(1,2)) ,
                  'priorMean.gam.h0'=seq(-8.5,0,len =ncol(BS_S)) , 'M.gam.h0'=M.gam.h0 ,
                  'a.gam.h0'=1 , 'b.gam.h0'=0.05)


Params_save=c('u_L','u_AP', #Random effects
              'beta_LYM','beta_ANC','beta_PLT',  #The betas
              'sig_LYM', 'sig_ANC' ,'sig_PLT', #Error variances
              'D_inv_L','D_inv_AP', #Random effect covariance matrices inverse
              'psi_LYM','psi_ANC','psi_PLT','theta','gam.h0','tau.gam.h0', #Elements of Hazard
              'pi','pi_star','clus', #Cluster properties
              'lik' #Cluster wise likelihood
)


Sys.time()
cov_index<-NULL
#Set covs
if(deg==1)
{
  cov_index<-c(1,2,4:15)
}
if(deg==2)
{
  cov_index<-c(1,2,3,4:15)
}
if(deg==3)
{
  cov_index<-c(1,2,3,3,4:15)
}

#Initial values for actual model
init_jm=function(jm_init,nCluster,set_vec)
{
  #set_vec=c(1,3,2,1)
  vec<-set_vec[1:nCluster]
  
  v=list( u_L = jm_init$mean$u_L[,,vec] , u_AP = jm_init$mean$u_AP,
          beta_LYM = jm_init$mean$beta_LYM[cov_index,vec],  beta_ANC = jm_init$mean$beta_ANC[cov_index], beta_PLT = jm_init$mean$beta_PLT[cov_index],
          prec_LYM = 1/jm_init$mean$sig_LYM[vec], prec_ANC = 1/jm_init$mean$sig_ANC, prec_PLT = 1/jm_init$mean$sig_PLT,
          # D_inv_L = jm_init$mean$D_inv_L[,,vec] ,D_inv_AP = jm_init$mean$D_inv_AP ,
          psi_LYM = jm_init$mean$psi_LYM[vec] , psi_ANC = jm_init$mean$psi_ANC[vec], psi_PLT = jm_init$mean$psi_PLT[vec],
          theta =jm_init$mean$theta , gam.h0 = jm_init$mean$gam.h0, tau.gam.h0 =jm_init$mean$tau.gam.h0,
          pi= jm_init$mean$pi[vec])
  
  return(v)
}

#--------------------------------------------------Trial 2-----------------------------------------------------#

jm_init<-readRDS(paste(path_out,paste("/JM_comp_Lymph_Cl_",nCluster,".rds",sep=""),sep=""))
init_clus<-function(){init_jm(jm_init,nCluster,c(1,2,3,1))}

Sys.time()
#Run Cluster JM
jagsfit_JM <- jags(data=dat_list_JM, 
                   inits=init_clus, 
                   parameters.to.save=Params_save,
                   model.file=paste(path_out,'JM_comp_Cl_Lymph.model.txt', sep="/"),parallel = TRUE,
                   n.chains=2, n.iter=20000, n.burnin=5000,
                   n.thin=50,n.adapt = 10000,
                   DIC=TRUE)
saveRDS(jagsfit_JM,paste(path_out,paste("/JM_cluster_Result_",deg,"_",nCluster,".rds",sep=""),sep=""))

#-------------------------------------------------------------------------------------------------------------#






sink(paste(path_out,"JM_comp_Cl_Lymph1.model.txt", sep="/"))
cat("
model {

for(i in 1:nId) 
{
 
 #Cluster
 #-------
 G[i] ~ dcat(pi[1:nCluster])
 
 #Random Effect
 #-------------
  u_AP[i,1:4] ~ dmnorm(ze[1:4],D_inv_AP[1:4,1:4]) #Global random intercept and slope of ANC , PLT
  
  for(cl in 1:nCluster)
  {
   u_L[i,1:2,cl] ~ dmnorm(ze[5:6],D_inv_L[1:2,1:2,cl]) #Cluster wise random intercept and slope  of LYM
  }

 #Longitudinal Process
 #---------------------
for( j in offset[i]:(offset[i+1]-1))
 {
  #Mean ANC 
  mu_ANC[j]= u_AP[i,1] + u_AP[i,3]*X[j,4] + inprod(beta_ANC[1:(ncX-2)], X[j,3:ncX])
  y[j,4] ~ dnorm(mu_ANC[j],prec_ANC) #Missing imputation of ANC
  
  #Mean PLT
  mu_PLT[j]= u_AP[i,2] + u_AP[i,4]*X[j,4]  + inprod(beta_PLT[1:(ncX-2)], X[j,3:ncX])
  y[j,5] ~ dnorm(mu_PLT[j],prec_PLT) #Missing imputation of PLT

  #Mean LYM
  mu_LYM[j]= u_L[i,1,G[i]] + u_L[i,2,G[i]]*X[j,4] + inprod(beta_LYM[1:(ncX-2),G[i]], X[j,3:ncX]) 
  y[j,6] ~dnorm (mu_LYM[j], prec_LYM[G[i]])  

  } #loop of j
  

   #mu_ANC(t) and mu_PLT(t) in hazard (Survival part 1) 
   #----------------------------------------------------
   mu_TANC[i] =u_AP[i,1] + u_AP[i,3]*XT[i,4]  + inprod(beta_ANC[1:(ncXT-2)], XT[i,3:ncXT])
   mu_TPLT[i] =u_AP[i,2] + u_AP[i,4]*XT[i,4]  + inprod(beta_PLT[1:(ncXT-2)], XT[i,3:ncXT])
   
   #mu_LYM(t) in hazard (Survival part 1) 
   #-------------------------------------
   mu_TLYM[i] =u_L[i,1,G[i]] + u_L[i,2,G[i]]*XT[i,4]  + inprod(beta_LYM[1:(ncXT-2),G[i]], XT[i,3:ncXT]) 
   
   #Net hazard (Survival Part 1)
   #-----------------------------
   log.hT[i]= (inprod(gam.h0[1:17,G[i]],BS_T[i,1:17]) #Baseline log hazard
                  + psi_LYM[G[i]]*mu_TLYM[i]
                  + psi_ANC[G[i]]*mu_TANC[i] + psi_PLT[G[i]]*mu_TPLT[i]
                  + inprod(theta[1:(ncF-2),G[i]],X_fix[i,3:ncF])
                  ) 
   
   for(l in offset_s[i]:(offset_s[i+1]-1) )
  {
   #mu_ANC(t) and mu_PLT(t) in Survival probability (Survival part 2) 
   #-----------------------------------------------------------------
    mu_SANC[l] =u_AP[i,1] + u_AP[i,3]*XS[l,4] + inprod(beta_ANC[1:(ncXS-2)], XS[l,3:ncXS])
    mu_SPLT[l] =u_AP[i,2] + u_AP[i,4]*XS[l,4] + inprod(beta_PLT[1:(ncXS-2)], XS[l,3:ncXS])
    
    #mu_LYM(t) in Survival probability (Survival part 2) 
    #-----------------------------------------------------
    mu_SLYM[l] =u_L[i,1,G[i]]  + u_L[i,2,G[i]]*XS[l,4] + inprod(beta_LYM[1:(ncXS-2),G[i]], XS[l,3:ncXS]) 
    
    
    #Net Survival probability at l^th time point (Survival Part 2)
    #--------------------------------------------------------------
    SurvLong[l] = wk_r[l]* ( exp( (inprod(gam.h0[1:17,G[i]],BS_S[l,1:17]) #Baseline log hazard
                                  + psi_LYM[G[i]]*mu_SLYM[l] + psi_ANC[G[i]]*mu_SANC[l] + psi_PLT[G[i]]*mu_SPLT[l]) ) )

  }

   #Baseline covariates in Survival part 2
   #--------------------------------------
   surv_base[i] = exp(inprod(theta[1:(ncF-2),G[i]],X_fix[i,3:ncF]))

   #The rest of the log-likelihood forcluster G[i]
   #----------------------------------------------
   loglik_r[i] <-  ((y_s[i,2]*log.hT[i])   #Hazard
                    -( P[i]*surv_base[i]*sum(SurvLong[offset_s[i]:(offset_s[i+1]-1)]) ) ) #Survival function

  #Zeros trick in JAGS for JM
  #---------------------------------
  phi_poi[i] = Const-loglik_r[i]
  zeros[i] ~ dpoi(phi_poi[i])
  
}#loop of i


#Prior
#-------------------------------------------------------------------------

#pi with dirichlet prior 
#------------------------------------------(I)
pi[1:nCluster] ~ ddirch(alpha[1:nCluster])


#Global priors for ANC, PLT 
#-------------------------------------------(II)

D_inv_AP[1:4,1:4] ~ dwish(Omega_AP[1:4,1:4], 5)
D_AP[1:4,1:4] = inverse(D_inv_AP[1:4,1:4])

#Betas of ANC and PLT
#--------------------
for(h in 1:(ncX-2))
  {
  beta_ANC[h]~dnorm(0,0.001)
  beta_PLT[h]~dnorm(0,0.001)
  }


#Error variance of ANC and PLT
#-----------------------------
prec_ANC ~ dgamma(0.001,0.001)
sig_ANC = 1/prec_ANC
prec_PLT ~ dgamma(0.001,0.001)
sig_PLT = 1/prec_PLT


#Cluster-specific priors
#---------------------------------------(III)

  for(c in 1:nCluster)
  {  

  #Error Variance
  #--------------
  prec_LYM[c] ~ dgamma(0.001,0.001)
  sig_LYM[c] = 1/prec_LYM[c]
  

  #Variance for Random Effect of LYM 
  #---------------------------------
  D_inv_L[1:2,1:2,c] ~ dwish(Omega[1:2,1:2], 3) 
  D_L[1:2,1:2,c]<-inverse(D_inv_L[1:2,1:2,c])
  
  
  #Beta
  #-----
  for(h in 1:(ncX-2))
  {
  beta_LYM[h,c]~dnorm(0,0.001)
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
dat_list_JM =list('nId'=nId,'nCluster'=nCluster, 'zeros'=rep(1,nId), 'ze'=rep(0,6) ,'Const'=10^9, 
                  'y'=y,'y_s'=y_s, 'offset'=offset,'offset_s'=offset_s , 'BS_S'=BS_S , 'BS_T'=BS_T ,
                  'X'=X,'ncX'=ncol(X),'XS'=XS,'ncXS'=ncol(XS),'XT'=XT,'ncXT'=ncol(XT),'X_fix'=X_fix, 'ncF'=ncol(X_fix), 
                  'P'=P , 'wk_r'=wk_r , 
                  'alpha'=c(1,1,1), 
                  'Omega_AP'=diag(rep(1,4)) ,'Omega'=diag(rep(1,2)) ,
                  'priorMean.gam.h0'=seq(-8.5,0,len =ncol(BS_S)) , 'M.gam.h0'=M.gam.h0 ,
                  'a.gam.h0'=1 , 'b.gam.h0'=0.05)


Params_save=c('u_L','u_AP', #Random effects
              'beta_LYM','beta_ANC','beta_PLT',  #The betas
              'sig_LYM', 'sig_ANC' ,'sig_PLT', #Error variances
              'D_L','D_AP', #Random effect covariance matrices 
              'psi_LYM','psi_ANC','psi_PLT','theta','gam.h0','tau.gam.h0', #Elements of Hazard
              'pi','G' #Cluster properties
)

Sys.time()
#Run Lagged JM
jagsfit_JM1 <- jags(data=dat_list_JM, 
                    inits=NULL, 
                    parameters.to.save=Params_save,
                    model.file=paste(path_out,'JM_comp_Cl_Lymph1.model.txt', sep="/"),parallel = TRUE,
                    n.chains=1, n.iter=25, n.burnin=5,
                    n.thin=1,n.adapt = 5,
                    DIC=TRUE)


saveRDS(jagsfit_JM1,paste(path_out,paste("/_New_JM_comp_Lymph_Cl_",nCluster,".rds",sep=""),sep=""))





#-----------------------------------------------------------------------------------------------------------------#
#Model selection of polynomial degree
#-----------------------------------------------------------------------------------------------------------------#

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
library(emdbook)


#-------------------------------------------------------------------------------------------#
#File Path
fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"

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
X$Week3<-X$Week^3
X<-X[,c("id_new","time_id","Intercept","Week","Week2","Week3","X6_mp","MTX")]
X<-merge(x=X, y=X_fix[,!colnames(X_fix)%in%c("Intercept")], by=c("id_new"), all.x=T)
X=X[with(X, order(id_new,time_id)), ]
#X<-X[,!colnames(X)%in%c(poly_removecov)]
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
XT$Time_to_event3=XT$Time_to_event^3
XT<-XT[,c("id_new","time_id","Intercept","Time_to_event","Time_to_event2","Time_to_event3","X6_mp","MTX")]
XT<-merge(x=XT, y=X_fix[,!colnames(X_fix)%in%c("Intercept")], by=c("id_new"), all.x=T)
XT=XT[with(XT, order(id_new,time_id)), ]
#XT<-XT[,!colnames(XT)%in%c(poly_removecovT)]
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
XS$Week3=XS$Week^3
XS<-XS[,c("id_new","time_id","Intercept","Week","Week2","Week3","X6_mp","MTX")]
XS<-merge(x=XS, y=X_fix[,!colnames(X_fix)%in%c("Intercept")], by=c("id_new"), all.x=T)
XS=XS[with(XS, order(id_new,time_id)), ]
#XS<-XS[,!colnames(XS)%in%c(poly_removecov)]
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





#Degree of Polynomial
List_of_criterion_poly=function(nCluster, deg)
{
  
  
  poly_removecov=c("Week","Week2","Week3")[-c(1:deg)]
  poly_removecovT=c("Time_to_event","Time_to_event2","Time_to_event3")[-c(1:deg)]
  
  X<-X[,!colnames(X)%in%c(poly_removecov)]
  XT<-XT[,!colnames(XT)%in%c(poly_removecovT)]
  XS<-XS[,!colnames(XS)%in%c(poly_removecov)]
  
  ncX<-ncol(X)
  ncXS<-ncol(XS)
  ncXT<-ncol(XT)
  ncF<-ncol(X_fix)
  
  #Cluster JM file
  jfit<-NULL
  jfit<-readRDS(paste(path_out,paste("/JM_cluster_Result_",deg,"_",nCluster,".rds",sep=""),sep=""))
  #likelihood for a any particular theta
  likelihood=function(jfit1,nCluster)
  {
    #Means
    mu_ANC<-rep(NA, nrow(X))
    mu_PLT<-rep(NA,nrow(X))
    mu_LYM<-matrix(rep(NA,nrow(X)*nCluster),ncol=nCluster)
    
    #Random Effect
    ld_AP_ran<-rep(NA,nId)
    ld_LYM_ran<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    
    #longitudinal log likelihood
    ld_AP<-rep(NA, nrow(X))
    ld_LYM<-matrix(rep(NA,nrow(X)*nCluster),ncol=nCluster)
    
    #Mean in hazard
    mu_TANC<-rep(NA,nId)
    mu_TPLT<-rep(NA,nId)
    mu_TLYM<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    
    #The hazard
    log.hT<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    
    #Means in Survival part
    mu_SANC<-rep(NA,nrow(XS))
    mu_SPLT<-rep(NA,nrow(XS))
    mu_SLYM<-matrix(rep(NA,nrow(XS)*nCluster),ncol=nCluster)
    
    #Survival probability components without baseline covs
    SurvLong<-matrix(rep(NA,nrow(XS)*nCluster),ncol=nCluster)
    
    #Baseline covariates in survival
    surv_base<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    
    #Cluster wise likelihood
    lik<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    
    #Total likelihood
    total_lik<-rep(NA,nId)
    
    #There is only 1 NA in plt at index 4644
    j0<-4644
    i0<-as.integer(y[4644,1])
    y[j0,"Platelet"]<-as.numeric(jfit1$u_AP[i0,2] + jfit1$u_AP[i0,4]*X[j0,4] + as.numeric((jfit1$beta_PLT[1:(ncX-2)]) %*% (X[j0,3:ncX])))
    
    
    for(i in 1:nId)
    {  
      #Random Effect of ANC, PLT
      ld_AP_ran[i]<-dmvnorm(x=jfit1$u_AP[i,1:4], mu=c(0,0,0,0), Sigma=solve(jfit1$D_inv_AP[1:4,1:4]), log = FALSE, tol = 1e-07)
      
      #Random Effect of LYM
      for(cl in 1:nCluster)
      {
        ld_LYM_ran[i,cl]<-dmvnorm(x=jfit1$u_L[i,1:2,cl], mu=c(0,0), Sigma=solve(jfit1$D_inv_L[1:2,1:2,cl]), log = TRUE, tol = 1e-07)
      }
      
      #Longitudinal part
      for(j in offset[i]:(offset[i+1]-1))
      {
        #Mean ANC
        mu_ANC[j]<- as.numeric(jfit1$u_AP[i,1] + jfit1$u_AP[i,3]*X[j,4] + as.numeric((jfit1$beta_ANC[1:(ncX-2)]) %*% (X[j,3:ncX])))
        #Mean PLT
        mu_PLT[j]<- as.numeric(jfit1$u_AP[i,2] + jfit1$u_AP[i,4]*X[j,4] + as.numeric((jfit1$beta_PLT[1:(ncX-2)]) %*% (X[j,3:ncX])))
        
        #Likelihood of ANC, PLT (Longitudinal)
        ld_AP[j]<-(dnorm(x=y[j,4],mean=mu_ANC[j],sd=sqrt(jfit1$sig_ANC), log = TRUE)+
                     dnorm(x=y[j,5],mean=mu_PLT[j],sd=sqrt(jfit1$sig_PLT), log = TRUE))
        
        for(cl in 1:nCluster)
        {
          #Mean LYM
          mu_LYM[j,cl]<- as.numeric(jfit1$u_L[i,1,cl] + jfit1$u_L[i,2,cl]*X[j,4] + as.numeric((jfit1$beta_LYM[1:(ncX-2),cl]) %*% (X[j,3:ncX])))
          #Likelihood of LYM (Longitudinal)
          ld_LYM[j,cl]<-dnorm(x=y[j,6],mean=mu_LYM[j,cl],sd=sqrt(jfit1$sig_LYM[cl]), log = TRUE)
        }
        
      }
      
      
      #Survival part 1
      #---------------#
      
      #mu_ANC(T_i) and mu_PLT(T_i) in hazard 
      #-------------------------------------#
      mu_TANC[i] =as.numeric(jfit1$u_AP[i,1] + jfit1$u_AP[i,3]*XT[i,4] + as.numeric((jfit1$beta_ANC[1:(ncXT-2)]) %*% (XT[i,3:ncXT])))
      mu_TPLT[i] =as.numeric(jfit1$u_AP[i,2] + jfit1$u_AP[i,4]*XT[i,4] + as.numeric((jfit1$beta_PLT[1:(ncXT-2)]) %*% (XT[i,3:ncXT])))
      
      #mu_LYM(T_i,cl) in hazard 
      #------------------------#
      for( cl in 1:nCluster)
      {
        mu_TLYM[i,cl] =as.numeric(jfit1$u_L[i,1,cl] + jfit1$u_L[i,2,cl]*XT[i,4] + as.numeric((jfit1$beta_LYM[1:(ncXT-2),cl]) %*% (XT[i,3:ncXT])))
        
        log.hT[i,cl]= (as.numeric((jfit1$gam.h0[1:17,cl])%*%(BS_T[i,1:17])) #Baseline log hazard
                       + jfit1$psi_LYM[cl]*mu_TLYM[i,cl]
                       + jfit1$psi_ANC[cl]*mu_TANC[i] 
                       + jfit1$psi_PLT[cl]*mu_TPLT[i]
                       + as.numeric((jfit1$theta[1:(ncF-2),cl]) %*% t(X_fix[i,3:ncF]))
        )
      }
      
      
      #Survival part 2
      #---------------#
      for(cl in 1:nCluster)
      {
        #Baseline covariates in Survival part 2
        #--------------------------------------#
        surv_base[i,cl] = exp(as.numeric((jfit1$theta[1:(ncF-2),cl]) %*% t(X_fix[i,3:ncF])))
      }
      for(l in offset_s[i]:(offset_s[i+1]-1) )
      {
        #mu_ANC(t) and mu_PLT(t) in Survival probability 
        #---------------------------------------------------#
        mu_SANC[l] =as.numeric(jfit1$u_AP[i,1] + jfit1$u_AP[i,3]*XS[l,4] + as.numeric((jfit1$beta_ANC[1:(ncXS-2)]) %*% (XS[l,3:ncXS])))
        mu_SPLT[l] =as.numeric(jfit1$u_AP[i,2] + jfit1$u_AP[i,4]*XS[l,4] + as.numeric((jfit1$beta_PLT[1:(ncXS-2)]) %*% (XS[l,3:ncXS])))
        
        
        #mu_LYM(t,cl) in Survival probability 
        #------------------------------------#
        for( cl in 1:nCluster)
        {
          mu_SLYM[l,cl] =as.numeric(jfit1$u_L[i,1,cl] + jfit1$u_L[i,2,cl]*XS[l,4] + as.numeric((jfit1$beta_LYM[1:(ncXS-2),cl]) %*% (XS[l,3:ncXS])))
          
          SurvLong[l,cl] = wk_r[l]* ( exp( ( as.numeric((jfit1$gam.h0[1:17,cl])%*%(BS_S[l,1:17])) #Baseline log hazard
                                             + jfit1$psi_LYM[cl]*mu_SLYM[l,cl]
                                             + jfit1$psi_ANC[cl]*mu_SANC[l] 
                                             + jfit1$psi_PLT[cl]*mu_SPLT[l]
                                             + log(surv_base[i,cl]) ) ) )
          
        }
      }
      
      for( cl in 1:nCluster)
      {
        
        #The likelihood clusterwise
        lik[i,cl] <-  exp((ld_LYM_ran[i,cl]                              #Random Effect of LYM
                           +sum(ld_LYM[offset[i]:(offset[i+1]-1),cl])    #Longitudinal LYM
                           +sum(ld_AP[offset[i]:(offset[i+1]-1)])        #Longitudinal ANC , PLT
                           +ld_AP_ran[i]                                 #Random Effect of ANC, PLT 
                           +(y_s[i,2]*log.hT[i,cl])   #Hazard
                           -( P[i]*surv_base[i,cl]*sum(SurvLong[offset_s[i]:(offset_s[i+1]-1),cl]) ) #Survival function
        )) 
      }
      
      total_lik[i]<-sum(lik[i,1:nCluster]*jfit$sims.list$pi[1:nCluster])
      
    }
    
    return(total_lik)
    
  }
  
  #likelihood iteration wise
  likelihood_computation=function(jfit_samp,ite,nCluster)
  {
    #Means
    mu_ANC<-rep(NA, nrow(X))
    mu_PLT<-rep(NA,nrow(X))
    
    #Random Effect
    ld_AP_ran<-rep(NA,nId)
    ld_LYM_ran<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    
    #longitudinal log likelihood
    ld_AP<-rep(NA, nrow(X))
    
    lik<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    total_lik<-rep(NA,nId)
    
    #There is only 1 NA in plt at index 4644
    j0<-4644
    i0<-as.integer(y[4644,1])
    y[j0,"Platelet"]<-as.numeric(jfit_samp$u_AP[ite,i0,2] + jfit_samp$u_AP[ite,i0,4]*X[j0,4] + as.numeric((jfit_samp$beta_PLT[ite,1:(ncX-2)]) %*% (X[j0,3:ncX])))
    
    
    for(i in 1:nId)
    {  
      #Random Effect of ANC, PLT
      ld_AP_ran[i]<-dmvnorm(x=jfit_samp$u_AP[ite,i,1:4], mu=c(0,0,0,0), Sigma=solve(jfit_samp$D_inv_AP[ite,1:4,1:4]), log = FALSE, tol = 1e-07)
      
      #Longitudinal part
      for(j in offset[i]:(offset[i+1]-1))
      {
        #Mean ANC
        mu_ANC[j]<- as.numeric(jfit_samp$u_AP[ite,i,1] + jfit_samp$u_AP[ite,i,3]*X[j,4] + as.numeric((jfit_samp$beta_ANC[ite,1:(ncX-2)]) %*% (X[j,3:ncX])))
        #Mean PLT
        mu_PLT[j]<- as.numeric(jfit_samp$u_AP[ite,i,2] + jfit_samp$u_AP[ite,i,4]*X[j,4] + as.numeric((jfit_samp$beta_PLT[ite,1:(ncX-2)]) %*% (X[j,3:ncX])))
        
        #Likelihood of ANC, PLT (Longitudinal)
        ld_AP[j]<-(dnorm(x=y[j,4],mean=mu_ANC[j],sd=sqrt(jfit_samp$sig_ANC[ite]), log = TRUE)+
                     dnorm(x=y[j,5],mean=mu_PLT[j],sd=sqrt(jfit_samp$sig_PLT[ite]), log = TRUE))
        
      }
      
      for( cl in 1:nCluster)
      {
        lik[i,cl] <- jfit_samp$lik[ite,i,cl] * jfit$sims.list$pi[ite,cl]      
      }
      
      total_lik[i]<-sum(lik[i,1:nCluster]) *exp((sum(ld_AP[offset[i]:(offset[i+1]-1)]) +ld_AP_ran[i]))
      
    }
    return(total_lik)
  }
  itera<-1:jfit$mcmc.info$n.samples
  
  llc<-matrix(rep(NA,nId*length(itera)),ncol=nId)
  d_theta<-rep(NA,length(itera))
  likli_inv<-matrix(rep(NA,nId*length(itera)),ncol=nId)
  likli<-matrix(rep(NA,nId*length(itera)),ncol=nId)
  
  for(g in 1:length(itera))
  {
    llc[g,]<-likelihood_computation(jfit_samp=jfit$sims.list,ite=itera[g],nCluster)
    d_theta[g]<- -2*sum(log(llc[g,]))
    likli_inv[g,]<- 1/llc[g,]
    likli[g,]<- llc[g,]
  }
  
  d_theta_hat<- -2*sum(log(likelihood(jfit1=jfit$mean,nCluster)))
  
  #log-likelihood
  log_liki<-sum(log(likelihood(jfit1=jfit$mean,nCluster)))
  #DIC
  p_D<-mean(d_theta)- d_theta_hat
  dic_comp<- mean(d_theta) + p_D
  
  paramater_vec<-c('u_L','u_AP', #Random effects
                   'beta_LYM','beta_ANC','beta_PLT',  #The betas
                   'sig_LYM', 'sig_ANC' ,'sig_PLT', #Error variances
                   'D_inv_L','D_inv_AP', #Random effect covariance matrices inverse
                   'psi_LYM','psi_ANC','psi_PLT','theta','gam.h0','tau.gam.h0' #Elements of Hazard
  )
  num_parameters<-sum(sapply(1:length(paramater_vec), function(i) length(jfit$mean[[paramater_vec[i]]])))
  
  #AIC
  aic_comp<- -2*sum(log(likelihood(jfit1=jfit$mean,nCluster))) + 2*num_parameters
  
  #BIC
  bic_comp<- -2*sum(log(likelihood(jfit1=jfit$mean,nCluster))) + num_parameters*log(nId)
  
  #LPML
  
  cpo<-sapply(1:nId,function(i) 1/mean(likli_inv[1:length(itera),i]))
  
  lpml_comp<- sum(log(cpo))
  
  #WAIC
  
  
  #WAIC_1
  p_waic1<- 2* ( sum(sapply(1:nId,function(i) log(mean(likli[1:length(itera),i])) ))
                 -mean(sapply(1:length(itera), function(g) sum(log(likli[g,1:nId])) ))
  )
  
  waic_1_comp<-mean(d_theta) + p_waic1
  
  #WAIC_2
  p_waic2<- sum(sapply(1:nId,function(i) var(log(likli[1:length(itera),i])) ))
  
  waic_2_comp<-mean(d_theta) + p_waic2
  
  #Minimum cluster percentage
  min_clust_comp<-round(min(table(sapply(1:nId,function(i)which.max(jfit$mean$pi_star[i,])))*100/184),3)
  #Criterion vector
  Criterion_vector=c(nCluster,deg,log_liki,dic_comp,aic_comp,bic_comp,lpml_comp,waic_1_comp,waic_2_comp,min_clust_comp)
  names(Criterion_vector)=c("Clucters_num","deg","log_likelihood","DIC","AIC","BIC","LPML","WAIC_1","WAIC_2","Minimum_Cluster_size")
  return(Criterion_vector)
  
}
jfit<-NULL
List_of_criterion_poly(nCluster=2, deg=1)
jfit<-NULL
List_of_criterion_poly(nCluster=2, deg=3)



#-----------------------------------------------------------------------------------------------------------------#

















#Diagnostics



#Total Clusters
nCluster=2

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
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)


#-------------------------------------------------------------------------------------------#
#File Path
fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"

#Input path
path=paste(fp,"/Cluster_JM/Exp_work/DATA",sep="")
#Output path
path_out=paste(fp,"/JAGS_Output/Cluster_JM",sep="")


get_all_jags=function(chain)
{
  jags_list=vector(mode="list", length=1)
  for(x in 1:1)
  {
    l=NULL
    l=readRDS(paste(path_out,paste("/JM_comp_Lymph_Cl_",nCluster,".rds",sep=""),sep=""))
    jags_list[[x]]<- l$samples[[chain]]
    x=x+1
  }
  return(jags_list)
}

#To view posterior density Chain wise
pd=function(chain, param)
{
  #Pooled sample
  pool<-get_all_jags(chain)
  #Extracting from the list named pool
  list_extract=function(v)
  {
    vec=c(as.vector(pool[[1]][,v]))
    return(vec)
  }
  vec=list_extract(v=param)
  #vec=vec[c(seq(1,20000,50))]
  df=data.frame(Iterations=1:length(vec), y=vec,y_cumu=cummean(vec))
  p1<-ggplot(data=df,aes(x=y))+
    geom_density()+theme_bw()+
    xlab("")+ylab(param)+
    geom_vline(xintercept = mean(vec), col="red")+
    theme(plot.title = element_text(hjust = 0.5))+ggtitle("Density Plot")
  p2<-ggplot(data=df,aes(x=Iterations,y=y))+
    geom_line()+theme_bw()+
    xlab("Iterations")+ylab(param)+
    geom_hline(yintercept = mean(vec), col="red")+
    theme(plot.title = element_text(hjust = 0.5))+ggtitle("Trace Plot")
  
  p3<-ggplot(data=df,aes(x=Iterations,y=y_cumu))+
    geom_line()+theme_bw()+
    xlab("Iterations")+ylab(param)+
    geom_hline(yintercept = mean(vec), col="red")+
    theme(plot.title = element_text(hjust = 0.5))+ggtitle("Cumulative Trace")
  
  bacf <- acf(df$y, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  conf.level <- 0.95
  ciline <- qnorm((1 - conf.level)/2)/sqrt(length(df$y))
  
  p4 <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = c(0)), col="black") +
    geom_hline(aes(yintercept = ciline), color = 'blue') +
    geom_hline(aes(yintercept = -ciline), color = 'blue')+
    geom_segment(mapping = aes(xend = lag, yend = 0))+theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+ggtitle("ACF")
  p5<-grid.arrange(p1,p2,p3,p4, ncol=2, nrow=2)
  return(p5)
}

#To view grand posterior density
pd_ALL=function( param,rows, cols)
{
  l=readRDS(paste(path_out,paste("/JM_comp_Lymph_Cl_",nCluster,".rds",sep=""),sep=""))
  vec1=l$sims.list[[param]]
  vec<-NULL
  param1<-NULL
  if(length(dim(vec1))>=3)
  {
    vec<-vec1[,rows,cols] 
    param1<-paste(param,"[",rows,",",cols,"]",sep="")
  }
  if(length(dim(vec1))<3)
  {
    vec<-vec1[,rows]
    param1<-paste(param,"[",rows,"]",sep="")
  }
  
  #vec=vec[c(seq(1,20000,50))]
  df=data.frame(Iterations=1:length(vec), y=vec,y_cumu=cummean(vec))
  p1<-ggplot(data=df,aes(x=y))+
    geom_density()+theme_bw()+
    xlab("")+ylab(param1)+
    geom_vline(xintercept = mean(vec), col="red")+
    theme(plot.title = element_text(hjust = 0.5))+ggtitle("Density Plot")
  p2<-ggplot(data=df,aes(x=Iterations,y=y))+
    geom_line()+theme_bw()+
    xlab("Iterations")+ylab(param1)+
    geom_hline(yintercept = mean(vec), col="red")+
    theme(plot.title = element_text(hjust = 0.5))+ggtitle("Trace Plot")
  
  p3<-ggplot(data=df,aes(x=Iterations,y=y_cumu))+
    geom_line()+theme_bw()+
    xlab("Iterations")+ylab(param1)+
    geom_hline(yintercept = mean(vec), col="red")+
    theme(plot.title = element_text(hjust = 0.5))+ggtitle("Cumulative Trace")
  
  bacf <- acf(df$y, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  conf.level <- 0.95
  ciline <- qnorm((1 - conf.level)/2)/sqrt(length(df$y))
  
  p4 <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = c(0)), col="black") +
    geom_hline(aes(yintercept = ciline), color = 'blue') +
    geom_hline(aes(yintercept = -ciline), color = 'blue')+
    geom_segment(mapping = aes(xend = lag, yend = 0))+theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+ggtitle("ACF")
  p5<-grid.arrange(p1,p2,p3,p4, ncol=2, nrow=2)
  return(p5)
}

#pd(chain=1,param="psi_ANC[2]")
pd_ALL(param="psi_LYM",rows=1,cols=NA)
pd_ALL(param="psi_ANC",rows=1,cols=NA)
pd_ALL(param="psi_PLT",rows=1,cols=NA)
pd_ALL(param="psi_LYM",rows=2,cols=NA)
pd_ALL(param="psi_ANC",rows=2,cols=NA)
pd_ALL(param="psi_PLT",rows=2,cols=NA)

pe=function(chain, statistic)
{
  
  #Pooled sample
  pool<-get_all_jags(chain)
  
  #Reference of size of parameters
  jagsfits<- readRDS(paste(path_out,paste("/JM_comp_Lymph_Cl_",nCluster,".rds",sep=""),sep=""))
  
  #Compute summary statistics
  summary_stat=function(vec, sum_stat)
  {
    vec1=vec#[c(seq(1,20000,50))]
    z=NULL
    if(sum_stat=="mean")
    {
      z=round(mean(vec1),3)
    }
    if(sum_stat=="median")
    {
      z=round(quantile(vec1, probs=c(0.5), names=F),3)
    }
    if(sum_stat=="sd")
    {
      z=round(sd(vec1),3)
    }
    if(sum_stat=="mode")
    {
      getmode <- function(v) 
      {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
      }
      
      z=getmode(vec1)
    }
    if(sum_stat=="significance")
    {
      z=ifelse(prod(round(quantile(vec1, probs=c(0.025,0.975),names=F),3))<=0,"No","Yes") 
    }
    if(sum_stat=="CI")
    {
      z1=round(quantile(vec1, probs=c(0.025,0.975),names=F),3)
      z=paste0("(",z1[1]," , ",z1[2],")",collapse="")
    }
    return(z)
  }
  
  #Extract summary values for parameters from matrix format
  Get_vals=function(pool,param)
  {
    arr=NULL
    dim_p=NULL
    dim_p=dim(jagsfits$mean[[param]])
    
    #Extracting from the list named pool
    list_extract=function(v)
    {
      vec=c(as.vector(pool[[1]][,v]))
      return(vec)
    }
    if(length(dim_p)<=0) # unit dimensional array
    {
      arr<-NA
      arr<-summary_stat(vec=list_extract(v=param),sum_stat=statistic)
    }
    
    if(length(dim_p)==1) #1d array
    {
      arr<-array(NA, dim=c(dim_p))
      for(i in 1:dim_p)
      {
        arr[i]<- summary_stat(vec=list_extract(v=paste0(param,"[",i,"]",collapse="")),sum_stat=statistic) 
        i=i+1
      }
    }
    
    if(length(dim_p)==2) #2d array
    {
      arr<-array(NA, dim=c(dim_p))
      for(i in 1:dim_p[1])
      {
        for(j in 1:dim_p[2])
        {
          arr[i,j]<- summary_stat(vec=list_extract(v=paste0(param,"[",i,",",j,"]",collapse="")),sum_stat=statistic) 
          j=j+1
        }
        i=i+1
      }
    }
    
    if(length(dim_p)>=3) #for 3d array
    {
      arr<-array(NA, dim = c(dim_p))
      for(i in 1:dim_p[1])
      {
        for(j in 1:dim_p[2])
        {
          for(k in 1:dim_p[3])
          {
            arr[i,j,k]<-  summary_stat(vec=list_extract(v=paste0(param,"[",i,",",j,",",k,"]",collapse="")),sum_stat=statistic)
            k=k+1
          }
          j=j+1
        }
        i=i+1
      }
      
      
    }
    return(arr)
  }
  params_interest<-NULL
  if(statistic!="mode")
  {
    params_interest<-jagsfits$parameters[!jagsfits$parameters%in%c("q_p","deviance","tau","D_inv","clus")]
    
  }
  if(statistic=="mode")
  {
    params_interest<-c("clus")
  }
  liz<-vector(mode="list",length=length(params_interest))
  names(liz)<-params_interest
  i=NA
  for(i in 1:length(params_interest))
  {
    liz[[params_interest[i]]]<- Get_vals(pool,params_interest[i])
    i=i+1
  }
  
  return(liz)
  
  
}

#Summary statistic
#95%CI
ci=pe(chain=2, statistic = "CI")
ci$psi_LYM
ci$psi_ANC
ci$psi_PLT

#Significance
signif=pe(chain=1, statistic ="significance" )
signif$psi_LYM
signif$psi_ANC
signif$psi_PLT

#Mean
pos_mean=pe(chain=2, statistic ="mean" )
pos_mean$pi
round(pos_mean$pi_star,2)

#Median
pos_median=pe(chain=1, statistic ="median" )

#Mode for discrete
pos_mode=pe(chain=1, statistic="mode")
table(pos_mode$clus)
#df=data.frame("id"=1:236,"Cluster"=pos_mode$clus,"subcluster"=pos_mode$subclus ) #subcluster is for DP
#round(table(df[df$Cluster==1,"subcluster"])*100/nrow(df[df$Cluster==1,]))

#sd
pos_sd=pe(chain=1, statistic ="sd" )

