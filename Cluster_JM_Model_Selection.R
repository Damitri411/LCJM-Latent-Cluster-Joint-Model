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

ncX<-ncol(X)
ncXS<-ncol(XS)
ncXT<-ncol(XT)
ncF<-ncol(X_fix)



#jfit1<-jfit$mean






#Number of cluster
List_of_criterion=function(nCluster)
{
  #Cluster JM file
  jfit<-NULL
  jfit<-readRDS(paste(path_out,paste("/JM_comp_Lymph_Cl_",nCluster,".rds",sep=""),sep=""))
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
  Criterion_vector=c(nCluster,log_liki,dic_comp,aic_comp,bic_comp,lpml_comp,waic_1_comp,waic_2_comp,min_clust_comp)
  names(Criterion_vector)=c("Clucters_num","log_likelihood","DIC","AIC","BIC","LPML","WAIC_1","WAIC_2","Minimum_Cluster_size")
  return(Criterion_vector)
  
}
jfit<-NULL
List_of_criterion(nCluster=2)
jfit<-NULL
List_of_criterion(nCluster=3)
