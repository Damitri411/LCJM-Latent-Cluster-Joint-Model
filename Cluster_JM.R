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

#Removing NAs in just ANC and PLT
data_org$ANC_PLT_na=data_org$ANC+data_org$Platelet
data_org=data_org[!is.na(data_org$ANC_PLT_na),]
data_org=data_org[,!colnames(data_org)%in%c("ANC_PLT_na")]


data_org$id_new=as.numeric(data_org$id_new)
data_org$time_id=as.numeric(data_org$time_id)
data_org=data_org[
  with(data_org, order(id_new, time_id)),
]


#Making every id to the same starting point.
min_w=data_org%>%
  group_by(Child.Id)%>%
  summarise(min_week=min(Week)*277)
data_org=merge(x=data_org, y=min_w,by=c("Child.Id"),all.x=T)

data_org$Week=data_org$Week*277-data_org$min_week+1
data_org$Last_follow_up_date=data_org$Last_follow_up_date*277-data_org$min_week+1
data_org$Time_to_event=data_org$Time_to_event*277-data_org$min_week+1

#Scaling the time (max week is now 271 weeks instead of previous 277)
max_week=max(max(data_org$Week),max(data_org$Last_follow_up_date),
             max(data_org$Time_to_event))
data_org$Week=data_org$Week/max_week
data_org$Time_to_event=data_org$Time_to_event/max_week
data_org$Last_follow_up_date=data_org$Last_follow_up_date/max_week

#Removing columns
data_org=data_org[,!colnames(data_org)%in%c("Height","Weight","BSA","min_week")]



#Removing id's with < 5 times points and updating the id names
id_set=unique(data_org$Child.Id)
dat=NULL
for(i in 1:length(id_set))
{
  d=subset(data_org, data_org$Child.Id==id_set[i])
  dat=rbind(dat, c(id_set[i],nrow(d)))
  i=i+1
}
dat=data.frame(dat)
rownames(dat)=1:nrow(dat)
dat$X2=as.numeric(dat$X2)
dat1=subset(dat,dat$X2<5)
data_org1=subset(data_org, !data_org$Child.Id%in%unique(dat1$X1))
unique(data_org1$id_new)

dat1=data_org1[,c("Child.Id","id_new")]
dat1=dat1[!duplicated(dat1), ]
dat1$id_new=sapply(dat1$Child.Id,function(x)as.numeric(strsplit(x,"_")[[1]][2]))
dat1=dat1[with(dat1, order(id_new)), ]
dat1$id_new=1:nrow(dat1)
data_org1=merge(x=data_org1[,!colnames(data_org1)%in%c("id_new")], y=dat1, by=c("Child.Id"), all.x=T)
data_org1=data_org1[,c(1,26,2:25)]
data_org=data_org1
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

data_org$WBC_new=exp(data_org$WBC)-exp(data_org$ANC)
data_org$WBC_new=log(data_org$WBC_new)
data_org$WBC=0
data_org$WBC=data_org$WBC_new

data_org=data_org[,!colnames(data_org)%in%c("WBC_new")]
data_org$Platelet=data_org$Platelet+log(192000)

#Centering to 0
#---------------#
data_org$WBC=data_org$WBC-median(data_org$WBC[!is.na(data_org$WBC)])
data_org$ANC=data_org$ANC-median(data_org$ANC[!is.na(data_org$ANC)])
data_org$Platelet=data_org$Platelet-median(data_org$Platelet[!is.na(data_org$Platelet)])

#Exclude variables Day 8 , Day 35 risk stratification and Status
#---------------#
data_org=data_org[,!colnames(data_org)%in%c("Day35","Status")]

#Getting intercept
data_org$Intercept=1

saveRDS(data_org,paste(path,"/Transformed_DATA.rds",sep=""))


#----------------------------------------------------------------------------------------------#

#number of ids 
#-----------------------------------------------------------------------------------------------------
nId <- length(unique(data_org$id_new))

#Response
#-----------------------------------------------------------------------------------------------------#
y <- as.matrix(data_org[,c("id_new","time_id","WBC","ANC","Platelet")])

#------------------------------------------------------------------------------------------------------#

#Fixed covariates
X_fix=data_org[,c("id_new","Intercept","Age","Gender","Lineage","WBC_high","NCI","Bulky","CNS","Risk","Risk_present","Day8","Morph_CR1")]
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



#-------------------------------------------------------------------------------------------#
#                                 Cluster JM Model 
#-------------------------------------------------------------------------------------------#

sink(paste(path_out,"Cluster_JM.model.txt", sep="/"))
cat("
model {

for(i in 1:nId) 
{
 #Cluster
  clus[i] ~ dcat(pi[i,1:nCluster])
  
 #Random Effect
 #----------
 for(cl in 1:nCluster)
  {
   u_L[i,cl] ~ dnorm(0,D_inv[cl]) #Cluster wise univariate random effect of Lymphocyte
  }
 u_AP[i,1:2] ~ dmnorm(ze[1:2],D_inv_AP[1:2,1:2]) #Random Effect of ANC and Platelet

 #Longitudinal Process
 #---------------------
for( j in offset[i]:(offset[i+1]-1))
 {
   
   mu_ANC[j] =u_AP[i,1] + inprod(beta_ANC[1:(ncX-2)], X[j,3:ncX]) #ANC mean
   mu_PLT[j] =u_AP[i,2] + inprod(beta_PLT[1:(ncX-2)], X[j,3:ncX]) #PLT mean
   
   y[j,4] ~ dnorm(mu_ANC[j], prec_ANC) #Model for ANC
   y[j,5] ~ dnorm(mu_PLT[j], prec_PLT) #Model for PLT
 
  #Mean
  mu_LYM[j]=(u_L[i,clus[i]]  + inprod(beta_LYM[1:(ncX-2),clus[i]], X[j,3:ncX]) 
                             + beta_LYM[(ncX-1),clus[i]]*mu_ANC[j]               #Fixed part of ANC
                             + beta_LYM[ncX,clus[i]]*mu_PLT[j] )                 #Fixed part of PLT
  
  y[j,3] ~ dnorm(mu_LYM[j], prec_LYM[clus[i]]) #Model for LYMPH
  
  } #loop of j
  
  #log probability (pi[i,c])
  #-------------------------
  log.pi[i]= log(pi[i,clus[i]])
  
  #Survival part 1
  #---------------
 
   mu_TANC[i] =u_AP[i,1] + inprod(beta_ANC[1:(ncXT-2)], XT[i,3:ncXT]) #ANC mean
   mu_TPLT[i] =u_AP[i,2] + inprod(beta_PLT[1:(ncXT-2)], XT[i,3:ncXT]) #PLT mean
   
   mu_TLYM[i] =(u_L[i,clus[i]] + inprod(beta_LYM[1:(ncXT-2),clus[i]], XT[i,3:ncXT]) 
                               + beta_LYM[(ncXT-1),clus[i]]*mu_TANC[i]            #Fixed part of ANC
                               + beta_LYM[ncXT,clus[i]]*mu_TPLT[i] )              #Fixed part of PLT
   
   log.hT[i]= (inprod(gam.h0[1:17,clus[i]],BS_T[i,1:17]) #Baseline log hazard
              + psi[clus[i]]*mu_TLYM[i]
              + inprod(theta[1:(ncF-2),clus[i]],X_fix[i,3:ncF]) ) #Baseline covariates
 
 
 #Survival part 2
 #---------------
 
 for(l in offset_s[i]:(offset_s[i+1]-1) )
  {
    mu_SANC[l] =u_AP[i,1] + inprod(beta_ANC[1:(ncXS-2)], XS[i,3:ncXS]) #ANC mean
    mu_SPLT[l] =u_AP[i,2] + inprod(beta_PLT[1:(ncXS-2)], XS[i,3:ncXS]) #PLT mean
    
    mu_SLYM[l] =(u[i,clus[i]]  + inprod(beta_LYM[1:(ncXS-2),clus[i]], XS[i,3:ncXS]) 
                               + beta[(ncXS-1),clus[i]]*mu_SANC[l]                #Fixed part of ANC
                               + beta[ncXS,clus[i]]*mu_SPLT[l] )                  #Fixed part of PLT
                                                        
    SurvLong[l] = wk_r[l]* ( exp( inprod(gam.h0[1:17,clus[i]],BS_S[l,1:17])+ #Baseline log hazard
                                  + psi[clus[i]]*mu_SLYM[l] ) )
  }
 
  #Zeros trick in JAGS for JM
  #---------------------------------
  phi_poi[i] = Const- ( log.pi[i] + (y_s[i,2]*log.hT[i]) -  ( P[i]*exp(inprod(theta[1:(ncF-2),clus[i]],X_fix[i,3:ncF]))*sum(SurvLong[offset_s[i]:(offset_s[i+1]-1)]) ) )
  zeros[i] ~ dpois(phi_poi[i])

}#loop of i


  
#Prior
#-------------------------------------

#Global parameters
#-----------------
#Betas
#-----
for( k in 1:(ncX-2))
{
 beta_ANC[k]~dnorm(0,0.001)
 beta_PLT[k]~dnorm(0,0.001)
}

#Precision of ANC and PLT
#------------------------
  prec_ANC ~ dgamma(0.001,0.001)
  sig_ANC = 1/prec_ANC
  prec_PLT ~ dgamma(0.001,0.001)
  sig_PLT =1/prec_PLT
  
#Random Effect var-cov matrix
#----------------------------
 D_inv_AP[1:2,1:2] ~ dwish(Omega[,], 3)
 D_AP[1:2,1:2] <-inverse(D_inv_AP[1:2,1:2])


#Construction of pi[i,c]
#-----------------------
  for( i in 1:nId)
  {
  
  for(c in 1:(nCluster-1))
     {
       pi_st[i,c]=exp((inprod(phi[c,1:(ncF-1)],X_fix[i,2:ncF])))
     }
    pi_st[i,nCluster]=1

    pi[i,1:nCluster]=pi_st[i,1:nCluster]/(sum(pi_st[i,1:nCluster]))
  }
 
 #The cluster parameters
 #-------------------------
 for(c in 1:(nCluster-1)) 
  {
    for(l in 1:(ncF-1))
    {
    phi[c,l] ~dnorm(0,0.001) #Coefficients of intercept and covariates
    }
  }

  for(c in 1:nCluster)
  {  
  
  #Error Variance
  #--------------
  prec_LYM[c] ~ dgamma(0.001,0.001)
  sig_LYM[c] = 1/prec_LYM[c]
  

  #Variance for Random Effect
  #--------------------------
  D_inv[c] ~ dgamma(0.001,0.001)
  D_L[c] <- 1/D_inv[c]
  
  
  #Beta
  #-----
  for(h in 1:ncX)
  {
  beta_LYM[h,c]~dnorm(0,0.001)
  }
  
  #Association parameter
  #----------------------
  psi[c] ~ dnorm(0,0.001)
  
  
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
  }
}
",fill = TRUE)
sink()



#Data list for Lagged JM
dat_list_Cluster_JM =list('nId'=nId,'nCluster'=nCluster, 'zeros'=zeros ,'Const'=10^7, 
                          'y'=y,'y_s'=y_s, 'offset'=offset,'offset_s'=offset_s , 'BS_S'=BS_S , 'BS_T'=BS_T ,
                          'X'=X,'ncX'=ncol(X),'XS'=XS,'ncXS'=ncol(XS),'XT'=XT,'ncXT'=ncol(XT),'X_fix'=X_fix, 'ncF'=ncol(X_fix), 
                          'P'=P , 'wk_r'=wk_r ,  
                          'ze'=c(0,0),'Omega' = diag(rep(1,2)) ,'priorMean.gam.h0'=seq(-8.5,0,len =ncol(BS_S)) , 'M.gam.h0'=M.gam.h0 ,
                          'a.gam.h0'=1 , 'b.gam.h0'=0.05)


Params_save=c('u_L','u_AP',
              'beta_LYM','beta_ANC','beta_PLT', 
              'sig_LYM','sig_ANC','sig_PLT' , 
              'D_L','D_AP', # var-cov
              'psi','theta','gam.h0', #Hazard
              'pi','clus','phi' #Cluster
)

Sys.time()
#Run Lagged JM
jagsfit_Cluster_JM <- jags(data=dat_list_Cluster_JM, 
                                  inits=NULL, 
                                  parameters.to.save=Params_save,
                                  model.file=paste(path_out,'Cluster_JM.model.txt', sep="/"),parallel = TRUE,
                                  n.chains=1, n.iter=12000, n.burnin=2000,
                                  n.thin=10,n.adapt = 2000,
                                  DIC=TRUE)


saveRDS(jagsfit_Cluster_JM,paste(path_out,"/Cluster_JM.rds",sep=""))



#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
#                                 Cluster JM EM Model 
#-------------------------------------------------------------------------------------------#
sink(paste(path_out,"Cluster_JM_EM.model.txt", sep="/"))
cat("
model {

for(i in 1:nId) 
{
 #Cluster
  clus[i] ~ dcat(pi[i,1:nCluster])
  
 #Random Effect
 #----------
 for(cl in 1:nCluster)
  {
   u_L[i,cl] ~ dnorm(0,D_inv[cl]) #Cluster wise univariate random effect of Lymphocyte
  }

 #Longitudinal Process
 #---------------------
for( j in offset[i]:(offset[i+1]-1))
 {
   
  #Mean
  for(cl in 1:nCluster)
  {
    mu_LYM[j,cl]=u_L[i,cl]  + inprod(beta_LYM[1:(ncX-2),cl], X[j,3:ncX]) 
  }
  y[j,3] ~ dnorm(mu_LYM[j,clus[i]], prec_LYM[clus[i]]) #Model for LYMPH
  
  #Log density cluster wise
  for(cl in 1:nCluster)
  {
  ld_LYM[j,cl] <- logdensity.norm(y[j,3], mu_LYM[j,cl], prec_LYM[cl])  
  }
  
  } #loop of j
   
   
   for(cl in 1:nCluster)
  {
  
  #Survival part 1
  #---------------
 
   mu_TLYM[i,cl] =u_L[i,cl] + inprod(beta_LYM[1:(ncXT-2),cl], XT[i,3:ncXT]) 
                               
   log.hT[i,cl]= (inprod(gam.h0[1:17,cl],BS_T[i,1:17]) #Baseline log hazard
                  + psi[cl]*mu_TLYM[i,cl]
                  + inprod(theta[1:(ncF-2),cl],X_fix[i,3:ncF]) ) #Baseline covariates
 
 
 #Survival part 2
 #---------------
 
 for(l in offset_s[i]:(offset_s[i+1]-1) )
  {
   
    mu_SLYM[l,cl] =u_L[i,cl]  + inprod(beta_LYM[1:(ncXS-2),cl], XS[i,3:ncXS]) 
                                             
    SurvLong[l,cl] = wk_r[l]* ( exp( (inprod(gam.h0[1:17,cl],BS_S[l,1:17]) #Baseline log hazard
                                  + psi[cl]*mu_SLYM[l,cl]) ) )
  }
   #Baseline covariates in Survival part 2
   #--------------------------------------
   surv_base[i,cl] = exp(inprod(theta[1:(ncF-2),cl],X_fix[i,3:ncF]))
   
    #The likelihood clusterwise
    lik[i,cl] <-  exp((logdensity.norm(u_L[i,cl], 0, D_inv[cl])     #Random Effect
                       +sum(ld_LYM[offset[i]:(offset[i+1]-1),cl])   #Longitudinal Lymphocyte
                       +(y_s[i,2]*log.hT[i,cl])   #Hazard
                       -( P[i]*surv_base[i,cl]*sum(SurvLong[offset_s[i]:(offset_s[i+1]-1),cl]) ) )) #Survival function
                     

  }
 
  #Zeros trick in JAGS for JM
  #---------------------------------
  phi_bern[i] = inprod(lik[i,1:nCluster],pi[i,1:nCluster])/Const
  ones[i] ~ dbern(phi_bern[i])

}#loop of i


  
#Prior
#-------------------------------------

#Construction of pi[i,c]
#-----------------------
  for( i in 1:nId)
  {
  
  for(c in 1:(nCluster-1))
     {
       pi_st[i,c]=exp((inprod(phi[c,1:(ncF-1)],X_fix[i,2:ncF])))
     }
    pi_st[i,nCluster]=1

    pi[i,1:nCluster]=pi_st[i,1:nCluster]/(sum(pi_st[i,1:nCluster]))
  }
 
 #The cluster parameters
 #-------------------------
 for(c in 1:(nCluster-1)) 
  {
    for(l in 1:(ncF-1))
    {
    phi[c,l] ~dnorm(0,0.001) #Coefficients of intercept and covariates
    }
  }

  for(c in 1:nCluster)
  {  
  
  #Error Variance
  #--------------
  prec_LYM[c] ~ dgamma(0.001,0.001)
  sig_LYM[c] = 1/prec_LYM[c]
  

  #Variance for Random Effect
  #--------------------------
  D_inv[c] ~ dgamma(0.001,0.001)
  D_L[c] <- 1/D_inv[c]
  
  
  #Beta
  #-----
  for(h in 1:(ncX-2))
  {
  beta_LYM[h,c]~dnorm(0,0.001)
  }
  
  #Association parameter
  #----------------------
  psi[c] ~ dnorm(0,0.001)
  
  
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
  }
}
",fill = TRUE)
sink()
#Data list for Lagged JM
dat_list_Cluster_JM =list('nId'=nId,'nCluster'=nCluster, 'ones'=rep(1,nId) ,'Const'=10^9, 
                          'y'=y,'y_s'=y_s, 'offset'=offset,'offset_s'=offset_s , 'BS_S'=BS_S , 'BS_T'=BS_T ,
                          'X'=X,'ncX'=ncol(X),'XS'=XS,'ncXS'=ncol(XS),'XT'=XT,'ncXT'=ncol(XT),'X_fix'=X_fix, 'ncF'=ncol(X_fix), 
                          'P'=P , 'wk_r'=wk_r ,  
                          'priorMean.gam.h0'=seq(-8.5,0,len =ncol(BS_S)) , 'M.gam.h0'=M.gam.h0 ,
                          'a.gam.h0'=1 , 'b.gam.h0'=0.05)


Params_save=c('u_L',
              'beta_LYM',
              'sig_LYM',
              'D_L',
              'psi','theta','gam.h0', #Hazard
              'pi','clus','phi' #Cluster
)

Sys.time()
#Run Lagged JM
jagsfit_Cluster_JM_EM <- jags(data=dat_list_Cluster_JM, 
                                  inits=NULL, 
                                  parameters.to.save=Params_save,
                                  model.file=paste(path_out,'Cluster_JM_EM.model.txt', sep="/"),parallel = TRUE,
                                  n.chains=1, n.iter=12000, n.burnin=2000,
                                  n.thin=10,n.adapt = 2000,
                                  DIC=TRUE)


saveRDS(jagsfit_Cluster_JM_EM,paste(path_out,"/Cluster_JM_EM.rds",sep=""))
#-------------------------------------------------------------------------------------------#





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

#Input path
path=paste(fp,"/Model_3/Exp_work/DATA",sep="")
#Output path
path_out=paste(fp,"/Model_3/Exp_work/JAGs_output",sep="")
#dir.create(path_out)


#Pool all JAGS output according to chain values

get_all_jags=function(chain)
{
  jags_list=vector(mode="list", length=1)
  for(x in 1:1)
  {
    l=NULL
    l=readRDS(paste(path_out,paste0("/Simple_Cluster_JM1.rds",collapse=""),sep=""))
    jags_list[[x]]<- l$samples[[chain]]
    x=x+1
  }
  return(jags_list)
}




chain=1
sum_stat="mean"
pe=function(chain, statistic)
{
  
  #Pooled sample
  pool<-get_all_jags(chain)
  
  #Reference of size of parameters
  jagsfits<- readRDS(paste(path_out,paste0("/Simple_Cluster_JM1.rds",collapse=""),sep=""))
  
  #Compute summary statistics
  summary_stat=function(vec, sum_stat)
  {
    z=NULL
    if(sum_stat=="mean")
    {
      z=round(mean(vec),3)
    }
    if(sum_stat=="median")
    {
      z=round(quantile(vec, probs=c(0.5), names=F),3)
    }
    if(sum_stat=="sd")
    {
      z=round(sd(vec),3)
    }
    if(sum_stat=="mode")
    {
      getmode <- function(v) 
      {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
      }
      
      z=getmode(vec)
    }
    if(sum_stat=="significance")
    {
      z=ifelse(prod(quantile(vec, probs=c(0.025,0.975),names=F))<0,"No","Yes") 
    }
    if(sum_stat=="CI")
    {
      z1=round(quantile(vec, probs=c(0.025,0.975),names=F),3)
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
ci=pe(chain=1, statistic = "CI")

#Significance
signif=pe(chain=2, statistic ="significance" )
#head(X_fix)
colnames(signif$phi)=colnames(X_fix[,c(-1)])
signif$phi

#Mean
pos_mean=pe(chain=1, statistic ="mean" )
pos_mean$phi
pos_mean$pi
pos_mean$beta

#Median
pos_median=pe(chain=2, statistic ="median" )

#Mode for discrete
pos_mode=pe(chain=2, statistic="mode")
table(pos_mode$clus)
#df=data.frame("id"=1:236,"Cluster"=pos_mode$clus,"subcluster"=pos_mode$subclus ) #subcluster is for DP
#round(table(df[df$Cluster==1,"subcluster"])*100/nrow(df[df$Cluster==1,]))

#sd
pos_sd=pe(chain=1, statistic ="sd" )


#To view posterior density
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

pd(chain=1,param="phi[1,12]")

