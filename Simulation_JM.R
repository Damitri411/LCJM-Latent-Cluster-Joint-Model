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
library(monomvn)
library(emdbook)
library(JMbayes)


#Parameters
n=250
n_i=10
nCluster=2
Clus_prob=c(0.6,0.4)


#betas
beta_y1<-matrix(rep(NA,5*nCluster),nrow=5)
rownames(beta_y1)<-c("Intercept","Week","Week2","z1","z2")
beta_y1[,1]<-c(1,2,-1,1.5,0.5)
beta_y1[,2]<-c(0.5,-1,1,-1,-1.5)

beta_y2<-c(-1,1,-1,0,1)
beta_y3<-c(1.2,-1,1,1,0)
names(beta_y2)=names(beta_y3)=rownames(beta_y1)

#Error sigma_sq
sig_y1<-c(1,2)
sig_y2<-1.6
sig_y3<-2

#Covariance matrix
D_y1<-array(rep(NA,2*2*nCluster), dim = c(2,2,nCluster))
D_y1[,,1]<-matrix(c(2,1,1,2),nrow=2)
D_y1[,,2]<-matrix(c(2,-1,-1,4),nrow=2)

D_y23<-rwish(5, diag(4))

#psi
psi_y1<-c(-1,1)
psi_y2<-c(0,1)
psi_y3<-c(1,0)

#Baseline haz
lambda0<-c(0.2,1)
#Baseline cov
theta<-matrix(rep(NA,2*nCluster),nrow=2)
rownames(theta)<-c("z1","z2")
theta[,1]=c(1,-1)
theta[,2]=c(-2,1)


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

#Generate Data
data_prod=function(trial)
{
  clus<-sapply(1:n, function(x) which.max(rmultinom(n=1, size=1, prob=Clus_prob)))
  #Censor and fixed covs
  Censor<-sample(x=4,size=n, replace=T)+60
  z1<-rnorm(n=n,mean=1,sd=1)
  z2<-rbinom(n=n, size=1,prob=0.45)
  z3<-runif(n=n, min=1,max=2)
  z4<-sapply(1:n, function(x) which.max(rmultinom(n=1, size=1, prob=c(0.2,0.5,0.3))))
  Week<-rep(NA,n_i*n)
  for(i in 1:n)
  {
    t1<-NULL
    t1<-cumsum(c(0,sapply(1:(n_i-1), function(x) which.max(rmultinom(n=1, size=1, prob=c(0.2,0.6,0.2)))))+1)
    Week[((i-1)*n_i +1):(n_i*i)]<-t1
    
  }
  max_Week<-max(Week,Censor)
  #Create data frame
  df<-data.frame(id=rep(1:n, each=n_i),time_id=rep(1:n_i,n),
                 Intercept=1,Week=Week/max_Week,Week2=(Week/max_Week)^2,
                 z1=rep(z1,each=n_i),
                 z2=rep(z2,each=n_i),
                 Censor=rep(Censor/max_Week,each=n_i),
                 clus=rep(clus,each=n_i),
                 y1=NA, y2=NA,y3=NA , Time=NA, Relapse=NA)
  offset<-(0:n)*n_i +1
  u_y1<-matrix(rep(NA,n*2), nrow=n)
  u_y23<-matrix(rep(NA,n*4), nrow=n)
  Surv<-matrix(rep(NA,max_Week*n), nrow=n)
  Time<-rep(NA,n)
  relapse<-rep(NA,n)
  for(i in 1:n)
  {
    #Random Effect
    u_y1[i,]<-mvrnorm(n=1,mu=c(0,0), Sigma=D_y1[,,clus[i]])
    u_y23[i,]<-mvrnorm(n=1,mu=c(0,0,0,0), Sigma=D_y23)
    
    #Longitudinal part
    for(j in offset[i]:(offset[i+1]-1))
    {
      df$y1[j]<-((u_y1[i,1] + u_y1[i,2]*df$Week[j]) + 
                   as.numeric((beta_y1[,clus[i]])%*%t(df[j,c(3:7)]))+
                   rnorm(n=1, mean=0, sd=sqrt(sig_y1[clus[i]])))
      
      df$y2[j]<-((u_y23[i,1] + u_y23[i,2]*df$Week[j]) + 
                   as.numeric(beta_y2%*%t(df[j,c(3:7)]))+
                   rnorm(n=1, mean=0, sd=sqrt(sig_y2)))
      
      df$y3[j]<-((u_y23[i,3] + u_y23[i,4]*df$Week[j]) + 
                   as.numeric(beta_y3%*%t(df[j,c(3:7)]))+
                   rnorm(n=1, mean=0, sd=sqrt(sig_y3)))
    }
    
    
    Survival_func=function(i, t)
    {
      t_gk<-((Gk_points$cuts+1)*(t-1)/2 +1)/max_Week
      
      haz=function(s)
      {
        mu_y1<-NA
        mu_y2<-NA
        mu_y3<-NA
        mu_y1<-(u_y1[i,1] + u_y1[i,2]*s + as.numeric(beta_y1[,clus[i]]%*%c(1,s,s^2,z1[i],z2[i])))
        mu_y2<-(u_y23[i,1] + u_y23[i,2]*s + as.numeric(beta_y2%*%c(1,s,s^2,z1[i],z2[i])))
        mu_y3<-(u_y23[i,3] + u_y23[i,4]*s + as.numeric(beta_y3%*%c(1,s,s^2,z1[i],z2[i])))
        
        z<-(mu_y1*psi_y1[clus[i]]+mu_y2*psi_y2[clus[i]] +mu_y3*psi_y3[clus[i]]+ as.numeric(theta[,clus[i]]%*%c(z1[i],z2[i])) )
        z1<-lambda0[clus[i]]*exp(z)
        return(z1)
      }
      sur<-sum(sapply(1:length(t_gk) , function(j) Gk_points$weights[j]*haz(t_gk[j])))*(t-1)/(max_Week*2)
      sur<-exp(-sur)
    }
    #1-Survival function at times 1:max_Week
    Surv[i,]<-1-c(1,sapply(2:max_Week, function(t) Survival_func(i, t)))
    un<-runif(n=1,min=0,max=1)
    Time[i]<-min((sum(ifelse(un>Surv[i,],1,0))+1),Censor[i])
    relapse[i]<-ifelse(Time[i]<Censor[i],1,0)
    df[offset[i]:(offset[i+1]-1),"Time"]<-rep(Time[i]/max_Week,n_i)
    df[offset[i]:(offset[i+1]-1),"Relapse"]<-rep(relapse[i],n_i)
    
  }
  
  X<-df[,c("id","time_id","Intercept","Week","Week2","z1","z2")]
  
  y<-df[,c("id","time_id","Week","y1","y2","y3")]
  offset<-(0:n)*n_i +1
  
  X_fix<-df[,c("id","time_id","z1","z2")]
  X_fix=X_fix[!duplicated(X_fix$id),]
  X_fix=X_fix[with(X_fix, order(id)), ]
  rownames(X_fix)=NULL
  
  y_s<-df[,c("id","Time","Relapse")]
  y_s=y_s[!duplicated(y_s$id),]
  y_s=y_s[with(y_s, order(id)), ]
  rownames(y_s)=NULL
  
  offset_s=(0:(n))*15+1
  
  P<-(y_s$Time -min(y$Week))*0.5
  wk <- Gk_points$weights
  wk_r <- rep(Gk_points$weights, n)
  
  XT<-data.frame(id=1:n,time_id=16,
                 Intercept=1, Week=NA , Week2=NA )
  XT$Week<-y_s$Time
  XT$Week2<-XT$Week^2
  XT<-merge(x=XT, y=X_fix[,c("id","z1","z2")], by=c("id"), all.x=T)
  rownames(XT)<-NULL
  
  XS<-data.frame(id=rep(1:n, each=nrow(Gk_points)),time_id=rep(1:nrow(Gk_points), n),
                 Intercept=1, Week=NA , Week2=NA )
  XS<-merge(x=XS, y=y_s[,c("id","Time")], by=c("id"), all.x=T)
  XS$cuts<-(Gk_points$cuts[XS$time_id]+1)*0.5
  XS$Week<-(XS$Time-min(y$Week))*XS$cuts +min(y$Week)
  XS$Week2<-XS$Week^2
  XS<-XS[,!colnames(XS)%in%c("cuts","Time")]
  XS<-merge(x=XS, y=X_fix[,c("id","z1","z2")], by=c("id"), all.x=T)
  rownames(XS)<-NULL
  
  True_Random_Effect<-list(u_y1, u_y23)
  #Data list 
  dat_list =list('n'=n,'nCluster'=NA, 'ones'=rep(1,n), 'zeros'=rep(0,6) ,'Const'=10^9, 'True_Random_Effect'=True_Random_Effect,'True_clus'=clus,'max_Week'=max_Week,
                 'y'=y,'y_s'=y_s, 'offset'=offset,'offset_s'=offset_s , 
                 'X'=X,'ncX'=ncol(X),'XS'=XS,'ncXS'=ncol(XS),'XT'=XT,'ncXT'=ncol(XT),'X_fix'=X_fix, 'ncF'=ncol(X_fix), 
                 'P'=P , 'wk_r'=wk_r , 
                 'alpha'=c(1,1,1,1,1), 
                 'Omega_y23'=diag(rep(1,4)) ,'Omega_y1'=diag(rep(1,2)) )
  return(dat_list)
}


df<-data_prod(1)
table(df$y_s[,2]/min(df$X[,"Week"]))
table(df$y_s[,3])
table(df$True_clus, df$y_s[,3])


#Input path
#File Path
fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"
path=paste(fp,"/Cluster_JM/Exp_work/DATA",sep="")
#Output path
path_out=paste(fp,"/JAGS_Output/Cluster_JM",sep="")
#---------------------------------------------------------------------------------------------------------#
# JM Cluster
#---------------------------------------------------------------------------------------------------------#
sink(paste(path_out,"Cluster1.model.txt", sep="/"))
cat("
model {

for(i in 1:n) 
{
 
 #Random Effect
 #-------------
  u_y23[i,1:4] ~ dmnorm(zeros[1:4],D_inv_y23[1:4,1:4]) #Global random intercept and slope of y2 , y3
  
 for(cl in 1:nCluster)
  {
   u_y1[i,1:2,cl] ~ dmnorm(zeros[5:6],D_inv_y1[1:2,1:2,cl]) #Cluster wise random intercept and slope  of y1
  }
  
 
 #Longitudinal Process
 #---------------------
for( j in offset[i]:(offset[i+1]-1))
 {
   
  mu_y2[j]= u_y23[i,1] + u_y23[i,2]*X[j,4] + inprod(beta_y2[1:(ncX-2)], X[j,3:ncX])
  y[j,5] ~ dnorm(mu_y2[j],prec_y2) 
  
  #Mean y3
  mu_y3[j]= u_y23[i,3] + u_y23[i,4]*X[j,4]  + inprod(beta_y3[1:(ncX-2)], X[j,3:ncX])
  y[j,6] ~ dnorm(mu_y3[j],prec_y3) 
 
  
  #Mean y1
  for(cl in 1:nCluster)
  {
    mu_y1[j,cl]= u_y1[i,1,cl] + u_y1[i,2,cl]*X[j,4] + inprod(beta_y1[1:(ncX-2),cl], X[j,3:ncX]) 
  
  #Log density cluster wise of y1
  ld_y1[j,cl] <- logdensity.norm(y[j,4], mu_y1[j,cl], prec_y1[cl])  
  }
  
  } #loop of j
  
  

   #mu_y2(t) and mu_y3(t) in hazard (Survival part 1) 
   #----------------------------------------------------
   mu_Ty2[i] =u_y23[i,1] + u_y23[i,2]*XT[i,4]  + inprod(beta_y2[1:(ncXT-2)], XT[i,3:ncXT])
   mu_Ty3[i] =u_y23[i,3] + u_y23[i,4]*XT[i,4]  + inprod(beta_y3[1:(ncXT-2)], XT[i,3:ncXT])
   
   #mu_y2(t) and mu_y3(t) in Survival probability (Survival part 2) 
   #-----------------------------------------------------------------
   for(l in offset_s[i]:(offset_s[i+1]-1) )
  {
    mu_Sy2[l] =u_y23[i,1] + u_y23[i,2]*XS[l,4] + inprod(beta_y2[1:(ncXS-2)], XS[l,3:ncXS])
    mu_Sy3[l] =u_y23[i,3] + u_y23[i,4]*XS[l,4] + inprod(beta_y3[1:(ncXS-2)], XS[l,3:ncXS])
  }
   
   for(cl in 1:nCluster)
  {
  
  #Survival part 1
  #---------------
 
   mu_Ty1[i,cl] =u_y1[i,1,cl] + u_y1[i,2,cl]*XT[i,4]  + inprod(beta_y1[1:(ncXT-2),cl], XT[i,3:ncXT]) 
                               
   log.hT[i,cl]= ( log.h0[cl] + psi_y1[cl]*mu_Ty1[i,cl]
                  + psi_y2[cl]*mu_Ty2[i] + psi_y3[cl]*mu_Ty3[i]
                  + inprod(theta[1:(ncF-2),cl],X_fix[i,3:ncF])
                  ) 
 
 
 #Survival part 2
 #---------------
 
 for(l in offset_s[i]:(offset_s[i+1]-1) )
  {
   
    mu_Sy1[l,cl] =u_y1[i,1,cl]  + u_y1[i,2,cl]*XS[l,4] + inprod(beta_y1[1:(ncXS-2),cl], XS[l,3:ncXS]) 
                                             
    SurvLong[l,cl] = wk_r[l]* ( exp( ( log.h0[cl] +psi_y1[cl]*mu_Sy1[l,cl] + psi_y2[cl]*mu_Sy2[l] + psi_y3[cl]*mu_Sy3[l]) ) )
  }
   #Baseline covariates in Survival part 2
   #--------------------------------------
   surv_base[i,cl] = exp(inprod(theta[1:(ncF-2),cl],X_fix[i,3:ncF]))
   
   
    #The likelihood clusterwise
    lik[i,cl] <-  exp((logdensity.mnorm(u_y1[i,1:2,cl], zeros[5:6], D_inv_y1[1:2,1:2,cl]) #Random Effect of y1
                       +sum(ld_y1[offset[i]:(offset[i+1]-1),cl])   #Longitudinal y1
                       +(y_s[i,3]*log.hT[i,cl])   #Hazard
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


#Global priors for y2, y3
#---------------------------
tau_y23[1:4,1]~dmnorm(zeros[1:4],Omega_y23)
tau_y23[1:4,2]~dmnorm(zeros[1:4],Omega_y23)
tau_y23[1:4,3]~dmnorm(zeros[1:4],Omega_y23)
tau_y23[1:4,4]~dmnorm(zeros[1:4],Omega_y23)
tau_y23[1:4,5]~dmnorm(zeros[1:4],Omega_y23)
D_inv_y23[1:4,1:4]<-tau_y23[1:4,1:5] %*% t(tau_y23[1:4,1:5]) + Omega_y23[1:4,1:4]*0.001

#Betas of y2 and y3
#--------------------

for(h in 1:(ncX-2))
  {
  beta_y2[h]~dnorm(0,0.001)
  beta_y3[h]~dnorm(0,0.001)
  }


#Error variance of y2 and y3
#-----------------------------
prec_y2 ~ dgamma(0.001,0.001)
sig_y2 = 1/prec_y2
prec_y3 ~ dgamma(0.001,0.001)
sig_y3 = 1/prec_y3


#pi with dirichlet prior
#-----------------------
pi[1:nCluster] ~ ddirch(alpha[1:nCluster])
 

  for(c in 1:nCluster)
  {  

  #Error Variance
  #--------------
  prec_y1[c] ~ dgamma(0.001,0.001)
  sig_y1[c] = 1/prec_y1[c]
  

  #Variance for Random Effect of y1 
  #---------------------------------
  tau_y1[1:2,1,c]~dmnorm(zeros[1:2],(Omega_y1[1:2,1:2]))
  tau_y1[1:2,2,c]~dmnorm(zeros[1:2],(Omega_y1[1:2,1:2]))
  tau_y1[1:2,3,c]~dmnorm(zeros[1:2],(Omega_y1[1:2,1:2]))
  D_inv_y1[1:2,1:2,c]<-tau_y1[1:2,1:3,c] %*% t(tau_y1[1:2,1:3,c]) + Omega_y1[1:2,1:2]*0.001

  
  #Beta
  #-----
  for(h in 1:(ncX-2))
  {
  beta_y1[h,c]~dnorm(0,0.001)
  }
  
  #Association parameter
  #----------------------
  psi_y1[c] ~ dnorm(0,0.001)
  psi_y2[c] ~ dnorm(0,0.001)
  psi_y3[c] ~ dnorm(0,0.001)
 
  #Theta
  #-----
  for(h in 1:(ncF-2))
  {
  theta[h,c]~dnorm(0,0.001)
  }
  
  #Prior for Baseline hazard
  #--------------------------
  log.h0[c] ~  dmnorm(0,0.001)
  
 
  
  }#End of cluster parameters
  
  
}
",fill = TRUE)
sink()


Params_save=c('u_y1','u_y23', #Random effects
              'beta_y1','beta_y2','beta_y3',  #The betas
              'sig_y1', 'sig_y2' ,'sig_y3', #Error variances
              'D_inv_y1','D_inv_y23', #Random effect covariance matrices inverse
              'psi_y1','psi_y2','psi_y3','theta','log.h0', #Elements of Hazard
              'pi','pi_star','clus', #Cluster properties
              'lik' #Cluster wise likelihood
)

#Number of cluster
List_of_criterion=function(jfit, dat_list)
{
  X<-dat_list$X
  XS<-dat_list$XS
  XT<-dat_list$XT
  X_fix<-dat_list$X_fix
  y<-dat_list$y
  y_s<-dat_list$y_s
  nCluster<-length(jfit$mean$pi)
  ncX<-ncol(X)
  ncXS<-ncol(XS)
  ncXT<-ncol(XT)
  ncF<-ncol(X_fix)
  nId<-dat_list$n
  wk_r<-dat_list$wk_r
  P<-dat_list$P
  offset<-dat_list$offset
  offset_s<-dat_list$offset_s
  
  
  #likelihood for a any particular theta
  likelihood=function(jfit1,nCluster)
  {
    #Means
    mu_y2<-rep(NA, nrow(X))
    mu_y3<-rep(NA,nrow(X))
    mu_y1<-matrix(rep(NA,nrow(X)*nCluster),ncol=nCluster)
    
    #Random Effect
    ld_y23_ran<-rep(NA,nId)
    ld_y1_ran<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    
    #longitudinal log likelihood
    ld_y23<-rep(NA, nrow(X))
    ld_y1<-matrix(rep(NA,nrow(X)*nCluster),ncol=nCluster)
    
    #Mean in hazard
    mu_Ty2<-rep(NA,nId)
    mu_Ty3<-rep(NA,nId)
    mu_Ty1<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    
    #The hazard
    log.hT<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    
    #Means in Survival part
    mu_Sy2<-rep(NA,nrow(XS))
    mu_Sy3<-rep(NA,nrow(XS))
    mu_Sy1<-matrix(rep(NA,nrow(XS)*nCluster),ncol=nCluster)
    
    #Survival probability components without baseline covs
    SurvLong<-matrix(rep(NA,nrow(XS)*nCluster),ncol=nCluster)
    
    #Baseline covariates in survival
    surv_base<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    
    #Cluster wise likelihood
    lik<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    
    #Total likelihood
    total_lik<-rep(NA,nId)
    
    
    for(i in 1:n) 
    {
      
      #Random Effect
      #-------------#
      ld_y23_ran[i]<-dmvnorm(x=jfit1$u_y23[i,1:4], mu=c(0,0,0,0), Sigma=solve(jfit1$D_inv_y23[1:4,1:4]), log = FALSE, tol = 1e-07)
      
      for(cl in 1:nCluster)
      {
        ld_y1_ran[i,cl]<-dmvnorm(x=jfit1$u_y1[i,1:2,cl], mu=c(0,0), Sigma=solve(jfit1$D_inv_y1[1:2,1:2,cl]), log = TRUE, tol = 1e-07)
      }
      
      
      #Longitudinal Process
      #---------------------#
      for( j in offset[i]:(offset[i+1]-1))
      {
        
        mu_y2[j]= jfit1$u_y23[i,1] + jfit1$u_y23[i,2]*X[j,4] + as.numeric(jfit1$beta_y2[1:(ncX-2)]%*% t(X[j,3:ncX]))
        mu_y3[j]= jfit1$u_y23[i,3] + jfit1$u_y23[i,4]*X[j,4] + as.numeric(jfit1$beta_y3[1:(ncX-2)]%*%t(X[j,3:ncX]))
        
        ld_y23[j]<-(dnorm(x=y[j,5],mean=mu_y2[j],sd=sqrt(jfit1$sig_y2), log = TRUE)+
                      dnorm(x=y[j,6],mean=mu_y3[j],sd=sqrt(jfit1$sig_y3), log = TRUE))
        
        #Mean y1
        for(cl in 1:nCluster)
        {
          mu_y1[j,cl]= jfit1$u_y1[i,1,cl] + jfit1$u_y1[i,2,cl]*X[j,4] + as.numeric(jfit1$beta_y1[1:(ncX-2),cl]%*%t(X[j,3:ncX]) )
          #Log density cluster wise of y1
          ld_y1[j,cl]<-dnorm(x=y[j,4],mean=mu_y1[j,cl],sd=sqrt(jfit1$sig_y1[cl]), log = TRUE)
        }
        
      } #loop of j
      
      
      
      #mu_y2(t) and mu_y3(t) in hazard (Survival part 1) 
      #----------------------------------------------------#
      mu_Ty2[i] =jfit1$u_y23[i,1] + jfit1$u_y23[i,2]*XT[i,4]  + as.numeric(jfit1$beta_y2[1:(ncXT-2)]%*%t(XT[i,3:ncXT]))
      mu_Ty3[i] =jfit1$u_y23[i,3] + jfit1$u_y23[i,4]*XT[i,4]  + as.numeric(jfit1$beta_y3[1:(ncXT-2)]%*%t(XT[i,3:ncXT]))
      
      #mu_y2(t) and mu_y3(t) in Survival probability (Survival part 2) 
      #-----------------------------------------------------------------#
      for(l in offset_s[i]:(offset_s[i+1]-1) )
      {
        mu_Sy2[l] =jfit1$u_y23[i,1] + jfit1$u_y23[i,2]*XS[l,4] + as.numeric(jfit1$beta_y2[1:(ncXS-2)]%*%t( XS[l,3:ncXS]))
        mu_Sy3[l] =jfit1$u_y23[i,3] + jfit1$u_y23[i,4]*XS[l,4] + as.numeric(jfit1$beta_y3[1:(ncXS-2)]%*%t( XS[l,3:ncXS]))
      }
      
      for(cl in 1:nCluster)
      {
        
        #Survival part 1
        #---------------#
        
        mu_Ty1[i,cl] =jfit1$u_y1[i,1,cl] + jfit1$u_y1[i,2,cl]*XT[i,4]  + as.numeric(jfit1$beta_y1[1:(ncXT-2),cl]%*%t( XT[i,3:ncXT]) )
        
        log.hT[i,cl]= ( jfit1$log.h0[cl] + jfit1$psi_y1[cl]*mu_Ty1[i,cl]
                        + jfit1$psi_y2[cl]*mu_Ty2[i] + jfit1$psi_y3[cl]*mu_Ty3[i]
                        + as.numeric(jfit1$theta[1:(ncF-2),cl]%*%t(X_fix[i,3:ncF]))
        ) 
        
        
        #Survival part 2
        #---------------#
        
        for(l in offset_s[i]:(offset_s[i+1]-1) )
        {
          
          mu_Sy1[l,cl] =jfit1$u_y1[i,1,cl]  + jfit1$u_y1[i,2,cl]*XS[l,4] + as.numeric(jfit1$beta_y1[1:(ncXS-2),cl]%*%t( XS[l,3:ncXS])) 
          
          SurvLong[l,cl] = wk_r[l]* ( exp( ( jfit1$log.h0[cl] +jfit1$psi_y1[cl]*mu_Sy1[l,cl] + jfit1$psi_y2[cl]*mu_Sy2[l] + jfit1$psi_y3[cl]*mu_Sy3[l]) ) )
        }
        #Baseline covariates in Survival part 2
        #--------------------------------------#
        surv_base[i,cl] = exp(as.numeric(theta[1:(ncF-2),cl]%*%t(X_fix[i,3:ncF])))
        
        
        #The likelihood clusterwise
        lik[i,cl] <-  exp((ld_y1_ran[i,cl]                              #Random Effect of y1
                           +sum(ld_y1[offset[i]:(offset[i+1]-1),cl])    #Longitudinal y1
                           +sum(ld_y23[offset[i]:(offset[i+1]-1)])        #Longitudinal y2 , y3
                           +ld_y23_ran[i]                                 #Random Effect of y2, y3 
                           +(y_s[i,3]*log.hT[i,cl])   #Hazard
                           -( P[i]*surv_base[i,cl]*sum(SurvLong[offset_s[i]:(offset_s[i+1]-1),cl]) ) )) #Survival function
        
        #Cluster wise likelihood times their cluster probs
        #-------------------------------------------------#
        total_lik[i]<-sum(lik[i,1:nCluster]*jfit$sims.list$pi[1:nCluster])
        
      }
      }
    return(total_lik)
  }
  
  #likelihood iteration wise
  likelihood_computation=function(jfit_samp,ite,nCluster)
  {
    
    #Means
    mu_y2<-rep(NA, nrow(X))
    mu_y3<-rep(NA,nrow(X))
    
    #Random Effect
    ld_y23_ran<-rep(NA,nId)
    ld_y1_ran<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    
    #longitudinal log likelihood
    ld_y23<-rep(NA, nrow(X))
    
    lik<-matrix(rep(NA,nId*nCluster),ncol=nCluster)
    total_lik<-rep(NA,nId)
    
    
    for(i in 1:nId)
    {  
      #Random Effect of y2, y3
      ld_y23_ran[i]<-dmvnorm(x=jfit_samp$u_y23[ite,i,1:4], mu=c(0,0,0,0), Sigma=solve(jfit_samp$D_inv_y23[ite,1:4,1:4]), log = FALSE, tol = 1e-07)
      
      #Longitudinal part
      for(j in offset[i]:(offset[i+1]-1))
      {
        mu_y2[j]<- as.numeric(jfit_samp$u_y23[ite,i,1] + jfit_samp$u_y23[ite,i,2]*X[j,4] + as.numeric((jfit_samp$beta_y2[ite,1:(ncX-2)]) %*% t(X[j,3:ncX])))
        mu_y3[j]<- as.numeric(jfit_samp$u_y23[ite,i,3] + jfit_samp$u_y23[ite,i,4]*X[j,4] + as.numeric((jfit_samp$beta_y3[ite,1:(ncX-2)]) %*% t(X[j,3:ncX])))
        
        #Likelihood of y2, y3 (Longitudinal)
        ld_y23[j]<-(dnorm(x=y[j,5],mean=mu_y2[j],sd=sqrt(jfit_samp$sig_y2[ite]), log = TRUE)+
                      dnorm(x=y[j,6],mean=mu_y3[j],sd=sqrt(jfit_samp$sig_y3[ite]), log = TRUE))
        
      }
      
      for( cl in 1:nCluster)
      {
        lik[i,cl] <- jfit_samp$lik[ite,i,cl] * jfit$sims.list$pi[ite,cl]      
      }
      
      total_lik[i]<-sum(lik[i,1:nCluster]) *exp((sum(ld_y23[offset[i]:(offset[i+1]-1)]) +ld_y23_ran[i]))
      
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
  
  paramater_vec<-c('u_y1','u_y23', #Random effects
                   'beta_y1','beta_y2','beta_y3',  #The betas
                   'sig_y1', 'sig_y2' ,'sig_y3', #Error variy2es
                   'D_inv_y1','D_inv_y23', #Random effect covariy2e matrices inverse
                   'psi_y1','psi_y2','psi_y3','theta','gam.h0','tau.gam.h0' #Elements of Hazard
  )
  num_parameters<-sum(sapply(1:length(paramater_vec), function(i) length(jfit$mean[[paramater_vec[i]]])))
  
  #AIC
  aic_comp<- -2*sum(log(likelihood(jfit1=jfit$mean,nCluster))) + 2*num_parameters
  
  #BIC
  bic_comp<- -2*sum(log(likelihood(jfit1=jfit$mean,nCluster))) + num_parameters*log(nId)
  
  #Misclassification rate
  misc<-NULL
  if(nCluster==2)
  {
    clus_pred<-sapply(1:nId,function(i)which.max(jfit$mean$pi_star[i,]))
    temp_tab<-as.matrix(table(dat_list$True_clus, clus_pred))
    misc=round(min((temp_tab[1,1]+temp_tab[2,2]),(temp_tab[1,2]+temp_tab[2,1]))*100/sum(temp_tab))

  }
  else
  {
    misc=NA
  }
  
  #Minimum cluster percentage
  min_clust_comp<-round(min(table(sapply(1:nId,function(i)which.max(jfit$mean$pi_star[i,])))*100/nId),3)
  #Criterion vector
  Criterion_vector=c(nCluster,log_liki,dic_comp,aic_comp,bic_comp,min_clust_comp,misc)
  names(Criterion_vector)=c("Clucters_num","log_likelihood","DIC","AIC","BIC","Minimum_Cluster_size","Misclassfication")
  return(Criterion_vector)
  
}

#Run Simulation for cluster =2,3,4
diagnostics=NULL
Basic_diagnostic=NULL
i=1
for(i in 1:200)
{
  #Generate data and data list
  dat_list<-NULL
  dat_list<-data_prod(i)
  
  #save data_list
  saveRDS(dat_list,paste(path_out,"/dat_list_",i,".rds",sep=""))
  
  #Run model for 2 clusters
  dat_list2<-dat_list[!names(dat_list)%in%c("True_Random_Effect","True_clus","max_Week")]
  dat_list2$nCluster<-2
  Sys.time()
  jagsfit_JM2 <- jags(data=dat_list2, 
                      inits=NULL, 
                      parameters.to.save=Params_save,
                      model.file=paste(path_out,'Cluster1.model.txt', sep="/"),parallel = TRUE,
                      n.chains=1, n.iter=7000, n.burnin=1000,
                      n.thin=10,n.adapt = 1000,
                      DIC=TRUE)
  
  #Run model for 3 clusters
  dat_list3<-dat_list2
  dat_list3$nCluster<-3
  Sys.time()
  jagsfit_JM3 <- jags(data=dat_list3, 
                      inits=NULL, 
                      parameters.to.save=Params_save,
                      model.file=paste(path_out,'Cluster1.model.txt', sep="/"),parallel = TRUE,
                      n.chains=1, n.iter=7000, n.burnin=1000,
                      n.thin=10,n.adapt = 1000,
                      DIC=TRUE)
  
  
  #Run model for 4 clusters
  dat_list4<-dat_list2
  dat_list4$nCluster<-4
  Sys.time()
  jagsfit_JM4 <- jags(data=dat_list4, 
                      inits=NULL, 
                      parameters.to.save=Params_save,
                      model.file=paste(path_out,'Cluster1.model.txt', sep="/"),parallel = TRUE,
                      n.chains=1, n.iter=7000, n.burnin=1000,
                      n.thin=10,n.adapt = 1000,
                      DIC=TRUE)
  
  #Model wise AIC, BIC, DIC  values
  diagn<-NULL
  diagn<-data.frame(rbind(c(i,List_of_criterion(jfit=jagsfit_JM2, dat_list=dat_list)),
                            c(i,List_of_criterion(jfit=jagsfit_JM3, dat_list=dat_list)),
                            c(i,List_of_criterion(jfit=jagsfit_JM4, dat_list=dat_list))))
  
  diagnostics<-rbind(diagnostics,diagn)
  
  
  #True parameter, bias, CI length, Coverage Probability for cluster=2 only
  jfit<-NULL
  jfit<-jagsfit_JM2
  clus_pred<-sapply(1:dat_list$n,function(i)which.max(jfit$mean$pi_star[i,]))
  temp_tab<-as.matrix(table(dat_list$True_clus, clus_pred))
  
  rev_ind=ifelse((temp_tab[1,1]+temp_tab[2,2])<(temp_tab[1,2]+temp_tab[2,1]), "reversed","original")
  
  #Did not include the random effect for variance covariance matrices
  all_clus<-NULL
  for(clus_p in 1:2)
  {
    clus<-ifelse(rev_ind=="reversed",3-clus_p,clus_p)
    
    bias_beta_y1=bias_sig_y1=bias_psi_y1=bias_psi_y1=bias_psi_y2=bias_theta=bias_pi=NULL
    
    #beta_y1
    bias_beta_y1<-data.frame(Param_name="beta_y1",
                             parameters=paste("beta_y1[",rownames(beta_y1)," , ",clus_p,"]",sep=""),
                             cluster=clus_p,
                             True_value=beta_y1[,clus_p],
                             bias=jfit$mean$beta_y1[,clus]-beta_y1[,clus_p],
                             CI_length=jfit$q97.5$beta_y1[,clus]-jfit$q2.5$beta_y1[,clus],
                             Coverage_ind=sapply(1:nrow(beta_y1), function(j) ifelse((beta_y1[j,clus_p]-jfit$q2.5$beta_y1[j,clus])*(jfit$q97.5$beta_y1[j,clus]-beta_y1[j,clus_p])>=0 ,1,0 )))
    #sig_y1
    bias_sig_y1<-data.frame(Param_name="sig_y1",
                            parameters=paste("sig_y1[",clus_p,"]",sep=""),
                            cluster=clus_p,
                            True_value=sig_y1[clus_p],
                            bias=jfit$mean$sig_y1[clus]-sig_y1[clus_p],
                            CI_length=jfit$q97.5$sig_y1[clus]-jfit$q2.5$sig_y1[clus],
                            Coverage_ind= ifelse((sig_y1[clus_p]-jfit$q2.5$sig_y1[clus])*(jfit$q97.5$sig_y1[clus]-sig_y1[clus_p])>=0 ,1,0 )
    )
    #psi_y1
    bias_psi_y1<-data.frame(Param_name="psi_y1",
                            parameters=paste("psi_y1[",clus_p,"]",sep=""),
                            cluster=clus_p,
                            True_value=psi_y1[clus_p],
                            bias=jfit$mean$psi_y1[clus]-psi_y1[clus_p],
                            CI_length=jfit$q97.5$psi_y1[clus]-jfit$q2.5$psi_y1[clus],
                            Coverage_ind= ifelse((psi_y1[clus_p]-jfit$q2.5$psi_y1[clus])*(jfit$q97.5$psi_y1[clus]-psi_y1[clus_p])>=0 ,1,0 )
    )
    
    bias_psi_y2<-data.frame(Param_name="psi_y2",
                            parameters=paste("psi_y2[",clus_p,"]",sep=""),
                            cluster=clus_p,
                            True_value=psi_y2[clus_p],
                            bias=jfit$mean$psi_y2[clus]-psi_y2[clus_p],
                            CI_length=jfit$q97.5$psi_y2[clus]-jfit$q2.5$psi_y2[clus],
                            Coverage_ind= ifelse((psi_y2[clus_p]-jfit$q2.5$psi_y2[clus])*(jfit$q97.5$psi_y2[clus]-psi_y2[clus_p])>=0 ,1,0 )
    )
    #psi_y3
    bias_psi_y3<-data.frame(Param_name="psi_y3",
                            parameters=paste("psi_y3[",clus_p,"]",sep=""),
                            cluster=clus_p,
                            True_value=psi_y3[clus_p],
                            bias=jfit$mean$psi_y3[clus]-psi_y3[clus_p],
                            CI_length=jfit$q97.5$psi_y3[clus]-jfit$q2.5$psi_y3[clus],
                            Coverage_ind= ifelse((psi_y3[clus_p]-jfit$q2.5$psi_y3[clus])*(jfit$q97.5$psi_y3[clus]-psi_y3[clus_p])>=0 ,1,0 )
    )
    #theta
    bias_theta<-data.frame(Param_name="theta",
                           parameters=paste("theta[",rownames(theta)," , ",clus_p,"]",sep=""),
                           cluster=clus_p,
                           True_value=theta[,clus_p],
                           bias=jfit$mean$theta[,clus]-theta[,clus_p],
                           CI_length=jfit$q97.5$theta[,clus]-jfit$q2.5$theta[,clus],
                           Coverage_ind=sapply(1:nrow(theta), function(j) ifelse((theta[j,clus_p]-jfit$q2.5$theta[j,clus])*(jfit$q97.5$theta[j,clus]-theta[j,clus_p])>=0 ,1,0 )))
    #D_y1
    #bias_D_y1[,,clus_p]<-as.vector(solve(jfit$mean$D_inv_y1[,,clus])-D_y1[,,clus_p])
    #pi
    bias_pi<-data.frame(Param_name="pi",
                        parameters=paste("pi[",clus_p,"]",sep=""),
                        cluster=clus_p,
                        True_value=Clus_prob[clus_p],
                        bias=jfit$mean$pi[clus]-Clus_prob[clus_p],
                        CI_length=jfit$q97.5$pi[clus]-jfit$q2.5$pi[clus],
                        Coverage_ind= ifelse((Clus_prob[clus_p]-jfit$q2.5$pi[clus])*(jfit$q97.5$pi[clus]-Clus_prob[clus_p])>=0 ,1,0 )
    )
    
    all_clus1<-NULL    
    all_clus1<-rbind(bias_beta_y1,bias_sig_y1,bias_psi_y1,bias_psi_y1,bias_psi_y2,bias_theta,bias_pi)
    all_clus<-rbind(all_clus,all_clus1)
  }
  bias_sig_y2=bias_sig_y3=bias_beta_y2=bias_beta_y3=NULL
  #sig_y2
  bias_sig_y2<-data.frame(Param_name="sig_y2",
                          parameters="sig_y2",
                          cluster=NA,
                          True_value=sig_y2,
                          bias=jfit$mean$sig_y2-sig_y2,
                          CI_length=jfit$q97.5$sig_y2-jfit$q2.5$sig_y2,
                          Coverage_ind= ifelse((sig_y2-jfit$q2.5$sig_y2)*(jfit$q97.5$sig_y2-sig_y2)>=0 ,1,0 ))
  
  #sig_y3
  bias_sig_y3<-data.frame(Param_name="sig_y3",
                          parameters="sig_y3",
                          cluster=NA,
                          True_value=sig_y3,
                          bias=jfit$mean$sig_y3-sig_y3,
                          CI_length=jfit$q97.5$sig_y3-jfit$q2.5$sig_y3,
                          Coverage_ind= ifelse((sig_y3-jfit$q2.5$sig_y3)*(jfit$q97.5$sig_y3-sig_y3)>=0 ,1,0 ))
  
  #beta_y2
  bias_beta_y2<-data.frame(Param_name="beta_y2",
                           parameters=paste("beta_y2[",rownames(beta_y1),"]",sep=""),
                           cluster=NA,
                           True_value=beta_y2,
                           bias=jfit$mean$beta_y2-beta_y2,
                           CI_length=jfit$q97.5$beta_y2-jfit$q2.5$beta_y2,
                           Coverage_ind=sapply(1:length(beta_y2), function(j) ifelse((beta_y2[j]-jfit$q2.5$beta_y2[j])*(jfit$q97.5$beta_y2[j]-beta_y2[j])>=0 ,1,0 )))
  
  
  #beta_y3
  bias_beta_y3<-data.frame(Param_name="beta_y3",
                           parameters=paste("beta_y3[",rownames(beta_y1),"]",sep=""),
                           cluster=NA,
                           True_value=beta_y3,
                           bias=jfit$mean$beta_y3-beta_y3,
                           CI_length=jfit$q97.5$beta_y3-jfit$q2.5$beta_y3,
                           Coverage_ind=sapply(1:length(beta_y3), function(j) ifelse((beta_y3[j]-jfit$q2.5$beta_y3[j])*(jfit$q97.5$beta_y3[j]-beta_y3[j])>=0 ,1,0 )))
  
  all_clus<-rbind(all_clus,bias_sig_y2,bias_sig_y3,bias_beta_y2,bias_beta_y3)
  all_clus<-data.frame(all_clus)
  rownames(all_clus)<-NULL
  all_clus$Iteration=i
  
  Basic_diagnostic<-rbind(Basic_diagnostic,all_clus)
  
  i=i+1
  
}
diagnostics<-data.frame(diagnostics)
rownames(diagnostics)<-NULL
colnames(diagnostics)<-c("Iterations","LogLik","DIC","AIC","BIC","Min_cluster_perc","Misclassfication")

Basic_diagnostic<-data.frame(Basic_diagnostic)
rownames(Basic_diagnostic)<-NULL
colnames(Basic_diagnostic)<-c("Param_name","parameters","cluster","True_value","bias","CI_length","Coverage_ind","Iteration")

