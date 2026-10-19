#Cluster number
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
library(latex2exp)


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


#Longitudinal trends
#Library

library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)
#File Path
fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"

#Input path
path=paste(fp,"/Cluster_JM/Exp_work/DATA",sep="")
#Output path
path_out=paste(fp,"/JAGS_Output/Cluster_JM",sep="")
#dir.create(path_out)

jfit<-readRDS(paste(path_out,paste("/JM_comp_Lymph_Cl_",nCluster,".rds",sep=""),sep=""))
time_set<-c(1,seq(10,110,10))/max_week
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
data_X=data.frame(id_new=rep(1:nId,each=length(time_set)),time_id=rep(1:length(time_set),nId),Intercept=1,Week=rep(time_set,nId))
data_X$Week2=data_X$Week^2
data_X$X6_mp=NA
data_X$MTX=NA
data_X=merge(x=data_X, y=X_fix[,!colnames(X_fix)%in%c("Intercept")], by=c("id_new"),all.x=T)

for(i in 1:nrow(data_X))
{
  v<-NULL
  v<-medi(id=data_X$id_new[i], t=data_X$Week[i])
  data_X$X6_mp[i]<-v[1]
  data_X$MTX[i]<-v[2]
  i=i+1
}
mu_lym=data.frame(id_new=rep(rep(1:nId,each=length(time_set)),nCluster),
                  time_id=rep(rep(1:length(time_set),nId),nCluster),
                  Cluster=rep(1:nCluster,each=(nId*(length(time_set)))),
                  Lymph=NA)
mu_lym$Week=time_set[mu_lym$time_id]*max_week


long_cl=function(id,t_id,Cluster)
{
  betas<- jfit$mean$beta_LYM[1:15,Cluster]
  d=subset(data_X, data_X$id_new==id )
  d=subset(d, d$time_id==as.numeric(t_id))
  v=jfit$mean$u_L[id,1,Cluster] + jfit$mean$u_L[id,2,Cluster]*d$Week +
    as.numeric(betas%*%t(d[,c(-1,-2)]))
  
  return(v)
    
}


mu_lym$Lymph=sapply(1:nrow(mu_lym),
function(i)long_cl(id=mu_lym$id_new[i],t_id=mu_lym$time_id[i],Cluster=mu_lym$Cluster[i]))

#Estimated cluster _ id combination and accordingly parameters for that cluster
mu_lym$st_cl=paste(mu_lym$id_new, mu_lym$Cluster, sep="_")
cluster_id_comb=sapply(1:nId, function(i) paste(i,which.max(jfit$mean$pi_star[i,]), sep="_"))
mu_lym<-subset(mu_lym, mu_lym$st_cl%in%cluster_id_comb)
  
dfit<-mu_lym%>%
  group_by(Week,Cluster)%>%
  summarise(Lymph=mean(Lymph))

dfit$Lymph<-dfit$Lymph+7.32449
dfit$Cluster=as.factor(dfit$Cluster)
ggplot(data=dfit,aes(x=Week, y=Lymph, color=Cluster))+
  geom_line(size=0.875)+ scale_color_manual(values=c("red","blue","green"))+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=12 ,color="black", face="bold"), axis.title.y =element_text(size=12 ,color="black", face="bold"))+
  xlab("Week")+ylab("log(LYMPH)")+
  geom_hline(yintercept = log(1664.95), col="black", linetype="dashed")+
  theme(plot.title = element_text(hjust = 0.5))

#ANC, PLT post analysis on the basis of pi_star groups

mu_anc_plt=data.frame(id_new=rep(1:nId,each=length(time_set)),
                      time_id=rep(1:length(time_set),nId),
                      ANC=NA, PLT=NA)
mu_anc_plt$Week=time_set[mu_anc_plt$time_id]*max_week


long_gen_a=function(id,t_id)
{
  betas_a<- jfit$mean$beta_ANC
  d=subset(data_X, data_X$id_new==id )
  d=subset(d, d$time_id==as.numeric(t_id))
  v=jfit$mean$u_AP[id,1] + jfit$mean$u_AP[id,3]*d$Week +
    as.numeric(betas_a%*%t(d[,c(-1,-2)]))
  
  return(v)
  
}
long_gen_p=function(id,t_id)
{
  betas_p<- jfit$mean$beta_PLT
  d=subset(data_X, data_X$id_new==id )
  d=subset(d, d$time_id==as.numeric(t_id))
  v=jfit$mean$u_AP[id,2] + jfit$mean$u_AP[id,4]*d$Week +
    as.numeric(betas_p%*%t(d[,c(-1,-2)]))
  return(v)
  
}

for( i in 1:nrow(mu_anc_plt))
{
  mu_anc_plt$ANC[i]<-long_gen_a(id=mu_anc_plt$id_new[i],t_id=mu_anc_plt$time_id[i])
  mu_anc_plt$PLT[i]<-long_gen_p(id=mu_anc_plt$id_new[i],t_id=mu_anc_plt$time_id[i])
  
  i=i+1
}

mu_anc_plt$Cluster=sapply(1:nrow(mu_anc_plt), function(i) which.max(jfit$mean$pi_star[mu_anc_plt$id_new[i],]))
mu_anc_plt$ANC=mu_anc_plt$ANC+log(1564)
mu_anc_plt$PLT=mu_anc_plt$PLT+log(192000)

lfit<-mu_anc_plt%>%
  group_by(Week,Cluster)%>%
  summarise(ANC=mean(ANC),PLT=mean(PLT))

lfit$Cluster=as.factor(lfit$Cluster)
p1<-ggplot(data=lfit,aes(x=Week, y=ANC, color=Cluster))+
  geom_line(size=0.875)+ scale_color_manual(values=c("red","blue","green"))+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=12 ,color="black", face="bold"), axis.title.y =element_text(size=12 ,color="black", face="bold"))+
  xlab("Week")+ylab("log(ANC)")+
  geom_hline(yintercept = log(1784.762), col="black", linetype="dashed")+
  theme(plot.title = element_text(hjust = 0.5))


p2<-ggplot(data=lfit,aes(x=Week, y=PLT, color=Cluster))+
  geom_line(size=0.875)+ scale_color_manual(values=c("red","blue","green"))+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=12 ,color="black", face="bold"), axis.title.y =element_text(size=12 ,color="black", face="bold"))+
  xlab("Week")+ylab("log(PLT)")+
  geom_hline(yintercept = log(204074.4), col="black", linetype="dashed")+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(p1, p2,  common.legend = TRUE,legend="right")

#-----------------------------------------------------------------------------------------------#
#Survival probabilities
time_st<-c(1,seq(10,210,10))/max_week
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


data_XS=data.frame(id_new=rep(rep(1:nId,each=length(time_st)),nCluster),
                   time_id=rep(rep(1:length(time_st),nId),nCluster),
                   Week=rep(rep(time_st,nId),nCluster),
                   Cluster=rep(1:nCluster,each=(nId*(length(time_st)))),
                   Surv_prob=NA)

Surv_cl=function(id,t,Cluster)
{
  gk_cut<-(Gk_points$cuts+1)*0.5*(t*max_week -1)/max_week + 1/max_week
  d=data.frame(id_new=rep(id,length(gk_cut)),time_id=1:length(gk_cut),
               Intercept=1 ,Week =gk_cut , Week2=gk_cut^2,
               X6_mp=NA, MTX=NA)
  
  d<-merge(x=d, y=X_fix[,!colnames(X_fix)%in%c("Intercept")], by=c("id_new"), all.x=T)
  
  for(i in 1:nrow(d))
  {
    v<-NULL
    v<-medi(id=d$id_new[i], t=d$Week[i])
    d$X6_mp[i]<-v[1]
    d$MTX[i]<-v[2]
    i=i+1
  }
  i=NULL
  betas_LYM<-jfit$mean$beta_LYM[,Cluster]
  betas_ANC<-jfit$mean$beta_ANC
  betas_PLT<-jfit$mean$beta_PLT
  
  
  base_haz<-jfit$mean$gam.h0[,Cluster]
  thetas<-jfit$mean$theta[,Cluster]
  
  assoc_LYM<-jfit$mean$psi_LYM[Cluster]
  assoc_ANC<-jfit$mean$psi_ANC[Cluster]
  assoc_PLT<-jfit$mean$psi_PLT[Cluster]
  
  #Expected Longitudinal value
  long=function(i)
  {
    v_lym<- jfit$mean$u_L[id,1,Cluster] + jfit$mean$u_L[id,2,Cluster]*d$Week[i] + as.numeric(betas_LYM%*%t(d[i,3:ncol(d)]))
    v_anc<- jfit$mean$u_AP[id,1] + jfit$mean$u_AP[id,3]*d$Week[i] + as.numeric(betas_ANC%*%t(d[i,3:ncol(d)]))
    v_plt<- jfit$mean$u_AP[id,2] + jfit$mean$u_AP[id,4]*d$Week[i] + as.numeric(betas_PLT%*%t(d[i,3:ncol(d)]))
    
    v<-assoc_LYM*v_lym + assoc_ANC*v_anc + assoc_PLT*v_plt
    
    return(v)
    
  }
  
  d$mu=sapply(1:nrow(d), function(i) long(i))
  
  #Log Baseline hazard 
  log_b_haz=function(i)
  {
    v<-as.numeric(base_haz%*%t(bbase(x=d$Week[i] ,xl=quantile(Time, probs =c(0.05),names=F) , xr=quantile(Time, probs =c(0.95),names=F) , ndx=14 , deg=3)))
    return(v)
  }
  
  d$log_base_hazard=sapply(1:nrow(d), function(i) log_b_haz(i))
  
  #Baseline Covariates
  d$baseline_cov= as.numeric(thetas%*%t(d[1,8:17]))
  
  #Weights
  d$Weights=Gk_points$weights
  
  #Multiplying factor for GK approximation of Surv
  P<-(t -1/max_week)/2
  
  #Survival Probability
  net_hazard<- exp(d$mu + d$baseline_cov + d$log_base_hazard)
  
  Sv<- as.numeric(exp(-(net_hazard%*%d$Weights)*P))
  
  return(Sv)
  
}


data_XS$Surv_prob=sapply(1:nrow(data_XS),
                         function(g) Surv_cl(id=data_XS$id_new[g],t=data_XS$Week[g],Cluster=data_XS$Cluster[g]))


data_XS1<-data_XS
data_XS1$Week=data_XS1$Week*max_week

#Estimated cluster _ id combination and accordingly parameters for that cluster
data_XS1$st_cl=paste(data_XS1$id_new, data_XS1$Cluster, sep="_")
cluster_id_comb=sapply(1:nId, function(i) paste(i,which.max(jfit$mean$pi_star[i,]), sep="_"))
data_XS1<-subset(data_XS1, data_XS1$st_cl%in%cluster_id_comb)

Sfit<-data_XS1%>%
  group_by(Week,Cluster)%>%
  summarise(Surv=mean(Surv_prob))
Sfit$Cluster=as.factor(Sfit$Cluster)

#Sfit1<-Sfit[Sfit$Cluster==1,]
#Sfit2<-Sfit[Sfit$Cluster==2 & Sfit$Week<=210,]
#Sfit<-data.frame(rbind(Sfit1,Sfit2))

d0<-NULL
for(cl in 1:nCluster)
{
  d_sub<-NULL
  d_sub=subset(Sfit,Sfit$Cluster==cl)
  
  lo<-loess(Surv~Week, data=d_sub)
  d_sub$Surv_new=predict(lo,d_sub$Week)
  d1<-d_sub
  
  d0<-rbind(d0,d1)
}

d0<-data.frame(d0)
rownames(d0)<-NULL
d0$Cluster=as.factor(d0$Cluster)

#Loess smoothing mass make some survival values in the initial weeks to be greater than 1 or can 
#make survival at time =0 to be <1
d0$Surv_new=ifelse(d0$Surv_new>1,1,d0$Surv_new)
d0[d0$Week==1,"Surv_new"]=1
#Smoothed survival
ggplot(data=d0,aes(x=Week, y=Surv_new, color=Cluster))+
  geom_line(size=0.875)+ scale_color_manual(values=c("red","blue","green"))+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=12 ,color="black", face="bold"), axis.title.y =element_text(size=12 ,color="black", face="bold"))+
  xlab("Week")+ylab("1-P(relapse)")+ ylim(c(0,1))+
  theme(plot.title = element_text(hjust = 0.5))

#-------------------------------------------------------------------------------------------------
#Polarity of Ids based on pi*

pi_stars<-data.frame(id_new=1:nId)
pi_stars$Cluster_1<-jfit$mean$pi_star[,1]
pi_stars$Cluster_2<-jfit$mean$pi_star[,2]
pi_stars$Dom=ifelse(pi_stars$Cluster_1>pi_stars$Cluster_2,pi_stars$Cluster_1,pi_stars$Cluster_2)
pi_stars$Cluster<-ifelse(pi_stars$Cluster_1>pi_stars$Cluster_2,"1","2")

ggplot(data=pi_stars,aes(x=id_new, y=Cluster_1, group=Cluster))+
  geom_point(aes(color=Cluster),size=1.5)+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9))+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=12 ,color="black", face="bold"), axis.title.y =element_text(size=12 ,color="black", face="bold"))+
  xlab("Week")+ylab("1-P(relapse)")+
  xlab("Id")+ylab(TeX("$P(G_i=1|rest)$"))+
  geom_hline(yintercept = c(0.45,0.5,0.55), col="black", linetype=c("dotted","dashed","dotted"))+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")+
  labs(color=NULL)+
  scale_color_manual(
    values=c("red","blue"),
    breaks=c("1","2"),
    labels=c(paste("Cluster_1 ","(",round(nrow(pi_stars[pi_stars$Cluster=="1",])*100/nrow(pi_stars),2),"%)", sep=""),
             paste("Cluster_2 ","(",round(nrow(pi_stars[pi_stars$Cluster=="2",])*100/nrow(pi_stars),2),"%)", sep=""))
    
  )

#------------------------------------------------------------------------------------------------------------------------#
#Kaplan Meier plot cluster-specific

library(survminer)
y_s1<-y_s
y_s1$Cluster<-sapply(1:nrow(y_s), function(i) which.max(jfit$mean$pi_star[y_s$id_new[i],]))
y_s1$Time_to_event=y_s1$Time_to_event*max_week
model_fit <- survfit(Surv(Time_to_event, Relapsed) ~ Cluster, data = y_s1)
ggsurvplot(model_fit, data = y_s1,
           ylab=c("1-P(relapse)"),
           xlab=c("Week"),
           break.time.by=c(50),
           xlim=c(5,205),
           # Change legends: title & labels
           legend.title = "Cluster",
           legend.labs = c("1", "2"),
           legend=c("right"),
           palette =c("red","blue"),
           censor.shape="|",
           ggtheme=theme_bw()+
             theme(axis.text.x = element_text(face="bold", 
                                              size=9))+
             theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
                   axis.title.x =element_text(size=12 ,color="black", face="bold"), axis.title.y =element_text(size=12 ,color="black", face="bold"))
)

#---------------------------------------------------------------------------------------------------------------------------#
#Confidence Intervals

Basic_df<-data.frame(Var=colnames(X)[c(-1,-2)],lci =NA,mean=NA, uci=NA, Response =NA, Cluster=NA, Parameter_type=NA)

#Beta
df_lym1<-Basic_df
df_lym1$Cluster=1
df_lym1$Response<-"LYM"
df_lym1$lci=jfit$q2.5$beta_LYM[,1]
df_lym1$mean=jfit$mean$beta_LYM[,1]
df_lym1$uci=jfit$q97.5$beta_LYM[,1]
df_lym1$Parameter_type="beta"


df_lym2<-Basic_df
df_lym2$Cluster=2
df_lym2$Response<-"LYM"
df_lym2$lci=jfit$q2.5$beta_LYM[,2]
df_lym2$mean=jfit$mean$beta_LYM[,2]
df_lym2$uci=jfit$q97.5$beta_LYM[,2]
df_lym2$Parameter_type="beta"

df_anc<-Basic_df
df_anc$Cluster=NA
df_anc$Response<-"ANC"
df_anc$lci=jfit$q2.5$beta_ANC
df_anc$mean=jfit$mean$beta_ANC
df_anc$uci=jfit$q97.5$beta_ANC
df_anc$Parameter_type="beta"

df_plt<-Basic_df
df_plt$Cluster=NA
df_plt$Response<-"PLT"
df_plt$lci=jfit$q2.5$beta_PLT
df_plt$mean=jfit$mean$beta_PLT
df_plt$uci=jfit$q97.5$beta_PLT
df_plt$Parameter_type="beta"


#theta
Basic_df1<-data.frame(Var=colnames(X)[c(-1:-7)],lci =NA,mean=NA, uci=NA, Response =NA, Cluster=NA, Parameter_type=NA)


df_theta1<-Basic_df1
df_theta1$Cluster=1
df_theta1$Response<-NA
df_theta1$lci=jfit$q2.5$theta[,1]
df_theta1$mean=jfit$mean$theta[,1]
df_theta1$uci=jfit$q97.5$theta[,1]
df_theta1$Parameter_type="theta"

df_theta2<-Basic_df1
df_theta2$Cluster=2
df_theta2$Response<-NA
df_theta2$lci=jfit$q2.5$theta[,2]
df_theta2$mean=jfit$mean$theta[,2]
df_theta2$uci=jfit$q97.5$theta[,2]
df_theta2$Parameter_type="theta"


#Association
Basic_df2<-data.frame(Var=NA,lci =NA,mean=NA, uci=NA, Response =NA, Cluster=1:nCluster, Parameter_type=NA)
df_lym_asso<-Basic_df2
df_lym_asso$Var="LYM"
df_lym_asso$Response="LYM"
df_lym_asso$lci=jfit$q2.5$psi_LYM
df_lym_asso$mean=jfit$mean$psi_LYM
df_lym_asso$uci=jfit$q97.5$psi_LYM
df_lym_asso$Parameter_type="Association"

df_anc_asso<-Basic_df2
df_anc_asso$Var="ANC"
df_anc_asso$Response="ANC"
df_anc_asso$lci=jfit$q2.5$psi_ANC
df_anc_asso$mean=jfit$mean$psi_ANC
df_anc_asso$uci=jfit$q97.5$psi_ANC
df_anc_asso$Parameter_type="Association"

df_plt_asso<-Basic_df2
df_plt_asso$Var="PLT"
df_plt_asso$Response="PLT"
df_plt_asso$lci=jfit$q2.5$psi_PLT
df_plt_asso$mean=jfit$mean$psi_PLT
df_plt_asso$uci=jfit$q97.5$psi_PLT
df_plt_asso$Parameter_type="Association"

df<-data.frame(rbind(df_lym1,df_lym2,df_anc, df_plt,df_theta1, df_theta2, df_lym_asso, df_anc_asso,df_plt_asso))





#Association
df_asso<-df[df$Parameter_type=="Association",]
df_asso$x<-rep(c(2,4,6),each=2)+rep(c(-0.25,0.25),3)
p_asso<- ggplot(df_asso, aes(x=mean, y=x ,color=as.factor(Cluster) )) + 
  scale_color_manual(name="Cluster",values=c("red","blue"))+
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(xmin=lci, xmax=uci), width=0.25, size=1,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9) ,legend.position="bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=12 ,color="black", face="bold"),
        axis.title.x =element_text(size=12 ,color="black", face="bold"), axis.title.y =element_text(size=12 ,color="black", face="bold"))+
  ylab("Response")+xlab("Association Estimates")+
  scale_y_continuous(labels=c("LYM","ANC","PLT"),
                     breaks=c(2,4,6))+
  geom_vline(xintercept = 0, linetype="dashed")
p_asso




#Baseline cov
df_base<-df[df$Parameter_type=="theta",]
df_base$x<-rep(c(1:(nrow(df_base)/2))*2,2)+rep(c(-0.25,0.25),each=(nrow(df_base)/2))
p_base<- ggplot(df_base, aes(x=mean, y=x ,color=as.factor(Cluster) )) + 
  scale_color_manual(name="Cluster",values=c("red","blue"))+
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(xmin=lci, xmax=uci), width=0.5, size=1,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9) ,legend.position="bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("Baseline Covariates")+xlab("Baseline Coefficients")+
  scale_y_continuous(labels=df_base$Var[1:(nrow(df_base)/2)],
                     breaks=c(1:(nrow(df_base)/2))*2)+
  geom_vline(xintercept = 0, linetype="dashed")
p_base

#Beta_LYM
df_lym<-df[df$Parameter_type=="beta"& df$Response=="LYM"& !df$Var%in%c("MTX","X6_mp","Week2","Week","Intercept"),]
df_lym$x<-rep(c(1:(nrow(df_lym)/2))*2,2)+rep(c(-0.25,0.25),each=(nrow(df_lym)/2))
p_lym<- ggplot(df_lym, aes(x=mean, y=x ,color=as.factor(Cluster) )) + 
  scale_color_manual(name="Cluster",values=c("red","blue"))+
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(xmin=lci, xmax=uci), width=0.5, size=1,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9) ,legend.position="bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("Fixed Covariates")+xlab("Beta_LYM")+
  scale_y_continuous(labels=df_lym$Var[1:(nrow(df_lym)/2)],
                     breaks=c(1:(nrow(df_lym)/2))*2)+
  geom_vline(xintercept = 0, linetype="dashed")
p_lym

df_lym_med<-df[df$Parameter_type=="beta"& df$Response=="LYM"& df$Var%in%c("MTX","X6_mp"),]
df_lym_med$x<-rep(c(1:(nrow(df_lym_med)/2))*2,2)+rep(c(-0.25,0.25),each=(nrow(df_lym_med)/2))
p_lym_med<- ggplot(df_lym_med, aes(x=mean, y=x ,color=as.factor(Cluster) )) + 
  scale_color_manual(name="Cluster",values=c("red","blue"))+
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(xmin=lci, xmax=uci), width=0.25, size=1,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9) ,legend.position="bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("Medicine")+xlab("Beta_LYM")+
  scale_y_continuous(labels=c("6MP","MTX"),
                     breaks=c(1:(nrow(df_lym_med)/2))*2)+
  geom_vline(xintercept = 0, linetype="dashed")
p_lym_med


#Beta_ANC, PLT
df_AP<-df[df$Parameter_type=="beta"& df$Response!="LYM"& !df$Var%in%c("MTX","X6_mp","Week2","Week","Intercept"),]
df_AP$x<-rep(c(1:(nrow(df_AP)/2))*2,2)+rep(c(-0.25,0.25),each=(nrow(df_AP)/2))
p_AP<- ggplot(df_AP, aes(x=mean, y=x ,color=as.factor(Response) )) + 
  scale_color_manual(name="Response",values=c("#DF536B","#61D04F"))+
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(xmin=lci, xmax=uci), width=0.5, size=1,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9) ,legend.position="bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("Fixed Covariates")+xlab("Beta")+
  scale_y_continuous(labels=df_AP$Var[1:(nrow(df_AP)/2)],
                     breaks=c(1:(nrow(df_AP)/2))*2)+
  geom_vline(xintercept = 0, linetype="dashed")
p_AP

df_AP_med<-df[df$Parameter_type=="beta"& df$Response!="LYM"& df$Var%in%c("MTX","X6_mp"),]
df_AP_med$x<-rep(c(1:(nrow(df_AP_med)/2))*2,2)+rep(c(-0.25,0.25),each=(nrow(df_AP_med)/2))
p_AP_med<- ggplot(df_AP_med, aes(x=mean, y=x ,color=as.factor(Response) )) + 
  scale_color_manual(name="Response",values=c("#DF536B","#61D04F"))+
  geom_point(size=2)+
  theme_bw()+
  geom_errorbar(aes(xmin=lci, xmax=uci), width=0.25, size=1,
                position=position_dodge(0.05))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=9) ,legend.position="bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(size=10 ,color="black", face="bold"),axis.text.y=element_text(size=10 ,color="black", face="bold"),
        axis.title.x =element_text(size=10 ,color="black", face="bold"), axis.title.y =element_text(size=10 ,color="black", face="bold"))+
  ylab("Medicine")+xlab("Beta")+
  scale_y_continuous(labels=c("6MP","MTX"),
                     breaks=c(1:(nrow(df_AP_med)/2))*2)+
  geom_vline(xintercept = 0, linetype="dashed")
p_AP_med


#Convergence plots
jagsfits<- readRDS(paste(path_out,paste("/JM_comp_Lymph_Cl_",nCluster,".rds",sep=""),sep=""))


#Association Convergence
response_name=c("LYM","ANC","PLT")
association=data.frame(Iterations=rep(1:600,nCluster*length(response_name)),
                       Response=rep(rep(response_name,each=600),nCluster),
                       Cluster=rep(1:nCluster,each=600*length(response_name)),
                       value=NA,
                       cumulative_value=NA)

for(i in 1:nCluster)
{
  
  association[association$Cluster==i & association$Response==response_name[1],"value"]=jagsfits$sims.list$psi_LYM[,i]
  association[association$Cluster==i & association$Response==response_name[2],"value"]=jagsfits$sims.list$psi_ANC[,i]
  association[association$Cluster==i & association$Response==response_name[3],"value"]=jagsfits$sims.list$psi_PLT[,i]
  association[association$Cluster==i & association$Response==response_name[1],"cumulative_value"]=cummean(jagsfits$sims.list$psi_LYM[,i])
  association[association$Cluster==i & association$Response==response_name[2],"cumulative_value"]=cummean(jagsfits$sims.list$psi_ANC[,i])
  association[association$Cluster==i & association$Response==response_name[3],"cumulative_value"]=cummean(jagsfits$sims.list$psi_PLT[,i])
  
}


p1<-ggplot(data=association[association$Cluster==1,],aes(x=value,color=Response))+
  geom_density(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Density")+
  geom_vline(xintercept = c(jagsfits$mean$psi_LYM[1],jagsfits$mean$psi_ANC[1],jagsfits$mean$psi_PLT[1])
             , color=c("blue","red","green"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p1<-p1+theme(legend.position="none")
p2<-ggplot(data=association[association$Cluster==1,],aes(x=Iterations,y=value, color=Response))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("Iterations")+ylab("Trace")+
  geom_hline(yintercept = c(jagsfits$mean$psi_LYM[1],jagsfits$mean$psi_ANC[1],jagsfits$mean$psi_PLT[1]),
              color=c("blue","red","green"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p2<-p2+theme(legend.position="none")
p3<-ggplot(data=association[association$Cluster==1,],aes(x=Iterations,y=cumulative_value, color=Response))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("Iterations")+ylab("Cumulative Trace")+
  geom_hline(yintercept = c(jagsfits$mean$psi_LYM[1],jagsfits$mean$psi_ANC[1],jagsfits$mean$psi_PLT[1]),
              color=c("blue","red","green"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))

p4<-grid.arrange(p1,p2,p3, ncol=3, nrow=1,widths=c(1,1,1.25),top="Cluster 1") 



p5<-ggplot(data=association[association$Cluster==2,],aes(x=value,color=Response))+
  geom_density(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Density")+
  geom_vline(xintercept = c(jagsfits$mean$psi_LYM[2],jagsfits$mean$psi_ANC[2],jagsfits$mean$psi_PLT[2]),
              color=c("blue","red","green"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p5<-p5+theme(legend.position="none")
p6<-ggplot(data=association[association$Cluster==2,],aes(x=Iterations,y=value, color=Response))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("Iterations")+ylab("Trace")+
  geom_hline(yintercept = c(jagsfits$mean$psi_LYM[2],jagsfits$mean$psi_ANC[2],jagsfits$mean$psi_PLT[2]),
              color=c("blue","red","green"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p6<-p6+theme(legend.position="none")
p7<-ggplot(data=association[association$Cluster==2,],aes(x=Iterations,y=cumulative_value, color=Response))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("Iterations")+ylab("Cumulative Trace")+
  geom_hline(yintercept = c(jagsfits$mean$psi_LYM[2],jagsfits$mean$psi_ANC[2],jagsfits$mean$psi_PLT[2]),
              color=c("blue","red","green"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))

p8<-grid.arrange(p5,p6,p7, ncol=3, nrow=1,widths=c(1,1,1.25),top="Cluster 2") 


grid.arrange(p4,p8, ncol=1, nrow=2)


#Beta Convergence
response_name=c("LYM","ANC","PLT")
beta_4=data.frame(Iterations=rep(1:600,(nCluster+length(response_name)-1)),
                  Response=rep(response_name[c(1,1,2,3)],each=600),
                  Cluster=rep(c(1,2,NA,NA),each=600),
                  value=NA,
                  cumulative_value=NA)
beta_7=beta_4

beta_4[beta_4$Response=="ANC","value"]=jagsfits$sims.list$beta_ANC[,4]
beta_4[beta_4$Response=="ANC","cumulative_value"]=cummean(jagsfits$sims.list$beta_ANC[,4])

beta_4[beta_4$Response=="PLT","value"]=jagsfits$sims.list$beta_PLT[,4]
beta_4[beta_4$Response=="PLT","cumulative_value"]=cummean(jagsfits$sims.list$beta_PLT[,4])

beta_4[beta_4$Response=="LYM" & beta_4$Cluster==1,"value"]=jagsfits$sims.list$beta_LYM[,4,1]
beta_4[beta_4$Response=="LYM"& beta_4$Cluster==1,"cumulative_value"]=cummean(jagsfits$sims.list$beta_LYM[,4,1])

beta_4[beta_4$Response=="LYM" & beta_4$Cluster==2,"value"]=jagsfits$sims.list$beta_LYM[,4,2]
beta_4[beta_4$Response=="LYM"& beta_4$Cluster==2,"cumulative_value"]=cummean(jagsfits$sims.list$beta_LYM[,4,2])

beta_4$Cov_name="6MP"
beta_7[beta_7$Response=="ANC","value"]=jagsfits$sims.list$beta_ANC[,7]
beta_7[beta_7$Response=="ANC","cumulative_value"]=cummean(jagsfits$sims.list$beta_ANC[,7])

beta_7[beta_7$Response=="PLT","value"]=jagsfits$sims.list$beta_PLT[,7]
beta_7[beta_7$Response=="PLT","cumulative_value"]=cummean(jagsfits$sims.list$beta_PLT[,7])

beta_7[beta_7$Response=="LYM" & beta_7$Cluster==1,"value"]=jagsfits$sims.list$beta_LYM[,7,1]
beta_7[beta_7$Response=="LYM"& beta_7$Cluster==1,"cumulative_value"]=cummean(jagsfits$sims.list$beta_LYM[,7,1])

beta_7[beta_7$Response=="LYM" & beta_7$Cluster==2,"value"]=jagsfits$sims.list$beta_LYM[,7,2]
beta_7[beta_7$Response=="LYM"& beta_7$Cluster==2,"cumulative_value"]=cummean(jagsfits$sims.list$beta_LYM[,7,2])
beta_7$Cov_name="Gender"

beta<-data.frame(rbind(beta_4,beta_7))
rownames(beta)<-NULL


#Beta 6MP
p1<-ggplot(data=beta[is.na(beta$Cluster) &beta$Cov_name=="6MP",],aes(x=value,color=Response))+
  geom_density(size=1.1)+theme_bw()+scale_color_manual(values=c("#DF536B","#61D04F"))+
  xlab("")+ylab("Density")+
  geom_vline(xintercept = c(jagsfits$mean$beta_ANC[4],jagsfits$mean$beta_PLT[4]),
              color=c("#DF536B","#61D04F"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p1<-p1+theme(legend.position="none")
p2<-ggplot(data=beta[is.na(beta$Cluster) &beta$Cov_name=="6MP",],aes(x=Iterations,y=value,color=Response))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("#DF536B","#61D04F"))+
  xlab("")+ylab("Trace")+
  geom_hline(yintercept = c(jagsfits$mean$beta_ANC[4],jagsfits$mean$beta_PLT[4]),
              color=c("#DF536B","#61D04F"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p2<-p2+theme(legend.position="none")
p3<-ggplot(data=beta[is.na(beta$Cluster) &beta$Cov_name=="6MP",],aes(x=Iterations,y=cumulative_value,color=Response))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("#DF536B","#61D04F"))+
  xlab("")+ylab("Cumulative Trace")+
  geom_hline(yintercept = c(jagsfits$mean$beta_ANC[4],jagsfits$mean$beta_PLT[4]),
              color=c("#DF536B","#61D04F"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))

p4<-grid.arrange(p1,p2,p3, ncol=3, nrow=1,widths=c(1,1,1.25),top="Effect of 6MP on ANC, PLT") 

beta$Cluster=as.factor(beta$Cluster)
p5<-ggplot(data=beta[!is.na(beta$Cluster) &beta$Cov_name=="6MP",],aes(x=value,color=Cluster))+
  geom_density(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Density")+
  geom_vline(xintercept = c(jagsfits$mean$beta_LYM[4,1],jagsfits$mean$beta_LYM[4,2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p5<-p5+theme(legend.position="none")
p6<-ggplot(data=beta[!is.na(beta$Cluster) &beta$Cov_name=="6MP",],aes(x=Iterations,y=value,color=Cluster))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Trace")+
  geom_hline(yintercept = c(jagsfits$mean$beta_LYM[4,1],jagsfits$mean$beta_LYM[4,2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p6<-p6+theme(legend.position="none")
p7<-ggplot(data=beta[!is.na(beta$Cluster) &beta$Cov_name=="6MP",],aes(x=Iterations,y=cumulative_value,color=Cluster))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Cumulative Trace")+
  geom_hline(yintercept = c(jagsfits$mean$beta_LYM[4,1],jagsfits$mean$beta_LYM[4,2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))

p8<-grid.arrange(p5,p6,p7, ncol=3, nrow=1,widths=c(1,1,1.25),top="Effect of 6MP on LYM") 
grid.arrange(p8,p4, ncol=1, nrow=2)



#Beta Gender
p1<-ggplot(data=beta[is.na(beta$Cluster) &beta$Cov_name=="Gender",],aes(x=value,color=Response))+
  geom_density(size=1.1)+theme_bw()+scale_color_manual(values=c("#DF536B","#61D04F"))+
  xlab("")+ylab("Density")+
  geom_vline(xintercept = c(jagsfits$mean$beta_ANC[7],jagsfits$mean$beta_PLT[7]),
              color=c("#DF536B","#61D04F"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p1<-p1+theme(legend.position="none")
p2<-ggplot(data=beta[is.na(beta$Cluster) &beta$Cov_name=="Gender",],aes(x=Iterations,y=value,color=Response))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("#DF536B","#61D04F"))+
  xlab("")+ylab("Trace")+
  geom_hline(yintercept = c(jagsfits$mean$beta_ANC[7],jagsfits$mean$beta_PLT[7]),
              color=c("#DF536B","#61D04F"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p2<-p2+theme(legend.position="none")
p3<-ggplot(data=beta[is.na(beta$Cluster) &beta$Cov_name=="Gender",],aes(x=Iterations,y=cumulative_value,color=Response))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("#DF536B","#61D04F"))+
  xlab("")+ylab("Cumulative Trace")+
  geom_hline(yintercept = c(jagsfits$mean$beta_ANC[7],jagsfits$mean$beta_PLT[7]),
              color=c("#DF536B","#61D04F"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))

p4<-grid.arrange(p1,p2,p3, ncol=3, nrow=1,widths=c(1,1,1.25),top="Effect of Gender on ANC, PLT") 

beta$Cluster=as.factor(beta$Cluster)
p5<-ggplot(data=beta[!is.na(beta$Cluster) &beta$Cov_name=="Gender",],aes(x=value,color=Cluster))+
  geom_density(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Density")+
  geom_vline(xintercept = c(jagsfits$mean$beta_LYM[7,1],jagsfits$mean$beta_LYM[7,2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p5<-p5+theme(legend.position="none")
p6<-ggplot(data=beta[!is.na(beta$Cluster) &beta$Cov_name=="Gender",],aes(x=Iterations,y=value,color=Cluster))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Trace")+
  geom_hline(yintercept = c(jagsfits$mean$beta_LYM[7,1],jagsfits$mean$beta_LYM[7,2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p6<-p6+theme(legend.position="none")
p7<-ggplot(data=beta[!is.na(beta$Cluster) &beta$Cov_name=="Gender",],aes(x=Iterations,y=cumulative_value,color=Cluster))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Cumulative Trace")+
  geom_hline(yintercept = c(jagsfits$mean$beta_LYM[7,1],jagsfits$mean$beta_LYM[7,2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))

p8<-grid.arrange(p5,p6,p7, ncol=3, nrow=1,widths=c(1,1,1.25),top="Effect of Gender on LYM") 
grid.arrange(p8,p4, ncol=1, nrow=2)



#Theta Convergence
Covariate_name=c("Age","WBC_high")
theta=data.frame(Iterations=rep(1:600,(nCluster*length(Covariate_name))),
                 Covariate=rep(rep(Covariate_name,each=nCluster),each=600),
                 Cluster=rep(rep(c(1,2),2),each=600),
                 value=NA,
                 cumulative_value=NA)


for(cl in 1:nCluster)
{
  theta[theta$Covariate=="Age" & theta$Cluster==cl,"value"]=jagsfits$sims.list$theta[,1,cl]
  theta[theta$Covariate=="Age"& theta$Cluster==cl,"cumulative_value"]=cummean(jagsfits$sims.list$theta[,1,cl])
  
  theta[theta$Covariate=="WBC_high" & theta$Cluster==cl,"value"]=jagsfits$sims.list$theta[,4,cl]
  theta[theta$Covariate=="WBC_high"& theta$Cluster==cl,"cumulative_value"]=cummean(jagsfits$sims.list$theta[,4,cl])
  
  
}


theta$Cluster=as.factor(theta$Cluster)

p1<-ggplot(data=theta[theta$Covariate=="Age",],aes(x=value,color=Cluster))+
  geom_density(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Density")+
  geom_vline(xintercept = c(jagsfits$mean$theta[1,1],jagsfits$mean$theta[1,2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p1<-p1+theme(legend.position="none")
p2<-ggplot(data=theta[theta$Covariate=="Age",],aes(x=Iterations,y=value,color=Cluster))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Trace")+
  geom_hline(yintercept = c(jagsfits$mean$theta[1,1],jagsfits$mean$theta[1,2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p2<-p2+theme(legend.position="none")
p3<-ggplot(data=theta[theta$Covariate=="Age",],aes(x=Iterations,y=cumulative_value,color=Cluster))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Cumulative Trace")+
  geom_hline(yintercept = c(jagsfits$mean$theta[1,1],jagsfits$mean$theta[1,2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))

p4<-grid.arrange(p1,p2,p3, ncol=3, nrow=1,widths=c(1,1,1.25),top="Baseline Coefficient of Age") 

p1<-ggplot(data=theta[theta$Covariate=="WBC_high",],aes(x=value,color=Cluster))+
  geom_density(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Density")+
  geom_vline(xintercept = c(jagsfits$mean$theta[4,1],jagsfits$mean$theta[4,2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p1<-p1+theme(legend.position="none")
p2<-ggplot(data=theta[theta$Covariate=="WBC_high",],aes(x=Iterations,y=value,color=Cluster))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Trace")+
  geom_hline(yintercept = c(jagsfits$mean$theta[4,1],jagsfits$mean$theta[4,2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p2<-p2+theme(legend.position="none")
p3<-ggplot(data=theta[theta$Covariate=="WBC_high",],aes(x=Iterations,y=cumulative_value,color=Cluster))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Cumulative Trace")+
  geom_hline(yintercept = c(jagsfits$mean$theta[4,1],jagsfits$mean$theta[4,2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))

p4<-grid.arrange(p1,p2,p3, ncol=3, nrow=1,widths=c(1,1,1.25),top="Baseline Coefficient of WBC_high") 



#Error variance
#sigma Convergence
response_name=c("LYM","ANC","PLT")
sigma=data.frame(Iterations=rep(1:600,(nCluster+length(response_name)-1)),
                 Response=rep(response_name[c(1,1,2,3)],each=600),
                 Cluster=rep(c(1,2,NA,NA),each=600),
                 value=NA,
                 cumulative_value=NA)

sigma[sigma$Response=="ANC","value"]=jagsfits$sims.list$sig_ANC
sigma[sigma$Response=="ANC","cumulative_value"]=cummean(jagsfits$sims.list$sig_ANC)

sigma[sigma$Response=="PLT","value"]=jagsfits$sims.list$sig_PLT
sigma[sigma$Response=="PLT","cumulative_value"]=cummean(jagsfits$sims.list$sig_PLT)

sigma[sigma$Response=="LYM" & sigma$Cluster==1,"value"]=jagsfits$sims.list$sig_LYM[,1]
sigma[sigma$Response=="LYM"& sigma$Cluster==1,"cumulative_value"]=cummean(jagsfits$sims.list$sig_LYM[,1])

sigma[sigma$Response=="LYM" & sigma$Cluster==2,"value"]=jagsfits$sims.list$sig_LYM[,2]
sigma[sigma$Response=="LYM"& sigma$Cluster==2,"cumulative_value"]=cummean(jagsfits$sims.list$sig_LYM[,2])

#Error variance
p1<-ggplot(data=sigma[is.na(sigma$Cluster) ,],aes(x=value,color=Response))+
  geom_density(size=1.1)+theme_bw()+scale_color_manual(values=c("#DF536B","#61D04F"))+
  xlab("")+ylab("Density")+
  geom_vline(xintercept = c(jagsfits$mean$sig_ANC,jagsfits$mean$sig_PLT),
              color=c("#DF536B","#61D04F"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p1<-p1+theme(legend.position="none")
p2<-ggplot(data=sigma[is.na(sigma$Cluster) ,],aes(x=Iterations,y=value,color=Response))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("#DF536B","#61D04F"))+
  xlab("")+ylab("Trace")+
  geom_hline(yintercept = c(jagsfits$mean$sig_ANC,jagsfits$mean$sig_PLT),
              color=c("#DF536B","#61D04F"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p2<-p2+theme(legend.position="none")
p3<-ggplot(data=sigma[is.na(sigma$Cluster) ,],aes(x=Iterations,y=cumulative_value,color=Response))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("#DF536B","#61D04F"))+
  xlab("")+ylab("Cumulative Trace")+
  geom_hline(yintercept = c(jagsfits$mean$sig_ANC,jagsfits$mean$sig_PLT),
              color=c("#DF536B","#61D04F"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))

p4<-grid.arrange(p1,p2,p3, ncol=3, nrow=1,widths=c(1,1,1.25),top="Error variance of ANC, PLT") 

sigma$Cluster=as.factor(sigma$Cluster)

p5<-ggplot(data=sigma[!is.na(sigma$Cluster) &sigma$Response=="LYM",],aes(x=value,color=Cluster))+
  geom_density(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Density")+
  geom_vline(xintercept = c(jagsfits$mean$sig_LYM[1],jagsfits$mean$sig_LYM[2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p5<-p5+theme(legend.position="none")
p6<-ggplot(data=sigma[!is.na(sigma$Cluster) &sigma$Response=="LYM",],aes(x=Iterations,y=value,color=Cluster))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Trace")+
  geom_hline(yintercept = c(jagsfits$mean$sig_LYM[1],jagsfits$mean$sig_LYM[2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))
p6<-p6+theme(legend.position="none")
p7<-ggplot(data=sigma[!is.na(sigma$Cluster) &sigma$Response=="LYM",],aes(x=Iterations,y=cumulative_value,color=Cluster))+
  geom_line(size=1.1)+theme_bw()+scale_color_manual(values=c("red","blue","green"))+
  xlab("")+ylab("Cumulative Trace")+
  geom_hline(yintercept = c(jagsfits$mean$sig_LYM[1],jagsfits$mean$sig_LYM[2]),
              color=c("red","blue"), size=1.1)+
  theme(plot.title = element_text(hjust = 0.5))

p8<-grid.arrange(p5,p6,p7, ncol=3, nrow=1,widths=c(1,1,1.25),top="Error variance of LYM") 
grid.arrange(p8,p4, ncol=1, nrow=2)



