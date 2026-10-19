#File Path
fp="C:/Users/Damitri/Desktop/Experiment JM/EXP_JAGS"

#Input path
path=paste(fp,"/Cluster_JM/Exp_work/DATA",sep="")
#Output path
path_out=paste(fp,"/JAGS_Output/Cluster_JM",sep="")
dir.create(path_out)


#Complete Data
data_org=readRDS(paste(path,"/data_0.rds",sep=""))

#The response
data_org$Lymph=exp(data_org$WBC+log(3200))-exp(data_org$ANC+log(1564))
data_org$ANC=exp(data_org$ANC+log(1564))
data_org$Platelet=exp(data_org$Platelet+log(192000))

#The medicine
data_org$X6_mp=data_org$X6_mp*275
data_org$MTX=data_org$MTX*12.5

#The covariates
data_org$WBC_high=exp(data_org$WBC_high)*21250
data_org$Age=data_org$Age*4.280008 + 6.589453


#The time
data_org$Week=data_org$Week*277
data_org$Time_to_event=data_org$Time_to_event*277
data_org$Last_follow_up_date=data_org$Last_follow_up_date*277



#Subset
data_org=data_org[,c("Child.Id","id_new","time_id","Week","Lymph","ANC","Platelet",
                     "X6_mp","MTX","Age","Gender","Lineage","WBC_high","NCI","Bulky","CNS","Risk","Risk_present","Morph_CR1",
                     "Relapsed","Last_follow_up_date","Time_to_event")]


data_org=subset(data_org, !is.na(data_org$Lymph))



#Keeping id's with > 5 time points  NA and updating the id names
id_set=unique(data_org$Child.Id)
dat=NULL
for(i in 1:length(id_set))
{
  d=subset(data_org, data_org$Child.Id==id_set[i])
  dat=rbind(dat, c(id_set[i],nrow(d)))
  i=i+1
}
dat=data.frame(dat)
dat$X2=as.numeric(dat$X2)
dat1=subset(dat, dat$X2>=5)
data_org1=subset(data_org, data_org$Child.Id%in%unique(dat1$X1))
unique(data_org1$id_new)
dat1=data_org1[,c("Child.Id","id_new")]
dat1=dat1[!duplicated(dat1), ]
dat1$id_new=sapply(dat1$Child.Id,function(x)as.numeric(strsplit(x,"_")[[1]][2]))
dat1=dat1[with(dat1, order(id_new)), ]
dat1$id_new=1:nrow(dat1)



data_org=merge(x=data_org1[,!colnames(data_org1)%in%c("id_new")], y=dat1, by=c("Child.Id"), all.x=T)
data_org=data_org[,c("Child.Id","id_new","time_id","Week","Lymph","ANC","Platelet",
                     "X6_mp","MTX","Age","Gender","Lineage","WBC_high","NCI","Bulky","CNS","Risk","Risk_present","Morph_CR1",
                     "Relapsed","Last_follow_up_date","Time_to_event")]

data_org=data_org[
  with(data_org, order(id_new, time_id)),
]

min_w=data_org%>%
  group_by(Child.Id)%>%
  summarise(min_week=min(Week))
data_org=merge(x=data_org, y=min_w,by=c("Child.Id"),all.x=T)
data_org$Week=data_org$Week-data_org$min_week+1
data_org$Time_to_event=data_org$Time_to_event-data_org$min_week+1
data_org$Last_follow_up_date=data_org$Last_follow_up_date-data_org$min_week+1



data_org=data_org[,c("Child.Id","id_new","time_id","Week","Lymph","ANC","Platelet",
                     "X6_mp","MTX","Age","Gender","Lineage","WBC_high","NCI","Bulky","CNS","Risk","Risk_present","Morph_CR1",
                     "Relapsed","Last_follow_up_date","Time_to_event")]


#Male :1 , Female :0
data_org$Gender=ifelse(data_org$Gender==1,"M","F")
#B :1 , T :2 , Other :0
data_org$Lineage=ifelse(data_org$Lineage==1,"B",ifelse(data_org$Lineage==2,"T","Other"))
#HR :1 , SR :2 , Un :0
data_org$NCI=ifelse(data_org$NCI==1,"HR",
                  ifelse(data_org$NCI==0,"Un","SR"))
#Y :1 , N :2 , Un :0
data_org$Bulky=ifelse(data_org$Bulky==1,"Y",
                    ifelse(data_org$Bulky==0,"Un","N"))
#Y :1 , N :2 , Un :0
data_org$CNS=ifelse(data_org$CNS==1,"Y",
                  ifelse(data_org$CNS==0,"Un","N"))
#Good :1 , Poor :2 , T cell :3 , Other :4 , Un :0
data_org$Risk=ifelse(data_org$Risk==1,"Good",ifelse(data_org$Risk==2,"Poor",
                                                ifelse(data_org$Risk==3,"T cell",
                                                       ifelse(data_org$Risk==0,"Un","Other"))))
#HR :1 , SR :2 , IR :3 , Un :0
data_org$Risk_present=ifelse(data_org$Risk_present==1,"HR",
                           ifelse(data_org$Risk_present==2,"SR",
                                  ifelse(data_org$Risk_present==0,"Un","IR")))


#Y :1 , N :2 , Un :0
data_org$Morph_CR1=ifelse(data_org$Morph_CR1==1,"Y",
                        ifelse(data_org$Morph_CR1==0,"Un","N"))

data_org=data_org[,!colnames(data_org)%in%c("time_id")]
write.csv(data_org, paste(path_out,paste("/data_paper3",".csv",sep=""),sep=""))
