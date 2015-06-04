
rm(list=ls(all=TRUE)) 
graphics.off()
#windows(record=T)
#setwd ("\\\\CCFHR-S-1534090\\popdyn1\\Purcell\\RedSnapIndex_SEDAR41")
library(MASS) 
library(doBy)
library(statmod)
library(Hmisc)
library(pscl)
library(lmtest)
library(ggplot2)
library(xtable)
library(stargazer)


# Read in Red Snapper Data
raw=read.csv("C:\\Users\\Kevin.Purcell\\Documents\\GitHub\\SEDAR41_RS_VI\\data\\redsnapper_SEDAR41_2015v2.csv")
names(raw)
str(raw)


# Full data subsetting: 
rs<-raw                                                #protect the original data
test <- rs[rs$Station_Type!="Recon",]                   # remove recon stations
rs <- subset(rs, rs$Station_Type !="Recon")
test <- rs[rs$A.Video.Readable == "Yes",]               # remove invalid videos
rs <- subset(rs, rs$A.Video.Readable == "Yes")
rs <- subset(rs, rs$Start_Depth > 0)                   # remove NA in depth
rs <- subset(rs, rs$Start_Depth < 100)                 # remove < 100 m deep
rs <- subset(rs, rs$LastOfTemp > 0)                    # remove blank water temps
rs <- subset(rs, rs$Turbidity != "Unknown")            # remove unknown turbidity values
rs <- subset(rs, rs$Turbidity!=0.5)                    # remove strange 0.5 value
rs <- subset(rs, rs$No.Readable.Frames ==41)           # remove zero readable frames
rs <- subset(rs, rs$Substrate !='Unknown')             # remove Substrate with unknowns
rs <- subset(rs, rs$Biotic_Density !='Unknown')        # remove biotic density with unknowns


#make the Substrate_Cat
summary(rs$Substrate)
length(unique(rs$Substrate))                           #90 unique values
temp=rep(NA,length(rs$Substrate))
temp[rs$Substrate=='Unknown']='Unknown'
rs$Substrate[rs$Substrate=='Unknown']=NA
rs$Substrate=as.numeric(levels(rs$Substrate)[rs$Substrate])
temp[rs$Substrate==0]='None'
temp[rs$Substrate>0 & rs$Substrate<10]='Low'           # 1-9
temp[rs$Substrate>=10 & rs$Substrate<40]='Moderate'    # 10-39
temp[rs$Substrate>=40]='High'                          # 40 & >
rs=cbind(rs,temp)
names(rs)=c(names(rs)[-dim(rs)[2]],'Substrate_Cat')
summary(rs$Substrate_Cat)



#make the Biotic_Density_Cat
summary(rs$Biotic_Density)
length(unique(rs$Biotic_Density))                     # 82 unique values
temp=rep(NA,length(rs$Biotic_Density))
temp[rs$Biotic_Density=='Unknown']='Unknown'
rs$Biotic_Density[rs$Biotic_Density=='Unknown']=NA
rs$Biotic_Density=as.numeric(levels(rs$Biotic_Density)[rs$Biotic_Density])
temp[rs$Biotic_Density==0]='None'                      # 0
temp[rs$Biotic_Density>0 & rs$Biotic_Density<10]='Low' # 1-9
temp[rs$Biotic_Density>=10 & rs$Biotic_Density<40]='Moderate'   # 10-39
temp[rs$Biotic_Density>=40]='High'                              # 40 and >
rs=cbind(rs,temp)                            
names(rs)=c(names(rs)[-dim(rs)[2]],'Biotic_Density_Cat')
summary(rs$Biotic_Density_Cat)                        



#rename startTime                   # Left from Lews code base IDK

# make sum count vector
# rs.frames <- rs[,49:89]
# rs$SumCount <- rowSums(rs.frames)

#Eliminate unnecessary columns
dat <- subset(rs,select=c(SumCount,No.Readable.Frames,Year,Turbidity,Current_Direction,
                       Current_Magnitude,Substrate_Cat,Relief,Size,Biotic_Density_Cat,
                       Biotic_Type,Biotic_Height,Start_Depth,Julian,Start_Latitude,
                       LastOfTemp,TOD))

# # Build Parameter Table
# orgnames=names(dat)
#rename to short names
names(dat)=c('SumCount','frames','y','wc','cd','cm','sc','sr','ss','bd','bt','bh','d','t','lat','temp','tod');head(dat)
shortnames=names(dat)
# 
# #create variable descriptions
# description <- c("Total count (sum of No.Readable.Frames)",
#                  "No.Readable.Frames",
#                  "Survey Year",
#                  "measure of Turbidity (how, units?)",
#                  "direction of current in reference to camera",
#                  "categorical rating for current magnitude",
#                  "categorical rating for bottom substrate",
#                  "categorical description of bottom relief",
#                  "Substrate size",
#                  "biotic density?",
#                  "biotic type",
#                  "biotic height",
#                  "Depth (m)",
#                  "Julian day of year?",
#                  "Latitude",
#                  "Bottom Temperature",
#                  "Time of day")
# 
# metadataBound <- cbind(shortnames, description)
# colnames(metadataBound) <- c("Parameter", "Description")
# metadataTable <- xtable::xtable(metadataBound, 
#                                 align="ccl",
#                                 caption="Data variables and description") 
# Metadata <- xtable::print.xtable(metadataTable, type= "latex" )



# # Table 1 data
stargazer(dat, type='text')
dat.2010<-subset(dat, dat$y==2010)
stargazer(dat.2010, type='text')
dat.2011<-subset(dat, dat$y==2011)
stargazer(dat.2011, type='text')
dat.2012<-subset(dat, dat$y==2012)
stargazer(dat.2012, type='text')
dat.2013<-subset(dat, dat$y==2013)
stargazer(dat.2013, type='text')
dat.2014<-subset(dat, dat$y==2014)
stargazer(dat.2013, type='text')
summaryBy(SumCount~y,data=dat,FUN=length)
#rescale frames
#dat$frames=dat$frames/41
#dat$frames=log(dat$frames)


# Figure for original continus variable distributions
par(mfrow=c(2,2))
#depth
hist(dat$d,breaks=seq(10,110,by=5),
     xlab="Depth (m)",
     main="Depth Distribution")
#latitude
hist(dat$lat,breaks=seq(27,36,by=0.25),
     main="Latitude Distribution",
     xlab="Latitude")
#day of year (t)
hist(dat$t,breaks=seq(110,305,by=5),
     xlab="Julian Day",
     main="Annual Sampling Distribution")
#water temperature (temp)
hist(dat$temp,breaks=seq(12,30, by=0.25),
     xlab="Temperature (C)",
     main="Temperature Distribution")


#depth
par(mfrow=c(1,2))
hist(dat$d,breaks=seq(10,110,by=5),
     xlab="Depth (m)",
     main="Depth Distribution")
#summary(dat$d)
#temp=cut(dat$d,breaks=c(14,25,41,52,115),labels=FALSE)#;temp;table(temp)
#temp=cut(dat$d,2,breaks=c(0,as.numeric(summary(dat$d))[-c(1,4)]),labels=FALSE);temp;table(temp)
temp=cut(dat$d,2,breaks=quantile(dat$d),labels=FALSE)#;temp;table(temp)
dat$d=temp
hist(dat$d,
     xlab="Depth (m)",
     main="Categorized Depth Distribution")

#latitude


par(mfrow=c(2,2))
hist(dat$lat[dat$y=='2010'],breaks=seq(27,36,by=0.25),
     xlab="Latitude",
     main="2010 Latitude Distribution")
hist(dat$lat[dat$y=='2011'],breaks=seq(27,36,by=0.25),
     xlab="Latitude",
     main="2011 Latitude Distribution")
hist(dat$lat[dat$y=='2012'],breaks=seq(27,36,by=0.25),
     xlab="Latitude",
     main="2012 Latitude Distribution")
hist(dat$lat[dat$y=='2013'],breaks=seq(27,36,by=0.25),
     xlab="Latitude",
     main="2013 Latitude Distribution")

par(mfrow=c(1,2))
hist(dat$lat,breaks=seq(27,36,by=0.25),
     main="Latitude Distribution",
     xlab="Latitude")
#summary(dat$lat)
#temp=cut(dat$lat,breaks=c(27,29.75,31.25,32.75,34,35.25),labels=FALSE)#;temp;table(temp)
#temp=cut(dat$lat,2,breaks=c(0,as.numeric(summary(dat$lat))[-c(1,4)]),labels=FALSE);temp;table(temp)
#temp=cut(dat$lat,2,breaks=quantile(dat$lat),labels=FALSE)#;temp;table(temp)
temp=cut(dat$lat,2,breaks=quantile(dat$lat,probs=seq(0,1,0.125)),labels=FALSE)#;temp;table(temp)
dat$lat=temp
# Categorical Result
hist(dat$lat,
     xlab="Latitude",
     main="Categorized Latitude Distribution")

#day of year (t)
par(mfrow=c(1,2))
hist(dat$t,breaks=seq(110,305,by=5),
     xlab="Julian Day",
     main="Annual Sampling Distribution")
#summary(dat$t)
#temp=cut(dat$t,breaks=c(##,##,##,##,##,##),labels=FALSE)#;temp;table(temp)
#temp=cut(dat$t,2,breaks=c(0,as.numeric(summary(dat$t))[-c(1,4)]),labels=FALSE);temp;table(temp)
#temp=cut(dat$t,2,breaks=quantile(dat$t),labels=FALSE)#;temp;table(temp)
temp=cut(dat$t,2,breaks=quantile(dat$t,probs=seq(0,1,0.125)),labels=FALSE)
dat$t=temp
hist(dat$t)


#water temperature (temp)
par(mfrow=c(1,2))
hist(dat$temp,breaks=seq(12.25,29.25,by=0.25),
     xlab="Temperature (C)",
     main="Temperature Distribution")
#summary(dat$temp)
#temp=cut(dat$t,breaks=c(##,##,##,##,##,##),labels=FALSE)#;temp;table(temp)
#temp=cut(dat$temp,2,breaks=c(0,as.numeric(summary(dat$temp))[-c(1,4)]),labels=FALSE);temp;table(temp)
quantile(dat$temp)  #give the divisions for the factors
temp=cut(dat$temp,2,breaks=quantile(dat$temp),labels=FALSE)#;temp;table(temp)
dat$temp=temp
hist(dat$temp)  # made into a factor of 1-4



#time of day (tod)
par(mfrow=c(1,2))
hist(dat$tod,breaks=seq(0.4,0.95,by=0.025),
     xlab="Time of Day",
     main="Time of Day Distribution")
#summary(dat$tod)
#temp=cut(dat$tod,breaks=c(##,##,##,##,##,##),labels=FALSE)#;temp;table(temp)
#temp=cut(dat$tod,2,breaks=c(0,as.numeric(summary(dat$tod))[-c(1,4)]),labels=FALSE);temp;table(temp)
quantile(dat$tod)
temp=cut(dat$tod,2,breaks=quantile(dat$tod),labels=FALSE)#;temp;table(temp)
dat$tod=temp
hist(dat$tod)

par(mfrow=c(2,2))
## DEPTH
hist(dat$d,
     main="Depth",
     xlab="Depth (m)")
## LATITUDE
hist(dat$lat,
     main="Latitude",
     xlab="Latitude (degrees)")
## DAY OF YEAR (t)
hist(dat$t,
     main="Day of the Year",
     xlab="Julian Day")
## WATER TEMP (temp)
hist(dat$temp,
     main="Temperature",
     xlab="Bottom Temperature (C)")


#Factorize variables
dat$y=factor(dat$y)
dat$d=factor(dat$d)
dat$lat=factor(dat$lat)
dat$t=factor(dat$t)
dat$temp=factor(dat$temp)
dat$tod=factor(dat$tod)

#get rid of unused factors for dat$wc
dat$wc=factor(dat$wc)
summary(dat)

plot.design(SumCount~y + wc + cd  ,data=dat)
plot.design(SumCount~ sc + bd ,data=dat)
plot.design(SumCount~d + t + lat + temp ,data=dat)


names(dat)


# So this bit evaluates the ZIP versus the ZINB **
  
#The ZINB is clearly preferred in the likelihood ratio test and fits the data better.  But you can see that neither model fit the data particularly well.


zipform=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp )
#zipform=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y  )
#zipform=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |1.  )
zipmod=zeroinfl(zipform,  dist = "poisson", link = "logit",data=dat);summary(zipmod)

nbform=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp )
#nbform=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y )
#nbform=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |1. )
nbmod=zeroinfl(nbform,  dist = "negbin", link = "logit",data=dat);summary(nbmod)

lrtest(zipmod,nbmod)




#```{r, echo=FALSE,warning=TRUE,error=TRUE, results='asis'}
# stargazer(zipmod,nbmod, type="html", title="Zero-Inflated Poisson versus Negative Binomial Models",
#           align=TRUE)
stargazer(zipmod,nbmod, type="text", 
          title="Zero-Inflated Poisson versus Negative Binomial Models",
          align=TRUE,
          column.labels=c("ZIP", "ZINB"),
          single.row=T,
          model.numbers=F)

par(mfrow=c(1,2))
hist(dat$SumCount,breaks=0:max(dat$SumCount),freq=T,right=TRUE,xlab='SumCount', main='ZIP',ylim=c(0,50), xlim=c(0,60))  
d=hist(predict(zipmod),breaks=0:max(dat$SumCount),plot=FALSE)
lines(seq(0.5,max(dat$SumCount),by=1),d$counts, col="blue",type='b') 

hist(dat$SumCount,breaks=0:max(dat$SumCount),freq=T,right=TRUE,xlab='SumCount', main='ZINB',ylim=c(0,50),xlim=c(0,60))  
d2=hist(predict(nbmod),breaks=0:max(dat$SumCount),plot=FALSE)
lines(seq(0.5,max(dat$SumCount),by=1),d2$counts, col="blue",type='b')   



## DIAGNOSTIC PLOTS
#windows(width=8,height=6,record=T)
resids=residuals(zipmod)
#cbind(fitted(zipmod),dat$SumCount,fitted(zipmod)-dat$SumCount,resids)
plot(fitted(zipmod),resids)
plot(dat$y,resids,xlab="Year",main="Residuals (zipmod)")
plot(dat$wc,resids,xlab="Water Clarity",main="Residuals (zipmod)")
plot(dat$cd,resids,xlab="Current Direction",main="Residuals (zipmod)")
plot(dat$bd,resids,xlab="Biotic Diversity",main="Residuals (zipmod)")
plot(dat$lat,resids,xlab="Latitude",main="Residuals (zipmod)")

plot(dat$SumCount,fitted(zipmod))

hist(dat$SumCount,breaks=0:max(dat$SumCount),freq=T,right=TRUE,
     xlab='Aggregate Fish Counted', 
     main='ZIP')  

d=hist(predict(zipmod),breaks=0:max(dat$SumCount),plot=FALSE)
lines(seq(0.5,max(dat$SumCount),by=1),d$counts, col="blue",type='l')      
hist(dat$SumCount,breaks=0:max(dat$SumCount),freq=T,right=TRUE,xlab='Aggregate Fish Counted', main='ZIP',ylim=c(0,50))  
lines(seq(0.5,max(dat$SumCount),by=1),d$counts, col="blue",type='b')      


#windows(width=8,height=6,record=T)
resids=residuals(nbmod)
#cbind(fitted(nbmod),dat$SumCount,fitted(nbmod)-dat$SumCount,resids)
plot(fitted(nbmod),resids)
plot(dat$SumCount,fitted(nbmod))

plot(dat$y,resids,xlab="Year",main="Residuals (nbmod)")
plot(dat$wc,resids,xlab="Water Clarity",main="Residuals (nbmod)")


plot(dat$cd,resids,xlab="Current Direction",main="Residuals (nbmod)")
plot(dat$bd,resids,xlab="Biotic Diversity",main="Residuals (nbmod)")
plot(dat$lat,resids,xlab="Latitude",main="Residuals (nbmod)")

hist(dat$SumCount,breaks=0:max(dat$SumCount),freq=T,right=TRUE,xlab='Aggregate Fish Counted', main='NB')  
d2=hist(predict(nbmod),breaks=0:max(dat$SumCount),plot=FALSE)
lines(seq(0.5,max(dat$SumCount),by=1),d2$counts, col="blue",type='l')      
hist(dat$SumCount,breaks=0:max(dat$SumCount),freq=T,right=TRUE,xlab='Aggregate Fish Counted', main='NB',ylim=c(0,50))  
lines(seq(0.5,max(dat$SumCount),by=1),d2$counts, col="blue",type='b')      

plot(dat$SumCount,fitted(nbmod))
points(dat$SumCount,fitted(zipmod),col='red',pch=19)



# So this bit allows variable selection within the ZINB**
  
#NULL formula Nbmod
#nbform=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp )

###Remove water clarity 
nbform1=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp)
#nbform1=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |1.)
nbmod1=zeroinfl(nbform1,  dist = "negbin", link = "logit",data=dat);summary(nbmod1)

###Remove cd 
nbform2=formula(SumCount~ y  + wc + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp)
#nbform2=formula(SumCount~ y  + wc + sc + bd + d + t + lat + temp |1.)
nbmod2=zeroinfl(nbform2,  dist = "negbin", link = "logit",data=dat);summary(nbmod2)

###Remove sc
nbform3=formula(SumCount~ y + wc + cd + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp )
#nbform3=formula(SumCount~ y + wc + cd + bd + d + t + lat + temp |1.)
nbmod3=zeroinfl(nbform3,  dist = "negbin", link = "logit",data=dat);summary(nbmod3)

###Remove bd
nbform4=formula(SumCount~ y + wc + cd + sc + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp )
#nbform4=formula(SumCount~ y + wc + cd + sc + d + t + lat + temp |1.)
nbmod4=zeroinfl(nbform4,  dist = "negbin", link = "logit",data=dat);summary(nbmod4)

###Remove d
nbform5=formula(SumCount~ y + wc + cd + sc + bd + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp)
#nbform5=formula(SumCount~ y + wc + cd + sc + bd + t + lat + temp |1.)
nbmod5=zeroinfl(nbform5,  dist = "negbin", link = "logit",data=dat);summary(nbmod5)

#### Remove t
nbform6=formula(SumCount~ y + wc + cd + sc + bd + d + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp)
#nbform6=formula(SumCount~ y + wc + cd + sc + bd + d + lat + temp |1.)
nbmod6=zeroinfl(nbform6,  dist = "negbin", link = "logit",data=dat);summary(nbmod6)

#### Remove Lat
nbform7=formula(SumCount~ y + wc + cd + sc + bd + d + t + temp |y + wc + cd + sc + bd + d + t + lat + temp)
#nbform7=formula(SumCount~ y + wc + cd + sc + bd + d + t + temp |1.)
nbmod7=zeroinfl(nbform7,  dist = "negbin", link = "logit",data=dat);summary(nbmod7)

#### Remove Temp
nbform8=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
#nbform8=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat  |1.)
nbmod8=zeroinfl(nbform8,  dist = "negbin", link = "logit",data=dat);summary(nbmod8)

#### Drop Year
nbform9=formula(SumCount~  wc + cd + sc + bd + d + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod9=zeroinfl(nbform9,  dist = "negbin", link = "logit",data=dat);summary(nbmod9)



# lr1 <- lrtest(nbmod1,nbmod)
# lr2 <- lrtest(nbmod2,nbmod)
# lr3 <- lrtest(nbmod3,nbmod)
# lr4 <- lrtest(nbmod4,nbmod)
# lr5 <- lrtest(nbmod5,nbmod)
# lr6 <- lrtest(nbmod6,nbmod)
# lr7 <- lrtest(nbmod7,nbmod)
# lr8 <- lrtest(nbmod8,nbmod)
# lr9 <- lrtest(nbmod9,nbmod)

lrtest(nbmod1,nbmod)
lrtest(nbmod2,nbmod)
lrtest(nbmod3,nbmod)
lrtest(nbmod4,nbmod)
lrtest(nbmod5,nbmod)
lrtest(nbmod6,nbmod)
lrtest(nbmod7,nbmod)
lrtest(nbmod8,nbmod)
lrtest(nbmod9,nbmod)

AIC(nbmod, nbmod1, nbmod2, nbmod3, nbmod4, 
    nbmod5, nbmod6,nbmod7, nbmod8, nbmod9)


AIC(nbmod)-AIC(nbmod1)
AIC(nbmod)-AIC(nbmod2)
AIC(nbmod)-AIC(nbmod3)
AIC(nbmod)-AIC(nbmod4)
AIC(nbmod)-AIC(nbmod5)
AIC(nbmod)-AIC(nbmod6)
AIC(nbmod)-AIC(nbmod7)
AIC(nbmod)-AIC(nbmod8)
AIC(nbmod)-AIC(nbmod9)

# # Following Zuur dropping logistic side parameters
# #NULL formula Nbmod
# #nbform=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp )
# 
# ###Remove wc 
# nbform10=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + t + lat + temp )
# nbmod10=zeroinfl(nbform10,  dist = "negbin", link = "logit",data=dat);summary(nbmod10)
# ###Remove cd
# nbform11=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y + wc + sc + bd + d + t + lat + temp )
# nbmod11=zeroinfl(nbform11,  dist = "negbin", link = "logit",data=dat);summary(nbmod11)
# ###Remove sc
# nbform12=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y + wc + cd + bd + d + t + lat + temp )
# nbmod12=zeroinfl(nbform12,  dist = "negbin", link = "logit",data=dat);summary(nbmod12)
# ###Remove bd
# nbform13=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + d + t + lat + temp )
# nbmod13=zeroinfl(nbform13,  dist = "negbin", link = "logit",data=dat);summary(nbmod13)
# ###Remove d
# nbform14=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd +  t + lat + temp )
# nbmod14=zeroinfl(nbform14,  dist = "negbin", link = "logit",data=dat);summary(nbmod14)
# ###Remove t
# nbform15=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d +  lat + temp )
# nbmod15=zeroinfl(nbform15,  dist = "negbin", link = "logit",data=dat);summary(nbmod15)
# ###Remove lat 
# nbform16=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + temp )
# nbmod16=zeroinfl(nbform16,  dist = "negbin", link = "logit",data=dat);summary(nbmod16)
# ###Remove temp
# nbform17=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat )
# nbmod17=zeroinfl(nbform17,  dist = "negbin", link = "logit",data=dat);summary(nbmod17)
# ###Remove y
# nbform18=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat + temp | wc + cd + sc + bd + d + t + lat + temp )
# nbmod18=zeroinfl(nbform18,  dist = "negbin", link = "logit",data=dat);summary(nbmod18)
# 
# 
# lr10 <- lrtest(nbmod11,nbmod)
# lr11 <- lrtest(nbmod12,nbmod)
# lr12 <- lrtest(nbmod13,nbmod)
# lr13 <- lrtest(nbmod14,nbmod)
# lr14 <- lrtest(nbmod15,nbmod)
# lr15 <- lrtest(nbmod16,nbmod)
# lr16 <- lrtest(nbmod17,nbmod)
# lr17 <- lrtest(nbmod18,nbmod)
# 
# 
# AIC(nbmod10, nbmod11, nbmod12, nbmod13, 
#     nbmod14, nbmod15, nbmod16,nbmod17, nbmod18)
# 
# 
# AIC(nbmod)-AIC(nbmod10)
# AIC(nbmod)-AIC(nbmod11)
# AIC(nbmod)-AIC(nbmod12)
# AIC(nbmod)-AIC(nbmod13)
# AIC(nbmod)-AIC(nbmod14)
# AIC(nbmod)-AIC(nbmod15)
# AIC(nbmod)-AIC(nbmod16)
# AIC(nbmod)-AIC(nbmod17)
# AIC(nbmod)-AIC(nbmod18)
# AIC(nbmod)-AIC(nbmod19)
# 
# DroppedTerms <- c("None", "Wc from Count", "cd from Count", "sc from Count", "bd from count",
#                  "d from count", "t from count", "lat from count", "temp from count", "y from count",
#                  "Wc from logistic", "cd from logistic", "sc from logistic", "bd from logistic",
#                  "d from logistic", "t from logistic", "lat from logistic", "temp from logistic", "y from logistic")
# aicValues <- AIC(nbmod, nbmod1, nbmod2, nbmod3, nbmod4, 
#                  nbmod5, nbmod6,nbmod7, nbmod8, nbmod9,
#                  nbmod10, nbmod11, nbmod12, nbmod13, 
#                  nbmod14, nbmod15, nbmod16,nbmod17, nbmod18)
# aicDF <- aicValues[1]
# aicVal <- aicValues[2]
# 
# chiSqr <- c("NA", lr1$Chisq[2], lr2$Chisq[2],lr3$Chisq[2],lr4$Chisq[2],lr5$Chisq[2],
#             lr6$Chisq[2],lr7$Chisq[2],lr8$Chisq[2],lr9$Chisq[2],lr10$Chisq[2],lr11$Chisq[2],
#             lr12$Chisq[2],lr13$Chisq[2],lr14$Chisq[2],lr15$Chisq[2],lr16$Chisq[2],lr17$Chisq[2])
# PrChiSqr <- c("NA", lr1$Pr[2], lr2$Pr[2],lr3$Pr[2],lr4$Pr[2],lr5$Pr[2],
#             lr6$Pr[2],lr7$Pr[2],lr8$Pr[2],lr9$Pr[2],lr10$Pr[2],lr11$Pr[2],
#             lr12$Pr[2],lr13$Pr[2],lr14$Pr[2],lr15$Pr[2],lr16$Pr[2],lr17$Pr[2])
# 
# modelOptTable <- cbind(aicDF, aicVal, chiSqr, PrChiSqr)

### Second Round ##### 
# Model nbmod1 is the best performing so step down to it
# nbform1=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp)

###Remove cd 
nbform1a=formula(SumCount~ y  + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp)
nbmod1a=zeroinfl(nbform1a,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a)
###Remove sc 
nbform1b=formula(SumCount~ y  + cd + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp)
nbmod1b=zeroinfl(nbform1b,  dist = "negbin", link = "logit",data=dat);summary(nbmod1b)
###Remove bd 
nbform1c=formula(SumCount~ y  + cd + sc + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp)
nbmod1c=zeroinfl(nbform1c,  dist = "negbin", link = "logit",data=dat);summary(nbmod1c)
###Remove d 
nbform1d=formula(SumCount~ y  + cd + sc + bd + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp)
nbmod1d=zeroinfl(nbform1d,  dist = "negbin", link = "logit",data=dat);summary(nbmod1d)
###Remove t 
nbform1e=formula(SumCount~ y  + cd + sc + bd + d + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp)
nbmod1e=zeroinfl(nbform1e,  dist = "negbin", link = "logit",data=dat);summary(nbmod1e)
###Remove lat 
nbform1f=formula(SumCount~ y  + cd + sc + bd + d + t + temp |y + wc + cd + sc + bd + d + t + lat + temp)
nbmod1f=zeroinfl(nbform1f,  dist = "negbin", link = "logit",data=dat);summary(nbmod1f)
###Remove temp 
nbform1g=formula(SumCount~ y  + cd + sc + bd + d + t + lat |y + wc + cd + sc + bd + d + t + lat + temp)
nbmod1g=zeroinfl(nbform1g,  dist = "negbin", link = "logit",data=dat);summary(nbmod1g)
###Remove y 
nbform1h=formula(SumCount~ cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp)
nbmod1h=zeroinfl(nbform1h,  dist = "negbin", link = "logit",data=dat);summary(nbmod1h)

lrtest(nbmod1a,nbmod1)
lrtest(nbmod1b,nbmod1)
lrtest(nbmod1c,nbmod1)
lrtest(nbmod1d,nbmod1)
lrtest(nbmod1e,nbmod1)
lrtest(nbmod1g,nbmod1)
lrtest(nbmod1g,nbmod1)
lrtest(nbmod1h,nbmod1)


AIC(nbmod1, nbmod1a, nbmod1b, nbmod1c, nbmod1d,
    nbmod1e, nbmod1f, nbmod1g, nbmod1h)


AIC(nbmod1)-AIC(nbmod1a)
AIC(nbmod1)-AIC(nbmod1b)
AIC(nbmod1)-AIC(nbmod1c)
AIC(nbmod1)-AIC(nbmod1d)
AIC(nbmod1)-AIC(nbmod1e)
AIC(nbmod1)-AIC(nbmod1f)
AIC(nbmod1)-AIC(nbmod1g)
AIC(nbmod1)-AIC(nbmod1h)


# Nothing dropped in the second round
# Moving to the logistic side of the model with nbmod1

### Third Round ##### 
# Model nbmod1 is the best performing so stay on it
# nbform1=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat + temp)

###Remove wc
nbform1a=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + t + lat + temp)
nbmod1a=zeroinfl(nbform1a,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a)
###Remove cd 
nbform1b=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + wc + sc + bd + d + t + lat + temp)
nbmod1b=zeroinfl(nbform1b,  dist = "negbin", link = "logit",data=dat);summary(nbmod1b)
###Remove sc 
nbform1c=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + wc + cd + bd + d + t + lat + temp)
nbmod1c=zeroinfl(nbform1c,  dist = "negbin", link = "logit",data=dat);summary(nbmod1c)
###Remove bd 
nbform1d=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + d + t + lat + temp)
nbmod1d=zeroinfl(nbform1d,  dist = "negbin", link = "logit",data=dat);summary(nbmod1d)
###Remove d 
nbform1e=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + t + lat + temp)
nbmod1e=zeroinfl(nbform1e,  dist = "negbin", link = "logit",data=dat);summary(nbmod1e)
###Remove t 
nbform1f=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + lat + temp)
nbmod1f=zeroinfl(nbform1f,  dist = "negbin", link = "logit",data=dat);summary(nbmod1f)
###Remove lat 
nbform1g=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + temp)
nbmod1g=zeroinfl(nbform1g,  dist = "negbin", link = "logit",data=dat);summary(nbmod1g)
###Remove temp 
nbform1h=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + wc + cd + sc + bd + d + t + lat )
nbmod1h=zeroinfl(nbform1h,  dist = "negbin", link = "logit",data=dat);summary(nbmod1h)
###Remove y 
nbform1i=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp | wc + cd + sc + bd + d + t + lat + temp)
nbmod1i=zeroinfl(nbform1i,  dist = "negbin", link = "logit",data=dat);summary(nbmod1i)

lrtest(nbmod1a,nbmod1)
lrtest(nbmod1b,nbmod1)
lrtest(nbmod1c,nbmod1)
lrtest(nbmod1d,nbmod1)
lrtest(nbmod1e,nbmod1)
lrtest(nbmod1f,nbmod1)
lrtest(nbmod1g,nbmod1)
lrtest(nbmod1h,nbmod1)
lrtest(nbmod1i,nbmod1)


AIC(nbmod1, nbmod1a, nbmod1b, nbmod1c, nbmod1d,
    nbmod1e, nbmod1f, nbmod1g, nbmod1h, nbmod1i)


AIC(nbmod1)-AIC(nbmod1a)
AIC(nbmod1)-AIC(nbmod1b)
AIC(nbmod1)-AIC(nbmod1c)
AIC(nbmod1)-AIC(nbmod1d)
AIC(nbmod1)-AIC(nbmod1e)
AIC(nbmod1)-AIC(nbmod1f)
AIC(nbmod1)-AIC(nbmod1g)
AIC(nbmod1)-AIC(nbmod1h)
AIC(nbmod1)-AIC(nbmod1i)



# After a round of logistic drops looks like the nbform1a is  best, so step down to that

### Forth Round ##### 
# Model nbmod1a is the best performing so stay on it
# nbform1a=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + t + lat + temp)

###Remove cd 
nbform1a1=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + sc + bd + d + t + lat + temp)
nbmod1a1=zeroinfl(nbform1a1,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a1)
###Remove sc
nbform1a2=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + bd + d + t + lat + temp)
nbmod1a2=zeroinfl(nbform1a2,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a2)
###Remove bd 
nbform1a3=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + d + t + lat + temp)
nbmod1a3=zeroinfl(nbform1a3,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a3)
###Remove d 
nbform1a4=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + t + lat + temp)
nbmod1a4=zeroinfl(nbform1a4,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a4)
###Remove t 
nbform1a5=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + lat + temp)
nbmod1a5=zeroinfl(nbform1a5,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a5)
###Remove lat 
nbform1a6=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + t +  temp)
nbmod1a6=zeroinfl(nbform1a6,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a6)
###Remove temp 
nbform1a7=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + t + lat )
nbmod1a7=zeroinfl(nbform1a7,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a7)
###Remove y 
nbform1a8=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp | cd + sc + bd + d + t + lat + temp)
nbmod1a8=zeroinfl(nbform1a8,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a8)

lrtest(nbmod1a1,nbmod1a)
lrtest(nbmod1a2,nbmod1a)
lrtest(nbmod1a3,nbmod1a)
lrtest(nbmod1a4,nbmod1a)
lrtest(nbmod1a5,nbmod1a)
lrtest(nbmod1a6,nbmod1a)
lrtest(nbmod1a7,nbmod1a)
lrtest(nbmod1a8,nbmod1a)



AIC(nbmod1a, nbmod1a1, nbmod1a2, nbmod1a3, nbmod1a4,
    nbmod1a5, nbmod1a6, nbmod1a7, nbmod1a8)


AIC(nbmod1a)-AIC(nbmod1a)
AIC(nbmod1a)-AIC(nbmod1b)
AIC(nbmod1a)-AIC(nbmod1c)
AIC(nbmod1a)-AIC(nbmod1d)
AIC(nbmod1a)-AIC(nbmod1e)
AIC(nbmod1a)-AIC(nbmod1f)
AIC(nbmod1a)-AIC(nbmod1g)
AIC(nbmod1a)-AIC(nbmod1h)
AIC(nbmod1a)-AIC(nbmod1i)


# After a round 4 of logistic drops looks like the nbform1a5 is  best, so step down to that
### Fifth Round ##### 
# Model nbmod1a is the best performing so stay on it
# nbform1a5=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + lat + temp)

###Remove cd 
nbform1a5a=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + sc + bd + d + lat + temp)
nbmod1a5a=zeroinfl(nbform1a5a,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a5a)
###Remove sc 
nbform1a5b=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + bd + d + lat + temp)
nbmod1a5b=zeroinfl(nbform1a5b,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a5b)
###Remove bd 
nbform1a5c=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + d + lat + temp)
nbmod1a5c=zeroinfl(nbform1a5c,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a5c)
###Remove d 
nbform1a5d=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + lat + temp)
nbmod1a5d=zeroinfl(nbform1a5d,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a5d)
###Remove lat 
nbform1a5e=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + temp)
nbmod1a5e=zeroinfl(nbform1a5e,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a5e)
###Remove temp 
nbform1a5f=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + lat )
nbmod1a5f=zeroinfl(nbform1a5f,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a5f)
###Remove y 
nbform1a5g=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |cd + sc + bd + d + lat + temp)
nbmod1a5g=zeroinfl(nbform1a5g,  dist = "negbin", link = "logit",data=dat);summary(nbmod1a5g)

lrtest(nbmod1a5a,nbmod1a5)
lrtest(nbmod1a5b,nbmod1a5)
lrtest(nbmod1a5c,nbmod1a5)
lrtest(nbmod1a5d,nbmod1a5)
lrtest(nbmod1a5e,nbmod1a5)
lrtest(nbmod1a5f,nbmod1a5)
lrtest(nbmod1a5g,nbmod1a5)

AIC(nbmod1a5, nbmod1a5a, nbmod1a5b, nbmod1a5c, nbmod1a5d,
    nbmod1a5e, nbmod1a5f, nbmod1a5g)

AIC(nbmod1a5)-AIC(nbmod1a5a)
AIC(nbmod1a5)-AIC(nbmod1a5b)
AIC(nbmod1a5)-AIC(nbmod1a5c)
AIC(nbmod1a5)-AIC(nbmod1a5d)
AIC(nbmod1a5)-AIC(nbmod1a5e)
AIC(nbmod1a5)-AIC(nbmod1a5f)
AIC(nbmod1a5)-AIC(nbmod1a5g)




nbbest=nbmod1a5
# nbform1a5=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + lat + temp)

resids=residuals(nbbest, type="pearson")



#cbind(fitted(nbbest),dat$SumCount,fitted(nbbest)-dat$SumCount,resids)
par(mfrow=c(1,2))
plot(fitted(nbbest),resids, ylab="Pearsons Residuals", xlab="Fitted Values") #pearson resids vs fitted values
plot(dat$SumCount,fitted(nbbest), ylab="Fitted Values", xlab="Original Values") #pearson values vs original data

# A single figure
par(mfrow=c(1,3))
plot(dat$y,resids,xlab="Year",main="Residuals (nbbest)")
plot(dat$t,resids,xlab="Season",main="Residuals (nbbest)")
plot(dat$lat,resids,xlab="Latitude",main="Residuals (nbbest)")

par(mfrow=c(3,2))
#plot(dat$wc,resids,xlab="Water Clarity",main="Residuals (nbbest)")
plot(dat$cd,resids,xlab="Current Direction",main="Residuals (nbbest)")
plot(dat$sc,resids,xlab="Substrate Composition",main="Residuals (nbbest)")
plot(dat$bd,resids,xlab="Biotic Diversity",main="Residuals (nbbest)")
plot(dat$d,resids,xlab="Depth",main="Residuals (nbbest)")
plot(dat$temp,resids,xlab="Temperature",main="Residuals (nbbest)")


par(mfrow=c(1,2))
#full histogram, commented out in favor of the range limited figure below
hist(dat$SumCount,breaks=0:max(dat$SumCount),freq=T,right=TRUE,xlab='Aggregate Fish Counted', main='NB')  
d3=hist(predict(nbbest),breaks=0:max(dat$SumCount),plot=FALSE)
lines(seq(0.5,max(dat$SumCount),by=1),d3$counts, col="blue",type='l')   

# This is where Lew limited the range.
hist(dat$SumCount,breaks=0:max(dat$SumCount),freq=T,right=TRUE,xlab='Aggregate Fish Counted', main='NB',ylim=c(0,100))  
lines(seq(0.5,max(dat$SumCount),by=1),d3$counts, col="blue",type='l')      
lines(seq(0.5,max(dat$SumCount),by=1),d2$counts, col="red",type='l')      

# Model Fit figure
par(mfrow=c(1,2))
hist(dat$SumCount,
     breaks=0:max(dat$SumCount),
     freq=T,right=TRUE,
     xlab='Sum Count',
     col ="lightgray",
     bor= "lightgray",
     main='NB')  
d3=hist(predict(nbbest),breaks=0:max(dat$SumCount),plot=FALSE)
lines(seq(0.5,max(dat$SumCount),by=1),d3$counts, col="blue",type='l')      
hist(dat$SumCount,
     breaks=0:max(dat$SumCount),
     freq=T,right=TRUE,
     xlab='Sum Count', 
     main='NB',
     ylim=c(0,50), xlim=c(0, 60))  
lines(seq(0.5,max(dat$SumCount),by=1),d3$counts, col="blue",type='b')      
#lines(seq(0.5,max(dat$SumCount),by=1),d2$counts, col="red",type='l')      

# final model form for FYI
# nbform1a5=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + lat + temp)

# Generating a data frame of all values
new.dat=expand.grid(y=levels(dat$y),
                    #wc=levels(dat$wc),
                    #cm=levels(dat$cm),
                    cd=levels(dat$cd),
                    sc=levels(dat$sc),
                    bd=levels(dat$bd),
                    #bt=levels(dat$bt),
                    #bh=levels(dat$bh),
                    d=levels(dat$d),
                    t=levels(dat$t),
                    lat=levels(dat$lat),
                    temp=levels(dat$temp)
                    #tod=levels(dat$tod),
                    #frames=1))
)
new.dat=cbind(new.dat,predict(nbbest,new.dat))
names(new.dat)[dim(new.dat)[2]]="Predicted"
resvec=summaryBy(Predicted~y,data=new.dat,FUN=mean)[,2]
index = resvec/mean(resvec)
#index = resvec/max(resvec)

plot(resvec,type='b',ylim=c(0,10))
plot(index,type='b',ylim=c(0,2))

resvec=summaryBy(Predicted~lat,data=new.dat,FUN=mean)[,2]
index = resvec/mean(resvec)
plot(index,type='b',ylim=c(0,4))

resvec=summaryBy(Predicted~lat,data=new.dat,FUN=mean)[,2]
plot(resvec,type='l',ylim=c(0,40),ylab='lat')
resvec=summaryBy(Predicted~lat,data=new.dat[new.dat$y=='2010',],FUN=mean)[,2]
lines(resvec,type='b',ylim=c(0,4),lty=2)
resvec=summaryBy(Predicted~lat,data=new.dat[new.dat$y=='2011',],FUN=mean)[,2]
lines(resvec,type='b',ylim=c(0,4),lty=3)
resvec=summaryBy(Predicted~lat,data=new.dat[new.dat$y=='2012',],FUN=mean)[,2]
lines(resvec,type='b',ylim=c(0,4),lty=4)
resvec=summaryBy(Predicted~lat,data=new.dat[new.dat$y=='2013',],FUN=mean)[,2]
lines(resvec,type='b',ylim=c(0,4),lty=5)

### Needs a legend

#set up data objects and specify number of bootstrap replications
ptm <- proc.time()
boots=500
# nbform1a5=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + lat + temp)

names(dat)
org.dat=dat[,c(1,3,5,7,10,13,14,15,16)];head(org.dat)
boot.dat=org.dat
numyrs=length(levels(org.dat$y));numyrs
index.boot=matrix(NA,nrow=boots,ncol=numyrs)
predmean.boot=index.boot
#yr.samp=summaryBy(~y,data=org.dat,FUN=length)[,2]   #Orginal
yr.samp=summaryBy(SumCount~y,data=org.dat,FUN=length)[,2]

#start bootstrap loop
for(boot in 1:boots){
  
  #Get bootstrap data for current replicate
  for(i in 1:numyrs){
    yr.rows=org.dat$y==levels(org.dat$y)[i]
    sub.dat=org.dat[yr.rows,]
    yr.boot=round(runif(yr.samp[i],1,yr.samp[i]),0)
    boot.dat[yr.rows,]=sub.dat[yr.boot,]
  }
  
  # Now Recalibrate 2010 Data 
  cal.prop=0.47
  #table(boot.dat$y,boot.dat$SumCount)
  #dat.2010=boot.dat[boot.dat$y==2010,]
  #cal.SumCount=rep(NA,length(dat.2010$SumCount))
  #for(i in 1:length(dat.2010$SumCount)){cal.SumCount[i]=rbinom(1,dat.2010$SumCount[i],cal.prop)}#;plot(cal.SumCount,dat.2010$SumCount,xlim=c(0,250),ylim=c(0,250))
  #  boot.dat[boot.dat$y==2010,1]= cal.SumCount 
  
  # nbform1a5=formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + lat + temp)
  #make a function to compute the index and return either a valid index or NA conditional on if the model converges
  getindex=function(){
    #define and fit model for current replicate.  Use the "try" function so that can continue to run if model does not converge on a particular replication  
    f2 =formula(SumCount~ y  + cd + sc + bd + d + t + lat + temp |y + cd + sc + bd + d + lat + temp)
    Nb2 = try(zeroinfl(f2, dist = "negbin", link = "logit", data = boot.dat)); Nb2
    
    #see if model converged, if it did, return the mean year effect over all covariate combinations
    # If it did not converge, return an NA vector for the year effects
    if (class(Nb2) != "try-error"){
      #Predict the year effect (index) by predicting for each covariate level and compute the mean
      new.dat=expand.grid(y=levels(dat$y),
                          #wc=levels(dat$wc),
                          #cm=levels(dat$cm),
                          cd=levels(dat$cd),
                          sc=levels(dat$sc),
                          bd=levels(dat$bd),
                          #bt=levels(dat$bt),
                          #bh=levels(dat$bh),
                          d=levels(dat$d),
                          t=levels(dat$t),
                          lat=levels(dat$lat),
                          temp=levels(dat$temp)
                          #tod=levels(dat$tod),
                          #frames=1)
      )
      new.dat=cbind(new.dat,predict(Nb2,new.dat))
      names(new.dat)[9]="Predicted"
      resvec=summaryBy(Predicted~y,data=new.dat,FUN=mean)[,2]
      return(resvec)
    }  
    resvec=rep(NA,numyrs)
    return(resvec)
  }
  
  #Call the getindex function for each bootstrap replicate
  predmean.boot[boot,]=getindex()
  
  #simple calibration method
  predmean.boot[boot,1]=predmean.boot[boot,1]*cal.prop
  
  index.boot[boot,]=predmean.boot[boot,]/mean(predmean.boot[boot,])
}
save.image("C:\\Users\\Kevin.Purcell\\Documents\\GitHub\\SEDAR41_RS_VI\\analysis\\SEDAR41_RS_VI.RData")

#windows(width=8,height=6,record=T)
medianidex=apply(predmean.boot[!is.na(predmean.boot[,1]),],2,quantile,c(0.5))
resvec=summaryBy(Predicted~y,data=new.dat,FUN=mean)[,2]
plot(resvec)
lines(medianidex)
#exclude runs that did not converge and compute the convergence rate
conv.index=index.boot[!is.na(index.boot[,1]),]
convrate=dim(conv.index)[1]/dim(index.boot)[1];convrate

#Calculate the CV by year

#Make a data frame from the boot strap output
cv.dat<-as.data.frame(index.boot)
cv.dat2<-cv.dat[complete.cases(cv.dat),]

# Make a function for coef vari
co.var<- function(x) (sd(x)/mean(x))

# Call co.var on each column of the bootstrap output df
cv.2010<-co.var(cv.dat2$V1)  
cv.2011<-co.var(cv.dat2$V2)
cv.2012<-co.var(cv.dat2$V3)
cv.2013<-co.var(cv.dat2$V4)

# Calculate the proportion positive raw values
temp<-subset(raw, raw$Year==2010)
PrP_2010<-length(temp$MeanCount[temp$MeanCount>0])/length(temp$MeanCount)
temp<-subset(raw, raw$Year==2011)
PrP_2011<-length(temp$MeanCount[temp$MeanCount>0])/length(temp$MeanCount)
temp<-subset(raw, raw$Year==2012)
PrP_2012<-length(temp$MeanCount[temp$MeanCount>0])/length(temp$MeanCount)
temp<-subset(raw, raw$Year==2013)
PrP_2013<-length(temp$MeanCount[temp$MeanCount>0])/length(temp$MeanCount)



(proc.time() - ptm)[3]/60/60
(proc.time() - ptm)[3]/60
(proc.time() - ptm)[3]

#conv.index
cI=apply(conv.index,2,quantile,c(0.025,0.5,0.975))
nomcpue<-summaryBy(SumCount ~ y, dat=dat, FUN=mean)
nomcpue.std <- nomcpue$SumCount.mean/mean(nomcpue$SumCount)
nomcpue.std <- as.data.frame(nomcpue.std)


matplot(t(cI),type='l',ylim=c(0,2.5),xaxt="n",col=c(1,2,1),lty=c(3,1,3),lwd=c(1,2,1),ylab='Relative CPUE')
axis(1,at=1:4,labels=2010:2013)
index = resvec/mean(resvec)
#lines(index,type='l',ylim=c(0,2),col='green',lwd=2)
lines(nomcpue.std$nomcpue.std, type='b', col='blue', pch=20, lty=1,lwd=1)
#approxgam=c(.37,.11,.21)
#lines(approxgam/mean(approxgam),lwd=2,col='blue')
legend('topright',legend=c('Standardized Index','Bootstrap CI','Nominal'),
       lty=c(1,3,1),col=c('red','black','blue'),lwd=c(2,1,1), bty="n")


sessionInfo()
save.image("C:\\Users\\Kevin.Purcell\\Documents\\GitHub\\SEDAR41_RS_VI\\analysis\\SEDAR41_RS_VI.RData")
