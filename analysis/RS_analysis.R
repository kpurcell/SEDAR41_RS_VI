
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
# stargazer(zipmod,nbmod, type="text", 
#           title="Zero-Inflated Poisson versus Negative Binomial Models",
#           align=TRUE,
#           column.labels=c("ZIP", "ZINB"),
#           single.row=T,
#           model.numbers=F)
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


lrtest(nbmod1,nbmod)
lrtest(nbmod2,nbmod)
lrtest(nbmod3,nbmod)
lrtest(nbmod4,nbmod)
lrtest(nbmod5,nbmod)
lrtest(nbmod6,nbmod)
lrtest(nbmod7,nbmod)
lrtest(nbmod8,nbmod)
lrtest(nbmod9,nbmod)

AIC(nbmod)
AIC(nbmod1)
AIC(nbmod2)
AIC(nbmod3)
AIC(nbmod4)
AIC(nbmod5)
AIC(nbmod6)
AIC(nbmod7)
AIC(nbmod8)
AIC(nbmod9)

AIC(nbmod)-AIC(nbmod1)
AIC(nbmod)-AIC(nbmod2)
AIC(nbmod)-AIC(nbmod3)
AIC(nbmod)-AIC(nbmod4)
AIC(nbmod)-AIC(nbmod5)
AIC(nbmod)-AIC(nbmod6)
AIC(nbmod)-AIC(nbmod7)
AIC(nbmod)-AIC(nbmod8)
AIC(nbmod)-AIC(nbmod9)
```


# Model 8 is the best performing so step down to it


```{r, cache=TRUE, echo=FALSE, warning=FALSE, error=FALSE, results='hide'}
#nbform8=formula(SumCount~ y + wc + cd + sc + bd + d + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )

###Remove water clarity and current variables from mean
nbform8a=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod10=zeroinfl(nbform8a,  dist = "negbin", link = "logit",data=dat);summary(nbmod10)

### Remove cd
nbform8b=formula(SumCount~ y + wc + sc + bd + d + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod11=zeroinfl(nbform8b,  dist = "negbin", link = "logit",data=dat);summary(nbmod11)

###Remove sc from mean
nbform8c=formula(SumCount~ y + wc + cd + bd + d + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod12=zeroinfl(nbform8c,  dist = "negbin", link = "logit",data=dat);summary(nbmod12)

### remove bd
nbform8d=formula(SumCount~ y + wc + cd + sc + d + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod13=zeroinfl(nbform8d,  dist = "negbin", link = "logit",data=dat);summary(nbmod13)

###Remove depth from mean
nbform8e=formula(SumCount~ y + wc + cd + sc + bd + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod14=zeroinfl(nbform8e,  dist = "negbin", link = "logit",data=dat);summary(nbmod14)

#### Remove Season From mean 
nbform8f=formula(SumCount~ y + wc + cd + sc + bd + d + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod15=zeroinfl(nbform8f,  dist = "negbin", link = "logit",data=dat);summary(nbmod15)

#### Remove Latitude From mean 
nbform8g=formula(SumCount~ y + wc + cd + sc + bd + d + t |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod16=zeroinfl(nbform8g,  dist = "negbin", link = "logit",data=dat);summary(nbmod16)

### remove year
nbform8h=formula(SumCount~ wc + cd + sc + bd + d + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod17=zeroinfl(nbform8h,  dist = "negbin", link = "logit",data=dat);summary(nbmod17)


lrtest(nbmod10,nbmod8)
lrtest(nbmod11,nbmod8)
lrtest(nbmod12,nbmod8)
lrtest(nbmod13,nbmod8)
lrtest(nbmod14,nbmod8)
lrtest(nbmod15,nbmod8)
lrtest(nbmod16,nbmod8)
lrtest(nbmod17,nbmod8)

AIC(nbmod10)
AIC(nbmod11)
AIC(nbmod12)
AIC(nbmod13)
AIC(nbmod14)
AIC(nbmod15)
AIC(nbmod16)
AIC(nbmod17)

AIC(nbmod8)-AIC(nbmod10)
AIC(nbmod8)-AIC(nbmod11)
AIC(nbmod8)-AIC(nbmod12)
AIC(nbmod8)-AIC(nbmod13)
AIC(nbmod8)-AIC(nbmod14)
AIC(nbmod8)-AIC(nbmod15)
AIC(nbmod8)-AIC(nbmod16)
AIC(nbmod8)-AIC(nbmod17)

```


# Model 8a is the best performing so step down to it

```{r, cache=TRUE, echo=FALSE, warning=FALSE, error=FALSE}
##nbform8a=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )


###Remove  and current variables from mean
nbform8a1=formula(SumCount~ y + sc + bd + d + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod18=zeroinfl(nbform8a1,  dist = "negbin", link = "logit",data=dat)
summary(nbmod18)

###Remove benthic variables from mean
nbform8a2=formula(SumCount~ y + cd + bd + d + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod19=zeroinfl(nbform8a2,  dist = "negbin", link = "logit",data=dat)
summary(nbmod19)

nbform8a3=formula(SumCount~ y + cd + sc + d + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod20=zeroinfl(nbform8a3,  dist = "negbin", link = "logit",data=dat)
summary(nbmod20)

###Remove depth from mean
nbform8a4=formula(SumCount~ y + cd + sc + bd + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod21=zeroinfl(nbform8a4,  dist = "negbin", link = "logit",data=dat)
summary(nbmod21)

#### Remove Season From mean 
nbform8a5=formula(SumCount~ y + cd + sc + bd + d + lat  |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod22=zeroinfl(nbform8a5,  dist = "negbin", link = "logit",data=dat)
summary(nbmod22)

#### Remove Latitude From mean 
nbform8a6=formula(SumCount~ y + cd + sc + bd + d + t   |y + wc + cd + sc + bd + d + t + lat + temp )
nbmod23=zeroinfl(nbform8a6,  dist = "negbin", link = "logit",data=dat)
summary(nbmod23)



lrtest(nbmod18,nbmod10)
lrtest(nbmod19,nbmod10)
lrtest(nbmod20,nbmod10)
lrtest(nbmod21,nbmod10)
lrtest(nbmod22,nbmod10)
lrtest(nbmod23,nbmod10)

AIC(nbmod18)
AIC(nbmod19)
AIC(nbmod20)
AIC(nbmod21)
AIC(nbmod22)
AIC(nbmod23)

AIC(nbmod10)-AIC(nbmod18)
AIC(nbmod10)-AIC(nbmod19)
AIC(nbmod10)-AIC(nbmod20)
AIC(nbmod10)-AIC(nbmod21)
AIC(nbmod10)-AIC(nbmod22)
AIC(nbmod10)-AIC(nbmod23)
```



# Model form 10 is the best performing for the positive part of the model

```{r, cache=TRUE, echo=FALSE, warning=FALSE, error=FALSE}
##nbmod10==>nbform8a=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + wc + cd + sc + bd + d + t + lat + temp )


##Now conduct variable selection on the binomial part
nbform8aLog1=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + cd + sc + bd + d + t + lat + temp)
nbmod24=zeroinfl(nbform8aLog1,  dist = "negbin", link = "logit",data=dat)
summary(nbmod24)

###Remove year from mean
nbform8aLog2=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + wc + sc + bd + d + t + lat + temp)
nbmod25=zeroinfl(nbform8aLog2,  dist = "negbin", link = "logit",data=dat)
summary(nbmod25)

###Remove current variables from mean
nbform8aLog3=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + wc + cd + bd + d + t + lat + temp)
nbmod26=zeroinfl(nbform8aLog3,  dist = "negbin", link = "logit",data=dat)
summary(nbmod26)

###Remove benthic variables from mean
nbform8aLog4=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + wc + cd + sc + d + t + lat + temp)
nbmod27=zeroinfl(nbform8aLog4,  dist = "negbin", link = "logit",data=dat)
summary(nbmod27)

nbform8aLog5=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + wc + cd + sc + bd + t + lat + temp)
nbmod28=zeroinfl(nbform8aLog5,  dist = "negbin", link = "logit",data=dat)
summary(nbmod28)

###Remove depth from mean
nbform8aLog6=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + wc + cd + sc + bd + d + lat + temp)
nbmod29=zeroinfl(nbform8aLog6,  dist = "negbin", link = "logit",data=dat)
summary(nbmod29)

#### Remove Season From mean 
nbform8aLog7=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + wc + cd + sc + bd + d + t + temp)
nbmod30=zeroinfl(nbform8aLog7,  dist = "negbin", link = "logit",data=dat)
summary(nbmod30)

#### Remove Latitude From mean 
nbform8aLog8=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + wc + cd + sc + bd + d + t + lat )
nbmod31=zeroinfl(nbform8aLog8,  dist = "negbin", link = "logit",data=dat)
summary(nbmod31)

nbform8aLog9=formula(SumCount ~ y + cd + sc + bd + d + t + lat | wc + cd + sc + bd + d + t + lat + temp)
nbmod32=zeroinfl(nbform8aLog9, dist = "negbin", link = "logit", data=dat)
summary(nbmod32)



lrtest(nbmod24,nbmod10)
lrtest(nbmod25,nbmod10)
lrtest(nbmod26,nbmod10)
lrtest(nbmod27,nbmod10)
lrtest(nbmod28,nbmod10)
lrtest(nbmod29,nbmod10)
lrtest(nbmod30,nbmod10)
lrtest(nbmod31,nbmod10)
lrtest(nbmod32,nbmod10)


AIC(nbmod24)
AIC(nbmod25)
AIC(nbmod26)
AIC(nbmod27)
AIC(nbmod28)
AIC(nbmod29)
AIC(nbmod30)
AIC(nbmod31)
AIC(nbmod32)

AIC(nbmod10)-AIC(nbmod24)
AIC(nbmod10)-AIC(nbmod25)
AIC(nbmod10)-AIC(nbmod26)
AIC(nbmod10)-AIC(nbmod27)
AIC(nbmod10)-AIC(nbmod28)
AIC(nbmod10)-AIC(nbmod29)
AIC(nbmod10)-AIC(nbmod30)
AIC(nbmod10)-AIC(nbmod31)
AIC(nbmod10)-AIC(nbmod32)



# nbmod24 new start one more round
#nbform8aLog1=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + cd + sc + bd + d + t + lat + temp)

nbform8aLog10=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + sc + bd + d + t + lat + temp)
nbmod33=zeroinfl(nbform8aLog10,  dist = "negbin", link = "logit",data=dat)
summary(nbmod33)

nbform8aLog11=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + cd + bd + d + t + lat + temp)
nbmod34=zeroinfl(nbform8aLog11,  dist = "negbin", link = "logit",data=dat)
summary(nbmod34)

nbform8aLog12=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + cd + sc + d + t + lat + temp)
nbmod35=zeroinfl(nbform8aLog12,  dist = "negbin", link = "logit",data=dat)
summary(nbmod35)

nbform8aLog13=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + cd + sc + bd + t + lat + temp)
nbmod36=zeroinfl(nbform8aLog13,  dist = "negbin", link = "logit",data=dat)
summary(nbmod36)

nbform8aLog14=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + cd + sc + bd + d + lat + temp)
nbmod37=zeroinfl(nbform8aLog14,  dist = "negbin", link = "logit",data=dat)
summary(nbmod37)

nbform8aLog15=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + cd + sc + bd + d + t + temp)
nbmod38=zeroinfl(nbform8aLog15,  dist = "negbin", link = "logit",data=dat)
summary(nbmod38)

nbform8aLog16=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + cd + sc + bd + d + t + lat)
nbmod39=zeroinfl(nbform8aLog16,  dist = "negbin", link = "logit",data=dat)
summary(nbmod39)



lrtest(nbmod33, nbmod24)
lrtest(nbmod34, nbmod24)
lrtest(nbmod35, nbmod24)
lrtest(nbmod36, nbmod24)
lrtest(nbmod37, nbmod24)
lrtest(nbmod38, nbmod24)
lrtest(nbmod39, nbmod24)

AIC(nbmod33)
AIC(nbmod34)
AIC(nbmod35)
AIC(nbmod36)
AIC(nbmod37)
AIC(nbmod38)
AIC(nbmod39)


AIC(nbmod24)-AIC(nbmod33)
AIC(nbmod24)-AIC(nbmod34)
AIC(nbmod24)-AIC(nbmod35)
AIC(nbmod24)-AIC(nbmod36)
AIC(nbmod24)-AIC(nbmod37)
AIC(nbmod24)-AIC(nbmod38)
AIC(nbmod24)-AIC(nbmod39)


#nbmod37
#nbform8aLog14=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + cd + sc + bd + d + lat + temp)
nbform8aLog19=formula(SumCount~ y + cd + sc + bd + d + t + lat  | y + sc + bd + d + lat + temp)
nbmod41=zeroinfl(nbform8aLog19,  dist = "negbin", link = "logit",data=dat)
summary(nbmod41)

nbform8aLog20=formula(SumCount~ y + cd + sc + bd + d + t + lat  | y + cd + bd + d + lat + temp)
nbmod42=zeroinfl(nbform8aLog20,  dist = "negbin", link = "logit",data=dat)
summary(nbmod42)

nbform8aLog21=formula(SumCount~ y + cd + sc + bd + d + t + lat  | y + cd + sc + d + lat + temp)
nbmod43=zeroinfl(nbform8aLog21,  dist = "negbin", link = "logit",data=dat)
summary(nbmod43)

nbform8aLog22=formula(SumCount~ y + cd + sc + bd + d + t + lat  | y + cd + sc + bd + lat + temp)
nbmod44=zeroinfl(nbform8aLog22,  dist = "negbin", link = "logit",data=dat)
summary(nbmod44)

nbform8aLog23=formula(SumCount~ y + cd + sc + bd + d + t + lat  | y + cd + sc + bd + d + temp)
nbmod45=zeroinfl(nbform8aLog23,  dist = "negbin", link = "logit",data=dat)
summary(nbmod45)

nbform8aLog24=formula(SumCount~ y + cd + sc + bd + d + t + lat  | y + cd + sc + bd + d + lat)
nbmod46=zeroinfl(nbform8aLog24,  dist = "negbin", link = "logit",data=dat)
summary(nbmod46)



lrtest(nbmod40, nbmod37)
lrtest(nbmod41, nbmod37)
lrtest(nbmod42, nbmod37)
lrtest(nbmod43, nbmod37)
lrtest(nbmod44, nbmod37)
lrtest(nbmod45, nbmod37)
lrtest(nbmod46, nbmod37)


AIC(nbmod40)
AIC(nbmod41)
AIC(nbmod42)
AIC(nbmod43)
AIC(nbmod44)
AIC(nbmod45)
AIC(nbmod46)


AIC(nbmod37)-AIC(nbmod40)
AIC(nbmod37)-AIC(nbmod41)
AIC(nbmod37)-AIC(nbmod42)
AIC(nbmod37)-AIC(nbmod43)
AIC(nbmod37)-AIC(nbmod44)
AIC(nbmod37)-AIC(nbmod45)
AIC(nbmod37)-AIC(nbmod46)




```

If we ignore the subtle difference between nbmod33 and nbmod42 and the effects of dropping y from the count process which I am not sure if a good idea what happens to the rest of the code.



```{r, echo=FALSE, warning=FALSE, error=FALSE}
nbbest=nbmod37
##nbform8aLog14=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + cd + sc + bd + d + lat + temp)

resids=residuals(nbbest, type="pearson")


```
```



```{r bestModPlots, echo=FALSE, warning=FALSE, error=FALSE}
#cbind(fitted(nbbest),dat$SumCount,fitted(nbbest)-dat$SumCount,resids)
par(mfrow=c(1,2))
plot(fitted(nbbest),resids, ylab="Pearsons Residuals", xlab="Fitted Values") #pearson resids vs fitted values
plot(dat$SumCount,fitted(nbbest), ylab="Fitted Values", xlab="Original Values") #pearson values vs original data
```

```{r, echo=FALSE, warning=FALSE, error=FALSE}
# A single figure
par(mfrow=c(1,3))
plot(dat$y,resids,xlab="Year",main="Residuals (nbbest)")
plot(dat$t,resids,xlab="Season",main="Residuals (nbbest)")
plot(dat$lat,resids,xlab="Latitude",main="Residuals (nbbest)")

par(mfrow=c(3,2))
plot(dat$wc,resids,xlab="Water Clarity",main="Residuals (nbbest)")
plot(dat$cd,resids,xlab="Current Direction",main="Residuals (nbbest)")
plot(dat$sc,resids,xlab="Substrate Composition",main="Residuals (nbbest)")
plot(dat$bd,resids,xlab="Biotic Diversity",main="Residuals (nbbest)")
plot(dat$d,resids,xlab="Depth",main="Residuals (nbbest)")
plot(dat$temp,resids,xlab="Temperature",main="Residuals (nbbest)")
```

```{r histNB, echo=FALSE, warning=FALSE, error=FALSE}
par(mfrow=c(1,2))
#full histogram, commented out in favor of the range limited figure below
hist(dat$SumCount,breaks=0:max(dat$SumCount),freq=T,right=TRUE,xlab='Aggregate Fish Counted', main='NB')  
d3=hist(predict(nbbest),breaks=0:max(dat$SumCount),plot=FALSE)
lines(seq(0.5,max(dat$SumCount),by=1),d3$counts, col="blue",type='l')   

# This is where Lew limited the range.
hist(dat$SumCount,breaks=0:max(dat$SumCount),freq=T,right=TRUE,xlab='Aggregate Fish Counted', main='NB',ylim=c(0,100))  
lines(seq(0.5,max(dat$SumCount),by=1),d3$counts, col="blue",type='l')      
lines(seq(0.5,max(dat$SumCount),by=1),d2$counts, col="red",type='l')      
```



```{r zinbFit, echo=FALSE, warning=FALSE, error=FALSE, fig.align='center'}
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
```




```{r, echo=FALSE, warning=FALSE, error=FALSE}
# final model form for FYI
#nbmod37
#nbform8aLog14=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + cd + sc + bd + d + lat + temp)

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
```



```{r bootStrap, cache=TRUE, echo=FALSE, warning=FALSE, error=FALSE}
#set up data objects and specify number of bootstrap replications
ptm <- proc.time()
boots=5000
#nbmod37
#nbform8aLog14=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + cd + sc + bd + d + lat + temp)
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
  
  #nbform8aLog14=formula(SumCount~ y + cd + sc + bd + d + t + lat  |y + cd + sc + bd + d + lat + temp)
  #make a function to compute the index and return either a valid index or NA conditional on if the model converges
  getindex=function(){
    #define and fit model for current replicate.  Use the "try" function so that can continue to run if model does not converge on a particular replication  
    f2 =formula(SumCount~ y + cd + sc + bd + d + t + lat  | y + cd + sc + bd + d + lat + temp)
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
save.image("\\\\CCFHR-S-1534090\\popdyn1\\Purcell\\RedSnapIndex_SEDAR41\\RedSnap_VideoIndex.RData")
```


```{r, echo=FALSE, warning=FALSE, error=FALSE, fig.align='center'}
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
lines(nomcpue.std$nomcpue.std, type='l', col='blue', lty=2,lwd=2)
#approxgam=c(.37,.11,.21)
#lines(approxgam/mean(approxgam),lwd=2,col='blue')
legend('topright',legend=c('Standardized Index','Bootstrap CI','Nominal'),lty=c(1,3,1),col=c('red','black','blue'),lwd=c(2,1,2), bty="n")
```


```{r sessionInfo, echo=FALSE, warning=FALSE, error=FALSE, results='hide'}
# produce R session info
sessionInfo()
save.image("\\\\CCFHR-S-1534090\\popdyn1\\Purcell\\RedSnapIndex_SEDAR41\\RedSnap_VideoIndex.RData")
```