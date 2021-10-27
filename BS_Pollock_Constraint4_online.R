#Script for estimating pollock habitat constraints metrics and relative figures

mydir<-'XXX'#Change this as appropirate
setwd(mydir) 
library(mgcv)
library(maps)
library(mapdata)
library(spacetime)
library(fields)
library(date)
library(colorRamps)  
library(itsadug)    
library(RColorBrewer)
library(viridis)
source('distance.function.R')

#DEPTH: import data from:
#http://www.ngdc.noaa.gov/mgg/gdas/gd_designagrid.html
bathy.dat<-read.table('/Users/lorenzociannelli/Documents/MyDocuments/Funding/NSF/Geoscience/data/BeringDepth.txt',sep='')
names(bathy.dat)<-c('lon','lat','depth')
bathy.dat$depth[bathy.dat$depth>0]<-NA#Avoid points above water
head(bathy.dat)
bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]

#Import pollock eggs, larvae, juvenile and adult data
subset.egg<-read.table("pollock_eggs.csv")
subset.larvae<-read.table("pollock_larvae.csv")
pk_adults_catch<-read.table("pollock_juvenile_adult_catch.csv")
pk_adults_length<-read.table("pollock_juvenile_adult_length.csv")

#Distribute size intervals based on cumulative distribution of catches as a function of size
length.seq<-seq(min(pk_adults_length$LENGTH),max(pk_adults_length$LENGTH),length=100)
length.cum<-length.seq*NA
for(i in 1:length(length.seq)){
length.cum[i]<-sum(pk_adults_length$cpuelt[pk_adults_length$LENGTH<=length.seq[i]])	
}

#Size bins = XX% of total numerical cacth
XX=0.125
perc.bin=seq(0,1,by=XX)
size.bin<-1:length(perc.bin)*NA
size.bin[1]<-min(pk_adults_length$LENGTH)
size.bin[length(size.bin)]<-max(pk_adults_length$LENGTH)
for(i in 2:(length(perc.bin)-1)){
size.bin[i]<-length.seq[(length.cum/sum(pk_adults_length$cpuelt))>perc.bin[i]][1]
}

plot(length.seq,length.cum/sum(pk_adults_length$cpuelt)*100,ylab="Cumulative biomass (%)")
abline(v=size.bin[2:length(size.bin)])
abline(h=(perc.bin*100)[2:length(perc.bin)])


#Build vector of sample sizes for each size bin and stage
sample.sizes<-1:(length(size.bin)+1)*NA
#Add eggs and larval sample sizes
sample.sizes[1:2]<-c(nrow(subset.egg[subset.egg$LARVALCATCHPER10M2>0,]),
                     nrow(subset.larvae[subset.larvae$LARVALCATCHPER10M2>0,]))

#CONSTRAINT ANALYSES: egg
#Grid for habitat extent calculations: note finer scale than that used for adults
nlat=40
nlon=60
latd=seq(min(c(subset.egg$lat,subset.larvae$lat,pk_adults_length$LAT)),
         max(c(subset.egg$lat,pk_adults_length$LAT)),length.out=nlat)
lond=seq(min(c(subset.egg$lon,subset.larvae$lat,pk_adults_length$LON)),
         max(c(subset.egg$lon,pk_adults_length$LON)),length.out=nlon)
grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

#Calculate distance of each grid point to closest 'observation': determining what the 'observation' is a critical juncture I am using 'observations' as location sampled and with positive catches in all data set, regardless of stage being considered
for(i in 1:nrow(grid.extent)){
	dist<-distance.function(grid.extent$lat[i],grid.extent$lon[i],
	                        c(subset.egg$lat,subset.larvae$lat,pk_adults_length$LAT),
	                        c(subset.egg$lon,subset.larvae$lat,pk_adults_length$LON))
	grid.extent$dist_outer[i]<-min(dist)
	dist<-distance.function(grid.extent$lat[i],grid.extent$lon[i],subset.egg$lat,subset.egg$lon)
	grid.extent$dist_inner[i]<-min(dist)
}


#Assign a within sample year and doy to the grid data frame and a vector to fill in cum biomass and habitat
grid.extent$year<-2005
grid.extent$doy<-median(subset.egg$doy)
grid.extent$cum.b<-NA
grid.extent$cum.h<-NA

#Define a cut off of cumulative biomass for the calculation of habitat extent
cut.off_b<-0.75 #(smallest habitat that contains 'cut.off_b%' of predicted biomass)

#Define vectors for output
rsq.ad<-NA*1:(length(size.bin)-1);all.var<-rsq.ad;all.h<-rsq.ad;

#Estimate habitat constraint and plot results on a map
quartz(height=4,width=10)
par(mfrow=c(2,5),mai=c(0.1,0.25,0.0,0.3),omi=c(0.25,0.1,0.50,0.2))

gam.pkegg<-gam(log(LARVALCATCHPER10M2)~factor(year)+s(lon,lat)+s(doy),data=subset.egg[subset.egg$LARVALCATCHPER10M2>0,])
summary(gam.pkegg)
#R-sq.(adj) =  0.352   Deviance explained = 36.9%
#GCV = 3.4621  Scale est. = 3.3699    n = 1955
gam.pkegg.base<-gam(log(LARVALCATCHPER10M2)~factor(year)+s(doy),data=subset.egg[subset.egg$LARVALCATCHPER10M2>0,])
summary(gam.pkegg.base)
#R-sq.(adj) =  0.124   Deviance explained = 13.5%
#GCV = 4.6141  Scale est. = 4.5538    n = 1955
#
gam.pkegg.pres<-gam(I(1*(LARVALCATCHPER10M2>0))~factor(year)+s(lon,lat)+s(doy),data=subset.egg,family=binomial)
summary(gam.pkegg.pres)
#R-sq.(adj) =  0.227   Deviance explained = 23.1%
#UBRE = -0.19901  Scale est. = 1         n = 2429
#

#MSE ratio
var_ratio.egg<-(summary(gam.pkegg.base)$scale-summary(gam.pkegg)$scale)/summary(gam.pkegg.base)$scale
rsq.egg<-summary(gam.pkegg)$r.sq

#Calculate habitat extent based on predictions over the standardized grid from habitat GAM
grid.extent$pred_ab<-exp(predict(gam.pkegg,newdata=grid.extent))
grid.extent$pred_ps<-predict(gam.pkegg.pres,newdata=grid.extent,type='response')
grid.extent$pred<-grid.extent$pred_ab*grid.extent$pred_ps
grid.extent$pred[grid.extent$dist_inner>30000]<-0
grid.extent$pred[grid.extent$dist_outer>30000]<-NA

#Plot predictions and observations
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),ncol=length(lond),byrow=T)),
           col=viridis(100),ylab=expression(paste("Latitude ("^0,'N)')),
           xlab=expression(paste("Longitude ("^0,'E)')),legend.width=2,
           xlim=range(pk_adults_length$LON),ylim=c(53,max(pk_adults_length$LAT)),axes=F)
axis(3,labels=T)
axis(2,labels=T)
symbols(subset.egg$lon[subset.egg$LARVALCATCHPER10M2>0],
        subset.egg$lat[subset.egg$LARVALCATCHPER10M2>0],
        circles=subset.egg$LARVALCATCHPER10M2[subset.egg$LARVALCATCHPER10M2>0],
        inches=0.1,bg=alpha('white',f=0.2),fg=alpha('black',f=0.1),add=T)

#Resume code  of habitat constraint calculation
grid.extent<-grid.extent[order(grid.extent$pred,decreasing=T),]
tmp1<-grid.extent[grid.extent$dist_outer<30000,]
tmp1$cum.b<-cumsum(tmp1$pred)/sum(tmp1$pred)
tmp1$cum.h<-1:nrow(tmp1)/nrow(tmp1)

all.h.egg<-max(tmp1$cum.h[tmp1$cum.b<=cut.off_b])
HC.egg<-var_ratio.egg*(1-all.h.egg)

#Add habitat extent
points(tmp1$lon[tmp1$cum.b<=cut.off_b],
       tmp1$lat[tmp1$cum.b<=cut.off_b],
       col=alpha('white',f=0.4),pch='+')
map("worldHires",fill=T,col="lightblue4",add=T)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200),labcex=0.4,col='black',add=T)
text(-174.5,54,labels='Eggs',cex=1.1)
text(-161.5,61.5,labels=expression(paste("# 10m"^-2)),cex=1.1,col="white")

#CONSTRAINT ANALYSES: LARVAE all sizes

#Grid for habitat extent calculations: note finer scale than that used for adults
nlat=40
nlon=60
latd=seq(min(c(subset.egg$lat,subset.larvae$lat,pk_adults_length$LAT)),
         max(c(subset.egg$lat,pk_adults_length$LAT)),length.out=nlat)
lond=seq(min(c(subset.egg$lon,subset.larvae$lat,pk_adults_length$LON)),
         max(c(subset.egg$lon,pk_adults_length$LON)),length.out=nlon)
grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')


#Calculate distance of each grid point to closest 'observation': determining what the 'observation' is a critical juncture I am using 'observations' as location sampled and with positive catches in all data set, regardless of stage being considered
for(i in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[i],grid.extent$lon[i],
                          c(subset.egg$lat,subset.larvae$lat,pk_adults_length$LAT),
                          c(subset.egg$lon,subset.larvae$lon,pk_adults_length$LON))
  grid.extent$dist_outer[i]<-min(dist)
  dist<-distance.function(grid.extent$lat[i],grid.extent$lon[i],subset.larvae$lat,subset.larvae$lon)
  grid.extent$dist_inner[i]<-min(dist)
}


#Assign a within sample year and doy to the grid data frame and a vector to fill in cum biomass and habitat
grid.extent$year<-2005
grid.extent$doy<-median(subset.larvae$doy)
grid.extent$cum.b<-NA
grid.extent$cum.h<-NA

gam.pklv<-gam(log(LARVALCATCHPER10M2)~factor(year)+
                s(lon,lat)+s(doy),data=subset.larvae[subset.larvae$LARVALCATCHPER10M2>0,])
summary(gam.pklv)
#R-sq.(adj) =  0.404   Deviance explained = 42.2%
#GCV = 2.6049  Scale est. = 2.5273    n = 1809
#
gam.pklv.base<-gam(log(LARVALCATCHPER10M2)~factor(year)+s(doy),data=subset.larvae[subset.larvae$LARVALCATCHPER10M2>0,])
summary(gam.pklv.base)
#R-sq.(adj) =  0.235   Deviance explained = 24.6%
#GCV = 3.2939  Scale est. = 3.2457    n = 1809
#
gam.pklv.pres<-gam(I(1*(LARVALCATCHPER10M2>0))~factor(year)+s(lon,lat)+s(doy),data=subset.larvae,family=binomial)
summary(gam.pklv.pres)
#R-sq.(adj) =  0.312   Deviance explained = 28.1%
#UBRE = -0.13929  Scale est. = 1         n = 2434
#

#MSE ratio
var_ratio.lv<-(summary(gam.pklv.base)$scale-summary(gam.pklv)$scale)/summary(gam.pklv.base)$scale
rsq.lv<-summary(gam.pklv)$r.sq

#Calculate habitat extent based on predictions over the standardized grid from habitat GAM
grid.extent$pred_ab<-exp(predict(gam.pklv,newdata=grid.extent))
grid.extent$pred_ps<-predict(gam.pklv.pres,newdata=grid.extent,type='response')
grid.extent$pred<-grid.extent$pred_ab*grid.extent$pred_ps
grid.extent$pred[grid.extent$dist_inner>30000]<-0
grid.extent$pred[grid.extent$dist_outer>30000]<-NA

#Plot predictions and observations
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),ncol=length(lond),byrow=T)),
           col=viridis(100),ylab="",xlab="",legend.width=2,
           xlim=range(pk_adults_length$LON),
           ylim=c(53,max(pk_adults_length$LAT)),axes=F)
symbols(subset.larvae$lon[subset.larvae$LARVALCATCHPER10M2>0],
        subset.larvae$lat[subset.larvae$LARVALCATCHPER10M2>0],
        circles=subset.larvae$LARVALCATCHPER10M2[subset.larvae$LARVALCATCHPER10M2>0],
        inches=0.1,bg=alpha('white',f=0.2),fg=alpha('black',f=0.1),add=T)


#Resume code  of habitat constraint calculation
grid.extent<-grid.extent[order(grid.extent$pred,decreasing=T),]
tmp1<-grid.extent[grid.extent$dist_outer<30000,]
tmp1$cum.b<-cumsum(tmp1$pred)/sum(tmp1$pred)
tmp1$cum.h<-1:nrow(tmp1)/nrow(tmp1)

all.h.lv<-max(tmp1$cum.h[tmp1$cum.b<=cut.off_b])
HC.lv<-var_ratio.lv*(1-all.h.lv)


#Add habitat extent
points(tmp1$lon[tmp1$cum.b<=cut.off_b],
       tmp1$lat[tmp1$cum.b<=cut.off_b],
       col=alpha('white',f=0.4),pch='+')
map("worldHires",fill=T,col="lightblue4",add=T)
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200),labcex=0.4,col='black',add=T)
text(-174.5,54,labels='Larvae',cex=1.1)
text(-161.5,61.5,labels=expression(paste("# 10m"^-2)),cex=1.1,col="white")

##### CONSTRAINT ANALYSES: Juveniles and Adults
#Grid for habitat extent calculations
nlat=20
nlon=30
latd=seq(min(c(subset.egg$lat,subset.larvae$lat,pk_adults_length$LAT)),
         max(c(subset.egg$lat,pk_adults_length$LAT)),length.out=nlat)
lond=seq(min(c(subset.egg$lon,subset.larvae$lat,pk_adults_length$LON)),
         max(c(subset.egg$lon,pk_adults_length$LON)),length.out=nlon)
grid.extent<-expand.grid(lond,latd)
names(grid.extent)<-c('lon','lat')

#Calculate distance of each grid point to closest 'observation': determining what the 'observation' is a critical juncture I am using 'observations' as location sampled and with positive catches in the data set where the stage considered came from (groundfish survey). Only stations close to positive catches are used for the observation.
for(i in 1:nrow(grid.extent)){
  dist<-distance.function(grid.extent$lat[i],grid.extent$lon[i],
                          c(subset.egg$lat,subset.larvae$lat,pk_adults_length$LAT),
                          c(subset.egg$lon,subset.larvae$lon,pk_adults_length$LON))  
  grid.extent$dist_outer[i]<-min(dist)
  dist<-distance.function(grid.extent$lat[i],grid.extent$lon[i],pk_adults_length$LAT,pk_adults_length$LON)
  grid.extent$dist_inner[i]<-min(dist)
}

#Assign a within sample year and doy to the grid data frame and a vector to fill in cum biomass and habitat
grid.extent$year<-median(pk_adults_length$YEAR)
grid.extent$doy<-median(pk_adults_length$doy)
grid.extent$cum.b<-NA
grid.extent$cum.h<-NA


for(i in 1:(length(size.bin)-1)){
tmp<-pk_adults_length[pk_adults_length$LENGTH<=size.bin[i+1]&pk_adults_length$LENGTH>size.bin[i],]
pk.l<-aggregate(tmp[,c("YEAR","LAT","LON","doy")],list(tmp$HAULJOIN),mean)
names(pk.l)<-c('HAULJOIN','year','lat','lon','doy')
pk.l$cpuelt<-aggregate(tmp$cpuelt,list(tmp$HAULJOIN),sum)$x

#Habitat GAM 
pk.gaml1<-gam(log(cpuelt)~factor(year)+s(lon,lat)+s(doy),data=pk.l)
#Base GAM
pk.gaml1_mean<-gam(log(cpuelt)~factor(year)+s(doy),data=pk.l)
#Pres/abs GAM
idx<-(1:nrow(pk_adults_catch))[!c(pk_adults_catch$HAULJOIN%in%pk.l$HAULJOIN)]
tmp1<-pk_adults_catch[idx,c('HAULJOIN','YEAR','START_LATITUDE','START_LONGITUDE','doy','NUMBER_FISH','dist')]
tmp1$NUMBER_FISH<-0
names(tmp1)<-c('HAULJOIN','year','lat','lon','doy','cpuelt','dist')
pk.pres<-rbind(pk.l,tmp1[tmp1$dist<30000,c('HAULJOIN','year','lat','lon','doy','cpuelt')])
gam.pres<-gam(I(1*(cpuelt>0))~factor(year)+s(lon,lat)+s(doy),data=pk.pres,family=binomial)


#MSE ratio and sample size
var_ratio<-(summary(pk.gaml1_mean)$scale-summary(pk.gaml1)$scale)/summary(pk.gaml1_mean)$scale
rsq<-summary(pk.gaml1)$r.sq
rsq.ad[i]<-rsq
all.var[i]<-var_ratio
sample.sizes[i+2]<-nrow(pk.l)

#Calculate habitat extent based on predictions over the standardized grid from habitat GAM
grid.extent$pred_ab<-exp(predict(pk.gaml1,newdata=grid.extent))
grid.extent$pred_ps<-predict(gam.pres,newdata=grid.extent,type='response')
grid.extent$pred<-grid.extent$pred_ab*grid.extent$pred_ps
grid.extent$pred[grid.extent$dist_inner>30000]<-0
grid.extent$pred[grid.extent$dist_outer>30000]<-NA

#Plot predictions and observations
image.plot(lond,latd,t(matrix(grid.extent$pred,nrow=length(latd),ncol=length(lond),byrow=T)),
           col=viridis(100),ylab="",xlab="",
           xlim=range(pk_adults_length$LON),legend.width=2,
           ylim=c(53,max(pk_adults_length$LAT)),axes=F)
symbols(pk.l$lon,pk.l$lat,circles=pk.l$cpuelt,
        inches=0.1,bg=alpha('white',f=0.1),
        fg=alpha('black',f=0.1),add=T)

#Resume: habitat extent
grid.extent2<-grid.extent[order(grid.extent$pred,decreasing=T),]
tmp1<-grid.extent2[grid.extent2$dist_outer<30000,]
tmp1$cum.b<-cumsum(tmp1$pred)/sum(tmp1$pred)
tmp1$cum.h<-1:nrow(tmp1)/nrow(tmp1)


all.h[i]<-max(tmp1$cum.h[tmp1$cum.b<=cut.off_b])
HC<-all.var*(1-all.h)

#Add habitat extent to predictions
text(-174.5,54,labels=paste(round(size.bin[i],0),"-",round(size.bin[i+1],0)),cex=1.1)
points(tmp1$lon[tmp1$cum.b<=cut.off_b],
       tmp1$lat[tmp1$cum.b<=cut.off_b],
       col=alpha('white',f=0.4),pch='+')
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-c(200),labcex=0.4,col='black',add=T)
map("worldHires",fill=T,col="lightblue4",add=T)
text(-161,61.5,labels=expression(paste("# ha"^-1)),cex=1.1,col="white")
}

a<-dev.cur()
dev.copy(jpeg,paste0(mydir,"/figures/",'pk_maps_lowres.jpg'),
         height=4,width=10,res=100,units='in')
dev.off()
dev.set(a)
dev.copy(jpeg,paste0(mydir,"/figures/",'pk_maps_hires.jpg'),
         height=4,width=10,res=600,units='in')
dev.off()


#Barplots of metrics of constraint
HC.all<-c(HC.egg,HC.lv,HC)
Rsq.all<-c(rsq.egg,rsq.lv,rsq.ad)
allages<-paste(round(size.bin[1:(length(size.bin)-1)],0),round(size.bin[2:length(size.bin)],0),sep="-")
allages<-c('Eggs','Lv',allages)

quartz(width=10,height=7)
par(mfcol=c(2,2))

barplot(Rsq.all,names.arg=allages,
        ylab=expression('R'^2),xlab='',las=2,
        main="Explained spatio-temporal variance",ylim=c(0,1))
box()
bar_data<-barplot(Rsq.all,plot=FALSE)
text(bar_data,0.05+Rsq.all,as.character(sample.sizes))

barplot(c(var_ratio.egg,var_ratio.lv,all.var),
        names.arg=allages,ylab='MSE ratio',
        xlab='Stage or size (mm)',main="Geographic consistency",ylim=c(0,1))
box()

barplot(c(all.h.egg,all.h.lv,all.h),
        names.arg=allages,ylab='Fraction of 75% total B',
        las=2,xlab='',main="Habitat Extent",ylim=c(0,1))
box()

barplot(HC.all,names.arg=allages,
        ylab='MSE ratio x (1- HE)',
        xlab='Stage or size (mm)',main="Habitat Constraint",
        ylim=c(0,1))
box()

a<-dev.cur()
dev.copy(jpeg,paste0(mydir,"/figures/",'pk_metrics_lowres.jpg'),
         height=7,width=10,res=100,units='in')
dev.off()
dev.set(a)
dev.copy(jpeg,paste0(mydir,"/figures/",'pk_metrics_hires.jpg'),
         height=7,width=10,res=600,units='in')
dev.off()

pk_metrics<-data.frame(lower.size=c(1,1,size.bin[1:(length(size.bin)-1)]),upper.size=c(1,2,size.bin[2:length(size.bin)]))
pk_metrics$Rsq<-Rsq.all
pk_metrics$var_ratio<-c(var_ratio.egg,var_ratio.lv,all.var)
pk_metrics$h_ext<-c(all.h.egg,all.h.lv,all.h)
pk_metrics$HC<-HC.all
pk_metrics$sample.size<-sample.sizes
write.table(pk_metrics,'pk_metrics.txt')
