#July 24th 2021

#Try with sgstat package following example of Roger Bivand (link below) for conditional gaussian random fields
#https://www.aspexit.com/simulating-spatial-datasets-with-known-spatial-variability/
#http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/
#Aug 1st, 2021: Add a fourth scenario, representing a conditional distribution, with center of abundance shifting position from year to year
mydir<-'xxx'
setwd(mydir) #Change this as appropriate
library(gstat)
library(sp)
library(mgcv)
library(viridis)
library(plotfunctions)
library(fields)

#SIMULATE DATA AND VISUALIZE SPATIAL FIELDS

#####1. CONDITIONAL DISTRIBUTION, highly constrained: 
#increase X and Y coeff, reduce nugget and range

# create spatial structure
xy <- expand.grid(1:20, 1:20)
names(xy) <- c("x","y")

# define the gstat object (spatial model)
Beta<-c(0,0.3,0.3)#First number is average magnitude
#Second is the linear correlation in X
#Third is linear corr in Y
Psill=1  ## Partial sill = Magnitude of variation
Range=0.1  ## Maximal distance of autocorrelation
Nugget=0.1  ## Small-scale variations
g.dummy <- gstat(formula=z~1+x+y, locations=~x+y, 
                 dummy=T, beta=Beta,
                 model=vgm(psill=Psill, range=Range, 
                           nugget=Nugget,model='Exp'), 
                 nmax=20)
# make n.rep simulations based on the stat object
n.rep=20
rm(yy_high_pred)
yy_high_pred_raw<- predict(g.dummy, newdata=xy, nsim=n.rep)

yy_high_pred<-yy_high_pred_raw

#Make data frame of simulated data for constraint analyses
sim.data_high<-cbind(xy,unlist(yy_high_pred[,3:(n.rep+2)]),year=rep(1:10,each=400))
names(sim.data_high)<-c("x","y","cpue","year")
head(sim.data_high)
tapply(sim.data_high$cpue,sim.data_high$year,mean)

# show first 10 simulations:
quartz(width=7,height=5)
par(mfrow=c(3,4),omi=c(0.18,0.18,0.1,0.1),mai=c(0.1,0.1,0,0))
for(i in 1:10){
        if(i==9){
                image(c(1:20),c(1:20),
                      matrix(yy_high_pred[,i+2],nrow=20,ncol=20,byrow=F),
                      xlab="X",ylab="Y",col=viridis(100),
                      zlim=range(yy_high_pred[,3:12]))
                text(5,19,paste("Sim",i),col="white",font=4)}else{
                              image(c(1:20),c(1:20),
                                    matrix(yy_high_pred[,i+2],
                                        nrow=20,ncol=20,byrow=F),
                                    xlab="",ylab="",col=viridis(100),axes=F,
                                    zlim=range(yy_high_pred[,3:12]))
                        text(5,19,paste("Sim",i),col="white",font=4)}}
par(oma=c( 0.9,0,0,24))
image.plot(legend.only=T,zlim=range(yy_high_pred[,3:12]),
           col=viridis(100),legend.width=5)

#dev.copy(jpeg,paste0(mydir,"/figures/",'sim_high.jpg'),height=5,width=7,res=300,units='in')
# dev.off()

######2. CONDITIONAL DISTRIBUTION: yearly shifting center of gravity
#keep X and Y coeff, nugget and range. 
#Randomly shift direction of X and Y coeff
n.rep<-20
# create spatial structure
rm(yy_high_pred_year_raw,xy)
xy <- expand.grid(1:20, 1:20)
names(xy) <- c("x","y")
yy_high_pred_year_raw<-xy

Beta_all<-matrix(c(
        c(0,0.3,0.3),
        c(0,-0.3,0.3),
        c(0,0.3,-0.3),
        c(0,-0.3,-0.3)),
        nrow=3,ncol=4)
Psill=1  ## Partial sill = Magnitude of variation
Range=0.1  ## Maximal distance of autocorrelation
Nugget=0.1  ## Small-scale variations

for(i in 1:n.rep){
# define the gstat object (spatial model)
Beta<-Beta_all[,sample(1:4,1)]
#Beta<-c(0,3,3) 
g.dummy <- gstat(formula=z~1+x+y, locations=~x+y, 
                 dummy=T, beta=Beta,
                 model=vgm(psill=Psill, range=Range, 
                           nugget=Nugget,model='Exp'), 
                 nmax=20)
# make n.rep simulations based on the stat object
tmp1<-predict(g.dummy, newdata=xy, nsim=n.rep)
yy_high_pred_year_raw<- cbind(yy_high_pred_year_raw,tmp1[,sample(3:(n.rep+2),1)])}

yy_high_pred_year<-yy_high_pred_year_raw

#Make data frame of simulated data for constraint analyses
sim.data_high_year<-cbind(xy,unlist(yy_high_pred_year[,3:(n.rep+2)]),year=rep(1:10,each=400))
names(sim.data_high_year)<-c("x","y","cpue","year")
head(sim.data_high_year)
tapply(sim.data_high_year$cpue,sim.data_high_year$year,mean)

# show first 10 simulations:
quartz(width=7,height=5)
par(mfrow=c(3,4),omi=c(0.18,0.18,0.1,0.1),
    mai=c(0.1,0.1,0,0))
for(i in 1:10){
        if(i==9){
                image(c(1:20),c(1:20),
                      matrix(yy_high_pred_year[,i+2],
                             nrow=20,ncol=20,byrow=F),
                      xlab="X",ylab="Y",col=viridis(100),
                      zlim=range(yy_high_pred_year[,3:12]))
                text(5,19,paste("Sim",i),col="white",font=4)}else{
                        image(c(1:20),c(1:20),
                              matrix(yy_high_pred_year[,i+2],
                                     nrow=20,ncol=20,byrow=F),
                              xlab="",ylab="",col=viridis(100),axes=F,
                              zlim=range(yy_high_pred_year[,3:12]))
                        text(5,19,paste("Sim",i),col="white",font=4)}}
par(oma=c( 0.9,0,0,24))
image.plot(legend.only=T,zlim=range(yy_high_pred_year[,3:12]),
           col=viridis(100),legend.width=5)

#dev.copy(jpeg,paste0(mydir,"/figures/",'sim_high_year.jpg'),height=5,width=7,res=300,units='in')
# dev.off()



######3. CONDITIONAL DISTRIBUTION: Reduce spatial constraint
#Reduce X and Y coeff, increase sill, nugget and range
# create structure
xy <- expand.grid(1:20, 1:20)
names(xy) <- c("x","y")
# define the gstat object (spatial model)
Beta<-c(0,0.15,0.15)#First number is average magnitude
#Second is the linear correlation in X
#Third is linear corr in Y
Psill=1  ## Partial sill = Magnitude of variation
Range=0.5  ## Maximal distance of autocorrelation
Nugget=0.5  ## Small-scale variations
g.dummy <- gstat(formula=z~1+x+y, locations=~x+y, 
                 dummy=T, beta=Beta,
                 model=vgm(psill=Psill, range=Range, 
                           nugget=Nugget,model='Exp'), 
                 nmax=20)

# make n.rep simulations based on the stat object
n.rep=20
rm(yy_med_pred)
yy_med_pred_raw<- predict(g.dummy, newdata=xy, nsim=n.rep)

yy_med_pred<-yy_med_pred_raw

#Make data frame of simulated data for constraint analyses
sim.data_med<-cbind(xy,unlist(yy_med_pred[,3:(n.rep+2)]),year=rep(1:10,each=400))
names(sim.data_med)<-c("x","y","cpue","year")
head(sim.data_med)
tapply(sim.data_med$cpue,sim.data_med$year,mean)


# show first 10 simulations:
quartz(width=7,height=5)
par(mfrow=c(3,4),omi=c(0.18,0.18,0.1,0.1),
    mai=c(0.1,0.1,0,0))
for(i in 1:10){
        if(i==9){
                image(c(1:20),c(1:20),
                      matrix(yy_med_pred[,i+2],
                             nrow=20,ncol=20,byrow=F),
                      xlab="X",ylab="Y",col=viridis(100),
                      zlim=range(yy_med_pred[,3:12]))
                text(5,19,paste("Sim",i),col="white",font=4)}else{
                        image(c(1:20),c(1:20),
                              matrix(yy_med_pred[,i+2],
                                     nrow=20,ncol=20,byrow=F),
                              xlab="",ylab="",col=viridis(100),axes=F,
                              zlim=range(yy_med_pred[,3:12]))
                        text(5,19,paste("Sim",i),col="white",font=4)}}
par(oma=c( 0.9,0,0,24))
image.plot(legend.only=T,zlim=range(yy_med_pred[,3:12]),
           col=viridis(100),legend.width=5)

#dev.copy(jpeg,paste0(mydir,"/figures/",'sim_med.jpg'),height=5,width=7,res=300,units='in')
#dev.off()




#4. Unconditional RANDOM FIELDS: remove coeff of X and Y, increase nugget and range
# create structure
xy <- expand.grid(1:20, 1:20)
names(xy) <- c("x","y")
# define the gstat object (spatial model)
Beta<-c(0)#First number is average magnitude
Psill=1  ## Partial sill = Magnitude of variation
Range=1  ## Maximal distance of autocorrelation
Nugget=1  ## Small-scale variations
g.dummy <- gstat(formula=z~1, locations=~x+y, 
                 dummy=T, beta=Beta,
                 model=vgm(psill=Psill, range=Range, 
                           nugget=Nugget,model='Exp'), 
                 nmax=20)

# make n.rep simulations based on the stat object
n.rep=20
rm(yy_low_pred)
yy_low_pred_raw<- predict(g.dummy, newdata=xy, nsim=n.rep)

yy_low_pred<-yy_low_pred_raw

#Make data frame of simulated data for constraint analyses
sim.data_low<-cbind(xy,unlist(yy_low_pred[,3:(n.rep+2)]),year=rep(1:10,each=400))
names(sim.data_low)<-c("x","y","cpue","year")
head(sim.data_low)
tapply(sim.data_low$cpue,sim.data_low$year,mean)

# show all 10 simulations:
quartz(width=7,height=5)
par(mfrow=c(3,4),omi=c(0.18,0.18,0.1,0.1),
    mai=c(0.1,0.1,0,0))
for(i in 1:10){
        if(i==9){
                image(c(1:20),c(1:20),
                      matrix(yy_low_pred[,i+2],
                             nrow=20,ncol=20,byrow=F),
                      xlab="X",ylab="Y",col=viridis(100),
                      zlim=range(yy_low_pred[,3:12]))
                text(5,19,paste("Sim",i),col="white",font=4)}else{
                        image(c(1:20),c(1:20),
                              matrix(yy_low_pred[,i+2],
                                     nrow=20,ncol=20,byrow=F),
                              xlab="",ylab="",col=viridis(100),axes=F,
                              zlim=range(yy_low_pred[,3:12]))
                        text(5,19,paste("Sim",i),col="white",font=4)}}
par(oma=c( 0.9,0,0,24))
image.plot(legend.only=T,zlim=range(yy_low_pred[,3:12]),
           col=viridis(100),legend.width=5)

#dev.copy(jpeg,paste0(mydir,"/figures/",'sim_low.jpg'),height=5,width=7,res=300,units='in')
#dev.off()

#START CONSTRAINT ANALYSES
#1. Conditional distribution 
gam.sim_high<-gam(cpue~factor(year)+s(x,y),data=sim.data_high)
summary(gam.sim_high)
#R-sq.(adj) =  0.845   Deviance explained = 84.6%
#GCV = 1.0962  Scale est. = 1.0929    n = 4000
gam.sim.base_high<-gam(cpue~factor(year),data=sim.data_high)
summary(gam.sim.base_high)
#R-sq.(adj) =  -0.00192   Deviance explained = 0.0333%
#GCV =  7.096  Scale est. = 7.0782    n = 4000
#MSE ratio
var_ratio.sim_high<-(summary(gam.sim.base_high)$scale-summary(gam.sim_high)$scale)/summary(gam.sim.base_high)$scale
rsq.sim_high<-summary(gam.sim_high)$r.sq
#Make predictions over the standardized grid from habitat GAM
cut.off_b<-0.75
grid.extent_high<-cbind(xy,1)
names(grid.extent_high)<-c("x","y","year")
grid.extent_high$pred<-exp(predict(gam.sim_high,newdata=grid.extent_high))

#2. Yearly variable distribution
#Start modeling 
gam.sim_high_year<-gam(cpue~factor(year)+s(x,y),data=sim.data_high_year)
summary(gam.sim_high_year)
#R-sq.(adj) =  0.0209   Deviance explained = 2.12%
#GCV = 375.84  Scale est. = 375.72    n = 40000
gam.sim.base_high_year<-gam(cpue~factor(year),data=sim.data_high_year)
summary(gam.sim.base_high_year)
#R-sq.(adj) =  0.00221   Deviance explained = 0.244%
#GCV = 382.99  Scale est. = 382.89    n = 40000
#MSE ratio
var_ratio.sim_high_year<-(summary(gam.sim.base_high_year)$scale-summary(gam.sim_high_year)$scale)/summary(gam.sim.base_high_year)$scale
rsq.sim_high_year<-summary(gam.sim_high_year)$r.sq
#Make predictions over the standardized grid from habitat GAM
cut.off_b<-0.75
grid.extent_high_year<-cbind(xy,3)
names(grid.extent_high_year)<-c("x","y","year")
grid.extent_high_year$pred<-exp(predict(gam.sim_high_year,newdata=grid.extent_high_year))

#3. Moderate conditional distribution
#Start modeling 
gam.sim_med<-gam(cpue~factor(year)+s(x,y),data=sim.data_med)
summary(gam.sim_med)
#R-sq.(adj) =  0.497   Deviance explained =   50%
#GCV = 1.5095  Scale est. = 1.5019    n = 4000
gam.sim.base_med<-gam(cpue~factor(year),data=sim.data_med)
summary(gam.sim.base_med)
#R-sq.(adj) =  5.78e-05   Deviance explained = 0.231%
#GCV = 2.9947  Scale est. = 2.9872    n = 4000
#MSE ratio
var_ratio.sim_med<-(summary(gam.sim.base_med)$scale-summary(gam.sim_med)$scale)/summary(gam.sim.base_med)$scale
rsq.sim_med<-summary(gam.sim_med)$r.sq
#Make predictions over the standardized grid from habitat GAM
cut.off_b<-0.75
grid.extent_med<-cbind(xy,1)
names(grid.extent_med)<-c("x","y","year")
grid.extent_med$pred<-exp(predict(gam.sim_med,newdata=grid.extent_med))
head(grid.extent_med)

#4. Random fields
#Start modeling 
gam.sim_random<-gam(cpue~factor(year)+s(x,y),data=sim.data_low)
summary(gam.sim_random)
#R-sq.(adj) =  0.0149   Deviance explained = 2.23%
#GCV = 1.9886  Scale est. = 1.9733    n = 4000
gam.sim.base_random<-gam(cpue~factor(year),data=sim.data_low)
summary(gam.sim.base_random)
#R-sq.(adj) =  0.00259   Deviance explained = 0.484%
#GCV =  2.003  Scale est. = 1.998     n = 4000
#MSE ratio
var_ratio.sim_random<-(summary(gam.sim.base_random)$scale-summary(gam.sim_random)$scale)/summary(gam.sim.base_random)$scale
rsq.sim_random<-summary(gam.sim_random)$r.sq
#Make predictions over the standardized grid from habitat GAM
cut.off_b<-0.75
grid.extent_random<-cbind(xy,1)
names(grid.extent_random)<-c("x","y","year")
grid.extent_random$pred<-exp(predict(gam.sim_random,newdata=grid.extent_random))
head(grid.extent_random)

#CONTINUE CONSTRAINT ANALYSES 
########Make Figures
quartz(width=11,height=7)
par(mfcol=c(2,4),omi=c(0.2,0.2,0.3,0.7),mai=c(0.5,0.5,0.2,0.0))
#1. High constraint
latd<-1:20
lond<-1:20
grid.extent_high$pred_scaled<-(grid.extent_high$pred-min(grid.extent_high$pred))/
  (max(grid.extent_high$pred)-min(grid.extent_high$pred))
range(grid.extent_high$pred_scaled)
  
image(lond,latd,
           t(matrix(grid.extent_high$pred_scaled,nrow=length(latd),ncol=length(lond),byrow=T)),
           col=viridis(100),ylab="Y",xlab="X",main="High constraint",zlim=c(0,1))
#Resume code  of habitat constraint calculation
tmp1<-grid.extent_high[order(grid.extent_high$pred_scaled,decreasing=T),]
tmp1$cum.b<-cumsum(tmp1$pred_scaled)/sum(tmp1$pred_scaled)
tmp1$cum.h<-1:nrow(tmp1)/nrow(tmp1)

all.h.sim_high<-max(tmp1$cum.h[tmp1$cum.b<=cut.off_b])
HC.sim_high<-var_ratio.sim_high*(1-all.h.sim_high)
#Add habitat extent
points(tmp1$x[tmp1$cum.b<=cut.off_b],
       tmp1$y[tmp1$cum.b<=cut.off_b],
       col=alpha('white',f=0.6),pch='+')

par(mfg=c(2,1))
barplot(c(rsq.sim_high,var_ratio.sim_high,all.h.sim_high,HC.sim_high),
          names.arg=c("Rsq","MSE_r","E","HC"),
        ylab="Habitat metrics",las=2,
        ylim=c(0,1))
box()


#2. Moderate constraint
par(mfg=c(1,2))
grid.extent_med$pred_scaled<-(grid.extent_med$pred-min(grid.extent_med$pred))/
  (max(grid.extent_med$pred)-min(grid.extent_med$pred))
range(grid.extent_med$pred_scaled)

image(lond,latd,
           t(matrix(grid.extent_med$pred_scaled,nrow=length(latd),ncol=length(lond),byrow=T)),
           col=viridis(100),ylab="",xlab="",main="Moderate constraint",
           zlim=c(0,1),axes=F)
#Resume code  of habitat constraint calculation
tmp1<-grid.extent_med[order(grid.extent_med$pred_scaled,decreasing=T),]
tmp1$cum.b<-cumsum(tmp1$pred_scaled)/sum(tmp1$pred_scaled)
tmp1$cum.h<-1:nrow(tmp1)/nrow(tmp1)

all.h.sim_med<-max(tmp1$cum.h[tmp1$cum.b<=cut.off_b])
HC.sim_med<-var_ratio.sim_med*(1-all.h.sim_med)

#Add habitat extent
points(tmp1$x[tmp1$cum.b<=cut.off_b],
       tmp1$y[tmp1$cum.b<=cut.off_b],
       col=alpha('white',f=0.4),pch='+')
par(mfg=c(2,2))
barplot(c(rsq.sim_med,var_ratio.sim_med,all.h.sim_med,HC.sim_med),
        names.arg=c("Rsq","MSE_r","E","HC"),
        ylab="",las=2,
        ylim=c(0,1))
box()

#3. Environmental constraint
par(mfg=c(1,3))
grid.extent_high_year$pred_scaled<-(grid.extent_high_year$pred-min(grid.extent_high_year$pred))/
  (max(grid.extent_high_year$pred)-min(grid.extent_high_year$pred))
range(grid.extent_high_year$pred_scaled)

image(lond,latd,
      t(matrix(grid.extent_high_year$pred_scaled,nrow=length(latd),ncol=length(lond),byrow=T)),
      col=viridis(100),ylab="",xlab="",main="Environmental constraint",
      zlim=c(0,1),axes=F)
#Resume code  of habitat constraint calculation
tmp1<-grid.extent_high_year[order(grid.extent_high_year$pred_scaled,decreasing=T),]
tmp1$cum.b<-cumsum(tmp1$pred_scaled)/sum(tmp1$pred_scaled)
tmp1$cum.h<-1:nrow(tmp1)/nrow(tmp1)

all.h.sim_high_year<-max(tmp1$cum.h[tmp1$cum.b<=cut.off_b])
HC.sim_high_year<-var_ratio.sim_high_year*(1-all.h.sim_high_year)

#Add habitat extent
points(tmp1$x[tmp1$cum.b<=cut.off_b],
       tmp1$y[tmp1$cum.b<=cut.off_b],
       col=alpha('white',f=0.6),pch='+')
par(mfg=c(2,3))
barplot(c(rsq.sim_high_year,var_ratio.sim_high_year,all.h.sim_high_year,HC.sim_high_year),
        names.arg=c("Rsq","MSE_r","E","HC"),
        ylab="",las=2,
        ylim=c(0,1))
box()

#4. Unconstrained
par(mfg=c(1,4))
grid.extent_random$pred_scaled<-(grid.extent_random$pred-min(grid.extent_random$pred))/
  (max(grid.extent_random$pred)-min(grid.extent_random$pred))
range(grid.extent_random$pred_scaled)

image(lond,latd,
           t(matrix(grid.extent_random$pred_scaled,nrow=length(latd),ncol=length(lond),byrow=T)),
           col=viridis(100),ylab="",xlab="",main="Unconstrained",
      zlim=c(0,1),axes=F)
#Resume code  of habitat constraint calculation
tmp1<-grid.extent_random[order(grid.extent_random$pred_scaled,decreasing=T),]
tmp1$cum.b<-cumsum(tmp1$pred_scaled)/sum(tmp1$pred_scaled)
tmp1$cum.h<-1:nrow(tmp1)/nrow(tmp1)

all.h.sim_random<-max(tmp1$cum.h[tmp1$cum.b<=cut.off_b])
HC.sim_random<-var_ratio.sim_random*(1-all.h.sim_random)

#Add habitat extent
points(tmp1$x[tmp1$cum.b<=cut.off_b],
       tmp1$y[tmp1$cum.b<=cut.off_b],
       col=alpha('white',f=0.6),pch='+')
par(mfg=c(2,4))
barplot(c(rsq.sim_random,var_ratio.sim_random,all.h.sim_random,HC.sim_random),
        names.arg=c("Rsq","MSE_r","E","HC"),
        ylab="",las=2,
        ylim=c(0,1))
box()
#par(oma=c( 25.9,31,0,2))

par(new=T,mfrow=c(1,4),mfg=c(1,4),omi=c(3.32,3.45,0.17,0.3))
image.plot(legend.only=T,zlim=c(0,1),
           col=viridis(100),legend.width=7,horizontal=F)

#dev.copy(jpeg,paste0(mydir,"/figures/",'Simulations.jpg'),height=7,width=10,res=300,units='in')
#dev.off()





