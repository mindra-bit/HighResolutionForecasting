#Article R-Code 
#High-Resolution Spatiotemporal Forecasting with Missing Observations with 
#Application to Daily PM2.5 Concentrations in Jakarta Province, #Indonesia
#I Gede Nyoman Mindra Jaya 1*, Henk Folmer2

#Setting directory 

setwd("/users/mindra/@MDPI2/Revised")


rm(list=ls())

library(geodata)
library(leaflet)
library(maps)
library(raster)
library(maptools); 
library(INLA); 
library(gtools)
library(sp); 
library(spdep); 
library(parallel)
library(ggplot2)
library(extrafont)
library(ggsn)
library(splancs)
library(brinla)
library(tidyverse)
require(mgcv)
require(dplyr)
require(DAAG)
require(reshape2)
library(gstat)
library(oce)
library(terra)
library(oce)

#Step 1: Prepare for the covariate including Altitude, Population Density, Precipitation 


#Covariate 1: Population density (We have to prepare high resolution data from city to grids)

#Generate Covariate Grid for Population Density-----------------------------------------------------------------------------------------------------------------

#Call Jakarta MAP
Indonesia1<-readRDS('gadm36_IDN_1_sp.rds')
JKT1<-Indonesia1[Indonesia1$NAME_1 == "Jakarta Raya",]

#Convert to UTM
UTM<-CRS("+proj=utm +zone=48 ellps=WGS84")
JKT1UTM <- spTransform(JKT1, UTM)

#Generate Grid for Jakarta Map
#1. Get Outline

OutlineUTM <- JKT1UTM@polygons[[1]]@Polygons[[91]]@coords  #Border Jakarta Area without Kepulauan Seribu
plot(OutlineUTM, type="l", lwd=2, col="blue") 


#2. Create Grid 100 x 100
grdUTM <- makegrid(JKT1UTM, n = 10000)
colnames(grdUTM) <- c('x','y')

#3. Use inout() for clipping the rectangular grid by outline
new_grdUTM <- grdUTM[inout(grdUTM,OutlineUTM), ]

#4. Visualize your clipped grid, which can be used for kriging!
ggplot(new_grdUTM) + geom_point(aes(x=x,y=y))
coordinates(new_grdUTM)<-~x+y
gridded(new_grdUTM) = TRUE
plot(new_grdUTM)


#Call Population Density from city area

PopDensity<-read.csv("PopDensity.csv", sep=";")

#Convert to UTM

UTM<-CRS("+proj=utm +zone=48 ellps=WGS84")
 
CoordPop<-PopDensity[,c(1,2)]
PopPoints <- SpatialPoints(CoordPop, proj4string = CRS("+proj=longlat")) 

CoordPopUTM <- spTransform(PopPoints, UTM)
CoordPopUTM<-as.data.frame(CoordPopUTM)

ValuePopDensity<-PopDensity$Density
DataPopDensity<-as.data.frame(ValuePopDensity)
DataPopDensity$x<-CoordPopUTM[,1]
DataPopDensity$y<-CoordPopUTM[,2]
coordinates(DataPopDensity)<-~x+y


###Kriging is used to get high-resolution data for population density 

Variogram = variogram(DataPopDensity$ValuePopDensity ~ 1, DataPopDensity)
plot(Variogram)
Variogram.fit = fit.variogram(Variogram, model = vgm(3.5*10^8, "Sph",7000, 2*10^8))
plot(Variogram, Variogram.fit)
Variogram.Kriging = krige(DataPopDensity$ValuePopDensity ~ 1, DataPopDensity, new_grdUTM, model = Variogram.fit)  
Pred.PopDensity<-Variogram.Kriging[,-c(1,2)]
Pred.PopDensity$PopDensity<-Variogram.Kriging$var1.pred 
class(Pred.PopDensity)



png("Figure.PopDensity.Kriging.png",  units="in", width=8, height=12, res=100) 

PopDensityInterpolation<-as.data.frame(Pred.PopDensity)  
ggplot(PopDensityInterpolation,aes(x=x, y=y)) + geom_tile(aes(fill=PopDensity)) + 
  geom_contour(aes(x=x, y=y, z= PopDensity, colour=stat(level)),size=0.5,colour="white",alpha=0.5)+
  coord_fixed(ratio = 1) + scale_fill_gradientn(colours = rev(rainbow(7))) + 
  theme_bw()+ylab("")+xlab("") + theme(legend.position="bottom", text = element_text(size=19))+ 
  guides(fill = guide_colorbar(title.position = "left", title.vjust = 1,  
                               frame.colour = "black",
                               barwidth = 20,
                               barheight = 1)) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  ) 


dev.off()


###IDW: Alternatively we can use IDW
PopDensity.IDW = krige(DataPopDensity$ValuePopDensity ~ 1, DataPopDensity, new_grdUTM)    
Pred.PopDensity.IDW<-PopDensity.IDW[,-c(1,2)]
Pred.PopDensity.IDW$PopDensity<-PopDensity.IDW$var1.pred
PopDensityInterpolation.IDW<-as.data.frame(Pred.PopDensity.IDW)  

png("Figure.PopDensity.IDW.png",  units="in", width=8, height=12, res=100) 


ggplot(PopDensityInterpolation.IDW,aes(x=x, y=y)) + geom_tile(aes(fill=PopDensity)) + 
  geom_contour(aes(x=x, y=y, z= PopDensity, colour=stat(level)),size=0.5,colour="white",alpha=0.5)+
  coord_fixed(ratio = 1) + scale_fill_gradientn(colours = rev(rainbow(7))) + 
  theme_bw()+ylab("")+xlab("") + theme(legend.position="bottom", text = element_text(size=19))+ 
  guides(fill = guide_colorbar(title.position = "left", title.vjust = 1,  
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1)) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  ) +labs(fill="Population density")


dev.off()


write.csv(PopDensityInterpolation.IDW, "PopDensityInterpolation.IDW.csv")


#Covariate 2: Altitude
#Get Altitude data from worldclim

alt <- getData("worldclim", var="alt", res=.5, lon=106, lat=-6.4)

#check wether the data is correct over Jakrta province

plot(alt)

PopDensityInterpolation.IDW<-read.csv("PopDensityInterpolation.IDW.csv", sep=",")
Coords<-utm2lonlat(PopDensityInterpolation.IDW$x, PopDensityInterpolation.IDW$y, zone = 48,hemisphere = "N", km = FALSE)
long=Coords$long
lat=Coords$lat

Coords<-data.frame(long,lat) 
ID<-c(1:nrow(Coords)) 

leaflet(data=Coords) %>% 
  addProviderTiles(providers$Esri.NatGeoWorldMap) %>% 
  addCircleMarkers(~long, ~lat, label=as.character(ID)) %>%
  addRectangles(
    lng1=min(long), lat1=min(lat),
    lng2=max(long), lat2=max(lat))

#Plot altitude data
plot(alt, xlab="Longitude", ylab="Latitude", 
     ylim=c(min(lat), max(lat)),
     xlim=c(min(long), max(long)))


#Covariate 3: Precipitation

r <- getData("worldclim", var="bio", res=0.5, lon=106, lat=-6.4)
r <- r[[c(1, 12)]] 

points <- SpatialPoints(Coords, proj4string = r@crs)

# we can name these two layers
names(r) <- c("Tmean", "Prec")

# the steps to extract values for the variables you want from the Coordinates:
points <- SpatialPoints(Coords, proj4string = r@crs)
# getting temp and precip for the points

clim <- extract(r, points)
# getting the 30s altitude for the points
altS <- extract(alt, points)
# bind it all into one dataframe
climate <- cbind.data.frame(Coords, altS, clim)
  
#save climate data 
write.csv(climate,"climate.csv")

#Note we have collect all the covariate called Covariates.csv

Covariates1<-read.csv("Covariates.csv", sep=";")
Covariates1<-Covariates1[Covariates1$Year==2023,]

 
png("Figure.Altitude.png",  units="in", width=8, height=12, res=100) 

 ggplot(Covariates1,aes(x=UTMX, y=UTMY)) + geom_tile(aes(fill=Altitude)) + 
  geom_contour(aes(x=UTMX, y=UTMY, z= Altitude, colour=stat(level)),size=0.5,colour="white",alpha=0.5)+
  coord_fixed(ratio = 1) + scale_fill_gradientn(colours = rev(rainbow(7))) + 
  theme_bw()+ylab("")+xlab("") + theme(legend.position="bottom", text = element_text(size=19))+ 
  guides(fill = guide_colorbar(title.position = "left", title.vjust = 1,  
                               frame.colour = "black",
                               barwidth = 15,
                               barheight = 1)) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  ) +labs(fill="Altitude")


dev.off()

#Note for forecasting precipitation R-code can be found in Forecasting Precipitation Data.R
#Step 2: Multivariate forecast of PM2.5 
#Forecast2023


DataF<-read.csv("Data13FinNA.csv", sep=";")

DataF$Time2<-DataF$Time1

#Use INLA with halfcauchy prior 

halfcauchy = "expression:
lambda = 25;
precision = exp(log_precision);
logdens = -1.5*log_precision-log(pi*lambda)-log(1+1/(precision*lambda^2));
log_jacobian = log_precision;
return(logdens+log_jacobian);"

hcprior = list(prec = list(prior = halfcauchy))

DataF$Time1<-DataF$Time
DataF$Time2<-DataF$Time


U <- 1

hyper.prec <- list(theta = list(prior = "pc.prec", param = c(U, 0.1)))


F1 <- PM2.5 ~ f(Station) + f(Time1, model = "rw2", scale.model = TRUE, hyper =hyper.prec, cyclic=TRUE)+
  f(Time2,model="seasonal", season.length=28, hyper =hyper.prec)


R1 <- inla(F1, control.compute = list(dic = TRUE), family = "gaussian",  data = DataF)


Fit1 <- R1$summary.fitted.values

DataF$mean.R1 <- Fit1$mean

DataF$lo.R1 <- Fit1$`0.025quant`

DataF$up.R1 <- Fit1$`0.975quant`


Time<-as.Date(format(seq(as.Date("2022-01-01"), as.Date("2023-01-31"), by = "1 day")))
Time<-rep(Time,13)

DataF$TimeLabel<-Time

DataF$Time2<-DataF$Time1

png("ForecastGaussian.png", units="in", width=16, height=8, res=100)
ggplot(data = DataF, aes(x = TimeLabel, y = PM2.5))+
  geom_point(size=0.2,col="red")+
  geom_line(aes(x = TimeLabel, y = mean.R1),size=0.2)+
  geom_ribbon(aes(x = TimeLabel, ymin = lo.R1, ymax = up.R1),  alpha = 0.2)+
  facet_wrap(~Station, scale="free")+
  xlab("Time")+
  ylab("Concentration of PM2.5")+
  theme_bw()+scale_x_date(date_breaks = "months" , date_labels = "%b")+theme_bw()

dev.off()
 

#Step 3: Get High-Resolution Forecast 2023 Using SPDE
 
##### SPDE============================================================##################################################

 
Indonesia1<-readRDS('gadm36_IDN_1_sp.rds')
JKT1<-Indonesia1[Indonesia1$NAME_1 == "Jakarta Raya",]
  
UTM<-CRS("+proj=utm +zone=48 ellps=WGS84")
JKT1UTM <- spTransform(JKT1, UTM)

#Generate Grid for Jakarta Map
#1. Get Outline
OutlineUTM <- JKT1UTM@polygons[[1]]@Polygons[[91]]@coords  #Border Jakarta Area without Kepulauan Seribu
plot(OutlineUTM, type="l", lwd=2, col="blue")

# Load the data for the 13 stations with forecast values
  
PM2.5_data <- read.csv("DataNew.csv",sep=";")

PM2.5_data<-PM2.5_data[PM2.5_data$Year==2023,]

Coordinates<-data.frame(long=PM2.5_data$long,lat=PM2.5_data$lat)
Coordinates <- SpatialPoints(Coordinates, proj4string = CRS("+proj=longlat")) 
Coordinates <- spTransform(Coordinates, UTM)
Coordinates<-as.data.frame(Coordinates)
colnames(Coordinates)<-c("UTMX","UTMY")

PM2.5_data$UTMX<-Coordinates$UTMX
PM2.5_data$UTMY<-Coordinates$UTMY


n_stations <- 13
n_data <- 31*13  
n_days <- 31

head(PM2.5_data)
mean_covariates <- apply(PM2.5_data[,c(8:10,15,16)],2,mean)
sd_covariates <- apply(PM2.5_data[,c(8:10,15,16)],2,sd)
PM2.5_data[,c(8:10,15,16)] <- scale(PM2.5_data[,c(8:10,15,16)],center=mean_covariates, scale=sd_covariates)
head(PM2.5_data)

#Load for covariate data 

Covariates<-read.csv("Covariates.csv", sep=";")
Covariates<-Covariates[Covariates$Year==2023,]

PM2.5_grid<-data.frame(UTMX=Covariates$UTMX,UTMY=Covariates$UTMY) 

mean_covariates <- apply(Covariates[,5:9],2,mean)
sd_covariates <- apply(Covariates[,5:9],2,sd)
Covariates[,5:9] <- scale(Covariates[,5:9],center=mean_covariates, scale=sd_covariates)

covariate_matrix_std <- data.frame(Covariates[,5:9])	

i_day <- 1
which_date <- unique(PM2.5_data$Time)[i_day]
print(paste("**---- You will get a prediction for ", which_date, "---**"))
 
plot(PM2.5_grid,col="grey",pch=18, asp=1, xlim=range(PM2.5_grid$UTMX))
lines(OutlineUTM, lwd=3, asp=1)
points(Coordinates$UTMX, Coordinates$UTMY, pch=20, cex=2)
  
PM2.5_mesh <- inla.mesh.2d(loc=cbind(Coordinates$UTMX,Coordinates$UTMY), 
                           loc.domain=OutlineUTM, offset=c(500, 7000), max.edge=c(2500, 50000))
summary(PM2.5_mesh)
 
plot(PM2.5_mesh,asp=1,main="")
lines(OutlineUTM, lwd=3)
points(Coordinates$UTMX, Coordinates$UTMY, pch=17, cex=1, col="red")

#median(dist(PM2.5_grid))

PM2.5_spde <- inla.spde2.pcmatern(
  mesh = PM2.5_mesh, alpha = 2, constr = TRUE,
  prior.range = c(15000, 0.1), # P(range < 5000) = 0.1
  prior.sigma = c(1, 0.1) # P(sigma > 1) = 0.1
)

 
A_est <- inla.spde.make.A(mesh=PM2.5_mesh,
                          loc=cbind(Coordinates$UTMX,Coordinates$UTMY),
                          group=PM2.5_data$Time-365,
                          n.group=n_days)
dim(A_est)



s_index <- inla.spde.make.index(name="spatial.field",
                                n.spde=PM2.5_spde$n.spde,
                                n.group=31)

names(s_index)


PM2.5_data$DT<-c(1:403)
Covariates$DT<-c(1:23312)


stack_est <- inla.stack(data=list(logPM2.5P=log(PM2.5_data$PM2.5P)),
                        A=list(A_est, 1),
                        effects=list(c(s_index,list(Intercept=1)), list(PM2.5_data[,c(8:10,16,17)],Time=PM2.5_data$Time-365, DT=PM2.5_data$DT)), tag="est")

A_pred <- inla.spde.make.A(mesh=PM2.5_mesh,
                           loc=cbind(PM2.5_grid$UTMX,PM2.5_grid$UTMY),
                           group=Covariates$Time-365,  #selected day for prediction
                           n.group=31)

dim(A_pred)

stack_pred <- inla.stack(data=list(logPM2.5P=NA),
                         A=list(A_pred,1),
                         effects=list(c(s_index,list(Intercept=1)), list(covariate_matrix_std,Time=Covariates$Time-365, DT=Covariates$DT)),
                         tag="pred")

stack <- inla.stack(stack_est, stack_pred)
 

control <- list(
  results = list(return.marginals.random = TRUE, return.marginals.predictor=TRUE),
  compute = list(hyperpar=TRUE, return.marginals.predictor=TRUE, return.marginals=TRUE, dic=TRUE, mlik = TRUE, cpo = TRUE, 
                 po = TRUE, waic=TRUE, graph=TRUE, openmp.strategy="huge")) 
 

U <- 1

hyper.prec <- list(theta = list(prior = "pc.prec", param = c(U, 0.1)))



U1 <- 0.8

hyper.prec1 <- list(theta = list(prior = "pc.prec", param = c(U1, 0.4)))


Fin <- logPM2.5P ~ -1 + Intercept + Altitude + PopDensity + Precipitation + 
  f(Time, model="rw1", scale.model = TRUE, hyper =hyper.prec, cyclic=TRUE)+ 
  f(spatial.field, model=PM2.5_spde,group=spatial.field.group, control.group=list(model="ar1"), hyper =hyper.prec1)
 


# ATTENTION: the run is computationally intensive!
RFin <- inla(Fin,
             data=inla.stack.data(stack, spde=PM2.5_spde),
             family="gaussian",
             control.compute = control$compute,
             control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),               
             control.predictor=list(A=inla.stack.A(stack), compute=TRUE), verbose=TRUE)              
  
             
             
  
Fin1 <- logPM2.5P ~ -1 + Intercept + Altitude + PopDensity + Precipitation + 
  f(Time, model="rw1", scale.model = TRUE, hyper =hyper.prec)+ 
  f(spatial.field, model=PM2.5_spde,group=spatial.field.group, control.group=list(model="ar1"))



# ATTENTION: the run is computationally intensive!
RFin1 <- inla(Fin1,
             data=inla.stack.data(stack, spde=PM2.5_spde),
             family="gaussian",
             control.compute = control$compute,
             control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),               
             control.predictor=list(A=inla.stack.A(stack), compute=TRUE), verbose=TRUE)              


output1<-RFin
summary(output1)

output1<-RFin
summary(output1)

plot(output1)

output1$summary.fitted.values[,1]


PM<-20

ExProb<-unlist(lapply(output1$marginals.fitted.values, function(X){
  1-inla.pmarginal(PM, exp(X))
}))



#######################----------------------------------------------------------------------------------------------------
index_pred <- inla.stack.index(stack,"pred")$data
index_est <- inla.stack.index(stack,"est")$data


Covariates1<-read.csv("Covariates.csv", sep=";")
Covariates1<-Covariates1[Covariates1$Year==2023,]


PM2.5_Pred<-data.frame(x=Covariates1$UTMX,y=Covariates1$UTMY,Time=Covariates1$Time-365) 
PM2.5_Est<-data.frame(x=Coordinates$UTMX,y=Coordinates$UTMY, Time=PM2.5_data$Time-365) 

head(PM2.5_data)
head(Covariates)




PM2.5_Pred$pred_mean <- exp(output1$summary.fitted.values[index_pred, "mean"])
PM2.5_Pred$pred_ll <- exp(output1$summary.fitted.values[index_pred, "0.025quant"])
PM2.5_Pred$pred_ul <- exp(output1$summary.fitted.values[index_pred, "0.975quant"])
summary(PM2.5_Pred$pred_mean)




PM2.5_Est$pred_mean <- exp(output1$summary.fitted.values[index_est, "mean"])
PM2.5_Est$pred_ll <- exp(output1$summary.fitted.values[index_est, "0.025quant"])
PM2.5_Est$pred_ul <- exp(output1$summary.fitted.values[index_est, "0.975quant"])
summary(PM2.5_Est$pred_mean)

PM2.5_Values<-rbind(PM2.5_Est,PM2.5_Pred)


PM2.5_Pred$pred_prob <- ExProb[index_pred]

summary((PM2.5_data$PM2.5P))


head(PM2.5_data)
write.csv(PM2.5_Pred,"PM2.5_Pred.csv")
write.csv(PM2.5_Est,"PM2.5_Est.csv")



Time<-as.Date(format(seq(as.Date("2023-01-01"), as.Date("2023-01-31"), by = "1 day")))
Time1<-(rep(Time,each=752))
head(Time1)

PM2.5_Pred$Time<-Time1
PM2.5_Pred$Time2<-rep(c(1:31),each=752)

head(PM2.5_Pred)


png("TimeRFin1.png", units="in", width=12, height=6, res=100)
ggplot(data = PM2.5_Pred)+ 
  geom_point(aes(x = Time, y = pred_mean),size=1)+
  geom_smooth(aes(x = Time, y = pred_mean))+
  xlab("Time")+
  ylab("PM2.5")+
  theme_bw()+scale_x_date(date_breaks = 'day', 
                          date_labels = '%d\n%a')+theme_bw()

dev.off()




plot(PM2.5_Est$pred_mean,PM2.5_data$PM2.5P,ylab="Predicted PM2.5", xlab="Observed PM2.5")
abline(0,1, col="red")

RESULT2<-data.frame(Predicted=PM2.5_Est$pred_mean,Obseved=PM2.5_data$PM2.5P)

write.csv(RESULT2, "RESULT2.csv")




head(PM2.5_Pred)

Day_names <- c( "1"="Sunday (1st)",
                "2"="Monday (2nd)",
                "3"="Tuesday (3rd)",
                "4"="Wednesday (4th)",
                "5"="Thursday (5th)",
                "6"="Friday (6th)",
                "7"="Saturday (7th)",
                "8"="Sunday (8th)",
                "9"="Monday (9th)",
                "10"="Tuesday (10th)",
                "11"="Wednesday (11th)",
                "12"="Thursday (12th)",
                "13"="Friday (13th)",
                "14"="Saturday (14th)",
                "15"="Sunday (15th)",
                "16"="Monday (16th)",
                "17"="Tuesday (17th)",
                "18"="Wednesday (18th)",
                "19"="Thursday (19th)",
                "20"="Friday (20th)",
                "21"="Saturday (21st)",
                "22"="Sunday (22nd)",
                "23"="Monday (23rd)",
                "24"="Tuesday (24th)",
                "25"="Wednesday (25th)",
                "26"="Thursday (26th)",
                "27"="Friday (27th)",
                "28"="Saturday (28th)",
                "29"="Sunday (29th)",
                "30"="Monday (30th)",
                "31"="Tuesday (31th)")


PM2.5_Est1<-PM2.5_Est

coordinates(PM2.5_Est1)<-~x+y

spplot(PM2.5_Est1,zcol=c("pred_mean"))


plot(PM2.5_Est$pred_mean)

png("Concentration1.png", units="in", width=16, height=14, res=100)

ggplot(PM2.5_Pred,aes(x=x, y=y)) + geom_tile(aes(fill=pred_mean)) + 
  geom_contour(aes(x=x, y=y, z= pred_mean, colour=stat(level)),size=0.5,colour="white",alpha=0.5)+
  coord_fixed(ratio = 1) + 
  scale_fill_gradientn(colours = rev(rainbow(2))) +
  #scale_fill_gradient(low = "green", high = "red3", na.value = NA)+
  facet_wrap(~ Time2, labeller=as_labeller(Day_names), ncol=7)+ 
  theme_bw()+ylab("")+xlab("") + theme(legend.position="bottom", text = element_text(size=19))+ 
  guides(fill = guide_colorbar(title.position = "left", title.vjust = 1,  
                               frame.colour = "black",
                               barwidth = 20,
                               barheight = 1.5)) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  )  + labs(fill = "PM2.5")

dev.off()

 

png("ConcentrationProb.png", units="in", width=16, height=14, res=100)

ggplot(PM2.5_Pred,aes(x=x, y=y)) + geom_tile(aes(fill=pred_prob)) + 
  geom_contour(aes(x=x, y=y, z= pred_mean, colour=stat(level)),size=0.5,colour="white",alpha=0.5)+
  coord_fixed(ratio = 1) + 
  #scale_fill_gradientn(colours = rev(rainbow(5)),limits=c(0,1),breaks=c(0.2, 0.4, 0.6, 0.8)) +
  scale_fill_gradient(low = "green", high = "red", na.value = NA)+
    facet_wrap(~ Time2, labeller=as_labeller(Day_names), ncol=7)+ 
  theme_bw()+ylab("")+xlab("") + theme(legend.position="bottom", text = element_text(size=19))+ 
  guides(fill = guide_colorbar(title.position = "left", title.vjust = 1,  
                               frame.colour = "black",
                               barwidth = 20,
                               barheight = 1.5)) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  )  + labs(fill = "Exceedance Probability of PM2.5")

dev.off()


