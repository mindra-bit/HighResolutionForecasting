##### Forecasting Precipitation Data
##### Data: Precipitation2.csv and Precip2
##### Reference https://tem11010.github.io/timeseries-inla/



rm(list=ls())

setwd("/users/mindra/MDPI2/NewData")


#call library
library(INLA); 
library(raster); 
library(maptools); 
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
require(INLA)
require(ggplot2)
require(dplyr)
require(DAAG)
require(reshape2)

#Call Data
Precipitation<-read.csv("Precipitation2.csv", sep=";")
Precipitation2022<-Precipitation[Precipitation$Year==2022,]





Time<-as.Date(format(seq(as.Date("2022-01-01"), as.Date("2022-12-31"), by = "1 day")))
Time<-rep(Time,2)
Precipitation2022$Time<-Time



png("Figure1.png", units="in", width=10, height=6, res=100)
ggplot(Precipitation2022)+
  geom_line(aes(x=Time, y=Precipitation, group=Name, col=Name)) +theme_bw()+ 
  xlab("Time")+ labs(colour = "Station")+xlab("Year 2022")+ylab("Precipitation (mm)")+
  scale_x_date(date_breaks = "months" , date_labels = "%b")     

dev.off()




 

#################################Test


###  SETUP INLA OUTPUT ###
control <- list(
  predictor = list(compute = TRUE,link=1),
  results = list(return.marginals.random = TRUE,  return.marginals.predictor=TRUE),
  compute = list(hyperpar=TRUE, return.marginals=TRUE, return.marginals.predictor=TRUE, dic=TRUE, mlik = TRUE, cpo = TRUE, 
                 po = TRUE, waic=TRUE))



U <- 50
hyper.prec <- list(theta = list(prior = "pc.prec", param = c(U, 0.01)))


Precipitation$Kemayoran <- as.numeric(Precipitation$Name == "Kemayoran")
Precipitation$Tanjungpriok <- as.numeric(Precipitation$Name == "Tanjungpriok") 
Precipitation$DayA1 <-Precipitation$Day
Precipitation$DayA2 <-Precipitation$Day
Precipitation$DayA3 <-Precipitation$Day
Precipitation$DayA4 <-Precipitation$Day
Precipitation$MonthA1 <-Precipitation$Month
Precipitation$MonthA2 <-Precipitation$Month

F1 <- Precipitation~ f(Station) + 
                     f(DayA1, Kemayoran, model = "rw1", scale.model = TRUE, hyper = hyper.prec) + 
                     f(DayA2, Tanjungpriok, model = "rw1", scale.model = TRUE, hyper = hyper.prec)+
                     f(MonthA1, Kemayoran, model = "rw1", scale.model = TRUE, hyper = hyper.prec)+
                     f(MonthA2, Tanjungpriok, model = "rw1", scale.model = TRUE, hyper = hyper.prec)+ 
                     f(DayA3, Kemayoran, model="seasonal", season.length=30, hyper = hyper.prec) + 
                     f(DayA4, Tanjungpriok, model="seasonal", season.length=30, hyper = hyper.prec) 
  
  
  
R1 <- inla(F1, control.compute = list(dic = TRUE), family = "gaussian", data = Precipitation, verbose = TRUE)

Precipitation$Forecast<-R1$summary.fitted.values[,1] 
Precipitation$ForecastL<-R1$summary.fitted.values[,3] 
Precipitation$ForecastU<-R1$summary.fitted.values[,5] 


Time<-as.Date(format(seq(as.Date("2022-01-01"), as.Date("2023-01-31"), by = "1 day")))
Time<-rep(Time,2)
Precipitation$Time<-Time





Precip2<-read.csv("Precip2.csv", sep=",")

head(Precip2)

U <- 10
hyper.prec <- list(theta = list(prior = "pc.prec", param = c(U, 0.01)))

Precip2$Day1<-Precip2$Day
Precip2$Day2<-Precip2$Day


####### Kemayoran

F2 <- Kemayoran~   
  f(Day, model = "rw2", scale.model = TRUE, hyper = hyper.prec) +
  f(Day1, model="seasonal", season.length=30, hyper = hyper.prec)+
  f(Day2, model = "iid", hyper = hyper.prec) 


R2 <- inla(F2, control.compute = list(dic = TRUE), family = "gaussian", data = Precip2, verbose = TRUE)

Precip2$Forecast<-R2$summary.fitted.values[,1] 
Precip2$ForecastL<-R2$summary.fitted.values[,3] 
Precip2$ForecastU<-R2$summary.fitted.values[,5] 

Time<-as.Date(format(seq(as.Date("2022-01-01"), as.Date("2023-01-31"), by = "1 day")))
Precip2$Time<-Time


png("Figure2a.png", units="in", width=16, height=8, res=100)

ggplot(data = Precip2)+
  geom_line(aes(x = Time, y = Forecast, col="blue"),size=0.2)+
  geom_line(aes(x = Time, y = Kemayoran, col="red"),size=0.2)+
  geom_ribbon(aes(x = Time, ymin = ForecastL, ymax = ForecastU),  alpha = 0.2)+ 
  xlab("Time")+
  ylab("Precipitation (mm)")+
  theme_bw()+scale_x_date(date_breaks = "months" , date_labels = "%b")+theme_bw()

########Tanjungpriok

head(Precip2)


F3 <- Tanjungpriok~   
  f(Day, model = "rw2", scale.model = TRUE, hyper = hyper.prec) +
  f(Day1, model="seasonal", season.length=30, hyper = hyper.prec)+
  f(Day2, model = "iid", hyper = hyper.prec) 


R3 <- inla(F3, control.compute = list(dic = TRUE), family = "gaussian", data = Precip2, verbose = TRUE)

Precip2$Forecast1<-R3$summary.fitted.values[,1] 
Precip2$ForecastL1<-R3$summary.fitted.values[,3] 
Precip2$ForecastU1<-R3$summary.fitted.values[,5] 

Time<-as.Date(format(seq(as.Date("2022-01-01"), as.Date("2023-01-31"), by = "1 day")))
Precip2$Time<-Time


png("Figure2b.png", units="in", width=16, height=8, res=100)

ggplot(data = Precip2)+
  geom_line(aes(x = Time, y = Forecast1, col="blue"),size=0.2)+
  geom_line(aes(x = Time, y = Tanjungpriok, col="red"),size=0.2)+
  geom_ribbon(aes(x = Time, ymin = ForecastL1, ymax = ForecastU1),  alpha = 0.2)+ 
  xlab("Time")+
  ylab("Precipitation (mm)")+
  theme_bw()+scale_x_date(date_breaks = "months" , date_labels = "%b")+theme_bw()

dev.off()




#######


write.csv(Precip2,"Precip2.csv")

######



png("Figure2b.png", units="in", width=16, height=8, res=100)

ggplot(data = Precip2)+
  geom_line(aes(x = Time, y = Forecast1, col="blue"),size=0.2)+
  geom_line(aes(x = Time, y = Tanjungpriok, col="red"),size=0.2)+
  geom_ribbon(aes(x = Time, ymin = ForecastL1, ymax = ForecastU1),  alpha = 0.2)+ 
  xlab("Time")+
  ylab("Precipitation (mm)")+
  theme_bw()+scale_x_date(date_breaks = "months" , date_labels = "%b")+theme_bw()

dev.off()


#Call Data
Precip3<-read.csv("Precip3.csv", sep=";")

Time<-as.Date(format(seq(as.Date("2022-01-01"), as.Date("2023-01-31"), by = "1 day")))
Time<-rep(Time,2)
Precip3$Time<-Time



png("Figure3.png", units="in", width=10, height=6, res=100)
ggplot(Precip3)+
  geom_point(aes(x=Time, y=Precipitation,col="Observed"))+ 
  geom_line(aes(x=Time, y=Forecast, col="Forecast"))+ 
  facet_wrap(~Name, scale="free")+
  theme_bw()+ 
  xlab("Time")+ labs(colour = "Precipitation")+xlab("Year")+ylab("Precipitation (mm)")+
  scale_x_date(date_breaks = "months" , date_labels = "%b")     

dev.off()






##### SPDE
setwd("/users/mindra/MDPI2/newdata")
Precipitation<-read.csv("Precipitation2.csv", sep=";")
Precipitation<-Precipitation[Precipitation$Year==2023,]


setwd("/users/mindra/MDPI2")
#Call Jakarta MAP
Indonesia1<-readRDS('gadm36_IDN_1_sp.rds')
JKT1<-Indonesia1[Indonesia1$NAME_1 == "Jakarta Raya",]
coordinates(JKT1)
#Convert to UTM

UTM<-CRS("+proj=utm +zone=48 ellps=WGS84")

#UTM <- CRS("+proj=robin +datum=WGS84") #Transform to meters

JKT1UTM <- spTransform(JKT1, UTM)

#Generate Grid for Jakarta Map
#1. Get Outline

OutlineUTM <- JKT1UTM@polygons[[1]]@Polygons[[117]]@coords  #Border Jakarta Area without Kepulauan Seribu
plot(OutlineUTM, type="l", lwd=2, col="blue")



# Load the data for the 24 stations and 182 days
Coordinates<-data.frame(long=Precipitation$x,lat=Precipitation$y)
Coordinates <- SpatialPoints(Coordinates, proj4string = CRS("+proj=longlat")) 
Coordinates <- spTransform(Coordinates, UTM)
Coordinates<-as.data.frame(Coordinates)
colnames(Coordinates)<-c("UTMX","UTMY")

Precipitation$UTMX<-Coordinates$UTMX
Precipitation$UTMY<-Coordinates$UTMY



head(Precipitation)
mean_covariates <- apply(Precipitation[,c(11:13)],2,mean)
sd_covariates <- apply(Precipitation[,c(11:13)],2,sd)
Precipitation[,c(11:13)] <- scale(Precipitation[,c(11:13)],center=mean_covariates, scale=sd_covariates)
summary(Precipitation$Forecast)



n_stations <- 13
n_data <- 31*13 #4368 space-time data
n_days <- 31

setwd("/users/mindra/MDPI2/newdata")
Covariates<-read.csv("Covariates.csv", sep=";")
Covariates<-Covariates[Covariates$Year==2023,]

Precipitation_grid<-data.frame(UTMX=Covariates$UTMX,UTMY=Covariates$UTMY) 


head(Covariates)

mean_covariates <- apply(Covariates[,c(5,8,9)],2,mean)
sd_covariates <- apply(Covariates[,c(5,8,9)],2,sd)
Covariates[,c(5,8,9)] <- scale(Covariates[,c(5,8,9)],center=mean_covariates, scale=sd_covariates)

covariate_matrix_std <- data.frame(Covariates[,c(5,8,9)])	



# *** Code for Figure 7.8 top
plot(Precipitation_grid,col="grey",pch=18, asp=1, xlim=range(Precipitation_grid$UTMX))
lines(OutlineUTM, lwd=3, asp=1)
points(Coordinates$UTMX, Coordinates$UTMY, pch=20, cex=1,col="red")
 
Precipitation_mesh <- inla.mesh.2d(loc=cbind(Coordinates$UTMX,Coordinates$UTMY), 
                           loc.domain=OutlineUTM, offset=c(500, 7000), max.edge=c(2500, 50000))


plot(Precipitation_mesh,asp=1,main="")
lines(OutlineUTM, lwd=3)
points(Coordinates$UTMX, Coordinates$UTMY, pch=20, cex=2,col="red")
# ***

#Precipitation_spde <- inla.spde2.matern(mesh=Precipitation_mesh, alpha=2)


Precipitation_spde <- inla.spde2.pcmatern(
  mesh = Precipitation_mesh, alpha = 2, constr = TRUE,
  prior.range = c(25000, 0.1), # P(range < 5000) = 0.1
  prior.sigma = c(1, 0.1) # P(sigma > 1) = 0.1
)

#max(dist(Precipitation_grid))
#head(Precipitation_grid)
#MatCoor<-as.matrix(Coordinates)
#class(MatCoor)


A_est <- inla.spde.make.A(mesh=Precipitation_mesh,
                          loc=cbind(Coordinates$UTMX,Coordinates$UTMY),
                          group=Precipitation$Day-365,
                          n.group=n_days)
dim(A_est)


head(Precipitation)
s_index <- inla.spde.make.index(name="spatial.field",
                                n.spde=Precipitation_spde$n.spde,
                                n.group=31)
 
stack_est <- inla.stack(data=list(Prec=log(Precipitation$Forecast+0.5)),
                        A=list(A_est, 1),
                        effects=list(c(s_index,list(Intercept=1)), list(Precipitation[,c(11,12,13)])), tag="est")

head(Precipitation[,c(11,12,13)])
head(Covariates)


A_pred <- inla.spde.make.A(mesh=Precipitation_mesh,
                           loc=cbind(Precipitation_grid$UTMX,Precipitation_grid$UTMY),
                           group=Covariates$Time-365,  #selected day for prediction
                           n.group=31)

dim(A_pred)

stack_pred <- inla.stack(data=list(Prec=NA),
                         A=list(A_pred,1),
                         effects=list(c(s_index,list(Intercept=1)), list(Covariates[,c(5,8,9)])),
                         tag="pred")

stack <- inla.stack(stack_est, stack_pred)


control <- list(
  results = list(return.marginals.random = TRUE, return.marginals.predictor=TRUE),
  compute = list(hyperpar=TRUE, return.marginals.predictor=TRUE, return.marginals=TRUE, dic=TRUE, mlik = TRUE, cpo = TRUE, 
                 po = TRUE, waic=TRUE, graph=TRUE, openmp.strategy="huge")) 


formula <- Prec ~ -1 + Intercept + Altitude + UTMX + UTMY + 
  f(spatial.field, model=Precipitation_spde,group=spatial.field.group, control.group=list(model="ar1"))




# ATTENTION: the run is computationally intensive!
output <- inla(formula,
               data=inla.stack.data(stack, spde=Precipitation_spde),
               family="gaussian",
               control.compute = control$compute,
               control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),               
               control.predictor=list(A=inla.stack.A(stack), compute=TRUE), verbose=TRUE)   
  
plot(output)


#######################----------------------------------------------------------------------------------------------------
index_pred <- inla.stack.index(stack,"pred")$data
index_est <- inla.stack.index(stack,"est")$data


Covariates1<-read.csv("Covariates.csv", sep=";")
Covariates1<-Covariates1[Covariates1$Year==2023,]


Precipitation_Pred<-data.frame(x=Covariates1$UTMX,y=Covariates1$UTMY,Time=Covariates1$Time-365)

Precipitation_Est<-data.frame(x=Coordinates$UTMX,y=Coordinates$UTMY, Time=Precipitation$Day) 

head(Precipitation)
head(Covariates)

Precipitation_Pred$pred_mean <- exp(output$summary.fitted.values[index_pred, "mean"])
Precipitation_Pred$pred_ll <- exp(output$summary.fitted.values[index_pred, "0.025quant"])
Precipitation_Pred$pred_ul <- exp(output$summary.fitted.values[index_pred, "0.975quant"])
summary(Precipitation_Pred$pred_mean)


write.csv(Precipitation_Pred,"Precipitation_Pred.csv")

Precipitation_Est$pred_mean <- exp(output$summary.fitted.values[index_est, "mean"])
Precipitation_Est$pred_ll <- exp(output$summary.fitted.values[index_est, "0.025quant"])
Precipitation_Est$pred_ul <- exp(output$summary.fitted.values[index_est, "0.975quant"])
summary(Precipitation_Est$pred_mean)

Precipitation_Values<-rbind(Precipitation_Est,Precipitation_Pred)

head(Precipitation_Pred)

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

 
png("Concentration1.png", units="in", width=18, height=12, res=100)
ggplot(Precipitation_Pred,aes(x=x, y=y)) + geom_tile(aes(fill=pred_mean)) + 
  geom_contour(aes(x=x, y=y, z= pred_mean, colour=stat(level)),size=0.5,colour="white",alpha=0.5)+
  coord_fixed(ratio = 1) + scale_fill_gradientn(colours = rev(rainbow(7))) +
  facet_wrap(~ Time, labeller=as_labeller(Day_names), ncol=7)+ 
  theme_bw()+ylab("")+xlab("") + theme(legend.position="bottom", text = element_text(size=19))+ 
  guides(fill = guide_colorbar(title.position = "left", title.vjust = 1,  
                               frame.colour = "black",
                               barwidth = 25,
                               barheight = 1.5)) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  ) + labs(fill="Precipitation (mm)")

dev.off()























