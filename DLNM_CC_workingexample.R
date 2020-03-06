
##Supplementairy material
##Epidemiology and infection
##Short-term associations between legionnaires' disease incidence and meteorological variables in Belgium, 2011-2019
##T. Braeye, F. Echahidi, A. Meghraoui, V. Laisnez, N. Hens

#########
## Combining the Case-Crossover Design with Distributed Lag Non-linear Models (DLNM)
#  Written by Toon Braeye, 18/02/2020 
#########

#! In this small simulation study, the event dates are a function of negative wind speed effects (day 4,5,6), positive rel. humidity (day 5), precipitation (day 6)
#! A total of 10000 cases are simulated
#! Two datasets (rescaled meteo-average by province and population totals by province by year) are necessary for this simulation

#Libraries
library(dlnm)
library(RcppBDT)
library(timeDate)
library(lubridate)
library(dplyr)
library(survival)
library(gnm)
library(splines)
library(cowplot)
library(ggplot2)
library(stringr)
library(gridGraphics)
library(data.table)


##Directory in which the input-data can be found and in which the results can be uploaded
directory.legionella<-"C:/..."


###
# Data input 
###

#read in the meteo.agg-data
meteo.agg<-readRDS(file.path(directory.legionella, "meteoProvS.Rdata"))

#read in the population data
pop.prov.year<-read.csv(file.path(directory.legionella, "pop_prov_year.csv"))

###
# Simulate a province and a date by the meteo.agg-dataset
###

#set probabilities per meteorological factor
meteo.agg$prob.wind4 <- -meteo.agg$wind.speed*0.2  # + as.numeric(format(all.dates.all.meteo$all.dates, "%u"))*0.1),
meteo.agg$prob.wind5 <- -meteo.agg$wind.speed*0.1 #+ as.numeric(format(all.dates.all.meteo$all.dates, "%u"))*0.1),
meteo.agg$prob.wind6 <- -meteo.agg$wind.speed*0.1 # + as.numeric(format(all.dates.all.meteo$all.dates, "%u"))*0.1),
meteo.agg$prob.rel.hum <- meteo.agg$rel.humidity*0.1 # + as.numeric(format(all.dates.all.meteo$all.dates, "%u"))*0.1),
meteo.agg$prob.precipitation <- meteo.agg$precipitation*0.06 # + as.numeric(format(all.dates.all.meteo$all.dates, "%u"))*0.1),

#combine the probabilities
meteo.agg[meteo.agg$date>=as.Date('2018-12-21'), c('temperature', 'rel.humidity', 'precipitation', 'wind.speed')]<-NA

meteo.agg$prob<-exp(meteo.agg$prob.wind4 + c(NA, meteo.agg$prob.rel.hum[1:(length(meteo.agg$prob.rel.hum)-1)]) + c(NA, NA, meteo.agg$prob.precipitation[1:(length(meteo.agg$prob.precipitation)-2)]) +
                    c(NA, meteo.agg$prob.wind5[1:(length(meteo.agg$prob.wind5)-1)]) + c(NA, NA, meteo.agg$prob.wind6[1:(length(meteo.agg$prob.wind6)-2)]))

#set negative probabilities to zero
meteo.agg$prob<-ifelse(meteo.agg$prob<0, 0, meteo.agg$prob)
#remove missing data
meteo.agg<-na.omit(meteo.agg)

sim.S<-data.table(meteo.agg[sample(c(1:nrow(meteo.agg)), size=10000, prob=meteo.agg$prob, replace=T), c('date', 'province')])

legionella.agg<-as.data.frame(table(sim.S$date+4, sim.S$province))
colnames(legionella.agg)<-c('date', 'province', 'Ncases')
legionella.agg$date<-as.Date(legionella.agg$date)

legionella.DLNM<-merge(legionella.agg, meteo.agg, by=c('date', 'province'), all.y=T)
legionella.DLNM$Ncases<-ifelse(is.na(legionella.DLNM$Ncases) & legionella.DLNM$date>=as.Date('01-01-2011', format='%d-%m-%Y') & legionella.DLNM$date<as.Date('31-12-2018', format='%d-%m-%Y'), 0, legionella.DLNM$Ncases)

#plot the number of simulated cases by day
ggplot(data=legionella.DLNM, aes(x=date, y=Ncases)) + geom_line()


###
# Data prep for analysis
###

#Because a delay with ten days is studied, we have to remove the first 10-days Ncases from the input-data
legionella.DLNM[legionella.DLNM$date>=as.Date('2019-08-21'), c('temperature', 'rel.humidity', 'precipitation', 'wind.speed')]<-NA
legionella.DLNM$doy<-as.numeric(format(legionella.DLNM$date, '%j'))
legionella.DLNM$doy.factor<-as.factor(legionella.DLNM$doy)
legionella.DLNM$year<-format(legionella.DLNM$date, "%Y")
legionella.DLNM$id<-as.character(paste(format(legionella.DLNM$date, '%m'), format(legionella.DLNM$date, '%d'), legionella.DLNM$province, sep="."))
legionella.DLNM$id.factor<-as.factor(legionella.DLNM$id)

legionella.DLNM<-merge(legionella.DLNM, pop.prov.year, by.x=c('year', 'province'), by.y=c('year', 'prov'))
legionella.DLNM$inc<-legionella.DLNM$Ncases/(legionella.DLNM$pop/100000)
legionella.DLNM<-legionella.DLNM[order(legionella.DLNM$province, legionella.DLNM$date),]


legionella.DLNM$temperature2<-as.factor(c(NA, NA, legionella.DLNM$temperatureF[1:(nrow(legionella.DLNM)-2)]))
legionella.DLNM$temperature3<-as.factor(c(NA, NA, NA, legionella.DLNM$temperatureF[1:(nrow(legionella.DLNM)-3)]))
legionella.DLNM$temperature4<-as.factor(c(NA, NA, NA, NA, legionella.DLNM$temperatureF[1:(nrow(legionella.DLNM)-4)]))
legionella.DLNM$temperature5<-as.factor(c(NA, NA, NA, NA, NA, legionella.DLNM$temperatureF[1:(nrow(legionella.DLNM)-5)]))
legionella.DLNM$temperature6<-as.factor(c(NA, NA, NA, NA, NA, NA, legionella.DLNM$temperatureF[1:(nrow(legionella.DLNM)-6)]))
legionella.DLNM$temperature7<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$temperatureF[1:(nrow(legionella.DLNM)-7)]))
legionella.DLNM$temperature8<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$temperatureF[1:(nrow(legionella.DLNM)-8)]))
legionella.DLNM$temperature9<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$temperatureF[1:(nrow(legionella.DLNM)-9)]))
legionella.DLNM$temperature10<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$temperatureF[1:(nrow(legionella.DLNM)-10)]))

legionella.DLNM$wind.speed2<-as.factor(c(NA, NA, legionella.DLNM$wind.speedF[1:(nrow(legionella.DLNM)-2)]))
legionella.DLNM$wind.speed3<-as.factor(c(NA, NA, NA, legionella.DLNM$wind.speedF[1:(nrow(legionella.DLNM)-3)]))
legionella.DLNM$wind.speed4<-as.factor(c(NA, NA, NA, NA, legionella.DLNM$wind.speedF[1:(nrow(legionella.DLNM)-4)]))
legionella.DLNM$wind.speed5<-as.factor(c(NA, NA, NA, NA, NA, legionella.DLNM$wind.speedF[1:(nrow(legionella.DLNM)-5)]))
legionella.DLNM$wind.speed6<-as.factor(c(NA, NA, NA, NA, NA, NA, legionella.DLNM$wind.speedF[1:(nrow(legionella.DLNM)-6)]))
legionella.DLNM$wind.speed7<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$wind.speedF[1:(nrow(legionella.DLNM)-7)]))
legionella.DLNM$wind.speed8<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$wind.speedF[1:(nrow(legionella.DLNM)-8)]))
legionella.DLNM$wind.speed9<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$wind.speedF[1:(nrow(legionella.DLNM)-9)]))
legionella.DLNM$wind.speed10<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$wind.speedF[1:(nrow(legionella.DLNM)-10)]))

legionella.DLNM$rel.humidity2<-as.factor(c(NA, NA, legionella.DLNM$rel.humidityF[1:(nrow(legionella.DLNM)-2)]))
legionella.DLNM$rel.humidity3<-as.factor(c(NA, NA, NA, legionella.DLNM$rel.humidityF[1:(nrow(legionella.DLNM)-3)]))
legionella.DLNM$rel.humidity4<-as.factor(c(NA, NA, NA, NA, legionella.DLNM$rel.humidityF[1:(nrow(legionella.DLNM)-4)]))
legionella.DLNM$rel.humidity5<-as.factor(c(NA, NA, NA, NA, NA, legionella.DLNM$rel.humidityF[1:(nrow(legionella.DLNM)-5)]))
legionella.DLNM$rel.humidity6<-as.factor(c(NA, NA, NA, NA, NA, NA, legionella.DLNM$rel.humidityF[1:(nrow(legionella.DLNM)-6)]))
legionella.DLNM$rel.humidity7<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$rel.humidityF[1:(nrow(legionella.DLNM)-7)]))
legionella.DLNM$rel.humidity8<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$rel.humidityF[1:(nrow(legionella.DLNM)-8)]))
legionella.DLNM$rel.humidity9<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$rel.humidityF[1:(nrow(legionella.DLNM)-9)]))
legionella.DLNM$rel.humidity10<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$rel.humidityF[1:(nrow(legionella.DLNM)-10)]))

legionella.DLNM$precipitation2<-as.factor(c(NA, NA, legionella.DLNM$precipitationF[1:(nrow(legionella.DLNM)-2)]))
legionella.DLNM$precipitation3<-as.factor(c(NA, NA, NA, legionella.DLNM$precipitationF[1:(nrow(legionella.DLNM)-3)]))
legionella.DLNM$precipitation4<-as.factor(c(NA, NA, NA, NA, legionella.DLNM$precipitationF[1:(nrow(legionella.DLNM)-4)]))
legionella.DLNM$precipitation5<-as.factor(c(NA, NA, NA, NA, NA, legionella.DLNM$precipitationF[1:(nrow(legionella.DLNM)-5)]))
legionella.DLNM$precipitation6<-as.factor(c(NA, NA, NA, NA, NA, NA, legionella.DLNM$precipitationF[1:(nrow(legionella.DLNM)-6)]))
legionella.DLNM$precipitation7<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$precipitationF[1:(nrow(legionella.DLNM)-7)]))
legionella.DLNM$precipitation8<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$precipitationF[1:(nrow(legionella.DLNM)-8)]))
legionella.DLNM$precipitation9<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$precipitationF[1:(nrow(legionella.DLNM)-9)]))
legionella.DLNM$precipitation10<-as.factor(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, legionella.DLNM$precipitationF[1:(nrow(legionella.DLNM)-10)]))


#create the cross-basis for meteo-factors
lagknots0 <- logknots(c(2,10), 6)
lagknots0 <- c(5,7,9)

lagknots <- c(3:9)

varknots.temperature<-equalknots(legionella.DLNM$temperature,fun="bs",df=2,degree=1)
varknots.precipitation<-equalknots(legionella.DLNM$precipitation,fun="bs",df=2,degree=1)
varknots.rel.humidity<-equalknots(legionella.DLNM$rel.humidity,fun="bs",df=2,degree=1)
varknots.wind.speed<-equalknots(legionella.DLNM$wind.speed,fun="bs",df=2,degree=1)

varknots.doy<-equalknots(legionella.DLNM$doy,fun="bs",df=20)

#Piecewise linear basis
basis.temperature <- crossbasis(legionella.DLNM$temperature, lag=c(2,10), argvar=list(fun="lin"), arglag=list(fun="ns",knots=lagknots))
basis.rel.humidity <- crossbasis(legionella.DLNM$rel.humidity, lag=c(2,10), argvar=list(fun="lin"), arglag=list(fun="ns",knots=lagknots))
basis.wind.speed <- crossbasis(legionella.DLNM$wind.speed, lag=c(2,10), argvar=list(fun="lin"), arglag=list(fun="ns",knots=lagknots))
basis.precipitation <- crossbasis(legionella.DLNM$precipitation, lag=c(2,10), argvar=list(fun="lin"), arglag=list(fun="ns",knots=lagknots))
basis.doy <- crossbasis(legionella.DLNM$doy, minlag=0, lag=0, argvar=list(fun="bs", degree=3, knots=varknots.doy))
#cubic spline
basis.temperature2 <- crossbasis(legionella.DLNM$temperature, lag=c(2,10), argvar=list(fun="bs"), arglag=list(fun="ns",knots=lagknots))
basis.rel.humidity2 <- crossbasis(legionella.DLNM$rel.humidity, lag=c(2,10), argvar=list(fun="bs"), arglag=list(fun="ns",knots=lagknots))
basis.wind.speed2 <- crossbasis(legionella.DLNM$wind.speed, lag=c(2,10), argvar=list(fun="bs"), arglag=list(fun="ns",knots=lagknots))
basis.precipitation2 <- crossbasis(legionella.DLNM$precipitation, lag=c(2,10), argvar=list(fun="bs"), arglag=list(fun="ns",knots=lagknots))

###
# Analysis
###

#Analysis by delay (single day models)
#Fit model for every day ----
CC.2              <- gnm(Ncases ~ offset(log(pop)) +  temperature2 + rel.humidity2 + precipitation2 + wind.speed2 + year , family=poisson(), data=legionella.DLNM, eliminate=factor(id.factor))
summary(CC.2)
CC.2.sum<-data.frame(summary(CC.2)$coef)
CC.2.sum$delay=2
CC.3              <- gnm(Ncases ~ offset(log(pop)) + temperature3 + rel.humidity3 + precipitation3 + wind.speed3 + year , family=poisson(), data=legionella.DLNM, eliminate=factor(id.factor))
summary(CC.3)
CC.3.sum<-data.frame(summary(CC.3)$coef)
CC.3.sum$delay=3
CC.4              <- gnm(Ncases ~  offset(log(pop)) + temperature4 + rel.humidity4 + precipitation4 + wind.speed4 + year , family=poisson(), data=legionella.DLNM, eliminate=factor(id.factor))
summary(CC.4)
CC.4.sum<-data.frame(summary(CC.4)$coef)
CC.4.sum$delay=4
CC.5              <- gnm(Ncases ~  offset(log(pop)) + temperature5 + rel.humidity5 + precipitation5 + wind.speed5 + year , family=poisson(), data=legionella.DLNM, eliminate=factor(id.factor))
summary(CC.5)
CC.5.sum<-data.frame(summary(CC.5)$coef)
CC.5.sum$delay=5
CC.6              <- gnm(Ncases ~  offset(log(pop)) + temperature6 + rel.humidity6 + precipitation6 + wind.speed6 + year , family=poisson(), data=legionella.DLNM, eliminate=factor(id.factor))
summary(CC.6)
CC.6.sum<-data.frame(summary(CC.6)$coef)
CC.6.sum$delay=6
CC.7              <- gnm(Ncases ~  offset(log(pop)) + temperature7 + rel.humidity7 + precipitation7 + wind.speed7 + year , family=poisson(), data=legionella.DLNM, eliminate=factor(id.factor))
summary(CC.7)
CC.7.sum<-data.frame(summary(CC.7)$coef)
CC.7.sum$delay=7
CC.8              <- gnm(Ncases ~  offset(log(pop)) + temperature8 + rel.humidity8 + precipitation8 + wind.speed8  + year , family=poisson(), data=legionella.DLNM, eliminate=factor(id.factor))
summary(CC.8)
CC.8.sum<-data.frame(summary(CC.8)$coef)
CC.8.sum$delay=8
CC.9              <- gnm(Ncases ~  offset(log(pop)) + temperature9 + rel.humidity9 + precipitation9 + wind.speed9  + year , family=poisson(), data=legionella.DLNM, eliminate=factor(id.factor))
summary(CC.9)
CC.9.sum<-data.frame(summary(CC.9)$coef)
CC.9.sum$delay=9
CC.10              <- gnm(Ncases ~  offset(log(pop)) + temperature10 + rel.humidity10 + precipitation10 + wind.speed10  + year , family=poisson(), data=legionella.DLNM, eliminate=factor(id.factor))
summary(CC.10)
CC.10.sum<-data.frame(summary(CC.10)$coef)
CC.10.sum$delay=10

CC<-rbind(CC.2.sum,CC.3.sum,CC.4.sum,CC.5.sum,CC.6.sum,CC.7.sum,CC.8.sum,CC.9.sum,CC.10.sum)
CC<-CC[!grepl('year', rownames(CC)),]
CC<-CC[!grepl('province', rownames(CC)),]

# Create a boxplot with the different significant effects. ====

RD<-CC[CC$Pr...z..<0.05,]

RD$LB<-RD$Estimate-1.96*RD$Std..Error
RD$UB<-RD$Estimate+1.96*RD$Std..Error

RD$var<-rownames(RD)
RD$mis.rel.hum<-grepl('rel.hum', RD$var)
RD$mis.temp<-grepl('temp', RD$var)
RD$mis.wind<-grepl('wind', RD$var)
RD$mis.prec<-grepl('prec', RD$var)
RD$mis.pres<-grepl('pres', RD$var)

RD$Label<-ifelse(RD$mis.pres==T, 'Atm Pressure', ifelse(RD$mis.wind==T, 'Wind Speed', ifelse(RD$mis.rel.hum==T, 'Rel Humidity', ifelse(RD$mis.temp==T, 'Temperature', 'Precipitation'))))


#Piecewise linear DLNM model after CC data collection 
CC.DLNM          <- gnm(Ncases ~  offset(log(pop)) + basis.rel.humidity + basis.temperature2 + year + 
                          basis.precipitation + basis.wind.speed, family=poisson(), data=legionella.DLNM, eliminate=factor(id.factor))
summary(CC.DLNM)

###
# Results Presentation
###

### single day models 

plot<- ggplot(data=RD, mapping=aes(x=delay, y=Estimate, ymin = LB, ymax = UB)) + 
  geom_hline(yintercept=0, colour='gray') + 
  geom_point(aes(group=Label, colour=Label, shape=Label), size=4) + 
  # scale_y_continuous(limits=c(0,ceiling(max(orthssi$Upper997CI)))) + 
  geom_linerange(mapping=aes(x=delay))   +
  ggtitle('Statistically significant coefficients') +
  xlab('Lag (days)') +
  ylab('Coefficient') + 
  theme_bw() +  theme(panel.grid.major.x = element_blank())

ggsave(file.path(directory.legionella,"fig3_sep_signcoef15022020.png"), plot)

### DLNM

#----
#Results Temperature
cen.temp<-round(median(legionella.DLNM$temperature, na.rm=T))
cen.temp
pred.temperature.cc <- crosspred(basis.temperature2, CC.DLNM, by=1, cen=cen.temp, lag=c(2,10))

overall.plot.temperature <- ~{plot(pred.temperature.cc,"overall",xlab="Temperature (°C)", ci="l", col=3, ylim=c(0.5,3.5), lwd=2,
                                   ci.arg=list(col=1,lty=3), main="Temperature", ylab='RR')}

lines.plot.temperature <- ~{plot(pred.temperature.cc, "slices", var=0, ci="l", col=1, ylim=c(0,3), lwd=1.5,
                                 main="Temperature", ylab='RR', xlab='Lag (days)')
  for(i in 1:2) lines(pred.temperature.cc, "slices", var=c(-1,1)[i], ci="l", col=i+1, lwd=1.5)
  legend("topright",paste("Temperature =",c(-1,0,1)), col=c(3,1,2), lwd=1.5)}


#Results precipitation
cen.prec<-round(median(legionella.DLNM$precipitation, na.rm=T))
cen.prec
pred.precipitation.CC <- crosspred(basis.precipitation, CC.DLNM, by=1, cen=cen.prec, lag=c(2,10))


overall.plot.precipitation <- ~{
  plot(pred.precipitation.CC,"overall",xlab="Precipitation (mm)", ci="l", col=3, ylim=c(0,4.5), lwd=2,
       ci.arg=list(col=1,lty=3), main="Precipitation", ylab='RR')}

lines.plot.precipitation<- ~{plot(pred.precipitation.CC, "slices", var=0, ci="l", col=1, ylim=c(0,3), lwd=1.5,
                                  main="Precipitation", ylab='RR', xlab='Lag (days)')
  for(i in 1:2) lines(pred.precipitation.CC, "slices", var=c(1,2)[i], ci="l", col=c('red2', 'red4')[i], lwd=1.5)
  legend("topright",paste("precipitation =",c(0,1,2)), col=c(1, 'red2', 'red4'), lwd=1.5)}


#Results Relative Humidity
cen.rh<-round(median(legionella.DLNM$rel.humidity, na.rm=T))
cen.rh
pred.rel.humidity.CC <- crosspred(basis.rel.humidity, CC.DLNM, by=1, cen=cen.rh, lag=c(2,10))

overall.plot.rel.humidity<- ~{plot(pred.rel.humidity.CC,"overall",xlab="Rel. Humidity (%)", ci="l", col=3, ylim=c(0,3.5), lwd=2,
                                   ci.arg=list(col=1,lty=3), main="Relative Humidity", ylab='RR')}

lines.plot.rel.humidity<- ~{plot(pred.rel.humidity.CC, "slices", var=0, ci="l", col=1, ylim=c(0,2), lwd=1.5,
                                 main="Relative Humidity", ylab='RR', xlab='Lag (days)')
  for(i in 1:2) lines(pred.rel.humidity.CC, "slices", var=c(-1,1)[i], ci="l", col=c('green', 'red')[i], lwd=1.5)
  legend("topright",paste("Rel.Humidity =",c(-1,0,1)), col=c('green', 1, 'red'), lwd=1.5)}


#Results Wind Speed
cen.ws<-round(median(legionella.DLNM$wind.speed, na.rm=T))
cen.ws
pred.wind.speed <- crosspred(basis.wind.speed, CC.DLNM, by=1, cen=cen.ws, lag=c(2,10))

overall.plot.wind.speed<- ~{plot(pred.wind.speed,"overall",xlab="Wind Speed (m/s)", ci="l", col=3, ylim=c(0,3), lwd=2,
                                 ci.arg=list(col=1,lty=3), main="Wind Speed", ylab='RR')}

lines.plot.wind.speed<- ~{plot(pred.wind.speed, "slices", var=0, ci="l", col=1, ylim=c(0,2), lwd=1.5,
                               main="Wind Speed", ylab='RR', xlab='Lag (days)')
  for(i in 1:2) lines(pred.wind.speed, "slices", var=c(-1,1)[i], ci="l", col=c('green', 'red')[i], lwd=1.5)
  legend("topright",paste("Wind Speed (m/s) =",c(-1,0,1)), col=c('green', 1, 'red'), lwd=1.5)}


##RR accumulated over lag
title.overall<-ggdraw() + draw_label("RR accumulated over lag (10-2 days)", fontface='bold')
overall.plots<-plot_grid(overall.plot.temperature, overall.plot.precipitation, overall.plot.rel.humidity, overall.plot.wind.speed, labels=c('A', 'B', 'C', 'D'))
overall.plots2<-plot_grid(title.overall, overall.plots, ncol=1, rel_heights = c(0.1,1))
##Lag-response curves
title.lines<-ggdraw() + draw_label("Lag-response curves", fontface='bold')
lines.plots<-plot_grid(lines.plot.temperature, lines.plot.precipitation, lines.plot.rel.humidity, lines.plot.wind.speed, labels=c('A', 'B', 'C', 'D'))
lines.plots2<-plot_grid(title.lines, lines.plots, ncol=1, rel_heights=c(0.1,1))


save_plot(filename = file.path(directory.legionella, "fig5a_DLNM_lines.png"), lines.plots2,
          ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 3, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)

save_plot(filename = file.path(directory.legionella, "fig4a_DLNM_overall.png"), overall.plots2,
          ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 3, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)


