##################### EHF FUNCTION ###############################
#https://gist.github.com/rensa/26e39b083fe3b5b0dd59
ehf <- function(tx, tn, t95)
{
  message('Calculating ehf')
  
  # use filter() to quickly calculate the moving averages required
  t3 = rowMeans(cbind(
    filter(tx, c(rep(1/3, 3), rep(0, 2)), method = 'convolution', sides = 2, circular = FALSE),
    filter(tn, c(rep(1/3, 3), rep(0, 4)), method = 'convolution', sides = 2, circular = FALSE)),
    na.rm = TRUE)
  t30 = rowMeans(cbind(
    filter(tx, c(rep(0, 31), rep(1/30, 30)), method = 'convolution', sides = 2, circular = FALSE),
    filter(tn, c(rep(0, 29), rep(1/30, 30)), method = 'convolution', sides = 2, circular = FALSE)),
    na.rm = TRUE)
  
  # bring it all together and return ehf
  ehi.sig = t3 - t95
  ehi.accl = t3 - t30
  ehf = ehi.sig * pmax(1, ehi.accl)
  
  message('Filling in missing ehf values')
  
  # which values are missing? (except for the edges that can't be done)
  missing.vals = which(is.na(ehf))
  missing.vals =
    missing.vals[! missing.vals %in% c(1:30, (length(tx) - 2):length(tx))]
  
  # fill in missing data manually
  for (i in missing.vals)
  {
    # get rolling tx, tn windows
    t3x = tx[i:(i + 2)]
    t3n = tn[(i + 1):(i + 3)]
    t30x = tx[(i - 30):(i - 1)]
    t30n = tn[(i - 29):i]
    
    # calc ehf if there's enough data
    if (length(which(is.na(t3x))) <= 1 ||
        length(which(is.na(t3n))) <= 1 ||
        length(which(is.na(t30x))) <= 5 ||
        length(which(is.na(t30n))) <= 5)
    {
      t3 = mean(c(
        mean(t3x, na.rm = TRUE),
        mean(t3n, na.rm = TRUE)))
      t30 = mean(c(
        mean(t30x, na.rm = TRUE),
        mean(t30n, na.rm = TRUE)))
      ehf[i] = (t3 - t95) * pmax(1, t3 - t30)
    }
  }
  return(ehf)
}

# FILL DATES B/W START AND END DATE 


x<- read.csv("C:/Users/jagad/Desktop/NC_Sur/NC_Dec2019/coastal_octf.csv", header = T,
                     fileEncoding="UTF-8-BOM")

x$F1<- as.Date(x$ï..F1, format="%m/%d/%Y")
x$F2<- as.Date(x$F2, format="%m/%d/%Y")


lst <- Map(function(x, y) seq(x,y, by = "1 day"), x$F1, x$F2)
i1 <- rep(1:nrow(x), lengths(lst)) 
y<- data.frame(x[i1,-3], dates = do.call("c", lst))


write.csv(y$dates, "C:\\Users\\jagad\\Desktop\\piedmont_99_HWedf\\attach\\Tmin_99_3.csv",
          row.names=FALSE)

#####################################################
#####################################################
###########LOAD -COASTAL DATA##########
coastal <- read.csv ("C:/Users/jagad/Desktop/NC_Sur/NC_Dec2019/coastal_octf.csv", header = T,
                     fileEncoding="UTF-8-BOM")

# COMPUTE EXCESS HEAT FACTOR VARIABLE
coastal$EHF<- ehf(coastal$tmax, coastal$tmin, t95 = 82.13)

### LOAD REQUIRED LIBS ######################
library(dplyr)
library(mgcv)
library(ggplot2)
library(lubridate)
library(caret)
library(weathermetrics)
library(ThermIndex)
library(corrplot)
library(naniar)
library(MASS)
library(rcompanion)
library(dlnm)
library(splines)
library(reshape2)
library(reshape2)
library(extrafont)
library(scales)   # to access breaks/formatting functions
library(gridExtra) # for arranging plots
library(grid)   # for arranging plots
loadfonts(device = "win")

#### CHECK FOR MISSING DATA ####
vis_miss(coastal) # If more than one variable is missing use option 2/3
#gg_miss_upset(coastal) #option -2
#gg_miss_upset(coastal, nsets = n_var_miss(coastal)) #option-3

#### Missing value imputation
coastal[is.na(coastal)] <- 0

###### DATA FORMATTING ######################
coastal$date <- as.Date(coastal$date, format = "%m/%d/%Y")
coastal$month<- as.factor(month(coastal$date))
coastal$year<- as.factor(format(coastal$date, '%Y'))
coastal$dow <- wday(as.Date(coastal$date, format = "%m/%d/%Y"))
weekdays1 <- c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday')
coastal$wDay <- factor((weekdays(coastal$date) %in% weekdays1), 
                       levels=c(FALSE, TRUE), labels=c('weekend', 'weekday'))

cat <- c("tavg_95_2_coas","tavg_98_2_coas","tavg_95_3_coas","Tavg_98_3_coas",
         "tavg_90_2_coas","tavg_90_3_coas","tmax_98_3_coas","tmax_98_2_coas",
         "tmax_90_2_coas","tmax_90_3_coas","tmax_95_2_coas","tmax_95_3_coas",
         "tmin_95_3_coas","tmin_98_2_coas","Tmin_98_3_coas",
         "tmin_90_2_coas","tmin_90_3_coas","tmin_95_2_coas", "dow", "wDay",
         "coas_tavg_99_2", "coas_tavg_99_3", "coas_tmax_99_2",
         "coas_tmax_99_3", "coas_tmin_99_2", "coas_tmin_99_3")
coastal[cat] <- lapply(coastal[cat], factor)

coastal$Rate_ER_visit <- as.numeric(coastal$imp_rate)
coastal$Log_rate_ER_visit<- log(coastal$Rate_ER_visit)
coastal$Count_ER_visit <- as.numeric(coastal$imp_count)
coastal$Max_temp <- as.numeric(coastal$tmax)
coastal$Min_temp <- as.numeric(coastal$tmin)
coastal$Avg_temp <- as.numeric(coastal$tavg)
coastal$DTR <- as.numeric(coastal$tmax - coastal$tmin)
coastal$MAT<- as.numeric(coastal$�..app_temp)
coastal$Dewpoint<- as.numeric(coastal$dewpoint)

##### COMPUTING ADDITIONAL VARIABLES ######
#Steadman HI
coastal$Steadman_HI <- heat.index(t = coastal$tmax,dp = coastal$dewpoint,
                                  temperature.metric = "fahrenheit")
# US NWS HI
coastal$NWS_HI<- heat.index.algorithm(t=coastal$tmax, rh=coastal$RH)

coastal$tavg_cel<- as.numeric(tempftoc(coastal$tavg)) #convert F to C
coastal$Humidex <- as.numeric(humidex(coastal$tavg_cel, coastal$RH))
#thermal discomfort index
coastal$TDI <- as.numeric(di(coastal$tavg_cel, coastal$RH))

#Steadman definitions
coastal$HI_26 <- as.factor(ifelse(coastal$MAT >= 95.85, 1,0)) #85th pct
coastal$HI_27 <- as.factor(ifelse(coastal$MAT >= 97.16, 1,0)) #90th pct
coastal$HI_28 <- as.factor(ifelse(coastal$MAT >= 98.99, 1,0)) #95th pct

# Tan et al HW definition
coastal$tmax_cel<- as.numeric(tempftoc(coastal$tmax)) #convert F to C
coastal$HI_17<- as.factor(ifelse(coastal$tmax_cel >= 35, 1,0))

#NWS HW definitions
coastal$HI_29 <- as.factor(ifelse(coastal$NWS_HI >= 105, 1,0))
coastal$HI_30 <- as.factor(ifelse(coastal$NWS_HI > 110, 1,0))


#REGRESSION COEFFICIENT
################## Continous variables ######################
myvars <- c("Rate_ER_visit","Count_ER_visit","Max_temp","Min_temp",
            "Avg_temp","DTR","Dewpoint","RH","MAT","Steadman_HI",
            "NWS_HI","Humidex","TDI","EHF")
dat<- coastal[myvars]

#descriptive stats
desc_stat<- stat.desc(dat, basic = F)
write.csv(desc_stat,
          "C:/Users/jagad/Desktop/work/1.csv", 
          row.names = T)

#categorical variable summary
y<- (coastal[cat [-c(19,20)]])
summary(y)

#Estimating rate of ER visits by year
y1<- aggregate(coastal$imp_count, by=list(Category=coastal$year), FUN=sum)
y2<- (y1$x/1222399)*100000

#	population
#Coastal	2741101
#Mountain	
#Piedmont	5571983

################## CORRELATION MATRIX ######################
M<-cor(dat)

##################### CORR MATRIX- MATRIX FUNCTION #############
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(dat)
head(p.mat[, 1:15])

#Correlation matrix
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method="color", col=col(200),  
         type="upper", order="original", addrect = 2,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "pch", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE)

##############Density plot ########################
ggplot(coastal, aes(x=NWS_HI)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=1,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")+  # Overlay with transparent density plot
  geom_vline(aes(xintercept = mean(NWS_HI, na.rm = T)),
             colour = "red", linetype ="longdash", size = .8)

#################################
#, by=tmin_95_2_coas +s(dpt,k=3, bs='cr')+s(dif,k=3, bs='cr')
#### GAM MODEL - COASTAL REGION #############
#Add day of week; month of year and year as factors in the model
m1<- gam(imp_rate ~ s(Max_temp,k=3, bs='cr')+dow+month+year,
         family=Gamma(link = log),
         method = "GCV.Cp",
         data=coastal)

summary(m1)
m1$aic
gam.check(m1)
RMSE_m1<- sqrt(mean(residuals.gam(m1,type="response")^2))

#computing RMSE values for model
accuracy(list(m1, m2), plotit=TRUE, digits=3)

m2<- gam(imp_rate ~ s(Max_temp,k=5, bs='cr')+dow+month+year,
         family=Gamma(link = log),
         method = "GCV.Cp",
         data=coastal)

anova(m1, m2, test = "Chisq")

## GAM plot
plot.gam(m1, pages=1,
         seWithMean = TRUE,
         too.far = 0.1,
         pch =25,
         cex.lab = 1.5,
         cex.axis = 1.5,
         las=1,
         ylim = c(-1.5,2.5),
         xlab = "Maximum temperature (�F)",
         ylab = "Residual - Log Rate of ER visits per 100,000",
         shade = TRUE, shade.col = "gainsboro")

abline(v=90, col="blue")
abline(v=c(100,104), col=c("blue", "blue"), lty=c(1,2), lwd=c(1, 3))

#Prediction for temperature value
p1<-data.frame(NWS_HI=seq(70,110, by=1)) #estimates of ER rate from NWS_HI
x<- predict(m2, p1, interval = "confidence", level = 0.95) #Prediction from m1/m2 model
x1<- round((exp(x)*1222399)/100000) #calculating the rate of ER visits back to count
#saving predictions
write.csv(coastal, "C:/Users/jagad/Desktop/work/1.csv", row.names = F)

#future RCP prediction (FOR PAPER 1)
x<- read.csv("C:/Users/jagad/Desktop/NC_manus/future_clim_dat/coas_densityplt.csv",
             header = T,fileEncoding="UTF-8-BOM")
x$Daily_Maximum_Temperature<- ((x$Tmax *(9/5))+32)

#correlation b/w CCSM4 and GFDL
x<- read.csv("C:\\Users\\jagad\\Desktop\\NC_manus\\future_clim_dat\\pied.csv", header = T,
             fileEncoding="UTF-8-BOM")
cor(x$ccsm4_45_tasmax_1116, x$gfdl_45_tasmax_1116, method = "pearson")

x$RCP<-as.factor(x$RCP)

#### TEMPERATURE DENSITY PLOT BY MODEL EMISSION SCENARIO
ggplot(x, aes(x=Daily_Maximum_Temperature, colour=RCP))+ geom_density()+
  facet_grid(Model~Time.period)+
  ggtitle("Projected maximum temperature - Coastal")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c('blue','firebrick1'))+
  theme(legend.position = 'bottom')+
  theme(text=element_text(size=16,  family="Arial Black"))

###PREDICTED HRI ER VISITS FOR FUTURE TEMPERATUE DATA
x$pred_ER<-exp(predict(m1, x))
write.csv(x, "C:/Users/jagad/Desktop/NC_manus/future_clim_dat/coas_densityplt_pred.csv")

###HRI ER VISIT PREDICTION PLOTS
x<- read.csv("C:\\Users\\jagad\\Desktop\\NC_manus\\future_clim_dat\\pied_densityplt_pred.csv", header = T,
             fileEncoding="UTF-8-BOM")
x$RCP<-as.factor(x$RCP)
x$date<- as.Date(x$date, "%m/%d/%y")

##Time series pred ER visits by RCP
plt<-ggplot(x, aes(date, pred_ER, group=RCP)) +
  geom_line(aes(color=RCP)) +
  facet_grid(Model~Time.period)+
  ggtitle("HRI ER Visits during climate change - Coastal")+
  labs(x="Date by time period", y="Estimated HRI ER visits (per 100,000)")+
  theme(axis.text.x = element_text(angle = 45))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c('blue','firebrick1'))+
  theme(legend.position = 'bottom')+
  theme(text=element_text(size=16,  family="Arial Black"))

### Boxplot  - Pred ER visits by timeframe
levels(x$RCP) <- c("RCP 4.5","RCP 8.5")

plt<-ggplot(x, aes(Time.period, pred_ER, group=Time.period)) +
  geom_boxplot(varwidth = TRUE, alpha=0.2, aes(fill=Time.period))+
  facet_grid(~RCP+Model, scales="fixed")+
  ggtitle("HRI ER Visits during climate change - Piedmont")+
  labs(x="", y="Estimated HRI ER visits (per 100,000)")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "bottom")+
  theme(plot.title = element_text(hjust = 0))+
  theme(text=element_text(size=20,  family="Arial Black"))+
  theme(axis.text = element_text(size = 20, family="Arial Black"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

##############################################
##############################################
##### GENERATING DLNM results
##############################################
##############################################

library(dlnm)
library(splines)
library(lubridate)
library(mgcv)
#=================================================================================
#coastal
dat <- read.csv ("C:/Users/jagad/Desktop/NC_Sur/NC_Dec2019/dat_octf.csv", header = T,
                 fileEncoding="UTF-8-BOM")
#piedmont
dat <- read.csv ("C:/Users/jagad/Desktop/NC_Sur/NC_Dec2019/pied_octf.csv", header = T,
                 fileEncoding="UTF-8-BOM")
dat$date <- dat$Date # only for piedmont
# ==================================================================================
dat$date <- as.Date(dat$date, format = "%m/%d/%Y")
dat$dow <- wday(as.Date(dat$date, format = "%m/%d/%Y"))
dat$month<- as.factor(month(dat$date))
dat$year<- as.factor(format(dat$date, '%Y'))
weekdays1 <- c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday')
dat$wDay <- factor((weekdays(dat$date) %in% weekdays1), 
                       levels=c(FALSE, TRUE), labels=c('weekend', 'weekday'))

cat <- c("month", "dow", "wDay", "year")
dat[cat] <- lapply(dat[cat], factor)
dat$Rate_ER_visit <- as.numeric(dat$imp_rate)
dat$Count_ER_visit <- as.numeric(dat$imp_count)
dat$Max_temp <- as.numeric(dat$tmax)
# ==================================================================================
varknots <- equalknots(dat$tmax,fun="bs",df=5,degree=2)
lagknots <- logknots(5, 1)
cb3.temp <- crossbasis(dat$Max_temp, lag=5, 
                       argvar=list(fun="bs",knots=varknots), arglag=list(knots=lagknots))
model3 <- gam(Count_ER_visit ~ cb3.temp+dow+month+year,
              family=quasipoisson(), dat)
# ==================================================================================
#how to set the cen value
pred3.temp <- crosspred(cb3.temp, model3,coef=NULL, vcov = NULL, cen=90, by=1, 
                        from = 64, to = 98)

#Print risk ratio matrix
print(pred3.temp$allRRfit)

#3d plot - mostly useless
plot(pred3.temp, xlab="Temperature", zlab="RR", theta=200, phi=40, lphi=30,
     main="3D graph of temperature effect")

#relative risk | Temperature | lag
plot(pred3.temp, "contour", xlab="Temperature", key.title=title("RR"),
     plot.title=title("Contour plot",xlab="Temperature",ylab="Lag"))

#Relative risk and lag days
plot(pred3.temp, "slices", var=95, ci="n", col=1, ylim=c(0.95,1.25), lwd=1.5,
     main="Lag-response curves")

#what happens at different lags and temp values [change var=c() for temp | lag=c() for lag days] 
plot(pred3.temp, "slices", var=c(85,95), lag=c(0,5), col=4,
     ci.arg=list(density=40,col=grey(0.7)))

#what happens when exposed to 95 deg F (change var for change in temp)
plot(pred3.temp, "slices", var=95, ci="bars", type="p", col=2, pch=19,
     ci.level=0.95, main="Lag-response a 10-unit increase above threshold (95CI)")


################## ANALYSIS _ EXTENSION #####################
### 10 FOLD CROSS VALIDATION ################################
#(NOT INCLUDED IN ANY PAPER)
set.seed(2)
ind <- sample(2, nrow(coastal), replace = TRUE, prob=c(0.9, 0.1))
trainset = coastal[ind == 1,]
testset = coastal[ind == 2,]
train_control = trainControl(method = "cv", number = 10)

model<- train(inc ~ tmax,
              data=trainset,
              method ="gamSpline",
              bs='cr', #penalized cubic regression spline
              k='6',
              trControl =train_control)

summary(model)
model$performances
print(model)

#check prediction for single temperature value
predict(m1, 90)

#prediction for testset
testset$pred.inc<- predict(model, testset) #attached predicted incidence
#, interval = "prediction"

write.csv(testset[c("pred.inc", "inc")],
          "C:/Users/jagad/Desktop/work/1.csv", 
          row.names = F) 

##### PRODUCTION Climate scenarios#####
pre<- read.csv ("C:/Users/jagad/Desktop/work/ped.csv", header = T)
pre$tmax<- as.numeric((pre$t*1.8)+32) #centigrade to fahrenheit 
pre$ped.inc<-predict(model,pre)
write.csv(pre, "C:/Users/jagad/Desktop/work/1.csv", row.names = F)


############### HW DEFINITION SENSITIVITY TESTING ####################
#PAPER 2
###################

#	population
#Coastal	2741101
#Mountain	
#Piedmont	5571983

coastal$pop<- as.numeric(2741101)

hw2<- glm.nb(imp_count~tmin_90_2_coas+offset(log(coastal$pop)),
          data=coastal)
(est <- cbind(Estimate = coef(hw2), confint(hw2)))
exp(est)

################# ESTIMATING ER visits - Non-anthropogenic emission scenario  ###
#Paper 3

# ==================================================================================
########################################################################
#######################################################################
#NATURAL SIMULATION PREDICTIONS
testset<- read.csv("C:\\Users\\jagad\\Desktop\\NC_Sur\\Nautal_scenario\\coas_Nat_act.csv", header = T,
                   fileEncoding="UTF-8-BOM")
testset$date <- as.Date(testset$Date, format = "%m/%d/%Y")
testset$month<- as.factor(month(testset$date))
testset$year<- as.factor(format(testset$date, '%Y'))
testset$dow <- as.factor (wday(as.Date(testset$date, format = "%m/%d/%Y")))
testset$Max_temp<- as.numeric(testset$tmax_nat)
#testset$NWS_HI<- heat.index.algorithm(t=testset$Max_temp, rh=testset$RH)
testset$tmax_preder_vis<-exp(predict(m1, testset, interval = "confidence", level = 0.95))
#testset$NWSHI_preder_vis<-predict(m2, testset)
#testset$tmax_pred_exp<- as.numeric(exp(testset$tmax_preder_vis))
testset$tmax_pred_count<- as.numeric(((testset$tmax_pred_exp)*2741101)/100000)


#for plot testing
p + geom_line(aes(y = lwr), color = "red", linetype = "dashed")+
    geom_line(aes(y = upr), color = "red", linetype = "dashed")