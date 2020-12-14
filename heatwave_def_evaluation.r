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


###########LOAD -COASTAL DATA##########
coastal <- read.csv ("C:/Users/jagad/Desktop/NC_Sur/NC_Dec2019/coastal_octf.csv", header = T)

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
library(pastecs)

#### CHECK FOR MISSING DATA ####
vis_miss(coastal) # If more than one variable is missing use option 2/3
#gg_miss_upset(coastal) #option -2
#gg_miss_upset(coastal, nsets = n_var_miss(coastal)) #option-3

#### Missing value imputation
coastal[is.na(coastal)] <- 0

###### DATA FORMATTING ######################
coastal$date <- as.Date(coastal$date, format = "%m/%d/%Y")
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
coastal$MAT<- as.numeric(coastal$ï..app_temp)
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

################## Continous variables ######################
myvars <- c("Avg_temp","Max_temp", "Min_temp", "DTR", "Dewpoint",
            "RH", "MAT","Steadman_HI", "NWS_HI", "Humidex", "TDI","EHF",
            "Count_ER_visit","Rate_ER_visit","Log_rate_ER_visit")
dat<- coastal[myvars]

#descriptive stats
desc_stat<- stat.desc(dat, basic = F)
write.csv(desc_stat,
          "C:/Users/jagad/Desktop/work/1.csv", 
          row.names = T)

#categorical variable summary
y<- (coastal[cat [-c(19,20)]])
summary(y)


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
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
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

#Estimating rate of ER visits by year
y1<- aggregate(coastal$imp_count, by=list(Category=coastal$year), FUN=sum)
y2<- (y1$x/1222399)*100000

#	population
#Coastal	2741101
#Mountain	
#Piedmont	5571983

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
         ylim = c(-1,4),
         xlab = "Maximum temperature (°F)",
         ylab = "Residual - Log Rate of ER visits per 100,000",
         shade = TRUE, shade.col = "gainsboro")

abline(v=90, col="blue")
abline(v=c(100,104), col=c("blue", "blue"), lty=c(1,2), lwd=c(1, 3))

#Prediction for temperature value
p1<-data.frame(NWS_HI=seq(70,110, by=1)) #estimates of ER rate from NWS_HI
x<- predict(m2, p1, interval = "confidence", level = 0.95) #Prediction from m1/m2 model
x1<- round((exp(x)*1222399)/100000) #calculating the rate of ER visits back to count
#Enter cost
x2<- (x1*c) #estimating the $ value based on the count

#saving predictions
write.csv(coastal, "C:/Users/jagad/Desktop/work/1.csv", row.names = F)


################## ANALYSIS _ EXTENSION #####################
### 10 FOLD CROSS VALIDATION ################################

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
coastal$pop<- as.numeric(1222399)

hw2<- glm.nb(imp_count~tmin_90_2_coas+offset(log(coastal$pop)),
          data=coastal)
(est <- cbind(Estimate = coef(hw1), confint(hw1)))
exp(est)

################# ESTIMATING ER visits - Non-anthropogenic emission scenario  ###

testset<- read.csv("C:\\Users\\jagad\\Desktop\\NC_Sur\\Nautal_scenario\\Nat_final\\X.csv", header=T)

testset$date <- as.Date(testset$ï..Date, format = "%m/%d/%Y")
testset$dow <- as.factor (wday(as.Date(testset$date, format = "%m/%d/%Y")))
testset$Max_temp<- as.numeric(testset$tmax)
testset$NWS_HI<- heat.index.algorithm(t=testset$Max_temp, rh=testset$RH)

testset$tmax_preder_vis<-predict(m1, testset)
testset$NWSHI_preder_vis<-predict(m2, testset)
write.csv(testset, "C:/Users/jagad/Desktop/NC_Sur/Nautal_scenario/Nat_final/piedmont_predictions_nat_sce", row.names = F)