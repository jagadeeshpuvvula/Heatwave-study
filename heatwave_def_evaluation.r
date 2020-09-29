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

###########LOAD -COASTAL DATA##########
coastal <- read.csv ("C:/Users/jagad/Desktop/NC_Sur/NC_Dec2019/coastal_sepf.csv", header = T)

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
         "tmin_90_2_coas","tmin_90_3_coas","tmin_95_2_coas", "dow", "wDay")
coastal[cat] <- lapply(coastal[cat], factor)

coastal$Rate_ER_visit <- as.numeric(coastal$imp_rate)
coastal$Log_rate_ER_visit<- log(coastal$Rate_ER_visit)
coastal$Count_ER_visit <- as.numeric(coastal$imp_count)
coastal$Max_temp <- as.numeric(coastal$tmax)
coastal$Min_temp <- as.numeric(coastal$tmin)
coastal$Avg_temp <- as.numeric(coastal$tavg)
coastal$doy <- as.numeric(format(coastal$date, "%d"))
coastal$month <- as.numeric(format(coastal$date, "%m"))
coastal$year <- as.numeric(format(coastal$date, "%Y"))
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

#for HI comparision
coastal$NWS_subset <- as.factor(ifelse(coastal$NWS_HI > 100, 1,0))

################## CORRELATION MATRIX ######################
myvars <- c("Avg_temp","Max_temp", "Min_temp", "DTR", "Dewpoint",
            "RH", "MAT","Steadman_HI", "NWS_HI", "Humidex", "TDI", "Rate_ER_visit", 
            "Count_ER_visit", "Log_rate_ER_visit", "EHF")
dat<- coastal[myvars]
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

#################################
#, by=tmin_95_2_coas +s(dpt,k=3, bs='cr')+s(dif,k=3, bs='cr')
#### GAM MODEL - COASTAL REGION #############
m1<- gam(Log_rate_ER_visit ~ s(NWS_HI,k=6, bs='cr')+wDay,
         family=gaussian,
         method = "GCV.Cp",
         data=coastal)

summary(m1)
m1$aic
gam.check(m1)


m2<- gam(Lg_rate_ER_visit ~ s(NWS_HI,k=6, bs='cr'),
         family=gaussian,
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
x<- predict(m2, p1) #Prediction from m1/m2 model
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
