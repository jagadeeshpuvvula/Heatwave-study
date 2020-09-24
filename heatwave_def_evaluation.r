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