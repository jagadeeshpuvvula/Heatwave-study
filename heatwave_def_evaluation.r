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
coastal <- read_csv ("/work/jessebell/puvvula/coastal_octf.csv")


# COMPUTE EXCESS HEAT FACTOR VARIABLE
coastal$EHF<- ehf(coastal$tmax, coastal$tmin, t95 = 82.13)

### LOAD REQUIRED LIBS ######################
library(tidyverse)
library(mgcv)
library(lubridate)
library(caret)
library(weathermetrics)
library(ThermIndex)
library(corrplot)
library(naniar)
library(MASS)
library(rcompanion)

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
coastal$temp_c <- ((coastal$tmax - 32) * 0.5555556)

################# GAM MODEL

m1<- gam(imp_rate ~ s(temp_c,k=3, bs='cr')+dow+month+year,
         family=Gamma(link = log),
         method = "GCV.Cp",
         data=coastal)

summary(m1)
m1$aic
gam.check(m1)

#computing RMSE values for model
accuracy(list(m1, m2),
         plotit=TRUE, digits=3)


m2<- gam(imp_rate ~ s(NWS_HI,k=3, bs='cr')+dow+month+year,
         family=Gamma(link = log),
         method = "GCV.Cp",
         data=coastal)

anova(m1, m2, test = "Chisq")
plot.gam(m1)
###############
#odds ratio table
or_gam(
  data = coastal, model = m1, pred = "Max_temp",
  percentage = 20, slice = TRUE)

#Generates odds ratio plots

plot_gam(
  model = m1,
  pred = "Max_temp",
  col_line = "blue",
  ci_line_col = "black",
  ci_line_type = "dashed",
  ci_fill = "grey",
  ci_alpha = 0.4,
  ci_line_size = 0.8,
  sm_fun_size = 1.1,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  limits_y = NULL, breaks_y = NULL
)
##########################
testset<- read.csv("C:\\Users\\jagad\\Desktop\\NC_Sur\\Nautal_scenario\\coas_Nat_act.csv", header = T,
                   fileEncoding="UTF-8-BOM")
testset$date <- as.Date(testset$Date, format = "%m/%d/%Y")
testset$month<- as.factor(month(testset$date))
testset$year<- as.factor(format(testset$date, '%Y'))
testset$dow <- as.factor (wday(as.Date(testset$date, format = "%m/%d/%Y")))
testset$Max_temp<- as.numeric(testset$tmax_nat)
testset$NWS_HI<- heat.index.algorithm(t=testset$Max_temp, rh=testset$RH)
testset$tmax_preder_vis<-predict(m1, testset)
testset$NWSHI_preder_vis<-predict(m2, testset)

testset$tmax_pred_exp<- as.numeric(exp(testset$tmax_preder_vis))
testset$tmax_pred_count<- as.numeric(((testset$tmax_pred_exp)*2741101)/100000)
write.csv(testset, "C:/Users/jagad/Desktop/NC_Sur/Nautal_scenario/Nat_final/coastal_predictions_nat_sce.csv", row.names = F)

plot.gam(m1, pages=1,
         seWithMean = TRUE,
         too.far = 0.1,
         pch =25,
         cex.lab = 1.5,
         cex.axis = 1.5,
         las=1,
         ylim = c(-1.5,3),
         xlim = c(65,100),
         xlab = "Maximum temperature  (°F)",
         ylab = "Residual - Rate of ER visits per 100,000",
         shade = T, shade.col = "gainsboro")

abline(v=95, col="blue")

##########################################################
##########################################################
##########################################################


#future RCP prediction

x<- read.csv("C:/Users/jagad/Desktop/NC_manus/future_clim_dat/coas_densityplt.csv",
             header = T,
             fileEncoding="UTF-8-BOM")
x$Max_temp<- ((x$Tmax *(9/5))+32)
x$Max_temp<- as.numeric(x$Max_temp)

x$pred_ER<-exp(predict(m1, x))

write.csv(x, "C:/Users/jagad/Desktop/NC_manus/future_clim_dat/coas_densityplt_pred.csv")


###########PIEDMONT DATA##########
#piedmont<- read.csv("C:/Users/jagad/Desktop/NC_Sur/NC_HW_new/Heat_NC/HeatWave_final/working_Set/piedmont.csv", header = TRUE) #old link
#piedmont[sapply(piedmont, is.integer)]<-lapply(piedmont[sapply(piedmont, is.integer)], as.factor)

piedmont <- read_csv ("/work/jessebell/puvvula/pied_octf.csv")

piedmont$EHF<- ehf(piedmont$tmax, piedmont$tmin, t95 = 80.55)

piedmont[is.na(piedmont)] <- 0

piedmont$date <- as.Date(piedmont$Date, format = "%m/%d/%Y")
piedmont$month<- as.factor(month(piedmont$date))
piedmont$year<- as.factor(format(piedmont$date, '%Y'))
piedmont$dat<- as.factor(format(piedmont$date, '%d'))
piedmont$dow <- wday(as.Date(piedmont$date, format = "%m/%d/%Y"))
weekdays2 <- c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday')
piedmont$wDay <- factor((weekdays(piedmont$date) %in% weekdays2), 
                        levels=c(FALSE, TRUE), labels=c('weekend', 'weekday'))
piedmont$temp_c <- ((piedmont$tmax - 32) * 0.5555556)


cat1 <- c("tavg_95_2_pie","tavg_98_2_pie","tavg_95_3_pie","tavg_98_3_pie","tavg_90_2_pie",
         "tavg_90_3_pie","tmax_95_3_pie","tmax_98_2_pie","tmax_98_3_pied","tmax_90_2_pie",
         "tmax_90_3_pie","tmax_95_2_pie","tmin_90_2_pie","tmin_90_3_pie","tmin_95_2_pie",
         "tmin_95_3_pie","tmin_98_2_pie","tmin_98_3_pie" ,"Tmin_99_2", "Tmin_99_3",
         "Tmax_99_2","Tmax_99_3", "Tavg_99_2","Tavg_99_3", "dow", "wDay")
piedmont[cat1] <- lapply(piedmont[cat1], factor)

piedmont$EHF<- as.numeric(piedmont$EHF)
piedmont$Rate_ER_visit <- as.numeric(piedmont$imp_rate)
piedmont$Log_rate_ER_visit<- log(piedmont$Rate_ER_visit)
piedmont$Count_ER_visit <- as.numeric(piedmont$imp_count)
piedmont$Max_temp <- as.numeric(piedmont$tmax)
piedmont$Min_temp <- as.numeric(piedmont$tmin)
piedmont$Avg_temp <- as.numeric(piedmont$TAVG)
piedmont$DTR <- as.numeric(piedmont$tmax-piedmont$tmin)
piedmont$MAT<- as.numeric(piedmont$max_appt)
piedmont$Dewpoint<- as.numeric(piedmont$dpt)
piedmont$RH<-as.numeric(piedmont$RH)

y1<- aggregate(piedmont$imp_count, by=list(Category=piedmont$year), FUN=sum)
y2<- (y1$x/5571983)*100000

#Steadman HI
piedmont$Steadman_HI <- heat.index(t = piedmont$tmax,dp = piedmont$Dewpoint,
                                  temperature.metric = "fahrenheit")
# US NWS HI
piedmont$NWS_HI<- heat.index.algorithm(t=piedmont$tmax, rh=piedmont$RH)

piedmont$tavg_cel<- as.numeric(tempftoc(piedmont$Avg_temp)) #convert F to C
piedmont$Humidex <- as.numeric(humidex(piedmont$tavg_cel, piedmont$RH))

#thermal discomfort index
piedmont$TDI <- as.numeric(di(piedmont$tavg_cel, piedmont$RH))


# Tan et al HW definition
piedmont$tmax_cel<- as.numeric(tempftoc(piedmont$tmax)) #convert F to C
piedmont$HI_17<- as.factor(ifelse(piedmont$tmax_cel >= 35, 1,0))


#Steadman definitions
piedmont$HI_26 <- as.factor(ifelse(piedmont$MAT >= 95.85, 1,0)) #85th pct
piedmont$HI_27 <- as.factor(ifelse(piedmont$MAT >= 97.16, 1,0)) #90th pct
piedmont$HI_28 <- as.factor(ifelse(piedmont$MAT >= 98.99, 1,0)) #95th pct

################## Continous variables ######################
myvars <- c("Rate_ER_visit","Count_ER_visit","Avg_temp","Max_temp", "Min_temp", "DTR", "Dewpoint",
            "RH", "MAT","Steadman_HI", "NWS_HI", "Humidex", "TDI","EHF",
            )
dat<- piedmont[myvars]

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

####
piedmont$pop <- as.numeric(5571983)

hw1<- glm.nb(imp_count~HI_28
             +offset(log(piedmont$pop)),
             data=piedmont)
(est <- cbind(Estimate = coef(hw1), confint(hw1)))
exp(est)

################# GAM MODEL

m2<- gam(imp_rate ~ s(temp_c,k=4, bs='cr')+dow+month+year,
         family=Gamma(link = log),
         method = "GCV.Cp",
         data=piedmont)

summary(m2)
m1$aic
gam.check(m1)

m4<- gam(imp_rate ~ s(temp_c,k=5, bs='cr')+dow+month+year,
         family=Gamma(link = log),
         method = "GCV.Cp",
         data=piedmont)

anova(m1, m2, test = "Chisq")
plot.gam(m1)

accuracy(list(m3, m4),
         plotit=TRUE, digits=3)

plot.gam(m1, pages=1,
         seWithMean = TRUE,
         too.far = 0.1,
         pch =25,
         cex.lab = 1.5,
         cex.axis = 1.5,
         las=1,
         ylim = c(-1.5,2.5),
         xlim = c(65,100),
         xlab = "Maximum temperature (°F)",
         ylab = "Residual - Rate of ER visits per 100,000",
         shade = TRUE, shade.col = "gainsboro")

#prediction production
#testset<- read.csv("C:\\Users\\jagad\\Desktop\\NC_Sur\\Nautal_scenario\\Nat_final\\pied_feb2020_natsim.csv", header=T)
testset<- read.csv("C:\\Users\\jagad\\Desktop\\NC_Sur\\Nautal_scenario\\Nat_final\\pied_feb2020_natsim.csv", header = T,
                   fileEncoding="UTF-8-BOM")
testset$date <- as.Date(testset$Date, format = "%m/%d/%Y")
testset$month<- as.factor(month(testset$date))
testset$year<- as.factor(format(testset$date, '%Y'))
testset$dat<- as.factor(format(testset$date, '%d'))
testset$dow <- as.factor(wday(as.Date(testset$date, format = "%m/%d/%Y")))
weekdays1 <- c('Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday')
testset$wDay <- factor((weekdays(testset$date) %in% weekdays1), 
                       levels=c(FALSE, TRUE), labels=c('weekend', 'weekday'))
testset$Max_temp<- as.numeric(testset$tmax)
testset$pied_nat<-predict(m1, testset)
testset$pied_nat_fin<- exp(testset$pied_nat)
testset$tmax_pred_count<- as.numeric(((testset$pied_nat_fin)*5571983)/100000)
write.csv(testset, "C:\\Users\\jagad\\Desktop\\work\\test_pied_predictions.csv")

####################################################
####################################################
# Manuscript figures
#Figure. 1
library(tidyverse)
dat<- read_csv("/work/jessebell/puvvula/attr_fig_one.csv")
dat$Temp_c <- ((dat$Temp - 32) * 0.5555556)

library(extrafont)
loadfonts()

ggplot(dat, aes(x=Temp_c, colour=Category, linetype=Category))+
  geom_density()+
  scale_linetype_manual(values=c(5,1))+
  scale_color_manual(values=c('blue','firebrick1'))+
  xlab("Maximum temperature (°C)")+
  ylab("Density")+
  facet_grid((Region~.), scales="fixed")+
  theme(axis.line = element_line(colour = "black"))+
  theme_bw()+
  theme(legend.position='bottom')+
  theme(text=element_text(size=12, family = "Helvetica"), 
        axis.text=element_text(size=10, family = "Helvetica"), 
        axis.title=element_text(size=12, family = "Helvetica"), 
        plot.title=element_text(size=12, family = "Helvetica"), 
        legend.text=element_text(size=12, family = "Helvetica"), 
        legend.title=element_blank()+
  strip.text.y = element_blank())

ggsave("fig_one.tiff", 
       width = 4,height = 9)

#Remove facet label        
#strip.text.y = element_blank())

#placing legend inside the plot
  #theme(strip.text.x = element_blank(),
   #     strip.background = element_rect(colour="transparent", 
    #                                    fill="transparent"),
     #   legend.position=c(0.1,0.95))


##################
##################
#Figure. 2
library(lubridate)
dat_two<- read_csv("/work/jessebell/puvvula/fig_two.csv")

dat_two$Date <- as.Date(dat_two$Date, format = "%m/%d/%Y")
dat_two$Month<- as.factor(months(dat_two$Date))
dat_two$Year<- as.factor(format(dat_two$Date, '%Y'))

#Reorder month variable
dat_two$Month <- factor(dat_two$Month, 
                        levels=c("May", "June", 
                                 "July", "August", "September"))

#percent increase in ER visit estimation
dat_two$anomaly<- ((dat_two$Observed-dat_two$Natural)/
                     dat_two$Natural)*100
#positive value - percent increase
#negative value - percent decrease

ggplot(dat_two, aes(Year, anomaly)) +
  geom_boxplot( alpha=0.2, aes(fill=Month), 
                outlier.size=1)+
  facet_grid(Region~., scales="fixed")+
  ylim(-100,1000)+
  labs(x="Year",
       y="Rate of HRI ER visits anomaly (%)")+
  geom_hline(yintercept = 0, linetype="dashed", color="grey")+
  theme_bw() + theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))+
  theme(
    legend.position='bottom',
    axis.line = element_line(colour = "black"),
    text=element_text(size=12, family = "Helvetica"),
    axis.text=element_text(size=10, family = "Helvetica"),
    axis.title=element_text(size=12, family = "Helvetica"),
    plot.title=element_text(size=12, family = "Helvetica"),
    legend.text=element_text(size=12, family = "Helvetica"),
    legend.title=element_blank())

ggsave("/work/jessebell/puvvula/fig_two.tiff", 
       width = 6,height = 6)


###########################
###########################
#Figure.3
dat_three<- read_csv("/work/jessebell/puvvula/fig_three.csv")

ggplot(dat_three, aes(Time_period, pred_ER, group=Time_period)) +
  geom_boxplot(varwidth = TRUE, alpha=0.2, aes(fill=Time_period))+
  facet_grid(region~RCP+Model, scales="fixed")+
  labs(x="", y="Estimated daily rate of HRI ER visits (per 100,000)")+
  theme_bw() + 
  theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))+
  theme(
    legend.position='bottom',
    axis.line = element_line(colour = "black"),
    text=element_text(size=12, family = "Helvetica"),
    axis.text=element_text(size=10, family = "Helvetica"),
    axis.title=element_text(size=12, family = "Helvetica"),
    plot.title=element_text(size=12, family = "Helvetica"),
    legend.text=element_text(size=12, family = "Helvetica"),
    legend.title=element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("/work/jessebell/puvvula/fig_three.tiff", 
       width = 8,height = 6)

####################################################
####################################################

#Supplement. 1
################## CORRELATION MATRIX ######################
myvars <- c("Rate_ER_visit","Count_ER_visit","Max_temp","Min_temp",
            "Avg_temp","DTR","Dewpoint","RH","MAT","Steadman_HI",
            "NWS_HI","Humidex","TDI","EHF")
dat<- coastal[myvars]
M<-cor(dat, method = "pearson")

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


#########################
#Supplement.2-A
library(visreg)

p.coas <- visreg(m1,  scale='response', "temp_c", line.par = list(col = 'red'), plot=FALSE)
p.pied <- visreg(m2,  scale='response', "temp_c", plot = FALSE)

dplyr::bind_rows(
  dplyr::mutate(p.coas$fit, plt = "Coastal"),
  dplyr::mutate(p.pied$fit, plt = "Piedmont")
) -> fits

ggplot() +
  geom_ribbon(
    data = fits, 
    aes(temp_c, ymin=visregLwr, ymax=visregUpr, group=plt), 
    fill="gray90") +
  geom_line(data = fits, aes(temp_c, visregFit, group=plt, color=plt)) +
  labs(x="Maximum temperature (°C)", 
       y="Rate of HRI ER visits (per 100,000)")+
  xlim(15,40)+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))+
  theme(
    legend.position='bottom',
    axis.line = element_line(colour = "black"),
    text=element_text(size=12, family = "Helvetica"),
    axis.text=element_text(size=10, family = "Helvetica"),
    axis.title=element_text(size=12, family = "Helvetica"),
    plot.title=element_text(size=12, family = "Helvetica"),
    legend.text=element_text(size=12, family = "Helvetica"),
    legend.title=element_blank())+
  theme(strip.text.x = element_blank(),
       strip.background = element_rect(colour="transparent", 
                                      fill="transparent"),
     legend.position=c(0.2,0.90))


ggsave("/work/jessebell/puvvula/supp_2a.tiff", 
       width = 4,height = 4)

#########################################
#Supplement. 6
##### GFDL CCSM correlation ####
x<- read.csv("C:\\Users\\jagad\\Desktop\\NC_manus\\future_clim_dat\\
             pied.csv", header = T,
             fileEncoding="UTF-8-BOM")
cor(x$ccsm4_45_tasmax_1116, x$gfdl_45_tasmax_1116, method = "pearson")

################################
#supplement. 7


x$Daily_Maximum_Temperature<- ((x$Tmax *(9/5))+32)
x$RCP<-as.factor(x$RCP)

ggplot(x, aes(x=Daily_Maximum_Temperature, colour=RCP))+ geom_density()+
  facet_grid(Model~Time.period)+
  ggtitle("Projected maximum temperature - Coastal")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c('blue','firebrick1'))+
  theme(legend.position = 'bottom')+
  theme(text=element_text(size=16,  family="Arial Black"))+
  



########################################
########################################
########################################

#Paper 2 analysis

cat <- c("tavg_95_2_coas","tavg_98_2_coas","tavg_95_3_coas","Tavg_98_3_coas",
         "tavg_90_2_coas","tavg_90_3_coas","tmax_98_3_coas","tmax_98_2_coas",
         "tmax_90_2_coas","tmax_90_3_coas","tmax_95_2_coas","tmax_95_3_coas",
         "tmin_95_3_coas","tmin_98_2_coas","Tmin_98_3_coas",
         "tmin_90_2_coas","tmin_90_3_coas","tmin_95_2_coas", "dow", "wDay",
         "coas_tavg_99_2", "coas_tavg_99_3", "coas_tmax_99_2",
         "coas_tmax_99_3", "coas_tmin_99_2", "coas_tmin_99_3")
coastal[cat] <- lapply(coastal[cat], factor)
coastal$Rate_ER_visit <- as.numeric(coastal$imp_rate)
#coastal$Log_rate_ER_visit<- log(coastal$Rate_ER_visit)
coastal$Count_ER_visit <- as.numeric(coastal$imp_count)
coastal$Max_temp <- as.numeric(coastal$tmax)
coastal$Min_temp <- as.numeric(coastal$tmin)
coastal$Avg_temp <- as.numeric(coastal$tavg)
coastal$doy <- as.numeric(format(coastal$date, "%d"))
#coastal$month <- as.numeric(format(coastal$date, "%m"))
#coastal$year <- as.numeric(format(coastal$date, "%Y"))
coastal$DTR <- as.numeric(coastal$tmax - coastal$tmin)
coastal$MAT<- as.numeric(coastal$app_temp)
coastal$Dewpoint<- as.numeric(coastal$dewpoint)
coastal$pop<- as.numeric(2741101)
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





##################### HW DEF SENSITIVITY ANALYSIS ############
ls(coastal)

hw1<- glm.nb(imp_count~tavg_90_3_coas,data=coastal)
(est <- cbind(Estimate = coef(hw1), confint(hw1)))
exp(est)


#paper 2 - main plot
################################################################3
library(tidyverse)
res<-read_csv("/work/jessebell/puvvula/heatwave_tbl.csv")

library(ggplot2)
cbbPalette <- c("#CC6666", "#9999CC", "#66CC99", "#000000")

res$Region<-as.factor(res$Region)
res$metric<-as.factor(res$metric)
levels(res$metric)
res$Metric<- factor(res$metric, levels = c("Mean daily temperature", 
                                           "Maximum daily temperature",
                                           "Minimum daily temperature",
                                           "Maximum daily apparent temperature"))

ggplot(res, aes(x = Heat_wave_definition, y = or, ymin = lcl, ymax = ucl,
                shape=Region)) + 
  geom_pointrange(aes(col = Metric), 
                  position=position_dodge(width=0.8),size = 1) + 
  ylab("Incidence rate ratio [95% CI]") +
  geom_hline(aes(yintercept = 1)) + 
  scale_colour_manual(values=cbbPalette) + 
  ggtitle("")+
  xlab("Heat wave definitions")+
  theme(legend.position = "bottom", legend.box = "vertical")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(size=15,  family="Arial Black"))+
  theme(axis.text = element_text(size = 15, family="Arial Black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(breaks = seq(0, 7, 1),limits=c(0, 7))+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  annotate("rect", xmin = 16.5, xmax = 17.5, ymin = 0, ymax = 7,
           alpha = .1,fill = "#666666")


library(tidyverse)

x<- read_csv("C:/Users/jagad/Desktop/NC_manus/HW_coas_pied_comp.csv")
x$HW_cos<- ifelse(x$coas_tmax>= 95, 1, 0)
x$HW_pied<- ifelse(x$pied_tmax>=95, 1, 0)
cor(x$HW_cos, x$HW_pied, method = "pearson")

#Matched/paired t-test
table(x$HW_cos, x$HW_pied) #2 by 2 table
prop.table(table(x$HW_cos, x$HW_pied),margin = 2)*100 #proportion table
mcnemar.test(x$HW_cos, x$HW_pied) # McNemar test (for binary data)
#t.test(x$HW_cos, x$HW_pied, paired = TRUE, alternative = "two.sided") #paired t test

#correlation between heatwave definition 17 and NWS database
x<- read_csv("C:/Users/jagad/Desktop/pied_cor.csv")
cor(x$hw_17, x$nws_hi, method = "pearson", use = "complete.obs")
prop.table(table(x$hw_17, x$nws_hi),margin = 1)*100 #margin =2 for column percentages
mcnemar.test(x$hw_17, x$nws_hi)


