#### script preparation ####


# load packages
library(chillR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(rms)
library(weathermetrics)
library(ggpubr)
library(rstudioapi)


# set working directory to source file location
setwd(dirname(getActiveDocumentContext()$path)) 
getwd()


#### load and modify data ####


# temperature data
dat.temperature <- read.csv("data_T-hourly.csv", stringsAsFactors = T)

# format temeperature data
# calculate daily temperature fluctuations for each treatment
dat.temperature2 <- dat.temperature %>% 
  mutate(Date.time = as.POSIXct(Date.time, format="%m/%d/%y %H:%H")) %>% 
  mutate(Start.date = as.Date(Date.time, format="%m/%d/%y")) %>% 
  group_by(Trt, Start.date) %>% 
  mutate(Tmax = max(Temp_C), 
         Tmin = min(Temp_C)) %>% 
  mutate(Tamp = Tmax-Tmin) %>% 
  ungroup()


# budbreak data 
dat.budbreak.ori <- read.csv("data_budbreak.csv", stringsAsFactors=TRUE)


# format budbreak data
dat.budbreak.ori <- dat.budbreak.ori %>% 
  #dplyr::select(-portions, -UT, -NC) %>%
  filter(Start.date != "10/8/19") %>% 
  mutate(Start.date = as.Date(Start.date, format="%m/%d/%y")) %>% 
  group_by(Year, Trt) %>% 
  mutate(days = as.numeric(Start.date-min(Start.date))) %>%
  ungroup()


# set a shared start data since Trt==fluctuating wasn't included in 2018_19
dat.flu.fill <- dat.budbreak.ori %>% 
  filter(days == 0 & Year == "2018_19" & Trt == "field") %>% 
  mutate(Trt = "fluctuating")


# merge
dat.budbreak.ori <- rbind(dat.budbreak.ori, dat.flu.fill)


# mean cold hardiness at AA collection
dat.budbreak.temp <- dat_budbreak.ori %>% 
  filter(days %in% c(0)) %>% 
  dplyr::select(Year, Cultivar, CH_start) %>% 
  group_by(Year, Cultivar) %>%
  dplyr::summarize(CH_zero = mean(CH_start)) %>% 
  ungroup() %>% 
  mutate(CH_zero_min = max(CH_zero))


# plot initial cold hardiness (i.e. AA collection) by Year and Cultivar
ggplot(data=dat.budbreak.temp, aes(x=as.factor(Year), y=CH_zero, group=Cultivar, color=Cultivar)) +
  geom_point()+
  geom_line()+
  theme_bw()+
  theme(aspect.ratio = 1)


# merge cold hardiness with budbreak dataset
# merge by: 
intersect(names(dat.budbreak.ori), names(dat.budbreak.temp))
dat.budbreak.ori2 <- inner_join(dat.budbreak.ori, dat.budbreak.temp, by=intersect(names(dat.budbreak.ori), names(dat.budbreak.temp)))


### budbreak adjustment
dat.budbreak.ori3 <- dat.budbreak.ori2 %>%
  mutate(kdeacc_est = (-CH_start)/budbreak) %>% 
  mutate(CH_diff = CH_start-CH_zero_min) %>% 
  mutate(CH_adj = -CH_diff/kdeacc_est) %>% 
  mutate(BB_adj = budbreak-CH_adj) %>% 
  mutate(BB_adj2 = (CH_start-CH_zero_min)/(CH_start/budbreak))


# summarize into mean BB
dat.BB50 <- dat.budbreak.ori3 %>%
  group_by(Year, Start.date, Trt, Cultivar) %>% 
  #na.omit() %>% 
  dplyr::select(everything()) %>% 
  dplyr::summarize(BB50 = mean(budbreak, na.rm=T), 
                   BB50_LPI = mean(budbreak_LPI),
                   adjBB50 = mean(BB_adj, na.rm=T),
                   adjBB50_LPI = mean(BB_adj, na.rm=T),
                   CH_zero = mean(CH_zero, na.rm = T),
                   CH_start = mean(CH_start, na.rm = T)) %>%
  ungroup()


### Create color palette for treatments
pal <- c("#e64e3c", "#78B7C5", "#f5d02c") #constant, #field, #fluctuating


# plot of temperatures
ggplot(data=dat.temperature2) +
  geom_line(aes(x=Start.date, y=Temp_C, color=Trt), alpha=0.5)+
  geom_line(aes(x=Start.date, y=Tamp), color="black", lty="dashed", alpha=0.5)+
  scale_x_date(breaks="2 month", date_labels = "%b")+
  scale_y_continuous(breaks=c(-30, -20, -10, 0,  10, 20, 30))+
  labs(x=element_blank(), y="Temperature (°C)")+
  scale_color_manual(values=pal)+
  facet_grid(Trt ~ Year, scales="free_x")+
  theme_bw()+
  theme(legend.position = "none")



#### create parameter combinations dataset ####

# all combinations of levels
param.combinations <- expand.grid(Chill.max = seq(0, 24, by=2),
                                  Chill.opt = seq(-6, 16, by=1),
                                  a = c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10, 15),
                                  Amp.slope = c(0.5, 1, 2, 4, 8, 16),
                                  Amp.inflection = seq(0, 8, by=2),
                                  d = seq(1, 2, by=0.5))


# select only combinations where Chill.max > Chill.opt
param.combinations <- param.combinations %>% 
  filter(Chill.max > Chill.opt)


# create placeholder columns for RMSE to be stored in subsequent loop
param.combinations$rmse.logmod <- NA
param.combinations$rmse.nls2 <- NA


# assign data for use in loop below
dat <- dat.temperature2
datout <- NULL


#### loop: parameter testing ####

# print time at start of loop
Sys.time()

# begin loop
for (i in 1:nrow(param.combinations)) {
  
  Chill.max = param.combinations$Chill.max[i]
  Chill.opt = param.combinations$Chill.opt[i]
  a = param.combinations$a[i]
  Amp.slope = param.combinations$Amp.slope[i]
  Amp.inflection = param.combinations$Amp.inflection[i]
  d = param.combinations$d[i]
  
  datout = NULL
  
  # separating Treatments
  for (j in levels(dat$Year)) {
    
    for (k in levels(dat$Trt)){
      
    sub <- dat %>% 
      filter(Year == j,
             Trt == k)
    
    
    sub$Temp = sub$Temp_C
    sub$Temp[sub$Temp_C>Chill.max] = Chill.max # ">" for standard chill curve
    #sub$Temp[sub$Temp_C<Chill.max] = Chill.max # "<" for inverted chill curve 
    
    #standard
    sub$chill = (((sub$Temp-Chill.max)/(Chill.opt-Chill.max))^a)*exp(a*(1-(sub$Temp-Chill.max)/(Chill.opt-Chill.max)))
    
    #amplitude enhancement
    sub$AmpEnh = 1+((d-1)/(1+exp(-Amp.slope*(sub$Tamp-Amp.inflection))))
    
    sub <- sub %>% 
      mutate(chill1 = ifelse(chill>0, chill, 0)) %>% 
      mutate(chill2 = chill1 * AmpEnh) %>% 
      mutate(chill3 = cumsum(chill2)) %>% 
      mutate(final.chill = chill3)
    
    sub2 <- sub %>% 
      dplyr::select(everything()) %>% 
      group_by(Start.date) %>%
      summarize_at(c("Temp_C", "Tmax", "Tmin", "Tamp", "final.chill"),
                   max,
                   na.rm=T) %>% 
      mutate(Year = j, Trt = k) %>% 
      ungroup()
    
    
    datout <- rbind(datout, 
                    sub2)

    }}
  
  
  # merge chill calculations with budbreak dataset
  
  # # merge based on full dataset
  #dat.merge <- merge(dat_budbreak.ori3[,c("Start.date", "Trt", "budbreak")], datout[,c("Start.date", "Trt", "final.chill")]) # this wasn't used in the analysis
  
  # merge based on mean budbreak dataset
  dat.merge <- merge(dat.BB50[,c("Start.date", "Trt", "BB50", "adjBB50")], datout[,c("Start.date", "Trt", "final.chill")])
  
  # rename variable for use in models below
  dat.merge$budbreak <- dat.merge$BB50
  dat.merge$budbreak2 <- dat.merge$adjBB50
  
  # remove NA budbreak (& adjusted budbreak) rows
  dat.merge2 <- dat.merge %>% 
    filter(!is.nan(budbreak))
  
  # log model for observed budbreak
  logmod <- lm(budbreak~final.chill+I(log(1+final.chill)), data=dat.merge2)
  
  # nls model for adjusted budbreak
  nls.mod2 <- nls(budbreak2 ~
                    (10.6+5) / (rate * (0.065 + (1-0.065) / (1 + exp(-b * (final.chill - e))))),
                  data=dat.merge2,
                  start = list(rate=1.85, b=0.01, e=500),
                  #lower = list(rate=0, b=0, e=0),
                  algorithm = "port",
                  control = list(maxiter = 50000, minFactor = 1/20000, warnOnly = TRUE))
  
  # alternative models (nls for obseved & log model for adjusted)
  # these weren't used in the analysis
  
  # # log model for adjusted budbreak
  # logmod2 <- lm(budbreak2~I(log((1+final.chill))), data=dat.merge2)
  
  
  # # nls model for observed budbreak
  # nls.mod <- nls(budbreak ~
  #                  (10.6+5) / (rate * (0.065 + (1-0.065) / (1 + exp(-b * (final.chill-e))))),
  #                data=dat.merge2,
  #                start = list(rate=1.05, b=0.01, e=500),
  #                #lower = list(rate=0, b=0, e=0),
  #                algorithm = "port",
  #                control = list(maxiter = 50000, minFactor = 1/20000, warnOnly = TRUE))

  
  # generate predictions for budbreak models
  
  # list chill quantities to predict on
  chilltemp <- data.frame(final.chill = seq(0, max(datout$final.chill), by=1))
  
  #chilltemp$preds.logmod <- predict(logmod, newdata=chilltemp)
  rmse.logmod = round(sqrt(sum(resid(logmod)^2)/nrow(na.omit(dat.merge2))), 2)
  
  # adjusted budbreak nls
  #chilltemp$preds.nls2 <- predict(nls.mod2, newdata=chilltemp)
  rmse.nls2 = round(sqrt(sum(resid(nls.mod2)^2)/nrow(na.omit(dat.merge2))), 2)
  
  param.combinations$rmse.logmod[i] <- rmse.logmod
  param.combinations$rmse.nls2[i] <- rmse.nls2

  print(i)

}

# print time at end of loop
Sys.time()


#### export parameters and RMSE ####

# assign st based on time at end of loop
st=format(Sys.time(), "%y%m%d_%I%M")

# export .csv named according to system time
write.csv(param.combinations, file = paste0(st, "_chill-parameters-export.csv", sep=""), row.names = F)


# <end>




