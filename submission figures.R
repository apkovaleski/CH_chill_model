#### script preparation ####

# load packages
library(rstudioapi)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(chillR)
library(weathermetrics)
library(rms)
library(ggpubr)
library(patchwork)
library(scales)


# set working directory to source file location
setwd(dirname(getActiveDocumentContext()$path)) 
getwd()

#### load and modify data ####

# temperature data
dat.temperature <- read.csv("data_T-hourly.csv", stringsAsFactors = T)

# format temperature data
dat.temperature2 <- dat.temperature %>% 
  mutate(Year.Trt = as.factor(paste(Year, Trt, sep="."))) %>% 
  mutate(Date.time = as.POSIXct(Date.time, format="%m/%d/%y %H:%H")) %>% 
  mutate(Start.date = as.Date(Date.time, format="%m/%d/%y")) %>% 
  group_by(Trt, Start.date) %>% 
  # calculate daily temperature fluctuations for each treatment
  mutate(Tmax = max(Temp_C), 
         Tmin = min(Temp_C),
         Tamp = Tmax-Tmin) %>% 
  ungroup()


# collection data
dat.collections <- read.csv("data_collection dates 3 yrs.csv", stringsAsFactors=TRUE)

dat.collections <- dat.collections %>% 
  mutate(Start.date = as.Date(Start.date, format="%m/%d/%y"))


# budbreak data
dat.budbreak.ori <- read.csv("data_budbreak.csv", stringsAsFactors=TRUE)

# format budbreak data
dat.budbreak.ori <- dat.budbreak.ori %>% 
  #dplyr::select(-portions, -UT, -NC) %>%
  filter(Start.date != "10/8/19") %>% 
  mutate(Start.date = as.Date(Start.date, format="%m/%d/%y")) %>% 
  group_by(Year, Trt) %>% 
  mutate(days = as.numeric(Start.date-min(Start.date))) %>%
  #filter(budbreak < 70) %>% 
  ungroup()

# set a shared start data since Trt==fluctuating wasn't included in 2018_19
dat.flu.fill <- dat.budbreak.ori %>% 
  filter(days == 0 & Year == "2018_19" & Trt == "field") %>% 
  mutate(Trt = "fluctuating")

# merge
dat.budbreak.ori <- rbind(dat.budbreak.ori, dat.flu.fill)

# mean cold hardiness at AA collection
dat.budbreak.temp <- dat.budbreak.ori %>% 
  filter(days %in% c(0)) %>% 
  dplyr::select(Year, Cultivar, CH_start) %>% 
  group_by(Year, Cultivar) %>%
  dplyr::summarize(CH_zero = mean(CH_start)) %>% 
  ungroup() %>% 
  mutate(CH_zero_min = max(CH_zero))
dat.budbreak.temp

# merge cold hardiness with budbreak dataset
# merge by:
intersect(names(dat.budbreak.ori), names(dat.budbreak.temp))
dat.budbreak.ori2 <- inner_join(dat.budbreak.ori, dat.budbreak.temp, by=intersect(names(dat.budbreak.ori), names(dat.budbreak.temp)))


# budbreak adjustment based on cold hardiness
dat.budbreak.ori3 <- dat.budbreak.ori2 %>%
  #estimated deacclimtion rate
  mutate(kdeacc_est = (-CH_start)/budbreak) %>% 
  #CH gained since start of chill treatments
  mutate(CH_diff = CH_start-CH_zero_min) %>% 
  #adjustment value
  mutate(CH_adj = -CH_diff/kdeacc_est) %>% 
  #subtracting adjustment
  mutate(BB_adj = budbreak-CH_adj) %>% 
  mutate(BB120d = ifelse(budbreak<=120, budbreak, NA))


# summarize into mean BB
dat.BB50 <- dat.budbreak.ori3 %>%
  group_by(Year, Start.date, Trt) %>% 
  #na.omit() %>% 
  dplyr::select(everything()) %>% 
  dplyr::summarize(BB50 = mean(budbreak, na.rm=T), 
                   BB20 = quantile(budbreak, 0.20, na.rm=T),
                   BB80 = quantile(budbreak, 0.80, na.rm=T),
                   BB50_LPI = mean(budbreak_LPI),
                   adjBB50 = mean(BB_adj, na.rm=T),
                   adjBB20 = quantile(BB_adj, 0.20, na.rm=T),
                   adjBB80 = quantile(BB_adj, 0.80, na.rm=T),
                   adjBB50_LPI = mean(BB_adj, na.rm=T),
                   SampleSize = n(),
                   countBB = sum(!is.na(BB120d)),
                   percBB = (countBB/SampleSize)*100,
                   CH_zero = mean(CH_zero, na.rm = T),
                   CH_start = mean(CH_start, na.rm = T),
                   CH_adj = mean(CH_adj, na.rm = T),
                   UT = mean(UT, na.rm = T),
                   NC = mean(NC, na.rm = T),
                   portions = mean(portions, na.rm=T),
                   days = mean(days)) %>%
  ungroup()


dat.BB50 <- dat.BB50 %>% 
  mutate(BB80B = ifelse(percBB<=30, NA, BB80),
         BB20B = ifelse(percBB<=30, NA, BB20),
         BBcensored = ifelse(days<=17 & Trt == "constant", 60.5,
                             ifelse(days<=12 & Trt == "field", 59.5,
                                    ifelse(days<=13 & Trt == "fluctuating", 60, NA))),
         BBcensored2 = ifelse(days==17 & Trt == "constant", BB20B,
                              ifelse(days==12 & Trt == "field", BB20B,
                                     ifelse(days==13 & Trt == "fluctuating", BB20B, NA))))



#### calculate and cumulatively sum chill units ####

## North Carolina model ##
# data frame for chill efficiency curve
data_NorthCarolina <- data.frame(NC_temperature = c(-1.1, 1.6, 7.2, 13, 16.5, 19, 20.7, 22.1, 23.3),
                                 NC_chill.units = c(0, 0.5, 1, 0.5, 0, -0.5, -1, -1.5, -2))
# create a lm with splines for temperature
spline.mod.NC <- lm(NC_chill.units ~ rcs(NC_temperature, c(-1.1, 1.6, 7.2, 13, 16.5, 19, 20.7, 22.1, 23.3)), data=data_NorthCarolina)
# dataframe for predictions and corresponding temperature values
data_NC_spline <- data.frame(NC_temperature=seq(-1,23.3, by=0.1))
# predict chill efficiency curve
data_NC_spline <- data_NC_spline %>% 
  mutate(spline.fit = predict(spline.mod.NC, newdata = data_NC_spline))
# ggplot()+
#   geom_line(data=NC_spline, aes(x=NC_temperature, y=spline.fit))+
#   geom_point(data=data_NorthCarolina, aes(x=NC_temperature, y=NC_chill.units))+
#   theme_bw()+
#   theme(aspect.ratio = 1)


## New Chill model ##
# parameters with min rmse (rmse == 3.98)
Chill.max = 16
Chill.opt = 13
a = 0.5
Amp.slope = 4
Amp.inflection = 4
d = 1.5

dat.ChillLoop <- dat.temperature2
datout.ChillLoop <- NULL

for (j in levels(dat.ChillLoop$Year.Trt)) {
  
  # subset by year and treatment
  sub <- dat.ChillLoop %>% 
    filter(Year.Trt == j)
  
  # assign temperatures based on subset
  sub$Temp = sub$Temp_C
  # avoid chill negation based on curve shape
  #  temperatures > chill.max are set equal to chill efficiency at chill.max (which is 0)
  sub$Temp[sub$Temp_C>Chill.max] = Chill.max
  
  # new chill calculation
  sub$chill = (((sub$Temp-Chill.max)/(Chill.opt-Chill.max))^a)*exp(a*(1-(sub$Temp-Chill.max)/(Chill.opt-Chill.max)))
  # new chill amplitude enhancement
  sub$AmpEnh = 1+((d-1)/(1+exp(-Amp.slope*(sub$Tamp-Amp.inflection))))
  
  # cumulative chill calculations for all chill models
  sub <- sub %>% 
    # New chill model
    # redundant step to avoid chill negation 
    mutate(chill1 = ifelse(chill>0, chill, 0)) %>% 
    # chill efficiency * amplitude enhancement
    mutate(chill2 = chill1 * AmpEnh) %>% 
    # cumulative summation
    mutate(chill3 = cumsum(chill2)) %>% 
    # store cumulative summation
    mutate(final.chill = chill3) %>% 
    # Dynamic chill model
    # uses function in chillR package
    mutate(portions = Dynamic_Model(Temp_C)) %>% 
    # Utah chill model
    # uses function in chillR package
    mutate(UT1 = Utah_Model(Temp_C)) %>% 
    # second cumulative summation step
    # this sets start date for cumulative chill based on Richardson et al. 1974
    mutate(UT = ifelse(row_number() < which.min(UT1), 0, abs(min(UT1))+UT1)) %>%
    # North Carolina chill model
    # predict according to spline model made above
    mutate(NC1 = ifelse(Temp_C >= -1.1 & Temp_C <= 23.3,
                        (predict(spline.mod.NC, data.frame(NC_temperature = Temp_C))), 
                        0)) %>% 
    # cumulative summation
    mutate(NC2 = cumsum(NC1)) %>% 
    # second cumulative summation step
    # this sets start date for cumulative chill based on Richardson et al. 1974
    mutate(NC = ifelse(row_number() < which.min(NC2), 0, abs(min(NC2))+NC2))
  
  # convert hourly calculations to daily values
  sub2 <- sub %>% 
    dplyr::select(everything()) %>% 
    group_by(Start.date) %>%
    summarise_at(c("Temp_C", "Tmax", "Tmin", "Tamp", "final.chill", "portions", "UT1", "UT", "NC1", "NC2", "NC"),
                 max,
                 na.rm=T) %>% 
    # store Year and Treatment id
    mutate(Year.Trt = j) %>%
    # also store Year and Treatment in individual columns
    separate_wider_delim(Year.Trt, ".", names=c("Year", "Trt")) %>% 
    ungroup()
  
  
  datout.ChillLoop <- rbind(datout.ChillLoop, 
                            sub2)
  
}





#### Figure 1 ####

# Title: Simplified conceptual method of evaluation for chill accumulation models though budbreak assays
# conceptual figure

t=seq(0,30,0.1)
BB1=10+(35/(exp(0.05*t)))
BB2=10+(35/(exp(0.2*t)))

datBB=data.frame(t,BB=BB1,Treat=5)
datBB2=data.frame(t,BB=BB2,Treat=12)

Temp=seq(-5,30,0.1)
modA=-0.0152*Temp*Temp+0.2437*Temp+0.0091
modA[modA<0]=NA

modB=-0.0152*Temp*Temp+0.3656*Temp-1.2095
modB[modB<0]=NA

datBB3=rbind(datBB,datBB2)

datBB3$thermtA=(-0.0152*datBB3$Treat*datBB3$Treat+0.2437*datBB3$Treat+0.0091)*datBB3$t
datBB3$thermtB=(-0.0152*datBB3$Treat*datBB3$Treat+0.365*datBB3$Treat-1.2095)*datBB3$t


Fig1A <- ggplot()+
  geom_line(data=datBB3, aes(x=t,y=BB, col=as.factor(Treat)), lwd=0.8)+
  scale_color_manual(values=c("#588398","#CF223B"))+
  scale_y_continuous(limits=c(0,47), expand = c(0,0),breaks=F)+
  scale_x_continuous(limits=c(0,32), expand = c(0,0),breaks=F)+
  xlab("Time under treatment")+ylab("Time to budbreak")+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_rect(color="#CCCCCC"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_blank(),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))

Fig1B <- ggplot()+
  geom_line(aes(x=Temp,y=modA), lty=3, lwd=0.3)+
  geom_line(aes(x=Temp,y=modB),lwd=0.3)+
  geom_point(aes(x=c(5,12),y=c(-0.0152*5*5+0.2437*5+0.0091,-0.0152*12*12+0.2437*12+0.0091)), col=c("#588398","#CF223B"), size=2)+
  geom_point(aes(x=c(5,12),y=c(-0.0152*5*5+0.3656*5-1.2095,-0.0152*12*12+0.3656*12-1.2095)), col=c("#588398","#CF223B"), size=2)+
  scale_x_continuous(limits=c(0,21.5), expand = c(0,0),breaks=F)+
  scale_y_continuous(limits=c(0,1.1),expand = c(0,0),breaks=F)+
  xlab("Temperature")+ylab("Efficiency")+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_rect(color="#CCCCCC"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_blank(),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))

Fig1C <- ggplot()+
  geom_line(data=datBB3, aes(x=thermtA,y=BB, col=as.factor(Treat)), lwd=0.8)+
  scale_color_manual(values=c("#588398","#CF223B"))+
  scale_x_continuous(limits=c(0,27), expand = c(0,0),breaks = F)+
  scale_y_continuous(limits=c(0,47), expand = c(0,0),breaks=F)+
  xlab("Thermal time A")+ylab("Time to budbreak")+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_rect(color="#CCCCCC"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_blank(),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))

Fig1D <- ggplot()+
  geom_line(data=datBB3, aes(x=thermtB,y=BB, col=as.factor(Treat)), lwd=0.8)+
  geom_line(data=subset(datBB3, Treat==5), aes(x=thermtB,y=BB, col=as.factor(Treat)), lwd=0.8)+
  scale_color_manual(values=c("#588398","#CF223B"))+
  scale_x_continuous(limits=c(0,32), expand = c(0,0),breaks = F)+
  scale_y_continuous(limits=c(0,47), expand = c(0,0),breaks=F)+
  xlab("Thermal time B")+ylab("Time to budbreak")+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_rect(color="#CCCCCC"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_blank(),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))


Fig1top <- Fig1A + plot_spacer() + Fig1B + plot_layout(width=c(1,0.03,1))
Fig1bot <- Fig1C + plot_spacer() + Fig1D + plot_layout(width=c(1,0.03,1))

Fig1 <- Fig1top / plot_spacer() / Fig1bot + plot_layout(height=c(1,0.03,1)) #Saved as 3.5 x 3
Fig1


#### Figure 2 ####

# Title: Potential effects of chilling treatments on cold hardiness and time to budbreak
# conceptual figure

deaccrate=0.6

t2=seq(0,30,0.1)
CH1=-25+(22.5/(exp(0.1*t2)))
CH2=-20+(17.5/(exp(0.043*t2)))

t2B=seq(30.1, -(((-25+(22.5/(exp(0.1*30)))+2.5)/deaccrate)-30) ,0.5)
t2B2=seq(30.1, -(((-20+(17.5/(exp(0.043*30)))+2.5)/deaccrate)-30) ,0.5)


CH1B=-25+(22.5/(exp(0.1*30)))+deaccrate*(t2B-30)
CH2B=-20+(17.5/(exp(0.043*30)))+deaccrate*(t2B2-30)



datCH=data.frame(t=t2,CH=CH1,Treat=5,Phase=as.factor("Acc"))
datCHB=data.frame(t=t2B,CH=CH1B, Treat=5,Phase=as.factor("Deacc"))

datCH2=data.frame(t=t2,CH=CH2,Treat=12,Phase=as.factor("Acc"))
datCH2B=data.frame(t=t2B2,CH=CH2B, Treat=12,Phase=as.factor("Deacc"))


datCH3=rbind(datCH,datCH2)
datCH3=rbind(datCH3,datCHB)
datCH3=rbind(datCH3,datCH2B)


Fig2AB=ggplot()+
  geom_ribbon(aes(x=c(0,30),ymax=c(-2,-2),ymin=c(-27,-27)),lty=1, fill="#074783", lwd=0.1, alpha=0.15 )+
  geom_ribbon(aes(x=c(30,75),ymax=c(-2,-2),ymin=c(-27,-27)),lty=1, fill="#F58216", lwd=0.1, alpha=0.15 )+
  geom_line(data=datCH3, aes(x=t,y=CH, col=as.factor(Treat),lty=Phase),lwd=0.8)+
  geom_point(aes(x=30,y=c((-25+(22.5/(exp(0.1*30)))),(-20+(17.5/(exp(0.043*30))))) ), col=c("#588398","#CF223B"), size=2.5)+
  scale_color_manual(values=c("#588398","#CF223B"))+
  scale_y_continuous(limits=c(-27,-2),expand = c(0,0),breaks=c(-27))+
  scale_x_continuous(limits=c(0,75), expand = c(0,0),breaks=F)+
  xlab("Time")+ylab("Cold hardiness")+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_rect(color="#CCCCCC"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_blank(),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))


Fig2 <- Fig2AB/Fig2AB #Saved as 3.45 x 3
Fig2


#### Figure 3 ####

# Title: Cold hardiness changes of grapevine buds during chilling

# figure dataframes
datFig3 <- dat.BB50

# creates subset dataframes
# Trt:Year for plotting geom_points separately
for (i in levels(datFig3$Trt)) {
  datFig3.b <- datFig3 %>% filter(Trt == "field" | Trt %in% c("constant", "fluctuating") & Year == "2021_22")
  dfName <- paste0("datFig3.",i)
  assign(dfName, subset(datFig3.b, Trt==i))
}

# color palette for treatments
pal <- c("#e64e3c", "#78B7C5", "#f5d02c") #constant, #field, #fluctuating

Fig3 <- ggplot(data=datFig3, aes(x=days/7, y=CH_start, color=Trt))+
  geom_abline(slope=0, intercept = -10.6, lty="dashed", alpha=0.5)+
  geom_line(data=datFig3.field, aes(x=days/7, y=CH_start, group=Year), color="#78B7C5", size=0.25)+
  geom_smooth(se=F, size=0.8)+
  geom_point(data=datFig3.field, aes(x=days/7, y=CH_start, group=Year), color="#78B7C5", size=1.5)+
  geom_point(data=datFig3.constant, aes(x=days/7, y=CH_start), color="#e64e3c", size=1.5)+
  geom_point(data=datFig3.fluctuating, aes(x=days/7, y=CH_start), color="#f5d02c", size=1.5)+
  labs(x="Chilling exposure (weeks)", y="Cold hardiness (°C)", color="Chill\ntreatment", fill="Chill\ntreatment")+ 
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  scale_y_continuous(limits=c(-35, -8), breaks=c(-35,-30, -25, -20, -15, -10), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.25, 27), breaks=c(-0.25, 2, 4, 6, 10, 14, 20, 26), labels=c(0, 2, 4, 6, 10, 14, 20, 26), expand=c(0,0))+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_rect(color="#CCCCCC"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))

Fig3

# export width=3.5" & height=3"


#### Figure 4 ####

# Title: Chill accumulation based on different chill models for three chill treatments


datFig4 <- datout.ChillLoop


datFig4.long <- pivot_longer(datFig4, cols=c("final.chill", "portions", "UT", "NC"), names_to="model", values_to="chill")


datFig4.long <- datFig4.long %>% 
  mutate(model = ordered(model, c("NC", "UT", "portions", "final.chill")))


dat.trt.limits <- dat.budbreak.ori3 %>%
  filter(Year == "2021_22") %>%
  group_by(Trt) %>%
  summarise(trt.limit = max(days))


datFig4.long2 <- datFig4.long %>% 
  filter(Year == "2021_22") %>% 
  mutate(model.name = recode(model, 
                             NC = "NC",
                             UT = "UT",
                             portions = "DM",
                             final.chill = "New")) %>% 
  mutate(initiation.date = ymd("2021-10-13")) %>% 
  mutate(days.in.trt = ifelse(Start.date<initiation.date, 0, Start.date-initiation.date)) %>% 
  left_join(dat.trt.limits, by="Trt") %>% 
  mutate(chill.limited = ifelse(days.in.trt < trt.limit, chill, NA))


# color palette for treatments
pal <- c("#e64e3c", "#78B7C5", "#f5d02c") #constant, #field, #fluctuating


# Chill accumulation vs treatment duration 
# y-axis scaled according to North Carolina and Utah chill models
ggplot(data=subset(datFig4.long2, model.name %in% c("NC", "UT", "DM")), aes(x=days.in.trt/7, y=chill.limited, color=Trt))+
  geom_line(size=0.8)+
  scale_x_continuous(limits=c(0, 26), breaks=c(0, 2, 4, 6, 10, 14, 20, 26), labels=c(0, 2, 4, 6, 10, 14, 20, 26), expand=c(0,0))+ #weeks
  scale_y_continuous(limits=c(0,2550), n.breaks=6, expand=c(0,0))+
  labs(x="Time under treatment (weeks)", y="Chill accumulation (units per model)")+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))+
  ggh4x::facet_grid2(.~model, scales = "free", independent="y")


# Chill accumulation vs treatment duration 
# y-axis scaled according to Dynamic chill model
ggplot(data=subset(datFig4.long2, model.name %in% c("NC", "UT", "DM")), aes(x=days.in.trt/7, y=chill.limited, color=Trt))+
  geom_line(size=0.8)+
  scale_x_continuous(limits=c(0, 26), breaks=c(0, 2, 4, 6, 10, 14, 20, 26), labels=c(0, 2, 4, 6, 10, 14, 20, 26), expand=c(0,0))+ #weeks
  scale_y_continuous(limits=c(0,85), n.breaks=6, expand=c(0,0))+
  labs(x="Time under treatment (weeks)", y="Chill accumulation (units per model)")+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))+
  ggh4x::facet_grid2(.~model.name, scales = "free", independent="y")


# we exported two plots above for use in powerpoint to manually crop and combine the three panels



# sample plots combined with gridExtra::grid.arrange

Fig4_NC <- ggplot(data=subset(datFig4.long2, model.name %in% c("NC")), aes(x=days.in.trt/7, y=chill.limited, color=Trt))+
  geom_line(size=0.8)+
  scale_x_continuous(limits=c(0, 26), breaks=c(0, 2, 4, 6, 10, 14, 20, 26), labels=c(0, 2, 4, 6, 10, 14, 20, 26), expand=c(0,0))+ #weeks
  scale_y_continuous(limits=c(0,2550), n.breaks=6, expand=c(0,0))+
  labs(x="Time under treatment (weeks)", y="Chill accumulation (units per model)")+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))+
  facet_grid(.~model.name)

Fig4_UT <- ggplot(data=subset(datFig4.long2, model.name %in% c("UT")), aes(x=days.in.trt/7, y=chill.limited, color=Trt))+
  geom_line(size=0.8)+
  scale_x_continuous(limits=c(0, 26), breaks=c(0, 2, 4, 6, 10, 14, 20, 26), labels=c(0, 2, 4, 6, 10, 14, 20, 26), expand=c(0,0))+ #weeks
  scale_y_continuous(limits=c(0,2550), n.breaks=6, expand=c(0,0))+
  labs(x="Time under treatment (weeks)", y=element_blank())+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))+
  facet_grid(.~model.name)

Fig4_DM <- ggplot(data=subset(datFig4.long2, model.name %in% c("DM")), aes(x=days.in.trt/7, y=chill.limited, color=Trt))+
  geom_line(size=0.8)+
  scale_x_continuous(limits=c(0, 26), breaks=c(0, 2, 4, 6, 10, 14, 20, 26), labels=c(0, 2, 4, 6, 10, 14, 20, 26), expand=c(0,0))+ #weeks
  scale_y_continuous(limits=c(0,85), n.breaks=6, expand=c(0,0))+
  labs(x="Time under treatment (weeks)", y=element_blank())+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))+
  facet_grid(.~model.name)

Fig4 <- Fig4_NC + plot_spacer() + Fig4_UT + plot_spacer() + Fig4_DM + plot_layout(width=c(1,0.03,1,0.03,1))
Fig4


#### Figure 5 ####

# Title: Effect of increasing time under chill treatment exposure on time to budbreak

datFig5 <- dat.BB50

Fig5A <- ggplot(data=datFig5, aes(x=days/7, y=BB50, color=Trt))+
  geom_point(na.rm = T, size=1.5)+
  geom_smooth(se=F, size=0.5)+
  #geom_point(stat="summary", fun="mean", na.rm = T)+
  labs(x="Time under treatment\n(weeks)", y=expression(bold(underline("Observed"))~plain(time~to~budbreak~(days))), color="Chill\ntreatment", fill="Chill\ntreatment")+ 
  scale_y_continuous(limits=c(0, 70), breaks=c(0, 15, 30,45, 60, 80), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.25, 27), breaks=c(-0.25, 2, 4, 6, 10, 14, 20, 26), labels=c(0, 2, 4, 6, 10, 14, 20, 26), expand=c(0,0,0.1,0.1))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8)+
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))


Fig5B <- ggplot(data=datFig5, aes(x=days/7, y=adjBB50, color=Trt))+
  geom_point(na.rm = T, size=1.5)+
  geom_smooth(se=F, size=0.5)+
  #geom_point(stat="summary", fun="mean", na.rm = T)+
  scale_y_continuous(limits=c(0, 70), breaks=c(0, 15, 30,45, 60, 80), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.25, 27), breaks=c(-0.25, 2, 4, 6, 10, 14, 20, 26), labels=c(0, 2, 4, 6, 10, 14, 20, 26), expand=c(0,0,0.1,0.1))+
  labs(x="Time under treatment\n(weeks)", y=expression(bold(underline("Adjusted"))~plain(time~to~budbreak~(days))), color="Chill\ntreatment", fill="Chill\ntreatment")+ 
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8)+
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))

Fig5 <- Fig5A+Fig5B
Fig5



#### Figure 6 ####

# Title: Effect of chill accumulation on ***observed*** time to budbreak


datFig6 <- datout.ChillLoop


datFig6.long <- pivot_longer(datFig6, cols=c("final.chill", "portions", "UT", "NC"), names_to="model", values_to="chill")


datFig6.long <- datFig6.long %>% 
  mutate(model = ordered(model, c("NC", "UT", "portions", "final.chill"))) %>% 
  mutate(model.name = recode(model, NC = "NC", UT = "UT", portions = "DM", final.chill = "New")) %>% 
  mutate(Year = as.factor(Year),
         Trt = as.factor(Trt)) %>% 
  select(-Temp_C, -Tmax, -Tmin, -Tamp, -UT1, -NC1, -NC2)

# reduce to chill at collections
datFig6.long2 <- inner_join(dat.collections, datFig6.long, by=intersect(names(datFig6.long), names(dat.collections)))

# combine chill at collections with budbreak dataset
datFig6.long3 <- full_join(dat.budbreak.ori3, datFig6.long2, by=intersect(names(datFig6.long2), names(dat.budbreak.ori3)), relationship="many-to-many")



# color palette for treatments
pal <- c("#e64e3c", "#78B7C5", "#f5d02c") #constant, #field, #fluctuating


Fig6 <- ggplot(data=subset(datFig6.long3, model.name %in% c("NC", "UT", "DM")), aes(x=chill, y=budbreak, color=Trt))+
  geom_point(stat="summary", fun="mean", na.rm = T, size=1.5)+
  geom_smooth(data=subset(datFig6.long3, model.name %in% c("NC", "UT", "DM")),
              aes(x=chill, y=budbreak),
              method="nls", 
              formula=y ~
                b1 + (b2 / exp(b3 * x)), # this is an nls argument
              method.args = list(start=c(b1=10, b2=42, b3=0.001),
                                 control = list(maxiter = 50000, minFactor = 1/2000000, warnOnly = TRUE)), # this too
              se=FALSE,  fullrange=T,
              color="black",
              size=0.8)+
  labs(x="Chill accumulation (units based on model)", y=expression(bold(underline("Observed"))~plain(time~to~budbreak~(days))), color="Chill\ntreatment", fill="Chill\ntreatment")+ 
  scale_y_continuous(limits=c(0, 70), breaks=c(0, 15, 30,45, 60, 80), expand=c(0,0))+
  scale_x_continuous(limits=c(0,NA),expand=c(0,0,0.1,0.1))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8)+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))+
  ggh4x::facet_grid2(.~model.name, scales = "free", independent="all")

Fig6


## calculate RMSE for plot panels ##
# !!!!!
# Select the model for rmse calculations here
# !!!!!
datFig6.rmse <- datFig6.long3 %>% 
  filter(model.name == "NC") %>% #options: NC, UT, DM
  filter(!is.na(budbreak))

dat.rmse.BBobs <- aggregate(budbreak~model.name+chill+Trt+Year, data=datFig6.rmse, mean)

mod.BBobs <- nls(budbreak ~
                       b1 + (b2 / exp(b3 * chill)),
                     data=dat.rmse.BBobs,
                     start = list(b1=10, b2=42, b3=0.001),
                     algorithm = "port",
                     control = list(maxiter = 50000, minFactor = 1/2000000, warnOnly = TRUE))


# summarizing rmse based on previous model selections
dat.rmse.BBobs1 <- dat.rmse.BBobs %>% 
  mutate(predmod.BBobs = predict(mod.BBobs)) %>% 
  mutate(residuals.BBobs = budbreak - predmod.BBobs) %>% 
  summarise(model.name = unique(model.name),
            Trt = "pooled",
            rmse.BBobs = sqrt(sum(na.omit(residuals.BBobs)^2)/nrow(data.frame(na.omit(budbreak)))))

dat.rmse.BBobs2 <- dat.rmse.BBobs %>% 
  mutate(predmod.BBobs = predict(mod.BBobs)) %>% 
  mutate(residuals.BBobs = budbreak - predmod.BBobs) %>% 
  filter(Trt == "field") %>% 
  summarise(model.name = unique(model.name),
            Trt = unique(Trt),
            rmse.BBobs = sqrt(sum(na.omit(residuals.BBobs)^2)/nrow(.)))

dat.rmse.BBobs3 <- dat.rmse.BBobs %>% 
  mutate(predmod.BBobs = predict(mod.BBobs)) %>% 
  mutate(residuals.BBobs = budbreak - predmod.BBobs) %>% 
  filter(Trt == "fluctuating") %>% 
  summarise(model.name = unique(model.name),
            Trt = unique(Trt),
            rmse.BBobs = sqrt(sum(na.omit(residuals.BBobs)^2)/nrow(.)))

dat.rmse.BBobs4 <- dat.rmse.BBobs %>% 
  mutate(predmod.BBobs = predict(mod.BBobs)) %>% 
  mutate(residuals.BBobs = budbreak - predmod.BBobs) %>% 
  filter(Trt == "constant") %>% 
  summarise(model.name = unique(model.name),
            Trt = unique(Trt),
            rmse.BBobs = sqrt(sum(na.omit(residuals.BBobs)^2)/nrow(.)))


# observed budbreak rmse
BBobs.rmse.table <- as.data.frame(rbind(dat.rmse.BBobs1, dat.rmse.BBobs2, dat.rmse.BBobs3, dat.rmse.BBobs4))
BBobs.rmse.table






#### Figure 7 ####

# Title: Effect of chill accumulation on ***cold hardiness-adjusted*** time to budbreak

datFig7 <- datout.ChillLoop

datFig7.long <- pivot_longer(datFig7, cols=c("final.chill", "portions", "UT", "NC"), names_to="model", values_to="chill")


datFig7.long <- datFig7.long %>% 
  mutate(model = ordered(model, c("NC", "UT", "portions", "final.chill"))) %>% 
  mutate(model.name = recode(model, NC = "NC", UT = "UT", portions = "DM", final.chill = "New")) %>% 
  mutate(Year = as.factor(Year),
         Trt = as.factor(Trt)) %>% 
  select(-Temp_C, -Tmax, -Tmin, -Tamp, -UT1, -NC1, -NC2)

# reduce to chill at collections
datFig7.long2 <- inner_join(dat.collections, datFig7.long, by=intersect(names(datFig7.long), names(dat.collections)))

# combine chill at collections with budbreak dataset
datFig7.long3 <- full_join(dat.budbreak.ori3, datFig7.long2, by=intersect(names(datFig7.long2), names(dat.budbreak.ori3)), relationship="many-to-many")



# color palette for treatments
pal <- c("#e64e3c", "#78B7C5", "#f5d02c") #constant, #field, #fluctuating


Fig7 <- ggplot(data=subset(datFig7.long3, model.name %in% c("NC", "UT", "DM")), aes(x=chill, y=BB_adj, color=Trt))+
  geom_point(stat="summary", fun="mean", na.rm = T, size=1.5)+
  geom_smooth(data=subset(datFig7.long3, model.name %in% c("NC", "UT", "DM")),
              aes(x=chill, y=BB_adj),
              method="nls", 
              formula=y ~
                b1 + (b2 / exp(b3 * x)), # this is an nls argument
              method.args = list(start=c(b1=10, b2=42, b3=0.001),
                                 control = list(maxiter = 50000, minFactor = 1/2000000, warnOnly = TRUE)), # this too
              se=FALSE,  fullrange=T,
              color="black",
              size=0.8)+
  labs(x="Chill accumulation (units based on model)", y=expression(bold(underline("Adjusted"))~plain(time~to~budbreak~(days))), color="Chill\ntreatment", fill="Chill\ntreatment")+ 
  scale_y_continuous(limits=c(0, 70), breaks=c(0, 15, 30,45, 60, 80), expand=c(0,0))+
  scale_x_continuous(limits=c(0,NA),expand=c(0,0,0.1,0.1))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8)+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))+
  ggh4x::facet_grid2(.~model.name, scales = "free", independent="all")

Fig7


## calculate RMSE for plot panels ##
# !!!!!
# Select the model for rmse calculations here
# !!!!!
datFig7.rmse <- datFig7.long3 %>% 
  filter(model.name == "NC") %>% #options: NC, UT, DM
  filter(!is.na(budbreak))

dat.rmse.BBadj <- aggregate(BB_adj~model.name+chill+Trt+Year, data=datFig7.rmse, mean)

mod.BBadj <- nls(BB_adj ~
                   b1 + (b2 / exp(b3 * chill)),
                 data=dat.rmse.BBadj,
                 start = list(b1=10, b2=42, b3=0.001),
                 algorithm = "port",
                 control = list(maxiter = 50000, minFactor = 1/2000000, warnOnly = TRUE))


# summarizing rmse based on previous model selections
dat.rmse.BBadj1 <- dat.rmse.BBadj %>% 
  mutate(predmod.BBadj = predict(mod.BBadj)) %>% 
  mutate(residuals.BBadj = BB_adj - predmod.BBadj) %>% 
  summarise(model.name = unique(model.name),
            Trt = "pooled",
            rmse.BBadj = sqrt(sum(na.omit(residuals.BBadj)^2)/nrow(data.frame(na.omit(BB_adj)))))

dat.rmse.BBadj2 <- dat.rmse.BBadj %>% 
  mutate(predmod.BBadj = predict(mod.BBadj)) %>% 
  mutate(residuals.BBadj = BB_adj - predmod.BBadj) %>% 
  filter(Trt == "field") %>% 
  summarise(model.name = unique(model.name),
            Trt = unique(Trt),
            rmse.BBadj = sqrt(sum(na.omit(residuals.BBadj)^2)/nrow(data.frame(na.omit(BB_adj)))))

dat.rmse.BBadj3 <- dat.rmse.BBadj %>% 
  mutate(predmod.BBadj = predict(mod.BBadj)) %>% 
  mutate(residuals.BBadj = BB_adj - predmod.BBadj) %>% 
  filter(Trt == "fluctuating") %>% 
  summarise(model.name = unique(model.name),
            Trt = unique(Trt),
            rmse.BBadj = sqrt(sum(na.omit(residuals.BBadj)^2)/nrow(data.frame(na.omit(BB_adj)))))

dat.rmse.BBadj4 <- dat.rmse.BBadj %>% 
  mutate(predmod.BBadj = predict(mod.BBadj)) %>% 
  mutate(residuals.BBadj = BB_adj - predmod.BBadj) %>% 
  filter(Trt == "constant") %>% 
  summarise(model.name = unique(model.name),
            Trt = unique(Trt),
            rmse.BBadj = sqrt(sum(na.omit(residuals.BBadj)^2)/nrow(data.frame(na.omit(BB_adj)))))

# adjusted budbreak rmse
BBadj.rmse.table <- as.data.frame(rbind(dat.rmse.BBadj1, dat.rmse.BBadj2, dat.rmse.BBadj3, dat.rmse.BBadj4))
BBadj.rmse.table



#### Figure 8 ####

# Title: Efficiency of chill accumulation as a function of temperature for different chill models
# conceptual figure

datport=NULL
for (i in seq(-20,25,0.5)) {
  
  HourTemp=rep(i, 24*100)
  port=max(Dynamic_Model(HourTemp))
  Utah=sign(max(Utah_Model(HourTemp)))*max(abs(Utah_Model(HourTemp)))
  
  
  datport=rbind(datport,c(i,port,Utah))
  
}

datport=data.frame(datport)

colnames(datport)=c("Temperature","Portions","Utah")

datport$Utah=datport$Utah/max(datport$Utah)
datport$Portions=datport$Portions/max(datport$Portions)

data_NorthCarolina <- data.frame(NC_temperature = c(-1.1, 1.6, 7.2, 13, 16.5, 19, 20.7, 22.1, 23.3),
                                 NC_chill.units = c(0, 0.5, 1, 0.5, 0, -0.5, -1, -1.5, -2))

#create a lm with splines for Temp; as well as vector with predicted values
spline.mod <- lm(NC_chill.units ~ rcs(NC_temperature, c(-1.1, 1.6, 7.2, 13, 16.5, 19, 20.7, 22.1, 23.3)), data=data_NorthCarolina)
NC_spline <- data.frame(NC_temperature=seq(-1,23.3, by=0.1))
NC_spline <- NC_spline %>%
  mutate(spline.fit = predict(spline.mod, newdata = NC_spline))

kmax=1
Tmax=16
Topt=13
a=0.5

datport2=data.frame(Temperature=seq(-20,25,0.1))

datport2$NewMod=kmax*(((Tmax-datport2$Temperature)/(Tmax-Topt))^a)*exp(a*(1-((Tmax-datport2$Temperature)/(Tmax-Topt))))
datport2$NewMod[datport2$NewMod<0]=0

Fig8A <- ggplot()+
  geom_line(data=NC_spline, aes(x=NC_temperature, y=spline.fit), lty=2, lwd=0.4)+
  geom_line(data=datport, aes(x=Temperature, y=Portions), lwd=0.4)+
  geom_line(data=datport, aes(x=Temperature, y=Utah), lty=3, lwd=0.4)+
  geom_line(data=datport2, aes(x=Temperature, y=NewMod), lwd=0.4, col="#1B9E77")+
  
  scale_x_continuous(limits=c(-10,22), expand=c(0,0), breaks=c(-10,-5,0,5,10,15,20))+
  scale_y_continuous(limits=c(-1.25,1.1), expand=c(0,0), breaks=c(-0.25,0,0.25,0.5,0.75,1))+
  
  xlab("Temperature (°C)")+ylab("Efficiency")+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_rect(color="#CCCCCC"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))+
  coord_cartesian(xlim = c(-10,22), ylim = c(-0.25, 1.1), expand = TRUE, default = FALSE, clip = "on")


datport3=data.frame(amplit=seq(0,10,0.1))

Tampinfl=4
Tampslope=4
Tampmax=1.5

datport3$Enhance=1+(Tampmax-1)/(1+exp(Tampslope*(Tampinfl-datport3$amplit)))

Fig8B <- ggplot()+
  geom_line(data=datport3, aes(x=amplit, y=Enhance), lwd=0.4, col="#1B9E77")+
  
  scale_x_continuous(limits=c(0,7), expand=c(0,0), breaks=c(0,2,4,6))+
  scale_y_continuous(limits=c(1,1.55), expand=c(0,0), breaks=c(1,1.25,1.5))+
  
  xlab("Daily temperature \n amplitude (°C)")+ylab("Enhancement")+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_rect(color="#CCCCCC"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))

Fig8AB = Fig8A + plot_spacer() + Fig8B + plot_layout(widths=c(1,0.03,0.3))

Fig8 <- Fig8AB/Fig8AB #Save at 3.5 x 4 
Fig8


#### Figure 9 ####

# Title: Chill accumulation and effects on time to budbreak of a new chill model

datFig9 <- datout.ChillLoop


datFig9.long <- pivot_longer(datFig9, cols=c("final.chill", "portions", "UT", "NC"), names_to="model", values_to="chill")

datFig9.long <- datFig9.long %>% 
  mutate(model = ordered(model, c("NC", "UT", "portions", "final.chill"))) %>% 
  mutate(model.name = recode(model, NC = "NC", UT = "UT", portions = "DM", final.chill = "New")) %>% 
  mutate(Year = as.factor(Year),
         Trt = as.factor(Trt)) %>% 
  select(-Temp_C, -Tmax, -Tmin, -Tamp, -UT1, -NC1, -NC2)


# reduce to chill at collections
datFig9.long2 <- inner_join(dat.collections, datFig9.long, by=intersect(names(datFig9.long), names(dat.collections)))


# combine chill at collections with budbreak dataset
datFig9.long3 <- full_join(dat.budbreak.ori3, datFig9.long2, by=intersect(names(datFig9.long2), names(dat.budbreak.ori3)), relationship="many-to-many")


datFig9.long <- datFig9.long %>% 
  mutate(model = ordered(model, c("NC", "UT", "portions", "final.chill")))


dat.trt.limits <- dat.budbreak.ori3 %>%
  filter(Year == "2021_22") %>%
  group_by(Trt) %>%
  summarise(trt.limit = max(days))

# for Fig9A
datFig9.longA <- datFig9.long %>% 
  filter(Year == "2021_22") %>% 
  mutate(model.name = recode(model, 
                             NC = "NC",
                             UT = "UT",
                             portions = "DM",
                             final.chill = "New")) %>% 
  mutate(initiation.date = ymd("2021-10-13")) %>% 
  mutate(days.in.trt = ifelse(Start.date<initiation.date, 0, Start.date-initiation.date)) %>% 
  left_join(dat.trt.limits, by="Trt") %>% 
  mutate(chill.limited = ifelse(days.in.trt < trt.limit, chill, NA))


# color palette for treatments
pal <- c("#e64e3c", "#78B7C5", "#f5d02c") #constant, #field, #fluctuating


# sample plots combined with gridExtra::grid.arrange

Fig9A <- ggplot(data=subset(datFig9.longA, model.name %in% c("New")), aes(x=days.in.trt/7, y=chill.limited, color=Trt))+
  geom_line(size=0.8)+
  scale_x_continuous(limits=c(0, 26), breaks=c(0, 2, 4, 6, 10, 14, 20, 26), labels=c(0, 2, 4, 6, 10, 14, 20, 26), expand=c(0,0))+ #weeks
  scale_y_continuous(limits=c(0,2550), n.breaks=6, expand=c(0,0))+
  labs(x="Time under treatment (weeks)", y="Chill accumulation (units per model)")+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))+
  facet_grid(.~model.name)

Fig9B <- ggplot(data=subset(datFig6.long3, model.name %in% c("New")), aes(x=chill, y=budbreak, color=Trt))+
  geom_point(stat="summary", fun="mean", na.rm = T, size=1.5)+
  geom_smooth(data=subset(datFig6.long3, model.name %in% c("New")),
              aes(x=chill, y=budbreak),
              method="nls", 
              formula=y ~
                b1 + (b2 / exp(b3 * x)), # this is an nls argument
              method.args = list(start=c(b1=10, b2=42, b3=0.001),
                                 control = list(maxiter = 50000, minFactor = 1/2000000, warnOnly = TRUE)), # this too
              se=FALSE,  fullrange=T,
              color="black",
              size=0.8)+
  labs(x="Chill accumulation\n(chill units)", y=expression(bold(underline("Observed"))~plain(time~to~budbreak~(days))), color="Chill\ntreatment", fill="Chill\ntreatment")+ 
  scale_y_continuous(limits=c(0, 70), breaks=c(0, 15, 30,45, 60, 80), expand=c(0,0))+
  scale_x_continuous(limits=c(0,NA),expand=c(0,0,0.1,0.1))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))+
  ggh4x::facet_grid2(.~model.name, scales = "free", independent="all")

Fig9C <- ggplot(data=subset(datFig7.long3, model.name %in% c("New")), aes(x=chill, y=BB_adj, color=Trt))+
  geom_point(stat="summary", fun="mean", na.rm = T, size=1.5)+
  geom_smooth(data=subset(datFig7.long3, model.name %in% c("New")),
              aes(x=chill, y=BB_adj),
              method="nls", 
              formula=y ~
                b1 + (b2 / exp(b3 * x)), # this is an nls argument
              method.args = list(start=c(b1=10, b2=42, b3=0.001),
                                 control = list(maxiter = 50000, minFactor = 1/2000000, warnOnly = TRUE)), # this too
              se=FALSE,  fullrange=T,
              color="black",
              size=0.8)+
  labs(x="Chill accumulation\n(chill units)", y=expression(bold(underline("Adjusted"))~plain(time~to~budbreak~(days))), color="Chill\ntreatment", fill="Chill\ntreatment")+ 
  scale_y_continuous(limits=c(0, 70), breaks=c(0, 15, 30,45, 60, 80), expand=c(0,0))+
  scale_x_continuous(limits=c(0,NA),expand=c(0,0,0.1,0.1))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))+
  ggh4x::facet_grid2(.~model.name, scales = "free", independent="all")


# print the figure with all 3 panels
Fig9 <- Fig9A+Fig9B+Fig9C 
Fig9



## calculate RMSE for plot panels ##
datFig9.rmse <- datFig9.long3 %>% 
  filter(model.name == "New") %>%
  filter(!is.na(budbreak))

datFig9.rmse.BBobs <- aggregate(budbreak~model.name+chill+Trt+Year, data=datFig6.rmse, mean)
modFig9.BBobs <- nls(budbreak ~
                   b1 + (b2 / exp(b3 * chill)),
                 data=datFig9.rmse.BBobs,
                 start = list(b1=10, b2=42, b3=0.001),
                 algorithm = "port",
                 control = list(maxiter = 50000, minFactor = 1/2000000, warnOnly = TRUE))

datFig9.rmse.BBadj <- aggregate(BB_adj~model.name+chill+Trt+Year, data=datFig7.rmse, mean)
modFig9.BBadj <- nls(BB_adj ~
                   b1 + (b2 / exp(b3 * chill)),
                 data=datFig9.rmse.BBadj,
                 start = list(b1=10, b2=42, b3=0.001),
                 algorithm = "port",
                 control = list(maxiter = 50000, minFactor = 1/2000000, warnOnly = TRUE))


# summarizing rmse based on previous model selections
dat.rmse.BBobs9.1 <- dat.rmse.BBobs %>% 
  mutate(predmod.BBobs = predict(modFig9.BBobs)) %>% 
  mutate(residuals.BBobs = budbreak - predmod.BBobs) %>% 
  summarise(model.name = unique(model.name),
            Trt = "pooled",
            rmse.BBobs = sqrt(sum(na.omit(residuals.BBobs)^2)/nrow(data.frame(na.omit(budbreak)))))

dat.rmse.BBobs9.2 <- dat.rmse.BBobs %>% 
  mutate(predmod.BBobs = predict(modFig9.BBobs)) %>% 
  mutate(residuals.BBobs = budbreak - predmod.BBobs) %>% 
  filter(Trt == "field") %>% 
  summarise(model.name = unique(model.name),
            Trt = unique(Trt),
            rmse.BBobs = sqrt(sum(na.omit(residuals.BBobs)^2)/nrow(.)))

dat.rmse.BBobs9.3 <- dat.rmse.BBobs %>% 
  mutate(predmod.BBobs = predict(modFig9.BBobs)) %>% 
  mutate(residuals.BBobs = budbreak - predmod.BBobs) %>% 
  filter(Trt == "fluctuating") %>% 
  summarise(model.name = unique(model.name),
            Trt = unique(Trt),
            rmse.BBobs = sqrt(sum(na.omit(residuals.BBobs)^2)/nrow(.)))

dat.rmse.BBobs9.4 <- dat.rmse.BBobs %>% 
  mutate(predmod.BBobs = predict(modFig9.BBobs)) %>% 
  mutate(residuals.BBobs = budbreak - predmod.BBobs) %>% 
  filter(Trt == "constant") %>% 
  summarise(model.name = unique(model.name),
            Trt = unique(Trt),
            rmse.BBobs = sqrt(sum(na.omit(residuals.BBobs)^2)/nrow(.)))

dat.rmse.BBadj9.1 <- dat.rmse.BBadj %>% 
  mutate(predmod.BBadj = predict(modFig9.BBadj)) %>% 
  mutate(residuals.BBadj = BB_adj - predmod.BBadj) %>% 
  summarise(model.name = unique(model.name),
            Trt = "pooled",
            rmse.BBadj = sqrt(sum(na.omit(residuals.BBadj)^2)/nrow(data.frame(na.omit(BB_adj)))))

dat.rmse.BBadj9.2 <- dat.rmse.BBadj %>% 
  mutate(predmod.BBadj = predict(modFig9.BBadj)) %>% 
  mutate(residuals.BBadj = BB_adj - predmod.BBadj) %>% 
  filter(Trt == "field") %>% 
  summarise(model.name = unique(model.name),
            Trt = unique(Trt),
            rmse.BBadj = sqrt(sum(na.omit(residuals.BBadj)^2)/nrow(data.frame(na.omit(BB_adj)))))

dat.rmse.BBadj9.3 <- dat.rmse.BBadj %>% 
  mutate(predmod.BBadj = predict(modFig9.BBadj)) %>% 
  mutate(residuals.BBadj = BB_adj - predmod.BBadj) %>% 
  filter(Trt == "fluctuating") %>% 
  summarise(model.name = unique(model.name),
            Trt = unique(Trt),
            rmse.BBadj = sqrt(sum(na.omit(residuals.BBadj)^2)/nrow(data.frame(na.omit(BB_adj)))))

dat.rmse.BBadj9.4 <- dat.rmse.BBadj %>% 
  mutate(predmod.BBadj = predict(modFig9.BBadj)) %>% 
  mutate(residuals.BBadj = BB_adj - predmod.BBadj) %>% 
  filter(Trt == "constant") %>% 
  summarise(model.name = unique(model.name),
            Trt = unique(Trt),
            rmse.BBadj = sqrt(sum(na.omit(residuals.BBadj)^2)/nrow(data.frame(na.omit(BB_adj)))))


# observed budbreak rmse
BBobs.rmse.table <- as.data.frame(rbind(dat.rmse.BBobs9.1, dat.rmse.BBobs9.2, dat.rmse.BBobs9.3, dat.rmse.BBobs9.4))
BBobs.rmse.table
# adjusted budbreak rmse
BBadj.rmse.table <- as.data.frame(rbind(dat.rmse.BBadj9.1, dat.rmse.BBadj9.2, dat.rmse.BBadj9.3, dat.rmse.BBadj9.4))
BBadj.rmse.table



#### Supplemental Figure 1 ####

datSuppFig1 <- dat.temperature2 %>% 
  mutate(year.name = recode(Year,
                            "2021_22" = "2021-2022",
                            "2019_20" = "2019-2020",
                            "2018_19" = "2018-2019"))

# color palette for treatments
pal2 <- c("#78B7C5", "#e64e3c", "#f5d02c") #constant, #field, #fluctuating

SuppFig1 <- ggplot()+
  geom_ribbon(data=subset(datSuppFig1, Trt %in% "field"), aes(x=Start.date, ymax=Tmax, ymin=Tmin, fill=Trt, color=Trt))+
  geom_line(data=subset(datSuppFig1, Trt %in% c("constant","fluctuating")), aes(x=Start.date, y=Tmin, color=Trt, fill=Trt))+
  geom_line(data=subset(datSuppFig1, Trt %in% c("constant","fluctuating")), aes(x=Start.date, y=Tmax, color=Trt, fill=Trt))+
  scale_x_date(breaks="2 month", date_labels = "%b")+
  scale_y_continuous(breaks=c(-30, -20, -10, 0,  10, 20, 30))+
  labs(y="Temperature °C", x=element_blank())+
  scale_color_manual(values=pal2)+
  scale_fill_manual(values=pal2)+
  theme_bw(base_size=8) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size="8", face = "bold"),
        strip.text.y = element_text(size="8", face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))+
  ggh4x::facet_grid2(Trt~year.name, axes="all", scales = "free_x")

SuppFig1

#### Supplemental Figure 2 ####

# Title: Chill accumulation based on four different chill models for three chill treatments and three seasons

datSuppFig2 <- datout.ChillLoop


datSuppFig2.long <- pivot_longer(datSuppFig2, cols=c("final.chill", "portions", "UT", "NC"), names_to="model", values_to="chill")


datSuppFig2.long <- datSuppFig2.long %>% 
  mutate(model = ordered(model, c("NC", "UT", "portions", "final.chill")))


dat.trt.limits <- dat.budbreak.ori3 %>%
  group_by(Year, Trt) %>%
  summarise(trt.limit = max(days))


datSuppFig2.long2 <- datSuppFig2.long %>% 
  #select(-Start.date) %>% 
  mutate(model.name = recode(model, 
                             NC = "NC",
                             UT = "UT",
                             portions = "DM",
                             final.chill = "New")) %>% 
  mutate(year.name = recode(Year,
                            "2021_22" = "2021-2022",
                            "2019_20" = "2019-2020",
                            "2018_19" = "2018-2019")) %>% 
  mutate(initiation.date = ifelse(Year == "2021_22", ymd("2021-10-13"),
                                  ifelse(Year == "2019_20", ymd("2019-10-15"),
                                         ifelse(Year == "2018_19", ymd("2018-10-11"), ymd("2000-01-01"))))) %>% 
  mutate(days.in.trt = ifelse(Start.date<initiation.date, 0, Start.date-initiation.date)) %>% 
  left_join(dat.trt.limits, by=c("Trt","Year")) %>% 
  mutate(chill.limited = ifelse(days.in.trt < trt.limit, chill, NA)) %>% 
  mutate(chill.limited2 = ifelse(Start.date < initiation.date, NA, chill.limited))



# color palette for treatments
pal <- c("#e64e3c", "#78B7C5", "#f5d02c") #constant, #field, #fluctuating


# Chill accumulation vs treatment duration 
# y-axis scaled according to North Carolina and Utah chill models
SuppFig2 <- ggplot(data=subset(datSuppFig2.long2, model.name %in% c("NC", "UT", "DM", "New")), aes(x=days.in.trt/7, y=chill.limited2, color=Trt))+
  geom_line(size=0.8)+
  scale_x_continuous(limits=c(0, 26), breaks=c(0, 2, 4, 6, 10, 14, 20, 26), labels=c(0, 2, 4, 6, 10, 14, 20, 26), expand=c(0,0))+ #weeks
  scale_y_continuous(limits=c(0,NA) ,n.breaks=6, expand=c(0,0))+
  labs(x="Time under treatment (weeks)", y="Chill accumulation (units per model)")+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8) +
  theme(aspect.ratio = 1,
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size="8", face = "bold"),
        strip.text.y = element_text(size="8", face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))+
  ggh4x::facet_grid2(model.name~year.name, scales = "free", axes="all", independent="y")

SuppFig2



#### Supplemental Figure 3 ####

datSuppFig3 <- dat.BB50


SuppFig3A <- ggplot(data=datSuppFig3, aes(x=days/7, y=percBB, color=Trt))+
  geom_point(stat="summary", fun="mean", na.rm = T, size=1.5)+
  geom_line(stat="summary", fun="mean", size=0.25)+
  geom_smooth(size=0.8,
              method = "nls", formula = y ~ (100/(1+exp((b)*(x-c)))),
              method.args = list(start = list(b=0.2, c=6)), se = F, fullrange=F)+
  labs(x="Time under treatment (weeks)", y="Percent budbreak (%)", color="Chill\ntreatment", fill="Chill\ntreatment")+ 
  scale_y_continuous(limits=c(0, 108), breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.25, 27), breaks=c(-0.25, 2, 4, 6, 10, 14, 20, 26), labels=c(0, 2, 4, 6, 10, 14, 20, 26), expand=c(0,0,0.1,0.1))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8)+
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))

SuppFig3A

SuppFig3B <- ggplot(data=datSuppFig3, aes(x=days/7, y=BB50, color=Trt))+
  geom_line(data=datSuppFig3, aes(x=days/7, y=BBcensored, color=Trt), lty="dashed", size=.5)+
  geom_segment(data=datSuppFig3, aes(x=days/7, xend=days/7, y=BBcensored, yend=BBcensored2, color=Trt), lty="dashed", size=.5)+
  geom_ribbon(data=datSuppFig3, aes(x=days/7, ymax=BB80B, ymin=BB20B, fill=Trt), na.rm=T, size=.5, alpha=0.15)+
  labs(x="Time under treatment (weeks)", y="Variability in observed time to budbreak\n(20 to 80 percentiles; days)", color="Chill\ntreatment", fill="Chill\ntreatment")+ 
  scale_y_continuous(limits=c(0, 62), breaks=c(0, 15, 30,45, 60, 80), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.25, 27), breaks=c(-0.25, 2, 4, 6, 10, 14, 20, 26), labels=c(0, 2, 4, 6, 10, 14, 20, 26), expand=c(0,0,0.1,0.1))+
  scale_fill_manual(values=pal)+
  scale_color_manual(values=pal)+
  theme_bw(base_size=8)+
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size=0.25,color="black"),
        axis.line.y = element_line(size=0.25,color="black"),
        axis.text = element_text(size=8, color="black"),
        axis.ticks = element_line(size=0.25,color="black"),
        axis.ticks.length = unit(0.3, 'lines'))

SuppFig3B


SuppFig3 <- SuppFig3A / SuppFig3B
SuppFig3




#### print all figures ####

Fig1 #Saved as 3.5 x 3
Fig2 #Saved as 3.45 x 3
Fig3
Fig4
Fig5
Fig6
Fig7
Fig8 #Save at 3.5 x 4 
Fig9
SuppFig1
SuppFig2
SuppFig3


# <end>






