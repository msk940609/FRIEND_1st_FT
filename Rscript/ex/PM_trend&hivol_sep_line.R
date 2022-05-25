library(extrafont)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)

getOption("digits")
options("digits" = 15)

pmcn=fread(file = "Datafile/PM/pm_cn.csv")
pmcn

pmcn=pmcn %>% 
  mutate(date2 = as.POSIXct(pmcn$Date, format = '%Y-%m-%d %H:%M'))
pmcn


cn_sampling=as.data.table(pmcn)

cn_sampling=cn_sampling %>% separate(Date,c("Year","Month","day2"),sep = "-")
cn_sampling=cn_sampling %>% separate(day2,c("Day","hour"),sep = " ")
cn_sampling=cn_sampling %>% separate(hour,c("Hour","min"),sep = ":")

#cn_sampling=cn_sampling %>% filter(Hour!="10")


cn_sampling_m=as.data.table(melt(cn_sampling[,c("Year","Month","Day","Hour","PM2.5")], id.vars=c("Year","Month","Day","Hour")))
cn_sampling_m$Date=paste0(cn_sampling_m$Year,"-",cn_sampling_m$Month,"-",cn_sampling_m$Day, " ", cn_sampling_m$Hour,":00")


cn_sampling_d=dcast.data.table(cn_sampling_m, Date~variable, mean)
cn_sampling_d
#fwrite(cn_sampling_d, file = "cn_PM_hour.csv")

pm2.5=fread("Datafile/PM/PM2.5_all_merge.csv")

pm2.5=pm2.5 %>% 
  mutate(date2 = as.POSIXct(pm2.5$Date, format = '%Y-%m-%d %H'))
pm2.5
#pm2.5$Site=factor(pm2.5$Site,levels = c("Seoul","Seosan"))

##Basic option====
lims <- as.POSIXct(strptime(c("2020-12-15 10:00","2021-01-15 23:00"), format = "%Y-%m-%d %H:%M"))

###Day average=====
#fwrite(event_day_d, file="day_avg_PM2.5_conc.csv")


###Hi vol sampling time based avarage=====

#write.csv(hivol_sampling,file = "1st_timest.csv")

hivol_sampling=fread(file = "Datafile/PM/1st_timest.csv")
hivol_sampling=hivol_sampling %>% 
  mutate(date2 = as.POSIXct(hivol_sampling$Date, format = '%Y-%m-%d %H'))

hivol_sampling$Hidate=hivol_sampling$date2-36000
hivol_sampling=hivol_sampling %>% mutate(day= as.Date(Hidate, format = '%Y-%m-%d',  tz="Asia/seoul"))
hivol_sampling=hivol_sampling %>% separate(Date,c("Year","Month","day2"),sep = "-")
hivol_sampling=hivol_sampling %>% separate(day2,c("Day","Hour"),sep = " ")
hivol_sampling=hivol_sampling %>% filter(Hour!="9:00")
hivol_sampling_m=as.data.table(melt(hivol_sampling[,c("PM2.5","Site","day")], id.vars=c("Site","day")))
hivol_sampling_d=dcast.data.table(hivol_sampling_m, Site~day, mean) %>% melt.data.table(id.vars = "Site") %>% `colnames<-`(c("Site","day","PM2.5"))

####
sample_day_raw=fread(file ="Datafile/PM/sampling_avg_PM2.5_conc.csv")
#sample_day_raw=sample_day_raw %>% filter(day!="NA") %>% filter(PM2.5!="NaN")

sample_day_raw=sample_day_raw %>% mutate(date2=as.POSIXct(day, format = '%Y-%m-%d')+36000) %>% 
  mutate(date3=as.POSIXct(date2, format = '%Y-%m-%d')+82800) %>% 
  mutate(label=as.Date(day,format = '%Y-%m-%d')) %>% 
  separate(label, c("Year","Month","Day"), sep = "-") %>% 
  mutate(label=paste(Month,Day,sep = "/"))


pm2.5$Site=factor(pm2.5$Site,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))
sample_day_raw$Site=factor(sample_day_raw$Site,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))

lims <- as.POSIXct(strptime(c("2020-12-15 00:00","2021-01-16 0:00"), format = "%Y-%m-%d %H:%M"))

#sample_day_raw=sample_day_raw %>% filter(day!=c("2021-01-16"))
subset(sample_day_raw,sample_day_raw$Event=="Non event"|sample_day_raw$Event=="Event")
subset(sample_day_raw,sample_day_raw$Event!="Fail")

ggplot()+
  geom_rect(data=subset(sample_day_raw,sample_day_raw$Event=="Non-event"|sample_day_raw$Event=="Event"),
            aes(xmin=date2,xmax=date2+86400, ymax=Inf, ymin=-Inf, group=Site,fill=Event),alpha=0.4)+
  geom_vline(data=sample_day_raw, aes(xintercept=date2), lty=2, col="grey50", size=0.7)+
  geom_line(data = pm2.5 ,aes(x=date2+36000, y=PM2.5, col=Site),size=0.7, na.rm = T)+
  geom_text(data = sample_day_raw, aes(label=label, x=date2+42600, 
                                       y=ifelse(Site=="Ulaanbaatar",550*0.85, ##550*0.92
                                                ifelse(Site=="Beijing",160*0.85,
                                                       ifelse(Site=="Seoul",160*0.85,
                                                              ifelse(Site=="Seosan",190*0.85,
                                                                     50*0.85))))), size=4)+
  geom_text(data = sample_day_raw, aes(label=label, x=date2+42600, 
                                       y=ifelse(Site=="Ulaanbaatar",550,
                                                ifelse(Site=="Beijing",160,
                                                       ifelse(Site=="Seoul",160,
                                                              ifelse(Site=="Seosan",190,
                                                                     50))))), size=4, col=NA)+
  geom_text(data = sample_day_raw, aes(label=paste0(sprintf("%0.1f", round(PM2.5, digits = 1))), x=date2+42600, 
  #geom_text(data = subset(sample_day_raw,sample_day_raw$Event!="Fail"), aes(label=paste0(sprintf("%0.1f", round(PM2.5, digits = 1))), x=date2+42600, 
                                       y=ifelse(Site=="Ulaanbaatar",550*0.75,
                                                ifelse(Site=="Beijing",160*0.75,
                                                       ifelse(Site=="Seoul",160*0.75,
                                                              ifelse(Site=="Seosan",190*0.75,
                                                                     50*0.75))))), size=3.5, col="black")+
  #geom_text(data = subset(sample_day_raw,sample_day_raw$Event=="Fail"), aes(label=Event, x=date2+42600, 
  #                                                                          y=ifelse(Site=="Ulaanbaatar",550*0.75,
  #                                                                                   ifelse(Site=="Beijing",160*0.75,
  #                                                                                          ifelse(Site=="Seoul",160*0.75,
  #                                                                                                 ifelse(Site=="Seosan",190*0.75,
  #                                                                                                        50*0.75))))), size=5, col="black")+
  #geom_text(data = sample_day_raw, aes(label=paste0("(",sprintf("%0.1f", round(PM2.5, digits = 1)),")"), x=date2+42600, y=125), size=2.5)+
  #geom_hline(yintercept = 35, col="grey70", size=1.5, lty="dashed")+
  facet_wrap(.~Site, ncol = 1, scales = "free")+
  scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  scale_fill_manual(values = c("orangered1","grey25"))+
  scale_x_datetime('',limits = lims,
                   date_breaks = '1 days',
                   date_labels = "%m/%d", expand = c(0.008,0.008))+
  scale_y_continuous(name = expression(bold("Mass conc."~"("*"\u03bcg/"*m^"3"*")")),expand = c(0.02,0.02))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.text = element_text(colour = "black", size = 16, face = "bold",margin = unit(c(0.3,0.2,0.2,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.text.x = element_text(size = 11,angle = 0,colour = "black", face = "bold",family = "Arial", vjust = 0.5, hjust = 0.5),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial"),
        axis.title.y.left = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = c(0.85,1.02))+
  xlab("")+
  guides(fill=guide_legend(order = 1,title =bquote(PM[2.5]~"Episode"),col=c(NA,NA,NA), linetype=c(1,1,1), alpha=0.4),
         col=F)+
  ggsave("Kr_PM_hivol_all_1.tiff",height = 30, width = 40, units = "cm", dpi = 300)

sample_day_raw


pm2.5

pm2.5$new=ifelse(pm2.5$Site=="Ulaanbaatar",pm2.5$PM2.5*0.25,pm2.5$PM2.5)

ggplot()+
  geom_rect(data=subset(sample_day_raw,sample_day_raw$Event=="Non event"|sample_day_raw$Event=="Event"|sample_day_raw$Event=="Intermediate"),
            aes(xmin=date2,xmax=date2+86400, ymax=Inf, ymin=-Inf, group=Site,fill=Event),alpha=0.4)+
  geom_vline(data=sample_day_raw, aes(xintercept=date2), lty=2, col="grey50", size=0.7)+
  geom_line(data = pm2.5 ,aes(x=date2, y=new, col=Site),size=0.7, na.rm = T)+
  #geom_text(data = sample_day_raw, aes(label=label, x=date2+42600, 
  #                                     y=ifelse(Site=="Ulaanbaatar",450,
  #                                              ifelse(Site=="Beijing",115,
  #                                                     ifelse(Site=="Seoul",100,
  #                                                            ifelse(Site=="Seosan",150,
  #                                                                   41))))), size=4)+
  #geom_text(data = sample_day_raw, aes(label=label, x=date2+42600, 
  #                                     y=ifelse(Site=="Ulaanbaatar",520,
  #                                              ifelse(Site=="Beijing",140,
  #                                                     ifelse(Site=="Seoul",120,
  #                                                            ifelse(Site=="Seosan",180,
  #                                                                   45))))), size=4, col=NA)+
  #geom_text(data = subset(sample_day_raw,sample_day_raw$day!="2021-01-15"), aes(label=label, x=date2+42600, 
  #                                                                              y=ifelse(Site=="Ulaanbaatar",430,
  #                                                                                       ifelse(Site=="Beijing",115,
  #                                                                                              ifelse(Site=="Seoul",103,
  #                                                                                                     ifelse(Site=="Seosan",155,
  #                                                                                                            39))))), size=4, col="black")+
  #geom_text(data = sample_day_raw, aes(label=paste0("(",sprintf("%0.1f", round(PM2.5, digits = 1)),")"), x=date2+42600, y=125), size=2.5)+
  #geom_hline(yintercept = 35, col="grey70", size=1.5, lty="dashed")+
  facet_wrap(.~Site, ncol = 1, scales = "fixed")+
  scale_color_manual(values = c("black","#754F44","#009E73","#0072B2","#631879B2"))+
  scale_fill_manual(values = c("orangered1","grey25"))+
  scale_x_datetime('',limits = lims,
                   date_breaks = '1 days',
                   date_labels = "%m/%d", expand = c(0.008,0.008))+
  scale_y_continuous(name = expression(bold("Mass conc."~"("*"\u03bcg/"*m^"3"*")")),expand = c(0.02,0.02))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.text = element_text(colour = "black", size = 16, face = "bold",margin = unit(c(0.3,0.2,0.2,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.text.x = element_text(size = 11,angle = 0,colour = "black", face = "bold",family = "Arial", vjust = 0.5, hjust = 0.5),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial"),
        axis.title.y.left = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = c(0.79,1.02))+
  xlab("")+
  guides(fill=guide_legend(order = 1,title =bquote(PM[2.5]~"Episode"),col=c(NA,NA,NA), linetype=c(1,1,1), alpha=0.4),
         col=F)+
  ggsave("Kr_PM_hivol_all_2.tiff",height = 30, width = 40, units = "cm", dpi = 300)


library(ggsci)
mypal = pal_aaas("default", alpha = 0.7)(9)

library("scales")
show_col(mypal)


