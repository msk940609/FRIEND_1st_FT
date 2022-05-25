library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(writexl)
library(xlsx)
library(vegan)
options(java.parameters = "- Xmx8192m")
library(extrafont)
loadfonts(device="win")
library(ggsci)
library(tidyr)
library("scales")
library(multcompView)
library(DescTools)
library(lemon)
getOption("digits")
options("digits" = 15)

source("Rscript/func_filename.R")
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")

###environmental analysis====
ft_envi=fread("Datafile/FRIEND_1st_envi_re.csv")
ft_envi

ft_envi_oc=ft_envi[,c("Group","No","Event","PM2.5","OC","EC","OC/EC","WSOC","WISOC","HULIS-C","Non HULIS-C",
                      "WSOC/OC","HULIS-C/WSOC","RH","Temp")]

ft_envi_oc_m=melt(ft_envi_oc,id.vars = c("Group","No","Event"))

ft_envi_oc_m
table(ft_envi_oc_m$Group)
ft_envi_oc_m$id=factor(ft_envi_oc_m$Group, levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"),
                       labels = c("UL","BJ","SUL","SS","NT"))

ft_envi_oc_m$Event=factor(ft_envi_oc_m$Event, levels = c("Event","Normal","Non-event"))


ggplot(ft_envi_oc_m, aes(x=id,y=value, fill=Event))+
  geom_boxplot()+
  facet_rep_wrap(.~variable,scales = "free", ncol = 4,repeat.tick.labels = T, strip.position = "left")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=22, colour = "black",face = "bold",family = "Arial",
                                  margin = unit(c(0.2,0.4,0.2,0.2),"cm")),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,-0.2,0.1,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.18,0.95),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("Org_distribution"), height = 60, width = 85, units = "cm", dpi = 300)

#environmental distribution=====  
ft_envi_air=ft_envi[,c("Group","No","Event","PM2.5","NH4+","NO3-","SO42-","CO","O3","NO","SO2")]
ft_envi_air_m=melt(ft_envi_air,id.vars = c("Group","No","Event"))

ft_envi_air_m
table(ft_envi_air_m$Group)
ft_envi_air_m$id=factor(ft_envi_air_m$Group, levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"),
                       labels = c("UL","BJ","SUL","SS","NT"))

ft_envi_air_m$Event=factor(ft_envi_air_m$Event, levels = c("Event","Normal","Non-event"))

ggplot(ft_envi_air_m, aes(x=id,y=value, fill=Event))+
  geom_boxplot()+
  facet_rep_wrap(.~variable,scales = "free", ncol = 4,repeat.tick.labels = T, strip.position = "left")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=22, colour = "black",face = "bold",family = "Arial",
                                  margin = unit(c(0.2,0.4,0.2,0.2),"cm")),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,-0.2,0.1,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.18,0.95),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("envi_distribution"), height = 40, width = 85, units = "cm", dpi = 300)

#environmental distribution (OC norm)=====  
ft_envi_ocair=ft_envi[,c("Group","No","Event","PM2.5","OC","NH4+","NO3-","SO42-","CO","O3","NO","SO2")]
ft_envi_ocair_m=melt(ft_envi_ocair,id.vars = c("Group","No","Event","OC"))

ft_envi_air_m
table(ft_envi_air_m$Group)
ft_envi_ocair_m$id=factor(ft_envi_ocair_m$Group, levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"),
                        labels = c("UL","BJ","SUL","SS","NT"))

ft_envi_ocair_m$Event=factor(ft_envi_ocair_m$Event, levels = c("Event","Normal","Non-event"))

ft_envi_ocair_m$ocnorm=ft_envi_ocair_m$value/ft_envi_ocair_m$OC
table(ft_envi_ocair_m)
ft_envi_ocair_m$variable2=paste0(ft_envi_ocair_m$variable,"/OC")

table(ft_envi_ocair_m$variable2)
ft_envi_ocair_m$variable2=factor(ft_envi_ocair_m$variable2,
                                 levels=c("PM2.5/OC","NH4+/OC","NO3-/OC","SO42-/OC",
                                          "CO/OC","O3/OC","NO/OC","SO2/OC"))

ggplot(ft_envi_ocair_m, aes(x=id,y=ocnorm, fill=Event))+
  geom_boxplot()+
  facet_rep_wrap(.~variable2,scales = "free", ncol = 4,repeat.tick.labels = T, strip.position = "left")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=22, colour = "black",face = "bold",family = "Arial",
                                  margin = unit(c(0.2,0.4,0.2,0.2),"cm")),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,-0.2,0.1,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.18,0.95),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("envi_distribution_ocnorm_re"), height = 40, width = 85, units = "cm", dpi = 300)

##SOC======
ft_envi=fread("Datafile/FRIEND_1st_envi_re.csv")

FT_envi=ft_envi
FT_envi$WSOCbb=2.94*FT_envi$Levoglucosan
FT_envi$WSOCnbb=FT_envi$WSOC-FT_envi$WSOCbb

FT_envi_sel=FT_envi[,c("Sample","Group","No","Event","Date","PM2.5","OC","WSOC","WISOC","WSOCbb","WSOCnbb")]
FT_envi_sel$POC=FT_envi_sel$WISOC+FT_envi_sel$WSOCbb
FT_envi_sel$SOC=FT_envi_sel$WSOCnbb

FT_envi_sel$POCp=FT_envi_sel$POC/FT_envi_sel$OC*100
FT_envi_sel$SOCp=FT_envi_sel$SOC/FT_envi_sel$OC*100

FT_envi_sel$OC=ifelse(FT_envi_sel$SOC<0,FT_envi_sel$OC-2*FT_envi_sel$SOC,FT_envi_sel$OC)
FT_envi_sel$SOC=ifelse(FT_envi_sel$SOC<0,-FT_envi_sel$SOC,FT_envi_sel$SOC)
FT_envi_sel

FT_soa_p=melt(FT_envi_sel[,c("Sample","Group","No","Event","Date","PM2.5","POCp","SOCp")],
              id.vars = c("Sample","Group","No","Event","Date","PM2.5"))

FT_soa_p$Group=factor(FT_soa_p$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))

ggplot()+
  geom_boxplot(data=FT_soa_p,aes(x=Group,y=value, fill=variable))

FT_envi_sel

FT_envi_trend=melt(FT_envi_sel[,c("Group","No","Event","Date","PM2.5","OC","POC","SOC")],
                   id.vars = c("Group","No","Event","Date","PM2.5","OC"),variable.name = "type",value.name = "conc")

FT_envi_trend2=melt(FT_envi_trend,
                   id.vars = c("Group","No","Event","Date","type","conc"), na.rm = F)

FT_envi_trend2

FT_envi_trend=FT_envi_trend %>% 
  mutate(date2 = as.POSIXct(FT_envi_trend$Date, format = '%Y-%m-%d %H:%M'))
FT_envi_trend

FT_envi_trend2=FT_envi_trend2 %>% 
  mutate(date2 = as.POSIXct(FT_envi_trend2$Date, format = '%Y-%m-%d %H:%M'))
FT_envi_trend2

FT_envi_trend2$val2=ifelse(FT_envi_trend2$variable=="PM2.5",FT_envi_trend2$value*0.2,FT_envi_trend2$value)
lims <- as.POSIXct(strptime(c("2020-12-15 8:00","2021-01-14 12:00"), format = "%Y-%m-%d %H:%M"))

FT_envi_trend$type=factor(FT_envi_trend$type, levels = c("POC","SOC"))
FT_envi_trend$Group=factor(FT_envi_trend$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))
FT_envi_trend2$Group=factor(FT_envi_trend2$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))

max(FT_envi_trend2$value)

ggplot(data=FT_envi_trend2)+
  geom_area(data = FT_envi_trend ,aes(x=date2, y=conc,fill=type),position = position_stack(reverse = T),na.rm = T)+
  geom_line(data = FT_envi_trend2 ,aes(x=date2, y=val2, lty=variable),size=0.7, na.rm = T)+
  #geom_point(data = FT_envi_trend2 ,aes(x=date2, y=value*1.2),size=NA, na.rm = F,col=NA)+
  scale_y_continuous(name = "OC",expand = c(0.02,0.02),
                     sec.axis = sec_axis(~.*5, name="PM2.5")
                    )+
   facet_wrap(.~Group, ncol = 1, scales = "free")+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  scale_x_datetime('',limits = lims,
                   date_breaks = '2 days',
                   date_labels = "%m/%d", expand = c(0.008,0.008))+
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
  ggsave("Kr_PM_hivol_all_0.tiff",height = 30, width = 40, units = "cm", dpi = 300)

#trend=====
##Ul====
FT_envi_trend_ul=subset(FT_envi_trend,FT_envi_trend$Group=="Ulaanbaatar")
FT_envi_trend2_ul=subset(FT_envi_trend2,FT_envi_trend2$Group=="Ulaanbaatar")

FT_envi_trend2_ul$val2=ifelse(FT_envi_trend2_ul$variable=="OC",FT_envi_trend2_ul$value*1.5,FT_envi_trend2_ul$value)

FT_envi_trend2_ul

FT_envi_trend_ul$conc2=FT_envi_trend_ul$conc*1.5

FT_envi_trend_ul_rect=subset(FT_envi_trend_ul,is.na(FT_envi_trend_ul$conc))
FT_envi_trend2_ul

FT_envi_trend2_ul

FT_envi_trend3_ul=subset(FT_envi_trend2_ul,!is.na(FT_envi_trend2_ul$value))

ggplot(data=FT_envi_trend2_ul)+
  geom_area(data = FT_envi_trend_ul ,aes(x=date2, y=conc2,fill=type),position = position_stack(reverse = T),na.rm = T)+
  geom_rect(data=FT_envi_trend_ul_rect,aes(xmin=date2-43200,xmax=date2+43200, ymax=Inf, ymin=-Inf), fill="white")+
  geom_line(data = FT_envi_trend2_ul ,aes(x=date2, y=val2, lty=variable),size=0.7, na.rm = F)+
  geom_line(data = FT_envi_trend3_ul ,aes(x=date2, y=val2, lty=variable),size=0.7, na.rm = F)+
  geom_point(data = FT_envi_trend2_ul ,aes(x=date2, y=val2*1.1),size=NA, na.rm = F,col=NA)+
  scale_y_continuous(name =expression(bold(PM[2.5]~"("*"\u03bcg/"*m^"3"*")")),expand = c(0.02,0.02),
                     sec.axis = sec_axis(~./1.5, name=expression(bold("OC"~"("*"\u03bcg/"*m^"3"*")")))
  )+
  scale_linetype_manual(values=c(1,2),labels=c(expression(bold(PM[2.5])),expression(bold(OC))))+
  scale_fill_manual(values = c("grey70","#FD7F20"))+
  facet_wrap(.~Group, ncol = 1, scales = "free")+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  scale_x_datetime('',limits = lims,
                   date_breaks = '2 days',
                   date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.text = element_text(colour = "black", size = 16, face = "bold",margin = unit(c(0.3,0.2,0.2,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.text.x = element_text(size = 16,angle = 0,colour = "black", face = "bold",family = "Arial", vjust = 0.5, hjust = 0.5),
        axis.ticks.length = unit(0.15,"cm"),
        axis.ticks = element_line(size = 1.2, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial"),
        axis.title.y.left = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial"),
        axis.title.y.right = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.0,0.1,0.4),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.position = c(0.79,1.07))+
  xlab("")+
  guides(fill=guide_legend(order = 2,title =""),
         lty=guide_legend(order=1, title = ""))+
  ggsave(filename("Kr_PM_hivol_all_ul"), height = 10, width = 35, units = "cm", dpi = 300)

##bj====
FT_envi_trend_bj=subset(FT_envi_trend,FT_envi_trend$Group=="Beijing")
FT_envi_trend2_bj=subset(FT_envi_trend2,FT_envi_trend2$Group=="Beijing")

FT_envi_trend2_bj$val2=ifelse(FT_envi_trend2_bj$variable=="OC",FT_envi_trend2_bj$value*4,FT_envi_trend2_bj$value)
FT_envi_trend_bj$conc2=FT_envi_trend_bj$conc*4

FT_envi_trend_bj_rect=subset(FT_envi_trend_bj,is.na(FT_envi_trend_bj$conc))

FT_envi_trend3_bj=subset(FT_envi_trend2_bj,!is.na(FT_envi_trend2_bj$value))


ggplot(data=FT_envi_trend2_bj)+
  geom_area(data = FT_envi_trend_bj ,aes(x=date2, y=conc2,fill=type),position = position_stack(reverse = T),na.rm = T)+
  geom_rect(data=FT_envi_trend_bj_rect,aes(xmin=date2-43200,xmax=date2+43200, ymax=Inf, ymin=-Inf), fill="white")+
  geom_line(data = FT_envi_trend2_bj ,aes(x=date2, y=val2, lty=variable),size=0.7, na.rm = F)+
  geom_line(data = FT_envi_trend3_bj ,aes(x=date2, y=val2, lty=variable),size=0.7, na.rm = F)+
  geom_point(data = FT_envi_trend2_bj ,aes(x=date2, y=val2*1.1),size=NA, na.rm = F,col=NA)+
  scale_y_continuous(name =expression(bold(PM[2.5]~"("*"\u03bcg/"*m^"3"*")")),expand = c(0.02,0.02),
                     sec.axis = sec_axis(~./4, name=expression(bold("OC"~"("*"\u03bcg/"*m^"3"*")")))
  )+
  scale_linetype_manual(values=c(1,2),labels=c(expression(bold(PM[2.5])),expression(bold(OC))))+
  scale_fill_manual(values = c("grey70","#FD7F20"))+
  facet_wrap(.~Group, ncol = 1, scales = "free")+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  scale_x_datetime('',limits = lims,
                   date_breaks = '2 days',
                   date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.text = element_text(colour = "black", size = 16, face = "bold",margin = unit(c(0.3,0.2,0.2,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.text.x = element_text(size = 16,angle = 0,colour = "black", face = "bold",family = "Arial", vjust = 0.5, hjust = 0.5),
        axis.ticks.length = unit(0.15,"cm"),
        axis.ticks = element_line(size = 1.2, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial"),
        axis.title.y.left = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial"),
        axis.title.y.right = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.0,0.1,0.4),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.position = c(0.79,1.07))+
  xlab("")+
  guides(fill=guide_legend(order = 2,title =""),
         lty=guide_legend(order=1, title = ""))+
  ggsave(filename("Kr_PM_hivol_all_bj"), height = 10, width = 35, units = "cm", dpi = 300)

##ss====
FT_envi_trend_ss=subset(FT_envi_trend,FT_envi_trend$Group=="Seosan")
FT_envi_trend_ss_sel=FT_envi_trend_ss


#FT_envi_trend_ss_sel$conc=ifelse(is.na(FT_envi_trend_ss_sel$conc),0,FT_envi_trend_ss_sel$conc)

FT_envi_trend2_ss=subset(FT_envi_trend2,FT_envi_trend2$Group=="Seosan")

FT_envi_trend2_ss$val2=ifelse(FT_envi_trend2_ss$variable=="OC",FT_envi_trend2_ss$value*4,FT_envi_trend2_ss$value)
FT_envi_trend_ss$conc2=FT_envi_trend_ss$conc*4

FT_envi_trend_ss_rect=subset(FT_envi_trend_ss,is.na(FT_envi_trend_ss$conc))

FT_envi_trend3_ss=subset(FT_envi_trend2_ss,!is.na(FT_envi_trend2_ss$value))

ggplot(data=FT_envi_trend2_ss)+
  geom_area(data = FT_envi_trend_ss ,aes(x=date2, y=conc2,fill=type),position = position_stack(reverse = T),na.rm = T)+
  geom_rect(data=FT_envi_trend_ss_rect,aes(xmin=date2-43200,xmax=date2+43200, ymax=Inf, ymin=-Inf), fill="white")+
  geom_line(data = FT_envi_trend2_ss ,aes(x=date2, y=val2, lty=variable),size=0.7, na.rm = F)+
  geom_line(data = FT_envi_trend3_ss ,aes(x=date2, y=val2, lty=variable),size=0.7, na.rm = F)+
  geom_point(data = FT_envi_trend2_ss ,aes(x=date2, y=val2*1.1),size=NA, na.rm = F,col=NA)+
  scale_y_continuous(name =expression(bold(PM[2.5]~"("*"\u03bcg/"*m^"3"*")")),expand = c(0.02,0.02),
                     sec.axis = sec_axis(~./4, name=expression(bold("OC"~"("*"\u03bcg/"*m^"3"*")")))
  )+
  scale_linetype_manual(values=c(1,2),labels=c(expression(bold(PM[2.5])),expression(bold(OC))))+
  scale_fill_manual(values = c("grey70","#FD7F20"))+
  facet_wrap(.~Group, ncol = 1, scales = "free")+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  scale_x_datetime('',limits = lims,
                   date_breaks = '2 days',
                   date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.text = element_text(colour = "black", size = 16, face = "bold",margin = unit(c(0.3,0.2,0.2,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.text.x = element_text(size = 16,angle = 0,colour = "black", face = "bold",family = "Arial", vjust = 0.5, hjust = 0.5),
        axis.ticks.length = unit(0.15,"cm"),
        axis.ticks = element_line(size = 1.2, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial"),
        axis.title.y.left = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial"),
        axis.title.y.right = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.0,0.1,0.4),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.position = c(0.79,1.07))+
  xlab("")+
  guides(fill=guide_legend(order = 2,title =""),
         lty=guide_legend(order=1, title = ""))+
  ggsave(filename("Kr_PM_hivol_all_ss_nofill"), height = 10, width = 35, units = "cm", dpi = 300)


##sul====
FT_envi_trend_sul=subset(FT_envi_trend,FT_envi_trend$Group=="Seoul")
FT_envi_trend2_sul=subset(FT_envi_trend2,FT_envi_trend2$Group=="Seoul")

FT_envi_trend2_sul$val2=ifelse(FT_envi_trend2_sul$variable=="OC",FT_envi_trend2_sul$value*4,FT_envi_trend2_sul$value)
FT_envi_trend_sul$conc2=FT_envi_trend_sul$conc*4

FT_envi_trend_sul_rect=subset(FT_envi_trend_sul,is.na(FT_envi_trend_sul$conc))

ggplot(data=FT_envi_trend2_sul)+
  geom_area(data = FT_envi_trend_sul ,aes(x=date2, y=conc2,fill=type),position = position_stack(reverse = T),na.rm = T)+
  #geom_rect(data=FT_envi_trend_sul_rect,aes(xmin=date2-43200,xmax=date2+43200, ymax=Inf, ymin=-Inf), fill="white")+
  geom_line(data = FT_envi_trend2_sul ,aes(x=date2, y=val2, lty=variable),size=0.7, na.rm = F)+
  geom_point(data = FT_envi_trend2_sul ,aes(x=date2, y=val2*1.1),size=NA, na.rm = F,col=NA)+
  scale_y_continuous(name =expression(bold(PM[2.5]~"("*"\u03bcg/"*m^"3"*")")),expand = c(0.02,0.02),
                     sec.axis = sec_axis(~./4, name=expression(bold("OC"~"("*"\u03bcg/"*m^"3"*")")))
  )+
  scale_linetype_manual(values=c(1,2),labels=c(expression(bold(PM[2.5])),expression(bold(OC))))+
  scale_fill_manual(values = c("grey70","#FD7F20"))+
  facet_wrap(.~Group, ncol = 1, scales = "free")+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  scale_x_datetime('',limits = lims,
                   date_breaks = '2 days',
                   date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.text = element_text(colour = "black", size = 16, face = "bold",margin = unit(c(0.3,0.2,0.2,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.text.x = element_text(size = 16,angle = 0,colour = "black", face = "bold",family = "Arial", vjust = 0.5, hjust = 0.5),
        axis.ticks.length = unit(0.15,"cm"),
        axis.ticks = element_line(size = 1.2, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial"),
        axis.title.y.left = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial"),
        axis.title.y.right = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.0,0.1,0.4),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.position = c(0.79,1.07))+
  xlab("")+
  guides(fill=guide_legend(order = 2,title =""),
         lty=guide_legend(order=1, title = ""))+
  ggsave(filename("Kr_PM_hivol_all_sul"), height = 10, width = 35, units = "cm", dpi = 300)


##nt====
FT_envi_trend_nt=subset(FT_envi_trend,FT_envi_trend$Group=="Noto")
FT_envi_trend2_nt=subset(FT_envi_trend2,FT_envi_trend2$Group=="Noto")

FT_envi_trend2_nt$val2=ifelse(FT_envi_trend2_nt$variable=="OC",FT_envi_trend2_nt$value*4,FT_envi_trend2_nt$value)
FT_envi_trend_nt$conc2=FT_envi_trend_nt$conc*4

FT_envi_trend_nt_rect=subset(FT_envi_trend_nt,is.na(FT_envi_trend_nt$conc))

FT_envi_trend3_nt=subset(FT_envi_trend2_nt,!is.na(FT_envi_trend2_nt$value))


ggplot(data=FT_envi_trend2_nt)+
  geom_area(data = FT_envi_trend_nt ,aes(x=date2, y=conc2,fill=type),position = position_stack(reverse = T),na.rm = T)+
  #geom_rect(data=FT_envi_trend_nt_rect,aes(xmin=date2-43200,xmax=date2+43200, ymax=Inf, ymin=-Inf), fill="white")+
  geom_line(data = FT_envi_trend2_nt ,aes(x=date2, y=val2, lty=variable),size=0.7, na.rm = F)+
  geom_line(data = FT_envi_trend3_nt ,aes(x=date2, y=val2, lty=variable),size=0.7, na.rm = F)+
  geom_point(data = FT_envi_trend2_nt ,aes(x=date2, y=val2*1.1),size=NA, na.rm = F,col=NA)+
  scale_y_continuous(name =expression(bold(PM[2.5]~"("*"\u03bcg/"*m^"3"*")")),expand = c(0.02,0.02),
                     sec.axis = sec_axis(~./4, name=expression(bold("OC"~"("*"\u03bcg/"*m^"3"*")")))
  )+
  scale_linetype_manual(values=c(1,2),labels=c(expression(bold(PM[2.5])),expression(bold(OC))))+
  scale_fill_manual(values = c("grey70","#FD7F20"))+
  facet_wrap(.~Group, ncol = 1, scales = "free")+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  scale_x_datetime('',limits = lims,
                   date_breaks = '2 days',
                   date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.text = element_text(colour = "black", size = 16, face = "bold",margin = unit(c(0.3,0.2,0.2,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.text.x = element_text(size = 16,angle = 0,colour = "black", face = "bold",family = "Arial", vjust = 0.5, hjust = 0.5),
        axis.ticks.length = unit(0.15,"cm"),
        axis.ticks = element_line(size = 1.2, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial"),
        axis.title.y.left = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial"),
        axis.title.y.right = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.0,0.1,0.4),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.position = c(0.79,1.07))+
  xlab("")+
  guides(fill=guide_legend(order = 2,title =""),
         lty=guide_legend(order=1, title = ""))+
  ggsave(filename("Kr_PM_hivol_all_nt"), height = 10, width = 35, units = "cm", dpi = 300)


###molecular_richness=====
FT_envi
ft_merge$Freq=1

ft_merge$Sample=paste(ft_merge$Group,ft_merge$No,sep = "_")
ft_merge

ft_merge_sel=ft_merge %>% left_join(fm_obs)
ft_merge_sel=ft_merge_sel %>% filter(cnt>1)

mol_rich=as.data.table(aggregate(ft_merge_sel$Freq,by=list(Sample=ft_merge_sel$Sample),sum))
mol_rich

mol_rich=as.data.table(aggregate(ft_merge_sel$Bromo.Inty*10000,by=list(Sample=ft_merge_sel$Sample),sum))
mol_rich

mol_rich_n=mol_rich %>% inner_join(FT_envi_sel[,c("Sample","Group","No","Date","PM2.5","OC","POC","SOC")])

mol_rich_n


mol_rich_n=mol_rich_n %>% 
  mutate(date2 = as.POSIXct(mol_rich_n$Date, format = '%Y-%m-%d %H:%M'))


mol_rich_n_m=melt(mol_rich_n, id.vars = c("Sample","Group","No","Date","date2","x"))
#fwrite(mol_rich_n_m,file = "cal_variation.csv")

mol_rich_nna=subset(mol_rich_n_m,!is.na(mol_rich_n_m$value))

mol_rich_min=as.data.table(aggregate(mol_rich_nna$value,
                                     by=list(Group=mol_rich_nna$Group,variable=mol_rich_nna$variable),min)) %>% 
  `colnames<-`(c("Group","variable","value"))

mol_rich_min2=mol_rich_min %>% left_join(mol_rich_nna[,c("Group","variable","value","Sample","x")])

mol_rich_n_m_sel=mol_rich_n_m
mol_rich_n_m_sel$Group=factor(mol_rich_n_m_sel$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))
mol_rich_n_m_sel$Group2=factor(mol_rich_n_m_sel$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                               labels = c(expression(bold("Ulaanbaatar")),
                                          expression(bold("Beijing")),
                                          expression(bold("Seosan")),
                                          expression(bold("Seoul")),
                                          expression(bold("Noto"))
                               ))

mol_rich_n_m_sel
mol_rich_n_m_sel$variable2=factor(mol_rich_n_m_sel$variable,levels = c("PM2.5","OC","POC","SOC"),
                              labels=c(expression(bold("PM"[2.5]~"("*"\u03bcg/"*m^"3"*")")),
                                       expression(bold("OC"~"("*"\u03bcg/"*m^"3"*")")),
                                       expression(bold("POC"~"("*"\u03bcg/"*m^"3"*")")),
                                       expression(bold("SOC"~"("*"\u03bcg/"*m^"3"*")"))
                                       ))

ggplot()+
  geom_line(data = mol_rich_n_m_sel ,aes(x=value, y=x, col=variable),size=0.7, na.rm = F)+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "loess",formula=y~log(x))+
  #geom_point(data = FT_envi_trend2 ,aes(x=date2, y=value*1.2),size=NA, na.rm = F,col=NA)+
  scale_y_continuous(name = "Sum of molecular intensities")+
  #scale_y_continuous(name = "Sum of molecular intensities")+
  #scale_y_continuous(name = "#AOM",expand = c(0.02,0.02),
  #                   sec.axis = sec_axis(~./fac, name="var2")
  #)+
  facet_wrap(Group2~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed)+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  #scale_x_datetime('',limits = lims,
  #                 date_breaks = '2 days',
  #                 date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 20, face=2 ,colour = "black",margin = unit(c(0.1,0.8,0.1,0.3),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = "NULL")+
  xlab("")+
  ggsave(filename("molvsOC_rich_all_0"), height = 50, width = 50, units = "cm", dpi = 300)


mol_rich_n_ul=subset(mol_rich_n_m_sel,mol_rich_n_m_sel$Group=="Ulaanbaatar")
mol_rich_n_ul$val2=mol_rich_n_ul$x/10000000

ggplot()+
  geom_smooth(data = mol_rich_n_ul ,aes(x=value, y=val2),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_rich_n_ul ,aes(x=value, y=val2, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  #scale_y_continuous(name = "Molecular richness")+
  scale_y_continuous(name = expression(bold("Sum of molecular intensities (x"~"10"^"7"~")")))+
  #scale_y_continuous(name = "#AOM",expand = c(0.02,0.02),
  #                   sec.axis = sec_axis(~./fac, name="var2")
  #)+
  facet_wrap(.~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed, strip.position = "bottom")+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  #scale_x_datetime('',limits = lims,
  #                 date_breaks = '2 days',
  #                 date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.placement = "outside",
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.3,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 16, face=2 ,colour = "black",margin = unit(c(0.1,0.8,0.1,0.3),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = "NULL")+
  xlab("")+
  ggtitle("Ulaanbaatar")+
  ggsave(filename("molvsOC_rich_ul"), height = 11, width = 50, units = "cm", dpi = 300)


mol_rich_n_bj=subset(mol_rich_n_m_sel,mol_rich_n_m_sel$Group=="Beijing")
ggplot()+
  geom_smooth(data = mol_rich_n_bj ,aes(x=value, y=x),col="black", lty=2,method = "loess",formbja=y~log(x))+
  geom_line(data = mol_rich_n_bj ,aes(x=value, y=x, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formbja=y~log(x))+
  #geom_point(data = FT_envi_trend2 ,aes(x=date2, y=value*1.2),size=NA, na.rm = F,col=NA)+
  #scale_y_continuous(name = "Molecular richness")+
  scale_y_continuous(name = "Sum of molecular intensities")+
  #scale_y_continuous(name = "#AOM",expand = c(0.02,0.02),
  #                   sec.axis = sec_axis(~./fac, name="var2")
  #)+
  facet_wrap(.~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed, strip.position = "bottom")+
  #scale_x_datetime('',limits = lims,
  #                 date_breaks = '2 days',
  #                 date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.placement = "outside",
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.3,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 20, face=2 ,colour = "black",margin = unit(c(0.1,0.8,0.1,0.3),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = "NULL")+
  xlab("")+
  ggtitle("Beijing")+
  ggsave(filename("molvsOC_rich_bj"), height = 11, width = 50, units = "cm", dpi = 300)


mol_rich_n_ss=subset(mol_rich_n_m_sel,mol_rich_n_m_sel$Group=="Seosan")
ggplot()+
  geom_smooth(data = mol_rich_n_ss ,aes(x=value, y=x),col="black", lty=2,method = "loess",formssa=y~log(x))+
  geom_line(data = mol_rich_n_ss ,aes(x=value, y=x, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formssa=y~log(x))+
  #geom_point(data = FT_envi_trend2 ,aes(x=date2, y=value*1.2),size=NA, na.rm = F,col=NA)+
  #scale_y_continuous(name = "Molecular richness")+
  scale_y_continuous(name = "Sum of molecular intensities")+
  #scale_y_continuous(name = "#AOM",expand = c(0.02,0.02),
  #                   sec.axis = sec_axis(~./fac, name="var2")
  #)+
  facet_wrap(.~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed, strip.position = "bottom")+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  #scale_x_datetime('',limits = lims,
  #                 date_breaks = '2 days',
  #                 date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.placement = "outside",
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.3,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 20, face=2 ,colour = "black",margin = unit(c(0.1,0.8,0.1,0.3),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = "NULL")+
  xlab("")+
  ggtitle("Seosan")+
  ggsave(filename("molvsOC_rich_ss"), height = 11, width = 50, units = "cm", dpi = 300)


mol_rich_n_sul=subset(mol_rich_n_m_sel,mol_rich_n_m_sel$Group=="Seoul")
ggplot()+
  geom_smooth(data = mol_rich_n_sul ,aes(x=value, y=x),col="black", lty=2,method = "loess",formsula=y~log(x))+
  geom_line(data = mol_rich_n_sul ,aes(x=value, y=x, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formsula=y~log(x))+
  #geom_point(data = FT_envi_trend2 ,aes(x=date2, y=value*1.2),size=NA, na.rm = F,col=NA)+
  #scale_y_continuous(name = "Molecular richness")+
  scale_y_continuous(name = "Sum of molecular intensities")+
  #scale_y_continuous(name = "#AOM",expand = c(0.02,0.02),
  #                   sec.axis = sec_axis(~./fac, name="var2")
  #)+
  facet_wrap(.~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed, strip.position = "bottom")+
  #scale_x_datetime('',limits = lims,
  #                 date_breaks = '2 days',
  #                 date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.placement = "outside",
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.3,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 20, face=2 ,colour = "black",margin = unit(c(0.1,0.8,0.1,0.3),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = "NULL")+
  xlab("")+
  ggtitle("Seoul")+
  ggsave(filename("molvsOC_rich_sul"), height = 11, width = 50, units = "cm", dpi = 300)

mol_rich_n_nt=subset(mol_rich_n_m_sel,mol_rich_n_m_sel$Group=="Noto")
ggplot()+
  geom_smooth(data = mol_rich_n_nt ,aes(x=value, y=x),col="black", lty=2,method = "loess",formnta=y~log(x))+
  geom_line(data = mol_rich_n_nt ,aes(x=value, y=x, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formnta=y~log(x))+
  #geom_point(data = FT_envi_trend2 ,aes(x=date2, y=value*1.2),size=NA, na.rm = F,col=NA)+
  #scale_y_continuous(name = "Molecular richness")+
  scale_y_continuous(name = "Sum of molecular intensities")+
  #scale_y_continuous(name = "#AOM",expand = c(0.02,0.02),
  #                   sec.axis = sec_axis(~./fac, name="var2")
  #)+
  facet_wrap(.~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed, strip.position = "bottom")+
  #scale_x_datetime('',limits = lims,
  #                 date_breaks = '2 days',
  #                 date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.placement = "outside",
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.3,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 20, face=2 ,colour = "black",margin = unit(c(0.1,0.8,0.1,0.3),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = "NULL")+
  xlab("")+
  ggtitle("Noto")+
  ggsave(filename("molvsOC_rich_nt"), height = 11, width = 50, units = "cm", dpi = 300)

###intensity weight=====
mol_rich_n_ul=subset(mol_rich_n_m_sel,mol_rich_n_m_sel$Group=="Ulaanbaatar")
mol_rich_n_ul$val2=mol_rich_n_ul$x/10000000

ggplot()+
  geom_smooth(data = mol_rich_n_ul ,aes(x=value, y=val2),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_rich_n_ul ,aes(x=value, y=val2, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  #scale_y_continuous(name = "Molecular richness")+
  scale_y_continuous(name = expression(bold("Sum of molecular intensities (x"~"10"^"7"~")")))+
  #scale_y_continuous(name = "#AOM",expand = c(0.02,0.02),
  #                   sec.axis = sec_axis(~./fac, name="var2")
  #)+
  facet_wrap(.~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed, strip.position = "bottom")+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  #scale_x_datetime('',limits = lims,
  #                 date_breaks = '2 days',
  #                 date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.placement = "outside",
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.3,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 16, face=2 ,colour = "black",margin = unit(c(0.1,0.8,0.1,0.3),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = "NULL")+
  xlab("")+
  ggtitle("Ulaanbaatar")+
  ggsave(filename("molvsOC_rich_ul"), height = 11, width = 50, units = "cm", dpi = 300)


mol_rich_n_bj=subset(mol_rich_n_m_sel,mol_rich_n_m_sel$Group=="Beijing")
mol_rich_n_bj$val2=mol_rich_n_bj$x/10000000

ggplot()+
  geom_smooth(data = mol_rich_n_bj ,aes(x=value, y=val2),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_rich_n_bj ,aes(x=value, y=val2, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  #scale_y_continuous(name = "Molecular richness")+
  scale_y_continuous(name = expression(bold("Sum of molecular intensities (x"~"10"^"7"~")")))+
  #scale_y_continuous(name = "#AOM",expand = c(0.02,0.02),
  #                   sec.axis = sec_axis(~./fac, name="var2")
  #)+
  facet_wrap(.~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed, strip.position = "bottom")+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  #scale_x_datetime('',limits = lims,
  #                 date_breaks = '2 days',
  #                 date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.placement = "outside",
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.3,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 16, face=2 ,colour = "black",margin = unit(c(0.1,0.8,0.1,0.3),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = "NULL")+
  xlab("")+
  ggtitle("Beijing")+
  ggsave(filename("molvsOC_rich_bj"), height = 11, width = 50, units = "cm", dpi = 300)

mol_rich_n_ss=subset(mol_rich_n_m_sel,mol_rich_n_m_sel$Group=="Seosan")
mol_rich_n_ss$val2=mol_rich_n_ss$x/10000000

ggplot()+
  geom_smooth(data = mol_rich_n_ss ,aes(x=value, y=val2),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_rich_n_ss ,aes(x=value, y=val2, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  #scale_y_continuous(name = "Molecular richness")+
  scale_y_continuous(name = expression(bold("Sum of molecular intensities (x"~"10"^"7"~")")))+
  #scale_y_continuous(name = "#AOM",expand = c(0.02,0.02),
  #                   sec.axis = sec_axis(~./fac, name="var2")
  #)+
  facet_wrap(.~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed, strip.position = "bottom")+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  #scale_x_datetime('',limits = lims,
  #                 date_breaks = '2 days',
  #                 date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.placement = "outside",
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.3,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 16, face=2 ,colour = "black",margin = unit(c(0.1,0.8,0.1,0.3),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = "NULL")+
  xlab("")+
  ggtitle("Seosan")+
  ggsave(filename("molvsOC_rich_ss"), height = 11, width = 50, units = "cm", dpi = 300)

mol_rich_n_sul=subset(mol_rich_n_m_sel,mol_rich_n_m_sel$Group=="Seoul")
mol_rich_n_sul$val2=mol_rich_n_sul$x/10000000

ggplot()+
  geom_smooth(data = mol_rich_n_sul ,aes(x=value, y=val2),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_rich_n_sul ,aes(x=value, y=val2, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  #scale_y_continuous(name = "Molecular richness")+
  scale_y_continuous(name = expression(bold("Sum of molecular intensities (x"~"10"^"7"~")")))+
  #scale_y_continuous(name = "#AOM",expand = c(0.02,0.02),
  #                   sec.axis = sec_axis(~./fac, name="var2")
  #)+
  facet_wrap(.~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed, strip.position = "bottom")+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  #scale_x_datetime('',limits = lims,
  #                 date_breaks = '2 days',
  #                 date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.placement = "outside",
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.3,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 16, face=2 ,colour = "black",margin = unit(c(0.1,0.8,0.1,0.3),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = "NULL")+
  xlab("")+
  ggtitle("Seoul")+
  ggsave(filename("molvsOC_rich_sul"), height = 11, width = 50, units = "cm", dpi = 300)

mol_rich_n_nt=subset(mol_rich_n_m_sel,mol_rich_n_m_sel$Group=="Noto")
mol_rich_n_nt$val2=mol_rich_n_nt$x/10000000

ggplot()+
  geom_smooth(data = mol_rich_n_nt ,aes(x=value, y=val2),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_rich_n_nt ,aes(x=value, y=val2, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  #scale_y_continuous(name = "Molecular richness")+
  scale_y_continuous(name = expression(bold("Sum of molecular intensities (x"~"10"^"7"~")")))+
  #scale_y_continuous(name = "#AOM",expand = c(0.02,0.02),
  #                   sec.axis = sec_axis(~./fac, name="var2")
  #)+
  facet_wrap(.~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed, strip.position = "bottom")+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  #scale_x_datetime('',limits = lims,
  #                 date_breaks = '2 days',
  #                 date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.placement = "outside",
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.3,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 16, face=2 ,colour = "black",margin = unit(c(0.1,0.8,0.1,0.3),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = "NULL")+
  xlab("")+
  ggtitle("Noto")+
  ggsave(filename("molvsOC_rich_nt"), height = 11, width = 50, units = "cm", dpi = 300)


##cal_variation=====
mol_rich_n_m_sel
mol_rich_min



mol_rich_n_m_sel=mol_rich_n_m_sel[,c("Group","No","Date","date2","x","variable","value")] %>% 
  `colnames<-`(c("Group","No","Date","date2","Richness","variable","conc"))

mol_var_min=mol_rich_n_m_sel %>% inner_join(mol_rich_min2)

mol_var_min$difx=mol_var_min$conc-mol_var_min$value
mol_var_min$dify=mol_var_min$Richness-mol_var_min$x

mol_var_min

mol_var_min$Group2=factor(mol_var_min$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                               labels = c(expression(bold("Ulaanbaatar")),
                                          expression(bold("Beijing")),
                                          expression(bold("Seosan")),
                                          expression(bold("Seoul")),
                                          expression(bold("Noto"))
                               ))

mol_var_min
mol_var_min$variable2=factor(mol_var_min$variable,levels = c("PM2.5","OC","POC","SOC"),
                                  labels=c(expression(bold("PM"[2.5]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("OC"~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("POC"~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("SOC"~"("*"\u03bcg/"*m^"3"*")"))
                                  ))

ggplot()+
  geom_line(data = mol_var_min ,aes(x=difx, y=dify, col=variable),size=0.7, na.rm = F)+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "loess",formula=y~log(x))+
  #geom_point(data = FT_envi_trend2 ,aes(x=date2, y=value*1.2),size=NA, na.rm = F,col=NA)+
  #scale_y_continuous(name = "Sum of molecular intensities")+
  #scale_y_continuous(name = "Sum of molecular intensities")+
  #scale_y_continuous(name = "#AOM",expand = c(0.02,0.02),
  #                   sec.axis = sec_axis(~./fac, name="var2")
  #)+
  facet_wrap(Group2~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed)+
  #scale_color_manual(values = c("black","#754F44","#0072B2","#009E73","#631879B2"))+
  #scale_fill_manual(values = c("orangered1","grey25"))+
  #scale_x_datetime('',limits = lims,
  #                 date_breaks = '2 days',
  #                 date_labels = "%m/%d", expand = c(0.008,0.008))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.2,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 20, face=2 ,colour = "black",margin = unit(c(0.1,0.8,0.1,0.3),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = "NULL")+
  xlab("")+
  ggsave(filename("rich_variation"), height = 50, width = 50, units = "cm", dpi = 300)


fwrite(mol_var_min,file = "Datafile/molecular_variation.csv")

