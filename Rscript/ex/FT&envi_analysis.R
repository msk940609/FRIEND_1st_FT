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
library(lawstat)
library(agricolae)
library(PMCMR)
library(PMCMRplus)
source("Rscript/func_filename.R")

library(ggrepel)

#0.Load datat====
FT_merge=fread("Datafile/FRIENDs_1st_FT.csv")
FT_merge$Sample=paste(FT_merge$Group,FT_merge$No,sep = "_")

FT_merge$C=ifelse(FT_merge$C.>0, "C","")
FT_merge$H=ifelse(FT_merge$H.>0, "H","")
FT_merge$O=ifelse(FT_merge$O.>0, "O","")
FT_merge$N=ifelse(FT_merge$N.>0, "N","")
FT_merge$S=ifelse(FT_merge$S.>0, "S","")

FT_merge=FT_merge %>% unite("Comp",c("C","H","O","N","S"),sep = "")
FT_merge$Comp=ifelse(FT_merge$O.==0, "Remainders",FT_merge$Comp)

FT_merge=subset(FT_merge,FT_merge$Comp!="Remainders") %>% droplevels()
FT_merge

FT_SOA=fread("Datafile/FRIENDs_1st_FT_SOA_S&B.csv")
FT_SOA$Sample=paste(FT_SOA$Group,FT_SOA$No,sep = "_")


ft_inty_all=as.data.table(aggregate(FT_merge$Bromo.Inty, by=list(`Sample`=FT_merge$Sample), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Tot"))
ft_inty=as.data.table(aggregate(FT_merge$Bromo.Inty, by=list(`Sample`=FT_merge$Sample, Comp=FT_merge$Comp), FUN=sum))
ft_soa=as.data.table(aggregate(FT_SOA$Bromo.Inty, by=list(`Sample`=FT_SOA$Sample, Comp=FT_SOA$SOA), FUN=sum))
ft_inty
ft_soa=subset(ft_soa,ft_soa$Comp=="SOA")

ft_inty=rbind(ft_inty,ft_soa)
  
ft_inty=ft_inty %>% inner_join(ft_inty_all,by = "Sample")
ft_inty$rel=round(ft_inty$x/ft_inty$Tot*100,1)

ft_inty=ft_inty %>% separate(Sample, c("Group","No"),sep = "_")
ft_inty

ft_inty$Comp=factor(ft_inty$Comp,levels = c("CHO","CHON","CHOS","CHONS","SOA"))

ft_inty$Group=factor(ft_inty$Group,levels = c("SUL","SS","B","NT","UL"),
                     labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
ft_inty

ft_inty_m=dcast(ft_inty,Group+No~Comp, sum, value.var = "rel")

ft_inty_m$Sample=paste(ft_inty_m$Group,ft_inty_m$No,sep = "_")

envi_all=fread("Datafile/FRIENDs_envi_S&B.csv")
envi_all

envi_all=envi_all %>% separate(Sample, c("Group","No"), sep = "_")

envi_all$Sample=paste(envi_all$Group,envi_all$No, sep = "_")

envi_all_sel=envi_all[, c("Sample","Group","No","OC","PM2.5","SO4","NO3","NH4","O3","CO","SO2","NO")]

envi_all_sel$`NO3/OC`=envi_all_sel$NO3/envi_all_sel$OC
envi_all_sel$`SO4/OC`=envi_all_sel$SO4/envi_all_sel$OC
envi_all_sel$`NO3/SO4`=envi_all_sel$NO3/envi_all_sel$SO4

ft_inty_env=ft_inty_m %>% inner_join(envi_all_sel)
wsoa=fread("Datafile/WOSCbb_S&B.csv")

ft_inty_env=ft_inty_env %>% inner_join(wsoa[,c("Sample","WSOCbb","WSOCnbb")])


ft_inty_env_m=melt(ft_inty_env, id.vars = c("Sample","Group","No","OC","CHO","CHON","CHOS","CHONS","SOA"), variable.name = "Envi", value.name = "Conc") %>% 
  melt(id.vars=c("Sample","Group","No","OC","Envi","Conc"),variable.name="Comp")
  
ft_inty_env_m$No=as.numeric(ft_inty_env_m$No)
ft_inty_env_m$mark="no"
ft_inty_env_m$mark=ifelse(ft_inty_env_m$No==5,"12/19",ft_inty_env_m$mark)
ft_inty_env_m$mark=ifelse(ft_inty_env_m$No==6,"12/20",ft_inty_env_m$mark)
ft_inty_env_m$mark=ifelse(ft_inty_env_m$No==7,"12/21",ft_inty_env_m$mark)
ft_inty_env_m$mark=ifelse(ft_inty_env_m$No==8,"12/22",ft_inty_env_m$mark)
ft_inty_env_m$mark=ifelse(ft_inty_env_m$No==9,"12/23",ft_inty_env_m$mark)

ft_inty_env_m
table(ft_inty_env_m$mark)

ft_inty_env_m_s=subset(ft_inty_env_m,ft_inty_env_m$Group=="Seoul")
ft_inty_env_m_s1=ft_inty_env_m_s %>% filter(Envi%in%c("PM2.5","NH4","NO3","SO4")) %>% droplevels()


ggplot()+
  geom_smooth(data=ft_inty_env_m_s1,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=ft_inty_env_m_s1,aes(x=value,y=Conc),col="royalblue2", size=3)+
  geom_blank(data=val_max_s1, aes(x=20,y=Conc*1.4))+
  geom_point(data=subset(ft_inty_env_m_s1,ft_inty_env_m_s1$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_inty_env_m_s1,ft_inty_env_m_s1$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envi~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of chemical composition (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsEnvi_test"),height = 30, width = 38, units = "cm", dpi = 300)


ft_inty_env_m_s1=ft_inty_env_m_s %>% filter(Envi%in%c("PM2.5","NH4","NO3","SO4")) %>% droplevels()
ft_inty_env_m_s1$Envilab=ft_inty_env_m_s1$Envi
ft_inty_env_m_s1$Envilab=factor(ft_inty_env_m_s1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                                labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))

ft_inty_env_m_s1$Complab=ft_inty_env_m_s1$Comp
ft_inty_env_m_s1$Complab=factor(ft_inty_env_m_s1$Complab, levels = c("CHO","CHON","CHOS","CHONS","SOA"),
                                labels = c(expression(bold("CHO")),expression(bold("CHON")),expression(bold("CHOS")),
                                           expression(bold("CHONS")),expression(bold("SOA"))))

val_max_s1=as.data.table(aggregate(ft_inty_env_m_s1$Conc, 
                                   by=list(`Comp`=ft_inty_env_m_s1$Comp,Envi=ft_inty_env_m_s1$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))
val_max_s1$Envilab=val_max_s1$Envi

val_max_s1$Envilab=factor(val_max_s1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                                labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))
val_max_s1

ft_inty_env_m_s1

ggplot()+
  geom_smooth(data=ft_inty_env_m_s1,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=ft_inty_env_m_s1,aes(x=value,y=Conc),col="royalblue2", size=3)+
  geom_blank(data=val_max_s1, aes(x=20,y=Conc*1.4))+
  geom_point(data=subset(ft_inty_env_m_s1,ft_inty_env_m_s1$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_inty_env_m_s1,ft_inty_env_m_s1$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of chemical composition (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsEnvi_sul1"),height = 30, width = 37.5, units = "cm", dpi = 300)


ft_inty_env_m_s2=ft_inty_env_m_s %>% filter(Envi%in%c("O3","CO","SO2","NO")) %>% droplevels()
ft_inty_env_m_s2$Envilab=ft_inty_env_m_s2$Envi

ft_inty_env_m_s2$Envilab=factor(ft_inty_env_m_s2$Envilab, levels = c("O3","CO","NO","SO2"),
                                labels = c(expression(bold("O"["3"]~"(ppb)")),
                                             expression(bold("CO (ppb)")),
                                             expression(bold("NO (ppb)")),
                                             expression(bold("SO"["2"]~"(ppb)"))))

ft_inty_env_m_s2$Complab=ft_inty_env_m_s2$Comp
ft_inty_env_m_s2$Complab=factor(ft_inty_env_m_s2$Complab, levels = c("CHO","CHON","CHOS","CHONS","SOA"),
                                labels = c(expression(bold("CHO")),expression(bold("CHON")),expression(bold("CHOS")),
                                           expression(bold("CHONS")),expression(bold("SOA"))))


val_max_s2=as.data.table(aggregate(ft_inty_env_m_s2$Conc, 
                                   by=list(`Comp`=ft_inty_env_m_s2$Comp,Envi=ft_inty_env_m_s2$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))

val_max_s2$Envilab=val_max_s2$Envi

val_max_s2$Envilab=factor(val_max_s2$Envilab, levels = c("O3","CO","NO","SO2"),
                                labels = c(expression(bold("O"["3"]~"(ppb)")),
                                           expression(bold("CO (ppb)")),
                                           expression(bold("NO (ppb)")),
                                           expression(bold("SO"["2"]~"(ppb)"))))

val_max_s2$Conc=ifelse(val_max_s2$Envi=="SO2",2.0,val_max_s2$Conc)

ggplot(ft_inty_env_m_s2, aes(x=value,y=Conc))+
  geom_smooth(formula = y~x, method = "lm", se = T, col="black", lty=2)+
  geom_point(col="royalblue2", size=3)+
  geom_blank(data=val_max_s2, aes(x=20,y=Conc*1.35))+
  geom_point(data=subset(ft_inty_env_m_s2,ft_inty_env_m_s2$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_inty_env_m_s2,ft_inty_env_m_s2$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of chemical composition (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsEnvi_sul2"),height = 30, width = 37.5, units = "cm", dpi = 300)

ft_inty_env_m_s3=ft_inty_env_m_s %>% filter(Envi%in%c("PM2.5","NO3/OC","SO4/OC","NO3/SO4")) %>% droplevels()

ft_inty_env_m_s3$Envilab=ft_inty_env_m_s3$Envi
ft_inty_env_m_s3$Envilab=factor(ft_inty_env_m_s3$Envilab, levels = c("PM2.5","NO3/OC","SO4/OC","NO3/SO4"),
                                labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NO"["3"]^{"-"}*"/OC")),
                                           expression(bold("SO"["4"]^{"2-"}*"/OC")),
                                           expression(bold("NO"["3"]^{"-"}*"/"*"SO"["4"]^{"2-"}))))
ft_inty_env_m_s3

ft_inty_env_m_s3$Complab=ft_inty_env_m_s3$Comp
ft_inty_env_m_s3$Complab=factor(ft_inty_env_m_s3$Complab, levels = c("CHO","CHON","CHOS","CHONS","SOA"),
                                labels = c(expression(bold("CHO")),expression(bold("CHON")),expression(bold("CHOS")),
                                           expression(bold("CHONS")),expression(bold("SOA"))))

val_max_s3=as.data.table(aggregate(ft_inty_env_m_s3$Conc, 
                                   by=list(`Comp`=ft_inty_env_m_s3$Comp,Envi=ft_inty_env_m_s3$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))
val_max_s3$Envilab=val_max_s3$Envi

val_max_s3$Envilab=factor(val_max_s3$Envilab, levels = c("PM2.5","NO3/OC","SO4/OC","NO3/SO4"),
                                labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NO"["3"]^{"-"}*"/OC")),
                                           expression(bold("SO"["4"]^{"2-"}*"/OC")),
                                           expression(bold("NO"["3"]^{"-"}*"/"*"SO"["4"]^{"2-"}))))

val_max_s3

ggplot()+
  geom_smooth(data=ft_inty_env_m_s3,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=ft_inty_env_m_s3,aes(x=value,y=Conc),col="royalblue2", size=3)+
  geom_blank(data=val_max_s3, aes(x=20,y=Conc*1.4))+
  geom_point(data=subset(ft_inty_env_m_s3,ft_inty_env_m_s3$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_inty_env_m_s3,ft_inty_env_m_s3$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_grid(Envi~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of chemical composition (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsEnvi_sul3"),height = 30, width = 37.5, units = "cm", dpi = 300)


ft_inty_env_m_b=subset(ft_inty_env_m,ft_inty_env_m$Group=="Beijing")
ft_inty_env_m_b1=ft_inty_env_m_b %>% filter(Envi%in%c("PM2.5","NH4","NO3","SO4")) %>% droplevels()

ft_inty_env_m_b1$Envilab=ft_inty_env_m_b1$Envi
ft_inty_env_m_b1$Envilab=factor(ft_inty_env_m_b1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                                labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))

ft_inty_env_m_b1$Complab=ft_inty_env_m_b1$Comp
ft_inty_env_m_b1$Complab=factor(ft_inty_env_m_b1$Complab, levels = c("CHO","CHON","CHOS","CHONS","SOA"),
                                labels = c(expression(bold("CHO")),expression(bold("CHON")),expression(bold("CHOS")),
                                           expression(bold("CHONS")),expression(bold("SOA"))))

val_max_b1=as.data.table(aggregate(ft_inty_env_m_b1$Conc, 
                                   by=list(`Comp`=ft_inty_env_m_b1$Comp,Envi=ft_inty_env_m_b1$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))
val_max_b1$Envilab=val_max_b1$Envi

val_max_b1$Envilab=factor(val_max_b1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                          labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))
val_max_b1

ggplot(ft_inty_env_m_b1, aes(x=value,y=Conc))+
  geom_smooth(formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(col="orangered3", size=3)+
  geom_blank(data=val_max_b1, aes(x=20,y=Conc*1.32))+
  geom_point(data=subset(ft_inty_env_m_b1,ft_inty_env_m_b1$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_inty_env_m_b1,ft_inty_env_m_b1$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.6,col="black", size=4)+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of chemical composition (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsEnvi_beigjing1"),height = 30, width = 37.5, units = "cm", dpi = 300)


ft_inty_env_m_b2=ft_inty_env_m_b %>% filter(Envi%in%c("O3","CO","SO2","NO")) %>% droplevels()
ft_inty_env_m_b2$Envilab=ft_inty_env_m_b2$Envi

ft_inty_env_m_b2$Envilab=factor(ft_inty_env_m_b2$Envilab, levels = c("O3","CO","NO","SO2"),
                                labels = c(expression(bold("O"["3"]~"(ppb)")),
                                           expression(bold("CO (ppb)")),
                                           expression(bold("NO (ppb)")),
                                           expression(bold("SO"["2"]~"(ppb)"))))

ft_inty_env_m_b2$Complab=ft_inty_env_m_b2$Comp
ft_inty_env_m_b2$Complab=factor(ft_inty_env_m_b2$Complab, levels = c("CHO","CHON","CHOS","CHONS","SOA"),
                                labels = c(expression(bold("CHO")),expression(bold("CHON")),expression(bold("CHOS")),
                                           expression(bold("CHONS")),expression(bold("SOA"))))

val_max_b2=as.data.table(aggregate(ft_inty_env_m_b2$Conc, 
                                   by=list(`Comp`=ft_inty_env_m_b2$Comp,Envi=ft_inty_env_m_b2$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))

val_max_b2$Envilab=val_max_b2$Envi

val_max_b2$Envilab=factor(val_max_b2$Envilab, levels = c("O3","CO","NO","SO2"),
                          labels = c(expression(bold("O"["3"]~"(ppb)")),
                                     expression(bold("CO (ppb)")),
                                     expression(bold("NO (ppb)")),
                                     expression(bold("SO"["2"]~"(ppb)"))))

val_max_b2$Conc=ifelse(val_max_b2$Envi=="CO",1050,val_max_b2$Conc)

ggplot(ft_inty_env_m_b2, aes(x=value,y=Conc))+
  geom_smooth(formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(col="orangered3", size=3)+
  geom_blank(data=val_max_b2, aes(x=20,y=Conc*1.32))+
  geom_point(data=subset(ft_inty_env_m_b2,ft_inty_env_m_b2$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_inty_env_m_b2,ft_inty_env_m_b2$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of chemical composition (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsEnvi_beijing2"),height = 30, width = 37.5, units = "cm", dpi = 300)


ft_inty_env_m_b3=ft_inty_env_m_b %>% filter(Envi%in%c("PM2.5","NO3/OC","SO4/OC","NO3/SO4")) %>% droplevels()

ft_inty_env_m_b3$Envilab=ft_inty_env_m_b3$Envi
ft_inty_env_m_b3$Envilab=factor(ft_inty_env_m_b3$Envilab, levels = c("PM2.5","NO3/OC","SO4/OC","NO3/SO4"),
                                labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NO"["3"]^{"-"}*"/OC")),
                                           expression(bold("SO"["4"]^{"2-"}*"/OC")),
                                           expression(bold("NO"["3"]^{"-"}*"/"*"SO"["4"]^{"2-"}))))
ft_inty_env_m_b3

ft_inty_env_m_b3$Complab=ft_inty_env_m_b3$Comp
ft_inty_env_m_b3$Complab=factor(ft_inty_env_m_b3$Complab, levels = c("CHO","CHON","CHOS","CHONS","SOA"),
                                labels = c(expression(bold("CHO")),expression(bold("CHON")),expression(bold("CHOS")),
                                           expression(bold("CHONS")),expression(bold("SOA"))))

val_max_b3=as.data.table(aggregate(ft_inty_env_m_b3$Conc, 
                                   by=list(`Comp`=ft_inty_env_m_b3$Comp,Envi=ft_inty_env_m_b3$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))
val_max_b3$Envilab=val_max_b3$Envi

val_max_b3$Envilab=factor(val_max_b3$Envilab, levels = c("PM2.5","NO3/OC","SO4/OC","NO3/SO4"),
                          labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NO"["3"]^{"-"}*"/OC")),
                                     expression(bold("SO"["4"]^{"2-"}*"/OC")),
                                     expression(bold("NO"["3"]^{"-"}*"/"*"SO"["4"]^{"2-"}))))

val_max_b3

val_max_b3$Conc=ifelse(val_max_b3$Envi=="NO3/OC",1.5,val_max_b3$Conc)
val_max_b3$Conc=ifelse(val_max_b3$Envi=="SO4/OC",0.55,val_max_b3$Conc)

ft_inty_env_m_b3

ggplot()+
  geom_smooth(data=ft_inty_env_m_b3,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=ft_inty_env_m_b3,aes(x=value,y=Conc),col="orangered3", size=3)+
  geom_blank(data=val_max_b3, aes(x=20,y=Conc*1.4))+
  geom_point(data=subset(ft_inty_env_m_b3,ft_inty_env_m_b3$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_inty_env_m_b3,ft_inty_env_m_b3$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_grid(Envi~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of chemical composition (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsEnvi_b3"),height = 30, width = 37.5, units = "cm", dpi = 300)

ft_inty_env_m_u=subset(ft_inty_env_m,ft_inty_env_m$Group=="Ulaanbaatar")
ft_inty_env_m_u1=ft_inty_env_m_u %>% filter(Envi%in%c("PM2.5","NH4","NO3","SO4")) %>% droplevels()

ft_inty_env_m_u1$Envilab=ft_inty_env_m_u1$Envi
ft_inty_env_m_u1$Envilab=factor(ft_inty_env_m_u1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                                labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))

ft_inty_env_m_u1$Complab=ft_inty_env_m_u1$Comp
ft_inty_env_m_u1$Complab=factor(ft_inty_env_m_u1$Complab, levels = c("CHO","CHON","CHOS","CHONS"),
                                labels = c(expression(bold("CHO")),expression(bold("CHON")),expression(bold("CHOS")),
                                           expression(bold("CHONS"))))

val_max_u1=as.data.table(aggregate(ft_inty_env_m_u1$Conc, 
                                   by=list(`Comp`=ft_inty_env_m_u1$Comp,Envi=ft_inty_env_m_u1$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))
val_max_u1$Envilab=val_max_u1$Envi

val_max_u1$Envilab=factor(val_max_u1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                          labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))
val_max_u1
ft_inty_env_m_u1

ggplot(ft_inty_env_m_u1, aes(x=value,y=Conc))+
  geom_smooth(formula = y~x, method = "lm", se = T, col="black")+
  geom_point(col="seagreen4", size=3)+
  geom_blank(data=val_max_u1, aes(x=ifelse(Comp=="CHON",15,
                                           ifelse(Comp=="CHOS",35,20)),y=Conc*1.32))+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of chemical composition (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsEnvi_ulan1"),height = 30, width = 30, units = "cm", dpi = 300)


ft_inty_env_m_u2=ft_inty_env_m_u %>% filter(Envi%in%c("CO","SO2","NO")) %>% droplevels()
ft_inty_env_m_u2$Envilab=ft_inty_env_m_u2$Envi

ft_inty_env_m_u2$Envilab=factor(ft_inty_env_m_u2$Envilab, levels = c("CO","NO","SO2"),
                                labels = c(expression(bold("CO (ppb)")),
                                           expression(bold("NO (ppb)")),
                                           expression(bold("SO"["2"]~"(ppb)"))))

ft_inty_env_m_u2$Complab=ft_inty_env_m_u2$Comp
ft_inty_env_m_u2$Complab=factor(ft_inty_env_m_u2$Complab, levels = c("CHO","CHON","CHOS","CHONS"),
                                labels = c(expression(bold("CHO")),expression(bold("CHON")),expression(bold("CHOS")),
                                           expression(bold("CHONS"))))




val_max_u2=as.data.table(aggregate(ft_inty_env_m_u2$Conc, 
                                   by=list(`Comp`=ft_inty_env_m_u2$Comp,Envi=ft_inty_env_m_u2$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))

val_max_u2$Envilab=val_max_u2$Envi

val_max_u2$Envilab=factor(val_max_u2$Envilab, levels = c("CO","NO","SO2"),
                          labels = c(expression(bold("CO (ppb)")),
                                     expression(bold("NO (ppb)")),
                                     expression(bold("SO"["2"]~"(ppb)"))))

#val_max_u2$Conc=ifelse(val_max_u2$Envi=="SO2",2.0,val_max_u2$Conc)

ggplot(ft_inty_env_m_u2, aes(x=value,y=Conc))+
  geom_smooth(formula = y~x, method = "lm", se = T, col="black")+
  geom_point(col="seagreen4", size=3)+
  geom_blank(data=val_max_u2, aes(x=ifelse(Comp=="CHON",15,
                                           ifelse(Comp=="CHOS",35,20)),y=Conc*1.32))+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of chemical composition (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsEnvi_ulan2"),height = 23, width = 30, units = "cm", dpi = 300)


##FT_envi correlation======
ft_inty_env_m

grp=unique(ft_inty_env_m$Group)
grp
comp=unique(ft_inty_env_m$Comp)
envi=unique(ft_inty_env_m$Envi)


dt=data.table()
for (i in 1:length(grp)) {
  temp=subset(ft_inty_env_m,ft_inty_env_m$Group==grp[i])
  
    for (k in 1:length(comp)) {
      tmp=subset(temp,temp$Comp==comp[k])
      
      for (l in 1:length(envi)) {
        tmp2=subset(tmp,tmp$Envi==envi[l])
        
        if(envi[l]=="O3"&grp[i]=="Ulaanbaatar") {
          tmp2$Conc=0
          new_row=data.table("Group"=grp[i],"Comp"=comp[k],"R2"=0,"Envi"=envi[l],"Cor"=0, "P"=1)
          
        } else {
          model=lm(Conc~value,tmp2)
          a=summary(model)
          #cor=cor.test(tmp2$value,tmp2$Conc, method = "pearson")
          cor=cor.test(tmp2$value,tmp2$Conc, method = "spearman", exact = F)
          #    cor
          
          new_row=data.table("Group"=grp[i],"Comp"=comp[k],"R2"=round(a$r.squared,8),"Envi"=envi[l],"Cor"=cor$estimate, "P"=round(cor$p.value,3))
          
        }
        
        dt=rbind(dt,new_row)
      }
    }
  }

dt

dt$sig=ifelse(dt$P<0.05,ifelse(dt$P<0.01,ifelse(dt$P<0.001,"***","**"),"*"),"")


fwrite(dt,"Datafile/211122FT&Envicor.csv")

dt$Envilab=dt$Envi
dt$Envilab=factor(dt$Envilab,levels = c("PM2.5","NH4","NO3","SO4","O3","CO","NO","SO2"),
                  labels = c(expression(bold("PM"["2.5"])),
                              expression(bold("NH"["4"]^{"+"})),
                              expression(bold("NO"["3"]^{"-"})),
                              expression(bold("SO"["4"]^{"2-"})),
                              expression(bold("O"["3"])),
                              expression(bold("CO")),
                              expression(bold("NO")),
                              expression(bold("SO"["2"]))
                             )
                  )

dt$Group=factor(dt$Group, levels = c("Ulaanbaatar","Beijing","Seoul"))

p=ggplot(dt, aes(x=Envilab,y=Cor,fill=Group))+
  geom_bar(stat="identity", position = "dodge", alpha=0.85)+
  scale_x_discrete("Environmental variables", labels=parse(text=levels(dt$Envilab)))+
  scale_y_continuous(name = "Spearman correlation coefficient (rho)                         Spearman correlation coefficient (rho)")+
  scale_fill_manual(values = c("seagreen4","orangered3","royalblue2"))+
  facet_wrap(.~Comp, scales = "free")+
  geom_text(aes(Envilab, y=ifelse(Cor>0,Cor*1.05,Cor*1.08), label=sig),size=10, position = position_dodge(width = 0.9))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.6,0,0.2),"cm")),
        strip.text = element_text(size = 20, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 26, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.key.width = unit(1.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.07,0.93),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )

p+ggsave(filename("compvsEnvi"),height = 40, width = 60, units = "cm", dpi = 300)

dt
dt$Envi=factor(dt$Envi,levels = rev(c("PM2.5","NH4","NO3","SO4","O3","CO","NO","SO2")))


matlab.like2(10)

p2=ggplot(dt)+
  geom_tile(data=subset(dt,dt$P<0.05), aes(x=Comp,y=Envi,fill=Cor),alpha=0.8)+
  facet_wrap(Group~.)+
  #scale_fill_gradient(low = "")
  #scale_fill_gradientn(colors=blue2yellow(2))+
  scale_fill_gradientn(colors=rev(matlab.like(20)))+
  scale_x_discrete(name="")+
  scale_y_discrete(name="",labels=rev(c(expression(bold("PM"["2.5"])),
                            expression(bold("NH"["4"]^{"+"})),
                            expression(bold("NO"["3"]^{"-"})),
                            expression(bold("SO"["4"]^{"2-"})),
                            expression(bold("O"["3"])),
                            expression(bold("CO")),
                            expression(bold("NO")),
                            expression(bold("SO"["2"])))))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.6,0,0.2),"cm")),
        strip.text = element_text(size = 20, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 26, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.key.width = unit(1.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )

p2+ggsave(filename("heat_map_compvsEnvi"),height = 30, width = 45, units = "cm", dpi = 300)

install.packages("colorRamps")
library(colorRamps)
matlab.like2(10)

YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
colorRampPalette(YlOrBr, space = "Lab",bias = 0.5)
#Event vs Non-event=====
##1.Chemical composition====
FT_merge

FT_envi=fread("Datafile/FRIENDs_envi.csv")

ft_inty_all=ft_inty
ft_inty_all$Sample=paste(ft_inty_all$Group,ft_inty_all$No,sep = "_")
ft_inty_all=ft_inty_all %>% inner_join(FT_envi[,c("Sample","Date","Event")]) 


ft_inty_ev=ft_inty_all %>% filter(Event%in%c("Event","Non-event"))
ft_inty_ev

ft_inty_ev$Group=factor(ft_inty_ev$Group,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))

ft_inty_ev$Complab=paste(ft_inty_ev$Comp,"(%)")
ft_inty_ev$Complab=factor(ft_inty_ev$Complab,levels = c("CHO (%)","CHON (%)","CHOS (%)","CHONS (%)"))

val_max_ev=as.data.table(aggregate(ft_inty_ev$rel, 
                                   by=list(`Comp`=ft_inty_ev$Comp,Ev=ft_inty_ev$Event), FUN=max)) %>% 
  `colnames<-`(c("Comp","Event","rel"))
val_max_ev
val_max_ev$Complab=paste(val_max_ev$Comp,"(%)")
val_max_ev$Complab=factor(val_max_ev$Complab,levels = c("CHO (%)","CHON (%)","CHOS (%)","CHONS (%)"))

p1=ggplot(data=ft_inty_ev, aes(x=Group, y=rel, fill=Event))+
  #geom_violin(position = "dodge")+
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.8), width=0.4)+
  geom_boxplot(position = position_dodge(width = 0.80), outlier.colour = NA)+
  facet_wrap(Complab~., strip.position = "left", scales = "free", ncol=2)+
  geom_blank(data=val_max_ev, aes(x=2.5,y=rel*1.12))+
  scale_x_discrete("")+
  scale_fill_manual(values = c("grey75","white"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.0,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.4,0.2,0.4),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 24, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.08,0.94),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )

p1+ggsave(filename("Eventcompare"),height = 30, width = 40, units = "cm", dpi = 300)

##Compare distribution (Event vsNon-event)=====
grp=as.vector(unique(ft_inty_ev$Group))
comp=as.vector(unique(ft_inty_ev$Comp))

data=ft_inty_ev
#data=subset(ft_inty_ev,ft_inty_ev$Sample!="Ulaanbaatar_8") %>% droplevels()

dt=data.table()
for (i in 1:length(comp)) {
#  i=3
  temp=subset(data,data$Comp==comp[i])
  for (j in 1:length(grp)) {
#    j=5
    temp2=subset(temp,temp$Group==grp[j])
    
    t_eve=subset(temp2,temp2$Event=="Event")
    t_nev=subset(temp2,temp2$Event=="Non-event")
    
    s1=shapiro.test(t_eve$rel)
    s2=shapiro.test(t_nev$rel)
    
    evtmp=levene.test(temp2$rel,temp2$Event)
    
    tt=t.test(t_eve$rel,t_nev$rel, paired = F,var.equal = F)
    tt2=wilcox.test(t_eve$rel,t_nev$rel, exact = F)
    
    new=data.table("Comp"=comp[i],"Group"=grp[j],"norm_eve"=round(s1$p.value,3),"norm_nev"=round(s2$p.value,3),
                   "Eq-var"=round(evtmp$p,3),"T-test"=round(tt$p.value,3),"non_ttest"=round(tt2$p.value,3))
      
    dt=rbind(dt,new)
    
  }
}
dt

#dt_non_equal=dt
dt
##4-2-2)Calculate label position=====
uf=data.table()
for (i in 1:length(comp)) {
  #i=4
  temp=subset(data,data$Comp==comp[i])
  for (j in 1:length(grp)) {
    #j=5
    temp2=subset(temp,temp$Group==grp[j])
    temp2
    
    t_eve=subset(temp2,temp2$Event=="Event")
    t_nev=subset(temp2,temp2$Event=="Non-event")
    
    iqr.df <- ddply(temp2, "Event", function (x) ((1.5*IQR(x[,"rel"])+quantile(x[,"rel"],0.75)))) %>% 
      `colnames<-`(c("Event","iqr"))
    
    t3=temp2 %>% inner_join(iqr.df)
    
    t4=subset(t3,t3$rel<t3$iqr)
    t4
    
    uf.ps <- ddply(t4, "Event", function (x) (max(fivenum(x[,"rel"]))))
    new=data.table("Comp"=comp[i],"Group"=grp[j],"uf_eve"=uf.ps[1,2],"uf_nev"=uf.ps[2,2] )
    uf=rbind(uf,new)
    
  }
}

uf$uf_pos=ifelse(uf$uf_eve>=uf$uf_nev,uf$uf_eve,uf$uf_nev)
uf$Group=factor(uf$Group,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))


uf$id=paste(uf$Comp,uf$Group,sep = "_")
dt$id=paste(dt$Comp,dt$Group,sep = "_")
#dt$Noneq_t=dt_non_equal$`T-test`
#dt$p=ifelse(dt$`T-test`<dt$Noneq_t,dt$`T-test`,dt$Noneq_t)
dt
uf=uf %>% inner_join(dt[,c("id","non_ttest")])
uf$sig=ifelse(uf$non_ttest<0.001,"***",
              ifelse(uf$non_ttest<0.01,"**",
                     ifelse(uf$non_ttest<0.05,"*","N.S")))

#uf$Group=factor(uf$Group,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))
uf$Comp=factor(uf$Comp,levels = c("CHO","CHON","CHOS","CHONS"))

uf$xs=ifelse(uf$Group=="Seoul",3,
             ifelse(uf$Group=="Seosan",4,
                    ifelse(uf$Group=="Beijing",2,
                           ifelse(uf$Group=="Noto",5,1))))

uf$Complab=paste(uf$Comp,"(%)")
uf$Complab=factor(uf$Complab,levels = c("CHO (%)","CHON (%)","CHOS (%)","CHONS (%)"))

ps=ggplot()+
  #geom_violin(position = "dodge")+
  stat_boxplot(data=ft_inty_ev, aes(x=Group, y=rel, fill=Event),geom = "errorbar",position = position_dodge(width = 0.8), width=0.4)+
  #stat_boxplot(data=subset(ft_inty_ev,ft_inty_ev$Sample!="Ulaanbaatar_8"), aes(x=Group, y=rel, fill=Event),geom = "errorbar",position = position_dodge(width = 0.8), width=0.4)+
  geom_boxplot(data=ft_inty_ev, aes(x=Group, y=rel, fill=Event),position = position_dodge(width = 0.80), outlier.colour = NA)+
  #geom_boxplot(data=subset(ft_inty_ev,ft_inty_ev$Sample!="Ulaanbaatar_8"), aes(x=Group, y=rel, fill=Event),position = position_dodge(width = 0.80), outlier.colour = NA)+
  facet_wrap(Complab~., strip.position = "left", scales = "free", ncol=2)+
  geom_blank(data=val_max_ev, aes(x=2.5,y=rel*1.12))+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), lty=2)+
  geom_segment(data=uf,aes(x=xs-0.2,xend=xs+0.2,y=uf_pos*1.03, yend=uf_pos*1.03))+
  #geom_segment(data=uf,aes(x=xs-0.2,xend=xs-0.2,y=uf_eve*1.01, yend=uf_pos*1.05))+
  #geom_segment(data=uf,aes(x=xs+0.2,xend=xs+0.2,y=uf_nev*1.01, yend=uf_pos*1.05))+
  #geom_text(data =comp_max, aes(x=1, y=x, label=Comp), col="white", size=0)+
  geom_text(data =subset(uf,uf$non_ttest<0.05), aes(x=Group, y=uf_pos*1.05, label=sig), col="black", size=6)+
  scale_x_discrete("")+
  scale_fill_manual(values = c("grey75","white"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.0,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.4,0.2,0.4),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 24, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.08,0.94),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_rect(fill="white")
  )

ps+ggsave(filename("Eventcompare_sig"),height = 30, width = 40, units = "cm", dpi = 300)

##molecular comparision======
FT_mole=fread("Datafile/FRIEND_molecule.csv")

FT_mole$Sample=paste(FT_mole$Group,FT_mole$No,sep = "_")

FT_mole_all=melt(FT_mole,id.vars = c("Sample","Group","No"))

FT_mole_all=FT_mole_all %>% inner_join(FT_envi[,c("Sample","Date","Event")]) 

FT_mole_all

FT_mole_ev=FT_mole_all %>% filter(Event%in%c("Event","Non-event"))
FT_mole_ev

FT_mole_ev$Group=factor(FT_mole_ev$Group,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))

FT_mole_ev$varlab=paste(FT_mole_ev$variable,"(%)")

unique(FT_mole_ev$varlab)


FT_mole_ev$varlab=factor(FT_mole_ev$varlab,levels = c("Lignins (%)","Condensed aromatics (%)","Tannins (%)", "Unsaturated hydrocarbons (%)",
                                                       "Carbohydrates (%)","Proteins (%)", "Lipids (%)"))

mol_max_ev=as.data.table(aggregate(FT_mole_ev$value, 
                                   by=list(`variable`=FT_mole_ev$variable,Ev=FT_mole_ev$Event), FUN=max)) %>% 
  `colnames<-`(c("variable","Event","value"))

mol_max_ev

mol_max_ev$varlab=paste(mol_max_ev$variable,"(%)")
mol_max_ev$varlab=factor(mol_max_ev$varlab,levels = c("Lignins (%)","Condensed aromatics (%)","Tannins (%)", "Unsaturated hydrocarbons (%)",
                                                      "Carbohydrates (%)","Proteins (%)", "Lipids (%)"))

p1=ggplot(data=FT_mole_ev, aes(x=Group, y=value, fill=Event))+
  #geom_violin(position = "dodge")+
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.8), width=0.4)+
  geom_boxplot(position = position_dodge(width = 0.80), outlier.colour = NA)+
  facet_wrap(varlab~., strip.position = "left", scales = "free", ncol=4)+
  geom_blank(data=mol_max_ev, aes(x=2.5,y=value*1.12))+
  scale_x_discrete("")+
  scale_fill_manual(values = c("grey75","white"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.0,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.4,0.2,0.4),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 24, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.82,0.94),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )

p1+ggsave(filename("Event_mole_compare"),height = 35, width = 75, units = "cm", dpi = 300)

##Compare distribution (Event vsNon-event)=====
grp=as.vector(unique(FT_mole_ev$Group))
comp=as.vector(unique(FT_mole_ev$variable))

data=FT_mole_ev
#data=subset(FT_mole_ev,FT_mole_ev$Sample!="Ulaanbaatar_8") %>% droplevels()

dt=data.table()
for (i in 1:length(comp)) {
   # i=1
  temp=subset(data,data$variable==comp[i])
  for (j in 1:length(grp)) {
    #    j=1
    temp2=subset(temp,temp$Group==grp[j])
    
    t_eve=subset(temp2,temp2$Event=="Event")
    t_nev=subset(temp2,temp2$Event=="Non-event")
    
    s1=shapiro.test(t_eve$value)
    s2=shapiro.test(t_nev$value)
    
    evtmp=levene.test(temp2$value,temp2$Event)
    
    tt=t.test(t_eve$value,t_nev$value, paired = F,var.equal = F)
    tt2=wilcox.test(t_eve$value,t_nev$value, exact = F)
    
    new=data.table("Comp"=comp[i],"Group"=grp[j],"norm_eve"=round(s1$p.value,3),"norm_nev"=round(s2$p.value,3),
                   "Eq-var"=round(evtmp$p,3),"T-test"=round(tt$p.value,3),"non_ttest"=round(tt2$p.value,3))
    
    dt=rbind(dt,new)
    
  }
}
dt

#dt_non_equal=dt
dt
##4-2-2)Calculate label position=====
ufm=data.table()
for (i in 1:length(comp)) {
  #i=1
  temp=subset(data,data$variable==comp[i])
  for (j in 1:length(grp)) {
    #j=1
    temp2=subset(temp,temp$Group==grp[j])
    temp2
    
    t_eve=subset(temp2,temp2$Event=="Event")
    t_nev=subset(temp2,temp2$Event=="Non-event")
    
    iqr.df <- ddply(temp2, "Event", function (x) ((1.5*IQR(x[,"value"])+quantile(x[,"value"],0.75)))) %>% 
      `colnames<-`(c("Event","iqr"))
    
    t3=temp2 %>% inner_join(iqr.df)
    
    t4=subset(t3,t3$value<t3$iqr)
    t4
    
    uf.ps_max <- ddply(t4, "Event", function (x) (max(fivenum(x[,"value"]))))
    uf.ps_min <- ddply(t4, "Event", function (x) (min(fivenum(x[,"value"]))))
    
    new=data.table("Comp"=comp[i],"Group"=grp[j],"uf_eve"=uf.ps_max[1,2],"uf_nev"=uf.ps_max[2,2],
                   "uf_eve_min"=uf.ps_min[1,2],"uf_nev_min"=uf.ps_min[2,2])
    ufm=rbind(ufm,new)
    
  }
}
ufm
#ufm$uf_eve=(ufm$uf_eve_max-ufm$uf_eve_min)/100*110+ufm$uf_eve_max
#ufm$uf_nev=(ufm$uf_nev_max-ufm$uf_nev_min)/100*110+ufm$uf_nev_max
ufm$uf_pos=ifelse(ufm$uf_eve>=ufm$uf_nev,ufm$uf_eve,ufm$uf_nev)


ufm$Group=factor(ufm$Group,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))


ufm$id=paste(ufm$Comp,uf$Group,sep = "_")
dt$id=paste(dt$Comp,dt$Group,sep = "_")
#dt$Noneq_t=dt_non_equal$`T-test`
#dt$p=ifelse(dt$`T-test`<dt$Noneq_t,dt$`T-test`,dt$Noneq_t)
dt
ufm=ufm %>% inner_join(dt[,c("id","non_ttest")])
ufm$sig=ifelse(ufm$non_ttest<0.001,"***",
              ifelse(ufm$non_ttest<0.01,"**",
                     ifelse(ufm$non_ttest<0.05,"*","N.S")))
ufm


ufm$xs=ifelse(ufm$Group=="Seoul",3,
             ifelse(ufm$Group=="Seosan",4,
                    ifelse(ufm$Group=="Beijing",2,
                           ifelse(ufm$Group=="Noto",5,1))))

ufm$varlab=paste(ufm$Comp,"(%)")

ufm$varlab=factor(ufm$varlab, levels = c("Lignins (%)","Condensed aromatics (%)","Tannins (%)", "Unsaturated hydrocarbons (%)",
                                                 "Carbohydrates (%)","Proteins (%)", "Lipids (%)"))

ufm$uf_pos=ifelse(ufm$Group=="Ulaanbaatar",ufm$uf_pos*1.1,ufm$uf_pos)
ufm$uf_pos=ifelse(ufm$id=="Unsaturated hydrocarbons_Seosan",ufm$uf_pos+0.1,ufm$uf_pos)
ufm$uf_pos=ifelse(ufm$id=="Lipids_Seosan",ufm$uf_pos+1.1,ufm$uf_pos)

unique(ufm$id)
ps=ggplot()+
  #geom_violin(position = "dodge")+
  stat_boxplot(data=FT_mole_ev, aes(x=Group, y=value, fill=Event),geom = "errorbar",position = position_dodge(width = 0.8), width=0.4)+
  #stat_boxplot(data=subset(FT_mole_ev,FT_mole_ev$Sample!="Ulaanbaatar_8"), aes(x=Group, y=rel, fill=Event),geom = "errorbar",position = position_dodge(width = 0.8), width=0.4)+
  geom_boxplot(data=FT_mole_ev, aes(x=Group, y=value, fill=Event),position = position_dodge(width = 0.80), outlier.colour = NA)+
  #geom_boxplot(data=subset(FT_mole_ev,FT_mole_ev$Sample!="Ulaanbaatar_8"), aes(x=Group, y=rel, fill=Event),position = position_dodge(width = 0.80), outlier.colour = NA)+
  facet_wrap(varlab~., strip.position = "left", scales = "free", ncol=4)+
  geom_blank(data=mol_max_ev, aes(x=2.5,y=value*1.12))+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), lty=2)+
  geom_segment(data=ufm,aes(x=xs-0.2,xend=xs+0.2,y=uf_pos*1.03, yend=uf_pos*1.03))+
  #geom_segment(data=uf,aes(x=xs-0.2,xend=xs-0.2,y=uf_eve*1.01, yend=uf_pos*1.05))+
  #geom_segment(data=uf,aes(x=xs+0.2,xend=xs+0.2,y=uf_nev*1.01, yend=uf_pos*1.05))+
  #geom_text(data =comp_max, aes(x=1, y=x, label=Comp), col="white", size=0)+
  geom_text(data =subset(ufm,ufm$non_ttest<0.05), aes(x=Group, y=uf_pos*1.05, label=sig), col="black", size=6)+
  scale_x_discrete("")+
  scale_fill_manual(values = c("grey75","white"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.0,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.4,0.2,0.4),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 24, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.82,0.9525),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_rect(fill="white")
  )

ps+ggsave(filename("Eventcompare_mol_sig"),height = 35, width = 75, units = "cm", dpi = 300)



##molecule====

ft_mole=fread("Datafile/FRIEND_molecule.csv")

ft_mole$Sample=paste(ft_mole$Group,ft_mole$No,sep = "_")
ft_mole

ft_mole_m=melt(ft_mole, id.vars = c("Sample","Group","No"))
ft_mole_m

ft_mole_envi=ft_mole %>% inner_join(envi_all_sel[,-c("Group","No")])
ft_mole_envi

ft_mole_envi=ft_mole_envi %>% filter(Group%in%c("Seoul","Beijing"))
ft_mole_envi

wsoa=fread("Datafile/WOSCbb_S&B.csv")

ft_mole_envi=ft_mole_envi %>% inner_join(wsoa[,c("Sample","WSOCbb","WSOCnbb")])



ft_mole_envi_m=melt(ft_mole_envi,id.vars = c("Sample","Group","No","Lignins","Condensed aromatics","OC",
                                             "Tannins","Carbohydrates","Unsaturated hydrocarbons","Lipids","Proteins"),variable.name = "Envi", value.name = "Conc") %>% 
  melt(id.vars=c("Sample","Group","No","OC","Envi","Conc"),variable.name="Comp")
  
ft_mole_envi_m

ft_mole_envi_m$No=as.numeric(ft_mole_envi_m$No)
ft_mole_envi_m$mark="no"
ft_mole_envi_m$mark=ifelse(ft_mole_envi_m$No==5,"12/19",ft_mole_envi_m$mark)
ft_mole_envi_m$mark=ifelse(ft_mole_envi_m$No==6,"12/20",ft_mole_envi_m$mark)
ft_mole_envi_m$mark=ifelse(ft_mole_envi_m$No==7,"12/21",ft_mole_envi_m$mark)
ft_mole_envi_m$mark=ifelse(ft_mole_envi_m$No==8,"12/22",ft_mole_envi_m$mark)
ft_mole_envi_m$mark=ifelse(ft_mole_envi_m$No==9,"12/23",ft_mole_envi_m$mark)

ft_mole_envi_m
table(ft_mole_envi_m$Envi)

ft_mole_envi_m_s=subset(ft_mole_envi_m,ft_mole_envi_m$Group=="Seoul")
ft_mole_envi_m_s1=ft_mole_envi_m_s %>% filter(Envi%in%c("PM2.5","NH4","NO3","SO4")) %>% droplevels()

table(ft_mole_envi_m_s1$Envi)

ft_mole_envi_m_s1$Envilab=ft_mole_envi_m_s1$Envi
ft_mole_envi_m_s1$Envilab=factor(ft_mole_envi_m_s1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                                labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))
                                           ))

unique(ft_mole_envi_m_s1$Comp)
ft_mole_envi_m_s1$Complab=ft_mole_envi_m_s1$Comp

ft_mole_envi_m_s1$Complab=factor(ft_mole_envi_m_s1$Complab, 
                                 levels = c("Lignins","Tannins","Condensed aromatics","Unsaturated hydrocarbons",
                                                       "Lipids","Proteins","Carbohydrates"),
                                 labels = c(expression(bold("Lignins")),expression(bold("Tannins")),
                                            expression(bold("Condensed \n aromatics")),
                                            expression(bold("Unsaturated\nhydrocarbons")),
                                            expression(bold("Lipids")),
                                            expression(bold("Proteins")),
                                            expression(bold("Carbohydrates"))))
                                            
mole_max_s1=as.data.table(aggregate(ft_mole_envi_m_s1$Conc, 
                                   by=list(`Comp`=ft_mole_envi_m_s1$Comp,Envi=ft_mole_envi_m_s1$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))

mole_max_s1$Envilab=mole_max_s1$Envi

mole_max_s1$Envilab=factor(mole_max_s1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                                 labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                            expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                            expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                            expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))
                                 ))

mole_max_s1


ft_mole_envi_m_s1

ggplot()+
  geom_smooth(data=ft_mole_envi_m_s1,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=ft_mole_envi_m_s1,aes(x=value,y=Conc),col="royalblue2", size=3)+
  geom_point(data=ft_mole_envi_m_s1,aes(x=value,y=Conc*1.5),col=NA, size=0)+
  #geom_blank(data=mole_max_s1, aes(x=2.5,y=Conc*1.4))+
  geom_point(data=subset(ft_mole_envi_m_s1,ft_mole_envi_m_s1$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_mole_envi_m_s1,ft_mole_envi_m_s1$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  facet_grid(Envilab~Complab, scales = "free",labeller = label_parsed, switch = "y")+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of molecular class (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20,vjust=0 ,hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(2.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("molevsEnvi_sul1"),height = 30, width = 52.5, units = "cm", dpi = 300)


ft_mole_envi_m_s3=ft_mole_envi_m_s %>% filter(Envi%in%c("WSOCbb","WSOCnbb")) %>% droplevels()

table(ft_mole_envi_m_s1$Envi)

ft_mole_envi_m_s3$Envilab=ft_mole_envi_m_s3$Envi
ft_mole_envi_m_s3$Envilab=factor(ft_mole_envi_m_s3$Envilab, levels = c("WSOCbb","WSOCnbb"),
                                 labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                            expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")"))
                                 ))

ft_mole_envi_m_s3$Complab=ft_mole_envi_m_s3$Comp

ft_mole_envi_m_s3$Complab=factor(ft_mole_envi_m_s3$Complab, 
                                 levels = c("Lignins","Tannins","Condensed aromatics","Unsaturated hydrocarbons",
                                            "Lipids","Proteins","Carbohydrates"),
                                 labels = c(expression(bold("Lignins")),expression(bold("Tannins")),
                                            expression(bold("Condensed \n aromatics")),
                                            expression(bold("Unsaturated\nhydrocarbons")),
                                            expression(bold("Lipids")),
                                            expression(bold("Proteins")),
                                            expression(bold("Carbohydrates"))))

mole_max_s3=as.data.table(aggregate(ft_mole_envi_m_s3$Conc, 
                                    by=list(`Comp`=ft_mole_envi_m_s3$Comp,Envi=ft_mole_envi_m_s3$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))

mole_max_s3$Envilab=mole_max_s3$Envi
mole_max_s3$Envilab=factor(mole_max_s3$Envilab, levels = c("WSOCbb","WSOCnbb"),
                                 labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                            expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")"))
                                 ))
mole_max_s3


ft_mole_envi_m_s3

ggplot()+
  geom_smooth(data=ft_mole_envi_m_s3,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=ft_mole_envi_m_s3,aes(x=value,y=Conc),col="royalblue2", size=3)+
  geom_point(data=ft_mole_envi_m_s3,aes(x=value,y=Conc*1.5),col=NA, size=0)+
  #geom_blank(data=mole_max_s1, aes(x=2.5,y=Conc*1.4))+
  geom_point(data=subset(ft_mole_envi_m_s3,ft_mole_envi_m_s3$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_mole_envi_m_s3,ft_mole_envi_m_s3$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  facet_grid(Envilab~Complab, scales = "free",labeller = label_parsed, switch = "y")+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of molecular class (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20,vjust=0 ,hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(2.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("molevsEnvi_sul3"),height = 20, width = 52.5, units = "cm", dpi = 300)





##Beijing====
ft_mole_envi_m_b=subset(ft_mole_envi_m,ft_mole_envi_m$Group=="Beijing")
ft_mole_envi_m_b1=ft_mole_envi_m_b %>% filter(Envi%in%c("PM2.5","NH4","NO3","SO4")) %>% droplevels()

table(ft_mole_envi_m_b1$Envi)

ft_mole_envi_m_b1$Envilab=ft_mole_envi_m_b1$Envi
ft_mole_envi_m_b1$Envilab=factor(ft_mole_envi_m_b1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                                 labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                            expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                            expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                            expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))
                                 ))

unique(ft_mole_envi_m_b1$Comp)
ft_mole_envi_m_b1$Complab=ft_mole_envi_m_b1$Comp

ft_mole_envi_m_b1$Complab=factor(ft_mole_envi_m_b1$Complab, 
                                 levels = c("Lignins","Tannins","Condensed aromatics","Unsaturated hydrocarbons",
                                            "Lipids","Proteins","Carbohydrates"),
                                 labels = c(expression(bold("Lignins")),expression(bold("Tannins")),
                                            expression(bold("Condensed \n aromatics")),
                                            expression(bold("Unsaturated\nhydrocarbons")),
                                            expression(bold("Lipids")),
                                            expression(bold("Proteins")),
                                            expression(bold("Carbohydrates"))))

mole_max_b1=as.data.table(aggregate(ft_mole_envi_m_b1$Conc, 
                                    by=list(`Comp`=ft_mole_envi_m_b1$Comp,Envi=ft_mole_envi_m_b1$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))

mole_max_b1$Envilab=mole_max_b1$Envi

mole_max_b1$Envilab=factor(mole_max_b1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                           labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                      expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                      expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                      expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))
                           ))

mole_max_b1


ft_mole_envi_m_b1

ggplot()+
  geom_smooth(data=ft_mole_envi_m_b1,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=ft_mole_envi_m_b1,aes(x=value,y=Conc),col="orangered3", size=3)+
  geom_point(data=ft_mole_envi_m_b1,aes(x=value,y=Conc*1.5),col=NA, size=0)+
  #geom_blank(data=mole_max_s1, aes(x=2.5,y=Conc*1.4))+
  geom_point(data=subset(ft_mole_envi_m_b1,ft_mole_envi_m_b1$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_mole_envi_m_b1,ft_mole_envi_m_b1$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  facet_grid(Envilab~Complab, scales = "free",labeller = label_parsed, switch = "y")+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of molecular class (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20,vjust=0 ,hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(2.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("molevsEnvi_b1"),height = 30, width = 52.5,units = "cm", dpi = 300)



    ft_mole_envi_m_b3=ft_mole_envi_m_b %>% filter(Envi%in%c("WSOCbb","WSOCnbb")) %>% droplevels()

table(ft_mole_envi_m_b3$Envi)

ft_mole_envi_m_b3$Envilab=ft_mole_envi_m_b3$Envi
ft_mole_envi_m_b3$Envilab=factor(ft_mole_envi_m_b3$Envilab, levels = c("WSOCbb","WSOCnbb"),
                                 labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                            expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")"))
                                 ))

ft_mole_envi_m_b3$Complab=ft_mole_envi_m_b3$Comp

ft_mole_envi_m_b3$Complab=factor(ft_mole_envi_m_b3$Complab, 
                                 levels = c("Lignins","Tannins","Condensed aromatics","Unsaturated hydrocarbons",
                                            "Lipids","Proteins","Carbohydrates"),
                                 labels = c(expression(bold("Lignins")),expression(bold("Tannins")),
                                            expression(bold("Condensed \n aromatics")),
                                            expression(bold("Unsaturated\nhydrocarbons")),
                                            expression(bold("Lipids")),
                                            expression(bold("Proteins")),
                                            expression(bold("Carbohydrates"))))

mole_max_b3=as.data.table(aggregate(ft_mole_envi_m_b3$Conc, 
                                    by=list(`Comp`=ft_mole_envi_m_b3$Comp,Envi=ft_mole_envi_m_b3$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))

mole_max_b3$Envilab=mole_max_b3$Envi
mole_max_b3$Envilab=factor(mole_max_b3$Envilab, levels = c("WSOCbb","WSOCnbb"),
                           labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                      expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")"))
                           ))
mole_max_b3


ft_mole_envi_m_b3

ggplot()+
  geom_smooth(data=ft_mole_envi_m_b3,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=ft_mole_envi_m_b3,aes(x=value,y=Conc),col="orangered3", size=3)+
  geom_point(data=ft_mole_envi_m_b3,aes(x=value,y=Conc*1.5),col=NA, size=0)+
  #geom_blank(data=mole_max_s1, aes(x=2.5,y=Conc*1.4))+
  geom_point(data=subset(ft_mole_envi_m_b3,ft_mole_envi_m_b3$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_mole_envi_m_b3,ft_mole_envi_m_b3$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  facet_grid(Envilab~Complab, scales = "free",labeller = label_parsed, switch = "y")+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of molecular class (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20,vjust=0 ,hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(2.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("molevsEnvi_b3"),height = 20, width = 52.5, units = "cm", dpi = 300)



ft_mole_envi_m

grp=unique(ft_mole_envi_m$Group)
grp
comp=unique(ft_mole_envi_m$Comp)
envi=unique(ft_mole_envi_m$Envi)


dt=data.table()
for (i in 1:length(grp)) {
  temp=subset(ft_mole_envi_m,ft_mole_envi_m$Group==grp[i])
  
  for (k in 1:length(comp)) {
    tmp=subset(temp,temp$Comp==comp[k])
    
    for (l in 1:length(envi)) {
      tmp2=subset(tmp,tmp$Envi==envi[l])
      
      if(envi[l]=="O3"&grp[i]=="Ulaanbaatar") {
        tmp2$Conc=0
        new_row=data.table("Group"=grp[i],"Comp"=comp[k],"R2"=0,"Envi"=envi[l],"Cor"=0, "P"=1)
        
      } else {
        model=lm(Conc~value,tmp2)
        a=summary(model)
        cor=cor.test(tmp2$value,tmp2$Conc, method = "pearson")
        #cor=cor.test(tmp2$value,tmp2$Conc, method = "spearman", exact = F)
        #    cor
        
        new_row=data.table("Group"=grp[i],"Comp"=comp[k],"R2"=round(a$r.squared,8),"Envi"=envi[l],"Cor"=cor$estimate, "P"=round(cor$p.value,3))
        
      }
      
      dt=rbind(dt,new_row)
    }
  }
}

dt

dt$sig=ifelse(dt$P<0.05,ifelse(dt$P<0.01,ifelse(dt$P<0.001,"***","**"),"*"),"")


fwrite(dt,"Datafile/211122molecule&Envicor_pearson.csv")


ft_mole_envi_m_3=ft_mole_envi_m %>% filter(Envi%in%c("WSOCbb","WSOCnbb")) %>% droplevels()

table(ft_mole_envi_m_s1$Envi)

ft_mole_envi_m_3$Envilab=ft_mole_envi_m_3$Envi
ft_mole_envi_m_3$Envilab=factor(ft_mole_envi_m_3$Envilab, levels = c("WSOCbb","WSOCnbb"),
                                 labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                            expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")"))
                                 ))

ft_mole_envi_m_3$Complab=ft_mole_envi_m_3$Comp

ft_mole_envi_m_3$Complab=factor(ft_mole_envi_m_3$Complab, 
                                 levels = c("Lignins","Tannins","Condensed aromatics","Unsaturated hydrocarbons",
                                            "Lipids","Proteins","Carbohydrates"),
                                 labels = c(expression(bold("Lignins")),expression(bold("Tannins")),
                                            expression(bold("Condensed \n aromatics")),
                                            expression(bold("Unsaturated\nhydrocarbons")),
                                            expression(bold("Lipids")),
                                            expression(bold("Proteins")),
                                            expression(bold("Carbohydrates"))))

ft_mole_envi_m_3

ft_mole_envi_m_3$Group=factor(ft_mole_envi_m_3$Group, levels = c("Seoul","Beijing"))



ggplot()+
  geom_smooth(data=ft_mole_envi_m_3,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=ft_mole_envi_m_3,aes(x=value,y=Conc,col=Group), size=3)+
  geom_point(data=ft_mole_envi_m_3,aes(x=value,y=Conc*1.5),col=NA, size=0)+
  #geom_blank(data=mole_max_s1, aes(x=2.5,y=Conc*1.4))+
  geom_point(data=subset(ft_mole_envi_m_3,ft_mole_envi_m_3$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_mole_envi_m_3,ft_mole_envi_m_3$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  facet_grid(Group*Envilab~Complab, scales = "free",labeller = label_parsed, switch = "y")+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of molecular class (%)")+
  scale_y_continuous(name = "")+
  scale_color_manual(values = c("royalblue2","orangered3"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20,vjust=0 ,hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(2.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.key.height = unit(1.25,"cm"),
        legend.key.width = unit(0.75,"cm")
  )+
  guides(col=guide_legend(override.aes = list(size=8)))+
  ggsave(filename("molevsEnvi_all3"),height = 40, width = 52.5, units = "cm", dpi = 300)





ft_mole_envi_m_b3=ft_mole_envi_m_b %>% filter(Envi%in%c("WSOCbb","WSOCnbb")) %>% droplevels()

table(ft_mole_envi_m_b3$Envi)

ft_mole_envi_m_b3$Envilab=ft_mole_envi_m_b3$Envi
ft_mole_envi_m_b3$Envilab=factor(ft_mole_envi_m_b3$Envilab, levels = c("WSOCbb","WSOCnbb"),
                                 labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                            expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")"))
                                 ))

ft_mole_envi_m_b3$Complab=ft_mole_envi_m_b3$Comp

ft_mole_envi_m_b3$Complab=factor(ft_mole_envi_m_b3$Complab, 
                                 levels = c("Lignins","Tannins","Condensed aromatics","Unsaturated hydrocarbons",
                                            "Lipids","Proteins","Carbohydrates"),
                                 labels = c(expression(bold("Lignins")),expression(bold("Tannins")),
                                            expression(bold("Condensed \n aromatics")),
                                            expression(bold("Unsaturated\nhydrocarbons")),
                                            expression(bold("Lipids")),
                                            expression(bold("Proteins")),
                                            expression(bold("Carbohydrates"))))

mole_max_b3=as.data.table(aggregate(ft_mole_envi_m_b3$Conc, 
                                    by=list(`Comp`=ft_mole_envi_m_b3$Comp,Envi=ft_mole_envi_m_b3$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))

mole_max_b3$Envilab=mole_max_b3$Envi
mole_max_b3$Envilab=factor(mole_max_b3$Envilab, levels = c("WSOCbb","WSOCnbb"),
                           labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                      expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")"))
                           ))
mole_max_b3


ft_mole_envi_m_b3

ggplot()+
  geom_smooth(data=ft_mole_envi_m_b3,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=ft_mole_envi_m_b3,aes(x=value,y=Conc),col="orangered3", size=3)+
  geom_point(data=ft_mole_envi_m_b3,aes(x=value,y=Conc*1.5),col=NA, size=0)+
  #geom_blank(data=mole_max_s1, aes(x=2.5,y=Conc*1.4))+
  geom_point(data=subset(ft_mole_envi_m_b3,ft_mole_envi_m_b3$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_mole_envi_m_b3,ft_mole_envi_m_b3$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  facet_grid(Envilab~Complab, scales = "free",labeller = label_parsed, switch = "y")+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of molecular class (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20,vjust=0 ,hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(2.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("molevsEnvi_b3"),height = 20, width = 52.5, units = "cm", dpi = 300)


##cor heatmap=====


heat=fread("Datafile/211122FT&Envicor_all.csv")
heat




ggplot(heat)+
  geom_tile(data=subset(heat,heat$P<0.1), aes(x=Comp,y=Envi,fill=Cor),alpha=0.8)+
  facet_wrap(Group~., ncol=1, scales = "free")+
  #scale_fill_gradient(low = "")
  #scale_fill_gradientn(colors=blue2yellow(2))+
  scale_fill_gradientn(colors=(topo.colors(40)[4:36]))+
  scale_x_discrete(name="")+
  #scale_y_discrete(name="",labels=rev(c(expression(bold("PM"["2.5"])),
  #                                      expression(bold("NH"["4"]^{"+"})),
  #                                      expression(bold("NO"["3"]^{"-"})),
  #                                      expression(bold("SO"["4"]^{"2-"})),
  #                                      expression(bold("O"["3"])),
  #                                      expression(bold("CO")),
  #                                      expression(bold("NO")),
  #                                      expression(bold("SO"["2"])))))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.6,0,0.2),"cm")),
        strip.text = element_text(size = 20, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 26, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.key.width = unit(1.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("heat_map_compvsEnvi"),height = 30, width = 45, units = "cm", dpi = 300)




heat$Envi=factor(heat$Envi, levels = rev(c("PM2.5","NH4","NO3","SO4","O3","CO","NO","SO2","NO3/OC","SO4/OC","NO3/SO4",
                                             "WSOCbb","WSOCnbb")))
heat$Complab=heat$Comp

heat$Complab=factor(heat$Complab, levels = c("SOA","CHO","CHON","CHOS","CHONS","Lignins","Tannins","Condensed aromatics","Unsaturated hydrocarbons",
                                       "Lipids","Proteins","Carbohydrates"),
                    labels=c("SOA","CHO","CHON","CHOS","CHONS","Lignins","Tannins","Condensed\naromatics","Unsaturated\nhydrocarbons",
                               "Lipids","Proteins","Carbohydrates"))
lab=c(expression(bold("PM"["2.5"])),
  expression(bold("NH"["4"]^{"+"})),
  expression(bold("NO"["3"]^{"-"})),
  expression(bold("SO"["4"]^{"2-"})),
  expression(bold("O"["3"]~"(ppb)")),
  expression(bold("CO (ppb)")),
  expression(bold("NO (ppb)")),
  expression(bold("SO"["2"]~"(ppb)")),
  expression(bold("NO"["3"]^{"-"}*"/OC")),
  expression(bold("SO"["4"]^{"2-"}*"/OC")),
  expression(bold("NO"["3"]^{"-"}*"/"*"SO"["4"]^{"2-"})),
  expression(bold("WSOC"["bb"])),
  expression(bold("WSOC"["nbb"])))

heat_s=subset(heat,heat$Group=="Seoul")
heat
ggplot()+
  geom_tile(data=heat_s, aes(x=Complab,y=Envi,fill=Cor),alpha=0.8)+
  geom_tile(data=subset(heat_s,heat_s$P>0.05), aes(x=Complab,y=Envi),fill="white",alpha=1)+
  facet_wrap(Group~., ncol=1, scales = "fixed")+
  #scale_fill_gradient(low = "")
  #scale_fill_gradientn(colors=blue2yellow(2))+
  scale_fill_gradientn(colors=(topo.colors(40)[6:34]))+
  scale_x_discrete(name="")+
  scale_y_discrete(name="",labels=rev(lab)
  )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16,hjust=0,vjust=1, colour = "black",face = "bold", family = "Arial", angle=315),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.6,0,0.2),"cm")),
        strip.text = element_text(size = 20, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 26, hjust=0.6,colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.4,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.key.width = unit(1.0,"cm"),
        legend.title = element_text(size = 26,hjust=0.2, colour = "black", family = "Arial",face = "bold"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.box.margin=margin(10,0,10,30),
  )+
  guides(fill=guide_colorbar(title = expression(bolditalic("ρ")),ticks.colour = "black", ticks.linewidth = 1.5,
                            barwidth = 2.5, barheight = 25))+
  ggsave(filename("heat_map_compvsEnvi_s"),height = 20, width = 35, units = "cm", dpi = 300)




heat_b=subset(heat,heat$Group=="Beijing")
heat
ggplot()+
  geom_tile(data=heat_b, aes(x=Complab,y=Envi,fill=Cor),alpha=0.8)+
  geom_tile(data=subset(heat_b,heat_b$P>0.05), aes(x=Complab,y=Envi),fill="white",alpha=1)+
  facet_wrap(Group~., ncol=1, scales = "fixed")+
  #scale_fill_gradient(low = "")
  #scale_fill_gradientn(colors=blue2yellow(2))+
  scale_fill_gradientn(colors=(topo.colors(40)[6:34]))+
  scale_x_discrete(name="")+
  scale_y_discrete(name="",labels=rev(lab)
  )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16,hjust=0,vjust=1, colour = "black",face = "bold", family = "Arial", angle=315),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.6,0,0.2),"cm")),
        strip.text = element_text(size = 20, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 26, hjust=0.6,colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        #legend.spacing = unit(0.4,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.key.width = unit(1.0,"cm"),
        legend.title = element_text(size = 26,hjust=0.2, colour = "black", family = "Arial",face = "bold"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box.margin=margin(10,0,10,30),
        legend.box = "horizontal",
        #legend.margin = unit(c(0.4,0.4,0.4,1.4),"cm"),
  )+
  guides(fill=guide_colorbar(title = expression(bolditalic("ρ")),ticks.colour = "black", ticks.linewidth = 1.5,
                             barwidth = 2.5, barheight = 25))+
  ggsave(filename("heat_map_compvsEnvi_b"),height = 20, width = 35, units = "cm", dpi = 300)
