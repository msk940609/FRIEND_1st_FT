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
library(ade4)
library(adespatial)
getOption("digits")
options("digits" = 15)

source("Rscript/func_filename.R")
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")
insert_minor_1 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_2 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}

###molecular_richness=====
FT_envi=fread("Datafile/FRIEND_1st_envi_re.csv")
ft_merge=fread("Datafile/FRIENDs_1st_FT.csv")
table(ft_merge$Group)
ft_merge$Group=factor(ft_merge$Group, levels =c("UL","B","SS","SUL","NT"),
                      labels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))
table(ft_merge$Group)

fm_obs=ft_merge[,c("Group","No","Formula","Freq")] %>% 
  dcast(Group~Formula,value.var = "Freq",sum) %>% 
  melt(id.vars=c("Group")) %>% `colnames<-`(c("Group","Formula","cnt"))

fm_obs_sel=subset(fm_obs,fm_obs$cnt>0)
fm_obs_sel

fm_obs_sel$Group2=factor(fm_obs_sel$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                               labels = c(expression(bold("Ulaanbaatar")),
                                          expression(bold("Beijing")),
                                          expression(bold("Seosan")),
                                          expression(bold("Seoul")),
                                          expression(bold("Noto"))
                               ))


ggplot(fm_obs_sel,aes(x=cnt))+
  geom_histogram(binwidth = 1, col="black",fill="grey70")+
  scale_y_continuous(name = "Occurrence", expand = c(0.05,0.05))+
  scale_x_continuous(name = "Number of days", expand = c(0.05,0.05), breaks = seq(1,30,5))+
  facet_rep_wrap(Group2~., scales = "free_x",ncol=5,dir = "h",repeat.tick.labels = T,labeller = label_parsed, strip.position = c("top"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.15,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.1,0.1,0.1),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 20, face=2 ,colour = "black",margin = unit(c(0.1,0.2,0.1,0.2),"cm"),family = "Arial"),
        axis.title.x = element_text(size = 20, face=2 ,colour = "black",margin = unit(c(0.1,0.2,0.1,0.2),"cm"),family = "Arial"),
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
  ggsave(filename("richnees_hist"), height = 10, width = 50, units = "cm", dpi = 300)

ft_tot=as.data.table(aggregate(ft_merge$Bromo.Inty,by=list(Sample=ft_merge$Sample), sum)) %>% `colnames<-`(c("Sample","tot"))
ft2=ft_merge %>% inner_join(ft_tot)
ft2$rel=ft2$Bromo.Inty/ft2$tot*100

fm_bromo=ft2[,c("Group","No","Formula","Bromo.Inty")] %>% 
  dcast(Group~Formula,value.var = "Bromo.Inty",mean) %>% 
  melt(id.vars=c("Group"),na.rm = T) %>% `colnames<-`(c("Group","Formula","Bromo.Inty"))

fm_bromo=ft2[,c("Group","No","Formula","rel")] %>% 
  dcast(Group~Formula,value.var = "rel",mean) %>% 
  melt(id.vars=c("Group"),na.rm = T) %>% `colnames<-`(c("Group","Formula","rel"))

fm_bromo[,1:3]
fm_obs

fm_bromo=fm_bromo %>% inner_join(fm_obs)

fm_bromo_sel=fm_bromo
fm_bromo_sel$Bromo2=fm_bromo_sel$Bromo.Inty*10000

fm_bromo_sel$Group2=factor(fm_bromo_sel$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                         labels = c(expression(bold("Ulaanbaatar")),
                                    expression(bold("Beijing")),
                                    expression(bold("Seosan")),
                                    expression(bold("Seoul")),
                                    expression(bold("Noto"))
                         ))

ggplot(fm_bromo_sel,aes(x=as.factor(cnt), y=log10(Bromo2)))+
  geom_boxplot(colour="black",outlier.shape = 21,outlier.fill = "white", outlier.size = 1)+
  #geom_histogram(binwidth = 1, col="black",fill="grey70")+
  scale_y_continuous(name = expression(bold(Log["10"]~"Mean intensity of AOM")), expand = c(0.05,0.05))+
  scale_x_discrete(name = "Number of days", expand = c(0.05,0.05), breaks = seq(1,30,5))+
  facet_rep_wrap(Group2~., scales = "free",ncol=5,dir = "h",repeat.tick.labels = T,labeller = label_parsed, strip.position = c("top"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.15,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.1,0.1,0.1),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 17, face=2 ,colour = "black",vjust = 0.5, hjust=0.6,margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial"),
        axis.title.x = element_text(size = 17, face=2 ,colour = "black",margin = unit(c(0.1,0.2,0.1,0.2),"cm"),family = "Arial"),
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
  #coord_cartesian(ylim = c(0,30))+
  ggsave(filename("Meanadundance_hist"), height = 10, width = 50, units = "cm", dpi = 300)




FT_envi
ft_merge$Freq=1

ft_merge$Sample=paste(ft_merge$Group,ft_merge$No,sep = "_")
ft_merge

ft_merge_sel=ft_merge %>% left_join(fm_obs)
ft_merge_sel=ft_merge_sel %>% filter(cnt>1)

#Frequency=====
mol_rich=as.data.table(aggregate(ft_merge_sel$Freq,by=list(Sample=ft_merge_sel$Sample),sum))
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

##cumulative curve=====
mol_rich
mol_rich_n

mol_acum=mol_rich_n[,c("Sample","Group","x","PM2.5")]

grp=unique(mol_acum$Group)
d_acu=data.table()
for (j in 1:length(grp)) {
  #j=2
  
  mol_acum_tmp=subset(mol_acum,mol_acum$Group==grp[j])
  mol_acum_tmp=mol_acum_tmp[order(mol_acum_tmp$PM2.5),]
  mol_acum_tmp$id=seq(1, nrow(mol_acum_tmp),by=1)
  
  ord=as.vector(mol_acum_tmp$Sample)
  dt=data.table()
  
  for (i in 1:length(ord)) {
    #i=2
    
    temp=subset(ft_merge,ft_merge$Sample==ord[i])
    temp_fm=temp[,c("Group","Formula")]
    dt=rbind(dt,temp_fm)
    
    dt=unique(dt)
    mol_temp=subset(mol_acum_tmp,mol_acum_tmp$Sample==ord[i])
    new=data.table(Sample=ord[i],Group=unique(dt$Group), mole_acu=dim(dt)[1], No=i, var=mol_temp$PM2.5)
    
    d_acu=rbind(d_acu, new)
    
  }
  
  
}

d_acu

table(d_acu$Group)


#d_acu2=d_acu %>% inner_join(FT_envi[,c("Sample","PM2.5")])
d_acu_oc=d_acu %>% inner_join(FT_envi[,c("Sample","PM2.5")])

#expression(bold("OC"~"("*"\u03bcg/"*m^"3"*")")),



ggplot()+
  geom_line(data = d_acu_oc ,aes(x=var   , y=mole_acu, col=Group),size=1, na.rm = F)+
  geom_point(data = d_acu_oc ,aes(x=var   , y=mole_acu, col=Group),size=2, na.rm = F)+
  geom_smooth(data = d_acu_oc ,aes(x=var   , y=mole_acu, group=Group),,col="black", lty=2, method = "lm",formula=y~log(x))+
  scale_y_continuous(name = "Cumulative number of AOM")+
  scale_x_continuous(name = expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")))+
  #scale_x_continuous(name = expression(bold("OC"~"("*"\u03bcg/"*m^"3"*")")))+
  #scale_x_continuous(name = "Number of days", expand = c(0.05,0.05), breaks = seq(1,30,5))+
  facet_wrap(Group~., scales = "free",ncol=5,dir = "h", strip.position = c("top"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, colour = "NA"),
        strip.text = element_text(colour = "black", size = 20, face = "bold",margin = unit(c(0.0,0.2,0.2,0.2),"cm")),
        plot.title= element_text(size = 24, colour = "black", face="bold",margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2, colour = "black"),
        axis.ticks.length = unit(0.15,"cm"),
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.1,0.1,0.1),"cm")),
        axis.text.x = element_text(size = 16, colour = "black" , face = "bold",family = "Arial",margin = unit(c(0.1,0.2,0.1,0.2),"cm")),
        axis.title.y.left = element_text(size = 17, face=2 ,colour = "black",vjust = 0.5, hjust=0.6,margin = unit(c(0.1,0.1,0.1,0.1),"cm"),family = "Arial"),
        axis.title.x = element_text(size = 17, face=2 ,colour = "black",margin = unit(c(0.1,0.2,0.1,0.2),"cm"),family = "Arial"),
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
  ggsave(filename("cumulaitive_richness_PM"), height = 11, width = 50, units = "cm", dpi = 300)

#richness olot=====
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
  #facet_wrap(Group2~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed)+
  facet_wrap(Group2~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed, strip.position = c("bottom"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
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

ggplot()+
  geom_smooth(data = mol_rich_n_ul ,aes(x=value, y=x),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_rich_n_ul ,aes(x=value, y=x, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_rich_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  scale_y_continuous(name = "Molecular richness")+
  #scale_y_continuous(name = expression(bold("Sum of molecular intensities (x"~"10"^"7"~")")))+
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
        axis.title.y.left = element_text(size = 20, face=2 ,colour = "black",margin = unit(c(0.1,0.3,0.1,0.3),"cm"),family = "Arial"),
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
  geom_smooth(data = mol_rich_n_bj ,aes(x=value, y=x),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_rich_n_bj ,aes(x=value, y=x, col=variable),size=1, na.rm = F)+
  scale_y_continuous(name = "Molecular richness")+
  #scale_y_continuous(name = "Sum of molecular intensities")+
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
        axis.title.y.left = element_text(size = 20, face=2 ,colour = "black",margin = unit(c(0.1,0.3,0.1,0.3),"cm"),family = "Arial"),
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
  geom_smooth(data = mol_rich_n_ss ,aes(x=value, y=x),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_rich_n_ss ,aes(x=value, y=x, col=variable),size=1, na.rm = F)+
  scale_y_continuous(name = "Molecular richness")+
  #scale_y_continuous(name = "Sum of molecular intensities")+
  facet_wrap(.~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed, strip.position = "bottom")+
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
        axis.title.y.left = element_text(size = 20, face=2 ,colour = "black",margin = unit(c(0.1,0.3,0.1,0.3),"cm"),family = "Arial"),
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
  geom_smooth(data = mol_rich_n_sul ,aes(x=value, y=x),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_rich_n_sul ,aes(x=value, y=x, col=variable),size=1, na.rm = F)+
  scale_y_continuous(name = "Molecular richness")+
  #scale_y_continuous(name = "Sum of molecular intensities")+
  facet_wrap(.~variable2, scales = "free",ncol=4,dir = "h",labeller = label_parsed, strip.position = "bottom")+
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
        axis.title.y.left = element_text(size = 20, face=2 ,colour = "black",margin = unit(c(0.1,0.3,0.1,0.3),"cm"),family = "Arial"),
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
  scale_y_continuous(name = "Molecular richness")+
  #scale_y_continuous(name = "Sum of molecular intensities")+
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
        axis.title.y.left = element_text(size = 20, face=2 ,colour = "black",margin = unit(c(0.1,0.3,0.1,0.3),"cm"),family = "Arial"),
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
mol_inty=as.data.table(aggregate(ft_merge_sel$Bromo.Inty*10000,by=list(Sample=ft_merge_sel$Sample),sum))
mol_inty

mol_inty_n=mol_inty %>% inner_join(FT_envi_sel[,c("Sample","Group","No","Date","PM2.5","OC","POC","SOC")])

mol_inty_n
mol_inty_n=mol_inty_n %>% 
  mutate(date2 = as.POSIXct(mol_inty_n$Date, format = '%Y-%m-%d %H:%M'))


mol_inty_n_m=melt(mol_inty_n, id.vars = c("Sample","Group","No","Date","date2","x"))
#fwrite(mol_inty_n_m,file = "cal_variation.csv")
mol_inty_nna=subset(mol_inty_n_m,!is.na(mol_inty_n_m$value))
mol_inty_min=as.data.table(aggregate(mol_inty_nna$value,
                                     by=list(Group=mol_inty_nna$Group,variable=mol_inty_nna$variable),min)) %>% 
  `colnames<-`(c("Group","variable","value"))
mol_inty_min2=mol_inty_min %>% left_join(mol_inty_nna[,c("Group","variable","value","Sample","x")])

mol_inty_n_m_sel=mol_inty_n_m
mol_inty_n_m_sel$Group=factor(mol_inty_n_m_sel$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))
mol_inty_n_m_sel$Group2=factor(mol_inty_n_m_sel$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                               labels = c(expression(bold("Ulaanbaatar")),
                                          expression(bold("Beijing")),
                                          expression(bold("Seosan")),
                                          expression(bold("Seoul")),
                                          expression(bold("Noto"))
                               ))

mol_inty_n_m_sel
mol_inty_n_m_sel$variable2=factor(mol_inty_n_m_sel$variable,levels = c("PM2.5","OC","POC","SOC"),
                                  labels=c(expression(bold("PM"[2.5]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("OC"~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("POC"~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("SOC"~"("*"\u03bcg/"*m^"3"*")"))
                                  ))





mol_inty_n_ul=subset(mol_inty_n_m_sel,mol_inty_n_m_sel$Group=="Ulaanbaatar")
mol_inty_n_ul$val2=mol_inty_n_ul$x/10000000

ggplot()+
  geom_smooth(data = mol_inty_n_ul ,aes(x=value, y=val2),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_inty_n_ul ,aes(x=value, y=val2, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_inty_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  #scale_y_continuous(name = "Molecular intyness")+
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
  ggsave(filename("molvsOC_inty_ul"), height = 11, width = 50, units = "cm", dpi = 300)


mol_inty_n_bj=subset(mol_inty_n_m_sel,mol_inty_n_m_sel$Group=="Beijing")
mol_inty_n_bj$val2=mol_inty_n_bj$x/10000000

ggplot()+
  geom_smooth(data = mol_inty_n_bj ,aes(x=value, y=val2),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_inty_n_bj ,aes(x=value, y=val2, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_inty_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  #scale_y_continuous(name = "Molecular intyness")+
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
  ggsave(filename("molvsOC_inty_bj"), height = 11, width = 50, units = "cm", dpi = 300)

mol_inty_n_ss=subset(mol_inty_n_m_sel,mol_inty_n_m_sel$Group=="Seosan")
mol_inty_n_ss$val2=mol_inty_n_ss$x/10000000

ggplot()+
  geom_smooth(data = mol_inty_n_ss ,aes(x=value, y=val2),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_inty_n_ss ,aes(x=value, y=val2, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_inty_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  #scale_y_continuous(name = "Molecular intyness")+
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
  ggsave(filename("molvsOC_inty_ss"), height = 11, width = 50, units = "cm", dpi = 300)

mol_inty_n_sul=subset(mol_inty_n_m_sel,mol_inty_n_m_sel$Group=="Seoul")
mol_inty_n_sul$val2=mol_inty_n_sul$x/10000000

ggplot()+
  geom_smooth(data = mol_inty_n_sul ,aes(x=value, y=val2),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_inty_n_sul ,aes(x=value, y=val2, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_inty_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  #scale_y_continuous(name = "Molecular intyness")+
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
  ggsave(filename("molvsOC_inty_sul"), height = 11, width = 50, units = "cm", dpi = 300)

mol_inty_n_nt=subset(mol_inty_n_m_sel,mol_inty_n_m_sel$Group=="Noto")
mol_inty_n_nt$val2=mol_inty_n_nt$x/10000000

ggplot()+
  geom_smooth(data = mol_inty_n_nt ,aes(x=value, y=val2),col="black", lty=2,method = "loess",formula=y~log(x))+
  geom_line(data = mol_inty_n_nt ,aes(x=value, y=val2, col=variable),size=1, na.rm = F)+
  #geom_smooth(data = mol_inty_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  #scale_y_continuous(name = "Molecular intyness")+
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
  ggsave(filename("molvsOC_inty_nt"), height = 11, width = 50, units = "cm", dpi = 300)

##cal_variation=====

mol_rich_n_m_sel=mol_rich_n_m
mol_rich_n_m_sel
mol_rich_min

mol_rich_n_m_sel=mol_rich_n_m_sel[,c("Group","No","Date","date2","x","variable","value")] %>% 
  `colnames<-`(c("Group","No","Date","date2","richness","variable","conc"))

mol_var_min=mol_rich_n_m_sel %>% inner_join(mol_rich_min2)

mol_var_min$difx=mol_var_min$conc-mol_var_min$value
mol_var_min$dify=mol_var_min$richness-mol_var_min$x

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
  #geom_smooth(data = mol_inty_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "lm",formula=y~log(x))+
  #geom_smooth(data = mol_inty_n_m_sel ,aes(x=value, y=x),col="black", lty=2,method = "loess",formula=y~log(x))+
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
  ggsave(filename("inty_variation"), height = 50, width = 50, units = "cm", dpi = 300)

fwrite(mol_var_min,file = "Datafile/molecular_variation.csv")

mol_var_min

####calculate molecular replacement and richness variation=====
ft_merge

#####ul=====
ft_merge_ul=subset(ft_merge,ft_merge$Group=="Ulaanbaatar")
ft_merge_ul

ft_ul=melt(ft_merge_ul[,c("Sample","Formula","Freq")],id.vars = c("Sample","Formula")) %>% 
  dcast(Sample~Formula, sum)

ft_ul_bet=as.data.frame(ft_ul[,-"Sample"]) %>% `row.names<-`(ft_ul$Sample)

ft_ul[,1:23]
beta_ul = beta.div.comp(ft_ul_bet, coef="J", quant=FALSE)

t1=as.data.frame(as.matrix(beta_ul$repl))
t1$Sample=row.names(t1)

repl_ul=t1[,c("Sample","Ulaanbaatar_1")] %>% `colnames<-`(c("Sample","replacement"))

t1=as.data.frame(as.matrix(beta_ul$rich))
t1$Sample=row.names(t1)

rich_ul=t1[,c("Sample","Ulaanbaatar_1")] %>% `colnames<-`(c("Sample","richness"))
var_ul=repl_ul %>% inner_join(rich_ul)
#####bj=====
ft_merge_bj=subset(ft_merge,ft_merge$Group=="Beijing")
ft_merge_bj

ft_bj=melt(ft_merge_bj[,c("Sample","Formula","Freq")],id.vars = c("Sample","Formula")) %>% 
  dcast(Sample~Formula, sum)

ft_bj_bet=as.data.frame(ft_bj[,-"Sample"]) %>% `row.names<-`(ft_bj$Sample)

ft_bj[,1:23]
beta_bj = beta.div.comp(ft_bj_bet, coef="J", quant=FALSE)

t1=as.data.frame(as.matrix(beta_bj$repl))
t1$Sample=row.names(t1)

repl_bj=t1[,c("Sample","Beijing_1")] %>% `colnames<-`(c("Sample","replacement"))

t1=as.data.frame(as.matrix(beta_bj$rich))
t1$Sample=row.names(t1)

rich_bj=t1[,c("Sample","Beijing_1")] %>% `colnames<-`(c("Sample","richness"))
var_bj=repl_bj %>% inner_join(rich_bj)

#####ss=====
ft_merge_ss=subset(ft_merge,ft_merge$Group=="Seosan")
ft_merge_ss

ft_ss=melt(ft_merge_ss[,c("Sample","Formula","Freq")],id.vars = c("Sample","Formula")) %>% 
  dcast(Sample~Formula, sum)

ft_ss_bet=as.data.frame(ft_ss[,-"Sample"]) %>% `row.names<-`(ft_ss$Sample)

ft_ss[,1:23]
beta_ss = beta.div.comp(ft_ss_bet, coef="J", quant=FALSE)

t1=as.data.frame(as.matrix(beta_ss$repl))
t1$Sample=row.names(t1)

repl_ss=t1[,c("Sample","Seosan_1")] %>% `colnames<-`(c("Sample","replacement"))

t1=as.data.frame(as.matrix(beta_ss$rich))
t1$Sample=row.names(t1)

rich_ss=t1[,c("Sample","Seosan_1")] %>% `colnames<-`(c("Sample","richness"))
var_ss=repl_ss %>% inner_join(rich_ss)

#####sul=====
ft_merge_sul=subset(ft_merge,ft_merge$Group=="Seoul")
ft_merge_sul

ft_sul=melt(ft_merge_sul[,c("Sample","Formula","Freq")],id.vars = c("Sample","Formula")) %>% 
  dcast(Sample~Formula, sum)

ft_sul_bet=as.data.frame(ft_sul[,-"Sample"]) %>% `row.names<-`(ft_sul$Sample)

ft_sul[,1:23]
beta_sul = beta.div.comp(ft_sul_bet, coef="J", quant=FALSE)

t1=as.data.frame(as.matrix(beta_sul$repl))
t1$Sample=row.names(t1)

repl_sul=t1[,c("Sample","Seoul_1")] %>% `colnames<-`(c("Sample","replacement"))

t1=as.data.frame(as.matrix(beta_sul$rich))
t1$Sample=row.names(t1)

rich_sul=t1[,c("Sample","Seoul_1")] %>% `colnames<-`(c("Sample","richness"))
var_sul=repl_sul %>% inner_join(rich_sul)

#####nt=====
ft_merge_nt=subset(ft_merge,ft_merge$Group=="Noto")
ft_merge_nt

ft_nt=melt(ft_merge_nt[,c("Sample","Formula","Freq")],id.vars = c("Sample","Formula")) %>% 
  dcast(Sample~Formula, sum)

ft_nt_bet=as.data.frame(ft_nt[,-"Sample"]) %>% `row.names<-`(ft_nt$Sample)

ft_nt[,1:23]
beta_nt = beta.div.comp(ft_nt_bet, coef="J", quant=FALSE)

t1=as.data.frame(as.matrix(beta_nt$repl))
t1$Sample=row.names(t1)

repl_nt=t1[,c("Sample","Noto_1")] %>% `colnames<-`(c("Sample","replacement"))

t1=as.data.frame(as.matrix(beta_nt$rich))
t1$Sample=row.names(t1)

rich_nt=t1[,c("Sample","Noto_1")] %>% `colnames<-`(c("Sample","richness"))
var_nt=repl_nt %>% inner_join(rich_nt)



var_ul
var_bj
var_ss
var_sul
var_nt

aom_var=rbind(var_ul,
              var_bj,
              var_ss,
              var_sul,
              var_nt)

aom_var=aom_var %>% inner_join(ft_envi[,c("Sample","Date")])

aom_var=aom_var %>%  mutate(date2 = as.POSIXct(aom_var$Date, format = '%Y-%m-%d %H:%M'))

aom_var_m=melt(aom_var,id.vars=c("Sample","Date","date2"))

aom_var_m=aom_var_m %>% separate("Sample",c("Group", "No"))

aom_var_m$Group=factor(aom_var_m$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))

aom_var_m$Group2=factor(aom_var_m$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                          labels = c(expression(bold("Ulaanbaatar")),
                                     expression(bold("Beijing")),
                                     expression(bold("Seosan")),
                                     expression(bold("Seoul")),
                                     expression(bold("Noto"))
                          ))

aom_var_m$variable2=factor(aom_var_m$variable, levels = c("replacement","richness"),
                          labels = c("Replacment","Richness difference"))

aom_var_m

ggplot(data=aom_var_m)+
  geom_line(data = aom_var_m ,aes(x=date2, y=value, col=variable2 ),size=0.7, na.rm = F)+
  geom_point(data = aom_var_m ,aes(x=date2, y=value,shape=variable2,col=variable2))+
  facet_wrap(.~Group, ncol = 1, scales = "free")+
  scale_y_continuous(name = "Podani's indice (Jaccard dissimilarity)")+
  #scale_y_continuous(name =expression(bold(PM[2.5]~"("*"\u03bcg/"*m^"3"*")")),expand = c(0.02,0.02),
  #                   sec.axis = sec_axis(~./1.5, name=expression(bold("OC"~"("*"\u03bcg/"*m^"3"*")")))
  #)+
  #scale_fill_manual(values = c("grey70","#FD7F20"))+
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
        axis.title.y.left = element_text(size = 20, colour = "black",face=2,margin = unit(c(0.1,0.6,0.1,0.1),"cm"),family = "Arial"),
        axis.title.y.right = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.0,0.1,0.4),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), size = 14,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.position = c(0.79,1.02))+
  xlab("")+
  guides(col=guide_legend(order = 1,title ="", override.aes = list(shape=c(17,18))),
         shape="none")+
  ggsave(filename("AOM_variation"), height = 30, width = 30, units = "cm", dpi = 300)


