library(ggmap)
options(java.parameters = "- Xmx8192m")
library(shiny)
library(plotly)
library(dplyr)
library(data.table)
options(shiny.maxRequestSize=300*1024^2)
getOption("digits")
options("digits" = 15)
library(data.table)
library(RColorBrewer)
options(scipen=1000)
library(rsconnect)
library(tidyverse)
library(bit64)
library(shinyFiles)
library(xlsx)
library(writexl)
library(ggplot2)
library(writexl)
library(xlsx)
library(vegan)
library(data.table)
library(extrafont)
loadfonts(device="win")
library(ggsci)
library(lawstat)
library(agricolae)
library(PMCMR)
library(PMCMRplus)
library(multcompView)
library(DescTools)
library(plyr)
library(grid)
library(cowplot)
library(ggrepel)
library(lmtest)
library(MASS)
library(nparcomp)
library(elevatr)
library(plotly)
library(raster)
library(rgdal)
library(rasterVis)
library(rgl)
library(htmlwidgets)
library(ggtext)
source("Rscript/func_filename.R")

register_google(key = 'AIzaSyC-1iYzvvPOeQd7zCBpmCxxjGI--ebjDxA')


##Slide 1. Project sampling point map======
mapmap=get_map(location = c(lon=122.7782691,lat=40.00161484),zoom=4, maptype = "terrain-background", color = "color")
ggmap(mapmap)

samp=fread("Datafile/Sampleing_point.csv")
ggmap(mapmap)+geom_point(data=samp, aes(x=long, y=lat, col=Location),size=1, alpha=0.9)+
  #geom_text_repel(data=tsub, aes(x=lon, y=lat, label=Sample),size=6,point.padding = unit(0.1,"cm"))+
  #annotate("text",x=96.0, y=40.4,label="G", size=6)+
  scale_x_continuous(name = "Longitude °",breaks = seq(100,140,10),labels =seq(100,140,10),limits = c(105,140), expand = c(0,0))+
  scale_y_continuous(name = "Latitude °", breaks = seq(30,50,5),labels =seq(30,50,5),limits = c(30,50), expand = c(0,0))+
  #scale_color_manual(values = c("#4131EB","#004200","#FF4400","sienna"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(0.2,1.0,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial", margin =unit(c(0.2,0.2,0.2,0.2),"cm") ),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin =unit(c(0.2,0.2,0.2,0.2),"cm")),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        legend.text = element_text(size = 20, colour = "black", family = "Arial",vjust = 0.6,margin = unit(c(0.0,-0.2,0.0,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.box.margin = margin(t=0,r=0,b=0,l=0, unit = "cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "NULL",
        legend.box.just = "left",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_blank()
  )+ ggsave(filename("map_all"), height = 20, width = 25, units = "cm", dpi = 700)


neAsia <- map_data(map = 'world',
                   region = c('South Korea', 'North Korea', 'China',
                              'Japan', 'Mongolia', 'Taiwan',"Russia"))

neAsia$PM=ifelse(neAsia$region=="Mongolia",93.082,
                 ifelse(neAsia$region=="China",30.041,
                        ifelse(neAsia$region=="South Korea",25.994,
                               ifelse(neAsia$region=="Japan",9.678, 28))))



ggplot() +
  geom_polygon(data = neAsia,mapping = aes(x = long, y = lat, group = group),fill = "White", color = 'black')+
  geom_point(data=samp, aes(x=long, y=lat, col=Location),size=1, alpha=0.9)+
  scale_x_continuous(name = "Longitude °",breaks = seq(100,140,10),labels =seq(100,140,10))+
  scale_y_continuous(name = "Latitude °", breaks = seq(30,50,5),labels =seq(30,50,5))+
  coord_map(xlim = c(105,140),ylim = c(30,48))+
  theme_bw()+
  #scale_fill_gradient2(low = "yellow",mid = "#FEB24C", midpoint = 0, high = "#F03B20")+
  #scale_fill_distiller(palette = "Spectral")+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "lightblue"),
        #plot.background = element_rect()
        plot.margin = unit(c(0.0,1.0,0.0,0.2),"cm"),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial", margin =unit(c(0.2,0.2,0.2,0.2),"cm") ),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin =unit(c(0.2,0.2,0.2,0.2),"cm")),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        legend.text = element_text(size = 20, colour = "black", family = "Arial",vjust = 0.6,margin = unit(c(0.0,-0.2,0.0,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.box.margin = margin(t=0,r=0,b=0,l=0, unit = "cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "NULL",
        legend.box.just = "left",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_blank()
  )+ ggsave(filename("map_all2"), height = 20, width = 25, units = "cm", dpi = 700)


ggplot() +
  geom_polygon(data = neAsia,mapping = aes(x = long, y = lat, group = group,fill = PM), color = 'black')+
  geom_point(data=samp, aes(x=long, y=lat, col=Location),size=1, alpha=0.9)+
  scale_x_continuous(name = "Longitude °",breaks = seq(100,140,10),labels =seq(100,140,10))+
  scale_y_continuous(name = "Latitude °", breaks = seq(30,50,5),labels =seq(30,50,5))+
  coord_map(xlim = c(105,140),ylim = c(30,48))+
  theme_bw()+
  scale_fill_gradient2(low = "yellow",mid = "#FEB24C", midpoint = 0, high = "#F03B20")+
  #scale_fill_distiller(palette = "Spectral")+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "lightblue"),
        #plot.background = element_rect()
        plot.margin = unit(c(0.0,1.0,0.0,0.2),"cm"),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial", margin =unit(c(0.2,0.2,0.2,0.2),"cm") ),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin =unit(c(0.2,0.2,0.2,0.2),"cm")),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        legend.text = element_text(size = 20, colour = "black", family = "Arial",vjust = 0.6,margin = unit(c(0.0,-0.2,0.0,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.box.margin = margin(t=0,r=0,b=0,l=0, unit = "cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "right",
        legend.box.just = "left",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_blank()
  )+ ggsave(filename("map_all3"), height = 20, width = 25, units = "cm", dpi = 700)

#Slide2-4 2D GC data====

##2-1 2D GC RT plot====.
gc2d_raw=fread("Datafile/gc_2d_merge.csv")
gc2d_raw
gc2d_raw$id=paste(gc2d_raw$RT1,gc2d_raw$RT2,sep = "_")

table(gc2d_raw$Sample)
gc2d_raw_sel=subset(gc2d_raw,gc2d_raw$Sample=="Seoul_Event_12_23")

gc2d_m_sel=melt(gc2d_raw_sel[,c("Group","Date","id","RT1","RT2","Area")], id.vars = c("Group","Date","id","RT1","RT2"))

timeHM_formatter <- function(x) {
  h <- floor(x/60)
  m <- floor(x %% 60)
  lab <- sprintf("%d:%02d", h, m) # Format the strings as HH:MM
  return(lab)
}
matlab.like(20)
?scale_size_area()
matlab.like(20)
tt=matlab.like(20)
ygobb(20)
gc2d_m_sel=gc2d_m_sel[order(gc2d_m_sel$value),]


ggplot()+
  geom_rect(aes(xmin=-Inf,xmax= 768, ymax=-Inf, ymin=4),fill="#EF423F",col="#EF423F",alpha=0.6)+
  geom_rect(aes(xmin=-Inf,xmax= 768, ymax=4, ymin=Inf),fill="orange1",col="orange1",alpha=0.6)+
  geom_rect(aes(xmin=768,xmax= 2892, ymax=Inf, ymin=4),fill="#00964C",col="#00964C",alpha=0.6)+
  geom_rect(aes(xmin=768,xmax= 2892, ymax=4, ymin=-Inf),fill="#FAE32C",col="#FAE32C",alpha=0.6)+
  geom_rect(aes(xmin=2892,xmax= Inf, ymax=Inf, ymin=4),fill="#503177",col="#503177",alpha=0.6)+
  geom_rect(aes(xmin=2892,xmax= Inf, ymax=4, ymin=-Inf),fill="steelblue1",col="steelblue1",alpha=0.6)+
  #geom_point(data=gc2d_m_sel,aes(x=RT1,y=RT2, size=value),fill="grey55",col="grey95", shape=21,alpha=0.8 )+
  geom_point(data=gc2d_m_sel,aes(x=RT1,y=RT2, size=value,fill=value),col="grey95", shape=21,alpha=0.8 )+
  #geom_point(data=gc2d_m_sel,aes(x=RT1,y=RT2, size=value,col=value),alpha=0.8 )+
  scale_fill_gradientn(colors=c("grey50",matlab.like(40)[6:25],matlab.like(40)[30:30]))+
  #scale_fill_gradientn(colors=ygobb(20)[5:10])+
  #scale_fill_gradientn(colors=ygobb(30)[12:20])+
  #scale_color_gradientn(colors=(cm.colors(40)[4:40]))+
  #scale_color_gradient(high="orange",low = "purple")+
  geom_hline(yintercept = 4, lty=2)+
  geom_vline(xintercept = 768,lty=2)+
  geom_vline(xintercept = 2892,lty=2)+
  scale_size(range = c(2, 15),
             breaks = 1000000 * c(2.5, 5.0, 7.50, 10.00),
             labels = c(expression(bold("2.5 x 10"^"6")), expression(bold("5.0 x 10"^"6")), expression(bold("7.5 x 10"^"6")),
                        expression(bold("10.0 x 10"^"6"))))+
  scale_x_continuous(name=expression(bold("1"^"st"~"RT (min:s)")),labels=timeHM_formatter)+
  scale_y_continuous(name=expression(bold("2"^"nd"~"RT (s)")),breaks=seq(0,8,2),labels=seq(0,8,2),limits = c(0,8.5))+
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
        axis.title.x = element_text(size = 20, colour = "black",margin = unit(c(0.2,0.1,0.1,0.1),"cm"),family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black" , face = "bold",family = "Arial"),
        axis.title.y.left = element_text(size = 20, colour = "black",margin = unit(c(0.1,0.2,0.1,0.1),"cm"),family = "Arial"),
        legend.text = element_text(size = 14, colour = "black",family = "Arial",margin = unit(c(0.1,0.1,0.2,0.1),"cm"), hjust = 0.0),
        legend.title = element_text(margin = unit(c(0.0,0.1,0.0,0.1),"cm"), face="bold",size = 18,family = "Arial"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(0.5,"cm"),
        legend.direction = "vertical",
        legend.background = element_blank(),
        #legend.position = c(0.08,1.02)
        legend.position = "right")+
  guides(size=guide_legend(title = "Intensity",title.hjust = 0.5, override.aes = list(shape=21,size=7, fill=c("grey55","#0055FF","#6DFFED","#FFAA00"),alpha=1)),
         fill="none")+
  ggsave(filename("2DGC_vocclass"), height = 20, width = 50, units = "cm", dpi = 300)


matlab.like(40)[6:25]
matlab.like(40)[30:30]

matlab.like(40)[25:40]
gc_merge=fread("Datafile/gc_2d_merge.csv")

gc_merge$Vol=ifelse(gc_merge$RT1<768,"VOC",
                    ifelse(gc_merge$RT1<2892,"SVOC","NVOC"))
table(gc_merge$Vol)

gc_merge$Pol=ifelse(gc_merge$RT2<4, "LP","HP")
gc_merge$class=paste(gc_merge$Pol,gc_merge$Vol,sep = "_")
table(gc_merge$class)
gc_merge

gc_class_all=as.data.table(aggregate(gc_merge$Area, by=list("Sample"=gc_merge$Sample,`Group`=gc_merge$Group), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Group","Tot"))
gc_class_all
gc_class=as.data.table(aggregate(gc_merge$Area, by=list("Sample"=gc_merge$Sample,`Group`=gc_merge$Group,
                                                             "Class"=gc_merge$class), FUN=sum))
gc_class
gc_class=gc_class %>% inner_join(gc_class_all)
gc_class$rel=round(gc_class$x/gc_class$Tot*100,1)

gc_class_m=dcast(gc_class,Sample~Class, value.var = "rel", fun.aggregate = sum)
gc_class_m


gc_merge$id=paste(gc_merge$RT1,gc_merge$RT2,sep = "_")

#####MDS chemical composition=====
MDS_2dgc=melt(gc_merge[,c("Sample","id","Norm.Area")], id.vars = c("Sample","id")) %>% 
  dcast(Sample~id, sum)

MDS_2dgc[,1:2]

NMDS_2dgc=metaMDS(MDS_2dgc[,-c(1)], k=2, distance = "bray", trymax = 20)
NMDS_2dgc ##Stress 0.177

gnmds_2dgc=as.data.table(NMDS_2dgc$points)
gnmds_2dgc$Sample=MDS_2dgc$Sample

gnmds_2dgc=gnmds_2dgc %>% tidyr::separate("Sample",c("Group","Event","Month","Day"),sep="_")
gnmds_2dgc

ggplot()+
  geom_point(data = gnmds_2dgc, aes(x=MDS1, y=MDS2, col=Group))


gnmds_2dgc$Group=factor(gnmds_2dgc$Group, levels = c("Seoul","Seosan","Beijing","Ulaanbaatar"))
gnmds_2dgc$Sample=paste(gnmds_2dgc$Group,gnmds_2dgc$Event,gnmds_2dgc$Month,gnmds_2dgc$Day,sep = "_")


ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_2dgc, aes(x=MDS1, y=MDS2, fill=Group, shape=Group), col="black", size=9)+
  scale_fill_npg()+
  #scale_fill_manual(values = c("#e97f02","#0080ff"))+
  scale_shape_manual(values = c(21,22,23,25))+
  #geom_segment(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=0, y=0, xend=MDS1*0.9, yend=MDS2*0.9),
  #             arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=MDS1*0.98, y=MDS2*0.98, label=variable), size=8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,-0.2,0.1,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.71, 0.14),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5,byrow = T , ncol = 2,
                           override.aes = list(size=7, shape=c(21,22,23,25))),
         shape="none")+
  ggsave(filename("MDS_2dgc"), height = 20, width = 20, units = "cm", dpi = 700)


vec_mds=gnmds_2dgc[,c(1,2)]
vec_mds

gc_class_m

vec_class=gnmds_2dgc %>% inner_join(gc_class_m, by=c("Sample"))

vec_class=vec_class[,c("HP_NVOC","HP_SVOC","HP_VOC","LP_NVOC","LP_SVOC","LP_VOC")]

vec_fit_class=envfit(vec_mds,vec_class, na.rm = T)
vec_fit_class$vectors$arrows
vec_fit_class$vectors$r
vec_fit_class$vectors$pval

arrow_4th_class=as.data.frame(scores(vec_fit_class, display = "vectors"))
arrow_4th_class$variable=row.names(arrow_4th_class)
arrow_4th_class$R=vec_fit_class$vectors$r
arrow_4th_class$p=vec_fit_class$vectors$pvals

arrow_4th_class
arrow_4th_class$variable=row.names(arrow_4th_class)

gnmds_1st

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_2dgc, aes(x=MDS1, y=MDS2, fill=Group), col="black", size=9, shape=21)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#3C5488FF"))+
  scale_x_continuous(name = "NMDS1", breaks =seq(-2,2,0.2), limits = c(-0.6,0.6))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-1.0,1.0,0.2), limits = c(-0.4,0.4))+
  #scale_shape_manual(values = c(24,21,22))+
  #scale_shape_manual(values = c(21,22,23,24,25))+
  geom_segment(data = subset(arrow_4th_class,arrow_4th_class$p<0.05), aes(x=0, y=0, xend=MDS1*0.50, yend=MDS2*0.50),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_class,arrow_4th_class$p<0.05), aes(x=MDS1*0.55, y=MDS2*0.55, label=variable), size=8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial", margin =unit(c(0.2,0.2,0.2,0.2),"cm") ),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        legend.text = element_text(size = 20, colour = "black", family = "Arial",vjust = 0.6,margin = unit(c(0.0,-0.2,0.0,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.box.margin = margin(t=0,r=0,b=0,l=-2.5, unit = "cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "bottom",
        legend.box.just = "left",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_blank()
  )+
  guides(fill=guide_legend(title = NULL,title.hjust = 0,title.position = "top",byrow = T , ncol = 5,
                           override.aes = list(size=7, shape=21), order = 1),
         shape=guide_legend(title = NULL,title.position = "top",override.aes = list(size=8, fill="black"), ncol=3))+
  ggsave(filename("MDS_1st_withclass_nolabel"), height = 24, width = 20, units = "cm", dpi = 700)

####
gc_class_d=dcast(gc_class,Group~Class, value.var = "rel", fun.aggregate = mean) %>% 
  melt(id.vars=c("Group"))

gc_class_d$Group=factor(gc_class_d$Group, levels = c("Ulaanbaatar","Beijing","Seoul","Seosan"))
gc_class_d$variable=factor(gc_class_d$variable, levels = c("HP_SVOC","LP_SVOC","HP_NVOC","HP_VOC","LP_NVOC",
                                                           "LP_VOC"),
                           labels = c("HP-SVOC","LP-SVOC","HP-NVOC","HP-VOC","LP-NVOC",
                                      "LP-VOC"))

ggplot()+
  geom_bar(data=gc_class_d, aes(x=Group, y=value, fill=variable),stat = "identity", position = position_stack(reverse = T))+
  scale_y_continuous(name = "Proportion of VOC class (%)", expand = c(0.008,0.008))+
  scale_fill_manual(values = c("#00964C","#FAE32C","#503177","#EF423F","steelblue1","orange1"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 16,colour = "black", angle = 0, hjust = 0.5,vjust = 1.0),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15,"cm"),
        axis.text.y = element_text(size = 16, colour = "black" ),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        axis.title.y.left = element_text(size = 18, colour = "black",margin = unit(c(0.2,0.5,0.2,0.2),"cm")),
        axis.title.y.right = element_text(size = 18, colour = "black",angle = 90,margin = unit(c(0.2,0.2,0.2,0.5),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.4,0.2,0.2,0.2),"cm")),
        strip.text.y = element_blank(),
        #legend.margin = unit(c(0.0,0.0,0.0,0.0),"cm"),
        legend.title = element_text(size = 22, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 18, colour = "black",margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        legend.key.width = unit(1.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.background = element_blank(),
        legend.direction = "vertical",
        legend.justification=c(0, 1),
        legend.position = "NULL")+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5, reverse = T))+
  ggsave(filename("Class_compare_bar_per"),height = 20, width = 20, units = "cm", dpi = 300)

gc_class_d2=dcast(gc_class,Group~Class, value.var = "x", fun.aggregate = mean) %>% 
  melt(id.vars=c("Group"))

gc_class_d2$Group=factor(gc_class_d2$Group, levels = c("Ulaanbaatar","Beijing","Seoul","Seosan"))
gc_class_d2$variable=factor(gc_class_d2$variable, levels = c("HP_SVOC","LP_SVOC","HP_NVOC","HP_VOC","LP_NVOC",
                                                           "LP_VOC"),
                           labels = c("HP-SVOC","LP-SVOC","HP-NVOC","HP-VOC","LP-NVOC",
                                      "LP-VOC"))

sci <- function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))}
  
ggplot()+
  geom_bar(data=gc_class_d2, aes(x=Group, y=value, fill=variable),stat = "identity", position = position_stack(reverse = T))+
  scale_y_continuous(name = "Distribution of VOC class ", expand = c(0.008,0.008),label = sci)+
  scale_fill_manual(values = c("#00964C","#FAE32C","#503177","#EF423F","steelblue1","orange1"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 16,colour = "black", angle = 0, hjust = 0.5,vjust = 1.0),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15,"cm"),
        axis.text.y = element_text(size = 16, colour = "black" ),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        axis.title.y.left = element_text(size = 18, colour = "black",margin = unit(c(0.2,0.5,0.2,0.2),"cm")),
        axis.title.y.right = element_text(size = 18, colour = "black",angle = 90,margin = unit(c(0.2,0.2,0.2,0.5),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.4,0.2,0.2,0.2),"cm")),
        strip.text.y = element_blank(),
        #legend.margin = unit(c(0.0,0.0,0.0,0.0),"cm"),
        legend.title = element_text(size = 22, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 14, colour = "black",margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        legend.key.width = unit(1.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.background = element_blank(),
        legend.direction = "vertical",
        legend.justification=c(0, 1),
        legend.position = c(0.77,1.02))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5,reverse=T))+
  ggsave(filename("Class_compare_bar_val"),height = 20, width = 22, units = "cm", dpi = 300)

##2-2 GC (PAH, N,Oxy PAH) NMDS vs 2D GC NMDS===== (valiadation with stress, adonis)
###2-2-1NMDS with PAH====
#pah=fread("Datafile/PAH_OCnorm_per_ewha.csv")
pah=fread("Datafile/PAH_OCnorm_ewha.csv")
#pah=fread("Datafile/PAH_merge_ewha.csv")
pah=pah[,c(1:42)]
pah_q_m=melt(pah[,-c("Sample Name_number","Sample Name_date")], id.vars = c("Group","No"))
pah_q_m$value=ifelse(pah_q_m$value<0,0,pah_q_m$value)

pah_q_m=pah_q_m %>% 
  mutate("Sample"=paste(Group,No,sep = "_")) %>% 
  dcast(Sample+Group+No~variable,sum)
pah_q_m

pah_q_sel=subset(pah_q_m,pah_q_m$Group!="Noto") %>% droplevels()


nmds_ewha=metaMDS(pah_q_sel[,-c(1:3)])
nmds_ewha##stress0.11

gnmds_ewha=as.data.table(nmds_ewha$points)
gnmds_ewha$Sample=pah_q_sel$Sample

gnmds_ewha=gnmds_ewha %>% tidyr::separate("Sample",c("Group","No"),sep="_")
gnmds_ewha

ggplot()+
  geom_point(data = gnmds_ewha, aes(x=MDS1, y=MDS2, col=Group))


gnmds_ewha$Group=factor(gnmds_ewha$Group, levels = c("Seoul","Seosan","Beijing","Ulaanbaatar"))
gnmds_ewha

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_ewha, aes(x=MDS1, y=MDS2, fill=Group),shape=21, col="black", size=7)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#3C5488FF"))+
  scale_x_continuous(name = "NMDS1", breaks =seq(-1.5,1.5,0.5), limits = c(-1.2,1.2))+
  #scale_y_continuous(name = "NMDS2", breaks =seq(-0.30,0.6,0.3), limits = c(-0.40,0.70))+
  #scale_shape_manual(values = c(24,21,22))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial", margin =unit(c(0.2,0.2,0.2,0.2),"cm") ),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        legend.text = element_text(size = 20, colour = "black", family = "Arial",vjust = 0.6,margin = unit(c(0.0,-0.2,0.0,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.box.margin = margin(t=0,r=0,b=0,l=-2.5, unit = "cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "bottom",
        legend.box.just = "left",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_blank()
  )+
  guides(fill=guide_legend(title = NULL,title.hjust = 0,title.position = "top",byrow = T , ncol = 5,
                           override.aes = list(size=7, shape=21), order = 1),
         shape=guide_legend(title = NULL,title.position = "top",override.aes = list(size=8, fill="black"), ncol=3))+
  ggsave(filename("MDS_ewhaGC"), height = 22, width = 20, units = "cm", dpi = 700)

adonis(pah_q_sel[,-c(1:3)]~Group,data=pah_q_sel) 
##R2:0.61***



###2-2-2 NMDS with 2D GC_PAH ====
pah_merge=fread("Datafile/pah_merge.csv")

pah_m=melt(pah_merge[,c("Group","Event","Formula","Normarea")], id.vars = c("Group","Event","Formula")) %>% 
  dcast(Group+Event~Formula, mean) %>% 
  melt(id.vars=c("Group","Event"), na.rm = T) %>% `colnames<-`(c("Group","Event","Formula","Norm Area"))

pah_m=melt(pah_merge[,c("Sample","Group","Name","Normarea")], id.vars = c("Sample","Group","Name")) %>% 
  dcast(Sample+Group~Name, value.var = "value",sum)

pah_m[,1:12]
dim(pah_m)

nmds_2d=metaMDS(pah_m[,-c(1,2)])
nmds_2d##stress"0.10

gnmds_2d=as.data.table(nmds_2d$points)
gnmds_2d$Sample=pah_m$Sample

gnmds_2d=gnmds_2d %>% tidyr::separate("Sample",c("Group","Date","Event"),sep="_")
gnmds_2d

ggplot()+
  geom_point(data = gnmds_2d, aes(x=MDS1, y=MDS2, col=Group))


gnmds_2d$Group=factor(gnmds_2d$Group, levels = c("Seoul","Seosan","Beijing","Ulaanbaatar"))
gnmds_2d

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_2d, aes(x=MDS1, y=MDS2, fill=Group),shape=21, col="black", size=9)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#3C5488FF"))+
  scale_x_continuous(name = "NMDS1", breaks =seq(-1.0,2.0,0.5), limits = c(-1.2,2.1))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-1.0,1.0,0.5), limits = c(-0.8,0.65))+
  #scale_shape_manual(values = c(24,21,22))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial", margin =unit(c(0.2,0.2,0.2,0.2),"cm") ),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        legend.text = element_text(size = 20, colour = "black", family = "Arial",vjust = 0.6,margin = unit(c(0.0,-0.2,0.0,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.box.margin = margin(t=0,r=0,b=0,l=-2.5, unit = "cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "bottom",
        legend.box.just = "left",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_blank()
  )+
  guides(fill=guide_legend(title = NULL,title.hjust = 0,title.position = "top",byrow = T , ncol = 5,
                           override.aes = list(size=7, shape=21), order = 1),
         shape=guide_legend(title = NULL,title.position = "top",override.aes = list(size=8, fill="black"), ncol=3))+
  ggsave(filename("MDS_2DGC"), height = 22, width = 20, units = "cm", dpi = 700)

adonis(pah_m[,-c(1,2)]~Group,data=pah_m) 
##R2:0.50***

##2-3 GC VOC class 

#Slide 5-8 FT ICR=====

##5-1 NMDS with comp* chp vector (all region)
##5-2 NMDS with comp* chp vector (Seoul, Seosan, Beijing)

##6-1 Overall comparision among Regions (Comp, Chp)=====



##7-1 Event vs Non-event comparision among Regions (Comp, Chp)====



##8-1 Correaltion heatmap between Comp + molecular class and Envi (Beijing, Mongol and Seoul) 
##9-1 Correaltion test between Comp + molecular class and Envi (Seoul)
##10-1 Correaltion test between Comp + molecular class and Envi (Beijing)
##11-1 Correaltion test between Comp + molecular class and Envi (Mongol)

# No slide====
##Pie chart averaged chemical composition according sampling point =====

#Slide 9~ SEM=====



