library(ggmap)
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
library(ggmap)
register_google(key = 'AIzaSyC-1iYzvvPOeQd7zCBpmCxxjGI--ebjDxA')
library(maps)
library(mapdata)
library(mapproj) 
library(raster)

samp=fread("Datafile/sampling_point.csv")

samp$Location=factor(samp$Location, levels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))

##Slide 1. Project sampling point map======
mapmap=get_map(location = c(lon=122.7782691,lat=40.00161484),zoom=4, maptype = "terrain-background", color = "color")
ggmap(mapmap)


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


neAsia <- map_data(map = 'world',                   region = c('South Korea', 'North Korea', 'China',
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
