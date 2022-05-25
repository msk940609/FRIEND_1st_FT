library(xlsx)
library(ggplot2)
library(writexl)
library(vegan)
options(java.parameters = "- Xmx8192m")
library(data.table)
library(tidyverse)
library(dplyr)
library(ggsci)
library(extrafont)
loadfonts(device="win")
source("Rscript/func_filename.R")

###file merge======
flist=list.files(pattern = ".xlsx") ###search csv file at working directory

dt=data.table() ###empty table

for (i in 1:length(flist)) {
  
  temp=as.data.table(read.xlsx2(flist[i], sheetIndex = 1))
  temp$Sample=paste(tools::file_path_sans_ext(flist[i]))
  dt=rbind(dt,temp)
} ###merge csv files

dt 

fwrite(dt, file = "Ulan_merge.csv") ###save merge files

###csv merge===
flist=list.files(pattern = ".csv") ###search csv file at working directory

dt=data.table() ###empty table

for (i in 1:length(flist)) {
  
  
  temp=fread(flist[i])
  dt=rbind(dt,temp)
} ###merge csv files

dt

gc_merge=dt

gc_merge=gc_merge %>% separate("Sample",c("Group","Event","Month","Day"), sep="_")
gc_merge$Date=paste(gc_merge$Month,gc_merge$Day,sep = "_")
gc_merge

gc_merge$Sample=paste(gc_merge$Group,gc_merge$Event,gc_merge$Date, sep = "_")

fwrite(gc_merge,file = "gc_2d_merge.csv")

table(gc_merge$Sample)

#data load====
gc_merge=fread("Datafile/gc_2d_merge.csv")

gc_merge$Vol=ifelse(gc_merge$RT1<768,"VOC",
                    ifelse(gc_merge$RT1<2892,"SVOC","NVOC"))
table(gc_merge$Vol)

gc_merge$Pol=ifelse(gc_merge$RT2<4, "LP","HP")
gc_merge$class=paste(gc_merge$Pol,gc_merge$Vol,sep = "_")
table(gc_merge$class)
gc_merge

gc_class_all=as.data.table(aggregate(gc_merge$Norm.Area, by=list("Sample"=gc_merge$Sample,`Group`=gc_merge$Group), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Group","Tot"))
gc_class_all
gc_class=as.data.table(aggregate(gc_merge$Norm.Area, by=list("Sample"=gc_merge$Sample,`Group`=gc_merge$Group,
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

##class distribution====
gc_class_m2=melt(gc_class_m,id.vars = c("Sample"))
gc_class_m2=gc_class_m2 %>% separate("Sample",c("Group","Event","Month","Day"))
gc_class_m2$Group=factor(gc_class_m2$Group, levels = c("Ulaanbaatar","Beijing","Seoul","Seosan"))
gc_class_m2$variable=factor(gc_class_m2$variable, levels = c("HP_VOC","HP_SVOC","HP_NVOC",
                                                             "LP_VOC","LP_SVOC","LP_NVOC"))

ggplot()+
  stat_boxplot(data=gc_class_m2, aes(x=Group, y=value, fill=Group),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=gc_class_m2, aes(x=Group, y=value, fill=Group),alpha=1, outlier.color = NA)+
  geom_boxplot()+
  scale_y_continuous(name = "Proportion of compound classes (%)")+
  facet_wrap(.~variable, scales = "free")+
  scale_fill_manual(values = c("#3C5488FF","#E64B35FF", "#4DBBD5FF","#00A087FF"))+
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
        legend.position = "none")+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("Class_compare_bon"),height = 30, width = 45, units = "cm", dpi = 300)


####
gc_class_d=dcast(gc_class,Group~Class, value.var = "x", fun.aggregate = mean) %>% 
  melt(id.vars=c("Group"))

gc_class_d$Group=factor(gc_class_d$Group, levels = c("Ulaanbaatar","Beijing","Seoul","Seosan"))
gc_class_d$variable=factor(gc_class_d$variable, levels = c("HP_VOC","LP_VOC","HP_SVOC","LP_SVOC",
                                                           "HP_NVOC","LP_NVOC"))

ggplot()+
  geom_bar(data=gc_class_d, aes(x=Group, y=value, fill=variable),stat = "identity", position = position_stack(reverse = T))+
  scale_y_continuous(name = "VOC compound class")+
  scale_fill_manual(values = c( "#AFABAB","#EF4B4B", "#F4B183", "#5B9BD5","#00B050","#7030A0"))+
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
        legend.position = "right")+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("Class_compare_bar"),height = 30, width = 45, units = "cm", dpi = 300)



