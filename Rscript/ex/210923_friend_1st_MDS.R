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

source("Rscript/func_filename.R")
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")

FT_b=fread("Datafile/FRIEND_B_merge.csv")
FT_NT=fread("Datafile/FRIEND_NT_merge.csv")
FT_ULb=fread("Datafile/FRIEND_ULb_merge.csv")
FT_SS=fread("Datafile/FRIEND_SS_merge.csv")
FT_SUL=fread("Datafile/FRIEND_SUL_merge.csv")

FT_b=FT_b %>% tidyr::separate(Sample, c("Group", "No"),sep="_", extra="drop")
FT_b$Sample=paste(FT_b$Group, FT_b$No, sep = "_")
#fwrite(FT_b, "Datafile/FRIEND_B_merge.csv")


FT_NT=FT_NT %>% tidyr::separate(Sample, c("Group", "No"),sep="_", extra="drop")
FT_NT$Sample=paste(FT_NT$Group, FT_NT$No, sep = "_")
#fwrite(FT_NT, "Datafile/FRIEND_NT_merge.csv")

FT_ULb=FT_ULb %>% tidyr::separate(Sample, c("Group", "No"),sep="_", extra="drop")
FT_ULb$No=as.numeric(FT_ULb$No)
FT_ULb$Sample=paste(FT_ULb$Group, FT_ULb$No, sep = "_")
#fwrite(FT_ULb, "Datafile/FRIEND_ULb_merge.csv")

FT_SS=FT_SS %>% tidyr::separate(Sample, c("Group", "No"),sep="_", extra="drop")
FT_SS$Sample=paste(FT_SS$Group, FT_SS$No, sep = "_")
#fwrite(FT_SS, "Datafile/FRIEND_SS_merge.csv")

FT_SUL=FT_SUL %>% tidyr::separate(Sample, c("Group", "No"),sep="_", extra="drop")
FT_SUL$Sample=paste(FT_SUL$Group, FT_SUL$No, sep = "_")
#fwrite(FT_SUL, "Datafile/FRIEND_SUL_merge.csv")

FT_merge=rbind(FT_SUL, FT_b, FT_SS, FT_NT, FT_ULb)

FT_merge$Sample=paste(FT_merge$Group, FT_merge$No, sep = "_")
FT_merge

#fwrite(FT_merge, "FRIENDs_1st_FT.csv")
#0.Load data=====
FT_merge=fread("Datafile/FRIENDs_1st_FT.csv")

FT_merge$C=ifelse(FT_merge$C.>0, "C","")
FT_merge$H=ifelse(FT_merge$H.>0, "H","")
FT_merge$O=ifelse(FT_merge$O.>0, "O","")
FT_merge$N=ifelse(FT_merge$N.>0, "N","")
FT_merge$S=ifelse(FT_merge$S.>0, "S","")

FT_merge=FT_merge %>% unite("Comp",c("C","H","O","N","S"),sep = "")
FT_merge$Comp=ifelse(FT_merge$O.==0, "Remainders",FT_merge$Comp)

FT_merge=subset(FT_merge,FT_merge$Comp!="Remainders") %>% droplevels()

##1.Multi dimension====
###1-1PCA=======
PC_1st=prcomp(MDS_1st_in[,-c(1)], center = T, scale. = F)
PC_1st 

gpca_1st=as.data.table(PC_1st$x)
gpca_1st$Sample=MDS_1st_in$Sample

gpca_1st=gpca_1st %>% tidyr::separate("Sample",c("Group","no"),sep="_")
gpca_1st

ggplot()+
  geom_point(data = gpca_1st, aes(x=PC1, y=PC2, col=Group))

###1-2NMDS======
MDS_1st=melt(FT_merge[,c("Sample","Formula","Bromo.Inty")], id.vars = c("Sample","Formula")) %>% 
  dcast(Sample~Formula, sum)

MDS_1st[,1:2]
table(MDS_1st_in$Sample)
MDS_1st_in=MDS_1st

NMDS_1st=metaMDS(MDS_1st_in[,-c(1)], k=5, distance = "bray", trymax = 20)
NMDS_1st ##Stress 0.08

gnmds_1st=as.data.table(NMDS_1st$points)
gnmds_1st$Sample=MDS_1st_in$Sample

gnmds_1st=gnmds_1st %>% tidyr::separate("Sample",c("Group","no"),sep="_")
gnmds_1st

ggplot()+
  geom_point(data = gnmds_1st, aes(x=MDS1, y=MDS2, col=Group))


gnmds_1st$Group=factor(gnmds_1st$Group, levels = c("SUL","SS","B","NT","UL"),
                       labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
gnmds_1st$Sample=paste(gnmds_1st$Group,gnmds_1st$no, sep = "_")

gnmds_1st$id=ifelse(gnmds_1st$Group=="Beijing",5,
                        ifelse(gnmds_1st$Group=="Seosan",4,
                               ifelse(gnmds_1st$Group=="seoul",3,2)))
gnmds_1st=gnmds_1st[order(gnmds_1st$id)]

gnmds_1st

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st, aes(x=MDS1, y=MDS2, fill=Group, shape=Group), col="black", size=9)+
  scale_fill_npg()+
  #scale_fill_manual(values = c("#e97f02","#0080ff"))+
  scale_shape_manual(values = c(21,22,23,24,25))+
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
                           override.aes = list(size=7, shape=c(21,22,23,24,25))),
         shape="none")+
  ggsave(filename("MDS_1st"), height = 20, width = 20, units = "cm", dpi = 700)


####1-2+ NMDS with vector=====
vec_mds=gnmds_1st[,c(1,2)]
vec_mds

pal=pal_npg("nrc")(5)

#####1-2-1 Enviromental data=====
gnmds_1st

envi_1st=fread("Datafile/envi_1st_sel.csv")
envi_1st
table(envi_1st$Group)

gnmds_1st=gnmds_1st %>% inner_join(envi_1st)

#gnmds_1st=gnmds_1st[order(gnmds_1st$id)]
gnmds_1st
table(gnmds_1st$Group)
gnmds_1st$Event=factor(gnmds_1st$Event,levels = c("Event","Normal","Non-event"))

vec_envi=gnmds_1st[,c("CO","SO2","NO","PM2.5")]
vec_envi

vec_fit_envi=envfit(vec_mds,vec_envi, na.rm = T)
vec_fit_envi$vectors$arrows
vec_fit_envi$vectors$r
vec_fit_envi$vectors$pval

arrow_4th_envi=as.data.frame(scores(vec_fit_envi, display = "vectors"))
arrow_4th_envi$variable=row.names(arrow_4th_envi)
arrow_4th_envi$R=vec_fit_envi$vectors$r
arrow_4th_envi$p=vec_fit_envi$vectors$pvals

arrow_4th_envi
arrow_4th_envi$variable=row.names(arrow_4th_envi)


gnmds_1st$Group=factor(gnmds_1st$Group, levels =  c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
gnmds_1st$Event=factor(gnmds_1st$Event,levels = c("Event","Normal","Non-event"))

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
  #scale_fill_npg()+
  scale_x_continuous(name = "NMDS1", breaks =seq(-2,2,0.5), limits = c(-1.0,2))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-1.5,1.5,0.5), limits = c(-1.2,1.2))+
  scale_shape_manual(values = c(24,21,22))+
  #scale_shape_manual(values = c(21,22,23,24,25))+
  geom_segment(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=0, y=0, xend=MDS1*1, yend=MDS2*1),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=MDS1*0.88, y=MDS2*0.88, label=variable), size=8)+
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
  ggsave(filename("MDS_1st_withenvi_nolabel2"), height = 24, width = 20, units = "cm", dpi = 700)


#####1-2-2 chemical composition=====
vec_mds=gnmds_1st[,c(1,2)]
vec_mds


chcomp=fread("Datafile/Chemicalcomposition_1st")
chcomp
chcomp$Sample=paste(chcomp$Group,chcomp$no,sep = "_")

vec_comp=dcast(chcomp, Sample+Group+no~Comp, sum, value.var = "rel")

vec_comp=gnmds_1st %>% inner_join(vec_comp[,-c("Group","no")], by=c("Sample"))

vec_comp=vec_comp[,c("CHO","CHON","CHOS","CHONS")]

vec_fit_comp=envfit(vec_mds,vec_comp, na.rm = T)
vec_fit_comp$vectors$arrows
vec_fit_comp$vectors$r
vec_fit_comp$vectors$pval

arrow_4th_comp=as.data.frame(scores(vec_fit_comp, display = "vectors"))
arrow_4th_comp$variable=row.names(arrow_4th_comp)
arrow_4th_comp$R=vec_fit_comp$vectors$r
arrow_4th_comp$p=vec_fit_comp$vectors$pvals

arrow_4th_comp
arrow_4th_comp$variable=row.names(arrow_4th_comp)

gnmds_1st

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
  scale_x_continuous(name = "NMDS1", breaks =seq(-2,2,0.5), limits = c(-1.0,2))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-1.5,1.5,0.5), limits = c(-1.2,1.2))+
  scale_shape_manual(values = c(24,21,22))+
  #scale_shape_manual(values = c(21,22,23,24,25))+
  geom_segment(data = subset(arrow_4th_comp,arrow_4th_comp$p<0.05), aes(x=0, y=0, xend=MDS1, yend=MDS2),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=MDS1*0.88, y=MDS2*0.88, label=variable), size=8)+
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
  ggsave(filename("MDS_1st_withComp_nolabel"), height = 24, width = 20, units = "cm", dpi = 700)


#####1-2-3 chemical composition=====
vec_mds=gnmds_1st[,c(1,2)]
vec_mds

chp_m=fread("Datafile/Chemicalprop_1st.csv")
chp_m
chp_m$Sample=paste(chp_m$Group,chp_m$no,sep = "_")

vec_chp=dcast(chp_m, Sample+Group+no~variable, sum, value.var = "val")
vec_chp

vec_chp=gnmds_1st %>% inner_join(vec_chp[,-c("Group","no")], by=c("Sample"))

vec_chp=vec_chp[,c("Mean AI","Mean DBE","Mean H/C","Mean N/C","Mean O/C","Mean S/C")]

vec_fit_chp=envfit(vec_mds,vec_chp, na.rm = T)
vec_fit_chp$vectors$arrows
vec_fit_chp$vectors$r
vec_fit_chp$vectors$pval

arrow_4th_chp=as.data.frame(scores(vec_fit_chp, display = "vectors"))
arrow_4th_chp$variable=row.names(arrow_4th_chp)
arrow_4th_chp$R=vec_fit_chp$vectors$r
arrow_4th_chp$p=vec_fit_chp$vectors$pvals

arrow_4th_chp
arrow_4th_chp$variable=row.names(arrow_4th_chp)

gnmds_1st

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
  scale_x_continuous(name = "NMDS1", breaks =seq(-2,2,0.5), limits = c(-1.0,2))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-1.5,1.5,0.5), limits = c(-1.2,1.2))+
  scale_shape_manual(values = c(24,21,22))+
  #scale_shape_manual(values = c(21,22,23,24,25))+
  geom_segment(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=0, y=0, xend=MDS1, yend=MDS2),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=MDS1*0.88, y=MDS2*0.88, label=variable), size=8)+
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
  ggsave(filename("MDS_1st_withChp_nolabel"), height = 24, width = 20, units = "cm", dpi = 700)


ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
  scale_x_continuous(name = "NMDS1", breaks =seq(-2,2,0.5), limits = c(-1.4,2))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-1.5,1.5,0.5), limits = c(-1.2,1.3))+
  scale_shape_manual(values = c(24,21,22))+
  #scale_shape_manual(values = c(21,22,23,24,25))+
  geom_segment(data = subset(arrow_4th_comp,arrow_4th_comp$p<0.05), aes(x=0, y=0, xend=MDS1*1.4, yend=MDS2*1.4),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  geom_segment(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=0, y=0, xend=MDS1*1.4, yend=MDS2*1.4),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=MDS1*1.4, y=MDS2*1.4, label=variable), size=8)+
  #geom_text(data = subset(arrow_4th_comp,arrow_4th_comp$p<0.05), aes(x=MDS1*1.4, y=MDS2*1.4, label=variable), size=8)+
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
  ggsave(filename("MDS_1st_withChp&Comp_nolabel"), height = 24, width = 20, units = "cm", dpi = 700)

##2. MDS select (Seoul, Seosan, Beijing)=====
###2-1NMDS======
table(FT_merge_sel$Group)

FT_merge_sel=FT_merge %>% filter(Group%in%c("SUL","SS","B"))
MDS_1st_sel=melt(FT_merge_sel[,c("Sample","Formula","Bromo.Inty")], id.vars = c("Sample","Formula")) %>% 
  dcast(Sample~Formula, sum)

MDS_1st_sel[,1:2]

MDS_1st_sel_in=MDS_1st_sel

NMDS_1st_sel=metaMDS(MDS_1st_sel_in[,-c(1)], k=5, distance = "bray", trymax = 100)
NMDS_1st_sel ##Stress 0.15

gnmds_1st_sel=as.data.table(NMDS_1st_sel$points)
gnmds_1st_sel$Sample=MDS_1st_sel_in$Sample

gnmds_1st_sel=gnmds_1st_sel %>% tidyr::separate("Sample",c("Group","no"),sep="_")
gnmds_1st_sel

ggplot()+
  geom_point(data = gnmds_1st_sel, aes(x=MDS1, y=MDS2, col=Group))


gnmds_1st_sel$Group=factor(gnmds_1st_sel$Group, levels = c("SUL","SS","B","NT","UL"),
                       labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
gnmds_1st_sel$Sample=paste(gnmds_1st_sel$Group,gnmds_1st_sel$no, sep = "_")

gnmds_1st_sel$id=ifelse(gnmds_1st_sel$Group=="Beijing",5,
                    ifelse(gnmds_1st_sel$Group=="Seosan",4,
                           ifelse(gnmds_1st_sel$Group=="seoul",3,2)))
gnmds_1st_sel=gnmds_1st_sel[order(gnmds_1st_sel$id)]

gnmds_1st_sel


ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st_sel, aes(x=MDS1, y=MDS2, fill=Group, shape=Group), col="black", size=9)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
  #scale_fill_npg()+
  scale_x_continuous(name = "NMDS1", breaks =seq(-1,1.0,0.5), limits = c(-1.0,1.0))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-0.6,0.6,0.3), limits = c(-0.65,0.65))+
  scale_shape_manual(values = c(24,21,22))+
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
  ggsave(filename("MDS_1st_sel"), height = 24, width = 20, units = "cm", dpi = 700)



p1=ggplot(data = gnmds_1st_sel, aes(x=MDS1, y=MDS2, fill=Group, shape=Group, id=Sample, Date=Date))+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(col="black", size=9)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
  #scale_fill_npg()+
  scale_x_continuous(name = "NMDS1", breaks =seq(-1,1.0,0.5), limits = c(-1.0,1.0))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-0.6,0.6,0.3), limits = c(-0.65,0.65))+
  scale_shape_manual(values = c(24,21,22))+
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
  )

ggplotly(p1)
####1-2+ NMDS with vector=====
vec_mds=gnmds_1st_sel[,c(1,2)]
vec_mds

pal=pal_npg("nrc")(5)

#####1-2-1 Enviromental data=====
gnmds_1st_sel


envi_1st

gnmds_1st_sel=gnmds_1st_sel %>% inner_join(envi_1st)

gnmds_1st_sel

gnmds_1st_sel$Event=factor(gnmds_1st_sel$Event,levels = c("Event","Normal","Non-event"))

vec_envi_sel=gnmds_1st_sel[,c("SO42-","NO3-","NH4+","CO","O3","SO2","NO","PM2.5")]
vec_envi_sel

vec_fit_envi=envfit(vec_mds,vec_envi_sel, na.rm = T)
vec_fit_envi$vectors$arrows
vec_fit_envi$vectors$r
vec_fit_envi$vectors$pval

arrow_4th_envi=as.data.frame(scores(vec_fit_envi, display = "vectors"))
arrow_4th_envi$variable=row.names(arrow_4th_envi)
arrow_4th_envi$R=vec_fit_envi$vectors$r
arrow_4th_envi$p=vec_fit_envi$vectors$pvals

arrow_4th_envi
arrow_4th_envi$variable=row.names(arrow_4th_envi)


gnmds_1st_sel$Group=factor(gnmds_1st_sel$Group, levels =  c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
gnmds_1st_sel$Event=factor(gnmds_1st_sel$Event,levels = c("Event","Normal","Non-event"))

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st_sel, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
  scale_x_continuous(name = "NMDS1", breaks =seq(-1,1.0,0.5), limits = c(-1.0,1.0))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-0.6,0.6,0.3), limits = c(-0.65,0.65))+
  scale_shape_manual(values = c(24,21,22))+
  #scale_shape_manual(values = c(21,22,23,24,25))+
  geom_segment(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=0, y=0, xend=MDS1*1, yend=MDS2*1),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=MDS1*0.88, y=MDS2*0.88, label=variable), size=8)+
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
  ggsave(filename("MDS_1st_sel_withenvi_nolabel2"), height = 24, width = 20, units = "cm", dpi = 300)

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st_sel, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
  scale_x_continuous(name = "NMDS1", breaks =seq(-1,1.0,0.5), limits = c(-1.0,1.0))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-0.6,0.6,0.3), limits = c(-0.65,0.65))+
  scale_shape_manual(values = c(24,21,22))+
  #scale_shape_manual(values = c(21,22,23,24,25))+
  geom_segment(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=0, y=0, xend=MDS1*1, yend=MDS2*1),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  geom_text(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=MDS1*0.88, y=MDS2*0.88, label=variable), size=8)+
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
  ggsave(filename("MDS_1st_sel_withenvi_label2"), height = 24, width = 20, units = "cm", dpi = 100)



#####2-2-2 chemical composition=====
vec_mds=gnmds_1st_sel[,c(1,2)]
vec_mds


chcomp
chcomp$Sample=paste(chcomp$Group,chcomp$no,sep = "_")

vec_comp=dcast(chcomp, Sample+Group+no~Comp, sum, value.var = "rel")

vec_comp=gnmds_1st_sel %>% inner_join(vec_comp[,-c("Group","no")], by=c("Sample"))

vec_comp=vec_comp[,c("CHO","CHON","CHOS","CHONS")]

vec_fit_comp=envfit(vec_mds,vec_comp, na.rm = T)
vec_fit_comp$vectors$arrows
vec_fit_comp$vectors$r
vec_fit_comp$vectors$pval

arrow_4th_comp=as.data.frame(scores(vec_fit_comp, display = "vectors"))
arrow_4th_comp$variable=row.names(arrow_4th_comp)
arrow_4th_comp$R=vec_fit_comp$vectors$r
arrow_4th_comp$p=vec_fit_comp$vectors$pvals

arrow_4th_comp
arrow_4th_comp$variable=row.names(arrow_4th_comp)

gnmds_1st_sel

arrow_4th_comp

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st_sel, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
  scale_x_continuous(name = "NMDS1", breaks =seq(-1,1.0,0.5), limits = c(-1.0,1.0))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-0.6,0.6,0.3), limits = c(-0.65,0.65))+
  scale_shape_manual(values = c(24,21,22))+
  #scale_shape_manual(values = c(21,22,23,24,25))+
  geom_segment(data = subset(arrow_4th_comp,arrow_4th_comp$p<0.05), aes(x=0, y=0, xend=MDS1, yend=MDS2),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=MDS1*0.88, y=MDS2*0.88, label=variable), size=8)+
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
  ggsave(filename("MDS_1st_sel_withComp_nolabel"), height = 24, width = 20, units = "cm", dpi = 700)


ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st_sel, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
  scale_x_continuous(name = "NMDS1", breaks =seq(-1,1.0,0.5), limits = c(-1.0,1.0))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-0.6,0.6,0.3), limits = c(-0.65,0.65))+
  scale_shape_manual(values = c(24,21,22))+
  #scale_shape_manual(values = c(21,22,23,24,25))+
  geom_segment(data = subset(arrow_4th_comp,arrow_4th_comp$p<0.05), aes(x=0, y=0, xend=MDS1, yend=MDS2),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  geom_text(data = subset(arrow_4th_comp,arrow_4th_comp$p<0.05), aes(x=MDS1*0.88, y=MDS2*0.88, label=variable), size=8)+
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
  ggsave(filename("MDS_1st_sel_withComp_label"), height = 24, width = 20, units = "cm", dpi = 100)

#####2-2-3 chemical composition=====
vec_mds=gnmds_1st_sel[,c(1,2)]
vec_mds

chp_m
vec_chp=dcast(chp_m, Sample+Group+no~variable, sum, value.var = "val")
vec_chp

vec_chp=gnmds_1st_sel %>% inner_join(vec_chp[,-c("Group","no")], by=c("Sample"))

vec_chp=vec_chp[,c("Mean AI","Mean DBE","Mean H/C","Mean N/C","Mean O/C","Mean S/C")]

vec_fit_chp=envfit(vec_mds,vec_chp, na.rm = T)
vec_fit_chp$vectors$arrows
vec_fit_chp$vectors$r
vec_fit_chp$vectors$pval

arrow_4th_chp=as.data.frame(scores(vec_fit_chp, display = "vectors"))
arrow_4th_chp$variable=row.names(arrow_4th_chp)
arrow_4th_chp$R=vec_fit_chp$vectors$r
arrow_4th_chp$p=vec_fit_chp$vectors$pvals

arrow_4th_chp
arrow_4th_chp$variable=row.names(arrow_4th_chp)

gnmds_1st_sel

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st_sel, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
  scale_x_continuous(name = "NMDS1", breaks =seq(-1,1.0,0.5), limits = c(-1.0,1.0))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-0.6,0.6,0.3), limits = c(-0.65,0.65))+
  scale_shape_manual(values = c(24,21,22))+
  #scale_shape_manual(values = c(21,22,23,24,25))+
  geom_segment(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=0, y=0, xend=MDS1*0.8, yend=MDS2*0.8),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=MDS1*0.88, y=MDS2*0.88, label=variable), size=8)+
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
  ggsave(filename("MDS_1st_sel_withChp_nolabel"), height = 24, width = 20, units = "cm", dpi = 700)




ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st_sel, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
  scale_x_continuous(name = "NMDS1", breaks =seq(-1,1.0,0.5), limits = c(-1.0,1.0))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-0.6,0.6,0.3), limits = c(-0.65,0.65))+
  scale_shape_manual(values = c(24,21,22))+
  geom_segment(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=0, y=0, xend=MDS1*0.8, yend=MDS2*0.8),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  geom_text(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=MDS1*0.88, y=MDS2*0.88, label=variable), size=8)+
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
  ggsave(filename("MDS_1st_sel_withChp_label"), height = 24, width = 20, units = "cm", dpi = 100)

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st_sel, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_manual(values = c("#E64B35FF", "#00A087FF","#4DBBD5FF","#F39B7FFF","#3C5488FF"))+
  scale_x_continuous(name = "NMDS1", breaks =seq(-1,1.0,0.5), limits = c(-1.0,1.0))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-0.6,0.6,0.3), limits = c(-0.70,0.70))+
  scale_shape_manual(values = c(24,21,22))+
  #scale_shape_manual(values = c(21,22,23,24,25))+
  geom_segment(data = subset(arrow_4th_comp,arrow_4th_comp$p<0.05), aes(x=0, y=0, xend=MDS1, yend=MDS2),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  geom_segment(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=0, y=0, xend=ifelse(variable=="Mean H/C",MDS1*0.8,
                                                                                            MDS1), 
                                                                      yend=ifelse(variable=="Mean H/C",MDS2*0.8,
                                                                                  MDS2)),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=MDS1*0.88, y=MDS2*0.88, label=variable), size=8)+
  #geom_text(data = subset(arrow_4th_comp,arrow_4th_comp$p<0.05), aes(x=MDS1*0.88, y=MDS2*0.88, label=variable), size=8)+
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
  ggsave(filename("MDS_1st_sel_withChp&Comp_nolabel"), height = 24, width = 20, units = "cm", dpi = 700)


