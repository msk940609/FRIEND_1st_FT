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
library(ggtext)
library(VennDiagram)
library("ggVennDiagram")

getOption("digits")
options("digits" = 15)

source("Rscript/func_filename.R")
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")

##NMDS=====
ft_merge$Sample=paste(ft_merge$Group,ft_merge$No, sep = "_")

###nmds with vector===========
MDS_1st=melt(ft_merge[,c("Sample","Formula","Bromo.Inty")], id.vars = c("Sample","Formula")) %>% 
  dcast(Sample~Formula, sum)

MDS_1st[,1:2]
MDS_1st_in=MDS_1st

NMDS_1st=metaMDS(MDS_1st_in[,-c(1)], k=2, distance = "bray", trymax = 20)
NMDS_1st ##Stress 0.03

gnmds_1st=as.data.table(NMDS_1st$points)
gnmds_1st$Sample=MDS_1st_in$Sample

gnmds_1st=gnmds_1st %>% tidyr::separate("Sample",c("Group","No"),sep="_")
gnmds_1st$No=as.numeric(gnmds_1st$No)

gnmds_1st$Group=factor(gnmds_1st$Group, 
                       levels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))

gnmds_1st$Sample=paste(gnmds_1st$Group,gnmds_1st$No, sep = "_")

envi_1st=fread("Datafile/FRIEND_1st_envi_re.csv")
envi_1st$Sample=paste(envi_1st$Group,envi_1st$No, sep = "_")


gnmds_1st=gnmds_1st %>% inner_join(envi_1st, by = c("Sample","Group","No"))
#gnmds_1st=gnmds_1st[order(gnmds_1st$id)]
gnmds_1st

gnmds_1st$Event=factor(gnmds_1st$Event,levels = c("Event","Normal","Non-event"))
                       

gnmds_1st_ord=gnmds_1st


gnmds_1st_ord$id=ifelse(gnmds_1st_ord$Group=="Beijing",5,
                        ifelse(gnmds_1st_ord$Group=="Seosan",4,
                               ifelse(gnmds_1st_ord$Group=="seoul",3,2)))
gnmds_1st$Group=factor(gnmds_1st$Group, 
                       levels = c("Seoul","Seosan","Beijing","Ulaanbaatar","Noto"))

gnmds_1st_ord=gnmds_1st[order(gnmds_1st$PM2.5),]

mypal = pal_npg("nrc", alpha = 1)(9)

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st_ord, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_npg()+
  scale_x_continuous(name = "NMDS1", breaks =seq(-2,2,0.5), limits = c(-1.0,2))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-1.5,1.5,0.5), limits = c(-1.2,1.2))+
  scale_shape_manual(values = c(24,21,22))+
  scale_fill_manual(values = c("#4DBBD5FF","#00A087FF","#E64B35FF","#3C5488FF","#F39B7FFF"))+
  #scale_shape_manual(values = c(21,22,23,24,25))+
  #geom_segment(data = subset(arrow_4th_chp,arrow_4th_chp$p<0.05), aes(x=0, y=0, xend=MDS1*0.8, yend=MDS2*0.8),
  #             arrow=arrow(length = unit(0.25, "cm")), size=1)+
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
        legend.text = element_text(size = 18, colour = "black", family = "Arial",vjust = 0.6,margin = unit(c(0.0,-0.2,0.0,-0.2),"cm")),
        legend.spacing = unit(0.2,"cm"),
        legend.box.margin = margin(t=0,r=0,b=0,l=0.0, unit = "cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "right",
        legend.box.just = c(0.5,0.5),
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_blank()
  )+
  guides(fill=guide_legend(title = NULL,title.hjust = 0,title.position = "top",byrow = T , ncol = 1,
                           override.aes = list(size=7, shape=21), order = 1),
         shape=guide_legend(title = NULL,title.position = "top",override.aes = list(size=8, fill="black"), ncol=1))+
  ggsave(filename("MDS_1st"), height = 20, width = 26, units = "cm", dpi = 700)

fm_obs=ft_merge[,c("Group","No","Formula","Freq")] %>% 
  dcast(Group~Formula,value.var = "Freq",sum) %>% 
  melt(id.vars=c("Group")) %>% `colnames<-`(c("Group","Formula","Freq"))

fm_obs_sel=fm_obs %>% filter(Freq>3)
fm_obs_sel

fm_obs_sel_u=subset(fm_obs_sel,fm_obs_sel$Group=="Ulaanbaatar")
fm_obs_sel_b=subset(fm_obs_sel,fm_obs_sel$Group=="Beijing")
fm_obs_sel_sul=subset(fm_obs_sel,fm_obs_sel$Group=="Seoul")
fm_obs_sel_ss=subset(fm_obs_sel,fm_obs_sel$Group=="Seosan")
fm_obs_sel_nt=subset(fm_obs_sel,fm_obs_sel$Group=="Noto")

fm_u=as.vector(fm_obs_sel_u$Formula)
fm_b=as.vector(fm_obs_sel_b$Formula)
fm_sul=as.vector(fm_obs_sel_sul$Formula)
fm_ss=as.vector(fm_obs_sel_ss$Formula)
fm_nt=as.vector(fm_obs_sel_nt$Formula)

tt <- list(
  Ulaanbaatar = fm_u,
  Beijing = fm_b,
  Seoul = fm_sul, 
  Seosan = fm_ss,
  Noto=fm_nt
)

display_venn <- function(x, fn=NULL,...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = fn, ...)
  grid.draw(venn_object)
}

venn.diagram(tt,
             category.names = c("Ulaanbaatar" , "Beijing" , "Seoul", "Seosan","Noto"),
             #fill =  c("#3C5488FF","#E64B35FF","#4DBBD5FF","#00A087FF","#F39B7FFF"),
             fill =  c("white","white","white","white","white"),
             cat.dist = c(0.1, 0.1, 0.1, 0.1,0.1),
             lty = 1,  lwd = 1,
             cat.cex=0,
             #cex=1.2,
             cat.fontfamily= "Arial",main.fontface = "bold",sub.fontface = "bold",
             fontfamily= "Arial",
             filename = "venn-5-dimensions.tiff")


venn=fread("Datafile/r_1st_venn.csv")
venn_m=melt(venn, id.vars = c("Group"), na.rm = T)
venn_m$Group=factor(venn_m$Group,levels = rev(c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto")))

ggplot(venn_m, aes(x=variable,y=Group, col=Group))+
  geom_line(aes(group=variable), col="black")+
  geom_point()+
  scale_x_discrete(name="",position="top")+
  scale_y_discrete(name="")+
  scale_color_manual(values = rev(c("#3C5488FF","#E64B35FF","#4DBBD5FF","#00A087FF","#F39B7FFF")))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.05, "cm"),
        axis.text.x = element_text(size = 10, colour = "black",face = "bold", family = "Arial", margin =unit(c(0.2,0.2,0.2,0.2),"cm") ),
        axis.text.y = element_text(size = 10, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        legend.text = element_text(size = 0, colour = "black", family = "Arial",vjust = 0.6,margin = unit(c(0.0,-0.2,0.0,-0.2),"cm")),
        legend.spacing = unit(0.2,"cm"),
        legend.box.margin = margin(t=0,r=0,b=0,l=0.0, unit = "cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "right",
        legend.box.just = c(0.5,0.5),
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_blank()
  )+
  guides(col="none")+
  ggsave(filename("ven_group"), height = 4, width = 20, units = "cm", dpi = 300)

fm_obs_sel
venn_fm=dcast(fm_obs_sel,Formula~Group,sum)

venn_fm$U=ifelse(venn_fm$Ulaanbaatar>0, "U","")
venn_fm$B=ifelse(venn_fm$Beijing>0, "B","")
venn_fm$S=ifelse(venn_fm$Seoul>0, "S","")
venn_fm$SS=ifelse(venn_fm$Seosan>0, "SS","")
venn_fm$N=ifelse(venn_fm$Noto>0, "N","")

venn_fm=venn_fm %>% unite( "ven_class",c("U","B","S","SS","N"),sep = "")

table(venn_fm$ven_class)

fwrite(venn_fm,file = "Datafile/venn_classlist.csv")

venn_class=fread("Datafile/venn_classlist_sel.csv")
fmlist

ft_vk_m=ft_vk_m %>% mutate(Comp = gsub('[[:digit:]]', '', Formula))
ft_vk_m$Comp=ifelse(ft_vk_m$`O`==0,"Remainders",ft_vk_m$Comp)
ft_vk_m$Comp=factor(ft_vk_m$Comp,levels = c("CHO","CHNO","CHOS","CHNOS","Remainders"),
                    labels = c("CHO","CHON","CHOS","CHONS","Remainders"))

venn_class=venn_class %>% mutate(Comp = gsub('[[:digit:]]', '', Formula))

venn_class$Comp=factor(venn_class$Comp,levels = c("CHO","CHNO","CHOS","CHNOS","CH","CHN","CHNS","CHS"),
                    labels = c("CHO","CHON","CHOS","CHONS","Remainders","Remainders","Remainders","Remainders"))

table(venn_class$Comp)

venn_class$Freq=1

venn_class_sum=as.data.table(aggregate(venn_class$Freq, by=list(variable=venn_class$variable, 
                                                                Comp=venn_class$Comp),sum))
venn_class_sum

ggplot(venn_class_sum, aes(x=variable,y=x, fill=Comp))+
  geom_bar(stat = "identity", position = position_stack(reverse = T))+
  scale_x_discrete(name="",position="bottom")+
  scale_y_continuous(name="The number of molecules", expand=c(0.01,0.01), limits = c(0,1700))+
  scale_fill_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF","grey50"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        plot.margin = unit(c(1.0,0.2,0.2,0.2),"cm"),
        axis.ticks.length = unit(0.05, "cm"),
        plot.title = element_text(size=0),
        axis.text.x = element_text(size = 10, colour = "black",face = "bold", family = "Arial", margin =unit(c(0.2,0.2,0.2,0.2),"cm") ),
        axis.text.y = element_text(size = 10, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 10, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 14, colour = "black",face = "bold", family = "Arial",
                                    margin =unit(c(0.2,0.4,0.2,0.2),"cm")),
        legend.text = element_text(size = 10, colour = "black",face = "bold", family = "Arial",hjust = 0.5,margin = unit(c(0.0,-0.2,0.0,-0.2),"cm")),
        legend.spacing = unit(0.2,"cm"),
        legend.box.margin = margin(t=0,r=0,b=0,l=0.0, unit = "cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.5,1.05),
        legend.box.just = c(0.5,0.5),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.background = element_blank()
  )+
  guides(fill=guide_legend(title = "", reverse = F))+
  ggsave(filename("ven_group_bar"), height = 10, width = 20, units = "cm", dpi = 300)
