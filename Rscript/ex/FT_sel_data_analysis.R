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
library(plotly)

source("Rscript/func_filename.R")

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

FT_merge=fread("Datafile/FRIENDs_1st_FT.csv")

FT_merge$C=ifelse(FT_merge$C.>0, "C","")
FT_merge$H=ifelse(FT_merge$H.>0, "H","")
FT_merge$O=ifelse(FT_merge$O.>0, "O","")
FT_merge$N=ifelse(FT_merge$N.>0, "N","")
FT_merge$S=ifelse(FT_merge$S.>0, "S","")

FT_merge=FT_merge %>% unite("Comp",c("C","H","O","N","S"),sep = "")
FT_merge$Comp=ifelse(FT_merge$O.==0, "Remainders",FT_merge$Comp)

FT_merge=subset(FT_merge,FT_merge$Comp!="Remainders") %>% droplevels()
FT_merge

####nmds_sel=====
table(FT_merge$Group)
FT_merge_sel=FT_merge %>% filter(Group%in%c("SUL","SS","B"))

table(FT_merge_sel$Group)

MDS_sel_1st=melt(FT_merge_sel[,c("Sample","Formula","Bromo.Inty")], id.vars = c("Sample","Formula")) %>% 
  dcast(Sample~Formula, sum)

MDS_sel_1st[,1:2]
MDS_sel_1st_in=MDS_sel_1st

NMDS_sel_1st=metaMDS(MDS_sel_1st_in[,-c(1)], k=3, distance = "bray", trymax = 100)
NMDS_sel_1st ##k=2 Stress 0.15 ,k=3 stress:0.08

gnMDS_sel_1st=as.data.table(NMDS_sel_1st$points)
gnMDS_sel_1st$Sample=MDS_sel_1st_in$Sample

gnMDS_sel_1st=gnMDS_sel_1st %>% tidyr::separate("Sample",c("Group","no"),sep="_")
gnMDS_sel_1st

ggplot()+
  geom_point(data = gnMDS_sel_1st, aes(x=MDS1, y=MDS2, col=Group))

gnMDS_sel_1st$Group=factor(gnMDS_sel_1st$Group, levels = c("SUL","SS","B"),
                           labels = c("Seoul","Seosan","Beijing"))
gnMDS_sel_1st$Sample=paste(gnMDS_sel_1st$Group,gnMDS_sel_1st$no,sep = "_")


gnMDS_sel_1st=gnMDS_sel_1st %>% inner_join(envi_1st, by=c("Sample","Group"))
gnMDS_sel_1st

gnMDS_sel_1st$PM2.5=round(gnMDS_sel_1st$PM2.5,1)

gnMDS_sel_1st$Event=factor(gnMDS_sel_1st$Event, levels = c("Event","Normal","Non-event"))

gnMDS_sel_1st$Group=factor(gnMDS_sel_1st$Group, levels=c("Seoul","Seosan","Beijing"))

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnMDS_sel_1st, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_npg()+
  scale_shape_manual(values = c(24,21,22))+
  #geom_segment(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=0, y=0, xend=MDS_sel1*0.9, yend=MDS_sel2*0.9),
  #             arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=MDS_sel1*0.98, y=MDS_sel2*0.98, label=variable), size=8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(0.2,0.6,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,-0.2,0.1,-0.2),"cm")),
        #legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "bottom",
        legend.spacing.y = unit(0.1,"cm"),
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_blank()
  )+
  xlab("NMDS1")+
  ylab("NMDS2")+
  guides(fill=guide_legend(title = NULL,title.hjust = 0,title.position = "top",byrow = T , ncol = 3,
                           override.aes = list(size=7, shape=21), order = 1),
         shape=guide_legend(title = NULL,title.position = "top",override.aes = list(size=8, fill="black"), ncol=3))+
  ggsave(filename("MDS_sel_1st"), height = 22, width = 20, units = "cm", dpi = 700)



gnMDS_sel_1st

envi_1st=fread("Datafile/envi_1st_sel.csv")

gnMDS_sel_1st$Sample=paste(gnMDS_sel_1st$Group,gnMDS_sel_1st$no,sep = "_")
gnMDS_sel_1st=gnMDS_sel_1st %>% inner_join(envi_1st)
gnMDS_sel_1st

gnMDS_sel_1st$PM2.5=round(gnMDS_sel_1st$PM2.5,1)

gnMDS_sel_1st$Group=factor(gnMDS_sel_1st$Group,levels = c("Seoul","Seosan","Beijing"))

p1=ggplot(data = gnMDS_sel_1st, aes(x=MDS1, y=MDS2, fill=Group, shape=Group,Sample=Sample,PM=PM2.5, Event=Event))+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(col="black", size=9)+
  scale_fill_npg()+
  #scale_fill_manual(values = c("#e97f02","#0080ff"))+
  scale_shape_manual(values = c(21,22,23))+
  #geom_segment(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=0, y=0, xend=MDS_sel1*0.9, yend=MDS_sel2*0.9),
  #             arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=MDS_sel1*0.98, y=MDS_sel2*0.98, label=variable), size=8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(0.2,0.6,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,-0.2,0.1,-0.2),"cm")),
        #legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.86, 0.14),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5,byrow = T , ncol = 1,
                           override.aes = list(size=7, shape=c(21,22,23))),
         shape="none")

ggplotly(p1,tooltip = c("x", "y","Formula", "Sample","PM","Event")) %>% 
  layout(legend = list(orientation = "h", x = 0.3, y = -0.2))


####p2 fill&shape logic error======
p2=ggplot(data = gnMDS_sel_1st, aes(x=MDS1, y=MDS2, fill=Group, shape=Event, Date=Date, PM2.5=PM2.5))+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(col="black", size=9)+
  scale_fill_npg()+
  scale_shape_manual(values = c(24,21,22))+
  #geom_segment(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=0, y=0, xend=MDS_sel1*0.9, yend=MDS_sel2*0.9),
  #             arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=MDS_sel1*0.98, y=MDS_sel2*0.98, label=variable), size=8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(0.2,0.6,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,-0.2,0.1,-0.2),"cm")),
        #legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "bottom",
        legend.spacing.y = unit(0.1,"cm"),
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_blank()
  )+
  guides(fill=guide_legend(title = NULL,title.hjust = 0,title.position = "top",byrow = T , ncol = 3,
                           override.aes = list(size=7, shape=21), order = 1),
         shape=guide_legend(title = NULL,title.position = "top",override.aes = list(size=8, fill="black"), ncol=3))

ggplotly(p2,tooltip = c("x", "y","Formula", "Sample","PM2.5","Event","Date")) %>% 
  layout(legend = list(orientation = "h", x = 0.3, y = -0.2))

###vkplot====

insert_minor <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_1 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_2 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}

ggplot(nmds_formula_spec, aes(x=`O.C`, y=`H.C`, col=Comp))+
  geom_point(size=3, alpha=0.4)+
  facet_wrap(.~Region, ncol = 3, scales = "fix")+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #scale_fill_manual(values = c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF","grey50"))+
  scale_color_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF"))+
  #scale_size_continuous(range = c(1,16))+
  #ggtitle("Beijing Event")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.text = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.22,0.60),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("Specific_vk"),height = 30, width = 45, units = "cm", dpi = 700)
  

FT_pm=fread("Datafile/FRIENDs_envi.csv")

FT_pm
gnmds_1st

gnmds_1st=gnmds_1st %>% inner_join(FT_pm, by = "Sample")
gnmds_1st

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st, aes(x=MDS1, y=MDS2, fill=PM2.5, shape=Group), col="black", size=9)+
  scale_shape_manual(values = c(21,22,23,24,25))+
  scico::scale_fill_scico(palette = "lajolla", begin = 0.05, end = 0.95)+
  #scale_fill_gradientn(colors = c("#e97f02","#e97f02","#e97f02","#0080ff","#0080ff"))+
  #geom_segment(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=0, y=0, xend=MDS1*0.9, yend=MDS2*0.8),
  #             arrow=arrow(length = unit(0.25, "cm")), size=1.0)+
  #geom_text(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=MDS1*0.75, y=MDS2*0.75, label=variable), size=6)+
  #geom_text(data = gnmds_4th, aes(x=MDS1, y=MDS2, label=Sample), size=6)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(0.2,0.3,0.2,0.2),"cm"),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.72, 0.16),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(fill=guide_colorbar(title = expression(bold(PM[2.5])), order = 2),shape=guide_legend(title = "Group",order = 1,override.aes = list(size=6,fill="black", col="black")))+
  ggsave(filename("MDS_1st_col_gradient_nonlabel"), height = 20, width = 20, units = "cm", dpi = 700)



#####Composition=====
ft_inty_all=as.data.table(aggregate(FT_merge$Bromo.Inty, by=list(`Sample`=FT_merge$Sample), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Tot"))
ft_inty=as.data.table(aggregate(FT_merge$Bromo.Inty, by=list(`Sample`=FT_merge$Sample, Comp=FT_merge$Comp), FUN=sum))
ft_inty

ft_inty=ft_inty %>% inner_join(ft_inty_all,by = "Sample")
ft_inty$rel=round(ft_inty$x/ft_inty$Tot*100,1)

ft_inty=ft_inty %>% separate(Sample, c("Group","no"),sep = "_")
ft_inty

ft_inty$Comp=factor(ft_inty$Comp,levels = c("CHO","CHON","CHOS","CHONS"))

ft_inty$Group=factor(ft_inty$Group,levels = c("SUL","SS","B","NT","UL"),
                     labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
ft_inty
stat

source("Rscript/func_upper_fence_label.R")

norm=data.table()
for (i in 1:length(comp)) {
  temp=subset(ft_inty,ft_inty$Comp==comp[i])
  grp=unique(temp$Group)
  for (j in 1:length(grp)) {
    tmp2=subset(temp,temp$Group==grp[j])
    stmp=shapiro.test(tmp2$rel)
    new=data.table("Comp"=comp[i],"Group"=grp[j],"p"=round(stmp$p.value,3))
    norm=rbind(norm,new)
  }
}

comp=unique(ft_inty$Comp)
stat=data.table()
for (i in 1:length(comp)) {
  temp=subset(ft_inty,ft_inty$Comp==comp[i])
  
  new=generate_label_df(temp, flev = "Group",value.var = "rel",offset = 4, pair.method="hsd")
  new$Comp=comp[i]
  stat=rbind.data.frame(stat, new)
  
}

stat=stat %>% `colnames<-`(c("Group","labels","rel","Comp"))
stat$Comp=factor(stat$Comp,levels = c("CHO","CHON","CHOS","CHONS"))
stat


comp_max=as.data.table(aggregate(ft_inty$rel, by=list(`Comp`=ft_inty$Comp), FUN=max))
comp_max$x=comp_max$x*1.2

ggplot()+
  stat_boxplot(data=ft_inty, aes(x=Group, y=rel, fill=Group),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=ft_inty, aes(x=Group, y=rel, fill=Group),alpha=1, outlier.color = NA)+
  geom_text(data =comp_max, aes(x=1, y=x, label=Comp), col="white", size=0)+
  geom_text(data =stat, aes(x=Group, y=rel, label=labels), 
            col="black", size=6)+
  facet_wrap(.~Comp, scales = "fix", nrow=2, dir = "h")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  scale_fill_npg()+
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
  ggsave(filename("Comp_compare_hsd"),height = 20, width = 30, units = "cm", dpi = 300)

####chemical properties====

FT_merge

chp_m=melt(FT_merge[,c("Sample","O.C","H.C","N.C","S.C","DBE","AI")], id.vars = c("Sample")) %>% 
  dcast(Sample~variable,mean) %>% `colnames<-`(c("Sample","Mean O/C","Mean H/C", "Mean N/C", "Mean S/C", "Mean DBE", "Mean AI")) %>% 
  melt(id.var=c("Sample"))

chp_m=chp_m %>% separate(Sample, c("Group","No"),sep = "_")
chp_m$Sample=paste(chp_m$Group,chp_m$No,sep = "_")

chp_m

chp_m$Group=factor(chp_m$Group,levels = c("SUL","SS","B","NT","UL"),
                     labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
chp_m=chp_m %>% `colnames<-`(c("Group","No","variable","val","Sample"))

norm=data.table()
chp=unique(chp_m$variable)
grp=unique(chp_m$Group)

for (i in 1:length(chp)) {
  temp=subset(chp_m,chp_m$variable==chp[i])
  grp=unique(temp$Group)
  for (j in 1:length(grp)) {
    tmp2=subset(temp,temp$Group==grp[j])
    stmp=shapiro.test(tmp2$val)
    new=data.table("Comp"=chp[i],"Group"=grp[j],"p"=round(stmp$p.value,3))
    norm=rbind(norm,new)
  }
}

norm

chp=unique(chp_m$variable)
stat=data.table()
for (i in 1:length(chp)) {
  temp=subset(chp_m,chp_m$variable==chp[i])
  new=generate_label_df(temp, flev = "Group",value.var = "val",offset = 0, pair.method="duncan")
  new$chp=chp[i]
  stat=rbind.data.frame(stat, new)
  
}

stat=stat %>% `colnames<-`(c("Group","labels","val","variable"))
#stat$variable=factor(stat$variable, levels = c("Mean O/C", ))

chp_max=as.data.table(aggregate(chp_m$val, by=list(`variable`=chp_m$variable), FUN=max))
chp_max$x=chp_max$x*1.1
chp_max$min=0

ggplot()+
  stat_boxplot(data=chp_m, aes(x=Group, y=val, fill=Group),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=chp_m, aes(x=Group, y=val, fill=Group), outlier.color = NA)+
  geom_text(data =chp_max, aes(x=1, y=x, label=variable), col="white", size=0)+
  #geom_text(data =chp_max, aes(x=1, y=min, label=variable), col="white", size=0)+
  geom_text(data =stat, aes(x=Group, y=ifelse(variable=="Mean O/C",val+0.02,
                                              ifelse(variable=="Mean H/C",val+0.05,
                                                     ifelse(variable=="Mean N/C",val+0.008,
                                                            ifelse(variable=="Mean S/C",val+0.006,
                                                                   ifelse(variable=="Mean DBE",val+0.8,val+0.02))))), label=labels), 
            col="black", size=6)+
  facet_wrap(.~variable, scales = "free", nrow=2, dir = "v")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  scale_fill_npg()+
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
  ggsave(filename("Chp_compare_duncan"),height = 20, width = 45, units = "cm", dpi = 300)

chp_m
chp_m=chp_m %>% `colnames<-`(c("Group","No","variable","val","Sample"))

chp=unique(chp_m$variable)
stat=data.table()
for (i in 1:length(chp)) {
  temp=subset(chp_m,chp_m$variable==chp[i])
  
  new=generate_label_df(temp, flev = "Group",value = "val",offset = 1.1, method="tukey")
  new$chp=chp[i]
  stat=rbind.data.frame(stat, new)
  
}

stat=stat %>% `colnames<-`(c("Group","labels","val","variable"))
#stat$variable=factor(stat$variable, levels = c("Mean O/C", ))

norm=data.table()
for (i in 1:length(comp)) {
  temp=subset(ft_inty,ft_inty$Comp==comp[i])
  grp=unique(temp$Group)
  for (j in 1:length(grp)) {
    tmp2=subset(temp,temp$Group==grp[j])
    stmp=shapiro.test(tmp2$rel)
    new=data.table("Comp"=comp[i],"Group"=grp[j],"p"=round(stmp$p.value,3))
    norm=rbind(norm,new)
  }
}

norm
library(DescTools)
####vk plot======

insert_minor <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_1 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_2 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}

table(FT_merge_sel$Group)

FT_merge_sel=FT_merge %>% filter(Group%in%c("SUL","SS","B"))

FT_merge_sel$Group=factor(FT_merge_sel$Group,levels = c("SUL","SS","B","NT","UL"),
                   labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))

FT_merge_sel=FT_merge_sel %>% droplevels()
FT_merge_sel$Sample=paste(FT_merge_sel$Group,FT_merge_sel$No, sep = "_")
envi_1st

FT_merge_sel_envi=FT_merge_sel %>% inner_join(envi_1st[,c("Sample","Group","No","Date")])

table(FT_merge_sel_envi$Date)

FT_red=FT_merge_sel_envi %>% filter(Date=="2020-12-23"|Date=="2020-12-28")
FT_red=FT_red %>% filter(Group%in%c("Seoul","Seosan"))

FT_blue=FT_merge_sel_envi %>% filter(Sample%in%c("Seoul_29","Seosan_30"))
FT_blue=FT_blue %>% filter(Group%in%c("Seoul","Seosan"))

FT_black=FT_merge_sel_envi %>% filter(Sample%in%c("Seosan_3","Seosan_6","Beijing_13","Beijing_8"))


ft_red_vk=melt(FT_red[,c("Sample","Group","Comp","Formula","O.C","H.C","Bromo.Inty")], 
           id.vars = c("Sample","Group","Comp","Formula","O.C","H.C")) %>% 
  dcast(Group+Sample~Formula, mean) %>% 
  melt(id.var=c("Group","Sample"), na.rm = T) %>% `colnames<-`(c("Group","Sample","Formula","Inty"))
ft_red_vk

fm=unique(FT_merge[,c("Formula","O.C","H.C","DBE","C.","Comp","Calc.m.z")])
fm

ft_red_vk=ft_red_vk %>% inner_join(fm, by="Formula")
ft_red_vk


ft_red_vk$Group=factor(ft_red_vk$Group,levels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
ft_red_vk

ft_red_vk$Sample=factor(ft_red_vk$Sample,levels = c("Seoul_9","Seosan_9","Seoul_14","Seosan_14"),
                        labels = c("Seoul 12/23","Seosan 12/23","Seoul 12/28","Seosan 12/28"))
ft_red_vk

ft_red_vk$Comp=factor(ft_red_vk$Comp, levels = c("CHO","CHON","CHOS","CHONS"))
ft_red_vk=ft_red_vk[order(ft_red_vk$Inty,decreasing=F),]

ggplot(ft_red_vk,aes(x=Inty))+
  geom_histogram(binwidth = 10)+coord_cartesian(ylim = c(0,100), xlim=c(0,100))

ft_red_vk$sf=ifelse(ft_red_vk$Inty<10,"0 - 10",
                    ifelse(ft_red_vk$Inty<100,"10 - 100",
                           ifelse(ft_red_vk$Inty<200,"100 - 200",
                                  ifelse(ft_red_vk$Inty<400,"200 - 400",
                                         ifelse(ft_red_vk$Inty<800,"400 - 800","800 <")))))

#ft_red_vk$sf=factor(ft_red_vk$sf,levels = c("0 - 100","100 - 300","300 - 600","600 - 900","900 - 1200",
#                                            "1200 <"))

ggplot(ft_red_vk, aes(x=`O.C`, y=`H.C`, col=Comp))+
  geom_point(aes(size=sf), alpha=0.4)+
  #geom_point(aes(size=sf), alpha=0.4)+
  facet_wrap(.~Sample, ncol = 5, scales = "fix")+
  #scale_size_continuous(range=c(1,30), breaks = seq(0,1000,100))+
  scale_size_manual(values=c(2,4,8,16,24,32))+
  #scale_size_continuous(range=c(1,30))+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #scale_fill_manual(values = c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF","grey50"))+
  scale_color_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF"))+
  #scale_size_continuous(range = c(1,16))+
  #ggtitle("Beijing Event")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.text = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "NULL",
        #legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)))+
  ggsave(filename("vk_red_noleg"),height = 18, width = 75, units = "cm", dpi = 700)


ft_blue_vk=melt(FT_blue[,c("Sample","Group","Comp","Formula","O.C","H.C","Bromo.Inty")], 
               id.vars = c("Sample","Group","Comp","Formula","O.C","H.C")) %>% 
  dcast(Group+Sample~Formula, mean) %>% 
  melt(id.var=c("Group","Sample"), na.rm = T) %>% `colnames<-`(c("Group","Sample","Formula","Inty"))
ft_blue_vk

fm=unique(FT_merge[,c("Formula","O.C","H.C","DBE","C.","Comp","Calc.m.z")])
fm

ft_blue_vk=ft_blue_vk %>% inner_join(fm, by="Formula")
ft_blue_vk


ft_blue_vk$Group=factor(ft_blue_vk$Group,levels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
ft_blue_vk
table(ft_blue_vk$Sample)

ft_blue_vk$Sample=factor(ft_blue_vk$Sample,levels = c("Seoul_29","Seosan_30"),
                        labels = c("Seoul 1/12","Seosan 1/13"))
ft_blue_vk

ft_blue_vk$Comp=factor(ft_blue_vk$Comp, levels = c("CHO","CHON","CHOS","CHONS"))
ft_blue_vk=ft_blue_vk[order(ft_blue_vk$Inty,decreasing=F),]

ggplot(ft_blue_vk,aes(x=Inty))+
  geom_histogram(binwidth = 100)+ylim(c(0,10))

ft_blue_vk$sf=ifelse(ft_blue_vk$Inty<10,"0 - 10",
                    ifelse(ft_blue_vk$Inty<100,"10 - 100",
                           ifelse(ft_blue_vk$Inty<200,"100 - 200",
                                  ifelse(ft_blue_vk$Inty<400,"200 - 400",
                                         ifelse(ft_blue_vk$Inty<800,"400 - 800","800 <")))))

#ft_blue_vk$sf=factor(ft_blue_vk$sf,levels = c("0 - 100","100 - 300","300 - 600","600 - 900","900 - 1200",
#                                            "1200 <"))


ggplot(ft_blue_vk, aes(x=`O.C`, y=`H.C`, col=Comp))+
  geom_point(aes(size=sf), alpha=0.4)+
  facet_wrap(.~Sample, ncol = 5, scales = "fix")+
  scale_size_manual(values=c(2,6,10,14,18,22))+
  #scale_size_continuous(range=c(1,30))+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #scale_fill_manual(values = c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF","grey50"))+
  scale_color_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF"))+
  #scale_size_continuous(range = c(1,16))+
  #ggtitle("Beijing Event")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.text = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "NULL",
        #legend.position = "bottom",
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("vk_blue_noleg"),height = 18, width = 37, units = "cm", dpi = 700)



ft_black_vk=melt(FT_black[,c("Sample","Group","Comp","Formula","O.C","H.C","Bromo.Inty")], 
                id.vars = c("Sample","Group","Comp","Formula","O.C","H.C")) %>% 
  dcast(Group+Sample~Formula, mean) %>% 
  melt(id.var=c("Group","Sample"), na.rm = T) %>% `colnames<-`(c("Group","Sample","Formula","Inty"))
ft_black_vk

fm=unique(FT_merge[,c("Formula","O.C","H.C","DBE","C.","Comp","Calc.m.z")])
fm

ft_black_vk=ft_black_vk %>% inner_join(fm, by="Formula")
ft_black_vk


ft_black_vk$Group=factor(ft_black_vk$Group,levels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
ft_black_vk
table(ft_black_vk$Sample)

ft_black_vk$Sample=factor(ft_black_vk$Sample,levels = c("Seosan_3","Beijing_13","Seosan_6","Beijing_8"),
                         labels = c("Seosan 12/17","Beijing 12/27","Seosan 12/20","Beijing 12/22"))
ft_black_vk

ft_black_vk$Comp=factor(ft_black_vk$Comp, levels = c("CHO","CHON","CHOS","CHONS"))
ft_black_vk=ft_black_vk[order(ft_black_vk$Inty,decreasing=F),]

ggplot(ft_black_vk,aes(x=Inty))+
  geom_histogram(binwidth = 100)+ylim(c(0,10))

ft_black_vk$sf=ifelse(ft_black_vk$Inty<10,"0 - 10",
                     ifelse(ft_black_vk$Inty<100,"10 - 100",
                            ifelse(ft_black_vk$Inty<200,"100 - 200",
                                   ifelse(ft_black_vk$Inty<400,"200 - 400",
                                          ifelse(ft_black_vk$Inty<800,"400 - 800","800 <")))))

#ft_black_vk$sf=factor(ft_black_vk$sf,levels = c("0 - 100","100 - 300","300 - 600","600 - 900","900 - 1200",
#                                              "1200 <"))


ggplot(ft_black_vk, aes(x=`O.C`, y=`H.C`, col=Comp))+
  geom_point(aes(size=sf), alpha=0.4)+
  facet_wrap(.~Sample, ncol = 5, scales = "fix")+
  #scale_size_continuous(range=c(1,30))+
  scale_size_manual(values=c(2,6,10,14,18,22))+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #scale_fill_manual(values = c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF","grey50"))+
  scale_color_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF"))+
  #scale_size_continuous(range = c(1,16))+
  #ggtitle("Beijing Event")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.text = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial",hjust=0.5, vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        #legend.position = c(0.18,0.14),
        #legend.position = "NULL",
        legend.position = "bottom",
        legend.direction = "vertical",
        #legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)),
  #       size=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7,alpha=1)))+
  ggsave(filename("vk_black_leg"),height = 18, width = 75, units = "cm", dpi = 700)



##vk all====

ft_red_vk$circle="Red"
ft_blue_vk$circle="Blue"
ft_black_vk$circle="Black"

ft_all_vk=rbind(ft_red_vk,ft_blue_vk,ft_black_vk)


ggplot(ft_all_vk, aes(x=`O.C`, y=`H.C`, col=Comp))+
  geom_point(aes(size=Inty), alpha=0.4)+
  facet_grid(circle~Sample, scales = "fix")+
  scale_size_continuous(range=c(1,30))+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.text = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        #legend.position = c(0.18,0.14),
        legend.position = "NULL",
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("vk_all"),height = 18, width = 75, units = "cm", dpi = 700)



ggplot(ft_vk, aes(x=`C.`, y=`Calc.m.z`, col=Comp))+
  geom_point(size=2, alpha=0.3)+
  facet_wrap(.~Group, ncol = 5, scales = "fix")+
  #scale_fill_manual(values = c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF","grey50"))+
  scale_color_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF"))+
  #scale_size_continuous(range = c(1,16))+
  #ggtitle("Beijing Event")+
  scale_y_continuous(name = "m/z")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.text = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "NULL",
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("C#mass"),height = 18, width = 75, units = "cm", dpi = 700)


###nmds with vector===========
envi_1st=fread("Datafile/envi_1st_sel.csv")
envi_1st

gnmds_1st$Sample=paste(gnmds_1st$Group,gnmds_1st$no, sep = "_")
gnmds_1st

gnmds_1st=gnmds_1st %>% inner_join(envi_1st)

gnmds_1st$id=ifelse(gnmds_1st$Group=="Beijing",5,
                    ifelse(gnmds_1st$Group=="Seosan",4,
                           ifelse(gnmds_1st$Group=="seoul",3,2)))
gnmds_1st=gnmds_1st[order(gnmds_1st$id)]
gnmds_1st

gnmds_1st$Event=factor(gnmds_1st$Event,levels = c("Event","Normal","Non-event"))
####chp=====
chp_m
gnmds_1st

chp_m$Sample=paste(chp_m$Group,chp_m$No, sep = "_")

vec_mds=gnmds_1st[,c(1,2)]

vec_chp=dcast(chp_m, Sample+Group+No~variable, sum, value.var = "val")

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

gnmds_1st$Group=factor(gnmds_1st$Group, levels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_npg()+
  scale_x_continuous(name = "NMDS1", breaks =seq(-2,2,0.5), limits = c(-1.0,2))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-1.5,1.5,0.5), limits = c(-1.2,1.2))+
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
  ggsave(filename("MDS_1st_withChp_nolabel"), height = 24, width = 20, units = "cm", dpi = 700)

####comp=====
ft_inty
gnmds_1st

ft_inty$Sample=paste(ft_inty$Group,ft_inty$no, sep = "_")

vec_mds=gnmds_1st[,c(1,2)]
vec_comp=dcast(ft_inty, Sample+Group+no~Comp, sum, value.var = "rel")

vec_comp$Group=factor(vec_comp$Group,levels = c("SUL","SS","B","NT","UL"),
                      labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))

vec_comp$Sample=paste(vec_comp$Group,vec_comp$no,sep = "_")

vec_comp=gnmds_1st %>% inner_join(vec_comp)
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


ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st, aes(x=MDS1, y=MDS2, fill=Group, shape=Event), col="black", size=9)+
  scale_fill_npg()+
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




####envi=====
envi_1st
gnmds_1st

gnmds_1st=gnmds_1st %>% inner_join(envi_1st)

vec_mds=gnmds_1st[,c(1,2)]

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
  scale_fill_npg()+
  scale_x_continuous(name = "NMDS1", breaks =seq(-2,2,0.5), limits = c(-1.0,2))+
  scale_y_continuous(name = "NMDS2", breaks =seq(-1.5,1.5,0.5), limits = c(-1.2,1.2))+
  scale_shape_manual(values = c(24,21,22))+
  #scale_shape_manual(values = c(21,22,23,24,25))+
  #geom_segment(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=0, y=0, xend=MDS1*1, yend=MDS2*1),
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

ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = gnmds_1st, aes(x=MDS1, y=MDS2, fill=Group, shape=Group), col="black", size=9)+
  scale_fill_npg()+
  #scale_fill_manual(values = c("#e97f02","#0080ff"))+
  scale_shape_manual(values = c(21,22,23,24,25))+
  geom_segment(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=0, y=0, xend=MDS1*0.8, yend=MDS2*0.8),
               arrow=arrow(length = unit(0.25, "cm")), size=1)+
  #geom_text(data = subset(arrow_4th_envi,arrow_4th_envi$p<0.05), aes(x=MDS1*0.88, y=MDS2*0.88, label=variable), size=8)+
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
  ggsave(filename("MDS_1st_withenvi_nolabel"), height = 20, width = 20, units = "cm", dpi = 700)

##spectrum====
FT_merge

FT_spec=melt(FT_merge[,c("Sample","Group","No","Calc.m.z","Bromo.Inty")],
             id.vars = c("Sample","Group","No","Calc.m.z")) %>% 
  dcast(Group~`Calc.m.z`, mean) %>% 
  melt(id.var="Group", na.rm = T) %>% `colnames<-`(c("Group","m/z","Relative intensity"))

FT_spec

FT_spec$`m/z`=as.numeric(FT_spec$`m/z`)

fwrite(FT_spec, file = "Datafile/mean_masslist.csv")
FT_spec=fread("Datafile/mean_masslist.csv")

FT_spec


FT_spec$Group=factor(FT_spec$Group,levels = c("SUL","SS","B","NT","UL"),
                   labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))

ggplot(FT_spec)+
  geom_segment(aes(x=`m/z`, xend=`m/z`, y=0.1, yend=`Relative intensity`), size=0.5)+
  facet_wrap(.~Group, scales = "free",nrow = 1)+
  #geom_line(aes(x=`m/z`, xend=`m/z`, y=`Inty Norm`,col=Assign), size=0.5)+
  #scale_color_manual(values = c("black","red"), labels=c("Assigned (1468)","Unassigned (933)"))+
  scale_x_continuous(expand = c(0.01,0.01), breaks = seq(250,1000,250), limits = c(150,990))+
  #scale_y_continuous(name = "Inty",expand = c(0.00,0.00), labels=scientific_10, limits = c(0.1,2890000000))+
  scale_y_continuous(name = "Relative intensity",expand = c(0.00,0.00))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial", hjust = 1),
        axis.title.y = element_text(size = 20, colour = "black",face = "bold", family = "Arial", margin = unit(c(0,0.5,0,0.2),"cm")),
        strip.background = element_blank(),
        strip.text = element_text(size = 18, colour = "black",face = "bold",family = "Arial"),
        legend.text = element_text(size = 20, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,0.2,0.2,0.2),"cm")),
        legend.title = element_text(size = 22, colour = "black", family = "Arial",face = "bold"),
        legend.box.background = element_blank(),
        legend.position = c(0.20, 0.84),
        legend.key.width = unit(2,"cm"),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = "", title.hjust = 0.0))+
  #coord_cartesian(ylim = c(0,100))+
  ggsave(filename("FT_mean_spec_full"), height = 15, width = 75, units = "cm", dpi = 700)


FT_spec$i2=ifelse(FT_spec$`Relative intensity`>100, 100,FT_spec$`Relative intensity`)

ggplot(FT_spec)+
  geom_segment(aes(x=`m/z`, xend=`m/z`, y=0.1, yend=i2), size=0.5)+
  facet_wrap(.~Group, scales = "free",nrow = 1)+
  #geom_line(aes(x=`m/z`, xend=`m/z`, y=`Inty Norm`,col=Assign), size=0.5)+
  #scale_color_manual(values = c("black","red"), labels=c("Assigned (1468)","Unassigned (933)"))+
  scale_x_continuous(expand = c(0.01,0.01), breaks = seq(250,1000,250), limits = c(150,990))+
  #scale_y_continuous(name = "Inty",expand = c(0.00,0.00), labels=scientific_10, limits = c(0.1,2890000000))+
  scale_y_continuous(name = "Relative intensity",expand = c(0.00,0.00))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial", hjust = 1),
        axis.title.y = element_text(size = 20, colour = "black",face = "bold", family = "Arial", margin = unit(c(0,0.5,0,0.2),"cm")),
        strip.background = element_blank(),
        strip.text = element_text(size = 18, colour = "black",face = "bold",family = "Arial"),
        legend.text = element_text(size = 20, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,0.2,0.2,0.2),"cm")),
        legend.title = element_text(size = 22, colour = "black", family = "Arial",face = "bold"),
        legend.box.background = element_blank(),
        legend.position = c(0.20, 0.84),
        legend.key.width = unit(2,"cm"),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = "", title.hjust = 0.0))+
  #coord_cartesian(ylim = c(0,100))+
  ggsave(filename("FT_mean_spec_zoom"), height = 15, width = 75, units = "cm", dpi = 700)




###294 vs PM conc====
tt=MDS_1st_in[,c("Sample","C10H17NO7S")]
tt=tt %>% tidyr::separate("Sample", c("Group","No"),sep="_")
tt$Group=factor(tt$Group,levels = c("SUL","SS","B","NT","UL"),
                     labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))

tt$Sample=paste(tt$Group,tt$No,sep = "_")
tt=tt[,c("Sample","Group","No","C10H17NO7S")]

tt=envi_1st %>% left_join(tt[,-c("Group","No")], by="Sample")

tt=tt %>% 
  mutate(date2 = as.POSIXct(tt$Date, format = '%Y-%m-%d %H'))

ggplot(tt, aes(x=PM2.5,y=C10H17NO7S))+
  geom_point()+
  facet_wrap(.~Group, scales="free")
lims <- as.POSIXct(strptime(c("2020-12-15 00:00","2021-01-16 0:00"), format = "%Y-%m-%d %H:%M"))

tt

tt$Group=factor(tt$Group, levels = c("Seoul","Seosan","Beijing","Noto",
                                     "Ulaanbaatar"))

tt$Pnorm=ifelse(tt$Group=="Beijing",tt$C10H17NO7S*0.05,
                ifelse(tt$Group=="Seosan",tt$C10H17NO7S*0.01,
                       ifelse(tt$Group=="Seoul",tt$C10H17NO7S*0.05,
                              ifelse(tt$Group=="Noto",tt$C10H17NO7S*0.03,tt$C10H17NO7S))))

tt_m=melt(tt[,c("Sample","Group","No","date2","Event","PM2.5","Pnorm")],
          id.vars = c("Sample","Group","No","date2","Event"), na.rm = T)

ggplot()+
  geom_line(data = tt_m ,aes(x=date2, y=value, col=variable),size=0.7, na.rm = T)+
  facet_wrap(.~Group, ncol = 1, scales = "free")+
  #scale_color_npg()+
  scale_x_datetime('',limits = lims,
                   date_breaks = '1 days',
                   date_labels = "%m/%d", expand = c(0.008,0.008))+
  #scale_y_continuous(name = expression(bold("Mass conc."~"("*"\u03bcg/"*m^"3"*")")),expand = c(0.02,0.02))+
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
        legend.position = "bottom")+
  xlab("")+
  #guides(fill=guide_legend(order = 1,title =bquote(PM[2.5]~"Episode"),col=c(NA,NA,NA), linetype=c(1,1,1), alpha=0.4),
  #       col=F)+
  ggsave(filename("294_comp"), height = 30, width = 40, units = "cm", dpi = 300)


