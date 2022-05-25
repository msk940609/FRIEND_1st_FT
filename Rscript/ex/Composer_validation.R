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

#####number of formula& intensity====

masslist=fread("Datafile/Masslist.csv")
masslist

masslist_m=melt(masslist[,c("Sample","Group","No","Masslist","Assigned formula","Ratio")],id.vars = c("Sample","Group","No"))

masslist_m

masslist_m$Group=factor(masslist_m$Group,levels = c("SUL","SS","B","NT","UL"),
                   labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))



ggplot()+
  stat_boxplot(data=masslist_m, aes(x=Group, y=value, fill=Group),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=masslist_m, aes(x=Group, y=value, fill=Group))+
  #geom_violin(data=masslist_m, aes(x=Group, y=value, fill=Group), trim=T)+
  #geom_boxplot(data=masslist_m, aes(x=Group, y=value, fill=Group),width=0.1,fill="white", outlier.color = "black")+
  #geom_text(data =value_max, aes(x=1, y=value, label=variable), col="white", size=0)+
  facet_wrap(.~variable, scales = "free", nrow=1, dir = "v")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  scale_fill_npg()+
  #geom_segment(data=sig_pos_m, aes(x=1, xend=2, y=ypos*1.01, yend=ypos*1.01))+
  #geom_segment(data=sig_pos, aes(x=Group, xend=Group, y=ypos.x*0.99, yend=ypos.y*1.01))+
  #geom_text(data =stat, aes(x=1.5, y=ifelse(stat$aster=="N.S",ifelse(stat$variable=="Mean DBE",11.5,ypos*1.07),ypos*1.02), label=aster), 
  #          col="black", size=ifelse(stat$aster=="N.S",6,8))+
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
  ggsave(filename("masslist_compare"),height = 15, width = 45, units = "cm", dpi = 300)

masslist_m$Sample=paste(masslist_m$Group,masslist_m$No,sep = "_")

pm2.5=fread(file = "Datafile/FRIENDs_envi.csv")

#scico_palette_show()

masslist_m=masslist_m %>% inner_join(pm2.5, by="Sample")

max(masslist_m$PM2.5)
library(scico)
masslist_m$exceed=ifelse(masslist_m$PM2.5>75,75,masslist_m$PM2.5)

ggplot()+
  stat_boxplot(data=masslist_m, aes(x=Group, y=value),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=masslist_m, aes(x=Group, y=value),fill="white", outlier.colour = NA)+
  geom_point(data=masslist_m, aes(x=Group, y=value, col=exceed), position=position_jitter(width = 0.3))+
  #geom_violin(data=masslist_m, aes(x=Group, y=value, fill=Group), trim=T)+
  #geom_boxplot(data=masslist_m, aes(x=Group, y=value, fill=Group),width=0.1,fill="white", outlier.color = "black")+
  #geom_text(data =value_max, aes(x=1, y=value, label=variable), col="white", size=0)+
  facet_wrap(.~variable, scales = "free", nrow=1, dir = "v")+
  #scico::scale_color_scico(palette = "lapaz", begin = 0.05, end = 0.95, limits=c(0,75), direction = -1)+
  scale_color_gradient2(low="orange",mid ="royalblue1", midpoint = 40 ,high="royalblue4", space ="Lab", limits=c(0,80))+
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
        legend.title = element_text(size = 18, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 14, colour = "black",margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        legend.key.width = unit(1.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.background = element_blank(),
        legend.direction = "vertical",
        legend.justification=c(0, 1),
        legend.position = c(0.008,0.42))+
  guides(color=guide_colorbar(title = expression(bold(PM[2.5])), order = 2,barheight = 7))+
  ggsave(filename("masslist_compare2"),height = 15, width = 45, units = "cm", dpi = 300)

masslist_md=dcast(masslist_m,Sample+Group+No+day+Event+PM2.5+exceed~variable, sum)
masslist_md

masslist_m=masslist_m %>% inner_join(masslist_md[,-c("Group","No","day","Event","PM2.5","exceed")], by="Sample")
masslist_m

pa=ggplot(masslist_m, aes(x=Group, y=value,`Mass list`=Masslist, `Assigned formula`=`Assigned formula`, `Assign ratio`=Ratio, PM2.5=PM2.5))+
  stat_boxplot(data=masslist_m, aes(x=Group, y=value),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=masslist_m, aes(x=Group, y=value),fill="white", outlier.colour = NA)+
  geom_point(data=masslist_m, aes(x=Group, y=value, col=exceed), position=position_jitter(width = 0.3))+
  #geom_violin(data=masslist_m, aes(x=Group, y=value, fill=Group), trim=T)+
  #geom_boxplot(data=masslist_m, aes(x=Group, y=value, fill=Group),width=0.1,fill="white", outlier.color = "black")+
  #geom_text(data =value_max, aes(x=1, y=value, label=variable), col="white", size=0)+
  facet_wrap(.~variable, scales = "free", nrow=1, dir = "v")+
  #scico::scale_color_scico(palette = "lapaz", begin = 0.05, end = 0.95, limits=c(0,75), direction = -1)+
  scale_color_gradient2(low="orange",mid ="royalblue1", midpoint = 40 ,high="royalblue4", space ="Lab", limits=c(0,80))+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
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
        legend.title = element_text(size = 18, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 14, colour = "black",margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        legend.key.width = unit(1.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.background = element_blank(),
        legend.direction = "vertical",
        legend.justification=c(0, 1),
        legend.position = c(0.008,0.42))+
  guides(color=guide_colorbar(title = expression(bold(PM[2.5])), order = 2,barheight = 7))

pa

ggplotly(pa,tooltip = c("x", "y","Mass list", "Assigned formula","Assign ratio","PM2.5")) %>% 
  layout(legend = list(orientation = "h", x = 0.3, y = -0.2))

pa=ggplot(masslist_m)+
 # stat_boxplot(geom='errorbar', linetype=1, width=0.25,position = position_dodge(width = 0.75))+
  geom_point(aes(x=Group, y=value,col=exceed,Sample=Sample,`Mass list`=Masslist, `Assigned formula`=`Assigned formula`, `Assign ratio`=Ratio, PM2.5=PM2.5),position=position_jitter(width = 0.3))+
  geom_boxplot(aes(x=Group, y=value),fill="white", outlier.color =  NA, outlier.size = 0, outlier.shape = NA)+
  facet_wrap(.~variable, scales = "free", nrow=1, dir = "v")+
  scale_color_gradient2(low="orange",mid ="royalblue1", midpoint = 40 ,high="royalblue4", space ="Lab", limits=c(0,80))+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        #axis.line.y.right = element_line(size = 1, color = "black"),
        #axis.line.x.top =  element_line(colour = "black"),
        #axis.text.x = element_text(size = 16,colour = "black", angle = 0, hjust = 0.5,vjust = 1.0),
        #axis.ticks.x = element_blank(),
        #axis.ticks.length.y = unit(0.15,"cm"),
        #axis.text.y = element_text(size = 16, colour = "black" ),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        axis.title.y.left = element_text(size = 18, colour = "black",margin = unit(c(0.0,0.0,0.0,0.0),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.4,0.2,0.2,0.2),"cm")),
        #strip.text.y = element_blank(),
        #legend.margin = unit(c(0.0,0.0,0.0,0.0),"cm"),
        #legend.title = element_text(size = 18, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        #legend.text = element_text(size = 14, colour = "black",margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        #legend.key.width = unit(1.0,"cm"),
        #legend.key.height = unit(1.0,"cm"),
        #legend.background = element_blank(),
        #legend.direction = "vertical",
        #legend.justification=c(0, 1),
        )+
  guides(color=guide_colorbar(title = expression(bold(PM[2.5])), order = 2,barheight = 7))

pa

ggplotly(pa,tooltip = c("x","Sample","Mass list", "Assigned formula","Assign ratio","PM2.5")) %>% 
  layout(legend = list(orientation = "v", x = 0.3, y = -0.2))


iris%>%
  group_by(Species) %>%
  do(p=plot_ly(., x = ~Sepal.Length, y = ~Sepal.Width, color = ~Species, type = "scatter")) %>%
  subplot(nrows = 1, shareX = TRUE, shareY = TRUE)



