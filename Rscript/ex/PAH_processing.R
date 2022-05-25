###start load library (MS default)=====
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

pah_b=fread("Beijing_merge.csv")
pah_sul=fread("Seoul_merge.csv")
pah_ss=fread("Seosan_merge.csv")
pah_ul=fread("Ulaanbaatar_merge.csv")

pah_b
pah_sul
pah_ss
pah_ul

pah_merge=rbind(pah_b,pah_sul,pah_ss,pah_ul)

pah_merge=pah_merge[,c("Peak #","Name","Expected Ion m/z","Formula","Area","Normarea","Comp",
             "#C","#H","#N","#O","#S","Sample")]

pah_merge=pah_merge %>% separate("Sample",c("Group","Date","Event"), sep = "_")
fwrite(pah_merge,"Datafile/pah_merge.csv")

pah_merge=fread("Datafile/pah_merge.csv")

pah_merge$`C#`=as.numeric(pah_merge$`C#`)
pah_merge$`H#`=as.numeric(pah_merge$`H#`)
pah_merge$`N#`=as.numeric(pah_merge$`N#`)
pah_merge$`O#`=as.numeric(pah_merge$`O#`)
pah_merge$`S#`=as.numeric(pah_merge$`S#`)

pah_merge$DBE=pah_merge$`C#`+1-pah_merge$`H#`/2-pah_merge$`N#`/2
pah_merge$`O/C`=pah_merge$`O#`/pah_merge$`C#`
pah_merge$`H/C`=pah_merge$`H#`/pah_merge$`C#`

pah_merge$CAI=pah_merge$`C#`-pah_merge$`N#`-pah_merge$`O#`/2-pah_merge$`S#`

pah_merge$DBEAI=1+pah_merge$`C#`-pah_merge$`O#`/2-pah_merge$`S#`-pah_merge$`H#`/2-pah_merge$`N#`/2

pah_merge$AI=ifelse(pah_merge$CAI<=0,0,ifelse(pah_merge$DBEAI<0,0,pah_merge$DBEAI/pah_merge$CAI))

pah_merge


pah_merge=fread("Datafile/pah_merge.csv")
pah_merge$Sample=paste(pah_merge$Group,pah_merge$Date,pah_merge$Event,sep = "_")

pah_m=melt(pah_merge[,c("Group","Event","Formula","Normarea")], id.vars = c("Group","Event","Formula")) %>% 
  dcast(Group+Event~Formula, mean) %>% 
  melt(id.vars=c("Group","Event"), na.rm = T) %>% `colnames<-`(c("Group","Event","Formula","Norm Area"))

pah_m

table(pah_m$Group)
pah_fm=pah_merge[,c("Formula","Comp","C#","H#","N#","O#","S#","O/C","H/C","DBE","AI")]
pah_fm=unique(pah_fm)

pah_m

insert_minor <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_1 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_2 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}


pah_m$Group=factor(pah_m$Group,levels = c("Seoul","Seosan","Beijing","Ulaanbaatar"),
                          labels = c("Seoul","Seosan","Beijing","Ulaanbaatar"))


pah_m=pah_m %>%inner_join(pah_fm)

table(pah_m$Group)

pah_m$id=ifelse(pah_m$Group=="Seoul","Seoul (141)",
                ifelse(pah_m$Group=="Seosan","Seosan (79)",
                       ifelse(pah_m$Group=="Beijing","Beijing (218)","Ulaanbaatar (228)")))

pah_m$id=factor(pah_m$id,levels = c("Seoul (141)","Seosan (79)","Beijing (218)","Ulaanbaatar (228)"),
                labels = c("Seoul (141)","Seosan (79)","Beijing (218)","Ulaanbaatar (228)"))

pah_m$Comp2=factor(pah_m$Comp,levels = c("CH","CHN","CHS","CHO","CHNO","CHOS"),
                  labels = c("CH","CHN","CHS","CHO","CHON","CHOS"))

pah_m=pah_m[order(pah_m$Comp2)]

pah_m$id2=paste(pah_m$Group,pah_m$Event, sep = "_")
table(pah_m$id2)
table(pah_m$Group)


mypal=c("#759BB2FF","#F17F42","grey45","#BC3C29FF","#EFC000FF","#008B45FF")
ggplot(pah_m, aes(x=`O/C`, y=`H/C`, col=Comp2))+
  geom_point(aes(size=`Norm Area`), alpha=0.4)+
  facet_grid(Event~Group, scales = "free")+
  #facet_wrap(.~id, scales = "free", ncol = 4)+
  scale_size_continuous(range=c(2,20))+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,0.35), breaks = round(seq(0,2.1,0.1),1))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01),limits = c(-0.01,1.40), breaks = round(seq(0,2.1,0.4),1))+
  #scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,0.4), breaks = round(seq(0,2.1,0.1),1),
  #                   labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  #scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
  #                   labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_manual(values = mypal)+
  #scale_color_igv()+
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
        legend.text = element_text(debug = F,size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.2,0,0,0),"cm")),
        legend.spacing = unit(0.2,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "bottom",
        legend.margin = margin(0.2,0.2,0.2,0.2, "cm"),
        legend.text.align=0,
        #legend.position = c(0.12,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(order=1,title.vjust = 0.8,title = "",ncol = 6,title.position = "left",byrow = T,override.aes = list(shape=21,size=7,fill=mypal,alpha=1)),
         size=guide_legend("Norm. area", byrow=T, ncol=3, title.position = "left"))+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("gc_vk_leg"),height = 30, width = 60, units = "cm", dpi = 700)
  #ggsave(filename("gc_vk_all"),height = 15, width = 60, units = "cm", dpi = 700)

pah_m$id=paste(pah_m$Group,pah_m$Event,sep = "_")
table(pah_m$id)

mypal=c("#759BB2FF","#F17F42","grey45","#BC3C29FF","#EFC000FF","#008B45FF")
ggplot(pah_m, aes(x=`C#`, y=`DBE`, col=Comp2))+
  geom_point(aes(size=`Norm Area`), alpha=0.4)+
  facet_grid(Event~Group, scales = "free")+
  #facet_wrap(.~id, scales = "free", ncol = 4)+
  scale_size_continuous(range=c(2,20), labels = seq(2,20,2), breaks = seq(2,20,2))+
  #scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,0.35), breaks = round(seq(0,2.1,0.1),1))+
  #scale_y_continuous(name = "H/C",expand = c(0.01,0.01),limits = c(-0.01,1.40), breaks = round(seq(0,2.1,0.4),1))+
  scale_color_manual(values = mypal)+
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
        legend.spacing = unit(-0.2,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "NULL",
        #legend.position = c(0.12,0.16),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = "",ncol = 3,byrow = T,override.aes = list(shape=21,size=7,fill=mypal,alpha=1)), size="none")+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("gc_dbe_all"),height = 30, width = 60, units = "cm", dpi = 700)


pah_m


table(pah_m$Comp2)
pah_merge

####Chemical composition======
pah_inty_all=as.data.table(aggregate(pah_merge$Normarea, by=list(`Sample`=pah_merge$Sample), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Tot"))

pah_inty=as.data.table(aggregate(pah_merge$Normarea, by=list(`Sample`=pah_merge$Sample,`Group`=pah_merge$Group,Event=pah_merge$Event,Comp=pah_merge$Comp), FUN=sum))
pah_inty

pah_inty=pah_inty %>% inner_join(pah_inty_all,by = "Sample")
pah_inty$rel=round(pah_inty$x/pah_inty$Tot*100,1)

table(pah_inty$Comp)



pah_inty$Comp2=factor(pah_inty$Comp,levels = c("CH","CHN","CHS","CHO","CHNO","CHOS"),
                   labels = c("CH","CHN","CHS","CHO","CHON","CHOS"))


pah_inty$Group=factor(pah_inty$Group,levels = c("Seoul","Seosan","Beijing","Ulaanbaatar"),
                     labels = c("Seoul","Seosan","Beijing","Ulaanbaatar"))
pah_inty


#comp_max=as.data.table(aggregate(ft_inty$rel, by=list(`Comp`=ft_inty$Comp), FUN=max))
#comp_max$x=comp_max$x*1.2

npg=pal_npg("nrc")(5)
npg_no_jp=npg[-4]

ggplot()+
  stat_boxplot(data=pah_inty, aes(x=Group, y=rel, fill=Event),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=pah_inty, aes(x=Group, y=rel, fill=Event),alpha=1, outlier.color = NA)+
  #geom_text(data =comp_max, aes(x=1, y=x, label=Comp), col="white", size=0)+
  #geom_text(data =stat, aes(x=Group, y=rel, label=labels), 
  #          col="black", size=6)+
  #facet_grid(Event~Comp2, scales = "free")+
  geom_vline(xintercept = c(1.5,2.5,3.5), lty=2)+
  facet_wrap(.~Comp2, scales = "free_y", nrow=2, dir = "h")+
  scale_x_discrete(name="")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  #scale_fill_manual(values = npg_no_jp)+
  scale_fill_manual(values = c("grey70", "#FFFFFF"))+
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
        legend.background = element_rect(fill = "white"),
        legend.direction = "vertical",
        legend.justification=c(0.5, 0.5),
        legend.position = c(0.06,0.052))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("pah_Comp_compare"),height = 30, width = 45, units = "cm", dpi = 300)


####Combustion prop======
pah_merge$`DBE/C`=pah_merge$DBE/pah_merge$`C#`
pah_merge$Comb=ifelse(pah_merge$`DBE/C`>0.7, "Combustion derived","organic")

pah_merge$Comb_inty=ifelse(pah_merge$Comb=="organic",0,pah_merge$Normarea)


table(pah_merge$Comp)

pah_inty_all

comb_inty=as.data.table(aggregate(pah_merge$Normarea, by=list(`Sample`=pah_merge$Sample,`Group`=pah_merge$Group,Event=pah_merge$Event,Comb=pah_merge$Comb), FUN=sum))
comb_inty

comb_inty=comb_inty %>% inner_join(pah_inty_all,by = "Sample")

comb_inty$rel=round(comb_inty$x/comb_inty$Tot*100,1)

comb_sel=subset(comb_inty,comb_inty$Comb=="Combustion derived")

comb_sel$Group=factor(comb_sel$Group,levels = c("Seoul","Seosan","Beijing","Ulaanbaatar"),
                     labels = c("Seoul","Seosan","Beijing","Ulaanbaatar"))


ggplot()+
  stat_boxplot(data=comb_sel, aes(x=Group, y=rel, fill=Event),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=comb_sel, aes(x=Group, y=rel, fill=Event),alpha=1, outlier.color = NA)+
  #geom_text(data =comp_max, aes(x=1, y=x, label=Comp), col="white", size=0)+
  #geom_text(data =stat, aes(x=Group, y=rel, label=labels), 
  #          col="black", size=6)+
  #facet_grid(Event~Comp2, scales = "free")+
  #facet_wrap(.~Comp2, scales = "free_y", nrow=2, dir = "h")+
  scale_x_discrete(name="")+
  geom_vline(xintercept = c(1.5,2.5,3.5), lty=2)+
  scale_y_continuous(name = "Proportion of combustion derived compounds (%)",
                     expand = c(0.02,0.02))+
  #scale_fill_manual(values = npg_no_jp)+
  scale_fill_manual(values = c("grey70", "#FFFFFF"))+
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
        legend.background = element_rect(fill = "white"),
        legend.direction = "vertical",
        legend.justification=c(0, 1),
        legend.position = c(0.02,0.15))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("pah_combustion_compare"),height = 20, width = 20, units = "cm", dpi = 300)

####Chemical prop======

pah_chp=melt(pah_merge[,c("Sample","Group","Event","AI","DBE","H/C","O/C")], id.vars=c("Sample","Group","Event")) %>% 
  dcast(Sample+Group+Event~variable, mean) %>% 
  melt(id.var=c("Sample","Group","Event"), na.rm = T)


pah_chp$Group=factor(pah_chp$Group,levels = c("Seoul","Seosan","Beijing","Ulaanbaatar"),
                      labels = c("Seoul","Seosan","Beijing","Ulaanbaatar"))

pah_chp$variable=factor(pah_chp$variable,levels = c("AI","DBE","H/C","O/C"),
                     labels = c("Mean AI","Mean DBE","Mean H/C","Mean O/C"))


ggplot(data=pah_chp, aes(x=Group,y=value,fill=Event))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(alpha=1, outlier.color = NA)+
  #geom_point(position = position_jitterdodge(jitter.width = 0.1))+
  #facet_grid(Event~Comp2, scales = "free")+
  facet_wrap(.~variable, scales = "free_y", nrow=2, dir = "h")+
  scale_x_discrete(name="")+
  geom_vline(xintercept = c(1.5,2.5,3.5), lty=2)+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  scale_fill_manual(values = c("grey70", "#FFFFFF"))+
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
        legend.background = element_rect(fill = "white"),
        legend.direction = "vertical",
        legend.justification=c(0, 1),
        legend.position = c(0.02,0.6275))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("pah_Chp_compare"),height = 30, width = 30, units = "cm", dpi = 300)




pah_chp_all=melt(pah_merge[,c("Sample","Group","Event","AI","DBE","H/C","O/C")], id.vars=c("Sample","Group","Event"))
pah_chp_all


ggplot(data=pah_chp_all, aes(x=Group,y=value,fill=Event))+
  #stat_boxplot(geom='errorbar', linetype=1, width=0.25,
  #             position = position_dodge(width = 0.75))+
  geom_violin()+
  #facet_grid(Event~Comp2, scales = "free")+
  facet_wrap(.~variable, scales = "free_y", nrow=2, dir = "h")+
  scale_x_discrete(name="")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  scale_fill_manual(values = c("orangered2","royalblue2"))+
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
        legend.position = c(0.02,0.63))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("pah_Chp_all_compare"),height = 30, width = 30, units = "cm", dpi = 300)


ggplot(data=pah_merge, aes(x=DBE))+
  geom_histogram(binwidth = 1)+
  facet_grid(Event~Group, scales = "free")

  
