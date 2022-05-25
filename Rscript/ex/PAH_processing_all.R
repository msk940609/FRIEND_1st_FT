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

source("Rscript/func_generate_lable.R")


#0.data load=====
#pah_b=fread("Beijing_merge.csv")
#pah_sul=fread("Seoul_merge.csv")
#pah_ss=fread("Seosan_merge.csv")
#pah_ul=fread("Ulaanbaatar_merge.csv")

#pah_merge=rbind(pah_b,pah_sul,pah_ss,pah_ul)

#pah_merge=pah_merge[,c("Peak #","Name","Expected Ion m/z","Formula","Area","Normarea","Comp",
#             "#C","#H","#N","#O","#S","Sample")]

#pah_merge=pah_merge %>% separate("Sample",c("Group","Date","Event"), sep = "_")
#fwrite(pah_merge,"Datafile/pah_merge.csv")

pah_merge=fread("Datafile/pah_merge.csv")

pah_merge$`C#`=as.numeric(pah_merge$`C#`)
pah_merge$`H#`=as.numeric(pah_merge$`H#`)
pah_merge$`N#`=as.numeric(pah_merge$`N#`)
pah_merge$`O#`=as.numeric(pah_merge$`O#`)
pah_merge$`S#`=as.numeric(pah_merge$`S#`)

pah_merge$DBE=pah_merge$`C#`+1-pah_merge$`H#`/2-pah_merge$`N#`/2-pah_merge$`O#`-pah_merge$`S#`
pah_merge$`O/C`=pah_merge$`O#`/pah_merge$`C#`
pah_merge$`H/C`=pah_merge$`H#`/pah_merge$`C#`

pah_merge$CAI=pah_merge$`C#`-pah_merge$`N#`-pah_merge$`O#`/2-pah_merge$`S#`

pah_merge$DBEAI=1+pah_merge$`C#`-pah_merge$`O#`/2-pah_merge$`S#`-pah_merge$`H#`/2-pah_merge$`N#`/2

pah_merge$AI=ifelse(pah_merge$CAI<=0,0,ifelse(pah_merge$DBEAI<0,0,pah_merge$DBEAI/pah_merge$CAI))

pah_merge
pah_merge$Sample=paste(pah_merge$Group,pah_merge$Date,pah_merge$Event,sep = "_")

#1.data processing======
pah_merge$Freq=1
pah_merge

pah_merge

pah_name_iso=melt(pah_merge[,c("Group","Formula","Name","Freq")], id.vars = c("Group","Name","Formula")) %>% 
  dcast(Group+Formula~Name, sum) %>% 
  melt(id.vars=c("Group","Formula")) %>% 
  dcast(Group~Formula, sum) %>% 
  melt(id.vars=c("Group")) %>% `colnames<-`(c("Group","Formula","Isomer"))


pah_name_iso
pah_name_iso_sel=subset(pah_name_iso,pah_name_iso$Isomer>0)  
pah_name_iso_sel


pah_m=melt(pah_merge[,c("Sample","Group","Formula","Normarea")], id.vars = c("Sample","Group","Formula")) %>% 
  #dcast(Sample+Group~Formula, sum) %>% 
  #melt(id.vars=c("Sample","Group")) %>% 
  dcast(Group~Formula, mean) %>% 
  melt(id.vars=c("Group"),na.rm = T) %>% `colnames<-`(c("Group","Formula","Normarea"))

dim(pah_m)
dim(pah_iso)
pah_m_iso=subset(pah_iso,pah_iso$Isomer>0)
unique()
#pah_m$Isomer=round(pah_m$Isomer,0)
max(pah_m$Normarea)

pah_m_sel=subset(pah_m,pah_m$Normarea!=0) %>% droplevels()
#pah_m_sel$Isomer=round(pah_m_sel$Isomer,0)
#pah_m_sel=subset(pah_m_sel,pah_m_sel$Isomer!=0) %>% droplevels()
pah_m_sel

table(pah_m_sel$Group)

pah_m_sel$Group=factor(pah_m_sel$Group,levels = c("Seoul","Seosan","Beijing","Ulaanbaatar"))

pah_fm=pah_merge[,c("Formula","Comp","C#","H#","N#","O#","S#","O/C","H/C","DBE","AI")]
pah_fm=unique(pah_fm)

pah_m_sel=pah_m_sel %>%inner_join(pah_fm)
#pah_m$id=ifelse(pah_m$Group=="Seoul","Seoul (141)",
#                ifelse(pah_m$Group=="Seosan","Seosan (79)",
#                       ifelse(pah_m$Group=="Beijing","Beijing (218)","Ulaanbaatar (228)")))

#pah_m$id=factor(pah_m$id,levels = c("Seoul (141)","Seosan (79)","Beijing (218)","Ulaanbaatar (228)"),
#                labels = c("Seoul (141)","Seosan (79)","Beijing (218)","Ulaanbaatar (228)"))

pah_m_sel$Comp2=factor(pah_m_sel$Comp,levels = c("CH","CHN","CHS","CHO","CHNO","CHOS"),
                  labels = c("CH","CHN","CHS","CHO","CHON","CHOS"))

pah_m_sel=pah_m_sel[order(pah_m_sel$Comp2)]

table(pah_m_sel$Group)

##1-1)vk plot=====
insert_minor <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_1 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_2 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}
?scale_size_area
mypal=c("#759BB2FF","#F17F42","grey45","#BC3C29FF","#EFC000FF","#008B45FF")

dim(pah_m_sel)
dim(pah_iso)


pah_m_sel=pah_m_sel %>% inner_join(pah_name_iso_sel, by=c("Group","Formula")) 

ggplot(pah_m_sel, aes(x=`O/C`, y=`H/C`, col=Comp2))+
  geom_point(aes(size=`Normarea`), alpha=0.4)+
  facet_grid(.~Group, scales = "free")+
  #facet_wrap(.~id, scales = "free", ncol = 4)+
  scale_size_continuous(range = c(1.5, 25))+
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
        legend.text = element_text(debug = F,size = 22, colour = "black", family = "Arial", vjust = 1.2,margin = unit(c(0.2,0,0,0),"cm")),
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
         size=guide_legend("Norm. area", byrow=T, ncol=5, title.position = "left"))+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("gc_vk_leg"),height = 20, width = 60, units = "cm", dpi = 300)


pah_m$id=paste(pah_m$Group,pah_m$Event,sep = "_")
table(pah_m$id)
##1-2)C#dbeplot=====
pah_m_sel=fread("Datafile/PAH_melt.csv")

pah_m_sel$Group=factor(pah_m_sel$Group,levels = c("Seoul","Seosan","Beijing","Ulaanbaatar"))
pah_m_sel
max(pah_m_sel$Isomer)
max(pah_m_sel$Isomer)
mypal=c("#759BB2FF","#F17F42","grey45","#BC3C29FF","#EFC000FF","#008B45FF")

ggplot(pah_m_sel,aes(x=Isomer))+
  geom_histogram(bins=20)+
  facet_wrap(.~Group)


table(pah_m_sel$Group)

#pah_m_sel$sf=ifelse(pah_m_sel$Isomer<4,"1 - 4",
#                    ifelse(pah_m_sel$Isomer<8,"5 - 8",
#                           ifelse(pah_m_sel$Isomer<12,"9 - 12",
#                                  ifelse(pah_m_sel$Isomer<16,"13 - 16","16 <"))))
pah_m_sel$sf=ifelse(pah_m_sel$Isomer<4,"1 - 4",
                    ifelse(pah_m_sel$Isomer<8,"5 - 8",
                           ifelse(pah_m_sel$Isomer<12,"9 - 12","12 <")))


pah_m_sel$sf=factor(pah_m_sel$sf,levels = c("1 - 4","5 - 8","9 - 12","12 <"))


fwrite(pah_m_sel,file = "pah_sel.csv")

pah_m_sel$Comp2=factor(pah_m_sel$Comp,levels = c("CH","CHN","CHS","CHO","CHNO","CHOS"),
                       labels = c("CH","CHN","CHS","CHO","CHON","CHOS"))

ggplot(pah_m_sel, aes(x=`C#`, y=`DBE`, col=Comp2))+
  geom_point(aes(size=`Normarea`), alpha=0.4)+
  facet_grid(.~Group, scales = "free")+
  scale_size_continuous(range = c(3, 25))+
  #scale_size_manual(values=c(4,10,16,22))+
  #scale_size_continuous(range=c(2,20), labels = seq(2,20,2), breaks = seq(2,20,2))+
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
        legend.text = element_text(size = 22, colour = "black", family = "Arial",vjust = 0.5,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(-0.5,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.050,0.78),
        #legend.position = "NULL",
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        #legend.position = c(0.12,0.16),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  #guides(col=guide_legend(title = "",ncol = 3,byrow = T,override.aes = list(shape=21,size=7,fill=mypal,alpha=1)), size="none")+
  guides(col="none",size=guide_legend(title = "Norm. area", title.hjust = 1.2))+
  ggsave(filename("gc_dbe_all_no_leg"),height = 18, width = 60, units = "cm", dpi = 300)

pah_m

table(pah_m$Comp2)
pah_merge

#3.Chemical composition======
#pah_inty_all=as.data.table(aggregate(pah_merge$Normarea, by=list(`Sample`=pah_merge$Sample), FUN=sum)) %>% 
#  `colnames<-`(c("Sample","Tot"))

#pah_inty=as.data.table(aggregate(pah_merge$Normarea, by=list(`Sample`=pah_merge$Sample,`Group`=pah_merge$Group,Event=pah_merge$Event,Comp=pah_merge$Comp), FUN=sum))
#pah_inty
#pah_inty=pah_inty %>% inner_join(pah_inty_all,by = "Sample")
#pah_inty$rel=round(pah_inty$x/pah_inty$Tot*100,1)

pah_merge

pah_inty_all=melt(pah_merge[,c("Sample","Group","Comp","Normarea")],id.vars = c("Sample","Group","Comp")) %>% 
  dcast(Sample+Group~Comp, sum)

pah_inty_all$tot=rowSums(pah_inty_all[,-c("Sample","Group")])

pah_inty_all



pah_inty=melt(pah_merge[,c("Sample","Group","Comp","Normarea")],id.vars = c("Sample","Group","Comp")) %>% 
  dcast(Sample+Group~Comp, sum) %>% 
  melt(id.vars=c("Sample","Group"))


pah_inty=pah_inty %>% inner_join(pah_inty_all[,c("Sample","tot")])
pah_inty

pah_inty=pah_inty %>% inner_join(pah_inty_all,by = "Sample")
pah_inty$rel=round(pah_inty$value/pah_inty$tot*100,1)
pah_inty=pah_inty %>% `colnames<-`(c("Sample","Group","Comp","value","tot","rel"))


table(pah_inty$Comp)

pah_inty$Comp2=factor(pah_inty$Comp,levels = c("CH","CHN","CHS","CHO","CHNO","CHOS"),
                   labels = c("CH","CHN","CHS","CHO","CHON","CHOS"))

pah_inty$Group=factor(pah_inty$Group,levels = c("Seoul","Seosan","Beijing","Ulaanbaatar"))
pah_inty

#comp_max=as.data.table(aggregate(ft_inty$rel, by=list(`Comp`=ft_inty$Comp), FUN=max))
#comp_max$x=comp_max$x*1.2

pah_inty

##3-1)normality, equal variace test====
norm=data.table()
comp=unique(pah_inty$Comp)

for (i in 1:length(comp)) {
  temp=subset(pah_inty,pah_inty$Comp==comp[i])
  evtmp=levene.test(temp$rel,temp$Group)
  grp=unique(temp$Group)
  for (j in 1:length(grp)) {
    tmp2=subset(temp,temp$Group==grp[j])
    
    if(sum(tmp2$rel)==0) {
      new=data.table("Comp"=comp[i],"Group"=grp[j],"p"=1, "EV"=round(evtmp$p.value,3))
    }else{
      stmp=shapiro.test(tmp2$rel)
      new=data.table("Comp"=comp[i],"Group"=grp[j],"p"=round(stmp$p.value,3), "EV"=round(evtmp$p.value,3))
      
    }
    norm=rbind(norm,new)
  }
}

norm


##3-2)statistical analysis (parametric)====
stat_par=data.table()
for (i in 1:length(comp)) {
#  i=1
  temp=subset(pah_inty,pah_inty$Comp==comp[i])
  
  tm=as.data.table(aggregate(temp$rel, by=list(`variable`=temp$Group), FUN=median))
  tm=tm[order(tm$x, decreasing = T),]
  
  temp$Group=factor(temp$Group,levels = tm$variable)
  #ans <- posthoc.kruskal.dunn.test(rel ~ Group,data=temp ,p.adjust="none") ##
  #ans.p <- get.pvalues(ans)
  #ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  ans=dunnettT3Test(rel ~ Group,data=temp, p.adjust.methods="fdr")
  ans.p <- get.pvalues(ans)
  ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  
  new=generate_label(temp, flev = "Group",value.var = "rel",offset = 4, cld_info = ans.mcV)
  new$Comp=comp[i]
  stat_par=rbind.data.frame(stat_par, new)
  
}


stat_par=stat_par %>% `colnames<-`(c("Group","labels","rel","Comp"))
stat_par


##3-3)statistical analysis (non parametric)====
stat_npr=data.table()
for (i in 1:length(comp)) {
  temp=subset(pah_inty,pah_inty$Comp==comp[i])
  
  tm=as.data.table(aggregate(temp$rel, by=list(`variable`=temp$Group), FUN=median))
  tm=tm[order(tm$x, decreasing = T),]
  
  temp$Group=factor(temp$Group,levels = tm$variable)
  ans <- posthoc.kruskal.dunn.test(rel ~ Group,data=temp ,p.adjust="bonferroni") ##
  ans.p <- get.pvalues(ans)
  ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  
  new=generate_label(temp, flev = "Group",value.var = "rel",offset = 0, cld_info = ans.mcV)
  new$Comp=comp[i]
  stat_npr=rbind.data.frame(stat_npr, new)
  
}

stat_npr=stat_npr %>% `colnames<-`(c("Group","labels","rel","Comp"))
stat_npr$Comp2=factor(stat_npr$Comp,levels = c("CH","CHN","CHS","CHO","CHNO","CHOS"),
                      labels = c("CH","CHN","CHS","CHO","CHON","CHOS"))

stat_npr

norm
#stat_par_sel=stat_par %>% filter(Comp%in%c("CHO","CHON"))
#stat_npr_sel=stat_npr %>% filter(Comp%in%c("CHOS","CHONS"))
stat_npr
stat=stat_npr

comp_max=as.data.table(aggregate(pah_inty$rel, by=list(`Comp`=pah_inty$Comp), FUN=max))
comp_max$x=comp_max$x*1.2
pah_inty_sel=subset(pah_inty,pah_inty$rel>0)
##3-4)plotting=====
ggplot()+
  stat_boxplot(data=pah_inty_sel, aes(x=Group, y=rel, fill=Group),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=pah_inty_sel, aes(x=Group, y=rel, fill=Group),alpha=1, outlier.color = NA)+
  #geom_text(data =comp_max, aes(x=1, y=rel, label=Comp), col="white", size=0)+
  geom_text(data =stat, aes(x=Group, y=ifelse(Comp2=="CH",rel+1,
                                              ifelse(Comp2=="CHN",rel+0.5,
                                                     ifelse(Comp2=="CHS",rel+0.25,
                                                            ifelse(Comp2=="CHO",rel+1,
                                                                   ifelse(Comp2=="CHON",rel+0.08,rel+0.15))))), label=labels), 
            col="black", size=6)+
  #facet_grid(Event~Comp2, scales = "free")+
  facet_wrap(.~Comp2, scales = "free_y", nrow=2, dir = "h")+
  scale_x_discrete(name="")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  #scale_fill_manual(values = npg_no_jp)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#3C5488FF"))+
#  scale_fill_manual(values = c("orangered2","royalblue2"))+
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
        #legend.position = c(0.02,0.1)
        legend.position = "Null"
        )+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("pah_Comp_compare"),height = 30, width = 45, units = "cm", dpi = 300)


#4.Chemical properties======
pah_chp=melt(pah_merge[,c("Sample","Group","AI","DBE","H/C","O/C")], id.vars=c("Sample","Group")) %>% 
  dcast(Sample+Group~variable, mean) %>% 
  melt(id.var=c("Sample","Group"), na.rm = T)

pah_chp$Group=factor(pah_chp$Group,levels = c("Seoul","Seosan","Beijing","Ulaanbaatar"))

pah_chp$variable=factor(pah_chp$variable,levels = c("AI","DBE","H/C","O/C"),
                     labels = c("Mean AI","Mean DBE","Mean H/C","Mean O/C"))


##4-1)normality, equal variace test====
norm=data.table()
comp=unique(pah_chp$variable)

for (i in 1:length(comp)) {
  temp=subset(pah_chp,pah_chp$variable==comp[i])
  evtmp=levene.test(temp$value,temp$Group)
  grp=unique(temp$Group)
  for (j in 1:length(grp)) {
    tmp2=subset(temp,temp$Group==grp[j])
    
    if(sum(tmp2$value)==0) {
      new=data.table("Comp"=comp[i],"Group"=grp[j],"p"=1, "EV"=round(evtmp$p.value,3))
    }else{
      stmp=shapiro.test(tmp2$value)
      new=data.table("Comp"=comp[i],"Group"=grp[j],"p"=round(stmp$p.value,3), "EV"=round(evtmp$p.value,3))
      
    }
    norm=rbind(norm,new)
  }
}

norm
##para: mean AI, DBE, H/C, O/C -> boneferroi, 
##unequal: O/C -> dunnett3

##4-2)statistical analysis (parametric)====
stat_par=data.table()
for (i in 1:length(comp)) {
  temp=subset(pah_chp,pah_chp$variable==comp[i])
  
  tm=as.data.table(aggregate(temp$value, by=list(`variable`=temp$Group), FUN=median))
  tm=tm[order(tm$x, decreasing = T),]
  
  temp$Group=factor(temp$Group,levels = tm$variable)
  #ans <- lsdTest(value ~ Group,data=temp ,p.adjust="fdr") ##
  #ans <- scheffeTest(value ~ Group,data=temp ,p.adjust="fdr") ##
  #ans.p <- get.pvalues(ans)
  #ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  ans <- dunnettT3Test(value ~ Group,data=temp ,p.adjust="fdr") ##unequal variance, normal distribution
  ans.p <- get.pvalues(ans)
  ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  
  new=generate_label(temp, flev = "Group",value.var = "value",offset = 0, cld_info = ans.mcV)
  new$Comp=comp[i]
  stat_par=rbind.data.frame(stat_par, new)
  
}


stat_par=stat_par %>% `colnames<-`(c("Group","labels","value","variable"))
stat_par


##4-3)statistical analysis (non parametric)====
stat_npr=data.table()
for (i in 1:length(comp)) {
  temp=subset(pah_chp,pah_chp$variable==comp[i])
  
  tm=as.data.table(aggregate(temp$value, by=list(`variable`=temp$Group), FUN=median))
  tm=tm[order(tm$x, decreasing = T),]
  
  temp$Group=factor(temp$Group,levels = tm$variable)
  ans <- posthoc.kruskal.dunn.test(rel ~ Group,data=temp ,p.adjust="bonferroni") ##
  ans.p <- get.pvalues(ans)
  ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  new=generate_label(temp, flev = "Group",value.var = "value",offset = 0, cld_info = ans.mcV)
  new$Comp=comp[i]
  stat_npr=rbind.data.frame(stat_npr, new)
  
}
stat_npr=stat_npr %>% `colnames<-`(c("Group","labels","rel","Comp"))
stat_npr$Comp2=factor(stat_npr$Comp,levels = c("CH","CHN","CHS","CHO","CHNO","CHOS"),
                      labels = c("CH","CHN","CHS","CHO","CHON","CHOS"))

stat_npr

norm

##4-4)plotting=====
stat=stat_par

ggplot(data=pah_chp, aes(x=Group,y=value,fill=Group))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(alpha=1, outlier.color = NA)+
  #geom_text(data =stat, aes(x=Group, y=ifelse(variable=="Mean AI",value+0.005,
  #                                            ifelse(variable=="Mean DBE",value+0.1,
  #                                                   ifelse(variable=="Mean H/C",value+0.01,value+0.005))), label=labels), 
  #          col="black", size=6)+
  facet_wrap(.~variable, scales = "free_y", nrow=2, dir = "h")+
  scale_x_discrete(name="")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#3C5488FF"))+
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
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("pah_Chp_compare_dunn"),height = 30, width = 30, units = "cm", dpi = 300)
  

#5.Combustion prop======
pah_merge$`DBE/C`=pah_merge$DBE/pah_merge$`C#`
pah_merge$Comb=ifelse(pah_merge$`DBE/C`>0.7, "Combustion derived","organic")

pah_merge$Comb_inty=ifelse(pah_merge$Comb=="organic",0,pah_merge$Normarea)


table(pah_merge$Comp)

pah_inty_all

comb_inty=as.data.table(aggregate(pah_merge$Normarea, by=list(`Sample`=pah_merge$Sample,`Group`=pah_merge$Group,Comb=pah_merge$Comb), FUN=sum))
comb_inty

comb_inty=comb_inty %>% inner_join(pah_inty_all[,c("Sample","tot")],by = "Sample")

comb_inty$rel=round(comb_inty$x/comb_inty$tot*100,1)

comb_sel=subset(comb_inty,comb_inty$Comb=="Combustion derived")

comb_sel$Group=factor(comb_sel$Group,levels = c("Seoul","Seosan","Beijing","Ulaanbaatar"),
                      labels = c("Seoul","Seosan","Beijing","Ulaanbaatar"))


##5-1)normality, equal variace test====
evtmp=levene.test(comb_sel$rel,comb_sel$Group)
stmp=shapiro.test(comb_sel$rel)

##5-2)statistical analysis (non parametric)====
stat_npr=data.table()
comb=unique(comb_sel$Comb)
for (i in 1:length(comb)) {
  temp=subset(comb_sel,comb_sel$Comb==comb[i])
  
  tm=as.data.table(aggregate(temp$rel, by=list(`variable`=temp$Group), FUN=median))
  tm=tm[order(tm$x, decreasing = T),]
  
  temp$Group=factor(temp$Group,levels = tm$variable)
  ans <- posthoc.kruskal.dunn.test(rel ~ Group,data=temp ,p.adjust="fdr") ##
  ans.p <- get.pvalues(ans)
  ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  new=generate_label(temp, flev = "Group",value.var = "rel",offset = 1.5, cld_info = ans.mcV)
  new$Comp=comb[i]
  stat_npr=rbind.data.frame(stat_npr, new)
  
}

stat_npr=stat_npr %>% `colnames<-`(c("Group","labels","rel","Comb"))

stat_npr
stat=stat_npr
norm

ggplot()+
  stat_boxplot(data=comb_sel, aes(x=Group, y=rel, fill=Group),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=comb_sel, aes(x=Group, y=rel, fill=Group),alpha=1, outlier.color = NA)+
  #geom_text(data =comp_max, aes(x=1, y=x, label=Comp), col="white", size=0)+
  geom_text(data =stat, aes(x=Group, y=rel, label=labels), 
            col="black", size=6)+
  #facet_grid(Event~Comp2, scales = "free")+
  #facet_wrap(.~Comp2, scales = "free_y", nrow=2, dir = "h")+
  scale_x_discrete(name="")+
  scale_y_continuous(name = "Proportion of combustion derived compounds (%)",
                     expand = c(0.02,0.02))+
  #scale_fill_manual(values = npg_no_jp)+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#3C5488FF"))+
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
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("pah_combustion_compare"),height = 20, width = 20, units = "cm", dpi = 300)

