library(ggplot2)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(PMCMRplus)
options(java.parameters = "- Xmx8192m")
library(extrafont)
loadfonts(device="win")
library(ggsci)
library(tidyr)
library("scales")
library(multcompView)
library(DescTools)
library(lawstat)
library(agricolae)
library(PMCMR)
library(plyr)


source("Rscript/func_filename.R")
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")
insert_minor_1 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_2 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}


#0.data loading====
FT_envi=fread("Datafile/FRIENDs_envi_sel.csv")
#1. Environmetal data analysis=====
FT_envi_sel=FT_envi %>% filter(Group%in%c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))

FT_envi_sel$WSOCbb=2.94*FT_envi_sel$Levoglucosan
FT_envi_sel$WSOCnbb=FT_envi_sel$WSOC-FT_envi_sel$WSOCbb
#FT_envi_sel$`SOC (%)`=FT_envi_sel$WSOCnbb
#FT_envi_sel$`POC (%)`=FT_envi_sel$OC-FT_envi_sel$`SOC`
FT_envi_sel

FT_envi_sel$`SOC (%)`=FT_envi_sel$WSOCnbb/FT_envi_sel$OC
FT_envi_sel$`POC (%)`=1-FT_envi_sel$`SOC (%)`

FT_envi_sel

SOA_m=melt(FT_envi_sel[,c("Sample","Group","No","Event","POC (%)","SOC (%)","WSOC","OC")],
               id.vars = c("Sample","Group","No","Event","WSOC","OC"), na.rm = T)
SOA_m

SOA_m$Event=factor(SOA_m$Event, levels = c("Event","Normal","Non-event"))

SOA_m=SOA_m %>% `colnames<-`(c("Sample","Group","No","Event","WSOC","OC","variable","val"))

ggplot(data=SOA_m, aes(x=variable, y=val, fill=Group))+
  geom_boxplot(outlier.color = NA)+
  facet_grid(.~Group,scales = "free")

  #geom_violin()+
  #geom_point(aes(shape=Group), position = position_dodge(width = 0.75))+
  #geom_text_repel(aes(label=No),position = position_dodge(width = 0.75), size=2)





stat_npr=data.table()
var=unique(SOA_m$variable)
grp=unique(SOA_m$Group)

for (a in 1:length(pd)) {
  
  #a=1
  data_in=subset(SOA_m,SOA_m$Pd==pd[a])
  
  for (i in 1:length(var)) {
    #i=1
    temp=subset(data_in,data_in$variable==var[i])
    
    tm=as.data.table(aggregate(temp$val, by=list(`Group`=temp$Group), FUN=median))
    tm=tm[order(tm$x, decreasing = T),]
    
    temp$Group=factor(temp$Group,levels = tm$Group)
    ans <- PMCMRplus::kwAllPairsDunnTest(val ~ Group,data=temp ,p.adjust="fdr") ##
    ans.p <- get.pvalues(ans)
    ans.mcV <- multcompLetters(ans.p, threshold=0.05)
    
    
    #new=generate_label(temp, flev = "variable",value.var = "val",offset = 0.0, cld_info = ans.mcV)
    
    flev="Group"
    value.var="val"
    cld_info=ans.mcV
    data_raw=as.data.frame(temp)
    
    plot.labels <- names(cld_info[['Letters']])
    # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
    # upper quantile and label placement
    iqr.df <- ddply(data_raw, flev, function (x) ((1.5*IQR(x[,value.var])+quantile(x[,value.var],0.75)))) %>% 
      `colnames<-`(c(flev,"iqr"))
    data_raw=data_raw %>% inner_join(iqr.df, by=flev)
    d2=subset(data_raw,data_raw[,value.var]<data_raw$iqr)
    
    boxplot.df <- ddply(d2, flev, function (x) (max(fivenum(x[,value.var])))+offset)
    #boxplot.df <- ddply(d2, flev, function (x) (2*max(fivenum(x[,value.var]))-min(fivenum(x[,value.var])))*offset)
    # Create a data frame out of the factor levels and Tukey's homogenous group letters
    plot.levels <- data.frame(plot.labels, labels = ans.mcV[['Letters']],
                              stringsAsFactors = FALSE)
    # Merge it with the labels
    labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = F)
    
    new=labels.df
    new$Group=var[i]
    new$Pd=pd[a]
    stat_npr=rbind.data.frame(stat_npr, new)
    
  }
}

stat_npr_grp=stat_npr
#stat_npr_var=stat_npr
stat_npr_grp
stat_npr_var


SOA_m$Pd=factor(SOA_m$Pd, levels = c("19s","19w","20w"))
ggplot()+
  stat_boxplot(data=SOA_m, aes(x=variable, y=val, fill=variable),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=SOA_m, aes(x=variable, y=val, fill=variable),alpha=1, outlier.color = NA)+
  #geom_point(data=SOA_m_ev, aes(x=Event, y=val, fill=variable), position = position_dodge(width = 0.75))+
  #geom_text_repel(data=SOA_m_ev, aes(x=Event, y=val, fill=variable, label=No), position = position_dodge(width = 0.75))+
  scale_y_continuous(name="",breaks = seq(0,1,0.25),labels = scales::percent_format(accuracy = 1), limits = c(0,1.2))+
  facet_wrap(.~Group, strip.position = "top", scales = "free", ncol=2, dir = "h")+
  #scale_fill_manual(values = c("grey70", "#FFFFFF"))+
  scale_fill_manual(values = c("#9DBCD4","#CB7723"))+
  #scale_fill_manual(values = c("#3C5488FF","#E64B35FF", "#4DBBD5FF", "#00A087FF","#F39B7FFF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.0),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 16,colour = "black", angle = 0, hjust = 0.5,vjust = 1.0),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15,"cm"),
        axis.text.y = element_text(size = 16, colour = "black" ),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        axis.title.y.left = element_text(size = 18, colour = "black",margin = unit(c(0.0,0.0,0.0,0.0),"cm")),
        axis.title.y.right = element_text(size = 18, colour = "black",angle = 90,margin = unit(c(0.2,0.2,0.2,0.5),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.2,0.4,0.2,0.2),"cm")),
        strip.text.y = element_blank(),
        strip.placement = "outside",
        #legend.margin = unit(c(0.0,0.0,0.0,0.0),"cm"),
        legend.title = element_text(size = 22, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 18, colour = "black",margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        legend.key.width = unit(1.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.background = element_blank(),
        legend.direction = "vertical",
        legend.justification=c(0.5, 0.5),
        legend.position = c(0.62,0.05))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("B&G_SOC_compare"),height = 30, width = 30, units = "cm", dpi = 300)



SOA_m$rel=SOA_m$val*100

grp=as.vector(unique(SOA_m$Group))
var=as.vector(unique(SOA_m$variable))
dt=data.table()
in_data=SOA_m

  for (i in 1:length(pd)) {
    #  i=1
    temp=subset(in_data,in_data$Pd==pd[i])
    for (j in 1:length(grp)) {
     #      j=1
      temp2=subset(temp,temp$Group==grp[j])
      t_poa=subset(temp2,temp2$variable=="POC (%)")
      t_soa=subset(temp2,temp2$variable=="SOC (%)")
      s1=shapiro.test(t_poa$rel)
      s2=shapiro.test(t_soa$rel)
      evtmp=levene.test(temp2$rel,temp2$variable)
      #tt=t.test(t_eve$value,t_nev$value, paired = F,var.equal = F)
      #tt=t.test(t_poa$rel,t_soa$rel, paired = T,var.equal = T)
      tt=t.test(t_poa$rel,t_soa$rel, paired = F,var.equal = T)
      #tt=wilcox.test(t_poa$rel,t_soa$rel, exact = F)
      new=data.table("Pd"=pd[i],"Group"=grp[j],"norm_poa"=round(s1$p.value,3),"norm_soa"=round(s2$p.value,3),
                     "Eq-var"=round(evtmp$p,3),"T-test"=round(tt$p.value,3))
      
      dt=rbind(dt,new)
      
    }
  }
  
dt  


##4-2-2)Calculate label position=====

SOA_m$rel=SOA_m$val*100
SOA_m_d=dcast(SOA_m,Group~variable, c(mean, sd), value.var = "rel")
SOA_m_d
#fwrite(SOA_m_d,file = "Datafile/SOA_table.csv")

SOA_m
SOA_m_ev=SOA_m %>% filter(Event%in%c("Event","Non-event"))
SOA_m_ev$Event2=ifelse(SOA_m_ev$Event=="Non-event","Nonevent",as.character(SOA_m_ev$Event))
SOA_m_ev
SOA_m_ev$val=ifelse(SOA_m_ev$val<0,(SOA_m_ev$val*-1),SOA_m_ev$val)


SOA_m_ev

SOA_m

stat_npr=data.table()
var=unique(SOA_m_ev$variable)
grp=unique(SOA_m_ev$Group)
ev=unique(SOA_m_ev$Event2)

data_in=SOA_m
for (i in 1:length(var)) {
    #i=1
    temp=subset(data_in,data_in$variable==var[i])
    
    tm=as.data.table(aggregate(temp$rel, by=list(`Group`=temp$Group), FUN=median))
    tm=tm[order(tm$x, decreasing = T),]
    
    temp$Group=factor(temp$Group,levels = tm$Group)
    ans <- PMCMRplus::kwAllPairsDunnTest(rel ~ Group,data=temp ,p.adjust="fdr") ##
    ans.p <- get.pvalues(ans)
    ans.mcV <- multcompLetters(ans.p, threshold=0.05)
    
    
    #new=generate_label(temp, flev = "variable",value.var = "val",offset = 0.0, cld_info = ans.mcV)
    
    flev="Group"
    value.var="rel"
    cld_info=ans.mcV
    data_raw=as.data.frame(temp)
    
    plot.labels <- names(cld_info[['Letters']])
    # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
    # upper quantile and label placement
    iqr.df_up <- ddply(data_raw, flev, function (x) ((1.5*IQR(x[,value.var])+quantile(x[,value.var],0.75)))) %>% 
      `colnames<-`(c(flev,"iqr_up"))
    
    iqr.df_down <- ddply(data_raw, flev, function (x) ((-1.5*IQR(x[,value.var])+quantile(x[,value.var],0.25)))) %>% 
      `colnames<-`(c(flev,"iqr_down"))
    
    
    data_raw=data_raw %>% inner_join(iqr.df_up, by=flev)
    d2=subset(data_raw,data_raw[,value.var]<data_raw$iqr_up)
    
    d2_raw=d2 %>% inner_join(iqr.df_down, by=flev)
    d2=subset(d2_raw,d2_raw[,value.var]>d2_raw$iqr_down)
    d2
    
    ans2 <- PMCMRplus::kwAllPairsDunnTest(rel ~ Group,data=d2 ,p.adjust="fdr") ##
    ans2.p <- get.pvalues(ans2)
    ans2.mcV <- multcompLetters(ans2.p, threshold=0.05)
    
    
    boxplot.df <- ddply(d2, flev, function (x) (max(fivenum(x[,value.var]))))
    #boxplot.df <- ddply(d2, flev, function (x) (2*max(fivenum(x[,value.var]))-min(fivenum(x[,value.var])))*offset)
    # Create a data frame out of the factor levels and Tukey's homogenous group letters
    plot.levels <- data.frame(plot.labels, labels = ans2.mcV[['Letters']],
                              stringsAsFactors = FALSE)
    
    # Merge it with the labels
    labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = F)
    
    new=labels.df
    new$variable=var[i]
    new$Group=grp[b]
    
    stat_npr=rbind.data.frame(stat_npr, new)
    
  }
  

stat_npr
stat_npr=stat_npr %>% `colnames<-`(c("Event","labels","val","variable","Group"))


for (b in 1:length(grp)) {
    #b=1
    data_in2=subset(data_in,data_in$Group==grp[b])
    
    for (i in 1:length(var)) {
      #i=2
      temp=subset(data_in2,data_in2$variable==var[i])
      
      tm=as.data.table(aggregate(temp$val, by=list(`Event2`=temp$Event2), FUN=median))
      tm=tm[order(tm$x, decreasing = T),]
      
      temp$Event2=factor(temp$Event2,levels = tm$Event2)
      ans <- PMCMRplus::kwAllPairsDunnTest(val ~ Event2,data=temp ,p.adjust="none") ##
      ans.p <- get.pvalues(ans)
      ans.mcV <- multcompLetters(ans.p, threshold=0.05)
      
      
      #new=generate_label(temp, flev = "variable",value.var = "val",offset = 0.0, cld_info = ans.mcV)
      
      flev="Event2"
      value.var="val"
      cld_info=ans.mcV
      data_raw=as.data.frame(temp)
      
      plot.labels <- names(cld_info[['Letters']])
      # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
      # upper quantile and label placement
      iqr.df_up <- ddply(data_raw, flev, function (x) ((1.5*IQR(x[,value.var])+quantile(x[,value.var],0.75)))) %>% 
        `colnames<-`(c(flev,"iqr_up"))
      
      iqr.df_down <- ddply(data_raw, flev, function (x) ((-1.5*IQR(x[,value.var])+quantile(x[,value.var],0.25)))) %>% 
        `colnames<-`(c(flev,"iqr_down"))
      
      
      data_raw=data_raw %>% inner_join(iqr.df_up, by=flev)
      d2=subset(data_raw,data_raw[,value.var]<data_raw$iqr_up)
      
      d2_raw=d2 %>% inner_join(iqr.df_down, by=flev)
      d2=subset(d2_raw,d2_raw[,value.var]>d2_raw$iqr_down)
      d2
      
      ans2 <- PMCMRplus::kwAllPairsDunnTest(val ~ Event2,data=d2 ,p.adjust="none") ##
      ans2.p <- get.pvalues(ans2)
      ans2.mcV <- multcompLetters(ans2.p, threshold=0.05)
      
      
      boxplot.df <- ddply(d2, flev, function (x) (max(fivenum(x[,value.var]))))
      #boxplot.df <- ddply(d2, flev, function (x) (2*max(fivenum(x[,value.var]))-min(fivenum(x[,value.var])))*offset)
      # Create a data frame out of the factor levels and Tukey's homogenous group letters
      plot.levels <- data.frame(plot.labels, labels = ans2.mcV[['Letters']],
                                stringsAsFactors = FALSE)
      # Merge it with the labels
      labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = F)
      
      new=labels.df
      new$variable=var[i]
      new$Group=grp[b]
      
      stat_npr=rbind.data.frame(stat_npr, new)
      
    }
    
}

stat_npr
stat_npr=stat_npr %>% `colnames<-`(c("Event","labels","val","variable","Group"))

SOA_m_ev
stat_npr$Event=ifelse(stat_npr$Event=="Nonevent","Non-event",as.character(stat_npr$Event))


stat_npr$Event=factor(stat_npr$Event,levels = c("Event","Normal","Non-event"))
stat_npr$labels2=ifelse(stat_npr$variable=="SOC (%)",toupper(stat_npr$labels),stat_npr$labels)
stat_npr[3,3]=0.4875612


stat_npr=data.table()
var=unique(SOA_m_ev$variable)
grp=unique(SOA_m_ev$Group)
ev=unique(SOA_m_ev$Event2)
data_in=SOA_m_ev

grp=as.vector(unique(SOA_m_ev$Group))
var=as.vector(unique(SOA_m_ev$variable))
dt=data.table()
in_data=subset(SOA_m_ev,SOA_m_ev$Sample!="Seoul_8")

dt=data.table()
for (i in 1:length(var)) {
  #  i=1
  temp=subset(in_data,in_data$var==var[i])
  for (j in 1:length(grp)) {
   #       j=1
    temp2=subset(temp,temp$Group==grp[j])
    t_ev=subset(temp2,temp2$Event=="Event")
    t_nev=subset(temp2,temp2$Event=="Non-event")
    s1=shapiro.test(t_ev$rel)
    s2=shapiro.test(t_nev$rel)
    evtmp=levene.test(temp2$rel,temp2$Event)
    #tt=t.test(t_eve$value,t_nev$value, paired = F,var.equal = F)
    #tt=t.test(t_poa$rel,t_soa$rel, paired = T,var.equal = T)
    tt=t.test(t_ev$rel,t_nev$rel, paired = F,var.equal = T)
    #tt=wilcox.test(t_poa$rel,t_soa$rel, exact = F)
    new=data.table("Variable"=var[i],"Group"=grp[j],"norm_eve"=round(s1$p.value,3),"norm_nev"=round(s2$p.value,3),
                   "Eq-var"=round(evtmp$p,3),"T-test"=round(tt$p.value,3))
    
    dt=rbind(dt,new)
    
  }
}
dt  

uf=data.table()
for (i in 1:length(var)) {
  #i=1
  temp=subset(in_data,in_data$var==var[i])
  for (j in 1:length(grp)) {
    #        j=1
    temp2=subset(temp,temp$Group==grp[j])
    temp2
    
    iqr.df <- ddply(temp2, "Event", function (x) ((1.5*IQR(x[,"val"])+quantile(x[,"val"],0.75)))) %>% 
      `colnames<-`(c("Event","iqr"))
    
    t3=temp2 %>% inner_join(iqr.df)
    
    t4=subset(t3,t3$val<t3$iqr)
    t4
    
    uf.ps <- ddply(t4, "Event", function (x) (max(fivenum(x[,"val"]))))
    new=data.table("variable"=var[i],"Group"=grp[j],"uf_ev"=uf.ps[1,2],"uf_nev"=uf.ps[2,2] )
    uf=rbind(uf,new)
    
  }
}
uf

uf$uf_pos=ifelse(uf$uf_ev>=uf$uf_nev,uf$uf_ev,uf$uf_nev)

uf=uf %>% inner_join(dt)
uf

uf$sig=ifelse(uf$`T-test`<0.001,"***",
              ifelse(uf$`T-test`<0.01,"**",
                     ifelse(uf$`T-test`<0.05,"*","N.S")))

uf$Group=factor(uf$Group,levels = c("Ulaanbaatar","Beijing","Seoul"))
uf$xs=ifelse(uf$variable=="POC (%)",1,2)

SOA_m_ev$Group=factor(SOA_m_ev$Group,levels = c("Ulaanbaatar","Beijing","Seoul"))

ggplot()+
  stat_boxplot(data=SOA_m_ev, aes(x=variable, y=val, fill=Event),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=SOA_m_ev, aes(x=variable, y=val, fill=Event),alpha=1, outlier.color = NA)+
  geom_segment(data=uf,aes(x=xs-0.2,xend=xs+0.2,
                           y=uf_pos+0.05 , yend=uf_pos+0.05 ))+
  #geom_segment(data=uf,aes(x=xs-0.2,xend=xs-0.2,
  #                         y=uf_pos+0.05 , yend=uf_ev+0.05))+
  #geom_segment(data=uf,aes(x=xs+0.2,xend=xs+0.2,
  #                         y=uf_pos+0.10 , yend=uf_SOC+0.05))+
  #geom_text(data =comp_max, aes(x=1, y=x, label=Comp), col="white", size=0)+
  geom_text(data =subset(uf,uf$sig!="N.S"), aes(x=variable, y=uf_pos*1.10, label=sig), col="black", size=6)+
  geom_text(data =subset(uf,uf$sig=="N.S"), aes(x=variable, y=uf_pos+0.1, label=sig), col="black", size=6)+
  #geom_text(data =stat_npr, aes(x=variable, y=val+0.10, label=labels2, Group=Event),
  #          position = position_dodge(width = 0.75),
  #          col="black", size=6)+
  #geom_point(data=SOA_m_ev, aes(x=Event, y=val, fill=variable), position = position_dodge(width = 0.75))+
  #geom_text_repel(data=SOA_m_ev, aes(x=Event, y=val, fill=variable, label=No), position = position_dodge(width = 0.75))+
  scale_y_continuous(name="",breaks = seq(0,1,0.25),labels = scales::percent_format(accuracy = 1), limits = c(0,1.2))+
  facet_wrap(.~Group, strip.position = "top", scales = "free", ncol=3, dir = "h")+
  scale_fill_manual(values = c("grey60", "#FFFFFF"))+
  #scale_fill_manual(values = c("#9DBCD4","#CB7723"))+
  #scale_fill_manual(values = c("#3C5488FF","#E64B35FF", "#4DBBD5FF", "#00A087FF","#F39B7FFF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.0),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 16,colour = "black", angle = 0, hjust = 0.5,vjust = 1.0),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15,"cm"),
        axis.text.y = element_text(size = 16, colour = "black" ),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        axis.title.y.left = element_text(size = 18, colour = "black",margin = unit(c(0.0,0.0,0.0,0.0),"cm")),
        axis.title.y.right = element_text(size = 18, colour = "black",angle = 90,margin = unit(c(0.2,0.2,0.2,0.5),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.3,0.4,0.3,0.2),"cm")),
        strip.text.y = element_blank(),
        strip.placement = "outside",
        #legend.margin = unit(c(0.0,0.0,0.0,0.0),"cm"),
        legend.title = element_text(size = 22, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 16, colour = "black",margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        legend.key.width = unit(1.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.background = element_blank(),
        legend.direction = "vertical",
        legend.justification=c(0.5, 0.5),
        legend.position = c(0.75,0.11))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("FRIEND_1st_SOC_event_compare"),height = 15, width = 45, units = "cm", dpi = 300)



SOA_m
FT_envi_sel

SOA_envi=FT_envi_sel[,c("Sample","Group","No","Event","POC (%)","SOC (%)","WSOCbb","WSOCnbb",
                        "NH4+","NO3-","SO42-","WSOC","WISOC","OC")]
SOA_envi


SOA_m_ev

SOA_m_ev2=subset(SOA_m_ev,SOA_m_ev$Sample!="Seoul_8") %>% droplevels()
tt=dcast(SOA_m_ev2, Group+Event~variable, value.var="rel", fun.aggregate = c(mean,sd))
tt
