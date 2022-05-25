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
library(multcompView)
library(lawstat)
library(agricolae)
library(PMCMR)
library(PMCMRplus)

getOption("digits")
options("digits" = 15)

source("Rscript/func_filename.R")
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")

ft_merge=fread("Datafile/FRIENDs_1st_FT.csv")
ft_merge=ft_merge %>% mutate(Comp = gsub('[[:digit:]]', '', Formula))
ft_merge$Comp=ifelse(ft_merge$O.==0,"Remainders",ft_merge$Comp)
ft_merge$Group=factor(ft_merge$Group, levels =c("UL","B","SS","SUL","NT"),
                      labels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))
ft_merge$Sample=paste(ft_merge$Group,ft_merge$No, sep = "_")

fm_obs

ft_merge=ft_merge %>% inner_join(fm_obs)

ft_merge_sel=subset(ft_merge,ft_merge$cnt>4)
##1-4)statistical analysis (non parametric)====

ft_inty_all=as.data.table(aggregate(ft_merge_sel$Bromo.Inty, by=list(`Sample`=ft_merge_sel$Sample), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Tot"))

ft_inty=as.data.table(aggregate(ft_merge_sel$Bromo.Inty, by=list(`Sample`=ft_merge_sel$Sample, Comp=ft_merge_sel$Comp), FUN=sum))
ft_inty

ft_inty=ft_inty %>% inner_join(ft_inty_all,by = "Sample")
ft_inty$rel=round(ft_inty$x/ft_inty$Tot*100,1)

ft_inty=ft_inty %>% separate(Sample, c("Group","no"),sep = "_")
ft_inty
ft_inty_sel=subset(ft_inty,ft_inty$Comp!="Remainders") %>% droplevels()

ft_inty_sel$Comp=factor(ft_inty_sel$Comp,levels = c("CHO","CHNO","CHOS","CHNOS"),
                        labels = c("CHO","CHON","CHOS","CHONS"))

ft_inty_sel$Group=factor(ft_inty_sel$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))

ft_inty_sel

comp=unique(ft_inty_sel$Comp)
data_in=ft_inty_sel
stat_npr=data.table()
for (i in 1:length(comp)) {
  
  temp=subset(data_in,data_in$Comp==comp[i])
  
  tm=as.data.table(aggregate(temp$rel, by=list(`variable`=temp$Group), FUN=median))
  tm=tm[order(tm$x, decreasing = T),]
  
  temp$Group=factor(temp$Group,levels = tm$variable)
  ans <- kwAllPairsDunnTest(rel ~ Group,data=temp ,p.adjust="fdr") ##
  ans.p <- get.pvalues(ans)
  ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  
  new=generate_label(temp, flev = "Group",value.var = "rel",offset = 4, cld_info = ans.mcV)
  new$Comp=comp[i]
  stat_npr=rbind.data.frame(stat_npr, new)
  
}

stat_npr=stat_npr %>% `colnames<-`(c("Group","labels","rel","Comp"))
stat_npr$Comp=factor(stat_npr$Comp,levels = c("CHO","CHON","CHOS","CHONS"))
stat_npr

stat=stat_npr

comp_max=as.data.table(aggregate(ft_inty_sel$rel, by=list(`Comp`=ft_inty_sel$Comp), FUN=max))
comp_max$x=comp_max$x*1.2

pal=pal_npg("nrc")(5)

ft_inty_sel$compvar=paste(ft_inty_sel$Comp, " (%)",sep = "")
ft_inty_sel$compvar=factor(ft_inty_sel$compvar,levels = c("CHO (%)","CHON (%)","CHOS (%)","CHONS (%)"))

stat$compvar=paste(stat$Comp, " (%)",sep = "")
stat$compvar=factor(stat$compvar,levels = c("CHO (%)","CHON (%)","CHOS (%)","CHONS (%)"))

##1-5)plotting=====

ggplot()+
  stat_boxplot(data=ft_inty_sel, aes(x=Group, y=rel, fill=Group),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=ft_inty_sel, aes(x=Group, y=rel, fill=Group),alpha=1, outlier.color = NA)+
  facet_rep_wrap(.~compvar, nrow=2, repeat.tick.labels = T, strip.position="left", scales="fix")
  

ggplot()+
  stat_boxplot(data=ft_inty_sel, aes(x=Group, y=rel, fill=Group),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=ft_inty_sel, aes(x=Group, y=rel, fill=Group),alpha=1, outlier.color = NA)+
  geom_text(data =ft_inty_sel, aes(x=1, y=rel*1.1, label=Comp), col=NA, size=0)+
  #geom_text(data =stat, aes(x=Group, y=rel, label=labels), 
  #          col="black", size=8)+
  geom_text(data =stat, aes(x=Group, y=ifelse(Comp=="CHO",rel-0.5,
                                                  ifelse(Comp=="CHON",rel-2.5,
                                                         ifelse(Comp=="CHOS",rel-2,
                                                                ifelse(Comp=="CHONS",rel+0,
                                                                       ifelse(Comp=="Mean DBE",rel+3,rel+0.012))))), label=labels), 
            col="black", size=8)+
  facet_rep_wrap(.~compvar, nrow=2, repeat.tick.labels = T, strip.position="left", scales="free")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  scale_fill_manual(values = c( "#3C5488FF","#E64B35FF", "#00A087FF","#4DBBD5FF","#F39B7FFF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.0),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 17,colour = "black", angle = 0, hjust = 0.5,vjust = 1.0),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15,"cm"),
        axis.text.y = element_text(size = 17, colour = "black", margin = unit(c(0.0,0.0,0.2,0.2),"cm")),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        #axis.title.y.left = element_text(size = 0, colour = "black",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        #axis.title.y.right = element_text(size = 0, colour = "black",angle = 90,margin = unit(c(0.2,0.2,0.2,0.5),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.placement = "outside",
        strip.text.y = element_text(size = 16, colour = "black",face = 2,margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        strip.text.x = element_blank(),
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
  ggsave(filename("Comp_compare_fdr"),height = 27, width = 36, units = "cm", dpi = 300)


ggplot()+
  #stat_boxplot(data=ft_inty_sel, aes(x=Group, y=rel, fill=Group),geom='errorbar', linetype=1, width=0.25,
  #             position = position_dodge(width = 0.75))+
  #geom_boxplot(data=ft_inty_sel, aes(x=Group, y=rel, fill=Group),alpha=1, outlier.color = NA)+
  geom_violin(data=ft_inty_sel, aes(x=Group, y=rel, fill=Group),alpha=1)+
  geom_text(data =ft_inty_sel, aes(x=1, y=rel*1.1, label=Comp), col=NA, size=0)+
  #geom_text(data =stat, aes(x=Group, y=rel, label=labels), 
  #          col="black", size=8)+
  geom_text(data =stat, aes(x=Group, y=ifelse(Comp=="CHO",rel-0.5,
                                              ifelse(Comp=="CHON",rel-2.5,
                                                     ifelse(Comp=="CHOS",rel-2,
                                                            ifelse(Comp=="CHONS",rel+0,
                                                                   ifelse(Comp=="Mean DBE",rel+3,rel+0.012))))), label=labels), 
            col="black", size=8)+
  facet_rep_wrap(.~compvar, nrow=2, repeat.tick.labels = T, strip.position="left", scales="free")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  scale_fill_manual(values = c( "#3C5488FF","#E64B35FF", "#00A087FF","#4DBBD5FF","#F39B7FFF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.0),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 17,colour = "black", angle = 0, hjust = 0.5,vjust = 1.0),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15,"cm"),
        axis.text.y = element_text(size = 17, colour = "black", margin = unit(c(0.0,0.0,0.2,0.2),"cm")),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        #axis.title.y.left = element_text(size = 0, colour = "black",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        #axis.title.y.right = element_text(size = 0, colour = "black",angle = 90,margin = unit(c(0.2,0.2,0.2,0.5),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.placement = "outside",
        strip.text.y = element_text(size = 16, colour = "black",face = 2,margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        strip.text.x = element_blank(),
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
  ggsave(filename("Comp_compare_violin"),height = 27, width = 36, units = "cm", dpi = 300)

##chp=====
ft_chp=melt(ft_merge_sel[,c("Group","No","AI","DBE","O.C","H.C")],id.vars = c("Group","No")) %>% 
  dcast(Group+No~variable,mean) %>% 
  melt(id.vars=c("Group","No"))

ft_chp=ft_chp %>% `colnames<-`(c("Group","No","variable","val"))

chp=unique(ft_chp$variable)
data_in=ft_chp
stat_npr=data.table()
for (i in 1:length(chp)) {
  temp=subset(data_in,data_in$variable==chp[i])
  
  tm=as.data.table(aggregate(temp$val, by=list(`variable`=temp$Group), FUN=median))
  tm=tm[order(tm$x, decreasing = T),]
  
  temp$Group=factor(temp$Group,levels = tm$variable)
  ans <- kwAllPairsDunnTest(val ~ Group,data=temp ,p.adjust="fdr") ##
  ans.p <- get.pvalues(ans)
  ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  new=generate_label(temp, flev = "Group",value.var = "val",offset = 0, cld_info = ans.mcV)
  new$chp=chp[i]
  stat_npr=rbind.data.frame(stat_npr, new)
  
}

stat_chp=stat_npr %>% `colnames<-`(c("Group","labels","val","variable"))
stat_chp$variable=factor(stat_chp$variable,levels = c("AI","DBE","O.C","H.C"),
                       labels = c("Mean AI","Mean DBE","Mean O/C", "Mean H/C"))
stat_chp

ft_chp
ft_chp$variable=factor(ft_chp$variable,levels = c("AI","DBE","O.C","H.C"),
                         labels = c("Mean AI","Mean DBE","Mean O/C", "Mean H/C"))
ft_chp$Group=factor(ft_chp$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))


ft_chp
stat_chp


ggplot()+
  stat_boxplot(data=ft_chp, aes(x=Group, y=val, fill=Group),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=ft_chp, aes(x=Group, y=val, fill=Group),alpha=1, outlier.color = NA)+
  #geom_text(data =chp_max, aes(x=1, y=x, label=variable), col="white", size=0)+
  geom_text(data =stat_chp, aes(x=Group, y=ifelse(variable=="Mean O/C",val+0.010,
                                              ifelse(variable=="Mean H/C",val+0.030,
                                                     ifelse(variable=="Mean N/C",val+0.006,
                                                            ifelse(variable=="Mean S/C",val+0.006,
                                                                   ifelse(variable=="Mean DBE",val+0.30,val+0.012))))), label=labels), 
            col="black", size=8)+
  
  geom_text(data =stat_chp, aes(x=Group, y=val*1.05), label=NA, col="white", size=0)+
  facet_rep_wrap(.~variable, nrow=2, repeat.tick.labels = T, strip.position="left", scales="free")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  scale_fill_manual(values = c( "#3C5488FF","#E64B35FF", "#00A087FF","#4DBBD5FF","#F39B7FFF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.0),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 17,colour = "black", angle = 0, hjust = 0.5,vjust = 1.0),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15,"cm"),
        axis.text.y = element_text(size = 17, colour = "black", margin = unit(c(0.0,0.0,0.2,0.2),"cm")),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        #axis.title.y.left = element_text(size = 0, colour = "black",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        #axis.title.y.right = element_text(size = 0, colour = "black",angle = 90,margin = unit(c(0.2,0.2,0.2,0.5),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.placement = "outside",
        strip.text.y = element_text(size = 16, colour = "black",face = 2,margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        strip.text.x = element_blank(),
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
  ggsave(filename("Chp_compare_fdr"),height = 27, width = 36, units = "cm", dpi = 300)

