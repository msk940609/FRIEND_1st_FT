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
library(lawstat)
library(agricolae)
library(PMCMR)
library(PMCMRplus)

source("Rscript/func_filename.R")
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")
#0.Load data====
FT_merge=fread("Datafile/FRIENDs_1st_FT.csv")
FT_merge$Sample=paste(FT_merge$Group,FT_merge$No,sep = "_")

FT_merge$C=ifelse(FT_merge$C.>0, "C","")
FT_merge$H=ifelse(FT_merge$H.>0, "H","")
FT_merge$O=ifelse(FT_merge$O.>0, "O","")
FT_merge$N=ifelse(FT_merge$N.>0, "N","")
FT_merge$S=ifelse(FT_merge$S.>0, "S","")

FT_merge=FT_merge %>% unite("Comp",c("C","H","O","N","S"),sep = "")
FT_merge$Comp=ifelse(FT_merge$O.==0, "Remainders",FT_merge$Comp)
FT_merge=subset(FT_merge,FT_merge$Comp!="Remainders") %>% droplevels()

#1.Chemical Composition distribution=====
##1-1)build chemical composition data====
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

#fwrite(ft_inty,"Datafile/Chemicalcomposition_1st")
source("Rscript/func_upper_fence_label.R")

##1-2)normality, equal variace test====
comp=unique(ft_inty$Comp)
stat_par=data.table()
norm=data.table()
for (i in 1:length(comp)) {
  temp=subset(ft_inty,ft_inty$Comp==comp[i])
  evtmp=levene.test(temp$rel,temp$Group)
  grp=unique(temp$Group)
  for (j in 1:length(grp)) {
    tmp2=subset(temp,temp$Group==grp[j])
    stmp=shapiro.test(tmp2$rel)
    new=data.table("Comp"=comp[i],"Group"=grp[j],"p"=round(stmp$p.value,3), "EV"=round(evtmp$p.value,3))
    norm=rbind(norm,new)
  }
}

norm

##1-3)statistical analysis (parametric)====
comp=unique(ft_inty$Comp)
stat_par=data.table()
for (i in 1:length(comp)) {
  temp=subset(ft_inty,ft_inty$Comp==comp[i])
  
  tm=as.data.table(aggregate(temp$rel, by=list(`variable`=temp$Group), FUN=median))
  tm=tm[order(tm$x, decreasing = T),]
  
  temp$Group=factor(temp$Group,levels = tm$variable)
  #ans <- posthoc.kruskal.dunn.test(rel ~ Group,data=temp ,p.adjust="none") ##
  #ans.p <- get.pvalues(ans)
  #ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  ans=dunnettT3Test(rel ~ Group,data=temp, p.adjust.methods="fdr")
  ans.p <- get.pvalues(ans)
  ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  new=generate_label(temp, flev = "Group",value.var = "rel",offset = 0, cld_info = ans.mcV)
  new$Comp=comp[i]
  stat_par=rbind.data.frame(stat_par, new)
  
}

stat_par=stat_par %>% `colnames<-`(c("Group","labels","rel","Comp"))
stat_par$Comp=factor(stat_par$Comp,levels = c("CHO","CHON","CHOS","CHONS"))
stat_par$Comp=factor(stat_par$Comp,levels = c("CHO","CHON","CHOS","CHONS"))
stat_par


##1-4)statistical analysis (non parametric)====
comp=unique(ft_inty$Comp)
stat_npr=data.table()
for (i in 1:length(comp)) {
  temp=subset(ft_inty,ft_inty$Comp==comp[i])
  
  tm=as.data.table(aggregate(temp$rel, by=list(`variable`=temp$Group), FUN=median))
  tm=tm[order(tm$x, decreasing = T),]
  
  temp$Group=factor(temp$Group,levels = tm$variable)
  ans <- kwAllPairsDunnTest(rel ~ Group,data=temp ,p.adjust="bonferroni") ##
  ans.p <- get.pvalues(ans)
  ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  
  new=generate_label(temp, flev = "Group",value.var = "rel",offset = 0, cld_info = ans.mcV)
  new$Comp=comp[i]
  stat_npr=rbind.data.frame(stat_npr, new)
  
}
stat_npr=stat_npr %>% `colnames<-`(c("Group","labels","rel","Comp"))
stat_npr$Comp=factor(stat_npr$Comp,levels = c("CHO","CHON","CHOS","CHONS"))
stat_npr$Comp=factor(stat_npr$Comp,levels = c("CHO","CHON","CHOS","CHONS"))
stat_npr

norm
stat_par_sel=stat_par %>% filter(Comp%in%c("CHO","CHON"))
stat_npr_sel=stat_npr %>% filter(Comp%in%c("CHOS","CHONS"))
stat_npr
stat=rbind(stat_par_sel,stat_npr_sel)
  
comp_max=as.data.table(aggregate(ft_inty$rel, by=list(`Comp`=ft_inty$Comp), FUN=max))
comp_max$x=comp_max$x*1.2

pal=pal_npg("nrc")(5)

##1-5)plotting=====


ft_inty$Complab=paste(ft_inty$Comp,"(%)")
ft_inty$Complab=factor(ft_inty$Complab,levels = c("CHO (%)","CHON (%)","CHOS (%)","CHONS (%)"))

comp_max$Complab=paste(comp_max$Comp,"(%)")
comp_max$Complab=factor(comp_max$Complab,levels = c("CHO (%)","CHON (%)","CHOS (%)","CHONS (%)"))
stat$Complab=paste(stat$Comp,"(%)")
stat$Complab=factor(stat$Complab,levels = c("CHO (%)","CHON (%)","CHOS (%)","CHONS (%)"))

ggplot()+
  stat_boxplot(data=ft_inty, aes(x=Group, y=rel, fill=Group),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=ft_inty, aes(x=Group, y=rel, fill=Group),alpha=1, outlier.color = NA)+
  geom_text(data =comp_max, aes(x=1, y=x, label=Comp), col="white", size=0)+
  geom_text(data =stat, aes(x=Group, y=ifelse(Complab=="CHO (%)",rel+4,
                                              ifelse(Complab=="CHON (%)",rel+2,
                                                     ifelse(Complab=="CHOS (%)",rel+3,rel+4))), label=labels), 
            col="black", size=6)+
  facet_wrap(Complab~., strip.position = "left", scales = "free", ncol=2)+
  #facet_wrap(.~Comp, scales = "fix", nrow=2, dir = "h")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
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
        strip.text.y = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.2,0.4,0.2,0.2),"cm")),
        strip.text.x = element_blank(),
        strip.placement = "outside",
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
  ggsave(filename("Comp_compare_bon"),height = 30, width = 31, units = "cm", dpi = 300)


#2.chemical properties====
##2-1)build chemical composition data====
FT_merge$Sample=paste(FT_merge$Group,FT_merge$No,sep = "_")

chp_m=melt(FT_merge[,c("Sample","O.C","H.C","N.C","S.C","DBE","AI")], id.vars = c("Sample"),
           measure.vars = c("O.C","H.C","N.C","S.C","DBE","AI")) %>% 
  dcast(Sample~variable,mean) %>% `colnames<-`(c("Sample","Mean O/C","Mean H/C", "Mean N/C", "Mean S/C", "Mean DBE", "Mean AI")) %>% 
  melt(id.var=c("Sample"))

chp_m=chp_m %>% separate(Sample, c("Group","no"),sep = "_")
chp_m$Sample=paste(chp_m$Group,chp_m$no,sep = "_")

chp_m$Group=factor(chp_m$Group,levels = c("SUL","SS","B","NT","UL"),
                   labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
chp_m=chp_m %>% `colnames<-`(c("Group","no","variable","val","Sample"))
chp_m$Sample=paste(chp_m$Group,chp_m$no,sep = "_")
chp_m

##fwrite(chp_m,"Datafile/Chemicalprop_1st.csv")

chp_m
##2-2)normality, equal variace test====
chp=unique(chp_m$variable)

norm=data.table()
for (i in 1:length(chp)) {
  temp=subset(chp_m,chp_m$variable==chp[i])
  evtmp=levene.test(temp$val,temp$Group)
  grp=unique(temp$Group)
  for (j in 1:length(grp)) {
    tmp2=subset(temp,temp$Group==grp[j])
    stmp=shapiro.test(tmp2$val)
    new=data.table("Comp"=chp[i],"Group"=grp[j],"p"=round(stmp$p.value,3), "EV"=round(evtmp$p.value,3))
    norm=rbind(norm,new)
  }
}

norm 
###Normality: O/C, S/C, DBE, AI (dunnettT3)
###DBE: equal variance: FSD or scheffe
###Non-normality: H/C, N/C (dunn)


##2-3)statistical analysis (parametric)====
stat_par=data.table()
for (i in 1:length(chp)) {
  temp=subset(chp_m,chp_m$variable==chp[i])
  
  tm=as.data.table(aggregate(temp$val, by=list(`variable`=temp$Group), FUN=median))
  tm=tm[order(tm$x, decreasing = T),]
  
  temp$Group=factor(temp$Group,levels = tm$variable)
  #ans <- posthoc.kruskal.dunn.test(rel ~ Group,data=temp ,p.adjust="none") ##
  #ans.p <- get.pvalues(ans)
  #ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  ans=dunnettT3Test(val ~ Group,data=temp, p.adjust.methods="fdr")
  #ans=lsdTest(val ~ Group,data=temp, p.adjust.methods="fdr")
  ans.p <- get.pvalues(ans)
  ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  new=generate_label(temp, flev = "Group",value.var = "val",offset = 0.0, cld_info = ans.mcV)
  new$chp=chp[i]
  stat_par=rbind.data.frame(stat_par, new)
  
}

stat_par=stat_par %>% `colnames<-`(c("Group","labels","value","variable"))
stat_par$variable=factor(stat_par$Comp,levels = c("Mean O/C","Mean H/C", "Mean N/C", "Mean S/C", "Mean DBE", "Mean AI"))
stat_par$Comp=factor(stat_par$Comp,levels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
stat_par


##2-4)statistical analysis (non parametric)====
stat_npr=data.table()
for (i in 1:length(chp)) {
  #i=1
  temp=subset(chp_m,chp_m$variable==chp[i])
  
  tm=as.data.table(aggregate(temp$val, by=list(`variable`=temp$Group), FUN=median))
  tm=tm[order(tm$x, decreasing = T),]
  
  temp$Group=factor(temp$Group,levels = tm$variable)
  ans <- PMCMRplus::kwAllPairsDunnTest(val ~ Group,data=temp ,p.adjust="fdr") ##
  ans.p <- get.pvalues(ans)
  ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  
  new=generate_label(temp, flev = "Group",value.var = "val",offset = 0.0, cld_info = ans.mcV)
  new$Comp=chp[i]
  stat_npr=rbind.data.frame(stat_npr, new)
  
}

stat_npr=stat_npr %>% `colnames<-`(c("Group","labels","value","variable"))
stat_npr$variable=factor(stat_npr$variable,levels = c("Mean O/C","Mean H/C", "Mean N/C", "Mean S/C", "Mean DBE", "Mean AI"))
stat_npr

norm
stat_par_sel=stat_par %>% filter(variable%in%c("Mean O/C","Mean S/C","Mean DBE","Mean AI"))
stat_npr_sel=stat_npr %>% filter(variable%in%c("Mean H/C","Mean N/C"))
stat_npr
stat=rbind(stat_par_sel,stat_npr_sel)

stat$variable=factor(stat$variable,levels = c("Mean O/C","Mean H/C", "Mean N/C", "Mean S/C", "Mean DBE", "Mean AI"))

chp_max=as.data.table(aggregate(chp_m$val, by=list(`variable`=chp_m$variable), FUN=max))
chp_max$x=chp_max$x*1.1
chp_max$min=0

pal=pal_npg("nrc")(5)

##2-5)plotting=====
ggplot()+
  stat_boxplot(data=chp_m, aes(x=Group, y=val, fill=Group),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=chp_m, aes(x=Group, y=val, fill=Group),alpha=1, outlier.color = NA)+
  #geom_text(data =chp_max, aes(x=1, y=x, label=variable), col="white", size=0)+
  geom_text(data =stat, aes(x=Group, y=ifelse(variable=="Mean O/C",value+0.015,
                                              ifelse(variable=="Mean H/C",value+0.035,
                                                     ifelse(variable=="Mean N/C",value+0.006,
                                                            ifelse(variable=="Mean S/C",value+0.006,
                                                                   ifelse(variable=="Mean DBE",value+0.5,value+0.015))))), label=labels), 
            col="black", size=6)+
  facet_wrap(variable~., strip.position = "left", scales = "free", ncol=3)+
  #facet_wrap(.~variable, scales = "free", nrow=2, dir = "v")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
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
        strip.placement = "outside",
        strip.text.y = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.2,0.4,0.2,0.2),"cm")),
        strip.text.x = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.4,0.2,0.2,0.2),"cm")),
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
  ggsave(filename("Chp_compare_fdr"),height = 30, width = 50, units = "cm", dpi = 300)

#3. vk plot & DBE vs #C plot=====
FT_merge

FT_merge

ft_vk=melt(FT_merge[,c("Sample","Group","Comp","Formula","O.C","H.C","Bromo.Inty")], 
           id.vars = c("Sample","Group","Comp","Formula","O.C","H.C")) %>% 
  dcast(Group~Formula, mean) %>% 
  melt(id.var="Group", na.rm = T) %>% `colnames<-`(c("Group","Formula","Inty"))

fm=unique(FT_merge[,c("Formula","O.C","H.C","DBE","C.","Comp","Calc.m.z")])
fm

ft_vk=ft_vk %>% inner_join(fm, by="Formula")
ft_vk


ft_vk$Group=factor(ft_vk$Group,levels = c("SUL","SS","B","NT","UL"),
                   labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))

table(ft_vk$Group)

ft_vk$id=ifelse(ft_vk$Group=="Beijing","Beijing (8095)",
                ifelse(ft_vk$Group=="Seoul","Seoul (9490)",
                       ifelse(ft_vk$Group=="Seosan","Seosan (9791)",
                              ifelse(ft_vk$Group=="Noto","Noto (5792)","Ulaanbaatar (7318)"))))

ft_vk$id=factor(ft_vk$id,levels = c("Seoul (9490)","Seosan (9791)","Beijing (8095)","Noto (5792)",
                                    "Ulaanbaatar (7318)"))

ggplot(ft_vk, aes(x=`O.C`, y=`H.C`, col=Comp))+
  geom_point(size=2, alpha=0.4)+
  facet_wrap(.~id, ncol = 5, scales = "fix")+
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
        legend.position = c(0.14,0.14),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("vk"),height = 18, width = 75, units = "cm", dpi = 700)




#4.Event vs non event=====
##4-1) vk& DBE plot=====
envi_1st=fread("Datafile/envi_1st_sel.csv")

FT_eve=FT_merge
FT_eve$Group=factor(FT_eve$Group,levels = c("SUL","SS","B","NT","UL"),
                    labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
FT_eve$Sample=paste(FT_eve$Group,FT_eve$No,sep = "_")
FT_eve=FT_eve %>% inner_join(envi_1st)

ft_vk_eve=melt(FT_eve[,c("Group","Event","Comp","Formula","O.C","H.C","Bromo.Inty")], 
               id.vars = c("Group","Event","Comp","Formula","O.C","H.C")) %>% 
  dcast(Group+Event~Formula, mean) %>% 
  melt(id.var=c("Group","Event"), na.rm = T) %>% `colnames<-`(c("Group","Event","Formula","Inty"))

fm=unique(FT_merge[,c("Formula","O.C","H.C","DBE","C.","Comp","Calc.m.z")])
fm

ft_vk_eve=ft_vk_eve %>% inner_join(fm, by="Formula")
ft_vk_eve

ft_vk_eve=subset(ft_vk_eve,ft_vk_eve$Event!="Normal") %>% droplevels()

ft_vk_eve$Group=factor(ft_vk_eve$Group,levels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
ft_vk_eve$Comp=factor(ft_vk_eve$Comp,levels = c("CHO","CHON","CHOS","CHONS"))

table(ft_vk_eve$Group)
ft_vk_eve=ft_vk_eve[order(ft_vk_eve$Inty),]

ggplot(ft_vk_eve, aes(x=`O.C`, y=`H.C`, col=Comp))+
  geom_point(aes(size=Inty), alpha=0.4)+
  #facet_wrap(.~id, ncol = 5, scales = "fix")+
  facet_grid(Event~Group, scales = "fix")+
  scale_size_continuous(range=c(2,20))+
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
        legend.position = c(0.13,0.08),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("vk_eve"),height = 30, width = 75, units = "cm", dpi = 700)


ggplot(ft_vk_eve, aes(x=`C.`, y=`DBE`, col=Comp))+
  geom_point(aes(size=Inty), alpha=0.3)+
  facet_grid(Event~Group, scales = "fix")+
  scale_size_continuous(range=c(2,15))+
  #facet_wrap(.~Group, ncol = 5, scales = "fix")+
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
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("DBE_eve"),height = 30, width = 75, units = "cm", dpi = 700)

##4-2) chemical composition====
ft_inty_all=as.data.table(aggregate(FT_merge$Bromo.Inty, by=list(`Sample`=FT_merge$Sample), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Tot"))
ft_inty=as.data.table(aggregate(FT_merge$Bromo.Inty, by=list(`Sample`=FT_merge$Sample, Comp=FT_merge$Comp), FUN=sum))

ft_inty=ft_inty %>% inner_join(ft_inty_all,by = "Sample")
ft_inty$rel=round(ft_inty$x/ft_inty$Tot*100,1)

ft_inty=ft_inty %>% separate(Sample, c("Group","no"),sep = "_")

ft_inty$Comp=factor(ft_inty$Comp,levels = c("CHO","CHON","CHOS","CHONS"))

ft_inty$Group=factor(ft_inty$Group,levels = c("SUL","SS","B","NT","UL"),
                     labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
ft_inty

ft_inty_eve=ft_inty

ft_inty_eve$Sample=paste(ft_inty_eve$Group,ft_inty_eve$no, sep = "_")

ft_inty_eve=ft_inty_eve %>% inner_join(envi_1st)
ft_inty_eve
table(ft_inty_eve$Event)

ft_inty_eve=subset(ft_inty_eve,ft_inty_eve$Event!="Normal") %>% droplevels()
ft_inty_eve_m=melt(ft_inty_eve[,c("Sample","Comp","Event","rel")],id.vars = c("Sample","Comp","Event"))
ft_inty_eve_m=ft_inty_eve_m %>% separate("Sample",c("Group","no"))

ft_inty_eve_m$Group=factor(ft_inty_eve_m$Group,levels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
ft_inty_eve_m

##4-2-1)Compare distribution (Event vsNon-event)=====
grp=as.vector(unique(ft_inty_eve_m$Group))
comp=as.vector(unique(ft_inty_eve_m$Comp))

dt=data.table()
for (i in 1:length(comp)) {
  temp=subset(ft_inty_eve_m,ft_inty_eve_m$Comp==comp[i])
  for (j in 1:length(grp)) {
    temp2=subset(temp,temp$Group==grp[j])
    
    t_eve=subset(temp2,temp2$Event=="Event")
    t_nev=subset(temp2,temp2$Event=="Non-event")
    if(grp[j]=="Ulaanbaatar"){
      
      new=data.table("Comp"=comp[i],"Group"=grp[j],"norm_eve"=1,"norm_nev"=1,
                     "Eq-var"=1,"T-test"=1)
    }else{
      s1=shapiro.test(t_eve$value)
      s2=shapiro.test(t_nev$value)
      evtmp=levene.test(temp2$value,temp2$Event)
      tt=t.test(t_eve$value,t_nev$value, paired = F,var.equal = F)
      new=data.table("Comp"=comp[i],"Group"=grp[j],"norm_eve"=round(s1$p.value,3),"norm_nev"=round(s2$p.value,3),
                     "Eq-var"=round(evtmp$p,3),"T-test"=round(tt$p.value,3))
      
    }
    dt=rbind(dt,new)
        
  }
}
dt

#dt_non_equal=dt
dt
##4-2-2)Calculate label position=====
uf=data.table()
for (i in 1:length(comp)) {
  #i=1
  temp=subset(ft_inty_eve_m,ft_inty_eve_m$Comp==comp[i])
  for (j in 1:length(grp)) {
    #    j=5
    temp2=subset(temp,temp$Group==grp[j])
    temp2
    
    t_eve=subset(temp2,temp2$Event=="Event")
    t_nev=subset(temp2,temp2$Event=="Non-event")
    
    iqr.df <- ddply(temp2, "Event", function (x) ((1.5*IQR(x[,"value"])+quantile(x[,"value"],0.75)))) %>% 
      `colnames<-`(c("Event","iqr"))
    
    t3=temp2 %>% inner_join(iqr.df)
    
    t4=subset(t3,t3$value<t3$iqr)
    t4
    
    uf.ps <- ddply(t4, "Event", function (x) (max(fivenum(x[,"value"]))))
    new=data.table("Comp"=comp[i],"Group"=grp[j],"uf_eve"=uf.ps[1,2],"uf_nev"=uf.ps[2,2] )
    uf=rbind(uf,new)
    
  }
}

uf$uf_pos=ifelse(uf$uf_eve>=uf$uf_nev,uf$uf_eve,uf$uf_nev)

uf$id=paste(uf$Comp,uf$Group,sep = "_")
dt$id=paste(dt$Comp,dt$Group,sep = "_")
dt$Noneq_t=dt_non_equal$`T-test`

dt$p=ifelse(dt$`T-test`<dt$Noneq_t,dt$`T-test`,dt$Noneq_t)

uf=uf %>% inner_join(dt[,c("id","p")])

uf
uf$sig=ifelse(uf$`p`<0.001,"***",
              ifelse(uf$`p`<0.01,"**",
                     ifelse(uf$`p`<0.05,"*","N.S")))
uf$Group=factor(uf$Group,levels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
uf$Comp=factor(uf$Comp,levels = c("CHO","CHON","CHOS","CHONS"))

uf$xs=ifelse(uf$Group=="Seoul",1,
             ifelse(uf$Group=="Seosan",2,
                    ifelse(uf$Group=="Beijing",3,
                           ifelse(uf$Group=="Noto",4,5))))
##4-2-3)plotting=====
ggplot()+
  stat_boxplot(data=ft_inty_eve_m, aes(x=Group, y=value, fill=Event),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=ft_inty_eve_m, aes(x=Group, y=value, fill=Event),alpha=1, outlier.color = NA)+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), lty=2)+
  geom_segment(data=subset(uf,uf$sig!="N.S"),aes(x=xs-0.2,xend=xs+0.2,y=uf_pos*1.03, yend=uf_pos*1.03))+
  #geom_text(data =comp_max, aes(x=1, y=x, label=Comp), col="white", size=0)+
  geom_text(data =subset(uf,uf$sig!="N.S"), aes(x=Group, y=uf_pos*1.05, label=sig), col="black", size=6)+
  geom_text(data =subset(uf,uf$sig!="N.S"), aes(x=Group, y=uf_pos*1.1, label=sig), col="white", size=6)+
  facet_wrap(.~Comp, scales = "free", nrow=2, dir = "h")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  #scale_fill_manual(values = c( "#4DBBD5FF", "#00A087FF", "#E64B35FF","#F39B7FFF","#3C5488FF"))+
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
        legend.position = c(0.01,0.4575))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("Comp_eve_"),height = 30, width = 30, units = "cm", dpi = 300)

##4-3) chemical properties====
FT_merge
table(FT_merge$Sample)
chp_m=melt(FT_merge[,c("Sample","O.C","H.C","N.C","S.C","DBE","AI")], id.vars = c("Sample"),
           measure.vars = c("O.C","H.C","N.C","S.C","DBE","AI")) %>% 
  dcast(Sample~variable,mean) %>% `colnames<-`(c("Sample","Mean O/C","Mean H/C", "Mean N/C", "Mean S/C", "Mean DBE", "Mean AI")) %>% 
  melt(id.var=c("Sample"))

chp_m

chp_m=chp_m %>% separate(Sample, c("Group","no"),sep = "_")
chp_m$Sample=paste(chp_m$Group,chp_m$no,sep = "_")


chp_m$Group=factor(chp_m$Group,levels = c("SUL","SS","B","NT","UL"),
                   labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
chp_m$Sample=paste(chp_m$Group,chp_m$no,sep = "_")

chp_m2=chp_m %>% inner_join(envi_1st, by=c("Sample","Group"))
#chp_m2=subset(chp_m2,chp_m2$Event!="Normal") %>% droplevels
table(chp_m2$Event)
chp_m2=chp_m2 %>% filter(Event%in%c("Event","Non-event"))

chp_m2$variable=factor(chp_m2$variable, levels=c("Mean O/C","Mean H/C","Mean N/C","Mean S/C","Mean DBE","Mean AI"))

chp_m2$Group=factor(chp_m2$Group, levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))


##4-3-1)Compare distribution (Event vsNon-event)=====
grp=as.vector(unique(chp_m2$Group))
var=as.vector(unique(chp_m2$variable))

dt=data.table()
for (i in 1:length(var)) {
#  i=1
  temp=subset(chp_m2,chp_m2$variable==var[i])
  for (j in 1:length(grp)) {
#    j=1
    temp2=subset(temp,temp$Group==grp[j])
    t_eve=subset(temp2,temp2$Event=="Event")
    t_nev=subset(temp2,temp2$Event=="Non-event")
    if(grp[j]=="Ulaanbaatar"){
      
      s1=shapiro.test(t_eve$value)
      s2=shapiro.test(t_nev$value)
      evtmp=levene.test(temp2$value,temp2$Event)
      #tt=t.test(t_eve$value,t_nev$value, paired = F,var.equal = F)
      tt=t.test(t_eve$value,t_nev$value, paired = F,var.equal = T)
      new=data.table("var"=var[i],"Group"=grp[j],"norm_eve"=round(s1$p.value,3),"norm_nev"=round(s2$p.value,3),
                     "Eq-var"=round(evtmp$p,3),"T-test"=round(tt$p.value,3))
    }else{
      s1=shapiro.test(t_eve$value)
      s2=shapiro.test(t_nev$value)
      evtmp=levene.test(temp2$value,temp2$Event)
      #tt=t.test(t_eve$value,t_nev$value, paired = F,var.equal = F)
      tt=t.test(t_eve$value,t_nev$value, paired = F,var.equal = T)
      new=data.table("var"=var[i],"Group"=grp[j],"norm_eve"=round(s1$p.value,3),"norm_nev"=round(s2$p.value,3),
                     "Eq-var"=round(evtmp$p,3),"T-test"=round(tt$p.value,3))
      
    }
    dt=rbind(dt,new)
    
  }
}
dt

#dt_non_equal=dt
dt$nont=dt_non_equal$`T-test`
##4-2-2)Calculate label position=====
uf=data.table()
for (i in 1:length(var)) {
#  i=1
  temp=subset(chp_m2,chp_m2$variable==var[i])
  for (j in 1:length(grp)) {
#        j=1
    temp2=subset(temp,temp$Group==grp[j])
    temp2
    
    t_eve=subset(temp2,temp2$Event=="Event")
    t_nev=subset(temp2,temp2$Event=="Non-event")
    
    iqr.df <- ddply(temp2, "Event", function (x) ((1.5*IQR(x[,"value"])+quantile(x[,"value"],0.75)))) %>% 
      `colnames<-`(c("Event","iqr"))
    
    t3=temp2 %>% inner_join(iqr.df)
    
    t4=subset(t3,t3$value<t3$iqr)
    t4
    
    uf.ps <- ddply(t4, "Event", function (x) (max(fivenum(x[,"value"]))))
    new=data.table("variable"=var[i],"Group"=grp[j],"uf_eve"=uf.ps[1,2],"uf_nev"=uf.ps[2,2] )
    uf=rbind(uf,new)
    
  }
}

uf$uf_pos=ifelse(uf$uf_eve>=uf$uf_nev,uf$uf_eve,uf$uf_nev)

uf$id=paste(uf$var,uf$Group,sep = "_")
dt$id=paste(dt$var,dt$Group,sep = "_")

#dt$p=ifelse(dt$`T-test`<dt$nont,dt$`T-test`,dt$nont)
dt$p=dt$`T-test`

uf=uf %>% inner_join(dt[,c("id","p")])

uf
uf$sig=ifelse(uf$`p`<0.001,"***",
              ifelse(uf$`p`<0.01,"**",
                     ifelse(uf$`p`<0.05,"*","N.S")))
uf$Group=factor(uf$Group,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))
uf$variable=factor(uf$variable,levels=c("Mean O/C","Mean H/C","Mean N/C","Mean S/C","Mean DBE","Mean AI"))

uf$xs=ifelse(uf$Group=="Seoul",3,
             ifelse(uf$Group=="Seosan",4,
                    ifelse(uf$Group=="Beijing",2,
                           ifelse(uf$Group=="Noto",5,1))))




ggplot()+
  stat_boxplot(data=chp_m2, aes(x=Group, y=value, fill=Event),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=chp_m2, aes(x=Group, y=value, fill=Event),alpha=1, outlier.color = NA)+
  facet_grid(Event~Comp, scales = "free")+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), lty=2)+
  geom_segment(data=subset(uf,uf$sig!="N.S"),aes(x=xs-0.2,xend=xs+0.2,
                                                 y=ifelse(variable=="Mean O/C",uf_pos+0.005,
                                                            ifelse(variable=="Mean H/C",uf_pos+0.008,
                                                                   ifelse(variable=="Mean N/C",uf_pos*1.03,
                                                                          ifelse(variable=="Mean S/C",uf_pos+0.006,
                                                                                 ifelse(variable=="Mean DBE",uf_pos+0.1,uf_pos+0.004))))), 
                                                 yend=ifelse(variable=="Mean O/C",uf_pos+0.005,
                                                             ifelse(variable=="Mean H/C",uf_pos+0.008,
                                                                    ifelse(variable=="Mean N/C",uf_pos*1.03,
                                                                           ifelse(variable=="Mean S/C",uf_pos+0.006,
                                                                                  ifelse(variable=="Mean DBE",uf_pos+0.1,uf_pos+0.004)))))))+
  #geom_text(data =comp_max, aes(x=1, y=x, label=Comp), col="white", size=0)+
  geom_text(data =subset(uf,uf$sig!="N.S"), aes(x=Group, y=ifelse(variable=="Mean O/C",uf_pos+0.008,
                                                                  ifelse(variable=="Mean H/C",uf_pos+0.015,
                                                                         ifelse(variable=="Mean N/C",uf_pos+0.003,
                                                                                ifelse(variable=="Mean S/C",uf_pos+0.010,
                                                                                       ifelse(variable=="Mean DBE",uf_pos+0.2,uf_pos+0.006))))), label=sig), col="black", size=6)+
  geom_text(data =subset(uf,uf$sig!="N.S"), aes(x=Group, y=uf_pos*1.15, label=sig), col="white", size=6)+
  #facet_wrap(.~variable, scales = "free", nrow=2, dir = "v")+
  facet_wrap(variable~., strip.position = "left", scales = "free", ncol=3, dir = "v")+
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
        axis.title.y.left = element_text(size = 18, colour = "black",margin = unit(c(0.2,0.0,0.2,0.0),"cm")),
        axis.title.y.right = element_text(size = 18, colour = "black",angle = 90,margin = unit(c(0.2,0.2,0.2,0.5),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.4,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.2,0.4,0.2,0.2),"cm")),
        strip.placement = "outside",
        #legend.margin = unit(c(0.0,0.0,0.0,0.0),"cm"),
        legend.title = element_text(size = 22, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 18, colour = "black",margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        legend.key.width = unit(1.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.background = element_rect(fill = "white"),
        legend.direction = "vertical",
        legend.justification=c(0, 1),
        legend.position = c(0.005,0.4820))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("Chp_compare_event"),height = 30, width = 50, units = "cm", dpi = 300)

#5.FT*Envi====
envi_all=fread("Datafile/FRIENDs_envi.csv")
envi_all

##intensity based correlation====
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


ft_inty
ft_comp=dcast(ft_inty,Group+no~Comp, value.var = c("rel","x"))
ft_comp$Sample=paste(ft_comp$Group,ft_comp$no,sep = "_")
ft_comp


ft_comp=ft_comp %>% inner_join(envi_all)
ft_comp
table(ft_comp$Group)

ft_comp_sel=ft_comp %>% filter(Group%in%c("Seoul","Beijing"))
ft_comp_sel

#fwrite(ft_comp_sel,file = "Datafile/ft_comp_sel.csv")

ft_comp_sel

ft_comp_sel_m=melt(ft_comp_sel[,c("Group","No","Sample","rel_CHON","rel_CHOS","rel_CHONS","x_CHON","x_CHOS","x_CHONS","SO42-","NO3-")],
                   id.vars = c("Group","No","Sample","SO42-","NO3-"), variable.name = "Comp")
ft_comp_sel_mm=melt(ft_comp_sel_m, id.vars = c("Group","No","Sample","Comp","value"), variable.name = "Envi",value.name = "Conc")
ft_comp_sel_mm=ft_comp_sel_mm %>% separate(Comp, c("type","Comp"),sep = "_")

ft_comp_sel_mm_s=subset(ft_comp_sel_mm,ft_comp_sel_mm$Group=="Seoul")
ft_comp_sel_mm_s_rel=subset(ft_comp_sel_mm_s,ft_comp_sel_mm_s$type=="rel")
ft_comp_sel_mm_s_val=subset(ft_comp_sel_mm_s,ft_comp_sel_mm_s$type=="x")

ft_comp_sel_mm_s_rel$Comp=factor(ft_comp_sel_mm_s_rel$Comp, levels = c("CHON","CHONS","CHOS"))
ft_comp_sel_mm_s_rel$Envi=factor(ft_comp_sel_mm_s_rel$Envi, levels = c("NO3-","SO42-"))

ft_comp_sel_mm_s
ft_comp_sel_mm_s_rel_max=as.data.table(aggregate(ft_comp_sel_mm_s$Conc, by=list(`Group`=ft_comp_sel_mm_s$Group,Envi=ft_comp_sel_mm_s$Envi), FUN=max))
ft_comp_sel_mm_s_rel_max

ggplot(ft_comp_sel_mm_s_rel, aes(x=value,y=Conc, col=Group))+
  #geom_smooth(method = "loess",se=F)+
  geom_blank(data=ft_comp_sel_mm_s_rel_max, aes(x=30,y=x*1.2,Group=Group))+
  geom_smooth(method = "lm",se=F,col="black", formula = y~log10(x))+
  geom_point(col="royalblue", size=3)+
  facet_grid(Envi~Comp, scales = "free")+
  scale_y_continuous(name = expression(bold("Conc."~"("*"\u03bcg/"*m^"3"*")")))+
  scale_x_continuous(name = "Proprotion of chemical composition (%)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.4,0,0.2),"cm")),
        strip.text = element_text(size = 20, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsEnvi_inty_rel_S"),height = 20, width = 30, units = "cm", dpi = 300)


ft_comp_sel_mm_s_val$Comp=factor(ft_comp_sel_mm_s_val$Comp, levels = c("CHON","CHONS","CHOS"))
ft_comp_sel_mm_s_val$Envi=factor(ft_comp_sel_mm_s_val$Envi, levels = c("NO3-","SO42-"))

ft_comp_sel_mm_s
ft_comp_sel_mm_s_val_max=as.data.table(aggregate(ft_comp_sel_mm_s$Conc, by=list(`Group`=ft_comp_sel_mm_s$Group,Envi=ft_comp_sel_mm_s$Envi), FUN=max))
ft_comp_sel_mm_s_val_max

ggplot(ft_comp_sel_mm_s_val, aes(x=value,y=Conc, col=Group))+
  geom_blank(data=ft_comp_sel_mm_s_val_max, aes(x=30,y=x*1.2,Group=Group))+
  geom_smooth(method = "lm",se=F,col="black", formula = y~log10(x))+
  geom_point(col="royalblue", size=3)+
  facet_grid(Envi~Comp, scales = "free")+
  scale_y_continuous(name = expression(bold("Conc."~"("*"\u03bcg/"*m^"3"*")")))+
  scale_x_continuous(name = "Total intensity of chemical composition")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.4,0,0.2),"cm")),
        strip.text = element_text(size = 20, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsEnvi_inty_val_S"),height = 20, width = 30, units = "cm", dpi = 300)
  

ft_comp_sel_mm_b=subset(ft_comp_sel_mm,ft_comp_sel_mm$Group=="Beijing")
ft_comp_sel_mm_b_rel=subset(ft_comp_sel_mm_b,ft_comp_sel_mm_b$type=="rel")
ft_comp_sel_mm_b_val=subset(ft_comp_sel_mm_b,ft_comp_sel_mm_b$type=="x")

ft_comp_sel_mm_b_rel$Comp=factor(ft_comp_sel_mm_b_rel$Comp, levels = c("CHON","CHONS","CHOS"))
ft_comp_sel_mm_b_rel$Envi=factor(ft_comp_sel_mm_b_rel$Envi, levels = c("NO3-","SO42-"))

ft_comp_sel_mm_b
ft_comp_sel_mm_b_rel_max=as.data.table(aggregate(ft_comp_sel_mm_b$Conc, by=list(`Group`=ft_comp_sel_mm_b$Group,Envi=ft_comp_sel_mm_b$Envi), FUN=max))
ft_comp_sel_mm_b_rel_max

ggplot(ft_comp_sel_mm_b_rel, aes(x=value,y=Conc, col=Group))+
  #geom_smooth(method = "loess",se=F)+
  geom_blank(data=ft_comp_sel_mm_b_rel_max, aes(x=30,y=x*1.2,Group=Group))+
  geom_smooth(method = "lm",se=F,col="black", formula = y~log10(x))+
  geom_point(col="orangered3", size=3)+
  facet_grid(Envi~Comp, scales = "free")+
  scale_y_continuous(name = expression(bold("Conc."~"("*"\u03bcg/"*m^"3"*")")))+
  scale_x_continuous(name = "Proprotion of chemical composition (%)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.4,0,0.2),"cm")),
        strip.text = element_text(size = 20, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsEnvi_inty_rel_B"),height = 20, width = 30, units = "cm", dpi = 300)


ft_comp_sel_mm_b_val$Comp=factor(ft_comp_sel_mm_b_val$Comp, levels = c("CHON","CHONS","CHOS"))
ft_comp_sel_mm_b_val$Envi=factor(ft_comp_sel_mm_b_val$Envi, levels = c("NO3-","SO42-"))

ft_comp_sel_mm_b
ft_comp_sel_mm_b_val_max=as.data.table(aggregate(ft_comp_sel_mm_b$Conc, by=list(`Group`=ft_comp_sel_mm_b$Group,Envi=ft_comp_sel_mm_b$Envi), FUN=max))
ft_comp_sel_mm_b_val_max

ggplot(ft_comp_sel_mm_b_val, aes(x=value,y=Conc, col=Group))+
  geom_blank(data=ft_comp_sel_mm_s_val_max, aes(x=30,y=x*1.2,Group=Group))+
  geom_smooth(method = "lm",se=F,col="black", formula = y~log10(x))+
  geom_point(col="orangered3", size=3)+
  facet_grid(Envi~Comp, scales = "free")+
  scale_y_continuous(name = expression(bold("Conc."~"("*"\u03bcg/"*m^"3"*")")))+
  scale_x_continuous(name = "Total intensity of chemical composition")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.4,0,0.2),"cm")),
        strip.text = element_text(size = 20, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsEnvi_inty_val_B"),height = 20, width = 30, units = "cm", dpi = 300)

##FT_envi correlation======
ft_comp_sel_mm$type=ifelse(ft_comp_sel_mm$type=="rel","rel","val")


grp=unique(ft_comp_sel_mm$Group)
grp
typ=unique(ft_comp_sel_mm$type)
comp=unique(ft_comp_sel_mm$Comp)
envi=unique(ft_comp_sel_mm$Envi)

dt=data.table()
for (i in 1:length(grp)) {
#  i=1
  temp=subset(ft_comp_sel_mm,ft_comp_sel_mm$Group==grp[i])
  
  for (j in 1:length(type)) {
 #   j=1
    temp2=subset(temp,temp$type==typ[j])
    temp2
    
    for (k in 1:length(comp)) {
  #    k=1
      tmp=subset(temp2,temp2$Comp==comp[k])
      
      for (l in 1:length(Envi)) {
   #     l=1
        tmp2=subset(tmp,tmp$Envi==envi[l])
        
        cor=cor.test(tmp2$value,tmp2$Conc, method = "spearman",exact = F)
    #    cor
        
        new_row=data.table("Group"=grp[i],"Type"=type[j],"Comp"=comp[k],"Envi"=envi[l],"Cor"=cor$estimate, "P"=round(cor$p.value,3))
        dt=rbind(dt,new_row)
      }
    }
  }
}


dt


##freqeuncy based correlation======
FT_merge$Freq=1

ft_freq_all=as.data.table(aggregate(FT_merge$Freq, by=list(`Sample`=FT_merge$Sample), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Tot"))

ft_freq=as.data.table(aggregate(FT_merge$Freq, by=list(`Sample`=FT_merge$Sample, Comp=FT_merge$Comp), FUN=sum))
ft_freq

ft_freq=ft_freq %>% inner_join(ft_freq_all,by = "Sample")
ft_freq$rel=round(ft_freq$x/ft_freq$Tot*100,1)


ft_freq=ft_freq %>% separate(Sample, c("Group","no"),sep = "_")
ft_freq

ft_freq$Comp=factor(ft_freq$Comp,levels = c("CHO","CHON","CHOS","CHONS"))

ft_freq$Group=factor(ft_freq$Group,levels = c("SUL","SS","B","NT","UL"),
                     labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
ft_freq

ft_comp_freq=dcast(ft_freq,Group+no~Comp, value.var = c("rel","x"))
ft_comp_freq$Sample=paste(ft_comp_freq$Group,ft_comp_freq$no,sep = "_")
ft_comp_freq


ft_comp_freq=ft_comp_freq %>% inner_join(envi_all)
ft_comp_freq
table(ft_comp_freq$Group)

ft_comp_freq_sel=ft_comp_freq %>% filter(Group%in%c("Seoul","Beijing"))
ft_comp_freq_sel
table(ft_comp_freq_sel$Group)
#fwrite(ft_comp_sel,file = "Datafile/ft_comp_sel.csv")

ft_comp_freq_sel

ft_comp_freq_sel_m=melt(ft_comp_freq_sel[,c("Group","No","Sample","rel_CHON","rel_CHOS","rel_CHONS","x_CHON","x_CHOS","x_CHONS","SO42-","NO3-","SO2","NO")],
                   id.vars = c("Group","No","Sample","SO42-","NO3-","SO2","NO"), variable.name = "Comp")
ft_comp_freq_sel_mm=melt(ft_comp_freq_sel_m, id.vars = c("Group","No","Sample","Comp","value"), variable.name = "Envi",value.name = "Conc")

ft_comp_freq_sel_mm=ft_comp_freq_sel_mm %>% separate(Comp, c("type","Comp"),sep = "_")

ft_comp_freq_sel_mm_s=subset(ft_comp_freq_sel_mm,ft_comp_freq_sel_mm$Group=="Seoul")
ft_comp_freq_sel_mm_s_rel=subset(ft_comp_freq_sel_mm_s,ft_comp_freq_sel_mm_s$type=="rel")
ft_comp_freq_sel_mm_s_s_val=subset(ft_comp_freq_sel_mm_s,ft_comp_freq_sel_mm_s$type=="x")

ggplot(ft_comp_freq_sel_mm_s_rel, aes(x=value,y=Conc, col=Group))+
  geom_smooth(method = "loess",se=F)+
  geom_point()+
  facet_grid(Envi~Comp, scales = "free")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggsave(filename("compvsEnvi_freq_rel"),height = 20, width = 30, units = "cm", dpi = 300)

ggplot(ft_comp_freq_sel_mm_s_s_val, aes(x=value,y=Conc, col=Group))+
  geom_smooth(method = "loess",se=F)+
  geom_point()+
  facet_grid(Envi~Comp, scales = "free")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggsave(filename("compvsEnvi_freq_val"),height = 20, width = 30, units = "cm", dpi = 300)


###pie chart=====
envi_1st=fread("Datafile/envi_1st_sel.csv")
ft_inty

#ft_inty$Sample=paste(ft_inty$Group,ft_inty$no,sep = "_")
ft_inty$Sample=paste(ft_inty$Group,ft_inty$no,sep = "_")
ft_inty

ft_inty_en=ft_inty %>% inner_join(envi_1st[,c("Sample","Event","PM2.5")])
ft_inty_en

ft_inty_en_d=dcast(ft_inty_en,Group~Comp, value.var = "rel",mean) %>% 
  melt(id.vars=c("Group"))
ft_inty_en_d_sel=ft_inty_en_d %>% filter(Group%in%c("Seoul","Beijing"))

ggplot()+
  geom_bar(data=ft_inty_en_d_sel,aes(x="",y=value, fill=variable), col="black",stat="identity", size=1.5,
           position = position_fill(reverse = T))+
  coord_polar("y")+
  scale_fill_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF"))+
  facet_wrap(.~Group, ncol = 1)+
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  theme_void()+
  theme(legend.text = element_text(size=10, hjust=0.5))+
  ggsave("Region_pie.png",height = 20, width = 10, units = "cm", dpi = 300)



ft_inty_ev=subset(ft_inty_en,ft_inty_en$Event=="Event") %>% 
  dcast(Sample+Group+no+PM2.5~Comp,value.var = "rel")

ft_inty_ev #Ulaanbaatar_16,Beijing_13, Seoul_9, Seosan_8, Noto_9

ft_pie_ev=ft_inty_ev %>% filter(Sample%in%c("Ulaanbaatar_16","Beijing_13", "Seoul_9", "Seosan_8", "Noto_9"))

ft_pie_ev_m=melt(ft_pie_ev,id.vars = c("Sample","Group","No","PM2.5"))


ggplot()+
  geom_bar(data=ft_pie_ev_m,aes(x="",y=value, fill=variable), col="black",stat="identity", size=2.5,
           position = position_fill(reverse = T))+
  coord_polar("y")+
  scale_fill_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF"))+
  facet_wrap(.~Group)+
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  theme_void()+
  theme(legend.text = element_text(size=10, hjust=0.5))+
  ggsave("Region_piechar_ev.png",height = 20, width = 30, units = "cm", dpi = 300)




ft_inty_nev=subset(ft_inty_en,ft_inty_en$Event=="Non-event") %>% 
  dcast(Sample+Group+No+PM2.5~Comp,value.var = "rel")

ft_inty_nev #Ulaanbaatar_15,Beijing_25, Seoul_1, Seosan_23, Noto_25

ft_pie_nev=ft_inty_nev %>% filter(Sample%in%c("Ulaanbaatar_15","Beijing_25", "Seoul_1", "Seosan_23", "Noto_25"))

ft_pie_nev_m=melt(ft_pie_nev,id.vars = c("Sample","Group","No","PM2.5"))

ggplot()+
  geom_bar(data=ft_pie_nev_m,aes(x="",y=value, fill=variable), col="black",stat="identity", size=2.5,
           position = position_fill(reverse = T))+
  coord_polar("y")+
  scale_fill_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF"))+
  facet_wrap(.~Group)+
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  theme_void()+
  theme(legend.text = element_text(size=10, hjust=0.5))+
  ggsave("Region_piechar_nev.png",height = 20, width = 30, units = "cm", dpi = 300)

###pie, eventvsnon-event WSOA cor==


sul_ev_cor=fread("Datafile/pie_sul_9_ev_sel.csv")
sul_nev_cor=fread("Datafile/pie_sul_1_nev_sel.csv")
b_ev_cor=fread("Datafile/pie_b_13_ev_sel.csv")
b_nev_cor=fread("Datafile/pie_b_25_nev_sel.csv")


pi_cor=rbind(sul_ev_cor,sul_nev_cor,b_ev_cor,b_nev_cor)
pi_cor

pi_cor$C=ifelse(pi_cor$C.>0, "C","")
pi_cor$H=ifelse(pi_cor$H.>0, "H","")
pi_cor$O=ifelse(pi_cor$O.>0, "O","")
pi_cor$N=ifelse(pi_cor$N.>0, "N","")
pi_cor$S=ifelse(pi_cor$S.>0, "S","")

pi_cor=pi_cor %>% unite("Comp",c("C","H","O","N","S"),sep = "")
pi_cor$Comp=ifelse(pi_cor$O.==0, "Remainders",pi_cor$Comp)
pi_cor$Sample=paste(pi_cor$Group,pi_cor$No,sep = "_")

pi_cor=subset(pi_cor,pi_cor$Comp!="Remainders") %>% droplevels()

pi_cor_inty_all=as.data.table(aggregate(pi_cor$Bromo.Inty, by=list(`Sample`=pi_cor$Sample), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Tot"))

pi_cor_inty_all

pi_cor_inty=as.data.table(aggregate(pi_cor$Bromo.Inty, by=list(`Sample`=pi_cor$Sample, Comp=pi_cor$Comp), FUN=sum))
pi_cor_inty

pi_cor_inty=pi_cor_inty %>% inner_join(pi_cor_inty_all,by = "Sample")
pi_cor_inty$rel=round(pi_cor_inty$x/pi_cor_inty$Tot*100,1)

pi_cor_inty$Event=ifelse(pi_cor_inty$Sample=="B_13","Event",
                         ifelse(pi_cor_inty$Sample=="B_25","Non-event",
                                ifelse(pi_cor_inty$Sample=="SUL_1","Non-event","Event")))


pi_cor_inty=pi_cor_inty %>% separate(Sample, c("Group","no"),sep = "_")

pi_cor_inty

pi_cor_inty$Comp=factor(pi_cor_inty$Comp,levels = c("CHO","CHON","CHOS","CHONS"))

ggplot()+
  geom_bar(data=pi_cor_inty,aes(x="",y=rel, fill=Comp), col="black",stat="identity", size=2.5,
           position = position_fill(reverse = T))+
  coord_polar("y")+
  scale_fill_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF"))+
  facet_wrap(Event~Group)+
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  theme_void()+
  theme(legend.text = element_text(size=10, hjust=0.5))+
  ggsave("Region_piechar_SOA.png",height = 20, width = 30, units = "cm", dpi = 300)

