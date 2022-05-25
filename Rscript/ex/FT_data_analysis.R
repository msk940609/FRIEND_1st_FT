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


###NMDS, PCA=======
PC_1st=prcomp(MDS_1st_in[,-c(1)], center = T, scale. = F)
PC_1st 

gpca_1st=as.data.table(PC_1st$x)
gpca_1st$Sample=MDS_1st_in$Sample

gpca_1st=gpca_1st %>% tidyr::separate("Sample",c("Group","no"),sep="_")
gpca_1st

ggplot()+
  geom_point(data = gpca_1st, aes(x=PC1, y=PC2, col=Group))


MDS_1st=melt(FT_merge[,c("Sample","Formula","Bromo.Inty")], id.vars = c("Sample","Formula")) %>% 
  dcast(Sample~Formula, sum)

MDS_1st[,1:2]

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


##nmds formula=====
MDS_1st[,1:2]

tmds=dcast(melt(MDS_1st, id.vars = "Sample"), variable ~ Sample)
tmds

fwrite(tmds, "Datafile/specific.csv")

spec=fread("Datafile/specific.csv")

NMDS_1st ##Stress 0.08
NMDS_1st$species
nmds_formula=as.data.table(NMDS_1st$species) %>% mutate("Formula"=row.names(NMDS_1st$species))
nmds_formula=nmds_formula %>% inner_join(spec)
nmds_formula


FT_formula=FT_merge[,c("Formula","C.","DBE","H.C","O.C","Comp")]
FT_formula=unique(FT_formula)

nmds_formula=nmds_formula %>% inner_join(FT_formula)
nmds_formula

nmds_formula

nmds_formula$Region=factor(nmds_formula$Region,
                           levels = c("Seoul","Seosan","Beijing","Noto","NUM","Common"),
                           labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar","Common"))


nmds_formula$order=ifelse(nmds_formula$Region=="Seoul",6,
                          ifelse(nmds_formula$Region=="Seosan",5,
                                ifelse(nmds_formula$Region=="Beijing",4,
                                      ifelse(nmds_formula$Region=="Noto",3,
                                            ifelse(nmds_formula$Region=="Ulaanbaatar",2,1)))))
nmds_formula=nmds_formula[order(nmds_formula$order),]
nmds_formula

pal=pal_jama("default")(5)


ggplot()+
  geom_vline(xintercept = 0, lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(data = nmds_formula, aes(x=MDS1, y=MDS2, col=Region), size=2)+
  scale_color_manual(values = c(pal,"grey70"))+
  #scale_color_jama()+
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
        legend.text = element_text(size = 16, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.2,0,0.2,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.76, 0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = NULL,title.hjust = 0.5, ncol = 2, override.aes = list(size=7)))+
  ggsave(filename("MDS_1st_formula"), height = 20, width = 20, units = "cm", dpi = 700)

table(nmds_formula$Region)


###Cvs DBE ====
nmds_formula_spec=subset(nmds_formula,nmds_formula$Specific=="Specific")

nmds_formula_spec$Comp=factor(nmds_formula_spec$Comp,levels = c("CHO","CHON","CHOS","CHONS"))

nmds_formula_spec$Region=factor(nmds_formula_spec$Region,
                           levels = c("Seoul","Seosan","Beijing","Noto","NUM"),
                           labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))


ggplot()+
  geom_point(data = nmds_formula_spec, aes(x=C., y=DBE, col=Comp), size=2,alpha=0.3)+
  facet_wrap(.~Region, ncol = 3, scales = "fix")+
  scale_x_continuous(name = "#C")+
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
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0,0.4,0,0.2),"cm")),
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
  guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)))+
  ggsave(filename("specific_CvsDBE_comp"),height = 20, width = 30, units = "cm", dpi = 300)


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

#fwrite(ft_inty,"Datafile/Chemicalcomposition_1st")

source("Rscript/func_upper_fence_label.R")

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


comp=unique(ft_inty$Comp)
stat=data.table()
for (i in 1:length(comp)) {
  temp=subset(ft_inty,ft_inty$Comp==comp[i])
  
  tm=as.data.table(aggregate(temp$rel, by=list(`variable`=temp$Group), FUN=median))
  tm=tm[order(tm$x, decreasing = T),]
  
  temp$Group=factor(temp$Group,levels = tm$variable)
  ans <- posthoc.kruskal.dunn.test(rel ~ Group,data=temp ,p.adjust="fdr")
  ans.p <- get.pvalues(ans)
  ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  
  new=generate_label(temp, flev = "Group",value.var = "rel",offset = 4, cld_info = ans.mcV)
  new$Comp=comp[i]
  stat=rbind.data.frame(stat, new)
  
}



stat=stat %>% `colnames<-`(c("Group","labels","rel","Comp"))
stat$Comp=factor(stat$Comp,levels = c("CHO","CHON","CHOS","CHONS"))
stat$Comp=factor(stat$Comp,levels = c("CHO","CHON","CHOS","CHONS"))
stat


comp_max=as.data.table(aggregate(ft_inty$rel, by=list(`Comp`=ft_inty$Comp), FUN=max))
comp_max$x=comp_max$x*1.2

pal=pal_npg("nrc")(5)


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
  ggsave(filename("Comp_compare_dunn"),height = 20, width = 30, units = "cm", dpi = 300)


ft_envi=fread("Datafile/envi_1st_sel.csv")
ft_inty$Sample=paste(ft_inty$Group,ft_inty$no, sep = "_")
ft_inty=ft_inty %>% inner_join(ft_envi)

ft_inty_eve=subset(ft_inty,ft_inty$Event!="Normal") %>% droplevels()

ft_inty_eve$Group=factor(ft_inty_eve$Group,levels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))

ggplot()+
  stat_boxplot(data=ft_inty_eve, aes(x=Group, y=rel, fill=Event),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=ft_inty_eve, aes(x=Group, y=rel, fill=Event),alpha=1, outlier.color = NA)+
  #geom_text(data =comp_max, aes(x=1, y=x, label=Comp), col="white", size=0)+
  #geom_text(data =stat, aes(x=Group, y=rel, label=labels), 
  #          col="black", size=6)+
  #facet_grid(Event~Comp, scales = "free")+
  facet_wrap(.~Comp, scales = "free_y", nrow=2, dir = "h")+
  scale_x_discrete(name="")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  #scale_fill_manual(values = npg_no_jp)+
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
        legend.position = c(0.30,0.63))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("Comp_compare_event"),height = 30, width = 30, units = "cm", dpi = 300)

####chemical properties====
FT_merge
?melt.data.table

chp_m=melt(FT_merge[,c("Sample","O.C","H.C","N.C","S.C","DBE","AI")], id.vars = c("Sample"),
           measure.vars = c("O.C","H.C","N.C","S.C","DBE","AI")) %>% 
  dcast(Sample~variable,mean) %>% `colnames<-`(c("Sample","Mean O/C","Mean H/C", "Mean N/C", "Mean S/C", "Mean DBE", "Mean AI")) %>% 
  melt(id.var=c("Sample"))

chp_m

chp_m=chp_m %>% separate(Sample, c("Group","no"),sep = "_")
chp_m$Sample=paste(chp_m$Group,chp_m$no,sep = "_")

chp_m

chp_m$Group=factor(chp_m$Group,levels = c("SUL","SS","B","NT","UL"),
                     labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
chp_m=chp_m %>% `colnames<-`(c("Group","no","variable","val","Sample"))
chp_m$Sample=paste(chp_m$Group,chp_m$no,sep = "_")
chp_m

fwrite(chp_m,"Datafile/Chemicalprop_1st.csv")


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

chp_m2=chp_m %>% inner_join(ft_envi, by=c("Sample","Group"))
chp_m2=subset(chp_m2,chp_m2$Event!="Normal") %>% droplevels
chp_m2

chp_m2$variable=factor(chp_m2$variable, levels=c("Mean O/C","Mean H/C","Mean N/C","Mean S/C","Mean DBE","Mean AI"))

chp_m2$Group=factor(chp_m2$Group, levels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))

ggplot()+
  stat_boxplot(data=chp_m2, aes(x=Group, y=val, fill=Event),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=chp_m2, aes(x=Group, y=val, fill=Event),alpha=1, outlier.color = NA)+
  #facet_grid(Event~Comp, scales = "free")+
  facet_wrap(.~variable, scales = "free_y", nrow=2, dir = "v")+
  scale_x_discrete(name="")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  #scale_fill_manual(values = npg_no_jp)+
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
        legend.position = c(0.01,0.63))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("Chp_compare_event"),height = 30, width = 45, units = "cm", dpi = 300)

####vk plot======
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


ggplot(ft_vk, aes(x=`C.`, y=`DBE`, col=Comp))+
  geom_point(size=3, alpha=0.3)+
  facet_wrap(.~Group, ncol = 5, scales = "fix")+
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
  ggsave(filename("DBE"),height = 18, width = 75, units = "cm", dpi = 700)



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


#####event vs non event=====
FT_merge

ft_envi
FT_eve=FT_merge
FT_eve$Group=factor(FT_eve$Group,levels = c("SUL","SS","B","NT","UL"),
                   labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
FT_eve$Sample=paste(FT_eve$Group,FT_eve$No,sep = "_")
FT_eve=FT_eve %>% inner_join(ft_envi)

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
MDS_1st=melt(FT_merge[,c("Sample","Formula","Bromo.Inty")], id.vars = c("Sample","Formula")) %>% 
  dcast(Sample~Formula, sum)

MDS_1st[,1:2]

MDS_1st_in=MDS_1st

NMDS_1st=metaMDS(MDS_1st_in[,-c(1)], k=5, distance = "bray", trymax = 20)
NMDS_1st ##Stress 0.08

gnmds_1st=as.data.table(NMDS_1st$points)
gnmds_1st$Sample=MDS_1st_in$Sample

gnmds_1st=gnmds_1st %>% tidyr::separate("Sample",c("Group","no"),sep="_")
gnmds_1st

gnmds_1st$Group=factor(gnmds_1st$Group, levels = c("SUL","SS","B","NT","UL"),
                       labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
gnmds_1st$Sample=paste(gnmds_1st$Group,gnmds_1st$no, sep = "_")

envi_1st=fread("Datafile/envi_1st_sel.csv")
envi_1st

gnmds_1st=gnmds_1st %>% inner_join(envi_1st)
#gnmds_1st=gnmds_1st[order(gnmds_1st$id)]
gnmds_1st

gnmds_1st$Event=factor(gnmds_1st$Event,levels = c("Event","Normal","Non-event"))

#gnmds_1st$id=ifelse(gnmds_1st$Group=="Beijing",5,
#                    ifelse(gnmds_1st$Group=="Seosan",4,
#                           ifelse(gnmds_1st$Group=="seoul",3,2)))
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
gnmds_1st

gnmds_1st_ord=gnmds_1st


gnmds_1st_ord$id=ifelse(gnmds_1st_ord$Group=="Beijing",5,
                    ifelse(gnmds_1st_ord$Group=="Seosan",4,
                           ifelse(gnmds_1st_ord$Group=="seoul",3,2)))

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


