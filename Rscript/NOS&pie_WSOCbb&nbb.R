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
library(lemon)
#devtools::install_github("cardiomoon/moonBook")
require(moonBook)
#install.packages("webr")
require(webr)
library(grid)
library(ggforce)

source("Rscript/func_filename.R")
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")
source("Rscript/piedonut_ms.R")
insert_minor_1 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_2 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}

ft_merge=fread("Datafile/FRIENDs_1st_FT.csv")
table(ft_merge$Group)
ft_merge$Group=factor(ft_merge$Group, levels =c("UL","B","SS","SUL","NT"),
                      labels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))

ft_merge=ft_merge %>% mutate(Comp = gsub('[[:digit:]]', '', Formula))
ft_merge

ft_merge$Comp=ifelse(ft_merge$`O.`==0,"Remainders",ft_merge$Comp)

ft_merge$Comp=factor(ft_merge$Comp,levels = c("CHO","CHNO","CHOS","CHNOS","Remainders"),
                     labels = c("CHO","CHON","CHOS","CHONS","Remainders"))


table(ft_merge$Comp)

ft_merge$NOS_val=ft_merge$O./(3*ft_merge$N.+3*ft_merge$S.)
ft_merge$NOS=ifelse(ft_merge$O./(3*ft_merge$N.+3*ft_merge$S.)>=1,"NOS","non")
ft_merge$NOS=ifelse(ft_merge$Comp=="CHON"&ft_merge$NOS=="NOS","ON",ft_merge$NOS)
ft_merge$NOS=ifelse(ft_merge$Comp=="CHOS"&ft_merge$NOS=="NOS","OS",ft_merge$NOS)
ft_merge$NOS=ifelse(ft_merge$Comp=="Remainders"|ft_merge$Comp=="CHO",as.character(ft_merge$Comp),ft_merge$NOS)
ft_merge

ft_merge_nos=ft_merge %>% filter(NOS%in%c("NOS","ON","OS"))



ft_merge_soa=ft_merge %>% left_join(ft_cor_pos[,c("Group","Formula","type")])

ft_merge_soa

ft_merge_soa$Freq=1

fwrite(ft_merge_soa,"ft_nos.csv")
ft_merge_soa$Sample=paste(ft_merge_soa$Group,ft_merge_soa$No, sep = "_")


temp=as.data.table(aggregate(ft_merge_soa$Freq,by=list(Sample=ft_merge_soa$Sample,Comp=ft_merge_soa$Comp, NOS=ft_merge_soa$NOS),
                             FUN=sum))



temp_tot=as.data.table(aggregate(ft_merge_soa$Freq,by=list(Sample=ft_merge_soa$Sample,Comp=ft_merge_soa$Comp),
                                 FUN=sum)) %>% `colnames<-`(c("Sample","Comp","Tot"))

temp=temp %>% inner_join(temp_tot)
temp$rel=temp$x/temp$Tot*100
temp$id=temp$Sample
temp=temp %>% separate(id,c("Group","No"))

temp_sel=temp %>% filter(!Comp%in%c("CHO","Remainders"))

table(temp_sel$Comp)
temp_sel$Comp=factor(temp_sel$Comp,levels = c("CHON","CHOS","CHONS"))
temp_sel$NOS2=factor(temp_sel$NOS,levels = c("ON","OS","NOS","non"),
                     labels = c("ON (O/N ≥ 3)","OS (O/S ≥ 3)","NOS (O ≥ 3N+3S)","Others"))
temp_sel
ggplot()+
  stat_boxplot(data=temp_sel, aes(x=Comp, y=rel/100, fill=NOS2),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=temp_sel, aes(x=Comp,y=rel/100,fill=NOS2))+
  scale_y_continuous(name="", labels = scales::percent_format(accuracy = 1L))+
  facet_rep_wrap(.~Group,repeat.tick.labels = T)+
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
        legend.position = "right")+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("1st_NOS_compare"),height = 30, width = 40, units = "cm", dpi = 300)


temp_sel_mean=melt(temp_sel[,c("Group","Comp","NOS2","rel")], 
                   id.vars = c("Group","Comp","NOS2")) %>% 
  dcast(Group+Comp~NOS2,mean) %>% melt(id.vars=c("Group","Comp"), na.rm = T) %>% 
  `colnames<-`(c("Group","Comp","NOS","mean"))
temp_sel_mean

temp_sel_sd=melt(temp_sel[,c("Group","Comp","NOS2","rel")], 
                 id.vars = c("Group","Comp","NOS2")) %>% 
  dcast(Group+Comp~NOS2,sd) %>% melt(id.vars=c("Group","Comp"), na.rm = T) %>% 
  `colnames<-`(c("Group","Comp","NOS","sd"))

temp_sel_mean
temp_sel_sd

temp_sel_avg=temp_sel_mean %>% inner_join(temp_sel_sd)
temp_sel_avg

temp_sel_stack=subset(temp_sel_avg,temp_sel_avg$NOS!="Others") %>% `colnames<-`(c("Group","Comp","NOS","error","std"))

temp_sel_stack

temp_sel_avg=temp_sel_avg %>% inner_join(temp_sel_stack[,c("Group","Comp","error")])
temp_sel_avg

temp_sel_avg$st_error=ifelse(temp_sel_avg$NOS=="Others",temp_sel_avg$mean+temp_sel_avg$error,temp_sel_avg$error)

temp_sel_avg

temp_sel_avg$Group=factor(temp_sel_avg$Group, levels =c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))

ggplot()+
  geom_bar(data=temp_sel_avg, aes(x=Comp, y=mean/100,fill=NOS), stat = "identity",position = position_stack(reverse = T))+
  geom_errorbar(data=temp_sel_avg,aes(x=Comp, ymin=(st_error-sd)/100,ymax=(st_error+sd)/100, fill=NOS), width=0.25)+
  scale_y_continuous(name="", labels = scales::percent_format(accuracy = 1L),limits = c(0,1.15), breaks=seq(0,1.1,0.25),expand = c(0.01,0.01))+
  scale_fill_manual(values = c("#EFC000FF","#008B45FF","#5F559BFF","grey50"))+
  facet_rep_wrap(Group~.,repeat.tick.labels = T, ncol=5)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.0),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 14,colour = "black", angle = 0, hjust = 0.5,vjust = 1.0),
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
        legend.direction = "horizontal",
        legend.justification=c(0.5, 0.5),
        legend.position = "bottom")+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("B&G_NOS_avg_distribution"),height = 15, width = 50, units = "cm", dpi = 300)

temp_sel$No=as.numeric(temp_sel$No)
temp_sel
temp_sel$Group=factor(temp_sel$Group, levels =c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))

ggplot()+
  geom_bar(data=temp_sel, aes(x=as.factor(No), y=rel/100,fill=NOS2), stat = "identity",position = position_stack(reverse = T))+
  scale_y_continuous(name="", labels = scales::percent_format(accuracy = 1L), expand = c(0.01,0.01))+
  scale_fill_manual(values=c("#EFC000FF","#008B45FF","#5F559BFF","grey50"))+
  #facet_rep_wrap(Group+Comp~.,repeat.tick.labels = T, ncol=3)+
  facet_rep_wrap(.~Comp+Group,repeat.tick.labels = T, ncol=5)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.0),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 14,colour = "black", angle = 0, hjust = 0.5,vjust = 1.0),
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
        legend.direction = "horizontal",
        legend.justification=c(0.5, 0.5),
        legend.position = "bottom")+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("B&G_NOS_distribution_winter"),height = 60, width = 100, units = "cm", dpi = 300)


nos=unique(ft_merge_soa[,c("Formula","NOS")])
nos


ft_cor_vk_sel2=ft_vk_m_filter_sel %>% inner_join(nos)
ft_cor_vk_sel2$Group=factor(ft_cor_vk_sel2$Group, levels =c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))

ft_cor_vk_sel2
ft_cor_vk_sel2_nos=subset(ft_cor_vk_sel2,ft_cor_vk_sel2$NOS!="non"|ft_cor_vk_sel2$NOS!="CHO"|
                            ft_cor_vk_sel2$NOS!="Remainders") %>% droplevels()

ft_cor_vk_sel2_nos=ft_cor_vk_sel2 %>% filter(NOS%in%c("ON","OS","NOS"))
ft_cor_vk_sel2_nos
ft_cor_vk_sel2_nos$NOS2=factor(ft_cor_vk_sel2_nos$NOS, levels = c("ON","OS","NOS"),
                               labels = c("ON (O/N ≥ 3)","OS (O/S ≥ 3)","NOS (O ≥ 3N+3S)"))

ft_cor_vk_sel2_non=ft_cor_vk_sel2 %>% filter(!NOS%in%c("ON","OS","NOS"))

ggplot()+
  geom_point(data=ft_cor_vk_sel2_non,aes(x=`O/C`, y=`H/C`),size=1.5,col="grey50", alpha=0.4)+
  geom_point(data=ft_cor_vk_sel2_nos,aes(x=`O/C`, y=`H/C`,col=NOS2),size=2, alpha=0.4)+
  facet_rep_grid(.~Group, repeat.tick.labels = "all",scales = "free")+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_manual(values = c("#EFC000FF","#008B45FF","#5F559BFF"))+
  scale_size_manual(values = c(0.1,1,2,4,10),
                    breaks =  c(3, 4, 5, 6,7),
                    labels = c(expression(bold("1 x 10"^"3")), expression(bold("1 x 10"^"4")), expression(bold("1 x 10"^"5")),
                               expression(bold("1 x 10"^"6")),expression(bold("1 x 10"^"7"))))+
  #scale_size(range = c(1,8),
  #           breaks = 1000 * c(1, 10, 100, 1000,10000),
  #           labels = c(expression(bold("1 x 10"^"3")), expression(bold("1 x 10"^"4")), expression(bold("1 x 10"^"5")),
  #                      expression(bold("1 x 10"^"6")),expression(bold("1 x 10"^"7"))))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.3,0,0.3,0),"cm")),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.3,0.0,0.3),"cm")),
        axis.title.x = element_text(size = 40, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0,0.0,0),"cm")),
        axis.title.y = element_text(size = 40, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.4,0.0,0.0),"cm")),
        strip.text.x = element_text(size = 44, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.3,0,0.5,0),"cm")),
        strip.text.y = element_text(size = 44, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.3,0.2,0.7),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 40, colour = "black", family = "Arial",hjust = 0.5,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
        #legend.position = "right",
        #legend.direction = "vertical",
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  guides(col=guide_legend(title = "",ncol = 5,override.aes = list(shape=21,size=14, fill=c("#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), 
         size="none")+
  ggsave(filename("vk_avg_nos"),height = 20, width = 100, units = "cm", dpi = 300)


ft_cor_vk_sel2_nos

table(ft_cor_vk_sel2_nos$variable)
ft_cor_vk_sel2$Freq=1

ft_vk_m_filter_sel

ft_cor_pos

fm_nos_list=ft_merge_soa[,c("Formula","Comp","NOS","Freq")]

fm_nos_list=unique(fm_nos_list)

ft_cor_pos2=ft_cor_pos %>% inner_join(fm_nos_list)

fm_wsoc=ft_cor_pos2[,c("Group","Formula","type","Comp","NOS","Freq")]


t=c("#BC3C29FF","#EFC000","#008B45","#5F559B","grey70","grey30")

fm_wsoc$variable2=factor(fm_wsoc$type, level=c("WSOCbb","WSOCnbb"),
                         labels=c("WSOC ","WSOC  "))


fm_wsoc$varlab=factor(fm_wsoc$type,levels = c("WSOCbb","WSOCnbb"),
                      labels = c(expression(bold("WSOC"["BB"])),
                                 expression(bold("WSOC"["NBB"]))))

table(fm_wsoc$NOS)
fm_wsoc$NOS2=factor(fm_wsoc$NOS, levels = c("CHO","ON","OS","NOS","non","Remainders"),
                    labels = c("CHO","ON","OS","NOS","NF","Remainders"))

fm_wsoc_b=subset(fm_wsoc,fm_wsoc$Group=="Beijing")

PieDonut_ms(fm_wsoc_b,aes(variable2 ,NOS2), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,showPieName = F,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill="manual",fillcol = c("#9DBCD4","#CB7723"),
            subfill = "manual", subfillcol=rep(t,2))


png("Figure/B_pie.png", width = 20, height = 20,units = "cm",res = 300, bg = "transparent")
PieDonut_ms(fm_wsoc_b,aes(variable2 ,NOS2), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,showPieName = F,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill="manual",fillcol = c("#9DBCD4","#CB7723"),
            subfill = "manual", subfillcol=rep(t,2))
dev.off()


fm_wsoc_u=subset(fm_wsoc,fm_wsoc$Group=="Ulaanbaatar")
fm_wsoc_u

png("Figure/UL_pie.png", width = 20, height = 20,units = "cm",res = 300, bg = "transparent")
PieDonut_ms(fm_wsoc_u,aes(variable2 ,NOS2), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,showPieName = F,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill="manual",fillcol = c("#9DBCD4","#CB7723"),
            subfill = "manual", subfillcol=rep(t,2))
dev.off()



fm_wsoc_s=subset(fm_wsoc,fm_wsoc$Group=="Seoul")

png("Figure/Sul_pie.png", width = 20, height = 20,units = "cm",res = 300,bg = "transparent")
PieDonut_ms(fm_wsoc_s,aes(variable2 ,NOS2), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,showPieName = F,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill="manual",fillcol = c("#9DBCD4","#CB7723"),
            subfill = "manual", subfillcol=rep(t,2))
dev.off()


fm_wsoc_ss=subset(fm_wsoc,fm_wsoc$Group=="Seosan")

png("Figure/220506ss_pie.png", width = 20, height = 20,units = "cm",res = 300,bg = "transparent")
PieDonut_ms(fm_wsoc_ss,aes(variable2 ,NOS2), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,showPieName = F,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill="manual",fillcol = c("#9DBCD4","#CB7723"),
            subfill = "manual", subfillcol=rep(t,2))
dev.off()

fm_wsoc_nt=subset(fm_wsoc,fm_wsoc$Group=="Noto")

png("Figure/220506nt_pie.png", width = 20, height = 20,units = "cm",res = 300,bg = "transparent")
PieDonut_ms(fm_wsoc_nt,aes(variable2 ,NOS2), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,showPieName = F,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill="manual",fillcol = c("#9DBCD4","#CB7723"),
            subfill = "manual", subfillcol=rep(t,2))
dev.off()
