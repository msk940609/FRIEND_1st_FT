library(ggplot2)
library(data.table)
library(dplyr)
library(tidyverse)
options(java.parameters = "- Xmx8192m")
library(extrafont)
loadfonts(device="win")
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


source("Rscript/func_filename.R")
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")
insert_minor_1 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_2 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}

FT_merge=fread("Datafile/FRIENDs_1st_FT.csv")
FT_merge$Sample=paste(FT_merge$Group,FT_merge$No,sep = "_")

FT_merge$C=ifelse(FT_merge$C.>0, "C","")
FT_merge$H=ifelse(FT_merge$H.>0, "H","")
FT_merge$O=ifelse(FT_merge$O.>0, "O","")
FT_merge$N=ifelse(FT_merge$N.>0, "N","")
FT_merge$S=ifelse(FT_merge$S.>0, "S","")

FT_merge=FT_merge %>% unite("Comp",c("C","H","O","N","S"),sep = "")
FT_merge$Comp=ifelse(FT_merge$O.==0, "Remainders",FT_merge$Comp)
#FT_merge=subset(FT_merge,FT_merge$Comp!="Remainders") %>% droplevels()

FT_bb=fread("Datafile/WOSCbb_S&B.csv")

FT_merge_sel=FT_merge %>% filter(Group%in%c("SUL","B"))
FT_merge_sel$Group=ifelse(FT_merge_sel$Group=="SUL","Seoul","Beijing")
FT_merge_sel$Sample=paste(FT_merge_sel$Group,FT_merge_sel$No,sep = "_")

FT_merge_mono=(aggregate(FT_merge_sel$Mono.Inty, by=list(`Sample`=FT_merge_sel$Sample), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Tot.mono"))

FT_merge_bromo=(aggregate(FT_merge_sel$Bromo.Inty, by=list(`Sample`=FT_merge_sel$Sample), FUN=sum))%>% 
  `colnames<-`(c("Sample","Tot.bromo"))

FT_merge_sel=FT_merge_sel %>% inner_join(FT_merge_mono)
FT_merge_sel=FT_merge_sel %>% inner_join(FT_merge_bromo)

FT_merge_sel$monorel=FT_merge_sel$Mono.Inty/FT_merge_sel$Tot.mono*100
FT_merge_sel$bromorel=FT_merge_sel$Bromo.Inty/FT_merge_sel$Tot.bromo*100
#tt=(aggregate(FT_merge_sel$rel, by=list(`Sample`=FT_merge_sel$Sample), FUN=sum))

FT_merge_sel=FT_merge_sel %>% inner_join(FT_bb[,c("Sample","WSOCbb","WSOCnbb","POA","SOA")])

FT_m=melt(FT_merge_sel[,c("Sample","Group","No","WSOCbb","WSOCnbb","Formula","Bromo.Inty")],
          id.vars = c("Sample","Group","No","WSOCbb","WSOCnbb","Formula")) %>% 
  dcast(Sample+Group+No+WSOCbb+WSOCnbb~Formula, value.var = "value",sum)

FT_m[,1:23]

fwrite(FT_m, file = "Datafile/FT_formula_bromo.csv")

cor_SOA=fread("Datafile/FT_formula_ob3_mono_rel_ex.csv")
#cor_SOA=fread("Datafile/FT_formula_ob3_mono_rel.csv")
cor_SOA[,1:23]

cor_temp=cor_SOA
grp=unique(cor_temp$Group)

df=data.frame()
for (i in 1:length(grp)) {
  #i=1
  temp=as.data.frame(subset(cor_temp,cor_temp$Group==grp[i]))
  for (j in 4:5) {
    #j=4
    for (k in 6:dim(temp)[2]) {
      #k=8  
      sum.j<-sum(temp[,j])
      sum.k<-sum(temp[,k])
      if(sum.j >0 & sum.k >0){
        test<-cor.test(temp[,j],temp[,k],method="spearman",na.action=na.rm,exact = F)
        rho<-test$estimate
        p.value<-test$p.value
      }
      
      if(sum.j <=0 | sum.k <= 0){
        rho<-0
        p.value<-1
      }	
      
      new=data.frame("Envi"=names(temp)[j],"Formula"=names(temp)[k],"rho"=round(test$estimate,3),
                     "p"=round(test$p.value,4), "Group"=grp[i])
      
      df<-rbind(df,new)			
      
    }
    print(i/length(grp))
    
  }
  fwrite(df, file = "cortest.csv")
}



#Cor specialization=====
#cor_vk=fread("cortest_bromo_rel.csv")
cor_vk=fread("cortest.csv")

#cor_vk=subset(cor_vk,cor_vk$p<0.05)
cor_vk

CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(cor_vk$Formula)

# extract numerical formulae
numericalFormula <- array(0, dim=c(length(molecularFormula), length(CHEMICAL_ELEMENTS)))
for (k in 1:length(molecularFormula)){
  #k=1
  formula <- molecularFormula[k]
  ge <- gregexpr("[A-Z]\\d*", formula, perl=TRUE)
  s_index <- ge[[1]]
  s_len <- attr(s_index, "match.length")
  for (i in 1:length(s_len)){
    #i=1
    token <- substr(formula, s_index[i], s_index[i] + s_len[i] - 1)
    element <- substr(token, 1, 1)
    if (grepl(element, "CHNOSP")) {
      idx = which(CHEMICAL_ELEMENTS %in% element)
      if (numericalFormula[k, idx] > 0) next  # for C13
      if (s_len[i] == 1) {
        numericalFormula[k, idx] = 1
      } else {
        numElement <- try(strtoi(substr(formula, s_index[i] + 1, s_index[i] + s_len[i] - 1)))
        if (class(numElement)=="integer"){
          numericalFormula[k, idx] = numElement
        } else {
          print(paste("[ERROR] an unknown chemical element found:", token, "in", formula))
        }
      }
    } else {
      print(paste("[ERROR] an unknown chemical element found:", element, "in", formula))
    }
  }
}

fm=as.data.table(numericalFormula) %>% `colnames<-`(c("C","H","N","O","S"))
fm$Formula=molecularFormula

fm

cor_vk_d=dcast(cor_vk, Group+Formula~Envi, sum, value.var = c("rho","p"))
cor_vk_d

cor_vk_d=cor_vk_d %>% inner_join(fm)
cor_vk_d$O.C=cor_vk_d$O/cor_vk_d$C
cor_vk_d$H.C=cor_vk_d$H/cor_vk_d$C

cor_vk_m=melt(cor_vk_d[,c("Group","Formula","O.C","H.C","rho_WSOCbb","rho_WSOCnbb","p_WSOCbb","p_WSOCnbb")],
                 id.vars=c("Group","Formula","O.C","H.C","rho_WSOCbb","rho_WSOCnbb"),
              variable.name = "p-var",value.name = "p") %>% 
  melt(id.vars=c("Group","Formula","O.C","H.C","p-var","p"))
cor_vk_m
cor_vk_m_sel=subset(cor_vk_m,abs(cor_vk_m$value)>0.001)

##specialization====
cor_vk_d$Type=ifelse(abs(cor_vk_d$rho_WSOCbb)>abs(cor_vk_d$rho_WSOCnbb),"BB","NBB")
#cor_vk_d$SOA=ifelse(abs(cor_vk_d$rho_SOA)>abs(cor_vk_d$rho_POA),"SOA","POA")
#cor_vk_ev$Type=ifelse(abs(abs(cor_vk_ev$rho_WSOCbb)-abs(cor_vk_ev$rho_WSOCnbb))<0.05,"Both",cor_vk_ev$Type)
table(cor_vk_d$Type)
#table(cor_vk_d$SOA)
cor_vk_d

cor_vk_d$variable=ifelse(cor_vk_d$Type=="BB","WSOCbb","WSOCnbb")
cor_vk_d$value=ifelse(cor_vk_d$Type=="BB",cor_vk_d$rho_WSOCbb,cor_vk_d$rho_WSOCnbb)
cor_vk_d$p=ifelse(cor_vk_d$Type=="BB",cor_vk_d$p_WSOCbb,cor_vk_d$p_WSOCnbb)
cor_vk_d

cor_vk_d_sel=cor_vk_d[,c("Group","Formula","O.C","H.C","Type","variable","value","p")]
cor_vk_d_sel=unique(cor_vk_d_sel)
cor_vk_d_sel
#cor_vk_d_sel=subset(cor_vk_d_sel,abs(cor_vk_d_sel$value)>0.05)

dim(cor_vk_d_sel)
cor_vk_d_sel_s=subset(cor_vk_d_sel,cor_vk_d_sel$Group=="Seoul")
cor_vk_d_sel_s

table(cor_vk_d_sel_s$Type)
table(cor_vk_d_sel_s$variable)

table(cor_vk_d_sel_s$Type)

cor_vk_sel_s1=subset(cor_vk_d_sel_s,abs(cor_vk_d_sel_s$value)>0.3)
#cor_vk_sel_s1=subset(cor_vk_d_sel_s,abs(cor_vk_d_sel_s$p)<0.05)
table(cor_vk_sel_s1$Type)

cor_vk_sel_s1$Grouplab=cor_vk_sel_s1$Group
cor_vk_sel_s1$Grouplab=factor(cor_vk_sel_s1$Grouplab,levels = c("Seoul","Beijing"),
                                 labels=c(expression(bold("Seoul")),
                                          expression(bold("Beijing"))))

cor_vk_sel_s1$varlab=cor_vk_sel_s1$variable
cor_vk_sel_s1$varlab=factor(cor_vk_sel_s1$varlab,levels = c("WSOCbb","WSOCnbb"),
                               labels = c(expression(bold("WSOC"["bb"])),
                                          expression(bold("WSOC"["nbb"]))))

cor_vk_sel_s2=subset(cor_vk_sel_s1,abs(cor_vk_sel_s1$value)>0.30)


tt=as.data.frame(table(cor_vk_sel_s2$Type)) 
tt$varlab=tt$Var1
tt$varlab=factor(tt$varlab,levels = c("BB","NBB"),
                 labels = c(expression(bold("WSOC"["bb"])),
                            expression(bold("WSOC"["nbb"]))))
tt$O.C=0.1
tt$H.C=0.15


ggplot()+
  geom_point(data=cor_vk_sel_s2, aes(x=O.C, y=H.C,col=value),size=2.7)+
  #geom_text(data=tt, aes(x=O.C, y=H.C,label=Freq), size=8)+
  geom_text(data=tt, aes(x=0.05, y=0.14,label=paste0(expression(italic("n=")))),parse=T, size=8)+
  geom_text(data=tt, aes(x=0.15, y=H.C,label=Freq),parse=T, size=8)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #annotate(x=0.1,y=0.25, geom = "text", label=32)+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  #scale_color_gradientn(colors=(topo.colors(20)[4:17]))+
  scale_color_gradientn(colors=(topo.colors(40)[4:40]))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.3,0,0.3,0),"cm")),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.3,0.0,0.3),"cm")),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.text.x = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.6,0,0.4,0),"cm")),
        strip.text.y = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial",hjust = 0.5,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26,hjust = 0.2 ,vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
        legend.position = "right",
        legend.direction = "vertical",
        #legend.position = "bottom",
        #legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  guides(col=guide_colorbar(title = expression(bolditalic("ρ")),ticks.colour = "black", ticks.linewidth = 1.5,
                            barwidth = 2.5, barheight = 25))+
  ggsave(filename("vk_cor_s_rel"),height = 20, width = 40, units = "cm", dpi = 300)




cor_vk_d_sel_b=subset(cor_vk_d_sel,cor_vk_d_sel$Group=="Beijing")
cor_vk_d_sel_b

cor_vk_d_sel_b=subset(cor_vk_d_sel_b,abs(cor_vk_d_sel_b$value)>0.0)
cor_vk_d_sel_b1=subset(cor_vk_d_sel_b,abs(cor_vk_d_sel_b$p)<0.05)

table(cor_vk_d_sel_b1$Type)

cor_vk_d_sel_b1$Grouplab=cor_vk_d_sel_b1$Group
cor_vk_d_sel_b1$Grouplab=factor(cor_vk_d_sel_b1$Grouplab,levels = c("Seoul","Beijing"),
                                 labels=c(expression(bold("Seoul")),
                                          expression(bold("Beijing"))))

cor_vk_d_sel_b1$varlab=cor_vk_d_sel_b1$variable
cor_vk_d_sel_b1$varlab=factor(cor_vk_d_sel_b1$varlab,levels = c("WSOCbb","WSOCnbb"),
                               labels = c(expression(bold("WSOC"["bb"])),
                                          expression(bold("WSOC"["nbb"]))))
cor_vk_d_sel_b1

cor_vk_d_sel_b2=subset(cor_vk_d_sel_b1,abs(cor_vk_d_sel_b1$value)>0.30)
cor_vk_d_sel_b2=subset(cor_vk_d_sel_b2,abs(cor_vk_d_sel_b2$p)<0.05)


tt2=as.data.frame(table(cor_vk_d_sel_b2$Type)) 
tt2$varlab=tt2$Var1
tt2$varlab=factor(tt2$varlab,levels = c("BB","NBB"),
                 labels = c(expression(bold("WSOC"["bb"])),
                            expression(bold("WSOC"["nbb"]))))
tt2$O.C=0.1
tt2$H.C=0.15


ggplot()+
  geom_point(data=cor_vk_d_sel_b2, aes(x=O.C, y=H.C,col=value),size=2.7)+
  geom_text(data=tt2, aes(x=0.05, y=0.14,label=paste0(expression(italic("n=")))),parse=T, size=8)+
  geom_text(data=tt2, aes(x=0.15, y=H.C,label=Freq),parse=T, size=8)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  scale_color_gradientn(colors=(topo.colors(30)[2:28]))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:32]))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.3,0,0.3,0),"cm")),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.3,0.0,0.3),"cm")),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.text.x = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.6,0,0.4,0),"cm")),
        strip.text.y = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial",hjust = 0.5,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26,hjust = 0.2 ,vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
        legend.position = "right",
        legend.direction = "vertical",
        #legend.position = "bottom",
        #legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  guides(col=guide_colorbar(title = expression(bolditalic("ρ")),ticks.colour = "black", ticks.linewidth = 1.5,
                            barwidth = 2.5, barheight = 25))+
  ggsave(filename("vk_cor_b_rel"),height = 20, width = 40, units = "cm", dpi = 300)


##display selected vk to event/non-event====
##Seoul====
FT_merge_sel
cor_vk_sel_s1

FT_merge_s_ev=subset(FT_merge_sel,FT_merge_sel$Sample=="Seoul_9")

cor_vk_sel_s1
cor_fm_s1=dcast(cor_vk_sel_s1[,c("Formula","variable","value")],Formula~variable, mean) %>% 
  `colnames<-`(c("Formula","rho_WSOCbb","rho_WSOCnbb"))
cor_fm_s1

FT_merge_s_ev=FT_merge_s_ev %>% left_join(cor_fm_s1,na_matches = "na" )

FT_merge_s_ev

FT_s_ev_vk=FT_merge_s_ev[,c("Formula","O.C","H.C","rho_WSOCbb","rho_WSOCnbb")]
FT_s_ev_vk_m=melt(FT_s_ev_vk,id.vars = c("Formula","O.C","H.C"))


FT_s_ev_vk_m=FT_s_ev_vk_m[rev(order(FT_s_ev_vk_m$value,decreasing = T)),]

ggplot()+
  geom_point(data=FT_s_ev_vk_m, aes(x=O.C, y=H.C,col=value),size=2.7)+
  #geom_text(data=tt, aes(x=O.C, y=H.C,label=Freq), size=8)+
  facet_grid(variable~., labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #annotate(x=0.1,y=0.25, geom = "text", label=32)+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  #scale_color_gradientn(colors=(topo.colors(20)[4:17]))+
  scale_color_gradientn(colors=(topo.colors(40)[4:40]),na.value = "grey85")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.3,0,0.3,0),"cm")),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.3,0.0,0.3),"cm")),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.text.x = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.6,0,0.4,0),"cm")),
        strip.text.y = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial",hjust = 0.5,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26,hjust = 0.2 ,vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
        legend.position = "right",
        legend.direction = "vertical",
        #legend.position = "bottom",
        #legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  guides(col=guide_colorbar(title = expression(bolditalic("ρ")),ticks.colour = "black", ticks.linewidth = 1.5,
                            barwidth = 2.5, barheight = 25))+
  ggsave(filename("vk_cor_s_ev"),height = 30, width = 25, units = "cm", dpi = 300)



FT_merge_s_ev=FT_merge_sel %>% filter(Sample%in%c("Seoul_25","Seoul_9"))

cor_vk_sel_s1
cor_fm_s1=dcast(cor_vk_sel_s1[,c("Formula","variable","value")],Formula~variable, mean) %>% 
  `colnames<-`(c("Formula","rho_WSOCbb","rho_WSOCnbb"))

FT_merge_s_ev=FT_merge_s_ev %>% left_join(cor_fm_s1,na_matches = "na" )

FT_merge_s_ev$event=ifelse(FT_merge_s_ev$Sample=="Seoul_25","Non-event","Event")
FT_s_ev_vk=FT_merge_s_ev[,c("Sample","Formula","O.C","H.C","rho_WSOCbb","rho_WSOCnbb","event")]

FT_s_ev_vk_m=melt(FT_s_ev_vk,id.vars = c("Sample","Formula","event","O.C","H.C"))

FT_s_ev_vk_m$varlab=FT_s_ev_vk_m$variable

FT_s_ev_vk_m=FT_s_ev_vk_m[rev(order(FT_s_ev_vk_m$value,decreasing = T)),]

FT_s_ev_vk_m$varlab=FT_s_ev_vk_m$variable
FT_s_ev_vk_m$varlab=factor(FT_s_ev_vk_m$varlab,levels = c("rho_WSOCbb","rho_WSOCnbb"),
                              labels = c(expression(bold("WSOC"["bb"])),
                                         expression(bold("WSOC"["nbb"]))))
FT_s_ev_vk_m$evlab=FT_s_ev_vk_m$event
FT_s_ev_vk_m$evlab=factor(FT_s_ev_vk_m$evlab,levels = c("Event","Non-event"),
                           labels = c(expression(bold("Event (12/23)" )),
                                      expression(bold("Non-event (1/5)"))))

fwrite(FT_s_ev_vk_m,file = "vk_s_ev&nev.csv")

ggplot()+
  geom_point(data=FT_s_ev_vk_m, aes(x=O.C, y=H.C,col=value),size=2.7)+
  #geom_text(data=tt, aes(x=O.C, y=H.C,label=Freq), size=8)+
  facet_grid(varlab~evlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #annotate(x=0.1,y=0.25, geom = "text", label=32)+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  #scale_color_gradientn(colors=(topo.colors(20)[4:17]))+
  scale_color_gradientn(colors=(topo.colors(40)[6:37]),na.value = "grey85")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.3,0,0.3,0),"cm")),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.3,0.0,0.3),"cm")),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.text.x = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.6,0,0.4,0),"cm")),
        strip.text.y = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial",hjust = 0.5,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26,hjust = 0.2 ,vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
        #legend.position = "right",
        #legend.direction = "vertical",
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  guides(col=guide_colorbar(title = expression(bolditalic("ρ")),ticks.colour = "black", ticks.linewidth = 1.5,
                            barwidth = 30, barheight = 2.5))+
  ggsave(filename("vk_cor_s_ev&nev"),height = 33, width = 30, units = "cm", dpi = 300)


##Beijing=====
cor_vk_d_sel_b1
FT_merge_b_ev=FT_merge_sel %>% filter(Sample%in%c("Beijing_13","Beijing_25"))

cor_vk_d_sel_b
cor_fm_b1=dcast(cor_vk_d_sel_b1[,c("Formula","variable","value")],Formula~variable, mean) %>% 
  `colnames<-`(c("Formula","rho_WSOCbb","rho_WSOCnbb"))
cor_fm_b1
FT_merge_b_ev=FT_merge_b_ev %>% left_join(cor_fm_b1,na_matches = "na" )

FT_merge_b_ev$event=ifelse(FT_merge_b_ev$Sample=="Beijing_25","Non-event","Event")
FT_b_ev_vk=FT_merge_b_ev[,c("Sample","Formula","O.C","H.C","rho_WSOCbb","rho_WSOCnbb","event")]

FT_b_ev_vk_m=melt(FT_b_ev_vk,id.vars = c("Sample","Formula","event","O.C","H.C"))

FT_b_ev_vk_m$varlab=FT_b_ev_vk_m$variable

FT_b_ev_vk_m=FT_b_ev_vk_m[rev(order(FT_b_ev_vk_m$value,decreasing = T)),]

FT_b_ev_vk_m$varlab=FT_b_ev_vk_m$variable
FT_b_ev_vk_m$varlab=factor(FT_b_ev_vk_m$varlab,levels = c("rho_WSOCbb","rho_WSOCnbb"),
                           labels = c(expression(bold("WSOC"["bb"])),
                                      expression(bold("WSOC"["nbb"]))))
FT_b_ev_vk_m$evlab=FT_b_ev_vk_m$event
FT_b_ev_vk_m$evlab=factor(FT_b_ev_vk_m$evlab,levels = c("Event","Non-event"),
                          labels = c(expression(bold("Event (12/27)")),
                                     expression(bold("Non-event (1/5)"))))

fwrite(FT_b_ev_vk_m,file = "vk_b_ev&nev.csv")

ggplot()+
  geom_point(data=FT_b_ev_vk_m, aes(x=O.C, y=H.C,col=value),size=2.7)+
  #geom_text(data=tt, aes(x=O.C, y=H.C,label=Freq), size=8)+
  facet_grid(varlab~evlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #annotate(x=0.1,y=0.25, geom = "text", label=32)+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  #scale_color_gradientn(colors=(topo.colors(20)[4:17]))+
  scale_color_gradientn(colors=topo.colors(30)[2:28],na.value = "grey85")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.3,0,0.3,0),"cm")),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.3,0.0,0.3),"cm")),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.text.x = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.6,0,0.4,0),"cm")),
        strip.text.y = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial",hjust = 0.5,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26,hjust = 0.2 ,vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
        #legend.position = "right",
        #legend.direction = "vertical",
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  guides(col=guide_colorbar(title = expression(bolditalic("ρ")),ticks.colour = "black", ticks.linewidth = 1.5,
                            barwidth = 30, barheight = 2.5))+
  ggsave(filename("vk_cor_b_ev&nev2"),height = 33, width = 30, units = "cm", dpi = 300)


###clustering=====
cor_vk_sel_s2
cor_vk_d_sel_b2

FT_merge_sel

cor_vk_sel_s2

cl_s=cor_vk_sel_s2[,c("Group","Formula","variable","value")]
cl_s$class=ifelse(cl_s$value>0,"pos","neg")
cl_s$Cl=paste(cl_s$variable,cl_s$class,sep = "_")

cl_s=cl_s %>% left_join(FT_merge_sel[,c("Formula","O.C","H.C","DBE","AI","OSc","Comp")])
cl_s=unique(cl_s)
cl_s$Freq=1

#fwrite(cl_s, file = "SUL_WSOC_corlist.csv")
#fwrite(cl_b, file = "beijing_WSOC_corlist.csv")


cl_s_tot=(aggregate(cl_s$Freq, by=list(`Cl`=cl_s$Cl), FUN=sum)) 

cl_s_d=dcast(cl_s[,c("Cl","Comp","Freq")],Cl~Comp,sum) %>% 
  melt(id.vars=c("Cl"))
cl_s_d=cl_s_d %>% inner_join(cl_s_tot)

cl_s_d$rel=round(cl_s_d$value/cl_s_d$x*100,1)

cl_s_d$variable=factor(cl_s_d$variable,levels = c("CHO","CHON","CHOS","CHONS"))

cl_s_d$Cl=factor(cl_s_d$Cl,levels = c("WSOCbb_pos","WSOCnbb_pos","WSOCbb_neg","WSOCnbb_neg"))


cl_s_d=subset(cl_s_d,cl_s_d$variable!="Remainders") %>% droplevels()

ggplot()+
  geom_bar(data=cl_s_d,aes(x="",y=value, fill=variable), col="black",stat="identity", size=1.5,
           position = position_fill(reverse = T))+
  coord_polar("y")+
  scale_fill_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF"))+
  facet_wrap(.~Cl, ncol = 4)+
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
  ggsave("Seoul_cluster_pie.png",height = 10, width = 40, units = "cm", dpi = 300)

cor_vk_d_sel_b2

cl_b=cor_vk_d_sel_b2[,c("Group","Formula","variable","value")]
cl_b$class=ifelse(cl_b$value>0,"pos","neg")
cl_b$Cl=paste(cl_b$variable,cl_b$class,sep = "_")

cl_b=cl_b %>% left_join(FT_merge_sel[,c("Formula","O.C","H.C","DBE","AI","OSc","Comp")])

cl_b=unique(cl_b)
cl_b$Freq=1

fwrite(cl_b,"test.csv")


cl_b_tot=(aggregate(cl_b$Freq, by=list(`Cl`=cl_b$Cl), FUN=sum)) 

cl_b_d=dcast(cl_b[,c("Cl","Comp","Freq")],Cl~Comp,sum) %>% 
  melt(id.vars=c("Cl"))


cl_b_d=cl_b_d %>% inner_join(cl_b_tot)

cl_b_d$rel=round(cl_b_d$value/cl_b_d$x*100,1)
cl_b_d
cl_b_d=subset(cl_b_d,cl_b_d$variable!="Remainders") %>% droplevels()

cl_b_d$variable=factor(cl_b_d$variable,levels = c("CHO","CHON","CHOS","CHONS"))

cl_b_d$Cl=factor(cl_b_d$Cl,levels = c("WSOCbb_pos","WSOCnbb_pos","WSOCbb_neg","WSOCnbb_neg"))

ggplot()+
  geom_bar(data=cl_b_d,aes(x="",y=value, fill=variable), col="black",stat="identity", size=0.5,
           position = position_fill(reverse = T))+
  coord_polar("y")+
  scale_fill_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF"))+
  facet_wrap(.~Cl, ncol = 4)+
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
  theme(legend.text = element_text(size=20, hjust=0.5),
        legend.position = "bottom",
        legend.direction = "horizontal")+
  ggsave("Beijing_cluster_pie.png",height = 10, width = 40, units = "cm", dpi = 300)

##chemical properties comparison=====

cl_c=rbind(cl_s,cl_b)
cl_c_sel=cl_c[,c("Group","Cl","O.C","H.C","DBE","AI")]
cl_c_sel_m=melt(cl_c_sel,id.vars = c("Group","Cl"))
  
cl_c_sel_m
cl_c_sel_m$Group=factor(cl_c_sel_m$Group, levels = c("Seoul","Beijing"))
cl_c_sel_m$variable=factor(cl_c_sel_m$variable, levels = c("DBE","AI","H.C","O.C"),
                           labels = c("DBE","AI","H/C","O/C"))

cl_c_sel_m$Cl=factor(cl_c_sel_m$Cl, levels = c("WSOCbb_pos","WSOCnbb_pos","WSOCbb_neg","WSOCnbb_neg"))

cl_c_sel_m$Clvar=cl_c_sel_m$Cl
cl_c_sel_m$Clvar=factor(cl_c_sel_m$Clvar, levels = c("WSOCbb_pos","WSOCnbb_pos","WSOCbb_neg","WSOCnbb_neg"),
                        labels = c(expression(bold("WSOC"["bb_pos"])),
                                   expression(bold("WSOC"["nbb_pos"])),
                                   expression(bold("WSOC"["bb_neg"])),
                                   expression(bold("WSOC"["nbb_neg"]))))

cl_max=as.data.table(aggregate(cl_c_sel_m$value, by=list("variable"=cl_c_sel_m$variable, Cl=
                                                           cl_c_sel_m$Cl), FUN=max))

cl_max$x=ifelse(cl_max$variable=="DBE",25,cl_max$x)

ggplot(cl_c_sel_m,aes(x=Group, y=value, fill=Cl))+
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.85), width=0.4)+
  geom_boxplot(position = position_dodge(width = 0.85), outlier.color = NA)+
  geom_blank(data=cl_max, aes(x=2.5,y=x*1.2))+
  scale_fill_lancet(labels=c(expression(bold("WSOC"["bb_pos"])),
                           expression(bold("WSOC"["nbb_pos"])),
                           expression(bold("WSOC"["bb_neg"])),
                           expression(bold("WSOC"["nbb_neg"]))))+
  facet_wrap(.~variable,scales = "free", strip.position = "left")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.3,0,0.3,0),"cm")),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.3,0.0,0.3),"cm")),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        strip.text.x = element_text(size = 25, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.6,0,0.4,0),"cm")),
        strip.text.y = element_text(size = 25, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.text = element_text(size = 22, colour = "black", family = "Arial",hjust = 0.5,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26,hjust = 0.2 ,vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
        legend.key.height = unit(2,"cm"),
        legend.key.width = unit(2,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  guides(fill=guide_legend(title = "", override.aes = list(size=1)))+
  ggsave(filename("cl_comp"),height = 30, width = 40, units = "cm", dpi = 300)

cl_c_sel_m_s=subset(cl_c_sel_m,cl_c_sel_m$Group=="Seoul")
cl_c_sel_m_b=subset(cl_c_sel_m,cl_c_sel_m$Group=="Beijing")


cl_c_sel_m_b$Region=cl_c_sel_m_b$Group
cl_c_sel_m_b$Group=cl_c_sel_m_b$Cl
chp=unique(cl_c_sel_m$variable)
stat_par=data.table()
for (i in 1:length(chp)) {
  temp=subset(cl_c_sel_m_b,cl_c_sel_m_b$variable==chp[i])
  
  tm=as.data.table(aggregate(temp$value, by=list(`variable`=temp$Cl), FUN=median))
  tm=tm[order(tm$x, decreasing = T),]
  
  temp$Cl=factor(temp$Cl,levels = tm$variable)
  #ans <- posthoc.kruskal.dunn.test(rel ~ Group,data=temp ,p.adjust="none") ##
  #ans.p <- get.pvalues(ans)
  #ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  ans=dunnettT3Test(value ~ Cl,data=temp, p.adjust.methods="fdr")
  #ans=lsdTest(val ~ Group,data=temp, p.adjust.methods="fdr")
  ans.p <- get.pvalues(ans)
  ans.mcV <- multcompLetters(ans.p, threshold=0.05)
  
  new=generate_label(temp, flev = "Group",value.var = "value",offset = 0.0, cld_info = ans.mcV)
  new$chp=chp[i]
  stat_par=rbind.data.frame(stat_par, new)
  
}

stat_sul

stat_b

###select pos & compare all=====
cl_c=rbind(cl_s,cl_b)
cl_c_sel=cl_c[,c("Group","Cl","O.C","H.C","DBE","AI")]

FT_merge_sel
tt=FT_merge_sel
tt$Cl="All"

tt=FT_merge_sel
tt$Cl="All"
tt_m=melt(tt[,c("Group","O.C","H.C","DBE","AI","Formula","Cl")],id.vars = c("Group","Formula","Cl")) %>% 
  dcast(Group+Formula+Cl~variable, value.var = "value",mean)
tt_m

cl_c_sel=cl_c[,c("Group","Cl","O.C","H.C","DBE","AI")]
cl_c_sel_p=cl_c_sel %>% filter(Cl%in%c("WSOCbb_pos","WSOCnbb_pos"))

cl_a_sel=rbind(cl_c_sel_p,tt_m[,-"Formula"])
cl_a_sel_m=melt(cl_a_sel,id.vars = c("Group","Cl"))

cl_a_sel_m
cl_a_sel_m$Group=factor(cl_a_sel_m$Group, levels = c("Seoul","Beijing"))
cl_a_sel_m$variable=factor(cl_a_sel_m$variable, levels = c("DBE","AI","H.C","O.C"),
                           labels = c("DBE","AI","H/C","O/C"))

cl_a_sel_m$Cl=factor(cl_a_sel_m$Cl, levels = c("All","WSOCbb_pos","WSOCnbb_pos"))

cl_a_sel_m$Clvar=cl_a_sel_m$Cl
cl_a_sel_m$Clvar=factor(cl_a_sel_m$Clvar, levels = c("All","WSOCbb_pos","WSOCnbb_pos"),
                        labels = c(expression(bold("All")),
                                   expression(bold("WSOC"["bb"])),
                                   expression(bold("WSOC"["nbb"]))))


cl_a_max=as.data.table(aggregate(cl_a_sel_m$value, by=list("variable"=cl_a_sel_m$variable, Cl=
                                                             cl_a_sel_m$Cl), FUN=max))

cl_a_max

cl_a_max$x=ifelse(cl_a_max$variable=="DBE",20,cl_a_max$x)

cl_a_sel_m

cl_a_sel_m2=subset(cl_a_sel_m,cl_a_sel_m$Cl!="All") %>% droplevels()

cl_a_max2=as.data.table(aggregate(cl_a_sel_m2$value, by=list("variable"=cl_a_sel_m2$variable, Cl=
                                                               cl_a_sel_m2$Cl), FUN=max))

ggplot(cl_a_sel_m2,aes(x=Group, y=value, fill=Cl))+
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.85), width=0.4, size=0.75)+
  geom_boxplot(position = position_dodge(width = 0.85), outlier.color = NA, size=0.75)+
  geom_blank(data=cl_a_max2, aes(x=2.5,y=x*1.2))+
  scale_fill_manual(values = c("#9DBCD4","#CB7723"),
                    labels=c(expression(bold("WSOC"["bb"])),
                             expression(bold("WSOC"["nbb"]))))+
  facet_wrap(.~variable,scales = "free", strip.position = "left")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.3,0,0.3,0),"cm")),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.3,0.0,0.3),"cm")),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        strip.text.x = element_text(size = 25, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.6,0,0.4,0),"cm")),
        strip.text.y = element_text(size = 25, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.text = element_text(size = 22, colour = "black", family = "Arial",hjust = 0.5,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26,hjust = 0.2 ,vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
        legend.key.height = unit(2,"cm"),
        legend.key.width = unit(2,"cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  guides(fill=guide_legend(title = "", override.aes = list(size=0.75)))+
  ggsave(filename("cl_comp_all"),height = 30, width = 40, units = "cm", dpi = 300)

cl_a_sel_m2

cl_a_sel_m2$Group
cl_a_sel_m2$rel=cl_a_sel_m2$value

grp=unique(cl_a_sel_m2$Group)
comp=unique(cl_a_sel_m2$variable)

uf=data.table()
data=cl_a_sel_m2
for (i in 1:length(comp)) {
  #i=1
  temp=subset(data,data$variable==comp[i])
  for (j in 1:length(grp)) {
    #j=1
    temp2=subset(temp,temp$Group==grp[j])
    temp2
    
    t_eve=subset(temp2,temp2$Cl=="WSOCbb_pos")
    t_nev=subset(temp2,temp2$Cl=="WSOCnbb_pos")
    
    iqr.df <- ddply(temp2, "Cl", function (x) ((1.5*IQR(x[,"rel"])+quantile(x[,"rel"],0.75)))) %>% 
      `colnames<-`(c("Cl","iqr"))
    
    t3=temp2 %>% inner_join(iqr.df)
    
    t4=subset(t3,t3$rel<t3$iqr)
    t4
    
    uf.ps <- ddply(t4, "Cl", function (x) (max(fivenum(x[,"rel"]))))
    new=data.table("Comp"=comp[i],"Group"=grp[j],"uf_eve"=uf.ps[1,2],"uf_nev"=uf.ps[2,2] )
    uf=rbind(uf,new)
    
  }
}




uf$uf_nev=ifelse(uf$Comp=="AI",0,uf$uf_nev)

uf$uf_pos=ifelse(uf$uf_eve>=uf$uf_nev,uf$uf_eve,uf$uf_nev)
uf$Group=factor(uf$Group,levels = c("Seoul","Beijing"))


uf$id=paste(uf$Comp,uf$Group,sep = "_")

cl_a_sel_m2

dt=data.table()
for (i in 1:length(comp)) {
  #i=1
  temp=subset(cl_a_sel_m2,cl_a_sel_m2$variable==comp[i])
  for (j in 1:length(grp)) {
    #j=1
    temp2=subset(temp,temp$Group==grp[j])
    
    t_eve=subset(temp2,temp2$Cl=="WSOCbb_pos")
    t_nev=subset(temp2,temp2$Cl=="WSOCnbb_pos")
    
    if(grp[j]=="Ulaanbaatar"){
      
      new=data.table("Comp"=comp[i],"Group"=grp[j],"norm_eve"=1,"norm_nev"=1,
                     "Eq-var"=1,"T-test"=1)
    }else{
      s1=shapiro.test(t_eve$value)
      s2=shapiro.test(t_nev$value)
      evtmp=levene.test(temp2$value,temp2$Cl)
      #tt=t.test(t_eve$value,t_nev$value, paired = F,var.equal = F)
      #tt=t.test(t_eve$value,t_nev$value, paired = F,var.equal = F)
      tt=wilcox.test(t_eve$value,t_nev$value,exact = F)
      new=data.table("Comp"=comp[i],"Group"=grp[j],"norm_eve"=round(s1$p.value,3),"norm_nev"=round(s2$p.value,3),
                     "Eq-var"=round(evtmp$p,3),"T-test"=round(tt$p.value,3))
      
    }
    dt=rbind(dt,new)
    
  }
}
dt


dt$id=paste(dt$Comp,dt$Group,sep = "_")
#dt$Noneq_t=dt_non_equal$`T-test`
#dt$p=ifelse(dt$`T-test`<dt$Noneq_t,dt$`T-test`,dt$Noneq_t)
dt
uf=uf %>% inner_join(dt[,c("id","T-test")])
uf$sig=ifelse(uf$`T-test`<0.001,"***",
              ifelse(uf$`T-test`<0.01,"**",
                     ifelse(uf$`T-test`<0.05,"*","N.S")))

#uf$Group=factor(uf$Group,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))
uf$Comp=factor(uf$Comp,levels = c("CHO","CHON","CHOS","CHONS"))
uf$variable=uf$Comp
uf$variable=factor(uf$variable, levels = c("DBE","AI","H/C","O/C"))


uf$xs=ifelse(uf$Group=="Seoul",1,
             ifelse(uf$Group=="Seosan",4,
                    ifelse(uf$Group=="Beijing",2,
                           ifelse(uf$Group=="Noto",5,3))))


uf


ggplot()+
  stat_boxplot(data=cl_a_sel_m2,aes(x=Group, y=value, fill=Cl),geom = "errorbar",position = position_dodge(width = 0.85), width=0.4, size=0.75)+
  geom_boxplot(data=cl_a_sel_m2,aes(x=Group, y=value, fill=Cl),position = position_dodge(width = 0.85), outlier.color = NA, size=0.75)+
  geom_segment(data=uf ,aes(x=xs-0.2,xend=xs+0.2,y=uf_pos*1.08, yend=uf_pos*1.08))+
  geom_text(data =uf, aes(x=Group, y=uf_pos*1.13, label=sig), col="black", size=6)+
  geom_blank(data=cl_a_max2, aes(x=2.5,y=x*1.2))+
  scale_fill_manual(values = c("#9DBCD4","#CB7723"),
                    labels=c(expression(bold("WSOC"["bb"])),
                             expression(bold("WSOC"["nbb"]))))+
  facet_wrap(.~variable,scales = "free", strip.position = "left")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.3,0,0.3,0),"cm")),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.3,0.0,0.3),"cm")),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        strip.text.x = element_text(size = 25, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.6,0,0.4,0),"cm")),
        strip.text.y = element_text(size = 25, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.text = element_text(size = 22, colour = "black", family = "Arial",hjust = 0.5,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26,hjust = 0.2 ,vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
        legend.key.height = unit(2,"cm"),
        legend.key.width = unit(2,"cm"),
        #legend.position = "bottom",
        legend.position = c(0.09,0.945),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  guides(fill=guide_legend(title = "", override.aes = list(size=0.75)))+
  ggsave(filename("cl_comp_all_sig"),height = 30, width = 40, units = "cm", dpi = 300)



#uf$Complab=paste(uf$Comp,"(%)")
#uf$Complab=factor(uf$Complab,levels = c("CHO (%)","CHON (%)","CHOS (%)","CHONS (%)"))


#stat_sul=stat_par

stat_b=stat_par


stat_sul

stat_b
