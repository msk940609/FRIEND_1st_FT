library(ggplot2)
library(data.table)
library(dplyr)
library(tidyverse)
options(java.parameters = "- Xmx8192m")
library(extrafont)
loadfonts(device="win")
require(ggplot2)
require(moonBook)
require(webr)
library(data.table)
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


FT_count=as.data.table(aggregate(FT_merge$Freq, by=list(`Group`=FT_merge$Group,`Formula`=FT_merge$Formula), FUN=sum))
FT_count_sel=FT_count %>% filter(Group%in%c("SUL","B","UL"))
FT_count_sel$Group=factor(FT_count_sel$Group,levels = c("SUL","B","UL"),
                          labels = c("Seoul","Beijing","Ulaanbaatar"))
cor_vk=fread("cortest.csv")
cor_vk
fwrite("Datafile/cortest_spearman.csv")
cor_vk=fread("Datafile/cortest_spearman.csv")

cor_vk=cor_vk %>% inner_join(FT_count_sel)
cor_vk

#fwrite(cor_vk,"cor_count.csv")
cor_vk=fread("cor_count.csv")

cor_vk=cor_vk %>% filter(x>2)

cor_vk=subset(cor_vk,cor_vk$p<0.05)
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


cor_vk=cor_vk %>% inner_join(fm)
cor_vk$O.C=cor_vk$O/cor_vk$C
cor_vk$H.C=cor_vk$H/cor_vk$C

cor_vk

cor_vk$C1=ifelse(cor_vk$C>0, "C","")
cor_vk$H1=ifelse(cor_vk$H>0, "H","")
cor_vk$O1=ifelse(cor_vk$O>0, "O","")
cor_vk$N1=ifelse(cor_vk$N>0, "N","")
cor_vk$S1=ifelse(cor_vk$S>0, "S","")

cor_vk=cor_vk %>% unite("Comp",c("C1","H1","O1","N1","S1"),sep = "")
cor_vk$Comp=ifelse(cor_vk$O==0, "Remainders",cor_vk$Comp)

table(cor_vk$Comp)
cor_vk$Comp=factor(cor_vk$Comp,levels = c("CHO","CHON","CHOS","CHONS","Remainders"))

cor_vk$varlab=cor_vk$Envi
cor_vk$varlab=factor(cor_vk$varlab,levels = c("WSOCbb","WSOCnbb"),
                           labels = c(expression(bold("WSOC"["bb"])),
                                      expression(bold("WSOC"["nbb"]))))

cor_vk$Grouplab=cor_vk$Group
cor_vk$Grouplab=factor(cor_vk$Grouplab,levels = c("Ulaanbaatar","Beijing","Seoul"),
                             labels=c(expression(bold("Ulaanbaatar")),
                                      expression(bold("Beijing")),
                                      expression(bold("Seoul"))))

cor_vk



ggplot(cor_vk, aes(x=O.C, y=H.C,col=rho))+
  geom_point(size=2.7)+
  facet_grid(varlab~Grouplab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]))+
  scale_color_gradientn(colors=(topo.colors(40)[4:35]))+
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
        legend.title = element_text(size = 26, vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
        #legend.position = "right",
        #legend.direction = "vertical",
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  guides(col=guide_colorbar(title = expression(bolditalic("ρ")),ticks.colour = "black", ticks.linewidth = 1.5,
                            barwidth = 35, barheight = 2.5))+
  ggsave(filename("vk_cor_raw"),height = 40, width = 60, units = "cm", dpi = 300)

#specialization=====
cor_vk

cor_vk_d=dcast(cor_vk, Group+Comp+Formula~Envi, sum, value.var = c("rho","p"))
cor_vk_d

cor_vk_d$Type=ifelse(abs(cor_vk_d$rho_WSOCbb)>abs(cor_vk_d$rho_WSOCnbb),"BB","NBB")
table(cor_vk_d$Type)
cor_vk_d

cor_vk_d$variable=ifelse(cor_vk_d$Type=="BB","WSOCbb","WSOCnbb")
cor_vk_d$value=ifelse(cor_vk_d$Type=="BB",cor_vk_d$rho_WSOCbb,cor_vk_d$rho_WSOCnbb)
cor_vk_d$p=ifelse(cor_vk_d$Type=="BB",cor_vk_d$p_WSOCbb,cor_vk_d$p_WSOCnbb)
cor_vk_d

cor_vk_d=cor_vk_d %>% inner_join(fm)
cor_vk_d$O.C=cor_vk_d$O/cor_vk_d$C
cor_vk_d$H.C=cor_vk_d$H/cor_vk_d$C

cor_vk_d_sel=cor_vk_d[,c("Group","Comp","Formula","O.C","H.C","Type","variable","value","p")]

cor_vk_d_sel=unique(cor_vk_d_sel)
cor_vk_d_sel
#cor_vk_d_sel=subset(cor_vk_d_sel,abs(cor_vk_d_sel$value)>0.05)

dim(cor_vk_d_sel)
cor_vk_d_sel=subset(cor_vk_d_sel,abs(cor_vk_d_sel$value)>0.3)
table(cor_vk_sel_s1$Type)

cor_vk_d_sel$Grouplab=cor_vk_d_sel$Group
cor_vk_d_sel$Grouplab=factor(cor_vk_d_sel$Grouplab,levels = c("Ulaanbaatar","Beijing","Seoul"),
                              labels=c(expression(bold("Ulaanbaatar")),
                                       expression(bold("Beijing")),
                                       expression(bold("Seoul"))))

cor_vk_d_sel$varlab=cor_vk_d_sel$variable
cor_vk_d_sel$varlab=factor(cor_vk_d_sel$varlab,levels = c("WSOCbb","WSOCnbb"),
                            labels = c(expression(bold("WSOC"["bb"])),
                                       expression(bold("WSOC"["nbb"]))))

##Seoul====
cor_vk_d_sel_s=subset(cor_vk_d_sel,cor_vk_d_sel$Group=="Seoul")
cor_vk_d_sel_s

table(cor_vk_d_sel_s$Type)
table(cor_vk_d_sel_s$variable)

cor_vk_sel_s1=subset(cor_vk_d_sel_s,abs(cor_vk_d_sel_s$value)>0.3)
#cor_vk_sel_s1=subset(cor_vk_d_sel_b,abs(cor_vk_d_sel_b$p)<0.05)
table(cor_vk_sel_s1$Type)

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
  #facet_grid(Grouplab~varlab, labeller = label_parsed )+
  facet_grid(varlab~Grouplab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #annotate(x=0.1,y=0.25, geom = "text", label=32)+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  #scale_color_gradientn(colors=(topo.colors(20)[4:17]))+
  scale_color_gradientn(colors=(topo.colors(40)[4:38]))+
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
                            barwidth = 25, barheight = 2.5))+
  ggsave(filename("vk_cor_s_rel"),height = 40, width = 20, units = "cm", dpi = 300)

##Beijing====
cor_vk_d_sel_b=subset(cor_vk_d_sel,cor_vk_d_sel$Group=="Beijing")
cor_vk_d_sel_b

table(cor_vk_d_sel_b$Type)
table(cor_vk_d_sel_b$variable)

cor_vk_d_sel_b1=subset(cor_vk_d_sel_b,abs(cor_vk_d_sel_b$value)>0.3)
#cor_vk_sel_s1=subset(cor_vk_d_sel_b,abs(cor_vk_d_sel_b$p)<0.05)
table(cor_vk_d_sel_b1$Type)

cor_vk_d_sel_b2=subset(cor_vk_d_sel_b1,abs(cor_vk_d_sel_b1$value)>0.30)


tt=as.data.frame(table(cor_vk_d_sel_b2$Type)) 
tt$varlab=tt$Var1
tt$varlab=factor(tt$varlab,levels = c("BB","NBB"),
                 labels = c(expression(bold("WSOC"["bb"])),
                            expression(bold("WSOC"["nbb"]))))
tt$O.C=0.1
tt$H.C=0.15

ggplot()+
  geom_point(data=cor_vk_d_sel_b2, aes(x=O.C, y=H.C,col=value),size=2.7)+
  #geom_text(data=tt, aes(x=O.C, y=H.C,label=Freq), size=8)+
  geom_text(data=tt, aes(x=0.05, y=0.14,label=paste0(expression(italic("n=")))),parse=T, size=8)+
  geom_text(data=tt, aes(x=0.15, y=H.C,label=Freq),parse=T, size=8)+
  #facet_grid(Grouplab~varlab, labeller = label_parsed )+
  facet_grid(varlab~Grouplab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #annotate(x=0.1,y=0.25, geom = "text", label=32)+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  #scale_color_gradientn(colors=(topo.colors(20)[4:17]))+
  scale_color_gradientn(colors=(topo.colors(40)[4:38]))+
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
                            barwidth = 25, barheight = 2.5))+
  ggsave(filename("vk_cor_b_rel"),height = 40, width = 20, units = "cm", dpi = 300)

##Ulan====
cor_vk_d
cor_vk_d_u=subset(cor_vk_d,cor_vk_d$Group=="Ulaanbaatar")

#cor_vk_d_u$Type=ifelse(abs(cor_vk_d_u$rho_WSOCnbb-cor_vk_d_u$rho_WSOCbb)<0.05,"BB",
#                       ifelse(abs(cor_vk_d_u$rho_WSOCbb)>abs(cor_vk_d_u$rho_WSOCnbb),
#                       "BB","NBB"))

cor_vk_d_u$Type=ifelse(abs(cor_vk_d_u$rho_WSOCbb)>abs(cor_vk_d_u$rho_WSOCnbb),"NBB","BB")

table(cor_vk_d_u$Type)

cor_vk_d

cor_vk_d_u$variable=ifelse(cor_vk_d_u$Type=="BB","WSOCbb","WSOCnbb")
cor_vk_d_u$value=ifelse(cor_vk_d_u$Type=="BB",cor_vk_d_u$rho_WSOCbb,cor_vk_d_u$rho_WSOCnbb)
cor_vk_d_u$p=ifelse(cor_vk_d_u$Type=="BB",cor_vk_d_u$p_WSOCbb,cor_vk_d_u$p_WSOCnbb)
cor_vk_d_u

cor_vk_d_u=cor_vk_d_u %>% inner_join(fm)
cor_vk_d_u$O.C=cor_vk_d_u$O/cor_vk_d_u$C
cor_vk_d_u$H.C=cor_vk_d_u$H/cor_vk_d_u$C

cor_vk_d_u_sel=cor_vk_d_u[,c("Group","Comp","Formula","O.C","H.C","Type","variable","value","p")]

cor_vk_d_u_sel=unique(cor_vk_d_u_sel)
cor_vk_d_u_sel
#cor_vk_d_sel=subset(cor_vk_d_sel,abs(cor_vk_d_sel$value)>0.05)

dim(cor_vk_d_u_sel)
cor_vk_d_u_sel=subset(cor_vk_d_u_sel,abs(cor_vk_d_u_sel$value)>0.01)
cor_vk_d_u_sel=cor_vk_d_u_sel
table(cor_vk_d_u_sel$Type)

cor_vk_d_u_sel$Grouplab=cor_vk_d_u_sel$Group
cor_vk_d_u_sel$Grouplab=factor(cor_vk_d_u_sel$Grouplab,levels = c("Ulaanbaatar","Beijing","Seoul"),
                             labels=c(expression(bold("Ulaanbaatar")),
                                      expression(bold("Beijing")),
                                      expression(bold("Seoul"))))

cor_vk_d_u_sel$varlab=cor_vk_d_u_sel$variable
cor_vk_d_u_sel$varlab=factor(cor_vk_d_u_sel$varlab,levels = c("WSOCbb","WSOCnbb"),
                           labels = c(expression(bold("WSOC"["bb"])),
                                      expression(bold("WSOC"["nbb"]))))

cor_vk_d_u_sel
cor_vk_d_sel_u=subset(cor_vk_d_u_sel,cor_vk_d_u_sel$Group=="Ulaanbaatar")
cor_vk_d_sel_u


table(cor_vk_d_sel_u$Type)
table(cor_vk_d_sel_u$variable)

cor_vk_d_sel_u1=subset(cor_vk_d_sel_u,abs(cor_vk_d_sel_u$value)>0.1)
#cor_vk_sel_s1=subset(cor_vk_d_sel_u,abs(cor_vk_d_sel_u$p)<0.05)
table(cor_vk_d_sel_u1$Type)

cor_vk_d_sel_u2=subset(cor_vk_d_sel_u1,abs(cor_vk_d_sel_u1$value)>0.30)
#cor_vk_d_sel_u2=subset(cor_vk_d_sel_u1,abs(cor_vk_d_sel_u1$p)<0.01)


tt=as.data.frame(table(cor_vk_d_sel_u2$Type)) 
tt$varlab=tt$Var1
tt$varlab=factor(tt$varlab,levels = c("BB","NBB"),
                 labels = c(expression(bold("WSOC"["bb"])),
                            expression(bold("WSOC"["nbb"]))))
tt$O.C=0.1
tt$H.C=0.15


ggplot()+
  geom_point(data=cor_vk_d_sel_u2, aes(x=O.C, y=H.C,col=value),size=2.7)+
  #geom_text(data=tt, aes(x=O.C, y=H.C,label=Freq), size=8)+
  geom_text(data=tt, aes(x=0.05, y=0.14,label=paste0(expression(italic("n=")))),parse=T, size=8)+
  geom_text(data=tt, aes(x=0.15, y=H.C,label=Freq),parse=T, size=8)+
  #facet_grid(Grouplab~varlab, labeller = label_parsed )+
  facet_grid(varlab~Grouplab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #annotate(x=0.1,y=0.25, geom = "text", label=32)+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  #scale_color_gradientn(colors=(topo.colors(20)[4:17]))+
  scale_color_gradientn(colors=(topo.colors(40)[4:38]))+
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
                            barwidth = 25, barheight = 2.5))+
  ggsave(filename("vk_cor_u_rel2"),height = 40, width = 20, units = "cm", dpi = 300)


###chemical composition piechart=====
cor_vk_sel_s2
cor_vk_d_sel_b2
cor_vk_d_sel_u2


t=c("#BC3C29","#EFC000","#008B45","#5F559B","grey50")

cor_vk_sel_s2_pos=subset(cor_vk_sel_s2,cor_vk_sel_s2$value>0)


PieDonut_ms(cor_vk_sel_s2_pos,aes(pies=variable,donuts=Comp),fill=c("#9DBCD4","#CB7723"),
            subfill = rep(t,2))


tiff("S_pie.tiff", width = 20, height = 20,units = "cm",res = 300)
PieDonut_ms(cor_vk_sel_s2_pos,aes(variable,Comp), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill=c("#9DBCD4","#CB7723"),
            subfill =rep(t,2))
dev.off()


cor_vk_d_sel_b2_pos=subset(cor_vk_d_sel_b2,cor_vk_d_sel_b2$value>0)


PieDonut_ms(cor_vk_d_sel_b2,aes(pies=variable,donuts=Comp),fill=c("#9DBCD4","#CB7723"),
            subfill = rep(t,2))

?PieDonut()
tiff("B_pie.tiff", width = 20, height = 20,units = "cm",res = 300)
PieDonut_ms(cor_vk_d_sel_b2_pos,aes(variable,Comp), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill=c("#9DBCD4","#CB7723"),
            subfill =rep(t,2))
dev.off()



cor_vk_d_sel_u2_pos=subset(cor_vk_d_sel_u2,cor_vk_d_sel_u2$value>0)
PieDonut_ms(cor_vk_d_sel_u2,aes(pies=variable,donuts=Comp),fill=c("#9DBCD4","#CB7723"),
            subfill = rep(t,2))

tiff("U_pie.tiff", width = 20, height = 20,units = "cm",res = 300)
PieDonut_ms(cor_vk_d_sel_u2_pos,aes(variable,Comp), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill=c("#9DBCD4","#CB7723"),
            subfill =rep(t,2))
dev.off()

cor_vk_sel_s2_pos
cor_vk_d_sel_b2_pos
cor_vk_d_sel_u2_pos




