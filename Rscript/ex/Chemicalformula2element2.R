library(ggplot2)
library(data.table)
library(dplyr)
library(tidyverse)
options(java.parameters = "- Xmx8192m")
library(extrafont)
loadfonts(device="win")

source("Rscript/func_filename.R")
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")
insert_minor_1 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_2 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}


cor_ev=fread("Datafile/WSOCnbb_cor_event.csv")
dim(cor_ev)

cor_nev=fread("Datafile/WSOCnbb_cor_nonevent.csv")
dim(cor_nev)

cor_temp=cor_nev

grp=unique(cor_temp$Group)

df=data.frame()
for (i in 1:length(grp)) {
  #i=1
  temp=as.data.frame(subset(cor_temp,cor_temp$Group==grp[i]))
  
  for (j in 2:3) {
    #j=4
    for (k in 4:dim(temp)[2]) {
      #k=5  
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
  fwrite(df, file = "cortest_nonevent.csv")
}


#Event=====
cor_vk_ev=fread("cortest_event.csv")
cor_vk_ev

cor_vk_ev=subset(cor_vk_ev,cor_vk_ev$p<0.05)

cor_vk_ev


CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(cor_vk_ev$Formula)

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

cor_vk_ev_m=melt(cor_vk_ev[,c("Group","Formula","O.C","H.C","rho_WSOCbb","rho_WSOCnbb","p_WSOCbb","p_WSOCnbb")],
                 id.vars=c("Group","Formula","O.C","H.C","rho_WSOCbb","rho_WSOCnbb"),variable.name = "p-var",
                 value.name = "p") %>%
  melt(id.vars=c("Group","Formula","O.C","H.C","p-var","p"))
cor_vk_ev_m_sel=subset(cor_vk_ev_m,abs(cor_vk_ev_m$value)>0.001)

cor_vk_ev_m_sel
cor_vk_ev_m_sel$Grouplab=cor_vk_ev_m_sel$Group
cor_vk_ev_m_sel$Grouplab=factor(cor_vk_ev_m_sel$Grouplab,levels = c("Seoul","Beijing"),
                                 labels=c(expression(bold("Seoul")),
                                          expression(bold("Beijing"))))

cor_vk_ev_m_sel$varlab=cor_vk_ev_m_sel$variable
cor_vk_ev_m_sel$varlab=factor(cor_vk_ev_m_sel$varlab,levels = c("rho_WSOCbb","rho_WSOCnbb"),
                               labels = c(expression(bold("WSOC"["bb"])),
                                          expression(bold("WSOC"["nbb"]))))
cor_vk_ev_m_sel=cor_vk_ev_m_sel[order(cor_vk_ev_m_sel$value),]

cor_vk_ev_m_sel_s=subset(cor_vk_ev_m_sel,cor_vk_ev_m_sel$Group=="Seoul")
cor_vk_ev_m_sel_s

cor_vk_ev_m_sel_s1=subset(cor_vk_ev_m_sel_s,cor_vk_ev_m_sel_s$p<0.01)

table(cor_vk_ev_m_sel_s1$variable)

ggplot(cor_vk_ev_m_sel_s1, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_gradientn(colors=(topo.colors(30)[4:24]),breaks = c(-0.9, 0.0, 0.9),limits = c(-1.00,1.00))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
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
        legend.title = element_text(size = 26,hjust = 0.15, vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
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
  ggsave(filename("vk_cor_s_ev"),height = 20, width = 40, units = "cm", dpi = 700)


cor_vk_ev_m_sel_b=subset(cor_vk_ev_m_sel,cor_vk_ev_m_sel$Group=="Beijing")
cor_vk_ev_m_sel_b

cor_vk_ev_m_sel_b1=subset(cor_vk_ev_m_sel_b,cor_vk_ev_m_sel_b$p<0.05)

table(cor_vk_ev_m_sel_b1$variable)

ggplot(cor_vk_ev_m_sel_b1, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_gradientn(colors=(topo.colors(30)[4:24]),breaks = c(-0.9, 0.0, 0.9),limits = c(-1.00,1.00))+
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
        legend.title = element_text(size = 26,hjust = 0.15, vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
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
  ggsave(filename("vk_cor_b_ev"),height = 20, width = 40, units = "cm", dpi = 700)

cor_vk_ev=cor_vk_ev %>% inner_join(fm)

cor_vk_ev$O.C=cor_vk_ev$O/cor_vk_ev$C
cor_vk_ev$H.C=cor_vk_ev$H/cor_vk_ev$C

cor_vk_ev=dcast(cor_vk_ev,Group+Formula+O.C+H.C~Envi,value.var = c("rho","p"),sum) %>%
  `colnames<-`(c("Group","Formula","O.C","H.C","rho_WSOCbb","rho_WSOCnbb","p_WSOCbb","p_WSOCnbb"))

cor_vk_ev

cor_vk_ev$Type=ifelse(abs(cor_vk_ev$rho_WSOCbb)>abs(cor_vk_ev$rho_WSOCnbb),"BB","NBB")
#cor_vk_ev$Type=ifelse(abs(abs(cor_vk_ev$rho_WSOCbb)-abs(cor_vk_ev$rho_WSOCnbb))<0.05,"Both",cor_vk_ev$Type)
table(cor_vk_ev$Type)

cor_vk_ev$variable=ifelse(cor_vk_ev$Type=="BB","WSOCbb","WSOCnbb")
cor_vk_ev$value=ifelse(cor_vk_ev$Type=="BB",cor_vk_ev$rho_WSOCbb,cor_vk_ev$rho_WSOCnbb)
cor_vk_ev$p=ifelse(cor_vk_ev$Type=="BB",cor_vk_ev$p_WSOCbb,cor_vk_ev$p_WSOCnbb)
cor_vk_ev

cor_vk_ev_sel=cor_vk_ev[,c("Group","Formula","O.C","H.C","Type","variable","value","p")]
cor_vk_ev_sel=unique(cor_vk_ev_sel)
cor_vk_ev_sel
cor_vk_ev_sel=subset(cor_vk_ev_sel,abs(cor_vk_ev_sel$value)>0.1)
cor_vk_ev_sel

cor_vk_ev_sel_s=subset(cor_vk_ev_sel,cor_vk_ev_sel$Group=="Seoul")
cor_vk_ev_sel_s

table(cor_vk_ev_sel_s$Type)
table(cor_vk_ev_sel_s$variable)

cor_vk_ev_sel_s1=subset(cor_vk_ev_sel_s,abs(cor_vk_ev_sel_s$value)>0.0)
cor_vk_ev_sel_s1=subset(cor_vk_ev_sel_s,abs(cor_vk_ev_sel_s$p)<0.05)

table(cor_vk_ev_sel_s1$Type)
cor_vk_ev_sel_s1

cor_vk_ev_sel_s1$Grouplab=cor_vk_ev_sel_s1$Group
cor_vk_ev_sel_s1$Grouplab=factor(cor_vk_ev_sel_s1$Grouplab,levels = c("Seoul","Beijing"),
                             labels=c(expression(bold("Seoul")),
                                      expression(bold("Beijing"))))

cor_vk_ev_sel_s1$varlab=cor_vk_ev_sel_s1$variable
cor_vk_ev_sel_s1$varlab=factor(cor_vk_ev_sel_s1$varlab,levels = c("WSOCbb","WSOCnbb"),
                           labels = c(expression(bold("WSOC"["bb"])),
                                      expression(bold("WSOC"["nbb"]))))

ggplot(cor_vk_ev_sel_s1, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  #scale_color_gradientn(colors=(topo.colors(20)[4:17]))+
  scale_color_gradientn(colors=(topo.colors(40)[4:38]),breaks = c(-0.4, 0.0,0.4, 0.8),limits = c(-0.4,0.7))+
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
  ggsave(filename("vk_cor_s_ev_sep"),height = 20, width = 40, units = "cm", dpi = 300)




cor_vk_ev_sel_b=subset(cor_vk_ev_sel,cor_vk_ev_sel$Group=="Beijing")
cor_vk_ev_sel_b

cor_vk_ev_sel_b1=subset(cor_vk_ev_sel_b,abs(cor_vk_ev_sel_b$value)>0.0)
cor_vk_ev_sel_b1=subset(cor_vk_ev_sel_b,abs(cor_vk_ev_sel_b$p)<0.05)

table(cor_vk_ev_sel_b1$Type)
cor_vk_ev_sel_b1

cor_vk_ev_sel_b1$Grouplab=cor_vk_ev_sel_b1$Group
cor_vk_ev_sel_b1$Grouplab=factor(cor_vk_ev_sel_b1$Grouplab,levels = c("Seoul","Beijing"),
                                 labels=c(expression(bold("Seoul")),
                                          expression(bold("Beijing"))))

cor_vk_ev_sel_b1$varlab=cor_vk_ev_sel_b1$variable
cor_vk_ev_sel_b1$varlab=factor(cor_vk_ev_sel_b1$varlab,levels = c("WSOCbb","WSOCnbb"),
                               labels = c(expression(bold("WSOC"["bb"])),
                                          expression(bold("WSOC"["nbb"]))))
cor_vk_ev_sel_b1
ggplot(cor_vk_ev_sel_b1, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  scale_color_gradientn(colors=(topo.colors(30)[4:24]),breaks = c(-0.9, 0.0, 0.9),limits = c(-1.00,1.00))+
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
  ggsave(filename("vk_cor_b_ev_sep"),height = 20, width = 40, units = "cm", dpi = 300)

#Non-event=====
cor_vk_nev=fread("cortest_nonevent.csv")
cor_vk_nev

cor_vk_nev=subset(cor_vk_nev,cor_vk_nev$p<0.05)

cor_vk_nev=dcast(cor_vk_nev,Group+Formula~Envi,value.var = c("rho","p"),sum) %>%
  `colnames<-`(c("Group","Formula","rho_WSOCbb","rho_WSOCnbb","p_WSOCbb","p_WSOCnbb"))

cor_vk_nev


CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(cor_vk_nev$Formula)

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

#cor_vk_nev_m_sel
#cor_vk_nev_m_sel$Grouplab=cor_vk_nev_m_sel$Group
#cor_vk_nev_m_sel$Grouplab=factor(cor_vk_nev_m_sel$Grouplab,lnevels = c("Seoul","Beijing"),
#                                labels=c(expression(bold("Seoul")),
#                                         expression(bold("Beijing"))))

#cor_vk_nev_m_sel$varlab=cor_vk_nev_m_sel$variable
#cor_vk_nev_m_sel$varlab=factor(cor_vk_nev_m_sel$varlab,lnevels = c("rho_WSOCbb","rho_WSOCnbb"),
#                              labels = c(expression(bold("WSOC"["bb"])),
#                                         expression(bold("WSOC"["nbb"]))))

#cor_vk_nev_m_sel_s=subset(cor_vk_nev_m_sel,cor_vk_nev_m_sel$Group=="Seoul")
#cor_vk_nev_m_sel_s

#cor_vk_nev_m_sel_s1=subset(cor_vk_nev_m_sel_s,cor_vk_nev_m_sel_s$p<0.01)



cor_vk_nev=cor_vk_nev %>% inner_join(fm,by="Formula")
cor_vk_nev$O.C=cor_vk_nev$O/cor_vk_nev$C
cor_vk_nev$H.C=cor_vk_nev$H/cor_vk_nev$C

cor_vk_nev$Type=ifelse(abs(cor_vk_nev$rho_WSOCbb)>abs(cor_vk_nev$rho_WSOCnbb),"BB","NBB")
cor_vk_nev$Type=ifelse(abs(abs(cor_vk_nev$rho_WSOCbb)-abs(cor_vk_nev$rho_WSOCnbb))<0.05,"Both",cor_vk_nev$Type)
table(cor_vk_nev$Type)

cor_vk_nev$variable=ifelse(cor_vk_nev$Type=="BB","WSOCbb","WSOCnbb")
cor_vk_nev$value=ifelse(cor_vk_nev$Type=="BB",cor_vk_nev$rho_WSOCbb,cor_vk_nev$rho_WSOCnbb)
cor_vk_nev$p=ifelse(cor_vk_nev$Type=="BB",cor_vk_nev$p_WSOCbb,cor_vk_nev$p_WSOCnbb)
cor_vk_nev

cor_vk_nev_sel=cor_vk_nev[,c("Group","Formula","O.C","H.C","Type","variable","value","p")]
cor_vk_nev_sel=unique(cor_vk_nev_sel)
cor_vk_nev_sel
cor_vk_nev_sel=subset(cor_vk_nev_sel,abs(cor_vk_nev_sel$value)>0.1)
cor_vk_nev_sel

cor_vk_nev_sel_s=subset(cor_vk_nev_sel,cor_vk_nev_sel$Group=="Seoul")
cor_vk_nev_sel_s

table(cor_vk_nev_sel_s$Type)
table(cor_vk_nev_sel_s$variable)

cor_vk_nev_sel_s1=subset(cor_vk_nev_sel_s,abs(cor_vk_nev_sel_s$value)>0.0)
cor_vk_nev_sel_s1=subset(cor_vk_nev_sel_s,abs(cor_vk_nev_sel_s$p)<0.05)

table(cor_vk_nev_sel_s1$Type)
cor_vk_nev_sel_s1

cor_vk_nev_sel_s1$Grouplab=cor_vk_nev_sel_s1$Group
cor_vk_nev_sel_s1$Grouplab=factor(cor_vk_nev_sel_s1$Grouplab,levels = c("Seoul","Beijing"),
                                 labels=c(expression(bold("Seoul")),
                                          expression(bold("Beijing"))))

cor_vk_nev_sel_s1$varlab=cor_vk_nev_sel_s1$variable
cor_vk_nev_sel_s1$varlab=factor(cor_vk_nev_sel_s1$varlab,levels = c("WSOCbb","WSOCnbb"),
                               labels = c(expression(bold("WSOC"["bb"])),
                                          expression(bold("WSOC"["nbb"]))))

ggplot(cor_vk_nev_sel_s1, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  #scale_color_gradientn(colors=(topo.colors(20)[4:17]))+
  scale_color_gradientn(colors=(topo.colors(40)[4:36]))+
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
  ggsave(filename("vk_cor_s_nev_sep"),height = 20, width = 40, units = "cm", dpi = 300)


cor_vk_nev_sel_b=subset(cor_vk_nev_sel,cor_vk_nev_sel$Group=="Beijing")
cor_vk_nev_sel_b

cor_vk_nev_sel_b1=subset(cor_vk_nev_sel_b,abs(cor_vk_nev_sel_b$value)>0.0)
cor_vk_nev_sel_b1=subset(cor_vk_nev_sel_b,abs(cor_vk_nev_sel_b$p)<0.05)

table(cor_vk_nev_sel_b1$Type)
cor_vk_nev_sel_b1

cor_vk_nev_sel_b1$Grouplab=cor_vk_nev_sel_b1$Group
cor_vk_nev_sel_b1$Grouplab=factor(cor_vk_nev_sel_b1$Grouplab,levels = c("Seoul","Beijing"),
                                 labels=c(expression(bold("Seoul")),
                                          expression(bold("Beijing"))))

cor_vk_nev_sel_b1$varlab=cor_vk_nev_sel_b1$variable
cor_vk_nev_sel_b1$varlab=factor(cor_vk_nev_sel_b1$varlab,levels = c("WSOCbb","WSOCnbb"),
                               labels = c(expression(bold("WSOC"["bb"])),
                                          expression(bold("WSOC"["nbb"]))))
cor_vk_nev_sel_b1
ggplot(cor_vk_nev_sel_b1, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  scale_color_gradientn(colors=(topo.colors(30)[4:24]),breaks = c(-0.9, 0.0, 0.9),limits = c(-1.00,1.00))+
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
  ggsave(filename("vk_cor_b_nev_sep"),height = 20, width = 40, units = "cm", dpi = 300)


