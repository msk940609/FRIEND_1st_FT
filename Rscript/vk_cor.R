library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
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

ft_merge=fread("Datafile/FRIENDs_1st_FT.csv")
ft_merge$Group=factor(ft_merge$Group, levels =c("UL","B","SS","SUL","NT"),
                      labels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))

ft_merge$Freq=1

fm_obs=ft_merge[,c("Group","No","Formula","Freq")] %>% 
  dcast(Group~Formula,value.var = "Freq",sum) %>% 
  melt(id.vars=c("Group")) %>% `colnames<-`(c("Group","Formula","Freq"))

ft_cor=fread("Datafile/molecular_cor.csv")
ft_cor=ft_cor %>% left_join(fm_obs)
ft_cor_sel=subset(ft_cor,ft_cor$cnt>4)

ft_cor_sig=subset(ft_cor_sel,ft_cor_sel$p<0.05)

ft_cor_sig

ft_cor_sel

ft_cor_insig=ft_cor_sel %>% filter(!Formula%in%c(ft_cor_sig$Formula))

CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(ft_cor_sel$Formula)
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

numericalFormula

insig_list=as.data.table(numericalFormula) %>% `colnames<-`(c("C","H","N","O","S"))
insig_list$Formula=molecularFormula
insig_list

insig_list$`O/C`=insig_list$O/insig_list$C
insig_list$`H/C`=insig_list$H/insig_list$C

ft_cor_insig=ft_cor_insig %>% inner_join(insig_list)
ft_cor_insig

ft_cor_insig$varlab=factor(ft_cor_insig$Envi,levels = c("WSOCbb","WSOCnbb"),
                              labels = c(expression(bold("WSOC"["BB"])),
                                         expression(bold("WSOC"["NBB"]))))

ft_cor_insig$Grouplab=factor(ft_cor_insig$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                                labels=c(expression(bold("Ulaanbaatar")),
                                         expression(bold("Beijing")),
                                         expression(bold("Seosan")),
                                         expression(bold("Seoul")),
                                         expression(bold("Noto"))))

ft_cor_insig$`O/C`=ft_cor_insig$O/ft_cor_insig$C
ft_cor_insig$`H/C`=ft_cor_insig$H/ft_cor_insig$C
ft_cor_insig


ft_cor_sig2=subset(ft_cor_sig,abs(ft_cor_sig$rho)>0.1)

ft_cor_d_rho=dcast(ft_cor_sig,Group+Formula~Envi,  value.var = c("rho"),fun.aggregate = sum) %>% 
  melt(id.vars=c("Group","Formula")) %>% `colnames<-`(c("Group","Formula","Type","rho"))
  
ft_cor_d_p=dcast(ft_cor_sig,Group+Formula~Envi,  value.var = c("p"),fun.aggregate = sum) %>% 
  melt(id.vars=c("Group","Formula")) %>% `colnames<-`(c("Group","Formula","Type","p"))

ft_cor_d_rho
ft_cor_d_p

ft_cor_d=ft_cor_d_rho %>% inner_join(ft_cor_d_p)

ft_cor_d$Group=factor(ft_cor_d$Group, 
                      levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))


CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(ft_cor_d$Formula)
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

numericalFormula

fmlist=as.data.table(numericalFormula) %>% `colnames<-`(c("C","H","N","O","S"))
fmlist$Formula=molecularFormula
fmlist

ft_cor_d=ft_cor_d %>% inner_join(fmlist)

ft_cor_d$`O/C`=ft_cor_d$O/ft_cor_d$C
ft_cor_d$`H/C`=ft_cor_d$H/ft_cor_d$C


ft_cor_d2=subset(ft_cor_d,abs(ft_cor_d$rho)>0.8)

ggplot(ft_cor_d2, aes(x=`O/C`,y=`H/C`))+
  geom_point(aes(col=rho))+
  facet_wrap(Type~Group, ncol = 5)+
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
  ggsave(filename("vk_cor_all_rho0.8"),height = 42, width = 100, units = "cm", dpi = 300)


ft_cor_d

##Seosan=====
ft_cor_d_ss=subset(ft_cor_d,ft_cor_d$Group=="Seosan")

ft_cor_d_ss

ft_cor_d_ss
ft_cor_d_ss=subset(ft_cor_d_ss,abs(ft_cor_d_ss$rho)>0.2)

ft_cor_d_ss

ft_cor_d_ss_spe=ft_cor_d_ss[,c("Group","Formula","Type","rho")] 

ft_cor_d_ss_spe=melt(ft_cor_d_ss[,c("Group","Formula","Type","rho")] ,id.vars = c("Group","Formula","Type"))
ft_cor_d_ss_spe=dcast(ft_cor_d_ss[,c("Group","Formula","Type","rho")] ,
                       Group+Formula~Type,value.var = "rho",fun.aggregate = sum)
ft_cor_d_ss_spe

ft_cor_d_ss_spe$type=ifelse(abs(ft_cor_d_ss_spe$WSOCbb)>abs(ft_cor_d_ss_spe$WSOCnbb),
                             "WSOCbb","WSOCnbb")

ft_cor_d_ss_spe$value=ifelse(ft_cor_d_ss_spe$type=="WSOCbb",ft_cor_d_ss_spe$WSOCbb,ft_cor_d_ss_spe$WSOCnbb)
ft_cor_d_ss_spe

ft_cor_d_ss_spe=ft_cor_d_ss_spe %>% inner_join(fmlist)
ft_cor_d_ss_spe$`O/C`=ft_cor_d_ss_spe$O/ft_cor_d_ss_spe$C
ft_cor_d_ss_spe$`H/C`=ft_cor_d_ss_spe$H/ft_cor_d_ss_spe$C

ft_cor_d_ss_spe

table(ft_cor_d_ss_spe$type)
ft_cor_d_ss_spe_pos=subset(ft_cor_d_ss_spe,ft_cor_d_ss_spe$value>0)
table(ft_cor_d_ss_spe_pos$type)

ft_cor_d_ss_spe$varlab=factor(ft_cor_d_ss_spe$type,levels = c("WSOCbb","WSOCnbb"),
                            labels = c(expression(bold("WSOC"["BB"])),
                                       expression(bold("WSOC"["NBB"]))))

ft_cor_d_ss_spe$Grouplab=factor(ft_cor_d_ss_spe$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                          labels=c(expression(bold("Ulaanbaatar")),
                                   expression(bold("Beijing")),
                                   expression(bold("Seosan")),
                                   expression(bold("Seoul")),
                                   expression(bold("Noto"))))

ft_cor_d_ss_spe

ft_cor_d_ss_spe2=subset(ft_cor_d_ss_spe,abs(ft_cor_d_ss_spe$value)>0.4)

table(ft_cor_d_ss_spe2$type)

ft_cor_sel_ss=subset(ft_cor_sel,ft_cor_sel$Group=="Seosan")

ft_cor_insig_ss=ft_cor_sel_ss %>% filter(!Formula%in%c(ft_cor_d_ss$Formula))
ft_cor_insig_ss=ft_cor_insig_ss %>% left_join(insig_list)

ggplot()+
  geom_point(data=ft_cor_insig_ss, aes(x=`O/C`, y=`H/C`),col="grey70",size=1.5)+
  geom_point(data=ft_cor_d_ss_spe2, aes(x=`O/C`, y=`H/C`,col=value),size=2.7)+
  facet_grid(varlab ~Grouplab, labeller = label_parsed )+
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
        strip.text.y = element_text(size = 40, colour = "white",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
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
  ggsave(filename("vk_cor_ss_rho0.4"),height = 40, width = 20, units = "cm", dpi = 300)


##sul====
ft_cor_d_sul=subset(ft_cor_d,ft_cor_d$Group=="Seoul")

ft_cor_d_sul=subset(ft_cor_d_sul,abs(ft_cor_d_sul$rho)>0.2)

ft_cor_d_sul

ft_cor_d_sul_spe=ft_cor_d_sul[,c("Group","Formula","Type","rho")] 

ft_cor_d_sul_spe=melt(ft_cor_d_sul[,c("Group","Formula","Type","rho")] ,id.vars = c("Group","Formula","Type"))
ft_cor_d_sul_spe=dcast(ft_cor_d_sul[,c("Group","Formula","Type","rho")] ,
                       Group+Formula~Type,value.var = "rho",fun.aggregate = sum)
ft_cor_d_sul_spe

ft_cor_d_sul_spe$type=ifelse(abs(ft_cor_d_sul_spe$WSOCbb)>abs(ft_cor_d_sul_spe$WSOCnbb),
                             "WSOCbb","WSOCnbb")

ft_cor_d_sul_spe$value=ifelse(ft_cor_d_sul_spe$type=="WSOCbb",ft_cor_d_sul_spe$WSOCbb,ft_cor_d_sul_spe$WSOCnbb)
ft_cor_d_sul_spe

ft_cor_d_sul_spe=ft_cor_d_sul_spe %>% inner_join(fmlist)
ft_cor_d_sul_spe$`O/C`=ft_cor_d_sul_spe$O/ft_cor_d_sul_spe$C
ft_cor_d_sul_spe$`H/C`=ft_cor_d_sul_spe$H/ft_cor_d_sul_spe$C

ft_cor_d_sul_spe

ft_cor_d_sul_spe$varlab=factor(ft_cor_d_sul_spe$type,levels = c("WSOCbb","WSOCnbb"),
                              labels = c(expression(bold("WSOC"["BB"])),
                                         expression(bold("WSOC"["NBB"]))))

ft_cor_d_sul_spe$Grouplab=factor(ft_cor_d_sul_spe$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                                labels=c(expression(bold("Ulaanbaatar")),
                                         expression(bold("Beijing")),
                                         expression(bold("Seosan")),
                                         expression(bold("Seoul")),
                                         expression(bold("Noto"))))

table(ft_cor_d_sul_spe$type)
ft_cor_d_sul_spe_pos=subset(ft_cor_d_sul_spe,ft_cor_d_sul_spe$value>0)
table(ft_cor_d_sul_spe_pos$type)

ft_cor_d_sul_spe2=subset(ft_cor_d_sul_spe,abs(ft_cor_d_sul_spe$value)>0.2)

ft_cor_sel_sul=subset(ft_cor_sel,ft_cor_sel$Group=="Seoul")
ft_cor_insig_sul=ft_cor_sel_sul %>% filter(!Formula%in%c(ft_cor_d_sul$Formula))
ft_cor_insig_sul=ft_cor_insig_sul %>% left_join(insig_list)

ggplot()+
  geom_point(data=ft_cor_insig_sul, aes(x=`O/C`, y=`H/C`),col="grey70",size=1.5)+
  geom_point(data=ft_cor_d_sul_spe2, aes(x=`O/C`, y=`H/C`,col=value),size=2.7)+
  facet_grid(varlab ~Grouplab, labeller = label_parsed )+
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
        strip.text.y = element_text(size = 40, colour = "white",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
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
  ggsave(filename("vk_cor_sul_spe_rho0.4"),height = 40, width = 20, units = "cm", dpi = 300)


##beijing====
ft_cor_d_bj=subset(ft_cor_d,ft_cor_d$Group=="Beijing")

ft_cor_d_bj=subset(ft_cor_d_bj,abs(ft_cor_d_bj$rho)>0.2)

ft_cor_d_bj

ft_cor_d_bj_spe=ft_cor_d_bj[,c("Group","Formula","Type","rho")] 

ft_cor_d_bj_spe=melt(ft_cor_d_bj[,c("Group","Formula","Type","rho")] ,id.vars = c("Group","Formula","Type"))
ft_cor_d_bj_spe=dcast(ft_cor_d_bj[,c("Group","Formula","Type","rho")] ,
                       Group+Formula~Type,value.var = "rho",fun.aggregate = sum)
ft_cor_d_bj_spe

ft_cor_d_bj_spe$type=ifelse(abs(ft_cor_d_bj_spe$WSOCbb)>abs(ft_cor_d_bj_spe$WSOCnbb),
                             "WSOCbb","WSOCnbb")

ft_cor_d_bj_spe$value=ifelse(ft_cor_d_bj_spe$type=="WSOCbb",ft_cor_d_bj_spe$WSOCbb,ft_cor_d_bj_spe$WSOCnbb)
ft_cor_d_bj_spe

ft_cor_d_bj_spe=ft_cor_d_bj_spe %>% inner_join(fmlist)
ft_cor_d_bj_spe$`O/C`=ft_cor_d_bj_spe$O/ft_cor_d_bj_spe$C
ft_cor_d_bj_spe$`H/C`=ft_cor_d_bj_spe$H/ft_cor_d_bj_spe$C

ft_cor_d_bj_spe

ft_cor_d_bj_spe$varlab=factor(ft_cor_d_bj_spe$type,levels = c("WSOCbb","WSOCnbb"),
                               labels = c(expression(bold("WSOC"["BB"])),
                                          expression(bold("WSOC"["NBB"]))))

ft_cor_d_bj_spe$Grouplab=factor(ft_cor_d_bj_spe$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                                 labels=c(expression(bold("Ulaanbaatar")),
                                          expression(bold("Beijing")),
                                          expression(bold("Seosan")),
                                          expression(bold("Seoul")),
                                          expression(bold("Noto"))))
table(ft_cor_d_bj_spe$type)
ft_cor_d_bj_spe_pos=subset(ft_cor_d_bj_spe,ft_cor_d_bj_spe$value>0)
table(ft_cor_d_bj_spe_pos$type)

ft_cor_d_bj_spe2=subset(ft_cor_d_bj_spe,abs(ft_cor_d_bj_spe$value)>0.4)

table(ft_cor_d_bj_spe2$type)


ft_cor_sel_bj=subset(ft_cor_sel,ft_cor_sel$Group=="Beijing")
ft_cor_insig_bj=ft_cor_sel_bj %>% filter(!Formula%in%c(ft_cor_d_bj$Formula))
ft_cor_insig_bj=ft_cor_insig_bj %>% left_join(insig_list)

ggplot()+
  geom_point(data=ft_cor_insig_bj, aes(x=`O/C`, y=`H/C`),col="grey70",size=1.5)+
  geom_point(data=ft_cor_d_bj_spe2, aes(x=`O/C`, y=`H/C`,col=value),size=2.7)+
  facet_grid(varlab ~Grouplab, labeller = label_parsed )+
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
        strip.text.y = element_text(size = 40, colour = "white",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
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
  ggsave(filename("vk_cor_bj_spe_rho0.4"),height = 40, width = 20, units = "cm", dpi = 300)

##Noto====
ft_cor_d_nt=subset(ft_cor_d,ft_cor_d$Group=="Noto")

ft_cor_d_nt=subset(ft_cor_d_nt,abs(ft_cor_d_nt$rho)>0.2)

ft_cor_d_nt

ft_cor_d_nt_spe=ft_cor_d_nt[,c("Group","Formula","Type","rho")] 

ft_cor_d_nt_spe=melt(ft_cor_d_nt[,c("Group","Formula","Type","rho")] ,id.vars = c("Group","Formula","Type"))
ft_cor_d_nt_spe=dcast(ft_cor_d_nt[,c("Group","Formula","Type","rho")] ,
                      Group+Formula~Type,value.var = "rho",fun.aggregate = sum)
ft_cor_d_nt_spe

ft_cor_d_nt_spe$type=ifelse(abs(ft_cor_d_nt_spe$WSOCbb)>abs(ft_cor_d_nt_spe$WSOCnbb),
                            "WSOCbb","WSOCnbb")

ft_ft_cor_d_nt_spe
ft_cor_d_nt_spe$varlab=factor(ft_cor_d_nt_spe$type,levels = c("WSOCbb","WSOCnbb"),
                              labels = c(expression(bold("WSOC"["BB"])),
                                         expression(bold("WSOC"["NBB"]))))

ft_cor_d_nt_spe$Grouplab=factor(ft_cor_d_nt_spe$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                                labels=c(expression(bold("Ulaanbaatar")),
                                         expression(bold("Beijing")),
                                         expression(bold("Seosan")),
                                         expression(bold("Seoul")),
                                         expression(bold("Noto"))))

cor_d_nt_spe$value=ifelse(ft_cor_d_nt_spe$type=="WSOCbb",ft_cor_d_nt_spe$WSOCbb,ft_cor_d_nt_spe$WSOCnbb)
ft_cor_d_nt_spe

ft_cor_d_nt_spe=ft_cor_d_nt_spe %>% inner_join(fmlist)
ft_cor_d_nt_spe$`O/C`=ft_cor_d_nt_spe$O/ft_cor_d_nt_spe$C
ft_cor_d_nt_spe$`H/C`=ft_cor_d_nt_spe$H/ft_cor_d_nt_spe$C


table(ft_cor_d_nt_spe$type)
ft_cor_d_nt_spe_pos=subset(ft_cor_d_nt_spe,ft_cor_d_nt_spe$value>0)
table(ft_cor_d_nt_spe_pos$type)

ft_cor_d_nt_spe2=subset(ft_cor_d_nt_spe,abs(ft_cor_d_nt_spe$value)>0.2)

table(ft_cor_d_nt_spe2$type)

ft_cor_sel_nt=subset(ft_cor_sel,ft_cor_sel$Group=="Noto")
ft_cor_insig_nt=ft_cor_sel_nt %>% filter(!Formula%in%c(ft_cor_d_nt$Formula))
ft_cor_insig_nt=ft_cor_insig_nt %>% left_join(insig_list)

ggplot()+
  geom_point(data=ft_cor_insig_nt, aes(x=`O/C`, y=`H/C`),col="grey70",size=1.5)+
  geom_point(data=ft_cor_d_nt_spe2, aes(x=`O/C`, y=`H/C`,col=value),size=2.7)+
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
  ggsave(filename("vk_cor_nt_spe_rho0.4"),height = 40, width = 20, units = "cm", dpi = 300)

##ulaanbaatar====
ft_cor_d_ul=subset(ft_cor_d,ft_cor_d$Group=="Ulaanbaatar")
ft_cor_d_ul=subset(ft_cor_d_ul,abs(ft_cor_d_ul$rho)>0.2)

ft_cor_d_ul
ft_cor_d_ul_spe=ft_cor_d_ul[,c("Group","Formula","Type","rho")] 

ft_cor_d_ul_spe=melt(ft_cor_d_ul[,c("Group","Formula","Type","rho")] ,id.vars = c("Group","Formula","Type"))
ft_cor_d_ul_spe=dcast(ft_cor_d_ul[,c("Group","Formula","Type","rho")] ,
                      Group+Formula~Type,value.var = "rho",fun.aggregate = sum)
ft_cor_d_ul_spe

ft_cor_d_ul_spe$type=ifelse(abs(ft_cor_d_ul_spe$WSOCbb)<abs(ft_cor_d_ul_spe$WSOCnbb),
                            "WSOCbb","WSOCnbb")

ft_cor_d_ul_spe$value=ifelse(ft_cor_d_ul_spe$type=="WSOCbb",ft_cor_d_ul_spe$WSOCbb,ft_cor_d_ul_spe$WSOCnbb)
ft_cor_d_ul_spe

ft_cor_d_ul_spe=ft_cor_d_ul_spe %>% inner_join(fmlist)
ft_cor_d_ul_spe$`O/C`=ft_cor_d_ul_spe$O/ft_cor_d_ul_spe$C
ft_cor_d_ul_spe$`H/C`=ft_cor_d_ul_spe$H/ft_cor_d_ul_spe$C

ft_cor_d_ul_spe

ft_cor_d_ul_spe$varlab=factor(ft_cor_d_ul_spe$type,levels = c("WSOCbb","WSOCnbb"),
                              labels = c(expression(bold("WSOC"["BB"])),
                                         expression(bold("WSOC"["NBB"]))))

ft_cor_d_ul_spe$Grouplab=factor(ft_cor_d_ul_spe$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                                labels=c(expression(bold("Ulaanbaatar")),
                                         expression(bold("Beijing")),
                                         expression(bold("Seosan")),
                                         expression(bold("Seoul")),
                                         expression(bold("Noto"))))

table(ft_cor_d_ul_spe$type)
ft_cor_d_ul_spe_pos=subset(ft_cor_d_ul_spe,ft_cor_d_ul_spe$value>0)
table(ft_cor_d_ul_spe_pos$type)

ft_cor_d_ul_spe2=subset(ft_cor_d_ul_spe,abs(ft_cor_d_ul_spe$value)>0.2)

table(ft_cor_d_ul_spe2$type)

ft_cor_sel_ul=subset(ft_cor_sel,ft_cor_sel$Group=="Ulaanbaatar")
ft_cor_insig_ul=ft_cor_sel_ul %>% filter(!Formula%in%c(ft_cor_d_ul$Formula))
ft_cor_insig_ul=ft_cor_insig_ul %>% left_join(insig_list)

ggplot()+
  geom_point(data=ft_cor_insig_ul, aes(x=`O/C`, y=`H/C`),col="grey70",size=1.5)+
  geom_point(data=ft_cor_d_ul_spe2, aes(x=`O/C`, y=`H/C`,col=value),size=2.7)+
  facet_grid(varlab~Grouplab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #annotate(x=0.1,y=0.25, geom = "text", label=32)+
  #scale_color_gradieuln(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  #scale_color_gradieuln(colors=(topo.colors(20)[4:17]))+
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
        strip.text.y = element_text(size = 40, colour = "white",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
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
  ggsave(filename("vk_cor_ul_spe_rho0.4"),height = 40, width = 20, units = "cm", dpi = 300)




ft_cor_d_ul_spe2
ft_cor_d_bj_spe2
ft_cor_d_sul_spe2
ft_cor_d_ss_spe2
ft_cor_d_nt_spe2

ft_cor=rbind(ft_cor_d_ul_spe2,
             ft_cor_d_bj_spe2,
             ft_cor_d_sul_spe2,
             ft_cor_d_ss_spe2,
             ft_cor_d_nt_spe2)

ft_cor

ft_cor_pos=subset(ft_cor,ft_cor$value>0)

ft_cor_pos$DBE=ft_cor_pos$C-ft_cor_pos$H/2+ft_cor_pos$N/2 +1
ft_cor_pos$CAI=ft_cor_pos$C-ft_cor_pos$N-ft_cor_pos$O-ft_cor_pos$S
ft_cor_pos$DBEAI=1+ft_cor_pos$`C`-ft_cor_pos$`O`-ft_cor_pos$`S`-ft_cor_pos$`H`/2-ft_cor_pos$`N`/2
ft_cor_pos$AI=ifelse(ft_cor_pos$CAI<=0,0,ifelse(ft_cor_pos$DBEAI<0,0,ft_cor_pos$DBEAI/ft_cor_pos$CAI))
ft_cor_pos$`N/C`=ft_cor_pos$N/ft_cor_pos$C
ft_cor_pos$`S/C`=ft_cor_pos$S/ft_cor_pos$C
ft_cor_pos$`N/O`=ft_cor_pos$N/ft_cor_pos$O
ft_cor_pos$`S/O`=ft_cor_pos$S/ft_cor_pos$O

ft_cor_pos$NOSC=4-(1+4*ft_cor_pos$C+ft_cor_pos$H-3*ft_cor_pos$N-2*ft_cor_pos$O-2*ft_cor_pos$S)/ft_cor_pos$C
ft_cor_pos$GCox=60.3-28.5*ft_cor_pos$NOSC

ft_cor_pos$Freq=1


ft_cor_pos$AIclass=ifelse(ft_cor_pos$AI==0,"Aliphatic",
                          ifelse(ft_cor_pos$AI<=0.5,"Olefinic","Aromatic"))


#fwrite(ft_cor_pos,file = "Datafile/cor_pos_chp.csv")


ft_cor_pos_m=melt(ft_cor_pos[,-c("WSOCbb","WSOCnbb","value","C","H","N","O","S","varlab","Grouplab","CAI","DBEAI"
                                 ,"N/C","S/C","N/O","S/O","Freq")],
                  id.vars = c("Group","Formula","type","AIclass"))
ft_cor_pos_m
ft_cor_pos_m$variable2=factor(ft_cor_pos_m$variable,levels = c("AI","DBE","H/C","O/C","NOSC","GCox"),
                             labels=c(expression(bold("AI")),
                                      expression(bold("DBE")),
                                      expression(bold("H/C")),
                                      expression(bold("O/C")),
                                      expression(bold("NOSC")),
                                      expression(bold(Delta*"G"[Cox]))))

ft_cor_pos_m


ft_cor_pos_m2=ft_cor_pos_m %>% filter(variable%in%c("AI","DBE","O/C","H/C"))

ggplot()+
  stat_boxplot(data=ft_cor_pos_m2, aes(x=Group, y=value, fill=type),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=ft_cor_pos_m2,aes(x=Group,y=value,fill=type), outlier.colour = NA)+
  geom_point(data=ft_cor_pos_m2,aes(x=Group,y=value*1.1),col=NA)+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  scale_fill_manual(values = c("#9DBCD4","#CB7723"),
                    labels=c(expression(bold("WSOC"["BB"])),
                             expression(bold("WSOC"["NBB"]))))+  
  facet_rep_wrap(.~variable2, nrow=2, repeat.tick.labels = T, strip.position="left", scales="free", dir="h",labeller = label_parsed)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.0),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 12,colour = "black",face=2, family = "Arial",angle = 0, hjust = 0.5,vjust = 1.0),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.12,"cm"),
        axis.ticks = element_line(size = 1),
        axis.text.y = element_text(size = 12, colour = "black",face=2, family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        axis.title.y.left = element_text(size = 0, face=2,colour = "black",margin = unit(c(0.0,0.0,0.0,0.0),"cm")),
        axis.title.y.right = element_text(size = 0, colour = "black",angle = 90,margin = unit(c(0.0,0.0,0.0,0.0),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        #strip.text.y = element_blank(),
        strip.placement = "outside",
        legend.title = element_text(size = 22, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 14,face=2, hjust=0,colour = "black",margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        legend.key.width = unit(0.8,"cm"),
        legend.key.height = unit(0.8,"cm"),
        legend.background = element_blank(),
        legend.direction = "vertical",
        legend.justification=c(0.5, 0.5),
        legend.position =c(0.92,0.96))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("bb&nbb_chp.tiff"),height = 22, width = 25, units = "cm", dpi = 300)


ft_cor_pos_m2

#wsocbb vs wsoc nbb=====
ft_cor_pos_m2

grp=unique(ft_cor_pos_m2$Group)
var=unique(ft_cor_pos_m2$variable)
data_in=ft_cor_pos_m2

dt=data.table()
for (i in 1:length(grp)) {
  #i=1
  
  temp=subset(data_in,data_in$Group==grp[i])
  
  for (j in 1:length(var)) {
    
    #j=1
    temp2=subset(temp,temp$variable==var[j])
    
    t1=subset(temp2,temp2$type=="WSOCbb")
    t2=subset(temp2,temp2$type=="WSOCnbb")
    
    tt=wilcox.test(t1$value,t2$value, exact = F)
    
    new=data.table(Group=grp[i],variable=var[j],stat=round(tt$statistic,3),p=round(tt$p.value,3))
    
    dt=rbind(dt,new)
  }
  
}

dt

ft_cor_pos_m2

bb=subset(ft_cor_pos_m2,ft_cor_pos_m2$type=="WSOCbb")
nbb=subset(ft_cor_pos_m2,ft_cor_pos_m2$type=="WSOCnbb")

kruskal(bb$value,bb$Group,console = T, p.adj = "fdr")
kruskal(nbb$value,nbb$Group,console = T, p.adj = "fdr")

ggplot()+
  #stat_boxplot(data=ft_cor_pos_m, aes(x=Group, y=value, fill=type),geom='errorbar', linetype=1, width=0.25,
  #             position = position_dodge(width = 0.75))+
  geom_violin(data=ft_cor_pos_m,aes(x=Group,y=value,fill=type), outlier.colour = NA)+
  geom_point(data=ft_cor_pos_m,aes(x=Group,y=value*1.1),col=NA)+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  scale_fill_manual(values = c("#9DBCD4","#CB7723"),
                    labels=c(expression(bold("WSOC"["BB"])),
                             expression(bold("WSOC"["NBB"]))))+  
  facet_rep_wrap(.~variable2, nrow=2, repeat.tick.labels = T, strip.position="left", scales="free", dir="h",labeller = label_parsed)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.0),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 17,colour = "black", face = 2,angle = 0, hjust = 0.5,vjust = 1.0),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15,"cm"),
        axis.text.y = element_text(size = 17, colour = "black",face = 2, margin = unit(c(0.0,0.0,0.2,0.2),"cm")),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        #axis.title.y.left = element_text(size = 0, colour = "black",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        #axis.title.y.right = element_text(size = 0, colour = "black",angle = 90,margin = unit(c(0.2,0.2,0.2,0.5),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.placement = "outside",
        strip.text.y = element_text(size = 22, colour = "black",face = 2,margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        strip.text.x = element_blank(),
        #legend.margin = unit(c(0.0,0.0,0.0,0.0),"cm"),
        legend.title = element_text(size = 22, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 26, colour = "black",face=2, hjust=0,margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        legend.key.width = unit(1.5,"cm"),
        legend.key.height = unit(1.5,"cm"),
        legend.position = c(0.58,0.96),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(fill=guide_legend(title=""))+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("bb&nbb_chp_violin.tiff"),height = 22, width = 25, units = "cm", dpi = 300)

ai=as.data.table(aggregate(ft_cor_pos$Freq, by=list(Group=ft_cor_pos$Group,
                                                       type=ft_cor_pos$type,
                                                       class=ft_cor_pos$AIclass), sum))
ai

ai$class=factor(ai$class,levels = c("Aromatic","Olefinic","Aliphatic"))
ai$varlab=factor(ai$type,levels = c("WSOCbb","WSOCnbb"),
                              labels = c(expression(bold("WSOC"["BB"])),
                                         expression(bold("WSOC"["NBB"]))))


ggplot(ai)+
  geom_bar(aes(x=Group,y=x,fill=type),stat = "identity", position = "dodge")+
  facet_wrap(.~class,scales = "free")+
  scale_y_continuous(name = "",
                     expand = c(0.01,0.01))+
  scale_fill_manual(values = c("#9DBCD4","#CB7723"),
                    labels=c(expression(bold("WSOC"["BB"])),
                             expression(bold("WSOC"["NBB"]))))+
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
        strip.text.y = element_text(size = 40, colour = "white",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
        strip.background = element_blank(),
        legend.key.height = unit(1.5,"cm"),
        legend.key.width = unit(1.5,"cm"),
        legend.text = element_text(size = 26, colour = "black",face=2, family = "Arial",hjust = 0.5,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26,hjust = 0.2 ,vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
        #legend.position = "right",
        #legend.direction = "vertical",
        legend.position = c(0.62,0.94),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  guides(fill=guide_legend(title=""))+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("bb&nbb_aiclass.tiff"),height = 40, width = 100, units = "cm", dpi = 300)


ggplot(ai)+
  geom_bar(aes(x=type,y=x,fill=class),stat = "identity", position = position_fill(reverse = T))+
  facet_wrap(.~Group,scales = "free", ncol=5)+
  scale_x_discrete(name="", expand = c(0.25,0.25),labels=c(expression(bold("WSOC"["BB"])),expression(bold("WSOC"["NBB"]))))+
  scale_y_continuous(name = "",labels = scales::percent,
                     expand = c(0.01,0.01))+
# scale_fill_manual(values = c("#9DBCD4","#CB7723"),
#                    labels=c(expression(bold("WSOC"["BB"])),
#                             expression(bold("WSOC"["NBB"]))))+
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
        strip.text.y = element_text(size = 40, colour = "white",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
        strip.background = element_blank(),
        legend.key.height = unit(1.5,"cm"),
        legend.key.width = unit(1.5,"cm"),
        legend.text = element_text(size = 26, colour = "black",face=2, family = "Arial",hjust = 0.5,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26,hjust = 0.2 ,vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
        #legend.position = "right",
        #legend.direction = "vertical",
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  guides(fill=guide_legend(title=""))+
  #guides(col=guide_legend(title = "",ncol = 2,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF"),alpha=1)), size="none")+
  ggsave(filename("bb&nbb_aiclass2.tiff"),height = 30, width = 100, units = "cm", dpi = 300)

ft_cor_pos

ft_cor_pos=ft_cor_pos %>% mutate(Comp = gsub('[[:digit:]]', '', Formula))
ft_cor_pos$Comp=ifelse(ft_cor_pos$`O`==0,"Remainders",ft_cor_pos$Comp)
ft_cor_pos$Comp=factor(ft_cor_pos$Comp,levels = c("CHO","CHNO","CHOS","CHNOS","Remainders"),
                     labels = c("CHO","CHON","CHOS","CHONS","Remainders"))
ft_cor_pos

soa_comp=as.data.table(aggregate(ft_cor_pos$Freq, by=list(Group=ft_cor_pos$Group,
                                                          type=ft_cor_pos$type,
                                                          Comp=ft_cor_pos$Comp),sum)) %>% 
  `colnames<-`(c("Group","Type","Comp","cnt"))

soa_comp_tot=as.data.table(aggregate(ft_cor_pos$Freq, by=list(Group=ft_cor_pos$Group,
                                                          type=ft_cor_pos$type),sum)) %>% 
  `colnames<-`(c("Group","Type","tot"))

soa_comp_tot

soa_comp=soa_comp %>% inner_join(soa_comp_tot)
soa_comp$rel=round(soa_comp$cnt/soa_comp$tot*100,1)

soa_comp_sel=subset(soa_comp,soa_comp$Comp!="Remainders") %>% droplevels()
soa_comp_sel


ggplot()+
  geom_boxplot(data=soa_comp_sel, aes(x=Group, y=rel, fill=Type),stat = "identity",position = "dodge")+
  facet_rep_wrap(.~Comp)

soa_comp_sel$compvar=paste(soa_comp_sel$Comp, " (%)",sep = "")
soa_comp_sel$compvar=factor(soa_comp_sel$compvar,levels = c("CHO (%)","CHON (%)","CHOS (%)","CHONS (%)"))


ggplot(soa_comp_sel, aes(x=Group,y=rel,fill=Type))+
  geom_bar(stat="identity",position = position_dodge(preserve = "single"))+
  geom_text(data=soa_comp_sel, aes(x=Group, y=rel*1.1),label="",size=0, color=NA)+
  facet_rep_wrap(.~compvar, nrow=2, repeat.tick.labels = T, strip.position="left", scales="free")+
  scale_y_continuous(name = "",
                     expand = c(0.03,0.03))+
  scale_fill_manual(values = c("#9DBCD4","#CB7723"),
                    labels=c(expression(bold("WSOC"["BB"])),
                             expression(bold("WSOC"["NBB"]))))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.0),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 12,colour = "black",face=2, family = "Arial",angle = 0, hjust = 0.5,vjust = 1.0),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.12,"cm"),
        axis.ticks = element_line(size = 1),
        axis.text.y = element_text(size = 12, colour = "black",face=2, family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        axis.title.y.left = element_text(size = 0, face=2,colour = "black",margin = unit(c(0.0,0.0,0.0,0.0),"cm")),
        axis.title.y.right = element_text(size = 0, colour = "black",angle = 90,margin = unit(c(0.0,0.0,0.0,0.0),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        #strip.text.y = element_blank(),
        strip.placement = "outside",
        legend.title = element_text(size = 22, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 14,face=2, hjust=0,colour = "black",margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        legend.key.width = unit(0.8,"cm"),
        legend.key.height = unit(0.8,"cm"),
        legend.background = element_blank(),
        legend.direction = "vertical",
        legend.justification=c(0.5, 0.5),
        legend.position =c(0.92,0.43))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("comp_comp_bb&nbb"),height = 22, width = 25, units = "cm", dpi = 300)

##chemical class distribution=====
ft_cor_pos
ft_data=ft_cor_pos

ft_data$PCA=ifelse(ft_data$AI>0.66, "Combustion-drived PCA","")

ft_data$aroma=ifelse(ft_data$AI>0.5,ifelse(ft_data$AI<=0.66,"Soil-derived polyphenols and PCAs with aliphatic chains",""),"")

ft_data$humic=ifelse(ft_data$AI<=0.5,
                     ifelse(ft_data$`H/C`<1.5,"Soil derived humics and highly unsaturated compounds",""),"")

ft_data$aliphatic=ifelse(ft_data$`H/C`<2.0&ft_data$`H/C`>=1.5,
                         ifelse(ft_data$AI<=0.5,"Unsaturated aliphatic compounds",""),"")

ft_data$fatty=ifelse(ft_data$`H/C`>=2.0,
                     ifelse(ft_data$AI<=0.5,"Saturated fatty and sulfonic acid,carbohyrates",""),"")

ft_data=ft_data %>%  unite("Molecularclass",c("PCA","aroma","humic","aliphatic","fatty"), sep = "")
ft_data$Molecularclass=ifelse(ft_data$Molecularclass=="","Unassigned", ft_data$Molecularclass) ###colum of molecular class
ft_data

ft_data$Molecularclass=factor(ft_data$Molecularclass,levels = rev(c("Saturated fatty and sulfonic acid,carbohyrates","Unsaturated aliphatic compounds",
                                                                    "Soil derived humics and highly unsaturated compounds",
                                                                    "Soil-derived polyphenols and PCAs with aliphatic chains",
                                                                    "Combustion-drived PCA")))

ft_data
ft_data$Freq=1


ft_data

ft_data$Class2=factor(ft_data$Molecularclass,levels = rev(c("Saturated fatty and sulfonic acid,carbohyrates","Unsaturated aliphatic compounds",
                                                            "Soil derived humics and highly unsaturated compounds",
                                                            "Soil-derived polyphenols and PCAs with aliphatic chains",
                                                            "Combustion-drived PCA")),
                      labels=c("Combustion-derived PCA","Polyphenols and PCAs with aliphatic chains",
                               "Humics and highly unsaturated compounds",
                               "Unsaturated aliphatic compounds",
                               "Saturated fatty and sulfonic acid, carbohyrates"
                      ))


pie_all=as.data.table(aggregate(ft_data$value, 
                                by=list(`Group`=ft_data$Group,type=ft_data$type,
                                        Class=ft_data$Molecularclass), FUN=sum))

pie_all_tot=as.data.table(aggregate(ft_data$value, 
                                    by=list(`Group`=ft_data$Group,type=ft_data$type), FUN=sum)) %>% `colnames<-`(c("Group","type","tot"))

pie_all=pie_all %>% inner_join(pie_all_tot)

pie_all$rel=pie_all$x/pie_all$tot*100
pie_all

pie_all$Class2=factor(pie_all$Class,levels = rev(c("Saturated fatty and sulfonic acid,carbohyrates","Unsaturated aliphatic compounds",
                                                   "Soil derived humics and highly unsaturated compounds",
                                                   "Soil-derived polyphenols and PCAs with aliphatic chains",
                                                   "Combustion-drived PCA")),
                      labels=c("Combustion-derived PCA","Polyphenols and PCAs with aliphatic chains",
                               "Humics and highly unsaturated compounds",
                               "Unsaturated aliphatic compounds",
                               "Saturated fatty and sulfonic acid, carbohyrates"
                      ))

#pie_inty=pie_all[,c(1,2,5)]
pie_frq=pie_all[,c(1,2,5)]

pie_frq$inty=round(pie_inty$rel,1)

pie_frq


ggplot()+
  geom_bar(data=pie_all,aes(x="",y=x, fill=Class2), col="black",stat="identity", size=0.5,
           position = position_fill(reverse = T))+
  coord_polar("y")+
  scale_fill_manual(values = rev(c("#586957","#FDB2AD","#4F94CD","#FAD692","#9D95CA")))+
  #scale_fill_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF","grey50"))+
  facet_wrap(type~Group, ncol = 5)+
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "transparent", color = NA, size = 0), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA, size=0), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    strip.text = element_text(size = 20, family = "Arial",face = "bold"),
    legend.title = element_text(size = 16, vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
    legend.text = element_text(size = 16, family = "Arial",face = "bold",colour = "black",hjust = 0.0,margin = unit(c(0,0.1,0,0.1),"cm")),
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent", color = NA, size=0), # get rid of legend panel bg
    legend.key.height = unit(0.9,"cm"),
    legend.key.width = unit(0.9,"cm"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  )+
  guides(fill=guide_legend(reverse = F, ncol = 3, byrow=T, title = NULL, title.position = "top"))+
  ggsave("220519pie_all_trfu_legend.png",height = 24, width = 50, units = "cm", dpi = 300)

tt=pie_all
tt$rel=round(tt$rel,1)
t2=subset(tt,tt$Group=="Noto")
t2[order(t2$type),1:6]

