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

getOption("digits")
options("digits" = 15)

source("Rscript/func_filename.R")
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")

insert_minor_1 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_2 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}


ft_merge=fread("Datafile/FRIENDs_1st_FT.csv")
ft_merge

ft_merge=ft_merge %>% mutate(Comp = gsub('[[:digit:]]', '', Formula))
ft_merge$Comp=ifelse(ft_merge$O.==0,"Remainders",ft_merge$Comp)

table(ft_merge$Group)
ft_merge$Group=factor(ft_merge$Group, levels =c("UL","B","SS","SUL","NT"),
                      labels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))
ft_merge$Sample=paste(ft_merge$Group,ft_merge$No, sep = "_")

ft_merge$Freq=1

fm_obs=ft_merge[,c("Group","No","Formula","Freq")] %>% 
  dcast(Group~Formula,value.var = "Freq",sum) %>% 
  melt(id.vars=c("Group")) %>% `colnames<-`(c("Group","Formula","cnt"))

fm_obs
ft_merge=ft_merge %>% inner_join(fm_obs)

ft_merge_sel=subset(ft_merge,ft_merge$cnt>3)
ft_merge_sel

###build cor data=====
ft_merge_cor=ft_merge

ft_inty_all=as.data.table(aggregate(ft_merge$Bromo.Inty, by=list(`Sample`=ft_merge$Sample), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Tot"))

ft_merge_cor=ft_merge_cor %>% inner_join(ft_inty_all,by = "Sample")
ft_merge_cor$rel=ft_merge_cor$Bromo.Inty/ft_merge_cor$Tot*100

ft_merge_cor_sel=subset(ft_merge_cor,ft_merge_cor$cnt>6)

ft_merge_cor_sel_m=melt(ft_merge_cor_sel[,c("Sample","Group","Formula","rel")],id.vars = c("Sample","Group","Formula")) %>% 
  dcast(Sample~Formula, value.var = "value",fun.aggregate = sum)
ft_merge_cor_sel_m[,1:12]

FT_envi=fread("Datafile/FRIEND_1st_envi_re.csv")
FT_envi$WSOCbb=2.94*FT_envi$Levoglucosan
FT_envi$WSOCnbb=FT_envi$WSOC-FT_envi$WSOCbb

FT_envi$Sample
ft_merge_cor_sel_m$Sample


FT_fm=FT_envi[,c("Sample","Group","No","WSOCbb","WSOCnbb")] %>% inner_join(ft_merge_cor_sel_m)

FT_fm[,1:12]

grp=as.vector(unique(FT_fm$Group))
df=data.frame()
for (i in 2:2) {
  #  i=1
  temp=as.data.frame(subset(FT_fm,FT_fm$Group==grp[i]))
  temp[,1:12]
  for (j in 4:5) {
    #j=4
    for (k in 6:dim(temp)[2]) {
      #k=10  
      sum.j<-sum(temp[,j], na.rm = T)
      sum.k<-sum(temp[,k])
      
      if(sum.j >0 & sum.k >0){
        test<-cor.test(temp[,j],temp[,k],method="spearman",na.action=na.rm,exact = F, )
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

df
table(df$Envi)
table(df$Group)

dim(FT_fm)
cor_vk_1st=fread("Datafile/corDB_WSOCBB&formula.csv")
fm_obs

cor_vk_1st_cnt=cor_vk_1st %>% inner_join(fm_obs)
cor_vk_1st
cor_vk_1st_cnt_sel=subset(cor_vk_1st_cnt,cor_vk_1st_cnt$cnt>3)
#fwrite(cor_vk_1st_cnt_sel,file = "Datafile/corDB_WSOCBB&formula_obs3.csv")

cor_vk_1st_cnt_sel




cor_vk_d=dcast(cor_vk_1st_cnt_sel, Group+Formula~Envi, sum, value.var = c("rho","p"))
cor_vk_d

cor_vk_d$Type=ifelse(abs(cor_vk_d$rho_WSOCbb)>abs(cor_vk_d$rho_WSOCnbb),"BB","NBB")
table(cor_vk_d$Type)
cor_vk_d

cor_vk_d$Type=ifelse(abs(cor_vk_d$rho_WSOCbb)>abs(cor_vk_d$rho_WSOCnbb),"BB","NBB")
table(cor_vk_d$Type)

cor_vk_d$variable=ifelse(cor_vk_d$Type=="BB","WSOCbb","WSOCnbb")
cor_vk_d$value=ifelse(cor_vk_d$Type=="BB",cor_vk_d$rho_WSOCbb,cor_vk_d$rho_WSOCnbb)
cor_vk_d$p=ifelse(cor_vk_d$Type=="BB",cor_vk_d$p_WSOCbb,cor_vk_d$p_WSOCnbb)
cor_vk_d

cor_vk_d$Grouplab=factor(cor_vk_d$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                         labels=c(expression(bold("Ulaanbaatar")),
                                  expression(bold("Beijing")),
                                  expression(bold("Seosan")),
                                  expression(bold("Seoul")),
                                  expression(bold("Noto"))))

cor_vk_d$varlab=factor(cor_vk_d$variable,levels = c("WSOCbb","WSOCnbb"),
                       labels = c(expression(bold("WSOC"["BB"])),
                                  expression(bold("WSOC"["NBB"]))))

cor_vk_d


#display as vkpot=====
ft_merge
ft_merge=ft_merge %>% mutate(Comp = gsub('[[:digit:]]', '', Formula))
ft_merge$Comp=ifelse(ft_merge$O.==0,"Remainders",ft_merge$Comp)

ft_vk=ft_merge[,c("Group","Formula","O.C","H.C","Bromo.Inty")]
ft_vk_m=dcast(ft_vk,Group~Formula, mean) %>% 
  melt(id.vars=c("Group"), na.rm = T) %>% 
  `colnames<-`(c("Group","Formula","inty"))

ft_vk_m
cor_vk_d
dim(ft_vk_m)

cor_vk_d_sel=subset(cor_vk_d,cor_vk_d$p<0.05)
cor_vk_d_sel

dim(ft_vk_m)

ft_cor_vk=ft_vk_m %>% left_join(cor_vk_d_sel[,c("Group","Formula","variable","value","p")])
ft_cor_vk$variable2=ifelse(is.na(ft_cor_vk$variable),"Non-sig",ft_cor_vk$variable)
ft_cor_vk

table(ft_cor_vk$variable2)

CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(ft_cor_vk$Formula)
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
        if (clasul(numElement)=="integer"){
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


ft_cor_vk
ft_cor_vk=ft_cor_vk %>% mutate(Comp = gsub('[[:digit:]]', '', Formula))
table(ft_cor_vk$Comp)

fmlist=as.data.table(numericalFormula) %>% `colnames<-`(c("C","H","N","O","S"))
ft_cor_vk0=ft_cor_vk

ft_cor_vk=ft_cor_vk %>% inner_join(fmlist)

ft_cor_vk$`O/C`=ft_cor_vk$O/ft_cor_vk$C
ft_cor_vk$`H/C`=ft_cor_vk$H/ft_cor_vk$C
ft_cor_vk

ft_cor_vk$Comp=ifelse(ft_cor_vk$`O`==0,"Remainders",ft_cor_vk$Comp)
ft_cor_vk$Comp=factor(ft_cor_vk$Comp,levels = c("CHO","CHNO","CHOS","CHNOS","Remainders"),
                      labels = c("CHO","CHON","CHOS","CHONS","Remainders"))

table(ft_cor_vk$Comp)

ft_cor_vk$Grouplab=factor(ft_cor_vk$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"),
                          labels=c(expression(bold("Ulaanbaatar")),
                                   expression(bold("Beijing")),
                                   expression(bold("Seosan")),
                                   expression(bold("Seoul")),
                                   expression(bold("Noto"))))

ft_cor_vk

ft_cor_vk
ft_cor_vk_sel=subset(ft_cor_vk,ft_cor_vk$variable2!="Non-sig") %>% droplevels()
ft_cor_vk_sel

ft_cor_vk
ft_cor_vk_sel$varlab=factor(ft_cor_vk_sel$variable2,levels = c("WSOCbb","WSOCnbb"),
                            labels = c(expression(bold("WSOC"["BB"])),
                                       expression(bold("WSOC"["NBB"]))))

ft_cor_vk_sel
ft_cor_vk_sel1=subset(ft_cor_vk_sel,ft_cor_vk_sel$p<0.05)


ft_cor_vk_non=subset(ft_cor_vk,ft_cor_vk$variable2=="Non-sig") %>% droplevels()

ft_cor_vk_non1=ft_cor_vk_non
ft_cor_vk_non1$variable2="WSOCbb"

ft_cor_vk_non2=ft_cor_vk_non
ft_cor_vk_non2$variable2="WSOCnbb"

ft_cor_vk_non1
ft_cor_vk_non2

ft_cor_vk_non=rbind.data.frame(ft_cor_vk_non1,ft_cor_vk_non2)

ft_cor_vk_non=ft_cor_vk_non %>% inner_join(fm_obs)

ft_cor_vk_non_sel=subset(ft_cor_vk_non,ft_cor_vk_non$cnt>4)

##Ul=====
ft_cor_ul=subset(ft_cor_vk_sel1,ft_cor_vk_sel1$Group=="Ulaanbaatar")

ft_cor_ul$variable2=factor(ft_cor_ul$variable,levels=c("WSOCnbb","WSOCbb"),
                           labels = c("WSOCbb","WSOCnbb"))

ft_cor_ul$varlab=factor(ft_cor_ul$variable,levels = c("WSOCnbb","WSOCbb"),
                        labels = c(expression(bold("WSOC"["BB"])),
                                   expression(bold("WSOC"["NBB"]))))
ft_cor_ul
table(ft_cor_ul$varlab)

ft_cor_ul_sel=subset(ft_cor_ul,ft_cor_ul$p<0.05)
ft_cor_ul_sel

table(ft_cor_ul_sel$varlab)
ft_cor_ul

tt=as.data.frame(table(ft_cor_ul$variable2))
tt$varlab=tt$Var1
tt$varlab=factor(tt$varlab,levels = c("WSOCbb","WSOCnbb"),
                 labels = c(expression(bold("WSOC"["BB"])),
                            expression(bold("WSOC"["NBB"]))))

tt$O.C=0.3
tt$H.C=0.15

ft_cor_vk_non

ft_cor_ul_non=subset(ft_cor_vk_non_sel,ft_cor_vk_non_sel$Group=="Ulaanbaatar")

ggplot()+
  geom_point(data=ft_cor_ul_non, aes(x=`O/C`, y=`H/C`),col="grey70",size=1.5)+
  geom_point(data=ft_cor_ul_sel, aes(x=`O/C`, y=`H/C`,col=value),size=2.7)+
  #geom_text(data=tt, aes(x=O.C, y=H.C,label=Freq), size=8)+
  geom_text(data=tt, aes(x=0.05, y=0.14,label=paste0(expression(italic("N =")))),parse=T, size=8)+
  geom_text(data=tt, aes(x=0.18, y=H.C,label=Freq),parse=T, size=8)+
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
  ggsave(filename("vk_cor_u_rel"),height = 40, width = 20, units = "cm", dpi = 300)

##Beijing=====
ft_cor_bj=subset(ft_cor_vk_sel1,ft_cor_vk_sel1$Group=="Beijing")

ft_cor_bj$varlab=factor(ft_cor_bj$variable,levels = c("WSOCbb","WSOCnbb"),
                        labels = c(expression(bold("WSOC"["BB"])),
                                   expression(bold("WSOC"["NBB"]))))
ft_cor_bj
table(ft_cor_bj$varlab)

ft_cor_bj_sel=subset(ft_cor_bj,ft_cor_bj$p<0.05)
ft_cor_bj_sel

table(ft_cor_bj_sel$varlab)
ft_cor_bj

tt_bj=as.data.frame(table(ft_cor_bj$variable2))
tt_bj$varlab=tt_bj$Var1
tt_bj$varlab=factor(tt_bj$varlab,levels = c("WSOCbb","WSOCnbb"),
                    labels = c(expression(bold("WSOC"["BB"])),
                               expression(bold("WSOC"["NBB"]))))

tt_bj$O.C=0.1
tt_bj$H.C=0.15

ft_cor_vk_non

ft_cor_bj_non=subset(ft_cor_vk_non_sel,ft_cor_vk_non_sel$Group=="Beijing")

ggplot()+
  geom_point(data=ft_cor_bj_non, aes(x=`O/C`, y=`H/C`),col="grey70",size=1.5)+
  geom_point(data=ft_cor_bj_sel, aes(x=`O/C`, y=`H/C`,col=value),size=2.7)+
  #geom_text(data=tt, aes(x=O.C, y=H.C,label=Freq), size=8)+
  geom_text(data=tt_bj, aes(x=0.05, y=0.14,label=paste0(expression(italic("N =")))),parse=T, size=8)+
  geom_text(data=tt_bj, aes(x=0.18, y=H.C,label=Freq),parse=T, size=8)+
  #facet_grid(Grouplab~varlab, labeller = label_parsed )+
  facet_grid(varlab~Grouplab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
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
        strip.text.x = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.6,0,0.2,0),"cm")),
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
  ggsave(filename("vk_cor_bj_rel"),height = 40, width = 20, units = "cm", dpi = 300)

##Seosan=====
ft_cor_ss=subset(ft_cor_vk_sel1,ft_cor_vk_sel1$Group=="Seosan")

ft_cor_ss$varlab=factor(ft_cor_ss$variable,levels = c("WSOCbb","WSOCnbb"),
                        labels = c(expression(bold("WSOC"["BB"])),
                                   expression(bold("WSOC"["NBB"]))))
ft_cor_ss

ft_cor_ss=ft_cor_ss[order(ft_cor_ss$value),]

table(ft_cor_ss$varlab)

ft_cor_ss_sel=subset(ft_cor_ss,ft_cor_ss$p<0.05)
ft_cor_ss_sel

table(ft_cor_ss_sel$varlab)
ft_cor_ss_sel

tt_ss=as.data.frame(table(ft_cor_ss_sel$variable2))
tt_ss$varlab=tt_ss$Var1
tt_ss$varlab=factor(tt_ss$varlab,levels = c("WSOCbb","WSOCnbb"),
                    labels = c(expression(bold("WSOC"["BB"])),
                               expression(bold("WSOC"["NBB"]))))

tt_ss$O.C=0.1
tt_ss$H.C=0.15

ft_cor_vk_non
ft_cor_ss_non=subset(ft_cor_vk_non_sel,ft_cor_vk_non_sel$Group=="Seosan")

ggplot()+
  geom_point(data=ft_cor_ss_non, aes(x=`O/C`, y=`H/C`),col="grey70",size=1.5)+
  geom_point(data=ft_cor_ss_sel, aes(x=`O/C`, y=`H/C`,col=value),size=2.7)+
  #geom_text(data=tt, aes(x=O.C, y=H.C,label=Freq), size=8)+
  geom_text(data=tt_ss, aes(x=0.05, y=0.14,label=paste0(expression(italic("N =")))),parse=T, size=8)+
  geom_text(data=tt_ss, aes(x=0.18, y=H.C,label=Freq),parse=T, size=8)+
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
  ggsave(filename("vk_cor_ss_rel"),height = 40, width = 20, units = "cm", dpi = 300)

#Seoul=====
ft_cor_sul=subset(ft_cor_vk_sel1,ft_cor_vk_sel1$Group=="Seoul")

ft_cor_sul$varlab=factor(ft_cor_sul$variable,levels = c("WSOCbb","WSOCnbb"),
                         labels = c(expression(bold("WSOC"["BB"])),
                                    expression(bold("WSOC"["NBB"]))))
ft_cor_sul
table(ft_cor_sul$varlab)

ft_cor_sul_sel=subset(ft_cor_sul,ft_cor_sul$p<0.05)
ft_cor_sul_sel

table(ft_cor_sul_sel$varlab)
ft_cor_sul

tt_sul=as.data.frame(table(ft_cor_sul$variable2))
tt_sul$varlab=tt_sul$Var1
tt_sul$varlab=factor(tt_sul$varlab,levels = c("WSOCbb","WSOCnbb"),
                     labels = c(expression(bold("WSOC"["BB"])),
                                expression(bold("WSOC"["NBB"]))))

tt_sul$O.C=0.1
tt_sul$H.C=0.15

ft_cor_vk_non
ft_cor_sul_non=subset(ft_cor_vk_non_sel,ft_cor_vk_non_sel$Group=="Seoul")

ggplot()+
  geom_point(data=ft_cor_sul_non, aes(x=`O/C`, y=`H/C`),col="grey70",size=1.5)+
  geom_point(data=ft_cor_sul_sel, aes(x=`O/C`, y=`H/C`,col=value),size=2.7)+
  #geom_text(data=tt, aes(x=O.C, y=H.C,label=Freq), size=8)+
  geom_text(data=tt_sul, aes(x=0.05, y=0.14,label=paste0(expression(italic("N =")))),parse=T, size=8)+
  geom_text(data=tt_sul, aes(x=0.18, y=H.C,label=Freq),parse=T, size=8)+
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
  ggsave(filename("vk_cor_sul_rel"),height = 40, width = 20, units = "cm", dpi = 300)


#noto=====
ft_cor_nt=subset(ft_cor_vk_sel1,ft_cor_vk_sel1$Group=="Noto")

ft_cor_nt$varlab=factor(ft_cor_nt$variable,levels = c("WSOCbb","WSOCnbb"),
                        labels = c(expression(bold("WSOC"["BB"])),
                                   expression(bold("WSOC"["NBB"]))))
ft_cor_nt
table(ft_cor_nt$varlab)

ft_cor_nt_sel=subset(ft_cor_nt,ft_cor_nt$p<0.05)
ft_cor_nt_sel

table(ft_cor_nt_sel$varlab)
ft_cor_nt

tt_nt=as.data.frame(table(ft_cor_nt$variable2))
tt_nt$varlab=tt_nt$Var1
tt_nt$varlab=factor(tt_nt$varlab,levels = c("WSOCbb","WSOCnbb"),
                    labels = c(expression(bold("WSOC"["BB"])),
                               expression(bold("WSOC"["NBB"]))))

tt_nt$O.C=0.1
tt_nt$H.C=0.15

ft_cor_vk_non
ft_cor_nt_non=subset(ft_cor_vk_non_sel,ft_cor_vk_non_sel$Group=="Noto")

ggplot()+
  geom_point(data=ft_cor_nt_non, aes(x=`O/C`, y=`H/C`),col="grey70",size=1.5)+
  geom_point(data=ft_cor_nt_sel, aes(x=`O/C`, y=`H/C`,col=value),size=2.7)+
  #geom_text(data=tt, aes(x=O.C, y=H.C,label=Freq), size=8)+
  geom_text(data=tt_nt, aes(x=0.05, y=0.14,label=paste0(expression(italic("N =")))),parse=T, size=8)+
  geom_text(data=tt_nt, aes(x=0.18, y=H.C,label=Freq),parse=T, size=8)+
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
  ggsave(filename("vk_cor_nt_rel"),height = 40, width = 20, units = "cm", dpi = 300)

###
ft_cor_ul$variable=ft_cor_ul$variable2
ft_cor_bj
ft_cor_ss
ft_cor_sul
ft_cor_nt

ft_cor_vk_m=rbind(ft_cor_ul,ft_cor_bj,ft_cor_ss,ft_cor_sul,ft_cor_nt)

ft_cor_vk_pos=ft_cor_vk_m %>% filter(value>0)
ft_cor_vk_pos$Freq=1
wsoc_bb_pos=as.data.table(aggregate(ft_cor_vk_pos$Freq, by=list(Group=ft_cor_vk_pos$Group,type=ft_cor_vk_pos$variable),sum))
wsoc_bb_pos

ft_cor_vk_neg=ft_cor_vk_m %>% filter(value<0)
ft_cor_vk_neg$Freq=1
wsoc_bb_neg=as.data.table(aggregate(ft_cor_vk_neg$Freq, by=list(Group=ft_cor_vk_neg$Group,type=ft_cor_vk_neg$variable),sum))

wsoc_bb_pos
wsoc_bb_neg

wsoc_bb_pos$neg=wsoc_bb_neg$x
wsoc_bb_pos$sum=wsoc_bb_pos$x+wsoc_bb_pos$neg
wsoc_bb_pos

ft_cor_vk_m

fm_chp=unique(ft_merge[,c("Formula","AI","DBE","O.C","H.C")])

ft_cor_vk_m=ft_cor_vk_m %>% inner_join(fm_chp)



ft_cor_chp=melt(ft_cor_vk_m[,c("Group","variable","AI","DBE","O.C","H.C")],id.vars = c("Group","variable")) %>% 
  `colnames<-`(c("Group","Type","variable","val"))
ft_cor_chp

ft_merge_soa=ft_merge %>% inner_join(ft_cor_vk_m[,c("Group","Formula","variable")])

table(ft_merge_soa$variable)
ft_merge_soa

ft_merge_soa_m=melt(ft_merge_soa[,c("Group","No","variable","AI","DBE","O.C","H.C")],
                    id.vars = c("Group","No","variable")) %>% 
  `colnames<-`(c("Group","No","Type","variable","val")) %>% 
  dcast(Group+No+Type~variable, mean) %>% 
  melt(id.vars=c("Group","No","Type")) %>% 
  `colnames<-`(c("Group","No","Type","variable","val"))

ft_merge_soa_m

ft_merge_soa_m$variable=factor(ft_merge_soa_m$variable,levels = c("AI","DBE","O.C","H.C"),
                               labels = c("AI","DBE","O/C", "H/C"))
ft_merge_soa_m$Group=factor(ft_merge_soa_m$Group,levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))

ft_merge_soa_m




ggplot()+
  stat_boxplot(data=ft_merge_soa_m, aes(x=Group, y=val, fill=Type),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=ft_merge_soa_m, aes(x=Group, y=val, fill=Type),alpha=1, outlier.color = NA)+
  geom_text(data =ft_merge_soa_m, aes(x=1, y=val*1.1, label=variable), col="white", size=0)+
  #geom_text(data =stat_chp, aes(x=Group, y=ifelse(variable=="Mean O/C",val+0.010,
  #                                                ifelse(variable=="Mean H/C",val+0.030,
  #                                                       ifelse(variable=="Mean N/C",val+0.006,
  #                                                              ifelse(variable=="Mean S/C",val+0.006,
  #                                                                     ifelse(variable=="Mean DBE",val+0.45,val+0.012))))), label=labels), 
  #          col="black", size=6)+
  facet_rep_wrap(.~variable, nrow=2, repeat.tick.labels = T, strip.position="left", scales="free")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  # scale_fill_manual(values = c("#9DBCD4","#CB7723"))+
  scale_fill_manual(values = c("#9DBCD4","#CB7723"),
                    labels=c(expression(bold("WSOC"["BB"])),
                             expression(bold("WSOC"["NBB"]))))+
  #  scale_fill_manual(values = c( "#3C5488FF","#E64B35FF","#4DBBD5FF", "#00A087FF","#F39B7FFF"))+
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
        axis.text.y = element_text(size = 14, colour = "black", margin = unit(c(0.0,0.0,0.2,0.2),"cm")),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        #axis.title.y.left = element_text(size = 0, colour = "black",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        #axis.title.y.right = element_text(size = 0, colour = "black",angle = 90,margin = unit(c(0.2,0.2,0.2,0.5),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.placement = "outside",
        strip.text.y = element_text(size = 16, colour = "black",face = 2,margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        strip.text.x = element_blank(),
        #legend.margin = unit(c(0.0,0.0,0.0,0.0),"cm"),
        legend.title = element_text(size = 22, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 18,hjust = 0, colour = "black",margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        legend.key.width = unit(1.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.background = element_blank(),
        legend.direction = "vertical",
        legend.justification=c(0, 1),
        legend.position = c(0.55,0.5))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("wsocbb_Chp_compare_fdr"),height = 20, width = 30, units = "cm", dpi = 300)

dt

ft_merge_soa_m


chp=unique(ft_merge_soa_m$variable)
grp=unique(ft_merge_soa_m$Group)
type=unique(ft_merge_soa_m$Type)

dt=data.table()
for (i in 1:length(chp)) {
  #i=1
  temp=subset(ft_merge_soa_m,ft_merge_soa_m$variable==chp[i])
  
  for (j in 1:length(grp)) {
    #j=1
    temp2=subset(temp,temp$Group==grp[j])
    
    tbb=subset(temp2,temp2$Type=="WSOCbb")
    tnbb=subset(temp2,temp2$Type=="WSOCnbb")
    
    tt=wilcox.test(tbb$val,tnbb$val,exact=F)
    
    new=data.table(Group=grp[j],var=chp[i],t=tt$statistic, p=round(tt$p.value,3))
    
    dt=rbind(dt,new)
  }
  
  
}
dt



chp2=unique(ft_cor_chp$variable)
grp2=unique(ft_cor_chp$Group)
type2=unique(ft_cor_chp$Type)

dt2=data.table()
for (i in 1:length(chp2)) {
  #i=1
  temp=subset(ft_cor_chp,ft_cor_chp$variable==chp2[i])
  
  for (j in 1:length(grp2)) {
    #j=1
    temp2=subset(temp,temp$Group==grp2[j])
    
    tbb=subset(temp2,temp2$Type=="WSOCbb")
    tnbb=subset(temp2,temp2$Type=="WSOCnbb")
    
    tt=wilcox.test(tbb$val,tnbb$val,exact=F)
    
    new=data.table(Group=grp2[j],var=chp2[i],t=tt$statistic, p=round(tt$p.value,3))
    
    dt2=rbind(dt2,new)
  }
  
  
}
dt2

ft_cor_chp
ft_cor_chp
ft_cor_chp$Group=factor(ft_cor_chp$Group,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))




ggplot()+
  stat_boxplot(data=ft_cor_chp, aes(x=Group, y=val, fill=Type),geom='errorbar', linetype=1, width=0.25,
               position = position_dodge(width = 0.75))+
  geom_boxplot(data=ft_cor_chp, aes(x=Group, y=val, fill=Type),alpha=1, outlier.color = NA)+
  geom_text(data =ft_cor_chp, aes(x=1, y=val*1.1, label=variable), col="white", size=0)+
  #geom_text(data =stat_chp, aes(x=Group, y=ifelse(variable=="Mean O/C",val+0.010,
  #                                                ifelse(variable=="Mean H/C",val+0.030,
  #                                                       ifelse(variable=="Mean N/C",val+0.006,
  #                                                              ifelse(variable=="Mean S/C",val+0.006,
  #                                                                     ifelse(variable=="Mean DBE",val+0.45,val+0.012))))), label=labels), 
  #          col="black", size=6)+
  facet_rep_wrap(.~variable, nrow=2, repeat.tick.labels = T, strip.position="left", scales="free")+
  scale_y_continuous(name = "",
                     expand = c(0.02,0.02))+
  # scale_fill_manual(values = c("#9DBCD4","#CB7723"))+
  scale_fill_manual(values = c("#9DBCD4","#CB7723"),
                    labels=c(expression(bold("WSOC"["BB"])),
                             expression(bold("WSOC"["NBB"]))))+
  #  scale_fill_manual(values = c( "#3C5488FF","#E64B35FF","#4DBBD5FF", "#00A087FF","#F39B7FFF"))+
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
        axis.text.y = element_text(size = 14, colour = "black", margin = unit(c(0.0,0.0,0.2,0.2),"cm")),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        #axis.title.y.left = element_text(size = 0, colour = "black",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        #axis.title.y.right = element_text(size = 0, colour = "black",angle = 90,margin = unit(c(0.2,0.2,0.2,0.5),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.placement = "outside",
        strip.text.y = element_text(size = 16, colour = "black",face = 2,margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        strip.text.x = element_blank(),
        #legend.margin = unit(c(0.0,0.0,0.0,0.0),"cm"),
        legend.title = element_text(size = 22, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 18,hjust = 0, colour = "black",margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        legend.key.width = unit(1.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.background = element_blank(),
        legend.direction = "vertical",
        legend.justification=c(0, 1),
        legend.position = c(0.55,0.5))+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("wsocbb_Chp_compare2"),height = 20, width = 30, units = "cm", dpi = 300)




