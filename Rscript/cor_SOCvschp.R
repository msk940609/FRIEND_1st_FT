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

ft_envi=fread("Datafile/FRIEND_1st_envi_re3.csv")
FT_envi_sel=ft_envi[,c("Sample","Group","No","Event","Date","PM2.5","OC","WSOC","WISOC","WSOCbb","WSOCnbb")]
FT_envi_sel$POC=FT_envi_sel$WISOC+FT_envi_sel$WSOCbb
FT_envi_sel$SOC=FT_envi_sel$WSOCnbb

FT_envi_sel$POCp=FT_envi_sel$POC/FT_envi_sel$OC*100
FT_envi_sel$SOCp=FT_envi_sel$SOC/FT_envi_sel$OC*100

ft_merge=fread("Datafile/FRIENDs_1st_FT.csv")
table(ft_merge$Group)
ft_merge$Group=factor(ft_merge$Group, levels =c("UL","B","SS","SUL","NT"),
                      labels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))

ft_merge=ft_merge %>% mutate(Comp = gsub('[[:digit:]]', '', Formula))
ft_merge

ft_merge$Sample=paste(ft_merge$Group,ft_merge$No, sep = "_")
ft_merge$Comp=ifelse(ft_merge$`O.`==0,"Remainders",ft_merge$Comp)

ft_merge$Comp=factor(ft_merge$Comp,levels = c("CHO","CHNO","CHOS","CHNOS","Remainders"),
                     labels = c("CHO","CHON","CHOS","CHONS","Remainders"))

ft_merge
ft_merge$Freq=1
fm_obs=ft_merge[,c("Group","No","Formula","Freq")] %>% 
  dcast(Group~Formula,value.var = "Freq",sum) %>% 
  melt(id.vars=c("Group")) %>% `colnames<-`(c("Group","Formula","cnt"))

ft_merge_filter=ft_merge %>% inner_join(fm_obs)

ft_merge_filter=subset(ft_merge_filter,ft_merge_filter$cnt>4)

ft_merge=ft_merge_filter

ft_merge$Freq=1


ft_cor_pos

ft_merge_soa=ft_merge %>% left_join(ft_cor_pos[,c("Group","Formula","type")])
ft_merge_soa

ft_merge_soa_nbb=ft_merge_soa %>% filter(type%in%c("WSOCnbb"))
ft_merge_soa_nbb=subset(ft_merge_soa_nbb,ft_merge_soa_nbb$Comp!="Remainders") %>% droplevels()

ft_merge_soa_nbb_chp=melt(ft_merge_soa_nbb[,c("Sample","Group","No","Comp","AI","DBE","O.C","H.C")],
                          id.vars = c("Sample","Group","No","Comp"),measure.vars = c("AI","DBE","O.C","H.C"))%>% 
  dcast(Sample+Group+No~Comp+variable,mean)

ft_merge_soa_nbb_chp



FT_envi_sel
ft_merge_soa_nbb_chp_soc=ft_merge_soa_nbb_chp %>% inner_join(FT_envi_sel[,c("Sample","WSOCnbb")])


source("Rscript/func_chartcor.R")
library(PerformanceAnalytics)
ft_comp=as.data.table(aggregate(ft_merge$Freq, by=list(Sample=ft_merge$Sample, Comp=ft_merge$Comp),sum)) %>% 
  `colnames<-`(c("Sample","Comp","cnt"))

ft_comp_tot=as.data.table(aggregate(ft_merge$Freq, by=list(Sample=ft_merge$Sample),sum)) %>% 
  `colnames<-`(c("Sample","tot"))

ft_comp=ft_comp %>% inner_join(ft_comp_tot)

ft_comp$rel=round(ft_comp$cnt/ft_comp$tot*100,1)
ft_comp_sel=ft_comp %>% filter(!Comp%in%c("Remainders"))
ft_comp_sel

ft_comp_sel_d=dcast(ft_comp_sel,Sample~Comp, value.var = "rel",sum)
ft_merge_soa_nbb_chp_soc=ft_merge_soa_nbb_chp_soc %>% inner_join(ft_comp_sel_d)

source("Rscript/func_chartcor.R")
library(PerformanceAnalytics)

##Seoul=====
ft_merge_soa_nbb_chp_soc$SOC=ft_merge_soa_nbb_chp_soc$WSOCnbb

ft_merge_soa_nbb_chp_soc_sul=subset(ft_merge_soa_nbb_chp_soc,ft_merge_soa_nbb_chp_soc$Group=="Seoul")
ft_merge_soa_nbb_chp_soc_sul_cho=ft_merge_soa_nbb_chp_soc_sul[,c("SOC","CHO","CHO_AI","CHO_DBE","CHO_O.C","CHO_H.C")]

tiff(filename = "sul_cor_chp_cho.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_sul_cho,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_sul_chon=ft_merge_soa_nbb_chp_soc_sul[,c("SOC","CHON","CHON_AI","CHON_DBE","CHON_O.C","CHON_H.C")]

tiff(filename = "sul_cor_chp_chon.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_sul_chon,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_sul_CHOS=ft_merge_soa_nbb_chp_soc_sul[,c("SOC","CHOS","CHOS_AI","CHOS_DBE","CHOS_O.C","CHOS_H.C")]
tiff(filename = "sul_cor_chp_CHOS.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_sul_CHOS,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_sul_CHONS=ft_merge_soa_nbb_chp_soc_sul[,c("SOC","CHONS","CHONS_AI","CHONS_DBE","CHONS_O.C","CHONS_H.C")]
tiff(filename = "sul_cor_chp_CHONS.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_sul_CHONS,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

###Ul=====
ft_merge_soa_nbb_chp_soc_ul=subset(ft_merge_soa_nbb_chp_soc,ft_merge_soa_nbb_chp_soc$Group=="Ulaanbaatar")
ft_merge_soa_nbb_chp_soc_ul_cho=ft_merge_soa_nbb_chp_soc_ul[,c("SOC","CHO","CHO_AI","CHO_DBE","CHO_O.C","CHO_H.C")]

tiff(filename = "ul_cor_chp_cho.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_ul_cho,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_ul_chon=ft_merge_soa_nbb_chp_soc_ul[,c("SOC","CHON","CHON_AI","CHON_DBE","CHON_O.C","CHON_H.C")]

tiff(filename = "ul_cor_chp_chon.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_ul_chon,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_ul_CHOS=ft_merge_soa_nbb_chp_soc_ul[,c("SOC","CHOS","CHOS_AI","CHOS_DBE","CHOS_O.C","CHOS_H.C")]
tiff(filename = "ul_cor_chp_CHOS.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_ul_CHOS,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_ul_CHONS=ft_merge_soa_nbb_chp_soc_ul[,c("SOC","CHONS","CHONS_AI","CHONS_DBE","CHONS_O.C","CHONS_H.C")]
tiff(filename = "ul_cor_chp_CHONS.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_ul_CHONS,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

###Bj=====
ft_merge_soa_nbb_chp_soc_bj=subset(ft_merge_soa_nbb_chp_soc,ft_merge_soa_nbb_chp_soc$Group=="Beijing")
ft_merge_soa_nbb_chp_soc_bj_cho=ft_merge_soa_nbb_chp_soc_bj[,c("SOC","CHO","CHO_AI","CHO_DBE","CHO_O.C","CHO_H.C")]

tiff(filename = "bj_cor_chp_cho.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_bj_cho,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_bj_chon=ft_merge_soa_nbb_chp_soc_bj[,c("SOC","CHON","CHON_AI","CHON_DBE","CHON_O.C","CHON_H.C")]

tiff(filename = "bj_cor_chp_chon.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_bj_chon,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_bj_CHOS=ft_merge_soa_nbb_chp_soc_bj[,c("SOC","CHOS","CHOS_AI","CHOS_DBE","CHOS_O.C","CHOS_H.C")]
tiff(filename = "bj_cor_chp_CHOS.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_bj_CHOS,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_bj_CHONS=ft_merge_soa_nbb_chp_soc_bj[,c("SOC","CHONS","CHONS_AI","CHONS_DBE","CHONS_O.C","CHONS_H.C")]
tiff(filename = "bj_cor_chp_CHONS.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_bj_CHONS,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

###SS=====
ft_merge_soa_nbb_chp_soc_ss=subset(ft_merge_soa_nbb_chp_soc,ft_merge_soa_nbb_chp_soc$Group=="Seosan")
ft_merge_soa_nbb_chp_soc_ss_cho=ft_merge_soa_nbb_chp_soc_ss[,c("SOC","CHO","CHO_AI","CHO_DBE","CHO_O.C","CHO_H.C")]

tiff(filename = "ss_cor_chp_cho.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_ss_cho,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_ss_chon=ft_merge_soa_nbb_chp_soc_ss[,c("SOC","CHON","CHON_AI","CHON_DBE","CHON_O.C","CHON_H.C")]

tiff(filename = "ss_cor_chp_chon.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_ss_chon,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_ss_CHOS=ft_merge_soa_nbb_chp_soc_ss[,c("SOC","CHOS","CHOS_AI","CHOS_DBE","CHOS_O.C","CHOS_H.C")]
tiff(filename = "ss_cor_chp_CHOS.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_ss_CHOS,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_ss_CHONS=ft_merge_soa_nbb_chp_soc_ss[,c("SOC","CHONS","CHONS_AI","CHONS_DBE","CHONS_O.C","CHONS_H.C")]
tiff(filename = "ss_cor_chp_CHONS.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_ss_CHONS,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

###nt=====
ft_merge_soa_nbb_chp_soc_nt=subset(ft_merge_soa_nbb_chp_soc,ft_merge_soa_nbb_chp_soc$Group=="Noto")
ft_merge_soa_nbb_chp_soc_nt_cho=ft_merge_soa_nbb_chp_soc_nt[,c("SOC","CHO","CHO_AI","CHO_DBE","CHO_O.C","CHO_H.C")]

tiff(filename = "nt_cor_chp_cho.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_nt_cho,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_nt_chon=ft_merge_soa_nbb_chp_soc_nt[,c("SOC","CHON","CHON_AI","CHON_DBE","CHON_O.C","CHON_H.C")]

tiff(filename = "nt_cor_chp_chon.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_nt_chon,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_nt_CHOS=ft_merge_soa_nbb_chp_soc_nt[,c("SOC","CHOS","CHOS_AI","CHOS_DBE","CHOS_O.C","CHOS_H.C")]
tiff(filename = "nt_cor_chp_CHOS.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_nt_CHOS,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()

ft_merge_soa_nbb_chp_soc_nt_CHONS=ft_merge_soa_nbb_chp_soc_nt[,c("SOC","CHONS","CHONS_AI","CHONS_DBE","CHONS_O.C","CHONS_H.C")]
tiff(filename = "nt_cor_chp_CHONS.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_merge_soa_nbb_chp_soc_nt_CHONS,
                   method = "spearman",exact=F,label.pos=0.65,blk=2.5)
dev.off()