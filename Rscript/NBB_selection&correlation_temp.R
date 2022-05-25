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

ft_envi=fread("Datafile/FRIEND_1st_envi_re.csv")

FT_envi=ft_envi
FT_envi$WSOCbb=2.94*FT_envi$Levoglucosan
FT_envi$WSOCnbb=FT_envi$WSOC-FT_envi$WSOCbb

FT_envi_sel=FT_envi[,c("Sample","Group","No","Event","Date","PM2.5","OC","WSOC","WISOC","WSOCbb","WSOCnbb")]
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
fm_obs=ft_merge[,c("Group","No","Formula","Freq")] %>% 
  dcast(Group~Formula,value.var = "Freq",sum) %>% 
  melt(id.vars=c("Group")) %>% `colnames<-`(c("Group","Formula","cnt"))

ft_merge_filter=ft_merge %>% inner_join(fm_obs)

ft_merge_filter=subset(ft_merge_filter,ft_merge_filter$cnt>4)

ft_merge=ft_merge_filter

ft_merge$Freq=1

ft_comp=as.data.table(aggregate(ft_merge$Freq, by=list(Sample=ft_merge$Sample, Comp=ft_merge$Comp),sum)) %>% 
  `colnames<-`(c("Sample","Comp","cnt"))

ft_comp_tot=as.data.table(aggregate(ft_merge$Freq, by=list(Sample=ft_merge$Sample),sum)) %>% 
  `colnames<-`(c("Sample","tot"))

ft_comp=ft_comp %>% inner_join(ft_comp_tot)

ft_comp$rel=round(ft_comp$cnt/ft_comp$tot*100,1)
ft_comp_sel=ft_comp %>% filter(!Comp%in%c("Remainders"))
ft_comp_sel

ft_comp_sel_d=dcast(ft_comp_sel,Sample~Comp, value.var = "cnt",sum)

ft_merge_soa=ft_merge %>% left_join(ft_cor_pos[,c("Group","Formula","type")])
table(ft_merge_soa$type)

ft_merge_soa_sel=ft_merge_soa %>% filter(type%in%c("WSOCbb","WSOCnbb"))
ft_comp_soa=as.data.table(aggregate(ft_merge_soa_sel$Bromo.Inty, by=list(Sample=ft_merge_soa_sel$Sample, Comp=ft_merge_soa_sel$Comp,
                                                                         Type=ft_merge_soa_sel$type),sum)) %>% 
  `colnames<-`(c("Sample","Comp","Type","cnt"))
ft_comp_soa_tot=as.data.table(aggregate(ft_merge_soa$Bromo.Inty, by=list(Sample=ft_merge_soa$Sample),sum)) %>% 
  `colnames<-`(c("Sample","tot"))
ft_comp_soa=ft_comp_soa %>% inner_join(ft_comp_soa_tot)
ft_comp_soa$rel=round(ft_comp_soa$cnt/ft_comp_soa$tot*100,1)
ft_comp_soa_sel=ft_comp_soa %>% filter(!Comp%in%c("Remainders")) %>% droplevels()
ft_comp_soa_sel

ft_comp_soa_sel$Comp2=paste0(ft_comp_soa_sel$Comp," (",ft_comp_soa_sel$Type,")")

FT_envi_sel
FT_envi_sel2=FT_envi_sel[,c("Sample","Group","No","WSOCbb","WSOCnbb")]

FT_envi_sel2

ft_comp_sel_d=dcast(ft_comp_sel,Sample~Comp, value.var = "rel",sum)
ft_comp_soa_sel_d=dcast(ft_comp_soa_sel,Sample~Comp2, value.var = "rel",sum)

ft_comp_sel_d
ft_comp_soa_sel_d

ft_comp_merge=ft_comp_sel_d %>% inner_join(ft_comp_soa_sel_d)

ft_comp_merge=ft_comp_merge %>% inner_join(FT_envi_sel2)

ft_comp_merge=ft_comp_merge[,c("Sample","Group","No","WSOCbb","WSOCnbb","CHO","CHON","CHOS","CHONS",
                               "CHO (WSOCbb)","CHON (WSOCbb)","CHOS (WSOCbb)","CHONS (WSOCbb)",
                               "CHO (WSOCnbb)","CHON (WSOCnbb)","CHOS (WSOCnbb)","CHONS (WSOCnbb)")]

ft_comp_merge


source("Rscript/func_chartcor.R")
library(PerformanceAnalytics)

conlab <- expression(atop(SOC),
                     atop(CHO),
                     atop(CHON),
                     atop(CHOS),
                     atop(CHONS))


ft_comp_merge_ss=subset(ft_comp_merge,ft_comp_merge$Group=="Seosan")
tiff(filename = "ss_cor_all_re.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_comp_merge_ss[,-c("Sample","Group","No","WSOCbb","CHO (WSOCbb)","CHON (WSOCbb)","CHOS (WSOCbb)",
                                        "CHONS (WSOCbb)","CHO (WSOCnbb)","CHON (WSOCnbb)","CHOS (WSOCnbb)",
                                        "CHONS (WSOCnbb)")],
                   method = "spearman",exact=F, labels=conlab,label.pos=0.65,blk=2.5)
dev.off()

tiff(filename = "ss_cor_SOC_re.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_comp_merge_ss[,-c("Sample","Group","No","WSOCbb","CHO (WSOCbb)","CHON (WSOCbb)","CHOS (WSOCbb)",
                                        "CHONS (WSOCbb)","CHO","CHON","CHOS",
                                        "CHONS")],
                   method = "spearman",exact=F, labels=conlab,label.pos=0.65,blk=2.5)
dev.off()

ft_comp_cor=ft_comp_merge[,-c("WSOCbb")]
ft_comp_cor_m=melt(ft_comp_cor,id.vars = c("Sample","Group","No","WSOCnbb"))
ft_comp_cor_m


grp=unique(ft_comp_cor_m$Group)
var=unique(ft_comp_cor_m$variable)
data_in=ft_comp_cor_m

dt=data.table()
for (i in 1:length(grp)) {
  #  i=1
  temp=subset(data_in,data_in$Group==grp[i])
  
  for (j in 1:length(var)) {
    #    j=1
    temp2=subset(temp,temp$variable==var[j])
    
    ct=cor.test(temp2$WSOCnbb,temp2$value, exact = F,method = "spearman")
    
    new=data.table("Group"=grp[i],"Var1"="WSOCnbb","Var2"=var[j],rho=round(ct$estimate,2), p=round(ct$p.value,3))
    
    dt=rbind(new,dt)
    
  }
  
}

dt0=dt
dt=dt %>% separate("Var2",c("Comp","type"),extra = "drop")

dt$Type=ifelse(is.na(dt$type),"All","SOC")
dt$mark=ifelse(dt$p<0.001,"***",
               ifelse(dt$p<0.01,"**",
                      ifelse(dt$p<0.05,"*","N.S")))

fwrite(dt,file = "cor_NBBvsComp_seosan.csv")