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

ft_envi=ft_envi
ft_envi$WSOCbb=2.94*ft_envi$Levoglucosan
ft_envi$WSOCnbb=ft_envi$WSOC-ft_envi$WSOCbb

ft_envi_sel=ft_envi[,c("Sample","Group","No","Event","Date","PM2.5","OC","WSOC","WISOC","WSOCbb","WSOCnbb")]
ft_envi_sel$POC=ft_envi_sel$WISOC+ft_envi_sel$WSOCbb
ft_envi_sel$SOC=ft_envi_sel$WSOCnbb

ft_envi_sel$POCp=ft_envi_sel$POC/ft_envi_sel$OC*100
ft_envi_sel$SOCp=ft_envi_sel$SOC/ft_envi_sel$OC*100


FT_envi=fread("Datafile/FRIEND_1st_envi_re.csv")
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

ft_merge$Freq=1

fm_obs=ft_merge[,c("Group","No","Formula","Freq")] %>% 
  dcast(Group~Formula,value.var = "Freq",sum) %>% 
  melt(id.vars=c("Group")) %>% `colnames<-`(c("Group","Formula","cnt"))

ft_merge_filter=ft_merge %>% inner_join(fm_obs)

ft_merge_filter=subset(ft_merge_filter,ft_merge_filter$cnt>4)

ft_merge=ft_merge_filter

ft_comp=as.data.table(aggregate(ft_merge$Freq, by=list(Sample=ft_merge$Sample, Comp=ft_merge$Comp),sum)) %>% 
  `colnames<-`(c("Sample","Comp","cnt"))

ft_comp_tot=as.data.table(aggregate(ft_merge$Freq, by=list(Sample=ft_merge$Sample),sum)) %>% 
  `colnames<-`(c("Sample","tot"))

ft_comp=ft_comp %>% inner_join(ft_comp_tot)

ft_comp$rel=round(ft_comp$cnt/ft_comp$tot*100,1)
ft_comp_sel=ft_comp %>% filter(!Comp%in%c("Remainders"))
ft_comp_sel

table(ft_merge$Comp)
ft_merge$NOS_val=ft_merge$O./(3*ft_merge$N.+4*ft_merge$S.)
ft_merge$NOS=ifelse(ft_merge$O./(3*ft_merge$N.+3*ft_merge$S.)>=1,"NOS","non")
ft_merge$NOS=ifelse(ft_merge$Comp=="CHON"&ft_merge$NOS=="NOS","ON",ft_merge$NOS)
ft_merge$NOS=ifelse(ft_merge$Comp=="CHOS"&ft_merge$NOS=="NOS","OS",ft_merge$NOS)
ft_merge$NOS=ifelse(ft_merge$Comp=="Remainders"|ft_merge$Comp=="CHO",as.character(ft_merge$Comp),ft_merge$NOS)
ft_merge

ft_merge_nos=ft_merge
ft_merge_soa_sel=ft_merge



table(ft_merge_soa_sel$NOS)
ft_comp_soa=as.data.table(aggregate(ft_merge_soa_sel$Freq, by=list(Sample=ft_merge_soa_sel$Sample, Comp=ft_merge_soa_sel$Comp,
                                                                         Type=ft_merge_soa_sel$NOS),sum)) %>% 
  `colnames<-`(c("Sample","Comp","Type","cnt"))

ft_comp_soa_tot=as.data.table(aggregate(ft_merge$Freq, by=list(Sample=ft_merge$Sample),sum)) %>% 
  `colnames<-`(c("Sample","tot"))

ft_comp_soa=ft_comp_soa %>% inner_join(ft_comp_soa_tot)
ft_comp_soa$rel=round(ft_comp_soa$cnt/ft_comp_soa$tot*100,1)
ft_comp_soa_sel=ft_comp_soa %>% filter(!Comp%in%c("Remainders")) %>% droplevels()
ft_comp_soa_sel

ft_comp_soa_sel$Comp2=paste0(ft_comp_soa_sel$Comp," (",ft_comp_soa_sel$Type,")")

FT_envi_sel

ft_envi_sel2=ft_envi_sel[,c("Sample","Group","No","WSOCbb","WSOCnbb")]
FT_envi_sel2=FT_envi_sel[,c("Sample","Group","No","WSOCbb","WSOCnbb")]

FT_envi_sel2
ft_comp_sel_d=dcast(ft_comp_sel,Sample~Comp, value.var = "rel",sum)
ft_comp_soa_sel_d=dcast(ft_comp_soa_sel,Sample~Comp2, value.var = "rel",sum)

ft_comp_sel_d
ft_comp_soa_sel_d

ft_comp_merge=ft_comp_sel_d %>% inner_join(ft_comp_soa_sel_d)

ft_comp_merge=ft_comp_merge %>% inner_join(FT_envi_sel2)

ft_comp_merge=ft_comp_merge[,c("Sample","Group","No","WSOCbb","WSOCnbb","CHO","CHON","CHOS","CHONS",
                               "CHON (ON)","CHOS (OS)","CHONS (NOS)")]
ft_comp_merge

source("Rscript/func_chartcor.R")
library(PerformanceAnalytics)

conlab <- expression(atop(SOC),
                     atop(CHO),
                     atop(CHON),
                     atop(CHOS),
                     atop(CHONS),
                     atop(CHON),
                     atop(CHOS),
                     atop(CHONS))

ft_comp_merge_ul=subset(ft_comp_merge,ft_comp_merge$Group=="Ulaanbaatar")
chart.Correlationx(ft_comp_merge_ul[,-c("Sample","Group","No","WSOCbb","CHO (WSOCbb)","CHON (WSOCbb)","CHOS (WSOCbb)",
                                        "CHONS (WSOCbb)")],
                   method = "spearman",exact=F,labels=conlab,label.pos=0.65,blk=2.5)

tiff(filename = "UL_cor_nbb.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_comp_merge_ul[,-c("Sample","Group","No","WSOCbb","CHO (WSOCbb)","CHON (WSOCbb)","CHOS (WSOCbb)",
                                        "CHONS (WSOCbb)")],
                   method = "spearman",exact=F, labels=conlab,label.pos=0.65,blk=2.5)
dev.off()

ft_comp_merge_bj=subset(ft_comp_merge,ft_comp_merge$Group=="Beijing")
tiff(filename = "Bj_cor_nbb.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_comp_merge_bj[,-c("Sample","Group","No","WSOCbb","CHO (WSOCbb)","CHON (WSOCbb)","CHOS (WSOCbb)",
                                        "CHONS (WSOCbb)")],
                   method = "spearman",exact=F,label.pos=0.65,labels=conlab,blk=2.5)
dev.off()

ft_comp_merge_ss=subset(ft_comp_merge,ft_comp_merge$Group=="Seosan")
tiff(filename = "SS_cor_nbb_re.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_comp_merge_ss[,-c("Sample","Group","No","WSOCbb","CHO (WSOCbb)","CHON (WSOCbb)","CHOS (WSOCbb)",
                                        "CHONS (WSOCbb)")],
                   method = "spearman",exact=F, label.pos=0.65,labels=conlab,blk=2.5)
dev.off()

ft_comp_merge_sul=subset(ft_comp_merge,ft_comp_merge$Group=="Seoul")
tiff(filename = "Sul_cor_nbb.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_comp_merge_sul[,-c("Sample","Group","No","WSOCbb","CHO (WSOCbb)","CHON (WSOCbb)","CHOS (WSOCbb)",
                                         "CHONS (WSOCbb)")],
                   method = "spearman",exact=F, label.pos=0.65,labels=conlab,blk=2.5)
dev.off()

ft_comp_merge_nt=subset(ft_comp_merge,ft_comp_merge$Group=="Noto")
tiff(filename = "nt_cor_nbb.tiff", width = 20,height = 20, units = "cm",res = 300)
chart.Correlationx(ft_comp_merge_nt[,-c("Sample","Group","No","WSOCbb","CHO (WSOCbb)","CHON (WSOCbb)","CHOS (WSOCbb)",
                                        "CHONS (WSOCbb)")],
                   method = "spearman",exact=F, label.pos=0.65,labels=conlab,blk=2.5)
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

dt$rho=ifelse(is.na(dt$rho),"0",dt$rho)
dt$p=ifelse(is.na(dt$p),"1",dt$p)

dt=dt %>% separate("Var2",c("Comp","type"),extra = "drop")

dt$Type=ifelse(is.na(dt$type),"All",
               ifelse(dt$type=="non","Non","SOA"))

dt
fwrite(dt,file = "cor_NBBvsNOS.csv")

dt=fread("Datafile/cor_NBBvsNOS.csv")
dt$mark=ifelse(dt$p<0.001,"***",
               ifelse(dt$p<0.01,"**",
                      ifelse(dt$p<0.05,"*",ifelse(dt$p==1,"","N.S"))))

dt$Group=factor(dt$Group, levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))
dt$Comp=factor(dt$Comp, levels = c("CHO","CHON","CHOS","CHONS"))

dt$rho=as.numeric(dt$rho)
dt

ft_cor_vk_sel2_nos$NOS2=factor(ft_cor_vk_sel2_nos$NOS, levels = c("ON","OS","NOS"),
                               labels = c("ON (O/N ≥ 3)","OS (O/S ≥ 3)","NOS (O ≥ 3N+3S)"))

dt$Type2=factor(dt$Type, levels=c("All","ON","OS","NOS"),
                labels = c("All","ON (O ≥ 3N)","OS (O ≥ 3S)","NOS (O ≥ 3N+3S)"))

ggplot(dt, aes(x=Comp,y=rho,fill=Type2))+
  geom_bar(stat="identity",position = position_dodge(preserve = "single"))+
  geom_text(data=dt, aes(x=Comp, y=ifelse(rho<0,rho-0.05,rho+0.03),label=ifelse(mark=="N.S","",mark)), position = position_dodge(width=0.85),size=7)+
  geom_text(data=dt, aes(x=Comp, y=ifelse(rho<0,rho-0.04,rho+0.04),label=ifelse(mark=="N.S",mark,"")), position = position_dodge(width=0.85),size=5)+
  facet_rep_wrap(.~Group,repeat.tick.labels = T,ncol=5)+
  scale_fill_manual(values = c("grey50","#EFC000FF","#008B45FF","#5F559BFF"))+
  #scale_fill_manual(values = c("grey70","#FD7F20"))+
  scale_y_continuous(name = expression(bold("Correlation coefficient (")*bolditalic("ρ")*bold(")")),
                     limits = c(-1.0,1.0))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.0),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 16,colour = "black",face=2, family = "Arial",angle = 0, hjust = 0.5,vjust = 1.0),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.12,"cm"),
        axis.ticks = element_line(size = 1),
        axis.text.y = element_text(size = 16, colour = "black",face=2, family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        axis.title.y.left = element_text(size = 18, face=2,colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        axis.title.y.right = element_text(size = 18, colour = "black",angle = 90,margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_blank(),
        strip.placement = "outside",
        legend.title = element_text(size = 22, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 16,face=2, hjust=0.5,colour = "black",margin = unit(c(0.0,0.0,0.0,0.0),"cm")),
        legend.key.width = unit(1.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.justification=c(0.5, 0.5),
        legend.position ="bottom")+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("rho_comp_nos_bar"),height = 15, width = 60, units = "cm", dpi = 300)


dt

ft_comp_cor_m

FT_envi_air=FT_envi[,c("Sample","Temp","RH","NH4+","NO3-","SO42-","CO","O3","SO2","NO")]

ft_comp_cor_m=ft_comp_cor_m %>% inner_join(FT_envi_air)
ft_comp_cor_m

ft_comp_cor_mm=melt(ft_comp_cor_m, id.vars = c("Sample","Group","No","WSOCnbb","variable","value"),
                    variable.name = "envi",value.name = "conc")
ft_comp_cor_mm




grp2=unique(ft_comp_cor_mm$Group)
env2=unique(ft_comp_cor_mm$envi)
var2=unique(ft_comp_cor_mm$variable)
data_in=ft_comp_cor_mm


data_in[is.na(data_in)] <- 0

d0=data.table()

for (i in 1:length(grp2)) {
  #  i=1
  temp=subset(data_in,data_in$Group==grp2[i])
  
  for (j in 1:length(var2)) {
    #  j=1
    temp2=subset(temp,temp$variable==var2[j])
    
    for (l in 1:length(env2)) {
      #l=7
      temp3=subset(temp2,temp2$envi==env2[l])
      sc=sum(temp3$conc)
      if(sc==0){
        new=data.table("Group"=grp2[i],"Var1"=env2[l],"Var2"=var2[j],rho=0, p=1)
        
      }
      else{
        ct2=cor.test(temp3$conc,temp3$value, exact = F,method = "spearman")
        new=data.table("Group"=grp2[i],"Var1"=env2[l],"Var2"=var2[j],rho=round(ct2$estimate,2), p=round(ct2$p.value,3))
        
      }
      
      d0=rbind(new,d0)
      
    }
    
  }
  
}

#d00=d0

d0=d0 %>% separate("Var2",c("Comp","type"),extra = "drop")

#fwrite(d0,file = "cor_nosvsenvi.csv")

d0=fread("Datafile/cor_nosvsenvi.csv")

d0$mark=ifelse(d0$p<0.001,"***",
               ifelse(d0$p<0.01,"**",
                      ifelse(d0$p<0.05,"*",ifelse(d0$p==1,"","N.S"))))

d0
d0$Group=factor(d0$Group, levels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))
d0$Comp=factor(d0$Comp, levels = c("CHO","CHON","CHOS","CHONS"))


d0$envi=factor(d0$Var1,levels = c("Sample","Temp","RH","NH4+","NO3-","SO42-","NO","SO2","CO","O3"))
d0=d0 %>% filter(!envi%in%c("Temp","RH"))

table(d0$Type)

d0$Type2=factor(d0$Type, levels=c("All","ON","OS","NOS"),
                labels = c("All","ON (O ≥ 3N)","OS (O ≥ 3S)","NOS (O ≥ 3N+3S)"))

d0

ggplot(d0, aes(x=envi,y=rho,fill=Type2))+
  geom_bar(stat="identity",position = position_dodge(preserve = "single"))+
  geom_text(data=d0, aes(x=envi, y=ifelse(rho<0,rho-0.10,rho+0.03),label=ifelse(mark=="N.S","",mark)), position = position_dodge(width=0.85),size=4)+
  #geom_text(data=d0, aes(x=Var1, y=ifelse(rho<0,rho-0.04,rho+0.04),label=ifelse(mark=="N.S",mark,"")), position = position_dodge(width=0.85),size=3)+
  facet_rep_grid(Comp~Group,repeat.tick.labels = T)+
  scale_fill_manual(values = c("grey50","#EFC000FF","#008B45FF","#5F559BFF"))+
  scale_x_discrete(labels=c(expression(bold(NH[4]^"+")),
                            expression(bold(NO[3]^"-")),
                            expression(bold(SO[4]^"2-")),
                            expression(bold(NO)),
                            expression(bold(SO[2])),
                            expression(bold(CO)),
                            expression(bold(O[3]))
  ))+
  scale_y_continuous(name = expression(bold("Spearman correlation (")*bolditalic("ρ")*bold(")")),
                     limits = c(-0.95,1.05))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0.2,0.2,0.2,0.0),"cm"),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.line.y.right = element_line(size = 1, color = "black"),
        axis.line.x.top =  element_line(colour = "black"),
        axis.text.x = element_text(size = 11,colour = "black",face=2, family = "Arial",angle = 0, hjust = 0.5,vjust = 1.0),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.12,"cm"),
        axis.ticks = element_line(size = 0.5),
        axis.text.y = element_text(size = 11, colour = "black",face=2, family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        axis.title.y.left = element_text(size = 18, face=2,colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        axis.title.y.right = element_text(size = 18, colour = "black",angle = 90,margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, colour = "black",face = 2,margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        #strip.text.y = element_blank(),
        strip.placement = "outside",
        legend.title = element_text(size = 22, colour = "black",margin = unit(c(0.0,0.2,0.0,0.2),"cm")),
        legend.text = element_text(size = 11,face=2, hjust=0.5,colour = "black",margin = unit(c(0.2,0.2,0.2,0.0),"cm")),
        legend.key.width = unit(0.8,"cm"),
        legend.key.height = unit(0.8,"cm"),
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.justification=c(0.5, 0.5),
        legend.position ="bottom")+
  guides(fill=guide_legend(title = NULL,title.hjust = 0.5))+
  ggsave(filename("rho_comp_nos_envi"),height = 24, width = 45, units = "cm", dpi = 300)


