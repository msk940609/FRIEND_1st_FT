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
library(plotly)
#0. Data loading=====
FT_merge=fread("Datafile/FRIENDs_1st_FT.csv")

FT_merge$C=ifelse(FT_merge$C.>0, "C","")
FT_merge$H=ifelse(FT_merge$H.>0, "H","")
FT_merge$O=ifelse(FT_merge$O.>0, "O","")
FT_merge$N=ifelse(FT_merge$N.>0, "N","")
FT_merge$S=ifelse(FT_merge$S.>0, "S","")

FT_merge=FT_merge %>% unite("Comp",c("C","H","O","N","S"),sep = "")
FT_merge$Comp=ifelse(FT_merge$O.==0, "Remainders",FT_merge$Comp)
FT_merge=subset(FT_merge,FT_merge$Comp!="Remainders") %>% droplevels()

##0-1)remove formula less than 3 times during the sampling period==== 
table(FT_merge$Group)
FT_merge_sel=FT_merge %>% filter(Group%in%c("SUL","B","UL"))
FT_merge_sel$Freq=1

fm_sel=melt(FT_merge_sel[,c("Sample","Group","Formula","Freq")],id.vars = c("Sample","Group","Formula")) %>% 
  dcast(Sample+Group~Formula, sum)

fm_sel[,1:23]

fm_sel_b=subset(fm_sel,fm_sel$Group=="B")
fm_sel_s=subset(fm_sel,fm_sel$Group=="SUL")
fm_sel_u=subset(fm_sel,fm_sel$Group=="UL")

tfm_b=as.data.frame(t(fm_sel_b[,-c(1,2)]))
tfm_b$Sum=rowSums(tfm_b)
tfm_b_sel=subset(tfm_b,tfm_b$Sum>3)
tfm_b_sel$Formula=rownames(tfm_b_sel)

tfm_s=as.data.frame(t(fm_sel_s[,-c(1,2)]))
tfm_s$Sum=rowSums(tfm_s)
tfm_s_sel=subset(tfm_s,tfm_s$Sum>3)
tfm_s_sel$Formula=rownames(tfm_s_sel)


tfm_u=as.data.frame(t(fm_sel_u[,-c(1,2)]))
tfm_u$Sum=rowSums(tfm_u)
tfm_u_sel=subset(tfm_u,tfm_u$Sum>3)
tfm_u_sel$Formula=rownames(tfm_u_sel)

tfm_u_sel
#fmlist=rbind.data.frame(tfm_b_sel[,c("Formula","Sum")],tfm_s_sel[,c("Formula","Sum")],
#                        tfm_u_sel[,c("Formula","Sum")])

fmlist=rbind.data.frame(tfm_b_sel[,c("Formula","Sum")],tfm_s_sel[,c("Formula","Sum")])

fmlist$type="cor"
fmlist=unique(fmlist[,c("Formula","type")])
fmlist

fm_cor=FT_merge_sel %>% left_join(fmlist)
fm_cor=fm_cor %>% filter(Group%in%c("SUL","B"))
fm_cor=subset(fm_cor,fm_cor$type=="cor")

##0-2)build formula table to correlation test====
fm_cor
ft_cor=melt(fm_cor[,c("Sample","Formula","Bromo.Inty")], id.vars = c("Sample","Formula")) %>% 
  dcast(Sample~Formula,sum)

ft_cor[,1:4]
#ft_cor=ft_cor %>% separate(Sample,c("Group"),sep = "_",extra = "drop")

envi_cor=fread("Datafile/envi_1st_S&B.csv")
envi_cor$Group=ifelse(envi_cor$Group=="Beijing","B","SUL")
envi_cor$Sample=paste(envi_cor$Group,envi_cor$No,sep = "_")
envi_cor

envi_cor_sel=envi_cor[,c("Sample","SO42-","NO3-")]

envi_cor_sel=envi_cor_sel %>% inner_join(ft_cor)
envi_cor_sel=envi_cor_sel %>% separate(Sample,c("Group"),sep = "_",extra = "drop")
envi_cor_sel[,1:7]


#1.Correlation between Envi vs Formula====
fwrite(envi_cor_sel,"Datafile/cor_S&B.csv")


grp=unique(envi_cor_sel$Group)
grp

cor_temp=envi_cor_sel
df=data.frame()
for (i in 1:length(grp)) {
  #i=1
  temp=as.data.frame(subset(cor_temp,cor_temp$Group==grp[i]))
  for (j in 2:3) {
    #j=2
    for (k in 4:dim(temp)[2]) {
    #k=4  
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
  
}

df


ft_cor=fread("Datafile/cortest.csv")
ft_cor_so4=subset(ft_cor,ft_cor$Envi=="SO42-")
ft_cor_no3=subset(ft_cor,ft_cor$Envi=="NO3-")

fm_cor_no3=fm_cor %>% inner_join(ft_cor_no3)
fm_cor_no3

fm_cor_so4=fm_cor %>% inner_join(ft_cor_so4)
fm_cor_so4

fm_cor_fin=rbind(fm_cor_no3,fm_cor_so4)
fm_cor_fin
fm_cor_vk=subset(fm_cor_fin,fm_cor_fin$p<0.05)
#fm_cor_vk=subset(fm_cor_vk,fm_cor_vk$Comp!="CHO") %>% droplevels()

pg=ggplot(fm_cor_vk, aes(x=O.C, y=H.C, col=rho,Formula=Formula))+
  geom_point(alpha=0.2, size=3)+
  facet_grid(Envi~Group)+
  scale_color_gradient2()+
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
        legend.text = element_text(size = 14, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank()
  )

pg+ggsave(filename("B&S_cor"), height = 44, width = 40, units = "cm", dpi = 300)


ggplotly(pg,tooltip = c("x", "y","Formula")) %>% 
  layout(legend = list(orientation = "h", x = 0.3, y = -0.2))

fm_cor_vk_neg=subset(fm_cor_vk,fm_cor_vk$rho<0)

fm_cor_vk_neg_b=subset(fm_cor_vk_neg,fm_cor_vk_neg$Group=="B")
fm_cor_vk_neg_s=subset(fm_cor_vk_neg,fm_cor_vk_neg$Group=="SUL")

table(fm_cor_vk_neg_b$Comp)
table(fm_cor_vk_neg_s$Comp)

