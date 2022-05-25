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

library(ggrepel)

#0.Load datat====
FT_merge=fread("Datafile/FRIENDs_1st_FT.csv")
FT_merge$Sample=paste(FT_merge$Group,FT_merge$No,sep = "_")

FT_merge$C=ifelse(FT_merge$C.>0, "C","")
FT_merge$H=ifelse(FT_merge$H.>0, "H","")
FT_merge$O=ifelse(FT_merge$O.>0, "O","")
FT_merge$N=ifelse(FT_merge$N.>0, "N","")
FT_merge$S=ifelse(FT_merge$S.>0, "S","")

FT_merge=FT_merge %>% unite("Comp",c("C","H","O","N","S"),sep = "")
FT_merge$Comp=ifelse(FT_merge$O.==0, "Remainders",FT_merge$Comp)

FT_merge=subset(FT_merge,FT_merge$Comp!="Remainders") %>% droplevels()
FT_merge

ft_inty_all=as.data.table(aggregate(FT_merge$Bromo.Inty, by=list(`Sample`=FT_merge$Sample), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Tot"))
ft_inty=as.data.table(aggregate(FT_merge$Bromo.Inty, by=list(`Sample`=FT_merge$Sample, Comp=FT_merge$Comp), FUN=sum))
ft_inty
FT_merge$Freq=1

#ft_inty_all=as.data.table(aggregate(FT_merge$Freq, by=list(`Sample`=FT_merge$Sample), FUN=sum)) %>% 
#  `colnames<-`(c("Sample","Tot"))
#ft_inty=as.data.table(aggregate(FT_merge$Freq, by=list(`Sample`=FT_merge$Sample, Comp=FT_merge$Comp), FUN=sum))
#ft_inty


ft_inty=ft_inty %>% inner_join(ft_inty_all,by = "Sample")
ft_inty$rel=round(ft_inty$x/ft_inty$Tot*100,1)

ft_inty=ft_inty %>% separate(Sample, c("Group","No"),sep = "_")
ft_inty

ft_inty$Comp=factor(ft_inty$Comp,levels = c("CHO","CHON","CHOS","CHONS"))

ft_inty$Group=factor(ft_inty$Group,levels = c("SUL","SS","B","NT","UL"),
                     labels = c("Seoul","Seosan","Beijing","Noto","Ulaanbaatar"))
ft_inty

ft_inty_m=dcast(ft_inty,Group+No~Comp, sum, value.var = "rel")
ft_inty_m$Sample=paste(ft_inty_m$Group,ft_inty_m$No,sep = "_")

FT_bb=fread("Datafile/WOSCbb_S&B.csv")

ft_inty_m_sel=ft_inty_m %>% filter(Group%in%c("Beijing","Seoul"))
ft_inty_m_sel

dim(FT_bb)
ft_WSOC_sel=FT_bb[,c("Sample","WSOCbb","WSOCnbb","POA","SOA")] %>% inner_join(ft_inty_m_sel)
ft_WSOC_sel

ft_WSOC_sel_m=melt(ft_WSOC_sel, id.vars = c("Sample","Group","No","WSOCbb","WSOCnbb","POA","SOA"),variable.name = "Comp") %>% 
  melt(id.vars=c("Sample","Group","No","Comp","value"),variable.name="Envi",value.name="Conc")

ft_WSOC_sel_m

ft_WSOC_sel_m$No=as.numeric(ft_WSOC_sel_m$No)
ft_WSOC_sel_m$mark="no"
ft_WSOC_sel_m$mark=ifelse(ft_WSOC_sel_m$No==5,"12/19",ft_WSOC_sel_m$mark)
ft_WSOC_sel_m$mark=ifelse(ft_WSOC_sel_m$No==6,"12/20",ft_WSOC_sel_m$mark)
ft_WSOC_sel_m$mark=ifelse(ft_WSOC_sel_m$No==7,"12/21",ft_WSOC_sel_m$mark)
ft_WSOC_sel_m$mark=ifelse(ft_WSOC_sel_m$No==8,"12/22",ft_WSOC_sel_m$mark)
ft_WSOC_sel_m$mark=ifelse(ft_WSOC_sel_m$No==9,"12/23",ft_WSOC_sel_m$mark)

ft_WSOC_sel_m
table(ft_WSOC_sel_m$mark)

ft_WSOC_sel_m$Envilab=ft_WSOC_sel_m$Envi
ft_WSOC_sel_m$Envilab=factor(ft_WSOC_sel_m$Envilab, levels = c("WSOCbb","WSOCnbb","POA","SOA"),
                                labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("POA"~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("SOA"~"("*"\u03bcg/"*m^"3"*")"))))

ft_WSOC_sel_m$Complab=ft_WSOC_sel_m$Comp
ft_WSOC_sel_m$Complab=factor(ft_WSOC_sel_m$Complab, levels = c("CHO","CHON","CHOS","CHONS"),
                                labels = c(expression(bold("CHO")),expression(bold("CHON")),expression(bold("CHOS")),
                                           expression(bold("CHONS"))))
ft_WSOC_sel_m

##Seoul====
ft_WSOC_sel_m_s=subset(ft_WSOC_sel_m,ft_WSOC_sel_m$Group=="Seoul")

val_max_s1=as.data.table(aggregate(ft_WSOC_sel_m_s$Conc, 
                                   by=list(`Comp`=ft_WSOC_sel_m_s$Comp,Envi=ft_WSOC_sel_m_s$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))
val_max_s1$Envilab=val_max_s1$Envi

val_max_s1$Envilab=val_max_s1$Envi
val_max_s1$Envilab=factor(val_max_s1$Envilab, levels = c("WSOCbb","WSOCnbb","POA","SOA"),
                             labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                        expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")")),
                                        expression(bold("POA"~"("*"\u03bcg/"*m^"3"*")")),
                                        expression(bold("SOA"~"("*"\u03bcg/"*m^"3"*")"))))

val_max_s1$Complab=val_max_s1$Comp
val_max_s1$Complab=factor(val_max_s1$Complab, levels = c("CHO","CHON","CHOS","CHONS"),
                             labels = c(expression(bold("CHO")),expression(bold("CHON")),expression(bold("CHOS")),
                                        expression(bold("CHONS"))))
val_max_s1

ggplot()+
  geom_smooth(data=ft_WSOC_sel_m_s,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=ft_WSOC_sel_m_s,aes(x=value,y=Conc),col="royalblue2", size=3)+
  geom_blank(data=val_max_s1, aes(x=20,y=Conc*1.2))+
  geom_point(data=subset(ft_WSOC_sel_m_s,ft_WSOC_sel_m_s$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_WSOC_sel_m_s,ft_WSOC_sel_m_s$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of chemical composition (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsWSOC&SOA_sul"),height = 30, width = 30, units = "cm", dpi = 300)


##Beijing====
ft_WSOC_sel_m_b=subset(ft_WSOC_sel_m,ft_WSOC_sel_m$Group=="Beijing")

val_max_b1=as.data.table(aggregate(ft_WSOC_sel_m_b$Conc, 
                                   by=list(`Comp`=ft_WSOC_sel_m_b$Comp,Envi=ft_WSOC_sel_m_b$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))
val_max_b1$Envilab=val_max_b1$Envi
val_max_b1$Envilab=factor(val_max_b1$Envilab, levels = c("WSOCbb","WSOCnbb","POA","SOA"),
                          labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("POA"~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("SOA"~"("*"\u03bcg/"*m^"3"*")"))))


val_max_b1$Complab=val_max_b1$Comp
val_max_b1$Complab=factor(val_max_b1$Complab, levels = c("CHO","CHON","CHOS","CHONS"),
                          labels = c(expression(bold("CHO")),expression(bold("CHON")),expression(bold("CHOS")),
                                     expression(bold("CHONS"))))
val_max_b1

ggplot()+
  geom_smooth(data=ft_WSOC_sel_m_b,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=ft_WSOC_sel_m_b,aes(x=value,y=Conc),col="orangered3", size=3)+
  geom_blank(data=val_max_b1, aes(x=30.5,y=Conc*1.2))+
  geom_point(data=subset(ft_WSOC_sel_m_b,ft_WSOC_sel_m_b$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(ft_WSOC_sel_m_b,ft_WSOC_sel_m_b$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "Proportion of chemical composition (%)")+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 18, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.85,0.10),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  ggsave(filename("compvsWSOCSOA_b"),height = 30, width = 30, units = "cm", dpi = 300)

##FT_envi correlation======
ft_WSOC_sel_m

grp=unique(ft_WSOC_sel_m$Group)
grp
comp=unique(ft_WSOC_sel_m$Comp)
envi=unique(ft_WSOC_sel_m$Envi)


dt=data.table()
for (i in 1:length(grp)) {
  temp=subset(ft_WSOC_sel_m,ft_WSOC_sel_m$Group==grp[i])
  
  for (k in 1:length(comp)) {
    tmp=subset(temp,temp$Comp==comp[k])
    
    for (l in 1:length(envi)) {
      tmp2=subset(tmp,tmp$Envi==envi[l])
      
      if(envi[l]=="O3"&grp[i]=="Ulaanbaatar") {
        tmp2$Conc=0
        new_row=data.table("Group"=grp[i],"Comp"=comp[k],"R2"=0,"Envi"=envi[l],"Cor"=0, "P"=1)
        
      } else {
        model=lm(Conc~value,tmp2)
        a=summary(model)
        cor=cor.test(tmp2$value,tmp2$Conc, method = "pearson")
        #cor=cor.test(tmp2$value,tmp2$Conc, method = "spearman", exact = F)
        #    cor
        
        new_row=data.table("Group"=grp[i],"Comp"=comp[k],"R2"=round(a$r.squared,8),"Envi"=envi[l],"Cor"=cor$estimate, "P"=round(cor$p.value,3))
        
      }
      
      dt=rbind(dt,new_row)
    }
  }
}

dt

dt$sig=ifelse(dt$P<0.05,ifelse(dt$P<0.01,ifelse(dt$P<0.001,"***","**"),"*"),"")


#fwrite(dt,"Datafile/211103FT&Envicor.csv")

dt$Envilab=dt$Envi
dt$Envilab=factor(dt$Envilab,levels = c("PM2.5","NH4","NO3","SO4","O3","CO","NO","SO2"),
                  labels = c(expression(bold("PM"["2.5"])),
                             expression(bold("NH"["4"]^{"+"})),
                             expression(bold("NO"["3"]^{"-"})),
                             expression(bold("SO"["4"]^{"2-"})),
                             expression(bold("O"["3"])),
                             expression(bold("CO")),
                             expression(bold("NO")),
                             expression(bold("SO"["2"]))
                  )
)

dt$Group=factor(dt$Group, levels = c("Ulaanbaatar","Beijing","Seoul"))

p=ggplot(dt, aes(x=Envilab,y=Cor,fill=Group))+
  geom_bar(stat="identity", position = "dodge", alpha=0.85)+
  scale_x_discrete("Environmental variables", labels=parse(text=levels(dt$Envilab)))+
  scale_y_continuous(name = "Spearman correlation coefficient (rho)                         Spearman correlation coefficient (rho)")+
  scale_fill_manual(values = c("seagreen4","orangered3","royalblue2"))+
  facet_wrap(.~Comp, scales = "free")+
  geom_text(aes(Envilab, y=ifelse(Cor>0,Cor*1.05,Cor*1.08), label=sig),size=10, position = position_dodge(width = 0.9))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.6,0,0.2),"cm")),
        strip.text = element_text(size = 20, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 26, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.key.width = unit(1.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.07,0.93),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )

p+ggsave(filename("compvsEnvi"),height = 40, width = 60, units = "cm", dpi = 300)

dt
dt$Envi=factor(dt$Envi,levels = rev(c("PM2.5","NH4","NO3","SO4","O3","CO","NO","SO2")))


matlab.like2(10)

p2=ggplot(dt)+
  geom_tile(data=subset(dt,dt$P<0.05), aes(x=Comp,y=Envi,fill=Cor),alpha=0.8)+
  facet_wrap(Group~.)+
  #scale_fill_gradient(low = "")
  #scale_fill_gradientn(colors=blue2yellow(2))+
  scale_fill_gradientn(colors=rev(matlab.like(20)))+
  scale_x_discrete(name="")+
  scale_y_discrete(name="",labels=rev(c(expression(bold("PM"["2.5"])),
                                        expression(bold("NH"["4"]^{"+"})),
                                        expression(bold("NO"["3"]^{"-"})),
                                        expression(bold("SO"["4"]^{"2-"})),
                                        expression(bold("O"["3"])),
                                        expression(bold("CO")),
                                        expression(bold("NO")),
                                        expression(bold("SO"["2"])))))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.6,0.4,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.6,0,0.2),"cm")),
        strip.text = element_text(size = 20, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 26, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.key.height = unit(1.0,"cm"),
        legend.key.width = unit(1.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )

p2+ggsave(filename("heat_map_compvsEnvi"),height = 30, width = 45, units = "cm", dpi = 300)

install.packages("colorRamps")
library(colorRamps)
matlab.like2(10)

YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
colorRampPalette(YlOrBr, space = "Lab",bias = 0.5)
#Event vs Non-event=====
##1.Chemical composition====
FT_merge

FT_envi=fread("Datafile/FRIENDs_envi.csv")

ft_inty_all=ft_inty
ft_inty_all$Sample=paste(ft_inty_all$Group,ft_inty_all$No,sep = "_")
ft_inty_all=ft_inty_all %>% inner_join(FT_envi[,c("Sample","Date","Event")]) 


ft_inty_ev=ft_inty_all %>% filter(Event%in%c("Event","Non-event"))
ft_inty_ev

ft_inty_ev$Group=factor(ft_inty_ev$Group,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))

ft_inty_ev$Complab=paste(ft_inty_ev$Comp,"(%)")
ft_inty_ev$Complab=factor(ft_inty_ev$Complab,levels = c("CHO (%)","CHON (%)","CHOS (%)","CHONS (%)"))

val_max_ev=as.data.table(aggregate(ft_inty_ev$rel, 
                                   by=list(`Comp`=ft_inty_ev$Comp,Ev=ft_inty_ev$Event), FUN=max)) %>% 
  `colnames<-`(c("Comp","Event","rel"))
val_max_ev
val_max_ev$Complab=paste(val_max_ev$Comp,"(%)")
val_max_ev$Complab=factor(val_max_ev$Complab,levels = c("CHO (%)","CHON (%)","CHOS (%)","CHONS (%)"))

p1=ggplot(data=ft_inty_ev, aes(x=Group, y=rel, fill=Event))+
  #geom_violin(position = "dodge")+
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.8), width=0.4)+
  geom_boxplot(position = position_dodge(width = 0.80), outlier.colour = NA)+
  facet_wrap(Complab~., strip.position = "left", scales = "free", ncol=2)+
  geom_blank(data=val_max_ev, aes(x=2.5,y=rel*1.12))+
  scale_x_discrete("")+
  scale_fill_manual(values = c("grey75","white"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.0,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.4,0.2,0.4),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 24, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.08,0.94),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )

p1+ggsave(filename("Eventcompare"),height = 30, width = 40, units = "cm", dpi = 300)

##Compare distribution (Event vsNon-event)=====
grp=as.vector(unique(ft_inty_ev$Group))
comp=as.vector(unique(ft_inty_ev$Comp))

data=ft_inty_ev
#data=subset(ft_inty_ev,ft_inty_ev$Sample!="Ulaanbaatar_8") %>% droplevels()

dt=data.table()
for (i in 1:length(comp)) {
  #  i=3
  temp=subset(data,data$Comp==comp[i])
  for (j in 1:length(grp)) {
    #    j=5
    temp2=subset(temp,temp$Group==grp[j])
    
    t_eve=subset(temp2,temp2$Event=="Event")
    t_nev=subset(temp2,temp2$Event=="Non-event")
    
    s1=shapiro.test(t_eve$rel)
    s2=shapiro.test(t_nev$rel)
    
    evtmp=levene.test(temp2$rel,temp2$Event)
    
    tt=t.test(t_eve$rel,t_nev$rel, paired = F,var.equal = F)
    tt2=wilcox.test(t_eve$rel,t_nev$rel, exact = F)
    
    new=data.table("Comp"=comp[i],"Group"=grp[j],"norm_eve"=round(s1$p.value,3),"norm_nev"=round(s2$p.value,3),
                   "Eq-var"=round(evtmp$p,3),"T-test"=round(tt$p.value,3),"non_ttest"=round(tt2$p.value,3))
    
    dt=rbind(dt,new)
    
  }
}
dt

#dt_non_equal=dt
dt
##4-2-2)Calculate label position=====
uf=data.table()
for (i in 1:length(comp)) {
  #i=4
  temp=subset(data,data$Comp==comp[i])
  for (j in 1:length(grp)) {
    #j=5
    temp2=subset(temp,temp$Group==grp[j])
    temp2
    
    t_eve=subset(temp2,temp2$Event=="Event")
    t_nev=subset(temp2,temp2$Event=="Non-event")
    
    iqr.df <- ddply(temp2, "Event", function (x) ((1.5*IQR(x[,"rel"])+quantile(x[,"rel"],0.75)))) %>% 
      `colnames<-`(c("Event","iqr"))
    
    t3=temp2 %>% inner_join(iqr.df)
    
    t4=subset(t3,t3$rel<t3$iqr)
    t4
    
    uf.ps <- ddply(t4, "Event", function (x) (max(fivenum(x[,"rel"]))))
    new=data.table("Comp"=comp[i],"Group"=grp[j],"uf_eve"=uf.ps[1,2],"uf_nev"=uf.ps[2,2] )
    uf=rbind(uf,new)
    
  }
}

uf$uf_pos=ifelse(uf$uf_eve>=uf$uf_nev,uf$uf_eve,uf$uf_nev)
uf$Group=factor(uf$Group,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))


uf$id=paste(uf$Comp,uf$Group,sep = "_")
dt$id=paste(dt$Comp,dt$Group,sep = "_")
#dt$Noneq_t=dt_non_equal$`T-test`
#dt$p=ifelse(dt$`T-test`<dt$Noneq_t,dt$`T-test`,dt$Noneq_t)
dt
uf=uf %>% inner_join(dt[,c("id","non_ttest")])
uf$sig=ifelse(uf$non_ttest<0.001,"***",
              ifelse(uf$non_ttest<0.01,"**",
                     ifelse(uf$non_ttest<0.05,"*","N.S")))

#uf$Group=factor(uf$Group,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))
uf$Comp=factor(uf$Comp,levels = c("CHO","CHON","CHOS","CHONS"))

uf$xs=ifelse(uf$Group=="Seoul",3,
             ifelse(uf$Group=="Seosan",4,
                    ifelse(uf$Group=="Beijing",2,
                           ifelse(uf$Group=="Noto",5,1))))

uf$Complab=paste(uf$Comp,"(%)")
uf$Complab=factor(uf$Complab,levels = c("CHO (%)","CHON (%)","CHOS (%)","CHONS (%)"))

ps=ggplot()+
  #geom_violin(position = "dodge")+
  stat_boxplot(data=ft_inty_ev, aes(x=Group, y=rel, fill=Event),geom = "errorbar",position = position_dodge(width = 0.8), width=0.4)+
  #stat_boxplot(data=subset(ft_inty_ev,ft_inty_ev$Sample!="Ulaanbaatar_8"), aes(x=Group, y=rel, fill=Event),geom = "errorbar",position = position_dodge(width = 0.8), width=0.4)+
  geom_boxplot(data=ft_inty_ev, aes(x=Group, y=rel, fill=Event),position = position_dodge(width = 0.80), outlier.colour = NA)+
  #geom_boxplot(data=subset(ft_inty_ev,ft_inty_ev$Sample!="Ulaanbaatar_8"), aes(x=Group, y=rel, fill=Event),position = position_dodge(width = 0.80), outlier.colour = NA)+
  facet_wrap(Complab~., strip.position = "left", scales = "free", ncol=2)+
  geom_blank(data=val_max_ev, aes(x=2.5,y=rel*1.12))+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), lty=2)+
  geom_segment(data=uf,aes(x=xs-0.2,xend=xs+0.2,y=uf_pos*1.03, yend=uf_pos*1.03))+
  #geom_segment(data=uf,aes(x=xs-0.2,xend=xs-0.2,y=uf_eve*1.01, yend=uf_pos*1.05))+
  #geom_segment(data=uf,aes(x=xs+0.2,xend=xs+0.2,y=uf_nev*1.01, yend=uf_pos*1.05))+
  #geom_text(data =comp_max, aes(x=1, y=x, label=Comp), col="white", size=0)+
  geom_text(data =subset(uf,uf$non_ttest<0.05), aes(x=Group, y=uf_pos*1.05, label=sig), col="black", size=6)+
  scale_x_discrete("")+
  scale_fill_manual(values = c("grey75","white"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.0,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.4,0.2,0.4),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 24, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.08,0.94),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_rect(fill="white")
  )

ps+ggsave(filename("Eventcompare_sig"),height = 30, width = 40, units = "cm", dpi = 300)

##molecular comparision======
FT_mole=fread("Datafile/FRIEND_molecule.csv")

FT_mole$Sample=paste(FT_mole$Group,FT_mole$No,sep = "_")

FT_mole_all=melt(FT_mole,id.vars = c("Sample","Group","No"))

FT_mole_all=FT_mole_all %>% inner_join(FT_envi[,c("Sample","Date","Event")]) 

FT_mole_all

FT_mole_ev=FT_mole_all %>% filter(Event%in%c("Event","Non-event"))
FT_mole_ev

FT_mole_ev$Group=factor(FT_mole_ev$Group,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))

FT_mole_ev$varlab=paste(FT_mole_ev$variable,"(%)")

unique(FT_mole_ev$varlab)


FT_mole_ev$varlab=factor(FT_mole_ev$varlab,levels = c("Lignins (%)","Condensed aromatics (%)","Tannins (%)", "Unsaturated hydrocarbons (%)",
                                                      "Carbohydrates (%)","Proteins (%)", "Lipids (%)"))

mol_max_ev=as.data.table(aggregate(FT_mole_ev$value, 
                                   by=list(`variable`=FT_mole_ev$variable,Ev=FT_mole_ev$Event), FUN=max)) %>% 
  `colnames<-`(c("variable","Event","value"))

mol_max_ev

mol_max_ev$varlab=paste(mol_max_ev$variable,"(%)")
mol_max_ev$varlab=factor(mol_max_ev$varlab,levels = c("Lignins (%)","Condensed aromatics (%)","Tannins (%)", "Unsaturated hydrocarbons (%)",
                                                      "Carbohydrates (%)","Proteins (%)", "Lipids (%)"))

p1=ggplot(data=FT_mole_ev, aes(x=Group, y=value, fill=Event))+
  #geom_violin(position = "dodge")+
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.8), width=0.4)+
  geom_boxplot(position = position_dodge(width = 0.80), outlier.colour = NA)+
  facet_wrap(varlab~., strip.position = "left", scales = "free", ncol=4)+
  geom_blank(data=mol_max_ev, aes(x=2.5,y=value*1.12))+
  scale_x_discrete("")+
  scale_fill_manual(values = c("grey75","white"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.0,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.4,0.2,0.4),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 24, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.82,0.94),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )

p1+ggsave(filename("Event_mole_compare"),height = 35, width = 75, units = "cm", dpi = 300)

##Compare distribution (Event vsNon-event)=====
grp=as.vector(unique(FT_mole_ev$Group))
comp=as.vector(unique(FT_mole_ev$variable))

data=FT_mole_ev
#data=subset(FT_mole_ev,FT_mole_ev$Sample!="Ulaanbaatar_8") %>% droplevels()

dt=data.table()
for (i in 1:length(comp)) {
  # i=1
  temp=subset(data,data$variable==comp[i])
  for (j in 1:length(grp)) {
    #    j=1
    temp2=subset(temp,temp$Group==grp[j])
    
    t_eve=subset(temp2,temp2$Event=="Event")
    t_nev=subset(temp2,temp2$Event=="Non-event")
    
    s1=shapiro.test(t_eve$value)
    s2=shapiro.test(t_nev$value)
    
    evtmp=levene.test(temp2$value,temp2$Event)
    
    tt=t.test(t_eve$value,t_nev$value, paired = F,var.equal = F)
    tt2=wilcox.test(t_eve$value,t_nev$value, exact = F)
    
    new=data.table("Comp"=comp[i],"Group"=grp[j],"norm_eve"=round(s1$p.value,3),"norm_nev"=round(s2$p.value,3),
                   "Eq-var"=round(evtmp$p,3),"T-test"=round(tt$p.value,3),"non_ttest"=round(tt2$p.value,3))
    
    dt=rbind(dt,new)
    
  }
}
dt

#dt_non_equal=dt
dt
##4-2-2)Calculate label position=====
ufm=data.table()
for (i in 1:length(comp)) {
  #i=1
  temp=subset(data,data$variable==comp[i])
  for (j in 1:length(grp)) {
    #j=1
    temp2=subset(temp,temp$Group==grp[j])
    temp2
    
    t_eve=subset(temp2,temp2$Event=="Event")
    t_nev=subset(temp2,temp2$Event=="Non-event")
    
    iqr.df <- ddply(temp2, "Event", function (x) ((1.5*IQR(x[,"value"])+quantile(x[,"value"],0.75)))) %>% 
      `colnames<-`(c("Event","iqr"))
    
    t3=temp2 %>% inner_join(iqr.df)
    
    t4=subset(t3,t3$value<t3$iqr)
    t4
    
    uf.ps_max <- ddply(t4, "Event", function (x) (max(fivenum(x[,"value"]))))
    uf.ps_min <- ddply(t4, "Event", function (x) (min(fivenum(x[,"value"]))))
    
    new=data.table("Comp"=comp[i],"Group"=grp[j],"uf_eve"=uf.ps_max[1,2],"uf_nev"=uf.ps_max[2,2],
                   "uf_eve_min"=uf.ps_min[1,2],"uf_nev_min"=uf.ps_min[2,2])
    ufm=rbind(ufm,new)
    
  }
}
ufm
#ufm$uf_eve=(ufm$uf_eve_max-ufm$uf_eve_min)/100*110+ufm$uf_eve_max
#ufm$uf_nev=(ufm$uf_nev_max-ufm$uf_nev_min)/100*110+ufm$uf_nev_max
ufm$uf_pos=ifelse(ufm$uf_eve>=ufm$uf_nev,ufm$uf_eve,ufm$uf_nev)


ufm$Group=factor(ufm$Group,levels = c("Ulaanbaatar","Beijing","Seoul","Seosan","Noto"))


ufm$id=paste(ufm$Comp,uf$Group,sep = "_")
dt$id=paste(dt$Comp,dt$Group,sep = "_")
#dt$Noneq_t=dt_non_equal$`T-test`
#dt$p=ifelse(dt$`T-test`<dt$Noneq_t,dt$`T-test`,dt$Noneq_t)
dt
ufm=ufm %>% inner_join(dt[,c("id","non_ttest")])
ufm$sig=ifelse(ufm$non_ttest<0.001,"***",
               ifelse(ufm$non_ttest<0.01,"**",
                      ifelse(ufm$non_ttest<0.05,"*","N.S")))
ufm


ufm$xs=ifelse(ufm$Group=="Seoul",3,
              ifelse(ufm$Group=="Seosan",4,
                     ifelse(ufm$Group=="Beijing",2,
                            ifelse(ufm$Group=="Noto",5,1))))

ufm$varlab=paste(ufm$Comp,"(%)")

ufm$varlab=factor(ufm$varlab, levels = c("Lignins (%)","Condensed aromatics (%)","Tannins (%)", "Unsaturated hydrocarbons (%)",
                                         "Carbohydrates (%)","Proteins (%)", "Lipids (%)"))

ufm$uf_pos=ifelse(ufm$Group=="Ulaanbaatar",ufm$uf_pos*1.1,ufm$uf_pos)
ufm$uf_pos=ifelse(ufm$id=="Unsaturated hydrocarbons_Seosan",ufm$uf_pos+0.1,ufm$uf_pos)
ufm$uf_pos=ifelse(ufm$id=="Lipids_Seosan",ufm$uf_pos+1.1,ufm$uf_pos)

unique(ufm$id)
ps=ggplot()+
  #geom_violin(position = "dodge")+
  stat_boxplot(data=FT_mole_ev, aes(x=Group, y=value, fill=Event),geom = "errorbar",position = position_dodge(width = 0.8), width=0.4)+
  #stat_boxplot(data=subset(FT_mole_ev,FT_mole_ev$Sample!="Ulaanbaatar_8"), aes(x=Group, y=rel, fill=Event),geom = "errorbar",position = position_dodge(width = 0.8), width=0.4)+
  geom_boxplot(data=FT_mole_ev, aes(x=Group, y=value, fill=Event),position = position_dodge(width = 0.80), outlier.colour = NA)+
  #geom_boxplot(data=subset(FT_mole_ev,FT_mole_ev$Sample!="Ulaanbaatar_8"), aes(x=Group, y=rel, fill=Event),position = position_dodge(width = 0.80), outlier.colour = NA)+
  facet_wrap(varlab~., strip.position = "left", scales = "free", ncol=4)+
  geom_blank(data=mol_max_ev, aes(x=2.5,y=value*1.12))+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), lty=2)+
  geom_segment(data=ufm,aes(x=xs-0.2,xend=xs+0.2,y=uf_pos*1.03, yend=uf_pos*1.03))+
  #geom_segment(data=uf,aes(x=xs-0.2,xend=xs-0.2,y=uf_eve*1.01, yend=uf_pos*1.05))+
  #geom_segment(data=uf,aes(x=xs+0.2,xend=xs+0.2,y=uf_nev*1.01, yend=uf_pos*1.05))+
  #geom_text(data =comp_max, aes(x=1, y=x, label=Comp), col="white", size=0)+
  geom_text(data =subset(ufm,ufm$non_ttest<0.05), aes(x=Group, y=uf_pos*1.05, label=sig), col="black", size=6)+
  scale_x_discrete("")+
  scale_fill_manual(values = c("grey75","white"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 16, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.0,0.0,0.2),"cm")),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0.0,0,0.0),"cm")),
        strip.placement = "outside",
        strip.text.x = element_text(size = 20, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.2,0.2,0.2),"cm")),
        strip.text.y = element_text(size = 18, hjust = 0.5,colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.4,0.2,0.4),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 24, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.82,0.9525),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_rect(fill="white")
  )

ps+ggsave(filename("Eventcompare_mol_sig"),height = 35, width = 75, units = "cm", dpi = 300)


