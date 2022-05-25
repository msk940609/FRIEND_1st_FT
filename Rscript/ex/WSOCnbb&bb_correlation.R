library(PerformanceAnalytics)


chart.Correlation(iris[1:4])

FT_bb=fread("Datafile/WOSCbb_S&B.csv")

FT_bb$WSOCbb=FT_bb$Levo*0.7*4.18*FT_bb$WSOC
FT_bb$WSOCnbb=FT_bb$WSOC-FT_bb$WSOCbb

FT_bb_s=subset(FT_bb,FT_bb$Group=="Seoul")
FT_bb_b=subset(FT_bb,FT_bb$Group=="Beijing")

chart.Correlation(FT_bb_s[,8:10], method = "pearson")

cor.test(FT_bb_s$WSOCbb,FT_bb_s$WSOCnbb, method = "spearman", exact = F)
cor.test(FT_bb_b$WSOCbb,FT_bb_b$WSOCnbb, method = "spearman", exact = F)

chart.Correlation(FT_bb_s[,8:10], method = "pearson")
chart.Correlation(FT_bb_b[,8:10], method = "pearson")

chart.Correlation(FT_bb_s[,8:10], method = "spearman")
chart.Correlation(FT_bb_b[,8:10], method = "spearman")

FT_bb
FT_bb_m=melt(FT_bb[,c(1,2,3,8,9,10:17)], id.vars = c("Sample","Group","No","WSOCbb","WSOCnbb"), variable.name = "Envi", value.name = "Conc") %>% 
  melt(id.vars=c("Sample","Group","No","Envi","Conc"),variable.name="WSOC", na.rm = T)
FT_bb_m

FT_bb_m$No=as.numeric(ft_inty_env_m$No)
FT_bb_m$mark="no"
FT_bb_m$mark=ifelse(FT_bb_m$No==5,"12/19",FT_bb_m$mark)
FT_bb_m$mark=ifelse(FT_bb_m$No==6,"12/20",FT_bb_m$mark)
FT_bb_m$mark=ifelse(FT_bb_m$No==7,"12/21",FT_bb_m$mark)
FT_bb_m$mark=ifelse(FT_bb_m$No==8,"12/22",FT_bb_m$mark)
FT_bb_m$mark=ifelse(FT_bb_m$No==9,"12/23",FT_bb_m$mark)

FT_bb_m

grp=unique(FT_bb_m$Group)
grp
comp=unique(FT_bb_m$WSOC)
envi=unique(FT_bb_m$Envi)


dt=data.table()
for (i in 1:length(grp)) {
  #i=1
  temp=subset(FT_bb_m,FT_bb_m$Group==grp[i])
  
  for (k in 1:length(comp)) {
    #k=2
    tmp=subset(temp,temp$WSOC==comp[k])
    
    for (l in 1:length(envi)) {
      tmp2=subset(tmp,tmp$Envi==envi[l])
      
      if(envi[l]=="O3"&grp[i]=="Ulaanbaatar") {
        tmp2$Conc=0
        new_row=data.table("Group"=grp[i],"Comp"=comp[k],"R2"=0,"Envi"=envi[l],"Cor"=0, "P"=1)
        
      } else {
        model=lm(Conc~value,tmp2)
        a=summary(model)
        #cor=cor.test(tmp2$value,tmp2$Conc, method = "pearson")
        cor=cor.test(tmp2$value,tmp2$Conc, method = "spearman", exact = F)
        #    cor
        
        new_row=data.table("Group"=grp[i],"Comp"=comp[k],"R2"=round(a$r.squared,8),"Envi"=envi[l],"Cor"=cor$estimate, "P"=round(cor$p.value,3))
        
      }
      
      dt=rbind(dt,new_row)
    }
  }
}

dt

dt$sig=ifelse(dt$P<0.05,ifelse(dt$P<0.01,ifelse(dt$P<0.001,"***","**"),"*"),"")

fwrite(dt,file="WSOCcor&Envi.csv")



FT_bb_m_s=subset(FT_bb_m,FT_bb_m$Group=="Seoul")
FT_bb_m_s1=FT_bb_m_s %>% filter(Envi%in%c("PM2.5","NH4","NO3","SO4")) %>% droplevels()

FT_bb_m_s1$Envilab=FT_bb_m_s1$Envi
FT_bb_m_s1$Envilab=factor(FT_bb_m_s1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                                labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))

FT_bb_m_s1$Complab=FT_bb_m_s1$WSOC

FT_bb_m_s1$Complab=factor(FT_bb_m_s1$Complab, levels = c("WSOCbb","WSOCnbb"),
                                labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")"))))


val_max_s1=as.data.table(aggregate(FT_bb_m_s1$Conc, 
                                   by=list(`Comp`=FT_bb_m_s1$WSOC,Envi=FT_bb_m_s1$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))
val_max_s1$Envilab=val_max_s1$Envi

val_max_s1$Envilab=factor(val_max_s1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                          labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))


ggplot(FT_bb_m_s1, aes(x=value,y=Conc))+
  geom_smooth(formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(col="royalblue2", size=3)+
  geom_blank(data=val_max_s1, aes(x=0.8,y=Conc*1.32))+
  geom_point(data=subset(FT_bb_m_s1,FT_bb_m_s1$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(FT_bb_m_s1,FT_bb_m_s1$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  #facet_wrap(.~Complab*Envilab, scales = "free",labeller = label_parsed,ncol=4)+
  facet_grid(Envilab~Complab, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "")+
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
  ggsave(filename("compvsEnvi_sul1"),height = 30, width = 20, units = "cm", dpi = 300)


FT_bb_m_s2=FT_bb_m_s %>% filter(Envi%in%c("O3","CO","SO2","NO")) %>% droplevels()
FT_bb_m_s2$Envilab=FT_bb_m_s2$Envi

FT_bb_m_s2$Envilab=factor(FT_bb_m_s2$Envilab, levels = c("O3","CO","NO","SO2"),
                                labels = c(expression(bold("O"["3"]~"(ppb)")),
                                           expression(bold("CO (ppb)")),
                                           expression(bold("NO (ppb)")),
                                           expression(bold("SO"["2"]~"(ppb)"))))
FT_bb_m_s2$Complab=FT_bb_m_s2$WSOC

FT_bb_m_s2$Complab=factor(FT_bb_m_s2$Complab, levels = c("WSOCbb","WSOCnbb"),
                          labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")"))))

val_max_s2=as.data.table(aggregate(FT_bb_m_s2$Conc, 
                                   by=list(`Comp`=FT_bb_m_s2$Comp,Envi=FT_bb_m_s2$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))

val_max_s2$Envilab=val_max_s2$Envi

val_max_s2$Envilab=factor(val_max_s2$Envilab, levels = c("O3","CO","NO","SO2"),
                          labels = c(expression(bold("O"["3"]~"(ppb)")),
                                     expression(bold("CO (ppb)")),
                                     expression(bold("NO (ppb)")),
                                     expression(bold("SO"["2"]~"(ppb)"))))

val_max_s2$Conc=ifelse(val_max_s2$Envi=="SO2",2.0,val_max_s2$Conc)

ggplot(FT_bb_m_s2, aes(x=value,y=Conc))+
  geom_smooth(formula = y~x, method = "lm", se = T, col="black")+
  geom_point(col="royalblue2", size=3)+
  geom_blank(data=val_max_s2, aes(x=0.8,y=Conc*1.32))+
  geom_point(data=subset(FT_bb_m_s2,FT_bb_m_s2$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(FT_bb_m_s2,FT_bb_m_s2$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_wrap(.~Complab*Envilab, scales = "free",labeller = label_parsed,ncol=4)+
  facet_grid(Envilab~Complab, scales = "free",labeller = label_parsed, switch = "y")+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  #scale_x_continuous(name = "Proportion of chemical composition (%)")+
  scale_x_continuous(name = "")+
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
  ggsave(filename("compvsEnvi_sul2"),height = 30, width = 20, units = "cm", dpi = 300)

##Beijing====
FT_bb_m_b=subset(FT_bb_m,FT_bb_m$Group=="Beijing")
FT_bb_m_b1=FT_bb_m_b %>% filter(Envi%in%c("PM2.5","NH4","NO3","SO4")) %>% droplevels()

FT_bb_m_b1$Envilab=FT_bb_m_b1$Envi
FT_bb_m_b1$Envilab=factor(FT_bb_m_b1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                          labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))

chart.Correlation(FT_bb_b[,c(7,8,9,10,13,12,11,14,15,16)], method = "spearman")

FT_bb_m_b1$Complab=FT_bb_m_b1$WSOC
FT_bb_m_b1$Complab=factor(FT_bb_m_b1$Complab, levels = c("WSOCbb","WSOCnbb"),
                          labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")"))))

val_max_b1=as.data.table(aggregate(FT_bb_m_b1$Conc, 
                                   by=list(`Comp`=FT_bb_m_b1$WSOC,Envi=FT_bb_m_b1$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))

val_max_b1$Envilab=val_max_b1$Envi

val_max_b1$Envilab=factor(val_max_b1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                          labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))


ggplot(FT_bb_m_b1, aes(x=value,y=Conc))+
  geom_smooth(formula = y~x, method = "lm", se = T, col="black",lyt=2)+
  geom_point(col="orangered3", size=3)+
  geom_blank(data=val_max_b1, aes(x=0.8,y=Conc*1.32))+
  geom_point(data=subset(FT_bb_m_b1,FT_bb_m_b1$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(FT_bb_m_b1,FT_bb_m_b1$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  #facet_wrap(.~Complab*Envilab, scales = "free",labeller = label_parsed,ncol=4)+
  facet_grid(Envilab~Complab, scales = "free",labeller = label_parsed, switch = "y")+
  scale_x_continuous(name = "")+
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
  ggsave(filename("compvsEnvi_b1"),height = 30, width = 20, units = "cm", dpi = 300)


FT_bb_m_b2=FT_bb_m_b %>% filter(Envi%in%c("O3","CO","SO2","NO")) %>% droplevels()
FT_bb_m_b2$Envilab=FT_bb_m_b2$Envi

FT_bb_m_b2$Envilab=factor(FT_bb_m_b2$Envilab, levels = c("O3","CO","NO","SO2"),
                          labels = c(expression(bold("O"["3"]~"(ppb)")),
                                     expression(bold("CO (ppb)")),
                                     expression(bold("NO (ppb)")),
                                     expression(bold("SO"["2"]~"(ppb)"))))
FT_bb_m_b2$Complab=FT_bb_m_b2$WSOC

FT_bb_m_b2$Complab=factor(FT_bb_m_b2$Complab, levels = c("WSOCbb","WSOCnbb"),
                          labels = c(expression(bold("WSOC"["bb"]~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("WSOC"["nbb"]~"("*"\u03bcg/"*m^"3"*")"))))

val_max_b2=as.data.table(aggregate(FT_bb_m_b2$Conc, 
                                   by=list(`Comp`=FT_bb_m_b2$Comp,Envi=FT_bb_m_b2$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))

val_max_b2$Envilab=val_max_b2$Envi

val_max_b2$Envilab=factor(val_max_b2$Envilab, levels = c("O3","CO","NO","SO2"),
                          labels = c(expression(bold("O"["3"]~"(ppb)")),
                                     expression(bold("CO (ppb)")),
                                     expression(bold("NO (ppb)")),
                                     expression(bold("SO"["2"]~"(ppb)"))))
val_max_b2$Complab=val_max_b2$Comp

val_max_b2$Conc=ifelse(val_max_b2$Envi=="CO",1050,val_max_b2$Conc)
#val_max_b2$Conc=ifelse(val_max_b2$Envi=="NO",55,val_max_s2$Conc)

ggplot(FT_bb_m_b2, aes(x=value,y=Conc))+
  geom_smooth(formula = y~x, method = "lm", se = T, col="black")+
  geom_point(col="orangered3", size=3)+
  geom_blank(data=val_max_b2, aes(x=0.3,y=Conc*1.3))+
  geom_point(data=subset(FT_bb_m_b2,FT_bb_m_b2$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(FT_bb_m_b2,FT_bb_m_b2$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_wrap(.~Complab*Envilab, scales = "free",labeller = label_parsed,ncol=4)+
  facet_grid(Envilab~Complab, scales = "free",labeller = label_parsed, switch = "y")+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  #scale_x_continuous(name = "Proportion of chemical composition (%)")+
  scale_x_continuous(name = "")+
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
  ggsave(filename("compvsEnvi_b2"),height = 30, width = 20, units = "cm", dpi = 300)


