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

FT_merge=fread("Datafile/FRIENDs_1st_FT.csv")
FT_merge$Sample=paste(FT_merge$Group,FT_merge$No,sep = "_")

FT_merge$C=ifelse(FT_merge$C.>0, "C","")
FT_merge$H=ifelse(FT_merge$H.>0, "H","")
FT_merge$O=ifelse(FT_merge$O.>0, "O","")
FT_merge$N=ifelse(FT_merge$N.>0, "N","")
FT_merge$S=ifelse(FT_merge$S.>0, "S","")

FT_merge=FT_merge %>% unite("Comp",c("C","H","O","N","S"),sep = "")
FT_merge$Comp=ifelse(FT_merge$O.==0, "Remainders",FT_merge$Comp)

FT_bb=fread("Datafile/WOSCbb_S&B.csv")

FT_merge_sel=FT_merge %>% filter(Group%in%c("SUL","B"))
FT_merge_sel$Group=ifelse(FT_merge_sel$Group=="SUL","Seoul","Beijing")
FT_merge_sel$Sample=paste(FT_merge_sel$Group,FT_merge_sel$No,sep = "_")

FT_merge_sel=FT_merge_sel %>% inner_join(FT_bb)

FT_m=melt(FT_merge_sel[,c("Sample","Group","No","WSOCbb","WSOCnbb","Formula","Bromo.Inty")],
          id.vars = c("Sample","Group","No","WSOCbb","WSOCnbb","Formula")) %>% 
  dcast(Sample+Group+No+WSOCbb+WSOCnbb~Formula, value.var = "value",sum)

FT_m[,1:23]

fwrite(FT_m, file = "Datafile/FT_formula.csv")

tt=fread("Datafile/WOSC_cortest.csv")
tt

##0-2)build formula table to correlation test====
tt

#1.Correlation between Envi vs Formula====
#fwrite(envi_cor_sel,"Datafile/cor_S&B.csv")

cor_SOA=fread("Datafile/FT_formula_ob3.csv")
cor_SOA[,1:23]

cor_temp=cor_SOA
grp=unique(cor_temp$Group)

df=data.frame()
for (i in 1:length(grp)) {
  #i=1
  temp=as.data.frame(subset(cor_temp,cor_temp$Group==grp[i]))
  for (j in 4:7) {
    #j=4
    for (k in 8:dim(temp)[2]) {
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
  fwrite(df, file = "cortest.csv")
}

df=fread("cortest.csv")

df_sel=subset(df,df$p<0.05)
df_sel_d=dcast(df_sel,Group+Formula~Envi,value.var = "rho",sum)

FT_cor=FT_merge_sel %>% inner_join(df_sel,by = c("Group","Formula"))
FT_cor

fwrite(FT_cor,file = "FT_cor.csv")

df$variable=ifelse(df$Envi=="WSOCbb","cor_WSOCbb","cor_WSOCnbb")

vk_cor=fread("Datafile/FT_cor_WSOCbb&nbb.csv")

vk_cor_m=melt(vk_cor[,c("Group","Formula","O.C","H.C","cor_WSOCbb","cor_WSOCnbb","NBB")],
              id.vars = c("Group","Formula","O.C","H.C","NBB"))

vk_cor_mp=vk_cor_m %>% inner_join(df[,-"Envi"])

vk_cor_mp$Group=factor(vk_cor_mp$Group,levels = c("Seoul","Beijing"))

vk_cor_mp$Grouplab=vk_cor_mp$Group
  
vk_cor_mp$Grouplab=factor(vk_cor_mp$Grouplab,levels = c("Seoul","Beijing"),
                       labels=c(expression(bold("Seoul")),
                                expression(bold("Beijing"))))

vk_cor_mp$varlab=vk_cor_mp$variable

vk_cor_mp$varlab=factor(vk_cor_mp$varlab,levels = c("cor_WSOCbb","cor_WSOCnbb"),
                          labels = c(expression(bold("WSOC"["bb"])),
                                     expression(bold("WSOC"["nbb"]))))

vk_cor_m_0=subset(vk_cor_mp,abs(vk_cor_mp$value)>0.05)
vk_cor_m_0=vk_cor_m_0[order(vk_cor_m_0$rho),]


vk_cor_m_r=subset(vk_cor_m_0,vk_cor_m_0$p<0.01)
vk_cor_m_r=subset(vk_cor_m_r,abs(vk_cor_m_r$value)>0.45)


library(colorRamps)

vk_cor_m_r=vk_cor_m_r[order(vk_cor_m_r$rho),]


ggplot(vk_cor_m_r, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(varlab~Grouplab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                    labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
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
  ggsave(filename("vk_cor_0.4_p0.01"),height = 40, width = 40, units = "cm", dpi = 700)

vk_cor_m_0_s=subset(vk_cor_m_0,vk_cor_m_0$Group=="Seoul")
vk_cor_m_1_s=subset(vk_cor_m_0_s,abs(vk_cor_m_0_s$rho)>0.45)
vk_cor_m_1_s=subset(vk_cor_m_1_s,vk_cor_m_1_s$p<0.005)

ggplot(vk_cor_m_1_s, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(varlab~Grouplab, labeller = label_parsed )+
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
                            barwidth = 15, barheight = 2.5))+
  ggsave(filename("vk_cor_sul"),height = 20, width = 20, units = "cm", dpi = 700)

vk_cor_m_0_s_bb=subset(vk_cor_m_0_s,vk_cor_m_0_s$NBB=="BB")
length(unique(vk_cor_m_0_s_bb$Formula))

vk_cor_m_0_s_nbb=subset(vk_cor_m_0_s,vk_cor_m_0_s$NBB=="NBB")
length(unique(vk_cor_m_0_s_nbb$Formula))

vk_cor_m_2_s=rbind.data.frame(vk_cor_m_0_s_bb,vk_cor_m_0_s_nbb)
length(unique(vk_cor_m_2_s$Formula))
vk_cor_m_2_s
table(vk_cor_m_2_s$NBB)

fwrite(vk_cor_m_2_s,file = "test.csv")

vk_cor_m_3_s=fread(file ="Datafile/WSOCbb&nbb_sul.csv")

vk_cor_m_3_s$Grouplab=vk_cor_m_3_s$Group
vk_cor_m_3_s$Grouplab=factor(vk_cor_m_3_s$Grouplab,levels = c("Seoul","Beijing"),
                          labels=c(expression(bold("Seoul")),
                                   expression(bold("Beijing"))))

vk_cor_m_3_s$varlab=vk_cor_m_3_s$variable
vk_cor_m_3_s$varlab=factor(vk_cor_m_3_s$varlab,levels = c("cor_WSOCbb","cor_WSOCnbb"),
                        labels = c(expression(bold("WSOC"["bb"])),
                                   expression(bold("WSOC"["nbb"]))))
table(vk_cor_m_3_s$variable)

vk_cor_m_3_s
vk_cor_m_4_s=subset(vk_cor_m_3_s,abs(vk_cor_m_3_s$rho)>0.400)


table(vk_cor_m_4_s$variable)

ggplot(vk_cor_m_4_s, aes(x=O.C, y=H.C,col=value))+
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
        legend.title = element_text(size = 26, vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
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
  ggsave(filename("vk_cor_sul2"),height = 20, width = 40, units = "cm", dpi = 300)



vk_cor_m_0_b=subset(vk_cor_m_0,vk_cor_m_0$Group=="Beijing")
vk_cor_m_1_b=subset(vk_cor_m_0_b,abs(vk_cor_m_0_b$rho)>0.45)
vk_cor_m_1_b=subset(vk_cor_m_1_b,vk_cor_m_1_b$p<0.001)
vk_cor_m_1_b

table(vk_cor_m_1_b$variable)
table(vk_cor_m_1_b$NBB)

ggplot(vk_cor_m_1_b, aes(x=O.C, y=H.C,col=value))+
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
        legend.title = element_text(size = 26, vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
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
  ggsave(filename("vk_cor_b"),height = 20, width = 40, units = "cm", dpi = 300)


vk_cor_m_0_b_bb=subset(vk_cor_m_0_b,vk_cor_m_0_b$NBB=="BB")
length(unique(vk_cor_m_0_b_bb$Formula))

vk_cor_m_0_b_nbb=subset(vk_cor_m_0_b,vk_cor_m_0_b$NBB=="NBB")
length(unique(vk_cor_m_0_s_nbb$Formula))

vk_cor_m_2_b=rbind.data.frame(vk_cor_m_0_b_bb,vk_cor_m_0_b_nbb)
length(unique(vk_cor_m_2_b$Formula))
vk_cor_m_2_b
table(vk_cor_m_2_b$NBB)
table(vk_cor_m_3_b$NBB)

fwrite(vk_cor_m_2_b,file = "WSOC_beijing_sel.csv")

vk_cor_m_3_b=fread(file ="Datafile/WSOCbb&nbb_b.csv")
#vk_cor_m_3_b=vk_cor_m_2_b
table(vk_cor_m_3_b$NBB)
table(vk_cor_m_3_b$variable)

vk_cor_m_3_b$Grouplab=vk_cor_m_3_b$Group
vk_cor_m_3_b$Grouplab=factor(vk_cor_m_3_b$Grouplab,levels = c("Seoul","Beijing"),
                             labels=c(expression(bold("Seoul")),
                                      expression(bold("Beijing"))))

vk_cor_m_3_b$varlab=vk_cor_m_3_b$variable
vk_cor_m_3_b$varlab=factor(vk_cor_m_3_b$varlab,levels = c("cor_WSOCbb","cor_WSOCnbb"),
                           labels = c(expression(bold("WSOC"["bb"])),
                                      expression(bold("WSOC"["nbb"]))))

vk_cor_m_3_b
table(vk_cor_m_3_b$NBB)

vk_cor_m_4_b=subset(vk_cor_m_3_b,abs(vk_cor_m_3_b$rho)>0.400)

table(vk_cor_m_4_b$NBB)

ggplot(vk_cor_m_4_b, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(Grouplab~varlab   , labeller = label_parsed )+
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
        legend.title = element_text(size = 26, vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
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
  ggsave(filename("vk_cor_b2"),height = 20, width = 40, units = "cm", dpi = 300)

