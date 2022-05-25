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
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")
insert_minor_raw <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", n_minor) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_raw(c("0","0.5","1.0","1.5","2.0","2.5"),8)


insert_minor_1 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_2 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}

#insert minor1
sapply(c(0,5,10,15,20,25), function(x) c(x, seq(1,4,1)))
sapply(c("0","0.5","1.0","1.5","2.0","2.5"), function(x) c(x, rep("", 4) ))
length(sapply(c("0","0.5","1.0","1.5","2.0","2.5"), function(x) c(x, rep("", 4) )))
sapply(c("0","0.5","1.0","1.5","2.0","2.5"), function(x) c(x, rep("", 4) ))[1:(30-4)]

#insert minor2
sapply(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"), 
       function(x) c(x, rep("", 1) ))
length(sapply(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"), 
              function(x) c(x, rep("", 1) )))
sapply(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"), 
       function(x) c(x, rep("", 1) ))[1:(24-1)]

insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),1)
insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),1)


#0.Load data====
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

table(FT_merge$Group)
FT_merge_sn=FT_merge %>% filter(Group%in%c("SS","NT"))
FT_merge_sn$Group=factor(FT_merge_sn$Group,levels = c("SS","NT"),labels = c("Seosan","Noto"))

#1.Chemical Composition distribution=====
##1-1)build chemical composition data====
ftsn_inty_all=as.data.table(aggregate(FT_merge_sn$Bromo.Inty, by=list(`Sample`=FT_merge_sn$Sample), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Tot"))


FT_merge_sn=FT_merge_sn %>% inner_join(ftsn_inty_all,by = "Sample")
FT_merge_sn$rel=FT_merge_sn$Bromo.Inty/FT_merge_sn$Tot*100
FT_merge_sn$Sample=paste(FT_merge_sn$Group,FT_merge_sn$No,sep = "_")
#aggregate(FT_merge_sn$rel, by=list(`Sample`=FT_merge_sn$Sample), FUN=sum)
FT_merge_sn$Freq=1
f_count=as.data.table(aggregate(FT_merge_sn$Freq, by=list(`Group`=FT_merge_sn$Group,
                                                          Formula=FT_merge_sn$Formula), FUN=sum))


FT_merge_sn_m=melt(FT_merge_sn[,c("Sample","Group","Formula","rel")],id.vars = c("Sample","Group","Formula"))

FT_merge_sn_m=FT_merge_sn_m %>% inner_join(f_count)
FT_merge_sn_m_sel=FT_merge_sn_m %>% filter(x>3)


FT_merge_sn_m_d=dcast(FT_merge_sn_m_sel,Sample~Formula, sum)
FT_merge_sn_m_d[,1:23]

FT_envi=fread("Datafile/FRIENDs_envi_sel.csv")
#1. Environmetal data analysis=====
FT_envi_sel=FT_envi %>% filter(Group%in%c("Seosan","Noto"))

FT_envi_sel$WSOCbb=2.94*FT_envi_sel$Levoglucosan
FT_envi_sel$WSOCnbb=FT_envi_sel$WSOC-FT_envi_sel$WSOCbb

FT_envi_sel

FT_fm=FT_envi_sel[,c("Sample","Group","No","WSOCbb","WSOCnbb")] %>% inner_join(FT_merge_sn_m_d)
FT_fm[,1:24]

grp=as.vector(unique(FT_fm$Group))
df=data.frame()
for (i in 1:length(grp)) {
  #i=1
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

cor_vk_sn=fread("cortest.csv")
#fwrite(cor_vk_sn,file = "cor_ss&nt.csv")

cor_vk_sn=subset(cor_vk_sn,cor_vk_sn$p<0.05)
cor_vk_sn

cor_vk_sn_d=dcast(cor_vk_sn, Group+Formula~Envi, sum, value.var = c("rho","p"))
cor_vk_sn_d

CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(cor_vk_sn_d$Formula)

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

cor_vk_sn_d=cor_vk_sn_d %>% inner_join(fm)
cor_vk_sn_d$O.C=cor_vk_sn_d$O/cor_vk_sn_d$C
cor_vk_sn_d$H.C=cor_vk_sn_d$H/cor_vk_sn_d$C

cor_vk_sn_d

cor_vk_sn_d$C1=ifelse(cor_vk_sn_d$C>0, "C","")
cor_vk_sn_d$H1=ifelse(cor_vk_sn_d$H>0, "H","")
cor_vk_sn_d$O1=ifelse(cor_vk_sn_d$O>0, "O","")
cor_vk_sn_d$N1=ifelse(cor_vk_sn_d$N>0, "N","")
cor_vk_sn_d$S1=ifelse(cor_vk_sn_d$S>0, "S","")

cor_vk_sn_d=cor_vk_sn_d %>% unite("Comp",c("C1","H1","O1","N1","S1"),sep = "")
cor_vk_sn_d$Comp=ifelse(cor_vk_sn_d$O==0, "Remainders",cor_vk_sn_d$Comp)


cor_vk_sn_d$Type=ifelse(abs(cor_vk_sn_d$rho_WSOCbb)>abs(cor_vk_sn_d$rho_WSOCnbb),"BB","NBB")
table(cor_vk_sn_d$Type)
cor_vk_sn_d

cor_vk_sn_d$variable=ifelse(cor_vk_sn_d$Type=="BB","WSOCbb","WSOCnbb")
cor_vk_sn_d$value=ifelse(cor_vk_sn_d$Type=="BB",cor_vk_sn_d$rho_WSOCbb,cor_vk_sn_d$rho_WSOCnbb)
cor_vk_sn_d$p=ifelse(cor_vk_sn_d$Type=="BB",cor_vk_sn_d$p_WSOCbb,cor_vk_sn_d$p_WSOCnbb)
cor_vk_sn_d

cor_vk_sn_d_sel=cor_vk_sn_d[,c("Group","Comp","Formula","O.C","H.C","Type","variable","value","p")]

cor_vk_sn_d_sel
#cor_vk_d_sel=subset(cor_vk_d_sel,abs(cor_vk_d_sel$value)>0.05)

dim(cor_vk_sn_d_sel)

cor_vk_sn_d_sel$Grouplab=factor(cor_vk_sn_d_sel$Group,levels = c("Seosan","Noto"),
                             labels=c(expression(bold("Seosan")),
                                      expression(bold("Noto"))))

cor_vk_sn_d_sel$varlab=factor(cor_vk_sn_d_sel$variable,levels = c("WSOCbb","WSOCnbb"),
                           labels = c(expression(bold("WSOC"["BB"])),
                                      expression(bold("WSOC"["NBB"]))))

cor_vk_sn_d_sel_ss=subset(cor_vk_sn_d_sel,cor_vk_sn_d_sel$Group=="Seosan")
cor_vk_sn_d_sel_ss=cor_vk_sn_d_sel_ss[order(cor_vk_sn_d_sel_ss$value),]


tt=as.data.frame(table(cor_vk_sn_d_sel_ss$Type)) 
tt$varlab=factor(tt$Var1,levels = c("BB","NBB"),
                 labels = c(expression(bold("WSOC"["BB"])),
                            expression(bold("WSOC"["NBB"]))))
tt$O.C=0.1
tt$H.C=0.15

ggplot()+
  geom_point(data=cor_vk_sn_d_sel_ss, aes(x=O.C, y=H.C,col=value),size=2.7)+
  #geom_text(data=tt, aes(x=O.C, y=H.C,label=Freq), size=8)+
  geom_text(data=tt, aes(x=0.05, y=0.14,label=paste0(expression(italic("n=")))),parse=T, size=8)+
  geom_text(data=tt, aes(x=0.15, y=H.C,label=Freq),parse=T, size=8)+
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

cor_vk_sn_d_sel_nt=subset(cor_vk_sn_d_sel,cor_vk_sn_d_sel$Group=="Noto")
cor_vk_sn_d_sel_nt=cor_vk_sn_d_sel_nt[order(cor_vk_sn_d_sel_nt$value),]


tt2=as.data.frame(table(cor_vk_sn_d_sel_nt$Type)) 
tt2$varlab=factor(tt$Var1,levels = c("BB","NBB"),
                 labels = c(expression(bold("WSOC"["BB"])),
                            expression(bold("WSOC"["NBB"]))))
tt2$O.C=0.1
tt2$H.C=0.15

ggplot()+
  geom_point(data=cor_vk_sn_d_sel_nt, aes(x=O.C, y=H.C,col=value),size=2.7)+
  geom_text(data=tt2, aes(x=0.05, y=0.14,label=paste0(expression(italic("n=")))),parse=T, size=8)+
  geom_text(data=tt2, aes(x=0.15, y=H.C,label=Freq),parse=T, size=8)+
  #facet_grid(Grouplab~varlab, labeller = label_parsed )+
  facet_grid(varlab~Grouplab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.2,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),1))+
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
  ggsave(filename("vk_cor_nt_rel_t"),height = 40, width = 20, units = "cm", dpi = 300)
''