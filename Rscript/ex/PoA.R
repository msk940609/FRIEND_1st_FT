library(ggplot2)
library(data.table)
library(dplyr)
library(tidyverse)
options(java.parameters = "- Xmx8192m")
library(extrafont)
loadfonts(device="win")
library(ggrepel)

source("Rscript/func_filename.R")
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")
insert_minor_1 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_2 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}

FT_soa=fread("Datafile/FRIENDs_1st_FT_SOA_S&B.csv")
FT_soa$Group=ifelse(FT_soa$Group=="SUL","Seoul","Beijing")
FT_soa$Sample=paste(FT_soa$Group,FT_soa$No,sep = "_")

FT_soa$C=ifelse(FT_soa$C.>0, "C","")
FT_soa$H=ifelse(FT_soa$H.>0, "H","")
FT_soa$O=ifelse(FT_soa$O.>0, "O","")
FT_soa$N=ifelse(FT_soa$N.>0, "N","")
FT_soa$S=ifelse(FT_soa$S.>0, "S","")

FT_soa=FT_soa %>% unite("Comp",c("C","H","O","N","S"),sep = "")
FT_soa$Comp=ifelse(FT_soa$O.==0, "Remainders",FT_soa$Comp)
FT_soa=subset(FT_soa,FT_soa$Comp!="Remainders") %>% droplevels()

wsoc=fread("Datafile/WOSCbb_S&B.csv")

FT_soa_mono=(aggregate(FT_soa$Mono.Inty, by=list(`Sample`=FT_soa$Sample), FUN=sum))%>% 
  `colnames<-`(c("Sample","Tot.mono"))

FT_soa=FT_soa %>% inner_join(FT_soa_mono)

FT_soa$rel=FT_soa$Mono.Inty/FT_soa$Tot.mono*100
#aggregate(FT_soa$rel, by=list(`Sample`=FT_soa$Sample), FUN=sum)

FT_soa=FT_soa %>% inner_join(wsoc[,c("Sample","WSOCbb","WSOCnbb","POA","SOA")], by="Sample")
FT_soa

FT_soa_m=melt(FT_soa[,c("Sample","Group","No","Formula","WSOCbb","WSOCnbb","POA","SOA","rel")], 
              id.vars = c("Sample","Group","No","Formula","WSOCbb","WSOCnbb","POA","SOA"))

FT_soa_m
#FT_soa$Freq=1
#FT_freq=(aggregate(FT_soa$Freq, by=list(`Group`=FT_soa$Group,Formula=FT_soa$Formula), FUN=sum))
#FT_freq_sel=FT_freq %>% filter(x>3)
#FT_freq_sel$ob3=1

#dim(FT_soa_m)
#dim(FT_soa_m2)
#dim(FT_freq_sel)

#FT_soa_m_sel=FT_soa_m %>% inner_join(FT_freq_sel)
#FT_soa_m_sel=dcast(FT_soa_m_sel, Sample+Group+No+WSOCbb+WSOCnbb+POA+SOA~Formula,value.var = "value",sum)
#FT_soa_m_sel[,1:23]

FT_soa_m
FT_soa_m_d=dcast(FT_soa_m, Sample+Group+No+WSOCbb+WSOCnbb+POA+SOA~Formula,value.var = "value",sum)
FT_soa_m_d[,1:23]

#fwrite(FT_soa_m_d,"Datafile/FT_soa_m.csv")

FT_soa_m_d=fread("Datafile/FT_soa_m_ob3.csv")
wsoc=fread("Datafile/WOSCbb_S&B.csv")
FT_soa_m_d$POA=wsoc$POA
FT_soa_m_d[,1:23]
cor_temp=FT_soa_m_d

cor_temp[,1:23]

grp=unique(cor_temp$Group)

df=data.frame()
for (i in 1:length(grp)) {
  #i=1
  temp=as.data.frame(subset(cor_temp,cor_temp$Group==grp[i]))
  for (j in 4:7) {
    #j=4
    for (k in 8:dim(temp)[2]) {
      #k=8  
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


df=fread("cortest_inr.csv")

unique(df$Envi)

###df=====

cor_vk_df=df
cor_vk_df

cor_vk_df=subset(cor_vk_df,cor_vk_df$p<0.05)

cor_vk_df


CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(cor_vk_df$Formula)

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
cor_vk_df=as.data.table(cor_vk_df)

cor_vk_d=dcast(cor_vk_df, Group+Formula~Envi, sum, value.var = c("rho","p"))
cor_vk_d

cor_vk_d=cor_vk_d %>% inner_join(fm)
cor_vk_d$O.C=cor_vk_d$O/cor_vk_d$C
cor_vk_d$H.C=cor_vk_d$H/cor_vk_d$C


cor_vk_d_m=melt(cor_vk_d[,c("Group","Formula","O.C","H.C","rho_WSOCbb","rho_WSOCnbb",
                            "p_WSOCbb","p_WSOCnbb")],
                 id.vars=c("Group","Formula","O.C","H.C","rho_WSOCbb","rho_WSOCnbb"),variable.name = "p-var",
                 value.name = "p") %>%
  melt(id.vars=c("Group","Formula","O.C","H.C","p-var","p"))

cor_vk_d_m
dim(cor_vk_d_m)
cor_vk_d_m_sel=subset(cor_vk_d_m,abs(cor_vk_d_m$value)>0.001)
dim(cor_vk_d_m_sel)


cor_vk_d_m_sel
cor_vk_d_m_sel$Grouplab=cor_vk_d_m_sel$Group
cor_vk_d_m_sel$Grouplab=factor(cor_vk_d_m_sel$Grouplab,levels = c("Seoul","Beijing"),
                                labels=c(expression(bold("Seoul")),
                                         expression(bold("Beijing"))))

cor_vk_d_m_sel$varlab=cor_vk_d_m_sel$variable
cor_vk_d_m_sel$varlab=factor(cor_vk_d_m_sel$varlab,levels = c("rho_WSOCbb","rho_WSOCnbb","rho_POA","rho_SOA"),
                              labels = c(expression(bold("WSOC"["bb"])),
                                         expression(bold("WSOC"["nbb"])),
                                         expression(bold("POA")),
                                         expression(bold("SOA"))
                                         ))

cor_vk_d_m_sel=cor_vk_d_m_sel[order(cor_vk_d_m_sel$value),]

cor_vk_d_m_sel_s=subset(cor_vk_d_m_sel,cor_vk_d_m_sel$Group=="Seoul")
cor_vk_d_m_sel_s

cor_vk_d_m_sel_s1=subset(cor_vk_d_m_sel_s,cor_vk_d_m_sel_s$p<0.05)

table(cor_vk_d_m_sel_s1$variable)

ggplot(cor_vk_d_m_sel_s1, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_gradientn(colors=(topo.colors(30)[4:24]),breaks = c(-0.9, 0.0, 0.9),limits = c(-1.00,1.00))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
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
        legend.title = element_text(size = 26,hjust = 0.15, vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
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
  ggsave(filename("vk_cor_s"),height = 20, width = 80, units = "cm", dpi = 700)


cor_vk_d_m_sel_b=subset(cor_vk_d_m_sel,cor_vk_d_m_sel$Group=="Beijing")
cor_vk_d_m_sel_b

cor_vk_d_m_sel_b1=subset(cor_vk_d_m_sel_b,cor_vk_d_m_sel_b$p<0.05)

table(cor_vk_d_m_sel_b1$variable)

ggplot(cor_vk_d_m_sel_b1, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_gradientn(colors=(topo.colors(30)[4:24]),breaks = c(-0.9, 0.0, 0.9),limits = c(-1.00,1.00))+
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
        legend.title = element_text(size = 26,hjust = 0.15, vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
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
  ggsave(filename("vk_cor_b"),height = 20, width = 80, units = "cm", dpi = 700)


###specialization=====

####nbb====
#cor_vk_d_nbb=cor_vk_d[,-c("rho_POA","rho_SOA","p_POA","p_SOA")]
cor_vk_d_nbb=cor_vk_d

cor_vk_d_nbb$Type=ifelse(abs(cor_vk_d_nbb$rho_WSOCbb)>abs(cor_vk_d_nbb$rho_WSOCnbb),"BB","NBB")

#cor_vk_ev$Type=ifelse(abs(abs(cor_vk_ev$rho_WSOCbb)-abs(cor_vk_ev$rho_WSOCnbb))<0.05,"Both",cor_vk_ev$Type)
table(cor_vk_d_nbb$Type)

cor_vk_d_nbb$variable=ifelse(cor_vk_d_nbb$Type=="BB","WSOCbb","WSOCnbb")
cor_vk_d_nbb$value=ifelse(cor_vk_d_nbb$Type=="BB",cor_vk_d_nbb$rho_WSOCbb,cor_vk_d_nbb$rho_WSOCnbb)
cor_vk_d_nbb$p=ifelse(cor_vk_d_nbb$Type=="BB",cor_vk_d_nbb$p_WSOCbb,cor_vk_d_nbb$p_WSOCnbb)
cor_vk_d_nbb

cor_vk_d_nbb_sel=cor_vk_d_nbb[,c("Group","Formula","O.C","H.C","Type","variable","value","p")]
cor_vk_d_nbb_sel=unique(cor_vk_d_nbb_sel)
cor_vk_d_nbb_sel
cor_vk_d_nbb_sel=subset(cor_vk_d_nbb_sel,abs(cor_vk_d_nbb_sel$value)>0.1)
cor_vk_d_nbb_sel

cor_vk_d_nbb_sel_s=subset(cor_vk_d_nbb_sel,cor_vk_d_nbb_sel$Group=="Seoul")
cor_vk_d_nbb_sel_s

table(cor_vk_d_nbb_sel_s$Type)
table(cor_vk_d_nbb_sel_s$variable)

cor_vk_d_nbb_sel_s1=subset(cor_vk_d_nbb_sel_s,abs(cor_vk_d_nbb_sel_s$value)>0.0)
cor_vk_d_nbb_sel_s1=subset(cor_vk_d_nbb_sel_s,abs(cor_vk_d_nbb_sel_s$p)<0.05)

table(cor_vk_d_nbb_sel_s1$Type)
cor_vk_d_nbb_sel_s1

cor_vk_d_nbb_sel_s1$Grouplab=cor_vk_d_nbb_sel_s1$Group
cor_vk_d_nbb_sel_s1$Grouplab=factor(cor_vk_d_nbb_sel_s1$Grouplab,levels = c("Seoul","Beijing"),
                                 labels=c(expression(bold("Seoul")),
                                          expression(bold("Beijing"))))

cor_vk_d_nbb_sel_s1$varlab=cor_vk_d_nbb_sel_s1$variable
cor_vk_d_nbb_sel_s1$varlab=factor(cor_vk_d_nbb_sel_s1$varlab,levels = c("WSOCbb","WSOCnbb"),
                               labels = c(expression(bold("WSOC"["bb"])),
                                          expression(bold("WSOC"["nbb"]))))
cor_vk_d_nbb_sel_s1_p=subset(cor_vk_d_nbb_sel_s1,cor_vk_d_nbb_sel_s1$value>0)

ggplot(cor_vk_d_nbb_sel_s1, aes(x=O.C, y=H.C,col=value))+
  #geom_point(size=2.7)+
  geom_point(data=cor_vk_d_nbb_sel_s1_p,aes(x=O.C, y=H.C,col=value),size=2.7)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  #scale_color_gradientn(colors=(topo.colors(20)[4:17]))+
  #scale_color_gradientn(colors=(topo.colors(40)[4:38]))+
  scale_color_gradientn(colors=(topo.colors(40)[28:38]))+
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
  ggsave(filename("vk_cor_s_nbb_special"),height = 20, width = 40, units = "cm", dpi = 300)

####beijing====

cor_vk_d_nbb_sel_b=subset(cor_vk_d_nbb_sel,cor_vk_d_nbb_sel$Group=="Beijing")
cor_vk_d_nbb_sel_b

table(cor_vk_d_nbb_sel_b$Type)
table(cor_vk_d_nbb_sel_b$variable)

cor_vk_d_nbb_sel_b1=subset(cor_vk_d_nbb_sel_b,abs(cor_vk_d_nbb_sel_b$value)>0.0)
cor_vk_d_nbb_sel_b1=subset(cor_vk_d_nbb_sel_b,abs(cor_vk_d_nbb_sel_b$p)<0.05)

table(cor_vk_d_nbb_sel_b1$Type)
cor_vk_d_nbb_sel_b1

cor_vk_d_nbb_sel_b1$Grouplab=cor_vk_d_nbb_sel_b1$Group
cor_vk_d_nbb_sel_b1$Grouplab=factor(cor_vk_d_nbb_sel_b1$Grouplab,levels = c("Seoul","Beijing"),
                                    labels=c(expression(bold("Seoul")),
                                             expression(bold("Beijing"))))

cor_vk_d_nbb_sel_b1$varlab=cor_vk_d_nbb_sel_b1$variable
cor_vk_d_nbb_sel_b1$varlab=factor(cor_vk_d_nbb_sel_b1$varlab,levels = c("WSOCbb","WSOCnbb"),
                                  labels = c(expression(bold("WSOC"["bb"])),
                                             expression(bold("WSOC"["nbb"]))))

ggplot(cor_vk_d_nbb_sel_b1, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  #scale_color_gradientn(colors=(topo.colors(20)[4:17]))+
  scale_color_gradientn(colors=(topo.colors(40)[4:38]))+
  #scale_color_gradientn(colors=(topo.colors(40)[28:38]))+
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
  ggsave(filename("vk_cor_b__nbb_special"),height = 20, width = 40, units = "cm", dpi = 300)


cor_vk_d_nbb_sel_b1_p=subset(cor_vk_d_nbb_sel_b1,cor_vk_d_nbb_sel_b1$value>0)
cor_vk_d_nbb_sel_b1_p

cor_vk_d_nbb_sel_b1_p_sel=cor_vk_d_nbb_sel_b1_p

cor_vk_d_nbb_sel_b1_p_sel$value=ifelse(cor_vk_d_nbb_sel_b1_p_sel$Type=="BB",
                                       ifelse(cor_vk_d_nbb_sel_b1_p_sel$O.C<0.04,cor_vk_d_nbb_sel_b1_p_sel$value*-1,
                                              cor_vk_d_nbb_sel_b1_p_sel$value),cor_vk_d_nbb_sel_b1_p_sel$value)
cor_vk_d_nbb_sel_b1_p2=subset(cor_vk_d_nbb_sel_b1_p_sel,cor_vk_d_nbb_sel_b1_p_sel$value>0)


ggplot(cor_vk_d_nbb_sel_b1_p2, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
  #scale_color_gradientn(colors=(topo.colors(20)[4:17]))+
  #scale_color_gradientn(colors=(topo.colors(40)[4:38]))+
  scale_color_gradientn(colors=(topo.colors(40)[28:38]))+
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
  ggsave(filename("vk_cor_b__nbb_special_sel"),height = 20, width = 40, units = "cm", dpi = 300)



####ec tracer====
cor_vk_d_ec=cor_vk_d[,-c("rho_WSOCbb","rho_WSOCnbb","p_WSOCbb","p_WSOCnbb")]

cor_vk_d_ec$Type=ifelse(abs(cor_vk_d_ec$rho_POA)>abs(cor_vk_d_ec$rho_SOA),"POA","SOA")

#cor_vk_ev$Type=ifelse(abs(abs(cor_vk_ev$rho_WSOCbb)-abs(cor_vk_ev$rho_WSOCnbb))<0.05,"Both",cor_vk_ev$Type)
table(cor_vk_d_ec$Type)

cor_vk_d_ec$variable=ifelse(cor_vk_d_ec$Type=="POA","POA","SOA")
cor_vk_d_ec$value=ifelse(cor_vk_d_ec$Type=="POA",cor_vk_d_ec$rho_POA,cor_vk_d_ec$rho_SOA)
cor_vk_d_ec$p=ifelse(cor_vk_d_ec$Type=="POA",cor_vk_d_ec$p_POA,cor_vk_d_ec$p_SOA)
cor_vk_d_ec

cor_vk_d_ec_sel=cor_vk_d_ec[,c("Group","Formula","O.C","H.C","Type","variable","value","p")]
cor_vk_d_ec_sel=unique(cor_vk_d_ec_sel)
cor_vk_d_ec_sel
cor_vk_d_ec_sel=subset(cor_vk_d_ec_sel,abs(cor_vk_d_ec_sel$value)>0.1)
cor_vk_d_ec_sel

cor_vk_d_ec_sel_s=subset(cor_vk_d_ec_sel,cor_vk_d_ec_sel$Group=="Seoul")
cor_vk_d_ec_sel_s

table(cor_vk_d_ec_sel_s$Type)
table(cor_vk_d_ec_sel_s$variable)

cor_vk_d_ec_sel_s1=subset(cor_vk_d_ec_sel_s,abs(cor_vk_d_ec_sel_s$value)>0.0)
cor_vk_d_ec_sel_s1=subset(cor_vk_d_ec_sel_s,abs(cor_vk_d_ec_sel_s$p)<0.05)

table(cor_vk_d_ec_sel_s1$Type)
cor_vk_d_ec_sel_s1

cor_vk_d_ec_sel_s1$Grouplab=cor_vk_d_ec_sel_s1$Group
cor_vk_d_ec_sel_s1$Grouplab=factor(cor_vk_d_ec_sel_s1$Grouplab,levels = c("Seoul","Beijing"),
                                    labels=c(expression(bold("Seoul")),
                                             expression(bold("Beijing"))))

cor_vk_d_ec_sel_s1$varlab=cor_vk_d_ec_sel_s1$variable
cor_vk_d_ec_sel_s1$varlab=factor(cor_vk_d_ec_sel_s1$varlab,levels = c("POA","SOA"),
                                  labels = c(expression(bold("POC")),
                                             expression(bold("SOC"))))

ggplot(cor_vk_d_ec_sel_s1, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
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
  ggsave(filename("vk_cor_s_ec_special"),height = 20, width = 40, units = "cm", dpi = 300)

####beijing====
cor_vk_d_ec_sel_b=subset(cor_vk_d_ec_sel,cor_vk_d_ec_sel$Group=="Beijing")
cor_vk_d_ec_sel_b

table(cor_vk_d_ec_sel_b$Type)
table(cor_vk_d_ec_sel_b$variable)

cor_vk_d_ec_sel_b1=subset(cor_vk_d_ec_sel_b,abs(cor_vk_d_ec_sel_b$value)>0.0)
cor_vk_d_ec_sel_b1=subset(cor_vk_d_ec_sel_b,abs(cor_vk_d_ec_sel_b$p)<0.05)

table(cor_vk_d_ec_sel_b1$Type)
cor_vk_d_ec_sel_b1

cor_vk_d_ec_sel_b1$Grouplab=cor_vk_d_ec_sel_b1$Group
cor_vk_d_ec_sel_b1$Grouplab=factor(cor_vk_d_ec_sel_b1$Grouplab,levels = c("Seoul","Beijing"),
                                   labels=c(expression(bold("Seoul")),
                                            expression(bold("Beijing"))))

cor_vk_d_ec_sel_b1$varlab=cor_vk_d_ec_sel_b1$variable
cor_vk_d_ec_sel_b1$varlab=factor(cor_vk_d_ec_sel_b1$varlab,levels = c("POA","SOA"),
                                 labels = c(expression(bold("POC")),
                                            expression(bold("SOC"))))

ggplot(cor_vk_d_ec_sel_b1, aes(x=O.C, y=H.C,col=value))+
  geom_point(size=2.7)+
  facet_grid(Grouplab~varlab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
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
  ggsave(filename("vk_cor_b_ec_special"),height = 20, width = 40, units = "cm", dpi = 300)

###cor with SOA related composition=====
wsoc=fread("Datafile/WOSCbb_S&B.csv")

FT_soa=FT_soa %>% inner_join(wsoc[,c("Sample","WSOCbb","WSOCnbb","POA","SOA")], by="Sample")
FT_soa
table(FT_soa$SOAtype)

FT_soa_sel=subset(FT_soa,FT_soa$SOAtype=="SOA")

soa_inty_all=as.data.table(aggregate(FT_soa$Bromo.Inty, by=list(`Sample`=FT_soa$Sample), FUN=sum)) %>% 
  `colnames<-`(c("Sample","Tot"))

soa_inty=as.data.table(aggregate(FT_soa_sel$Bromo.Inty, by=list(`Sample`=FT_soa_sel$Sample, Comp=FT_soa_sel$Comp), FUN=sum))
soa_inty

soa_inty=soa_inty %>% inner_join(soa_inty_all,by = "Sample")
soa_inty$rel=round(soa_inty$x/soa_inty$Tot*100,1)

tt=aggregate(soa_inty$rel, by=list(`Sample`=soa_inty$Sample), FUN=sum)
tt=tt %>% separate(Sample, c("Group","No"))

ggplot(tt,aes(x=Group, y=x, fill=Group))+
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 0.85), width=0.4)+
  geom_boxplot(position = position_dodge(width = 0.85), outlier.color = NA)+
  geom_blank(aes(x=1.5, y=50))+
  scale_fill_manual(values = c("orangered2", "royalblue2"))+
  scale_y_continuous(name = "Proportion of SOA (%)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.3,0,0.3,0),"cm")),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.3,0.0,0.3),"cm")),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 20, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.8,0.0,0.3),"cm")),
        #strip.text.x = element_text(size = 25, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.6,0,0.4,0),"cm")),
        #strip.text.y = element_text(size = 25, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0.6,0.4,0.2),"cm")),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.text = element_text(size = 16, colour = "black", family = "Arial",hjust = 0.5,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26,hjust = 0.2 ,vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
        legend.key.height = unit(2,"cm"),
        legend.key.width = unit(2,"cm"),
        #legend.position = c(0.15,0.9),
        legend.position = "Null",
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  guides(fill=guide_legend(title = "", override.aes = list(size=1)))+
  ggsave(filename("SOA_comp"),height = 20, width = 20, units = "cm", dpi = 300)


tt_b=subset(tt,tt$Group=="Beijing")
tt_s=subset(tt,tt$Group=="Seoul")

wilcox.test(tt_b$x,tt_s$x,exact = F)

soa_inty

soa_inty=soa_inty %>% separate(Sample, c("Group","No"),sep = "_")
soa_inty

soa_inty$Comp=factor(soa_inty$Comp,levels = c("CHO","CHON","CHOS","CHONS","SOA"))

soa_inty

soa_inty_m=dcast(soa_inty,Group+No~Comp, sum, value.var = "rel")

soa_inty_m$Sample=paste(soa_inty_m$Group,soa_inty_m$No,sep = "_")

wsoc
wsoc_sel=wsoc[, c("Sample","Group","No","OC","PM2.5","SO4","NO3","NH4","O3","CO","SO2","NO")]

wsoc_sel$`NO3/OC`=wsoc_sel$NO3/wsoc_sel$OC
wsoc_sel$`SO4/OC`=wsoc_sel$SO4/wsoc_sel$OC
wsoc_sel$`NO3/SO4`=wsoc_sel$NO3/wsoc_sel$SO4

soa_inty_m$No=as.numeric(soa_inty_m$No)
soa_inty_envi=soa_inty_m %>% inner_join(wsoc_sel, by=c("Sample","Group","No"))
soa_inty_envi

soa_inty_envi_m=melt(soa_inty_envi, id.vars = c("Sample","Group","No","OC","CHO","CHON","CHOS","CHONS"), variable.name = "Envi", value.name = "Conc") %>% 
  melt(id.vars=c("Sample","Group","No","OC","Envi","Conc"),variable.name="Comp")

soa_inty_envi_m

soa_inty_envi_m$No=as.numeric(soa_inty_envi_m$No)
soa_inty_envi_m$mark="no"
soa_inty_envi_m$mark=ifelse(soa_inty_envi_m$No==5,"12/19",soa_inty_envi_m$mark)
soa_inty_envi_m$mark=ifelse(soa_inty_envi_m$No==6,"12/20",soa_inty_envi_m$mark)
soa_inty_envi_m$mark=ifelse(soa_inty_envi_m$No==7,"12/21",soa_inty_envi_m$mark)
soa_inty_envi_m$mark=ifelse(soa_inty_envi_m$No==8,"12/22",soa_inty_envi_m$mark)
soa_inty_envi_m$mark=ifelse(soa_inty_envi_m$No==9,"12/23",soa_inty_envi_m$mark)

soa_inty_envi_m
table(soa_inty_envi_m$mark)

soa_inty_envi_m_s=subset(soa_inty_envi_m,soa_inty_envi_m$Group=="Seoul")
soa_inty_envi_m_s1=soa_inty_envi_m_s %>% filter(Envi%in%c("PM2.5","NH4","NO3","SO4")) %>% droplevels()


ggplot()+
  geom_smooth(data=soa_inty_envi_m_s1,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=soa_inty_envi_m_s1,aes(x=value,y=Conc),col="royalblue2", size=3)+
  geom_point(data=subset(soa_inty_envi_m_s1,soa_inty_envi_m_s1$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(soa_inty_envi_m_s1,soa_inty_envi_m_s1$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envi~Comp, scales = "free",labeller = label_parsed, switch = "y")+
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
  ggsave(filename("compvsEnvi_test"),height = 30, width = 38, units = "cm", dpi = 300)


soa_inty_envi_m_s1=soa_inty_envi_m_s %>% filter(Envi%in%c("PM2.5","NH4","NO3","SO4")) %>% droplevels()
soa_inty_envi_m_s1$Envilab=soa_inty_envi_m_s1$Envi
soa_inty_envi_m_s1$Envilab=factor(soa_inty_envi_m_s1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                                labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                           expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))

soa_inty_envi_m_s1$Complab=soa_inty_envi_m_s1$Comp
soa_inty_envi_m_s1$Complab=factor(soa_inty_envi_m_s1$Complab, levels = c("CHO","CHON","CHOS","CHONS","SOA"),
                                labels = c(expression(bold("CHO (SOA)")),expression(bold("CHON (SOA)")),expression(bold("CHOS (SOA)")),
                                           expression(bold("CHONS (SOA)")),expression(bold("SOA"))))

val_max_s1=as.data.table(aggregate(soa_inty_envi_m_s1$Conc, 
                                   by=list(`Comp`=soa_inty_envi_m_s1$Comp,Envi=soa_inty_envi_m_s1$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))
val_max_s1$Envilab=val_max_s1$Envi
val_max_s1$Envilab=factor(val_max_s1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                          labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))
val_max_s1

soa_inty_envi_m_s1

ggplot()+
  geom_smooth(data=soa_inty_envi_m_s1,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=soa_inty_envi_m_s1,aes(x=value,y=Conc),col="royalblue2", size=3)+
  geom_point(data=soa_inty_envi_m_s1,aes(x=value,y=Conc*1.5),col=NA, size=3)+
  #geom_blank(data=val_max_s1, aes(x=20,y=Conc*1.4))+
  geom_point(data=subset(soa_inty_envi_m_s1,soa_inty_envi_m_s1$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(soa_inty_envi_m_s1,soa_inty_envi_m_s1$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envilab~Complab, scales = "free",labeller = label_parsed, switch = "y")+
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
  ggsave(filename("compvsEnvi_SOA_sul1"),height = 30, width = 30, units = "cm", dpi = 300)


soa_inty_envi_m_b=subset(soa_inty_envi_m,soa_inty_envi_m$Group=="Beijing")
soa_inty_envi_m_b1=soa_inty_envi_m_b %>% filter(Envi%in%c("PM2.5","NH4","NO3","SO4")) %>% droplevels()


soa_inty_envi_m_b1=soa_inty_envi_m_b %>% filter(Envi%in%c("PM2.5","NH4","NO3","SO4")) %>% droplevels()
soa_inty_envi_m_b1$Envilab=soa_inty_envi_m_b1$Envi
soa_inty_envi_m_b1$Envilab=factor(soa_inty_envi_m_b1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                                  labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                             expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                             expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                             expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))

soa_inty_envi_m_b1$Complab=soa_inty_envi_m_b1$Comp
soa_inty_envi_m_b1$Complab=factor(soa_inty_envi_m_b1$Complab, levels = c("CHO","CHON","CHOS","CHONS","SOA"),
                                  labels = c(expression(bold("CHO (SOA)")),expression(bold("CHON (SOA)")),expression(bold("CHOS (SOA)")),
                                             expression(bold("CHONS (SOA)")),expression(bold("SOA"))))

val_max_b1=as.data.table(aggregate(soa_inty_envi_m_b1$Conc, 
                                   by=list(`Comp`=soa_inty_envi_m_b1$Comp,Envi=soa_inty_envi_m_b1$Envi), FUN=max)) %>% 
  `colnames<-`(c("Comp","Envi","Conc"))
val_max_b1$Envilab=val_max_b1$Envi

val_max_b1$Envilab=factor(val_max_b1$Envilab, levels = c("PM2.5","NH4","NO3","SO4"),
                          labels = c(expression(bold("PM"["2.5"]~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NH"["4"]^{"+"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("NO"["3"]^{"-"}~"("*"\u03bcg/"*m^"3"*")")),
                                     expression(bold("SO"["4"]^{"2-"}~"("*"\u03bcg/"*m^"3"*")"))))
val_max_b1

soa_inty_envi_m_b1

ggplot()+
  geom_smooth(data=soa_inty_envi_m_b1,aes(x=value,y=Conc),formula = y~x, method = "lm", se = T, col="black",lty=2)+
  geom_point(data=soa_inty_envi_m_b1,aes(x=value,y=Conc),col="orangered3", size=3)+
  geom_point(data=soa_inty_envi_m_b1,aes(x=value,y=Conc*1.5),col=NA, size=3)+
  #geom_blank(data=val_max_s1, aes(x=20,y=Conc*1.4))+
  geom_point(data=subset(soa_inty_envi_m_b1,soa_inty_envi_m_b1$mark!="no"),aes(x=value,y=Conc),col="black", size=3)+
  geom_text_repel(data=subset(soa_inty_envi_m_b1,soa_inty_envi_m_b1$mark!="no"),aes(x=value,y=Conc, label=mark),
                  min.segment.length = 0.2,box.padding = 0.4,col="black", size=4)+
  #facet_grid(Envilab~Comp, scales = "free",labeller = label_parsed, switch = "y")+
  facet_grid(Envilab~Complab, scales = "free",labeller = label_parsed, switch = "y")+
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
  ggsave(filename("compvsEnvi_SOA_b1"),height = 30, width = 30, units = "cm", dpi = 300)

##correlation test====

cor_temp=soa_inty_envi_m

grp=unique(soa_inty_envi_m$Group)
grp
comp=unique(soa_inty_envi_m$Comp)
envi=unique(soa_inty_envi_m$Envi)


dt=data.table()
for (i in 1:length(grp)) {
  #i=1
  temp=subset(soa_inty_envi_m,soa_inty_envi_m$Group==grp[i])
  
  for (k in 1:length(comp)) {
    #k=2
    tmp=subset(temp,temp$Comp==comp[k])
    
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

fwrite(dt,file="SOA_CompvsEnvi.csv")



###Correlation with levoglucosan=====
ft_levo=fread("Datafile/FT_formula_ob3_rel.csv")
ft_levo[,1:23]

cor_temp=ft_levo

cor_temp[,1:23]

grp=unique(cor_temp$Group)

df=data.frame()
for (i in 1:length(grp)) {
  #i=1
  temp=as.data.frame(subset(cor_temp,cor_temp$Group==grp[i]))
  for (j in 4:6) {
    #j=4
    for (k in 7:dim(temp)[2]) {
      #k=8  
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


df

table(df$Envi)

cor_levo=subset(df,df$Envi=="Levo")
cor_levo


CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(cor_levo$Formula)

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

cor_levo=cor_levo %>% inner_join(fm)
cor_levo$O.C=cor_levo$O/cor_levo$C
cor_levo$H.C=cor_levo$H/cor_levo$C


cor_levo_sel=subset(cor_levo,cor_levo$p<0.05)
cor_levo_sel=subset(cor_levo_sel,abs(cor_levo_sel$rho)>0.001)
cor_levo_sel=cor_levo_sel[order(cor_levo_sel$rho),]


cor_levo_sel
cor_levo_sel$Grouplab=cor_levo_sel$Group
cor_levo_sel$Grouplab=factor(cor_levo_sel$Grouplab,levels = c("Seoul","Beijing"),
                                labels=c(expression(bold("Seoul")),
                                         expression(bold("Beijing"))))


cor_levo_sel_s=subset(cor_levo_sel,cor_levo_sel$Group=="Seoul")
cor_levo_sel_s1=subset(cor_levo_sel_s,abs(cor_levo_sel_s$rho)>0.4)
cor_levo_sel_s2=subset(cor_levo_sel_s1,cor_levo_sel_s1$p<0.05)

ggplot(cor_levo_sel_s2, aes(x=O.C, y=H.C,col=rho))+
  geom_point(size=2.7)+
  facet_grid(.~Grouplab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_gradientn(colors=(topo.colors(30)[4:28]))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
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
        legend.title = element_text(size = 26,hjust = 0.15, vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
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
  ggsave(filename("vk_levo_cor_s"),height = 20, width = 25, units = "cm", dpi = 700)



cor_levo_sel_b=subset(cor_levo_sel,cor_levo_sel$Group=="Beijing")
cor_levo_sel_b1=subset(cor_levo_sel_b,abs(cor_levo_sel_b$rho)>0.4)
cor_levo_sel_b2=subset(cor_levo_sel_b1,cor_levo_sel_b1$p<0.001)

ggplot(cor_levo_sel_b2, aes(x=O.C, y=H.C,col=rho))+
  geom_point(size=2.7)+
  facet_grid(.~Grouplab, labeller = label_parsed )+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_gradientn(colors=(topo.colors(30)[4:26]))+
  #scale_color_gradientn(colors=(matlab.like(40)[6:34]),breaks = c(-0.5, 0, 0.5), limits = c(-0.75,0.93))+
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
        legend.title = element_text(size = 26,hjust = 0.15, vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
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
  ggsave(filename("vk_levo_cor_b"),height = 20, width = 25, units = "cm", dpi = 700)
