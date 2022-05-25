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

ggplot(cor_vk_d_nbb_sel_s1, aes(x=O.C, y=H.C,col=value))+
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


cor_vk_d_nbb_sel_p=subset(cor_vk_d_nbb_sel,cor_vk_d_nbb_sel$value>0)
cor_vk_d_nbb_sel_p

fwrite(cor_vk_d_nbb_sel_p,file = "nbb_formula_list.csv")

##
fm_list=fread("Datafile/nbb_formula_list_fin.csv")
fm_list

CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(fm_list$Formula)

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


fm_list=fm_list %>% inner_join(fm)

fm_list
fm_list$DBE=fm_list$C-fm_list$H/2+fm_list$N/2+1
fm_list$`O/C`=fm_list$`O`/fm_list$`C`
fm_list$`H/C`=fm_list$`H`/fm_list$`C`
fm_list$`S/C`=fm_list$`S`/fm_list$`C`
fm_list$`N/C`=fm_list$`N`/fm_list$`C`

fm_list$CAI=fm_list$`C`-fm_list$`N`-fm_list$`O`-fm_list$`S`
fm_list$DBEAI=1+fm_list$`C`-fm_list$`O`-fm_list$`S`-fm_list$`H`/2-fm_list$`N`/2
fm_list$AI=ifelse(fm_list$CAI<=0,0,ifelse(fm_list$DBEAI<0,0,fm_list$DBEAI/fm_list$CAI))
fm_list$NOSC=round(2*fm_list$`O/C`-fm_list$`H/C`+3*fm_list$`N/C`+2*fm_list$`S/C`,4)

fm_list

library(party)
set.seed(123)
library(tidyverse)
library(caret)
library(rpart)

?train
model <- train(
  Type ~., data = fm_list[,c("Type","O/C","H/C")], method = "ctree2",
  trControl = trainControl("cv", number = 10),
  tuneGrid = expand.grid(maxdepth = 3, mincriterion = 0.95 )
)
?plot
plot(model$finalModel)


install.packages("partykit", repos = "http://R-Forge.R-project.org")

library(rpart)
install.packages("rpart.plot")
library(rpart.plot)

fit <- rpart(Type~., data = fm_list[,c("Type","O/C","H/C")], method = 'class')

rpart.plot(fit, extra = 106)

fm_list$Type=as.factor(fm_list$Type)
output.tree <- ctree(
  Type ~ ., 
  data = fm_list[,c("Type","O/C","H/C")])
plot(output.tree)

library("partykit")

ct <- ctree(Species ~ ., data = iris[,c("Species","Sepal.Length","Sepal.Width")])
plot(ct, tp_args = list(text = TRUE))

ft <- ctree(Type ~ ., data =  fm_list[,c("Type","O/C","H/C","DBE")])
plot(ft, tp_args = list(text = TRUE))

ft <- ctree(Type ~ ., data =  fm_list[,c("Type","O/C","H/C")],
              control = ctree_control(mincriterion = 0.999,
                                      minsplit = 20,
                                      minbucket = 10
              ))

plot(ft, tp_args = list(text = TRUE))




cor_vk_d_nbb_sel_p
fm_s=subset(cor_vk_d_nbb_sel_p,cor_vk_d_nbb_sel_p$Group=="Seoul")

CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(fm_s$Formula)

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


fm_s=fm_s %>% inner_join(fm)

fm_s
fm_s$DBE=fm_s$C-fm_s$H/2+fm_s$N/2+1
fm_s$`O/C`=fm_s$`O`/fm_s$`C`
fm_s$`H/C`=fm_s$`H`/fm_s$`C`
fm_s$`S/C`=fm_s$`S`/fm_s$`C`
fm_s$`N/C`=fm_s$`N`/fm_s$`C`

fm_s$CAI=fm_s$`C`-fm_s$`N`-fm_s$`O`-fm_s$`S`
fm_s$DBEAI=1+fm_s$`C`-fm_s$`O`-fm_s$`S`-fm_s$`H`/2-fm_s$`N`/2
fm_s$AI=ifelse(fm_s$CAI<=0,0,ifelse(fm_s$DBEAI<0,0,fm_s$DBEAI/fm_s$CAI))
fm_s$NOSC=round(2*fm_s$`O/C`-fm_s$`H/C`+3*fm_s$`N/C`+2*fm_s$`S/C`,4)

fm_s
fm_s$Type=as.factor(fm_s$Type)

ft_s <- ctree(Type ~ ., data =  fm_s[,c("Type","O/C","H/C")],
            control = ctree_control(mincriterion = 0.999,
                                    minsplit = 20,
                                    minbucket = 10
            ))

ft_s <- ctree(Type ~ ., data =  fm_s[,c("Type","O/C","H/C")],
              control = ctree_control(mincriterion = 0.999,
                                      minsplit = 20,
                                      minbucket = 5,
                                      minprob = 0
              ))

plot(ft_s, tp_args = list(text = TRUE))


##beijng decision tree===
fm_b=subset(cor_vk_d_nbb_sel_p,cor_vk_d_nbb_sel_p$Group=="Beijing")

CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(fm_b$Formula)

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


fm_b=fm_b %>% inner_join(fm)

fm_b
fm_b$DBE=fm_b$C-fm_b$H/2+fm_b$N/2+1
fm_b$`O/C`=fm_b$`O`/fm_b$`C`
fm_b$`H/C`=fm_b$`H`/fm_b$`C`
fm_b$`S/C`=fm_b$`S`/fm_b$`C`
fm_b$`N/C`=fm_b$`N`/fm_b$`C`

fm_b$CAI=fm_b$`C`-fm_b$`N`-fm_b$`O`-fm_b$`S`
fm_b$DBEAI=1+fm_b$`C`-fm_b$`O`-fm_b$`S`-fm_b$`H`/2-fm_b$`N`/2
fm_b$AI=ifelse(fm_b$CAI<=0,0,ifelse(fm_b$DBEAI<0,0,fm_b$DBEAI/fm_b$CAI))
fm_b$NOSC=round(2*fm_b$`O/C`-fm_b$`H/C`+3*fm_b$`N/C`+2*fm_b$`S/C`,4)

fm_b
fm_b$Type=as.factor(fm_b$Type)

ft_b <- ctree(Type ~ ., data =  fm_b[,c("Type","O/C","H/C")])
plot(ft_b, tp_args = list(text = TRUE))






