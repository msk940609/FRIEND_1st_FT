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
library(ggtext)
library(lemon)

getOption("digits")
options("digits" = 15)

source("Rscript/func_filename.R")
source("Rscript/func_generate_lable.R")
source("Rscript/func_upper_fence_label.R")

ft_merge=fread("Datafile/FRIENDs_1st_FT.csv")
table(ft_merge$Group)
ft_merge$Group=factor(ft_merge$Group, levels =c("UL","B","SS","SUL","NT"),
                      labels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))
table(ft_merge$Group)
insert_minor_1 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

insert_minor_2 <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", 1) ) ) )
labs[1:(length(labs)-n_minor)]}

fancy_scientific <- function(l) {
  
  l <- format(l, scientific = TRUE)
  
  ifelse(l=="0e+00", "0", 
         
         l <- parse(text=gsub("(.*)e(\\+?)(\\-?[0-9]+)", "\\1 * 10^(\\3)", l))
  )
}
###

ft_merge=fread("Datafile/FRIENDs_1st_FT.csv")
table(ft_merge$Group)
ft_merge$Group=factor(ft_merge$Group, levels =c("UL","B","SS","SUL","NT"),
                      labels = c("Ulaanbaatar","Beijing","Seosan","Seoul","Noto"))
table(ft_merge$Group)
ft_merge$Freq=1

ft_merge$Sample=paste(ft_merge$Group,ft_merge$No,sep = "_")
ft_merge

mol_rich=as.data.table(aggregate(ft_merge$Freq,by=list(Sample=ft_merge$Sample),sum))
mol_rich


##reconstructed m/z distribution=====
ft_mz=ft_merge[,c("Group","Formula","Calc.m.z","Mono.Inty")]

ft_mz=melt(ft_merge[,c("Group","Formula","Calc.m.z","Mono.Inty")], id.vars = c("Group","Formula","Calc.m.z")) %>% 
  dcast(Group~`Calc.m.z`, mean) %>% 
  melt(id.vars=c("Group"), na.rm = T)

ft_mz
ft_mz$variable2=as.numeric(as.character(ft_mz$variable))
ft_mz$value=as.numeric(ft_mz$value)

max(ft_mz$variable2)

ft_mz_tot=as.data.table(aggregate(ft_mz$value, by=list(Group=ft_mz$Group),FUN=sum))
ft_mz=ft_mz %>% inner_join(ft_mz_tot)
ft_mz$rel=ft_mz$value/ft_mz$x*100

ggplot(ft_mz)+
  geom_segment(aes(x=variable2,xend=variable2,y=0,yend=rel))+
  facet_wrap(.~Group, scales = "free_y",ncol = 1)+ 
  scale_x_continuous(name = "", breaks = seq(0,1000,50))+
  scale_y_continuous(name = "")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=22, colour = "black",face = "bold",family = "Arial",
                                  margin = unit(c(0.2,0.4,0.2,0.2),"cm")),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,-0.2,0.1,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.18,0.95),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  coord_cartesian(ylim = c(0,1.5))+
  ggsave(filename("mz distribution"), height = 10, width = 40, units = "cm", dpi = 300)

ft_mz_ul=subset(ft_mz,ft_mz$Group=="Ulaanbaatar")

ggplot()+
  geom_segment(data=ft_mz_ul,aes(x=variable2,xend=variable2,y=0,yend=value), size=0.3)+
  facet_wrap(.~Group, scales = "free_y",ncol = 1)+ 
  scale_x_continuous(name = "", breaks = seq(200,950,100),expand = c(0.02,0.02))+
  scale_y_continuous(name = "",label = fancy_scientific,expand = c(0.02,0.02))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        #axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_markdown(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=22, colour = "black",face = "bold",family = "Arial",
                                  margin = unit(c(0.2,0.4,0.2,0.2),"cm")),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,-0.2,0.1,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.18,0.95),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  #coord_cartesian(ylim = c(0,784000000))+
  ggsave(filename("mz distribution_ul"), height = 8, width = 40, units = "cm", dpi = 300)


ft_mz_bj=subset(ft_mz,ft_mz$Group=="Beijing")

ggplot()+
  geom_segment(data=ft_mz_bj,aes(x=variable2,xend=variable2,y=0,yend=value),size=0.3)+
  facet_wrap(.~Group, scales = "free_y",ncol = 1)+ 
  scale_x_continuous(name = "", breaks = seq(200,950,100),expand = c(0.02,0.02))+
  scale_y_continuous(name = "",label = fancy_scientific,expand = c(0.02,0.02))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        #axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_markdown(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=22, colour = "black",face = "bold",family = "Arial",
                                  margin = unit(c(0.2,0.4,0.2,0.2),"cm")),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,-0.2,0.1,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.18,0.95),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  #coord_cartesian(ylim = c(0,1174000000))+
  ggsave(filename("mz distribution_bj"), height = 8, width = 40, units = "cm", dpi = 300)


ft_mz_ss=subset(ft_mz,ft_mz$Group=="Seosan")

ggplot()+
  geom_segment(data=ft_mz_ss,aes(x=variable2,xend=variable2,y=0,yend=value),size=0.3)+
  facet_wrap(.~Group, scales = "free_y",ncol = 1)+ 
  scale_x_continuous(name = "", breaks = seq(200,950,100),expand = c(0.02,0.02))+
  scale_y_continuous(name = "",label = fancy_scientific,expand = c(0.02,0.02))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        #axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_markdown(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=22, colour = "black",face = "bold",family = "Arial",
                                  margin = unit(c(0.2,0.4,0.2,0.2),"cm")),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,-0.2,0.1,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.18,0.95),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  #coord_cartesian(ylim = c(0,2941000000))+
  ggsave(filename("mz distribution_ss"), height = 8, width = 40, units = "cm", dpi = 300)

ft_mz_sul=subset(ft_mz,ft_mz$Group=="Seoul")

ggplot()+
  geom_segment(data=ft_mz_sul,aes(x=variable2,xend=variable2,y=0,yend=value),size=0.3)+
  facet_wrap(.~Group, scales = "free_y",ncol = 1)+ 
  scale_x_continuous(name = "", breaks = seq(200,950,100),expand = c(0.02,0.02))+
  scale_y_continuous(name = "",label = fancy_scientific,expand = c(0.02,0.02))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        #axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_markdown(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=22, colour = "black",face = "bold",family = "Arial",
                                  margin = unit(c(0.2,0.4,0.2,0.2),"cm")),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,-0.2,0.1,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.18,0.95),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  #coord_cartesian(ylim = c(0,1470000000))+
  ggsave(filename("mz distribution_sul"), height = 8, width = 40, units = "cm", dpi = 300)

ft_mz_nt=subset(ft_mz,ft_mz$Group=="Noto")

ggplot()+
  geom_segment(data=ft_mz_nt,aes(x=variable2,xend=variable2,y=0,yend=value),size=0.3)+
  facet_wrap(.~Group, scales = "free_y",ncol = 1)+ 
  scale_x_continuous(name = "", breaks = seq(200,950,100),expand = c(0.02,0.02))+
  scale_y_continuous(name = "",label = fancy_scientific,expand = c(0.02,0.02))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        #axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_markdown(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 0, colour = "black",face = "bold", family = "Arial"),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size=22, colour = "black",face = "bold",family = "Arial",
                                  margin = unit(c(0.2,0.4,0.2,0.2),"cm")),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0.1,-0.2,0.1,-0.2),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        legend.position = c(0.18,0.95),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  #coord_cartesian(ylim = c(0,588000000))+
  ggsave(filename("mz distribution_NT"), height = 8, width = 40, units = "cm", dpi = 300)

##reconstructed vk plot distribution=====
ft_merge

ft_merge=ft_merge %>% mutate(Comp = gsub('[[:digit:]]', '', Formula))

ft_merge

ft_merge$Comp=ifelse(ft_merge$O.==0,"Remainders",ft_merge$Comp)

ft_vk=ft_merge[,c("Group","Formula","O.C","H.C","Bromo.Inty")]


ft_vk_m=dcast(ft_vk,Group~Formula, mean) %>% 
  melt(id.vars=c("Group"), na.rm = T)

ft_vk_m

CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(ft_vk_m$variable)
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

numericalFormula
ft_vk_m

fmlist=as.data.table(numericalFormula) %>% `colnames<-`(c("C","H","N","O","S"))
fmlist$Formula=molecularFormula

ft_vk_m0=ft_vk_m

ft_vk_m=ft_vk_m %>% `colnames<-`(c("Group","Formula","value"))

ft_vk_m=ft_vk_m %>% inner_join(fmlist)

ft_vk_m$`O/C`=ft_vk_m$O/ft_vk_m$C
ft_vk_m$`H/C`=ft_vk_m$H/ft_vk_m$C

ft_vk_m
ft_vk_m=ft_vk_m %>% mutate(Comp = gsub('[[:digit:]]', '', Formula))

ft_vk_m$Comp=ifelse(ft_vk_m$`O`==0,"Remainders",ft_vk_m$Comp)
ft_vk_m$Comp=factor(ft_vk_m$Comp,levels = c("CHO","CHNO","CHOS","CHNOS","Remainders"),
                    labels = c("CHO","CHON","CHOS","CHONS","Remainders"))

table(ft_vk_m$Comp)

ggplot(ft_vk_m, aes(x=`O/C`, y=`H/C`,col=Comp))+
  geom_abline(slope = -1, intercept = 1.1,lty=2,size=1)+
  geom_abline(slope = -0.76, intercept = 0.75, lty=2,size=1)+
  geom_hline(yintercept = 1.5,lty=2, size=1)+
  geom_hline(yintercept = 2.0,lty=2,size=1)+
  geom_point(size=2, alpha=0.9)+
  facet_rep_wrap(.~Group, repeat.tick.labels = T, ncol=5, scales = "free")+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_manual(values = c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF","grey50"))+
  #scale_fill_manual(values = rev(c("#F2C9C3","#B8CBB5","#9DB9D1","#F9D38A","#958CC5")))+
  #scale_color_manual(values = c("mediumpurple2","olivedrab4","wheat3","steelblue3","rosybrown3"))+
  #scale_color_manual(values = rev(c("#FDA9A3","#B8CBB5","#9DB9D1","#F9D38A","#958CC5")))+
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
        #legend.position = "bottom",
        #legend.direction = "horizontal",
        legend.position = "NULL",
        legend.box = "horizontal",
        legend.background = element_blank(),
  )+
  #guides(col=guide_legend(title = "",ncol = 5,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF","grey50"),alpha=1)), 
  #       size="none")+
  ggsave(filename("vk_class"),height = 20, width = 100, units = "cm", dpi = 300)

ft_vk_m

ggplot(ft_vk_m, aes(x=`O/C`, y=`H/C`,col=Comp))+
  geom_abline(slope = -1, intercept = 1.1,lty=2,size=1)+
  geom_abline(slope = -0.76, intercept = 0.75, lty=2,size=1)+
  geom_hline(yintercept = 1.5,lty=2, size=1)+
  geom_hline(yintercept = 2.0,lty=2,size=1)+
  geom_point(aes(size=value), alpha=0.4)+
  scale_size_continuous(range=c(1.5,30))+
  facet_rep_wrap(.~Group, repeat.tick.labels = T, ncol=5, scales = "free")+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF","grey50"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        strip.text = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        #legend.position = c(0.18,0.14),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = "",title.position = "left",
                          ncol = 5,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF","grey50"),alpha=1)), size="none")+
  ggsave(filename("vk_all"),height = 20, width = 75, units = "cm", dpi = 700)


##observation filtering
ft_merge$Freq=1

fm_obs=ft_merge[,c("Group","No","Formula","Freq")] %>% 
  dcast(Group~Formula,value.var = "Freq",sum) %>% 
  melt(id.vars=c("Group")) %>% `colnames<-`(c("Group","Formula","cnt"))

fm_obs

ft_vk_m_filter
table(ft_vk_m_filter$Group)

ft_vk_m_filter=ft_vk_m %>% inner_join(fm_obs)

ft_vk_m_filter_sel=ft_vk_m_filter %>% filter(cnt>4)
table(ft_vk_m_filter_sel$Group)

ggplot(ft_vk_m_filter_sel, aes(x=`O/C`, y=`H/C`,col=Comp))+
  geom_abline(slope = -1, intercept = 1.1,lty=2,size=1)+
  geom_abline(slope = -0.76, intercept = 0.75, lty=2,size=1)+
  geom_hline(yintercept = 1.5,lty=2, size=1)+
  geom_hline(yintercept = 2.0,lty=2,size=1)+
  geom_point(aes(size=value), alpha=0.4)+
  scale_size_continuous(range=c(1.5,30))+
  facet_rep_wrap(.~Group, repeat.tick.labels = T, ncol=5, scales = "free")+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF","grey50"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0,0.0,0),"cm")),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.4,0.0,0.0),"cm")),
        strip.text = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        #legend.position = c(0.18,0.14),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = "",title.position = "left",
                          ncol = 5,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF","grey50"),alpha=1)), size="none")+
  ggsave(filename("vk_all_obs6"),height = 20, width = 75, units = "cm", dpi = 700)

fm_chp=ft_merge[,c("Formula","AI","DBE","O.C","H.C")]
fm_chp=fm_chp %>% `colnames<-`(c("Formula","AI","DBE","O/C","H/C"))
fm_chp=unique(fm_chp)

ft_vk_m_filter_sel

ft_m=ft_vk_m_filter_sel[,-c("O/C","H/C")] %>% left_join(fm_chp)

ft_m

ft_data=ft_m

ft_data$PCA=ifelse(ft_data$AI>0.66, "Combustion-drived PCA","")

ft_data$aroma=ifelse(ft_data$AI>0.5,ifelse(ft_data$AI<=0.66,"Soil-derived polyphenols and PCAs with aliphatic chains",""),"")

ft_data$humic=ifelse(ft_data$AI<=0.5,
                     ifelse(ft_data$`H/C`<1.5,"Soil derived humics and highly unsaturated compounds",""),"")

ft_data$aliphatic=ifelse(ft_data$`H/C`<2.0&ft_data$`H/C`>=1.5,
                         ifelse(ft_data$AI<=0.5,"Unsaturated aliphatic compounds",""),"")

ft_data$fatty=ifelse(ft_data$`H/C`>=2.0,
                     ifelse(ft_data$AI<=0.5,"Saturated fatty and sulfonic acid,carbohyrates",""),"")

ft_data=ft_data %>%  unite("Molecularclass",c("PCA","aroma","humic","aliphatic","fatty"), sep = "")
ft_data$Molecularclass=ifelse(ft_data$Molecularclass=="","Unassigned", ft_data$Molecularclass) ###colum of molecular class
ft_data

ft_data$Molecularclass=factor(ft_data$Molecularclass,levels = rev(c("Saturated fatty and sulfonic acid,carbohyrates","Unsaturated aliphatic compounds",
                                                                      "Soil derived humics and highly unsaturated compounds",
                                                                      "Soil-derived polyphenols and PCAs with aliphatic chains",
                                                                      "Combustion-drived PCA")))

ft_data
ft_data$Freq=1


ft_data

ft_data$Class2=factor(ft_data$Molecularclass,levels = rev(c("Saturated fatty and sulfonic acid,carbohyrates","Unsaturated aliphatic compounds",
                                                   "Soil derived humics and highly unsaturated compounds",
                                                   "Soil-derived polyphenols and PCAs with aliphatic chains",
                                                   "Combustion-drived PCA")),
                      labels=c("Combustion-derived PCA","Polyphenols and PCAs with aliphatic chains",
                               "Humics and highly unsaturated compounds",
                               "Unsaturated aliphatic compounds",
                               "Saturated fatty and sulfonic acid, carbohyrates"
                      ))

ft_tt=ft_data[,c("Formula","O/C","H/C","Class2")]

ft_tt=unique(ft_tt)

ggplot(ft_tt, aes(x=`O/C`, y=`H/C`,col=Class2))+
  geom_abline(slope = -1, intercept = 1.1,lty=2,size=1)+
  geom_abline(slope = -0.76, intercept = 0.75, lty=2,size=1)+
  geom_hline(yintercept = 1.5,lty=2, size=1)+
  geom_hline(yintercept = 2.0,lty=2,size=1)+
  geom_point(size=1, alpha=0.4)+
  scale_size_continuous(range=c(1.5,30))+
  scale_x_continuous(name = "O/C",expand = c(0.01,0.01),limits = c(-0.01,1.0), breaks = round(seq(0,2.1,0.1),1),
                     labels =insert_minor_2(c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0","2.2"),2))+
  scale_y_continuous(name = "H/C",expand = c(0.01,0.01), limits = c(0.0,2.2),breaks = round(seq(0.0,2.5,0.1),1),
                     labels =insert_minor_1(c("0","0.5","1.0","1.5","2.0","2.5"),4))+
  scale_color_manual(values = rev(c("#586957","#FDB2AD","#4F94CD","#FAD692","#9D95CA")))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 2),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"),
        axis.ticks = element_line(size = 1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.4,0,0.0,0),"cm")),
        axis.title.y = element_text(size = 22, colour = "black",face = "bold", family = "Arial",margin = unit(c(0.0,0.4,0.0,0.0),"cm")),
        strip.text = element_text(size = 40, colour = "black",face = "bold",family = "Arial",margin = unit(c(0.2,0,0.4,0),"cm")),
        strip.background = element_blank(),
        legend.text = element_text(size = 22, colour = "black", family = "Arial", vjust = 0.6,margin = unit(c(0,0,0,0),"cm")),
        legend.spacing = unit(0.0,"cm"),
        legend.title = element_text(size = 26, colour = "black", family = "Arial",face = "bold"),
        #legend.position = c(0.18,0.14),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  guides(col=guide_legend(title = "",title.position = "left",
                          ncol = 5,override.aes = list(shape=21,size=7, fill=c("#BC3C29FF","#EFC000FF","#008B45FF","#7876B1FF","grey50"),alpha=1)), size="none")+
  ggsave(filename("vk_ex"),height = 20, width = 20, units = "cm", dpi = 700)


pie_all=as.data.table(aggregate(ft_data$value, 
                               by=list(`Group`=ft_data$Group,
                                       Class=ft_data$Molecularclass), FUN=sum))

pie_all_tot=as.data.table(aggregate(ft_data$value, 
                                   by=list(`Group`=ft_data$Group), FUN=sum)) %>% `colnames<-`(c("Group","tot"))

pie_all=pie_all %>% inner_join(pie_all_tot)

pie_all$rel=pie_all$x/pie_all$tot*100
pie_all

pie_all$Class2=factor(pie_all$Class,levels = rev(c("Saturated fatty and sulfonic acid,carbohyrates","Unsaturated aliphatic compounds",
                                                "Soil derived humics and highly unsaturated compounds",
                                                "Soil-derived polyphenols and PCAs with aliphatic chains",
                                                "Combustion-drived PCA")),
                     labels=c("Combustion-derived PCA","Polyphenols and PCAs with aliphatic chains",
                              "Humics and highly unsaturated compounds",
                              "Unsaturated aliphatic compounds",
                              "Saturated fatty and sulfonic acid, carbohyrates"
                              ))




#pie_inty=pie_all[,c(1,2,5)]
pie_frq=pie_all[,c(1,2,5)]

pie_frq$inty=round(pie_inty$rel,1)

pie_frq


ggplot()+
  geom_bar(data=pie_all,aes(x="",y=x, fill=Class2), col="black",stat="identity", size=0.5,
           position = position_fill(reverse = T))+
  coord_polar("y")+
  scale_fill_manual(values = rev(c("#586957","#FDB2AD","#4F94CD","#FAD692","#9D95CA")))+
  #scale_fill_manual(values =  c("#BC3C29FF","#EFC000FF","#008B45FF","#5F559BFF","grey50"))+
  facet_wrap(.~Group, ncol = 5)+
  theme_bw()+
  theme(
    panel.background = element_rect(fill = "transparent", color = NA, size = 0), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA, size=0), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    strip.text = element_text(size = 20, family = "Arial",face = "bold"),
    legend.title = element_text(size = 16, vjust = 0.75,colour = "black", family = "Arial",face = "bold"),
    legend.text = element_text(size = 16, family = "Arial",face = "bold",colour = "black",hjust = 0.0,margin = unit(c(0,0.1,0,0.1),"cm")),
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent", color = NA, size=0), # get rid of legend panel bg
    legend.key.height = unit(0.9,"cm"),
    legend.key.width = unit(0.9,"cm"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    legend.position = "bottom"
  )+
  guides(fill=guide_legend(reverse = F, ncol = 3, byrow=T, title = NULL, title.position = "top"))+
  ggsave("220506pie_all_trfu_legend.png",height = 12, width = 50, units = "cm", dpi = 300)



###BB&NBB piechart=====

ft_m

ft_m_soa=ft_m %>% left_join(ft_cor_pos[,c("Group","Formula","type")])
ft_m_soa

ft_m_soa_bb=subset(ft_m_soa,ft_m_soa$type=="WSOCbb")
ft_m_soa_nbb=subset(ft_m_soa,ft_m_soa$type=="WSOCnbb")

ft_m_soa_bb
table(ft_m_soa_bb$Group)
ft_m_soa_nbb
table(ft_m_soa_nbb$Group)

ft_m_soa_merge=rbind(ft_m_soa_bb,ft_m_soa_nbb)
ft_m_soa_merge


ft_m_soa_merge$variable2=factor(ft_m_soa_merge$type, level=c("WSOCbb","WSOCnbb"),
                         labels=c("WSOC ","WSOC  "))


table(ft_m_soa_merge$Comp)

ft_m_soa_merge$Comp2=factor(ft_m_soa_merge$Comp, levels = c("CHO","CHON","CHOS","CHONS","Remainders"))

ft_m_soa_merge_b=subset(ft_m_soa_merge,ft_m_soa_merge$Group=="Beijing")

tt=c("#BC3C29FF","#EFC000","#008B45","#5F559B","grey30")

PieDonut_ms(ft_m_soa_merge_b,aes(variable2 ,Comp2), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,showPieName = F,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill="manual",fillcol = c("#9DBCD4","#CB7723"),
            subfill = "manual", subfillcol=rep(tt,2))

png("Figure/B_pie_comp.png", width = 20, height = 20,units = "cm",res = 300, bg = "transparent")
PieDonut_ms(ft_m_soa_merge_b,aes(variable2 ,Comp2), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,showPieName = F,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill="manual",fillcol = c("#9DBCD4","#CB7723"),
            subfill = "manual", subfillcol=rep(tt,2))
dev.off()

ft_m_soa_merge_ul=subset(ft_m_soa_merge,ft_m_soa_merge$Group=="Ulaanbaatar")
png("Figure/UL_pie_comp.png", width = 20, height = 20,units = "cm",res = 300, bg = "transparent")
PieDonut_ms(ft_m_soa_merge_ul,aes(variable2 ,Comp2), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,showPieName = F,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill="manual",fillcol = c("#9DBCD4","#CB7723"),
            subfill = "manual", subfillcol=rep(tt,2))
dev.off()

ft_m_soa_merge_sul=subset(ft_m_soa_merge,ft_m_soa_merge$Group=="Seoul")
png("Figure/Sul_pie_comp.png", width = 20, height = 20,units = "cm",res = 300, bg = "transparent")
PieDonut_ms(ft_m_soa_merge_sul,aes(variable2 ,Comp2), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,showPieName = F,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill="manual",fillcol = c("#9DBCD4","#CB7723"),
            subfill = "manual", subfillcol=rep(tt,2))
dev.off()


ft_m_soa_merge_ss=subset(ft_m_soa_merge,ft_m_soa_merge$Group=="Seosan")
png("Figure/ss_pie_comp.png", width = 20, height = 20,units = "cm",res = 300, bg = "transparent")
PieDonut_ms(ft_m_soa_merge_ss,aes(variable2 ,Comp2), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,showPieName = F,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill="manual",fillcol = c("#9DBCD4","#CB7723"),
            subfill = "manual", subfillcol=rep(tt,2))
dev.off()

ft_m_soa_merge_nt=subset(ft_m_soa_merge,ft_m_soa_merge$Group=="Noto")
png("Figure/nt_pie_comp.png", width = 20, height = 20,units = "cm",res = 300, bg = "transparent")
PieDonut_ms(ft_m_soa_merge_nt,aes(variable2 ,Comp2), pieLabelSize=8,donutLabelSize = 6,
            pieAlpha = 0.8, donutAlpha = 0.8,showPieName = F,
            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.01),
            labelpositionThreshold=0.05,labelposition=2,fill="manual",fillcol = c("#9DBCD4","#CB7723"),
            subfill = "manual", subfillcol=rep(tt,2))
dev.off()



ft_chons=subset(ft_merge,ft_merge$Comp=="CHONS")
ft_chons

table(ft_chons$N.)
table(ft_chons$S.)


ft_chons$nostype=paste0("CHO","N",ft_chons$N.,"S",ft_chons$S.)
table(ft_chons$nostype)


ft_chons$Freq=1
ft_chons=ft_chons %>% left_join(fm_obs)
ft_chons_sel=subset(ft_chons,ft_chons$cnt>4)
ft_chons_sel

ft_chons_fm=ft_chons_sel


chons_table=as.data.table(ft_chons$Fre1,by=list(S))