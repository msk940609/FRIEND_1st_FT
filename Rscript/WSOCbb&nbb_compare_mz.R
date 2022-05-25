


ft_cor_vk_m_pos=subset(ft_cor_vk_m,ft_cor_vk_m$value>0)


fwrite(ft_cor_vk_m_pos, file = "Datafile/cor_pos_list.csv")
ft_cor_vk_m


ft_cor_sel=ft_cor_vk_m_pos[,c("Group","Formula","variable2")]

venn_fm=dcast(fm_obs_sel,Formula~Group,sum)
fm_obs_sel

ft_cor_bb=subset(ft_cor_sel,ft_cor_sel$variable2=="WSOCbb")
ft_cor_nbb=subset(ft_cor_sel,ft_cor_sel$variable2=="WSOCnbb")

ft_cor_bb$Freq=1
venn_fm_bb=dcast(ft_cor_bb,Formula~Group,sum,value.var = "Freq")

ft_cor_nbb$Freq=1
venn_fm_nbb=dcast(ft_cor_nbb,Formula~Group,sum,value.var = "Freq")


venn_fm_bb
ft_cor_bb_ul=subset(ft_cor_bb,ft_cor_bb$Group=="Ulaanbaatar")
ft_cor_bb_bj=subset(ft_cor_bb,ft_cor_bb$Group=="Beijing")
ft_cor_bb_sul=subset(ft_cor_bb,ft_cor_bb$Group=="Seoul")
ft_cor_bb_ss=subset(ft_cor_bb,ft_cor_bb$Group=="Seosan")
ft_cor_bb_nt=subset(ft_cor_bb,ft_cor_bb$Group=="Noto")

fm_u=as.vector(ft_cor_bb_ul$Formula)
fm_b=as.vector(ft_cor_bb_bj$Formula)
fm_sul=as.vector(ft_cor_bb_sul$Formula)
fm_ss=as.vector(ft_cor_bb_ss$Formula)
fm_nt=as.vector(ft_cor_bb_nt$Formula)

tt <- list(
  Ulaanbaatar = fm_u,
  Beijing = fm_b,
  Seoul = fm_sul, 
  Seosan = fm_ss,
  Noto=fm_nt
)

display_venn <- function(x, fn=NULL,...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = fn, ...)
  grid.draw(venn_object)
}

venn.diagram(tt,
             category.names = c("Ulaanbaatar" , "Beijing" , "Seoul", "Seosan","Noto"),
             #fill =  c("#3C5488FF","#E64B35FF","#4DBBD5FF","#00A087FF","#F39B7FFF"),
             fill =  c("white","white","white","white","white"),
             cat.dist = c(0.1, 0.1, 0.1, 0.1,0.1),
             lty = 1,  lwd = 1,
             cat.cex=1,
             #cex=1.2,
             cat.fontfamily= "Arial",main.fontface = "bold",sub.fontface = "bold",
             fontfamily= "Arial",
             filename = "bb_venn-5-dimensions.tiff")


venn_fm_nbb


venn_fm_nbb
ft_cor_nbb_ul=subset(ft_cor_nbb,ft_cor_nbb$Group=="Ulaanbaatar")
ft_cor_nbb_bj=subset(ft_cor_nbb,ft_cor_nbb$Group=="Beijing")
ft_cor_nbb_sul=subset(ft_cor_nbb,ft_cor_nbb$Group=="Seoul")
ft_cor_nbb_ss=subset(ft_cor_nbb,ft_cor_nbb$Group=="Seosan")
ft_cor_nbb_nt=subset(ft_cor_nbb,ft_cor_nbb$Group=="Noto")

fm_u=as.vector(ft_cor_nbb_ul$Formula)
fm_b=as.vector(ft_cor_nbb_bj$Formula)
fm_sul=as.vector(ft_cor_nbb_sul$Formula)
fm_ss=as.vector(ft_cor_nbb_ss$Formula)
fm_nt=as.vector(ft_cor_nbb_nt$Formula)

tt <- list(
  Ulaanbaatar = fm_u,
  Beijing = fm_b,
  Seoul = fm_sul, 
  Seosan = fm_ss,
  Noto=fm_nt
)


venn.diagram(tt,
             category.names = c("Ulaanbaatar" , "Beijing" , "Seoul", "Seosan","Noto"),
             #fill =  c("#3C5488FF","#E64B35FF","#4DnbbD5FF","#00A087FF","#F39B7FFF"),
             fill =  c("white","white","white","white","white"),
             cat.dist = c(0.1, 0.1, 0.1, 0.1,0.1),
             lty = 1,  lwd = 1,
             cat.cex=1,
             #cex=1.2,
             cat.fontfamily= "Arial",main.fontface = "bold",sub.fontface = "bold",
             fontfamily= "Arial",
             filename = "nbb_venn-5-dimensions.tiff")


venn_fm_bb$U=ifelse(venn_fm_bb$Ulaanbaatar>0, "U",NA)
venn_fm_bb$B=ifelse(venn_fm_bb$Beijing>0, "B",NA)
venn_fm_bb$S=ifelse(venn_fm_bb$Seoul>0, "S",NA)
venn_fm_bb$SS=ifelse(venn_fm_bb$Seosan>0, "SS",NA)
venn_fm_bb$N=ifelse(venn_fm_bb$Noto>0, "N",NA)

venn_fm_bb=venn_fm_bb %>% unite( "ven_class",c("U","B","S","SS","N"),sep = "&",na.rm=TRUE)
venn_fm_bb

venn_fm_nbb$U=ifelse(venn_fm_nbb$Ulaanbaatar>0, "U",NA)
venn_fm_nbb$B=ifelse(venn_fm_nbb$Beijing>0, "B",NA)
venn_fm_nbb$S=ifelse(venn_fm_nbb$Seoul>0, "S",NA)
venn_fm_nbb$SS=ifelse(venn_fm_nbb$Seosan>0, "SS",NA)
venn_fm_nbb$N=ifelse(venn_fm_nbb$Noto>0, "N",NA)

venn_fm_nbb=venn_fm_nbb %>% unite( "ven_class",c("U","B","S","SS","N"),sep = "&",na.rm=TRUE)
venn_fm_nbb


table(venn_fm_bb$ven_class)
bb_all=subset(venn_fm_bb,venn_fm_bb$ven_class=="U&B&S&SS&N")
nbb_all=subset(venn_fm_nbb,venn_fm_nbb$ven_class=="U&B&S&SS&N")


venn_fm_bb

CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(venn_fm_bb$Formula)
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

bblist=as.data.table(numericalFormula) %>% `colnames<-`(c("C","H","N","O","S"))
bblist$Formula=molecularFormula
bblist


venn_fm_bb=venn_fm_bb %>% inner_join(bblist)

CHEMICAL_ELEMENTS = c("C","H","N","O","S")
molecularFormula <- unique(venn_fm_nbb$Formula)
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

nbblist=as.data.table(numericalFormula) %>% `colnames<-`(c("C","H","N","O","S"))
nbblist$Formula=molecularFormula
nbblist

venn_fm_nbb=venn_fm_nbb %>% inner_join(nbblist)

fwrite(venn_fm_bb, file = "bb_all.csv")
fwrite(venn_fm_nbb, file = "nbb_all.csv")
ft_merge

mz_list=ft_merge[,c("Formula","Calc.m.z")]


mz_list=unique(mz_list)
ft_cor_bb_mz=ft_cor_bb %>% left_join(mz_list)

ft_cor_nbb_mz=ft_cor_nbb %>% inner_join(mz_list)

ft_cor_mz=rbind(ft_cor_bb_mz,ft_cor_nbb_mz)

ft_mz=melt(ft_merge[,c("Group","Formula","Calc.m.z","Mono.Inty")], id.vars = c("Group","Formula","Calc.m.z")) %>% 
  dcast(Group~`Calc.m.z`, mean) %>% 
  melt(id.vars=c("Group"), na.rm = T) %>% 
  `colnames<-`(c("Group","Calc.m.z","value"))

ft_mz$Calc.m.z=as.numeric(as.character(ft_mz$Calc.m.z))


ft_cor_mz=ft_cor_mz %>% left_join(ft_mz, by = c("Group","Calc.m.z"))
ft_cor_mz$value=as.numeric(as.character(ft_cor_mz$value))

ft_cor_mz$inty=ifelse(ft_cor_mz$variable2=="WSOCbb",ft_cor_mz$value,ft_cor_mz$value*-1)

ft_merge

ft_comp=ft_merge[,c("Formula","Comp")]

ft_cor_mz=ft_cor_mz %>% inner_join(ft_comp)

ggplot(ft_cor_mz)+
  geom_hline(yintercept = 0)+
  geom_segment(aes(x=`Calc.m.z`,xend=`Calc.m.z`,y=0,yend=inty,col=Comp))+
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
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_blank()
  )+
  coord_cartesian(ylim = c(-500000000,500000000))+
  ggsave(filename("mz distribution"), height = 40, width = 40, units = "cm", dpi = 300)


