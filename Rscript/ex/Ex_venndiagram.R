library("ggVennDiagram")
genes <- paste("gene",1:1000,sep="")
genes

x <- list(
  A = sample(genes,300), 
  B = sample(genes,525), 
  C = sample(genes,440),
  D = sample(genes,350)
)

ggVennDiagram(
  x, label_alpha = 0,
  category.names = c("Stage 1","Stage 2","Stage 3", "Stage4")
) +
  ggplot2::scale_fill_gradient(low="blue",high = "yellow")


pah_merge


pah_name=melt(pah_merge[,c("Sample","Name","Freq")], id.vars = c("Sample","Name")) %>% 
  dcast(Sample~Name, sum) %>% 
  melt(id.vars=c("Sample"))

pah_name=pah_name %>% separate("Sample",c("Group","Date","Event"),sep = "_")
pah_name

table(pah_name_sel$Group)

pah_name_sel=subset(pah_name,pah_name$value>0)

pah_name_b=subset(pah_name,pah_name$Group=="Beijing"&pah_name$value>0)
pah_b=as.vector(unique(pah_name_b$variable))
length(pah_b)
#fwrite(as.data.table(pah_b),file = "Datafile/pah_b.csv")

pah_b=fread("Datafile/pah_b.csv", fill = T)
pah_b=as.vector(pah_b$pah_b)

pah_name_sul=subset(pah_name,pah_name$Group=="Seoul"&pah_name$value>0)
pah_sul=as.vector(unique(pah_name_sul$variable))
length(pah_sul)
#fwrite(as.data.table(pah_sul),file = "Datafile/pah_sul.csv")

pah_sul=fread("Datafile/pah_sul.csv", fill = T)
pah_sul=as.vector(pah_sul$pah_sul)

pah_name_ss=subset(pah_name,pah_name$Group=="Seosan"&pah_name$value>0)
pah_ss=as.vector(unique(pah_name_ss$variable))
#fwrite(as.data.table(pah_ss),file = "Datafile/pah_ss.csv")

pah_ss=fread("Datafile/pah_ss.csv", fill = T)
pah_ss=as.vector(pah_ss$pah_ss)

pah_name_ul=subset(pah_name,pah_name$Group=="Ulaanbaatar"&pah_name$value>0)
pah_ul=as.vector(unique(pah_name_ul$variable))
length(pah_ul)
#fwrite(as.data.table(pah_ul),file = "Datafile/pah_ul.csv")

pah_ul=fread("Datafile/pah_ul.csv", fill = T)
pah_ul=as.vector(pah_ul$pah_ul)



x <- list(
  Seoul = pah_sul, 
  Seosan = pah_ss, 
  Beijing = pah_b,
  Ulaanbaatar = pah_ul
)


library(ggVennDiagram)

ggVennDiagram(
  x, label_alpha = 0,
  #category.names = c("Stage 1","Stage 2","Stage 3", "Stage4")
) +
  ggplot2::scale_fill_gradient(low="yellow",high = "royalblue2")



library(eulerr)

wilkinson2012 <-  c(A = 4, B = 6, C = 3, D = 2, E = 7, F = 3,
                    "A&B" = 2, "A&F" = 2, "B&C" = 2, "B&D" = 1,
                    "B&F" = 2, "C&D" = 1, "D&E" = 1, "E&F" = 1,
                    "A&B&F" = 1, "B&C&D" = 1)
fit3 <- euler(wilkinson2012, shape = "ellipse")
plot(fit3)


pah_ven <- c(A=52,B=37,C=150, D=220,
             "A&B"=9,"A&B&C"=10,"A&B&C&D"=62,"A&B&D"=5,
             "A&C"=54,"A&C&D"=78,
             "A&D"=19,
             "B&C"=6,"B&C&D"=10,
             "B&D"=5,"C&D"=86)
pah_fit <- euler(pah_ven, shape = "ellipse", input = "disjoint")
pah_fit <- euler(pah_ven, shape = "circle", input = "disjoint")
#pah_fit <- euler(pah_ven, shape = "ellipse", input = "union")
#pah_fit <- euler(pah_ven, shape = "circle")

plot(pah_fit,labels = list(font = 4))

display_venn(
  x,
  category.names = c("Set 1" , "Set 2 " , "Set 3", "Set 4"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
  # Numbers
  cex = .9,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.055, 0.055, 0.1, 0.1)
)

install.packages("VennDiagram")
library(VennDiagram)
venn.diagram(x, filename = "venn-4-dimensions.png")

display_venn <- function(x, fn=NULL,...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = fn, ...)
  grid.draw(venn_object)
}
set.seed(20190708)
genes <- paste("gene",1:1000,sep="")
x <- list(
  A = sample(genes,300), 
  B = sample(genes,525), 
  C = sample(genes,440),
  D = sample(genes,350)
)

display_venn(x)




tt <- list(
  Ulaanbaatar = pah_ul,
  Beijing = pah_b,
  Seoul = pah_sul, 
  Seosan = pah_ss 
)


display_venn(
  tt,
  #fn = "ven_test.tiff",
  category.names = c("Ulaanbaatar" , "Beijing" , "Seoul", "Seosan"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c( "#3C5488FF","#E64B35FF","#4DBBD5FF", "#00A087FF"),
  # Numbers
  cex = .9,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.055, 0.055, 0.1, 0.1)
)
?venn.diagram

venn.diagram(tt,
             category.names = c("Ulaanbaatar" , "Beijing" , "Seoul", "Seosan"),
             fill = c( "#3C5488FF","#E64B35FF","#4DBBD5FF", "#00A087FF"),
             lty = 'blank',  lwd = 2,
             filename = "san.tiff")


venn.diagram(tt,
             category.names = c("Ulaanbaatar" , "Beijing" , "Seoul", "Seosan"),
             fill = c( "#3C5488FF","#E64B35FF","#4DBBD5FF", "#00A087FF"),
             lty = 'blank',  lwd = 2,
             cat.fontfamily= "Arial",
             fontfamily= "Arial",
             filename = "venn-4-dimensions.tiff")

  
