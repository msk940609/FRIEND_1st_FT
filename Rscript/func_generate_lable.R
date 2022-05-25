generate_label <- function(data, flev, value.var,cld_info, offset=1.2){
  
  
  
  
  
  data_raw=as.data.frame(data)
  
  plot.labels <- names(cld_info[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  iqr.df <- ddply(data_raw, flev, function (x) ((1.5*IQR(x[,value.var])+quantile(x[,value.var],0.75)))) %>% 
    `colnames<-`(c(flev,"iqr"))
  data_raw=data_raw %>% inner_join(iqr.df, by=flev)
  d2=subset(data_raw,data_raw[,value.var]<data_raw$iqr)
  
  boxplot.df <- ddply(d2, flev, function (x) (max(fivenum(x[,value.var])))+offset)
  #boxplot.df <- ddply(d2, flev, function (x) (2*max(fivenum(x[,value.var]))-min(fivenum(x[,value.var])))*offset)
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = ans.mcV[['Letters']],
                            stringsAsFactors = FALSE)
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = F)
  
  return(labels.df)
}

#data=subset(ft_inty,ft_inty$Comp=="CHONS")

#dm=as.data.table(aggregate(data$rel, by=list(`variable`=data$Group), FUN=median))
#dm=dm[order(dm$x, decreasing = T),]

#data$Group=factor(data$Group,levels = dm$variable)
#ans <- posthoc.kruskal.dunn.test(rel ~ Group,data=data ,p.adjust="fdr")
#ans.p <- get.pvalues(ans)
#ans.mcV <- multcompLetters(ans.p, threshold=0.05)

#tt=generate_label(data=data, flev="Group",value.var = "rel", offset = 1, cld_info = ans.mcV )
#tt


