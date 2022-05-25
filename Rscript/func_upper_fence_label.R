generate_label_df <- function(data, flev, value.var, offset=1.2, pair.method="tukey"){
  #data=subset(ft_inty,ft_inty$Comp=="CHONS")
  #flev = "Group"
  #value.var = "rel"
  #offset = 1
  #pair.method="scheffe"
  
  data=as.data.frame(data)
  
  a <- aov(data[,value.var]~data[,flev], data=data)
  
  PHT=PostHocTest(a, method = pair.method, ordered = T)
  names(PHT) <- (flev)
  
  PHTb = PHT$Group
  Group2 = gsub( "-.*", "", rownames(PHTb))
  Group1 = gsub( ".*-", "", rownames(PHTb))
  GroupComp = paste0(Group1, "-", Group2)
  rownames(PHTb) = GroupComp
  PHTb[,'diff'] = -1 * PHTb[,'diff']
  PHTb[,'lwr.ci'] = -1 * PHTb[,'lwr.ci']
  PHTb[,'upr.ci'] = -1 * PHTb[,'upr.ci']
  ### Now we're happy
  PTb2 <- PHTb[,'pval']
  
  
  Tukey.labels <- multcompLetters(PTb2, reversed = T)['Letters']
  
  plot.labels <- names(Tukey.labels[['Letters']])
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  iqr.df <- ddply(data, flev, function (x) ((1.5*IQR(x[,value.var])+quantile(x[,value.var],0.75)))) %>% 
    `colnames<-`(c("Group","iqr"))
  data=data %>% inner_join(iqr.df, by="Group")
  d2=subset(data,data[,value.var]<data$iqr)
  
  boxplot.df <- ddply(d2, flev, function (x) (max(fivenum(x[,value.var])))+offset)
  #boxplot.df <- ddply(d2, flev, function (x) (2*max(fivenum(x[,value.var]))-min(fivenum(x[,value.var])))*offset)
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = F)
  
  return(labels.df)
}
