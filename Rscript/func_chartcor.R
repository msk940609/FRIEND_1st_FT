chart.Correlationx <-
  function (R, histogram = TRUE, method=c("pearson", "kendall", "spearman"), pch=1,blk=0.5, ...)
  { # @author R Development Core Team
    # @author modified by Peter Carl
    # Visualization of a Correlation Matrix. On top the (absolute) value of the
    # correlation plus the result of the cor.test as stars. On botttom, the
    # bivariate scatterplots, with a fitted line
    
    x = checkData(R, method="matrix")
    
    if(missing(method)) method=method[1] #only use one
    cormeth <- method
    
    # Published at http://addictedtor.free.fr/graphiques/sources/source_137.R
    panel.cor <- function(x, y, digits=2, prefix="", use="pairwise.complete.obs", method=cormeth, cex.cor, ...)
    {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- cor(x, y, use=use, method=method) # MG: remove abs here
      txt <- format(c(r, 0.123456789), digits=digits)[1]
      txt <- paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
      
      test <- cor.test(as.numeric(x),as.numeric(y), method=method)
      # borrowed from printCoefmat
      Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("***", "**", "*", ".", " "))
      # MG: add abs here and also include a 30% buffer for small numbers
      text(0.5, 0.5, txt, cex = cex * (abs(r) + .3) / 1.3)
      text(.8, .8, Signif, cex=cex, col=2)
    }
    f <- function(t) {
      dnorm(t, mean=mean(x), sd=sd.xts(x) )
    }
    
    #remove method from dotargs
    dotargs <- list(...)
    dotargs$method <- NULL
    rm(method)
    
    hist.panel = function (x, bf=blk,...=NULL ) {
      par(new = TRUE)
      hist(x,
           col = "light gray",
           probability = TRUE,
           axes = FALSE,
           main = "",
           breaks = "FD",
           ylim = c(0,max(density(x,na.rm=TRUE)$y)*bf)
           )
      lines(density(x, na.rm=TRUE),
            col = "red",
            lwd = 1)
      #lines(f, col="blue", lwd=1, lty=1) how to add gaussian normal overlay?
      rug(x)
    }
    
    # Draw the chart
    if(histogram)
      pairs(x, gap=0, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=hist.panel ,pch=pch, ...)
    else
      pairs(x, gap=0, lower.panel=panel.smooth, upper.panel=panel.cor, pch=pch, ...)
  }



