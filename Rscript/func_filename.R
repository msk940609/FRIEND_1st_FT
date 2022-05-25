filename <- function(filename){
  # check that ftmsObj is of the correct class #
  
  filename=paste0("Figure/",filename,"_",format(Sys.time(), "%Y_%m_%d"),".tiff")
  
  return(filename)
}



