library(data.table)
library(ggplot2)
library(dplyr)
library(writexl)
library(xlsx)
options(java.parameters = "- Xmx8192m")


####set wd bejing=====
flist=list.files(path = ".",pattern = ".xlsx")

dt=data.table()
for (i in 1:length(flist)) {
  temp=as.data.table(read.xlsx2(flist[i], sheetIndex = 1))
  temp$Sample=flist[i]
  
  dt=rbind(dt, temp,fill=TRUE)
  
}
dt
table(dt$Sample)

fwrite(dt,"FRIEND_beijing_merge.csv")




