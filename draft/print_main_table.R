## read in tables number 3 and 4.

rm(list=ls())
library(xtable)
my_path<-"/net/holyparkesec/data/tata/doubleml4bounds"

### print main text table (Table 1) ####

Table_lasso<-read.csv(paste0(my_path,"/Tables/Table_lasso1.csv"),row.names=1)
Table_series<-read.csv(paste0(my_path,"/Tables/Table_series.csv"),row.names=1)

## finaltable
finaltable = matrix(0,12,8)
for (k in 1:3) {
  finaltable[4*(k-1)+1:4 ,1:4]<-as.matrix(Table_series[,4*(k-1)+1:4])
  finaltable[4*(k-1)+1:4,5:8]<-as.matrix(Table_lasso[,4*(k-1)+1:4])
  
  
}

write.csv( finaltable,paste0(my_path,"/Tables/MainTable.csv"))
write.table(print(xtable(finaltable, type="latex",digits=2),include.rownames =FALSE ),paste0(my_path,"/Tables/MainTable.txt"),append=FALSE)

