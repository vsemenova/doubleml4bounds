rm(list=ls())
library(xtable)
my_path<-"/net/holyparkesec/data/tata/doubleml4bounds"


### print Appendix table

Table_oracle_std<-read.csv(paste0(my_path,"/Tables/Table_oracle_nonortho.csv"),row.names=1)
Table_oracle_prop<-read.csv(paste0(my_path,"/Tables/Table_oracle_ortho.csv"),row.names=1)

## finaltable
finaltable = matrix(0,12,8)
for (k in 1:3) {
  finaltable[4*(k-1)+1:4 ,1:4]<-as.matrix(Table_oracle_std[,4*(k-1)+1:4])
  finaltable[4*(k-1)+1:4,5:8]<-as.matrix(Table_oracle_prop[,4*(k-1)+1:4])
  
  
}

write.csv( finaltable,paste0(my_path,"/Tables/AppendixTable.csv"))
write.table(print(xtable(finaltable, type="latex",digits=2),include.rownames =FALSE ),paste0(my_path,"/Tables/AppendixTable.txt"),append=FALSE)


### print intermediate table


Table_lasso_nonortho<-read.csv(paste0(my_path,"/Tables/Table_lasso_nonortho.csv"),row.names=1)
Table_series_ortho<-read.csv(paste0(my_path,"/Tables/Table_series_ortho.csv"),row.names=1)

## finaltable
finaltable = matrix(0,12,8)
for (k in 1:3) {
  finaltable[4*(k-1)+1:4 ,1:4]<-as.matrix(Table_lasso_nonortho[,4*(k-1)+1:4])
  finaltable[4*(k-1)+1:4,5:8]<-as.matrix(Table_series_ortho[,4*(k-1)+1:4])
  
  
}

write.csv( finaltable,paste0(my_path,"/Tables/AppendixTable2.csv"))
write.table(print(xtable(finaltable, type="latex",digits=2),include.rownames =FALSE ),paste0(my_path,"/Tables/AppendixTable2.txt"),append=FALSE)
