rm(list=ls())
library(xtable)
my_path<-"/net/holyparkesec/data/tata/doubleml4bounds/empirical_application/"

### print main text table (Table 1) ####

Table_lasso<-read.csv(paste0(my_path,"/Tables/Lasso.csv"),row.names=1)
Table_forest<-read.csv(paste0(my_path,"/Tables/Forest.csv"),row.names=1)

## finaltable
finaltable = matrix(0,4,8)
finaltable[1:4 ,1:4]<-as.matrix(Table_lasso)
finaltable[1:4,5:8]<-as.matrix(Table_forest)

write.csv( finaltable,paste0(my_path,"/Tables/MainTable.csv"))
write.table(print(xtable(finaltable, type="latex",digits=3),include.rownames =FALSE ),paste0(my_path,"/Tables/MainTable.txt"),append=FALSE)
