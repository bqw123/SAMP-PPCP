install.packages("tidyr")
library(tidyr)
setwd("D:\\分组_comb2\\TCGA-UCEC")
express = read.table("D:\\分组_comb2\\TCGA-UCEC\\female_white_grade.txt",header=T,row.names = NULL)
filew = "female_white_grade_comb.txt"

dim(express)

res3 = t(combn(express[,1],2))
nrow(res3)

exp = matrix(nrow = nrow(res3) ,ncol = ncol(express))
dim(exp)
colnames(exp) = colnames(express)

for (i in 1:ncol(express)) {
  res = t(combn(express[,i],2))
  res2 = paste(res[,1],res[,2],sep=",")
  res2 = matrix(res2,nrow(res),1)
  exp[,i] = res2[,1]
  print(i)
}

write.table(exp,file=filew,sep="\t",quote=F,row.names=F)

exp = c(1)
express1 = c(1)
express2 = c(1)


express1 = read.table("D:\\分组\\TCGA-PRAD\\male_white_grade_1_comb.txt",header=T)
express2 = read.table("D:\\分组\\TCGA-PRAD\\male_white_grade_2_comb.txt",header=T)

exp = cbind(express1,express2)
write.table(exp,file="male_white_grade_comb.txt",sep="\t",quote=F,row.names=F)
