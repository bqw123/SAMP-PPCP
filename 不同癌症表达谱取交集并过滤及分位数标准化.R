setwd("D:\\TCGA express")

#> 获得存储路径下所有的csv格式文件的文件名

files<-dir(path = "D:\\TCGA express",
           full.names = T,
           pattern = ".csv")

files

#> 批量将数据读入

library(tidyverse)
df<-map(files,read.csv)
class(df) 

#> 数据取交集

df1<-reduce(df,inner_join)

rownames(df1) = df1[,1]

df1 = df1[,-1]

#> 表达矩阵的过滤 过滤掉在所有样本表达量为0的mirna


k = apply(df1, 1, function(x) sum(x > 0) > 0);table(k)

res = df1[k,]

#> 把文件写入为txt

write.table(res,file="expressfilter.txt",sep="\t",quote=F,row.names=T)

#> 进行分位数标准化

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("preprocessCore")
library("preprocessCore")

expressfilterquanor = normalize.quantiles(x=as.matrix(res))

rownames(expressfilterquanor) = rownames(res)

colnames(expressfilterquanor) = colnames(res)

#> 把文件写入为txt

write.table(expressfilterquanor,file="expressfilterquanor.txt",sep="\t",quote=F,row.names=T)




