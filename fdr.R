setwd("D:\\分组_comb2\\TCGA-UCEC\\female_white")
file = read.table("D:\\分组_comb2\\TCGA-UCEC\\female_white\\UCEC_female_white_result.txt",sep = "\t",header = F)

pvalue = file[,3]

pvalue1 = pvalue[-1]

fdr = p.adjust(pvalue1,method="fdr",n=length(pvalue1))

a = c("FDR",fdr)

file2 = cbind(file,a)


write.table(file2,file="UCEC_female_white_result_fdr.txt",sep="\t",quote=F,row.names=F,col.names = F)

colnames(file2)<-file2[1,]
file2<-file2[-1,]

file2[,5] = as.numeric(file2[,5])

file2 = file2[order(file2[,5]),]


for (i in 1:nrow(file2)) {
  if(file2[i + 1,"FDR"] >= 0.05){
    print(i)
  break;
  }
}

file3 = file2[1:114898,]

write.table(file3,file="UCEC_female_white_result_fdr_0.05.txt",sep="\t",quote=F,row.names=F,col.names = T)
