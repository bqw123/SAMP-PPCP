library("survival")
library("survminer")

setwd("/home/shijiangcheng/wubaoqin/cox")
result = read.table("/home/shijiangcheng/wubaoqin/cox/BRCA_result_fdr_0.05_samn_filter.txt",sep = "\t",header=T)
comb = read.table("/home/shijiangcheng/wubaoqin/cox/BRCA_grade_comb.txt",sep = "\t",header=F)
cli = read.table("/home/shijiangcheng/wubaoqin/cox/BRCA_clinical2.txt",sep = "\t",header = T )
colnames(comb)<-comb[1,]
row.names(comb)<-comb[,1]

for (i in 1:nrow(result)) {
  
  for (j in 1:nrow(cli)) {
    if(comb[result[i,"pairname"],cli[j,"id"]] == result[i,"combi"]){
      cli[j,"belong"] = 1
    }else{
      cli[j,"belong"] = 0
    }
    print(i)
    print(j)
  }
  with(cli,{
  my_surv <<- Surv(OS_time,deadORlive == '0')
  cox <<- coxph(my_surv ~ belong+age+gender+race,data=cli)
})
  a = summary(cox)
  coxpvalue = format(a$coefficients[1,5],scientific=TRUE,digit=3)
  result[i,"coxpvalue"] = coxpvalue
  
}

write.table(result,file= "BRCA_result_fdr_0.05_samn_filter_cox.txt",sep="\t",quote=F,row.names=F,col.names = T)





