
library('survival')
library('survminer')
setwd("D:\\分组_comb\\BRCA_female_white_grade_comb")
file = read.table("D:\\分组_comb\\BRCA_female_white_grade_comb\\BRCA_female_white_grade_comb(25).txt",sep = "\t",header = F)
cli = read.table("D:\\TCGA\\TCGA-BRCA mature\\BRCA_clinical.txt",sep = "\t",header = T )
row.names(cli)<-cli[,1]
cli<-cli[,-1]


ncol(file)
#> 生成有pairmatch分组临床样本
cli2 = matrix(nrow = (ncol(file)-1),ncol = 14,dimnames=list(c(),c("sampleid","classfication_of_tumor","tumor_stage"
                                                     ,"gender","year_to_birth","year_to_death"
                                                     ,"year_to_diagnosis","OS_time","age","deadORlive"
                                                     ,"race","alcohol","smoked","pairmatch")))
#> 生成空结果文件
result = matrix(nrow =14265810 ,ncol = 4,dimnames=list(c(),c("pairname","combi"
                                                                    ,"pval","pairmatch")))
count = 1
#> 进行一个1585090个pair的循环
print(system.time({
for (i in 2:50001) {
  #> 分别循环9种组合
  for (combi in c("low,low","low,mid","low,high",
                  "mid,low","mid,mid","mid,high",
                  "high,low","high,mid","high,high")) {
  #> 按列进行样本二分
    for (j in 2:ncol(file)) {
      if(file[i,j] == combi){
        cli2[(j-1),"classfication_of_tumor"] = cli[file[1,j],"classfication_of_tumor"]
        cli2[(j-1),"tumor_stage"] = cli[file[1,j],"tumor_stage"]
        cli2[(j-1),"gender"] = cli[file[1,j],"gender"]
        cli2[(j-1),"year_to_birth"] = cli[file[1,j],"year_to_birth"]
        cli2[(j-1),"year_to_death"] = cli[file[1,j],"year_to_death"]
        cli2[(j-1),"year_to_diagnosis"] = cli[file[1,j],"year_to_diagnosis"]
        cli2[(j-1),"OS_time"] = cli[file[1,j],"OS_time"]
        cli2[(j-1),"age"] = cli[file[1,j],"age"]
        cli2[(j-1),"deadORlive"] = cli[file[1,j],"deadORlive"]
        cli2[(j-1),"race"] = cli[file[1,j],"race"]
        cli2[(j-1),"alcohol"] = cli[file[1,j],"alcohol"]
        cli2[(j-1),"smoked"] = cli[file[1,j],"smoked"]
        cli2[(j-1),"pairmatch"] = "match"
        cli2[(j-1),"sampleid"] = file[1,j]
      }
      if(file[i,j] != combi){
        cli2[(j-1),"classfication_of_tumor"] = cli[file[1,j],"classfication_of_tumor"]
        cli2[(j-1),"tumor_stage"] = cli[file[1,j],"tumor_stage"]
        cli2[(j-1),"gender"] = cli[file[1,j],"gender"]
        cli2[(j-1),"year_to_birth"] = cli[file[1,j],"year_to_birth"]
        cli2[(j-1),"year_to_death"] = cli[file[1,j],"year_to_death"]
        cli2[(j-1),"year_to_diagnosis"] = cli[file[1,j],"year_to_diagnosis"]
        cli2[(j-1),"OS_time"] = cli[file[1,j],"OS_time"]
        cli2[(j-1),"age"] = cli[file[1,j],"age"]
        cli2[(j-1),"deadORlive"] = cli[file[1,j],"deadORlive"]
        cli2[(j-1),"race"] = cli[file[1,j],"race"]
        cli2[(j-1),"alcohol"] = cli[file[1,j],"alcohol"]
        cli2[(j-1),"smoked"] = cli[file[1,j],"smoked"]
        cli2[(j-1),"pairmatch"] = "mismatch"
        cli2[(j-1),"sampleid"] = file[1,j]
      }
    }
  #> 进行km生存分析 
    
    cli2 <- as.data.frame(cli2)
    cli2[,"OS_time"] = as.numeric(cli2[,"OS_time"])
    if(length(unique(cli2[,"pairmatch"])) != 1){
      attach(cli2)
      my_surv = Surv(OS_time,deadORlive == '0')
      fit = survfit(my_surv ~ pairmatch,data=cli2)
      b = surv_pvalue(fit, method = "Log-rank",data=cli2)
      p_val = format(b[1,"pval"],scientific=TRUE,digit=3)
      a = surv_summary(fit)
      surv1 = 0
      surv2 = 0
      count1 = 0
      count2 = 0
      for (k in 1:nrow(a)) {
        if(a[k,'strata'] == 'pairmatch=match'){
          surv1 = surv1 + a[k,'surv']
          count1 = count1 + 1
        }
        if(a[k,'strata'] == 'pairmatch=mismatch'){
          surv2 = surv2 + a[k,'surv']
          count2 = count2 + 1
        }
      }
      ave_surv1 = surv1/count1
      ave_surv2 = surv2/count2
    #> 写入结果文件 
      result[count,"pairname"] = file[i,1]
      result[count,"combi"] = combi
      result[count,"pval"] = p_val
    #> 记录哪一组生存率较高
      if(ave_surv1>ave_surv2){
        result[count,"pairmatch"] = "match"
      }else{
        result[count,"pairmatch"] = "mismatch"
      }
      count = count + 1
      survdiff(my_surv ~ pairmatch,data=cli2)
      detach(cli2)
    }
    
  }
    print(i*100/50001)
}
}))

write.table(result,file="BRCA_female_white_result(25).txt",sep="\t",quote=F,row.names=F)
