setwd("D:\\分组\\TCGA-UCEC")

express = read.table("D:\\分组\\TCGA-UCEC\\female_white.txt",header=T)

filew = "female_white_grade.txt"

rownames(express) = express[,1]
express = express[,-1]

dim(express)

nrow(express)
ncol(express)

for (i in 1:nrow(express)) {
  for (j in 1:ncol(express)) {
    if(as.numeric(express[i,j])>=0 && as.numeric(express[i,j])<3.254326){
      express[i,j]="low"
    }else if(as.numeric(express[i,j])>=3.254326 && as.numeric(express[i,j])<41.70731 ){
      express[i,j]="mid"
    }else{
      express[i,j]="high"
    }
  }
}

write.table(express,file=filew,sep="\t",quote=F,row.names=T)

