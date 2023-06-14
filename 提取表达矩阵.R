add = "D:\\TCGA\\TCGA-UCEC mature"

#> 整理表达矩阵
#> 先任选两个counts文件读取，并观察geneid的顺序是否一致。

add2 = "D:\\TCGA\\TCGA-UCEC mature\\"
fileadd = "D:\\TCGA\\TCGA-UCEC mature\\fe9231a6-59bc-4ce4-addb-301b1456ff41\\16a35678-23e1-44f6-b58b-f63f8c1fd3ff.mirbase21.isoforms.quantification.txt"
metadatename = "metadata.cart.2021-09-07.json"
resultname = "UCEC express2.txt"
setwd(add)
library(tidyverse)
x = read.table(fileadd,sep = "\t",header = T) 
head(x)

#> 接下来需要：
#> 只要成熟体id和count两列
#> 按照成熟体分组求和，得出每个成熟体的count之和

x2 = x %>%
  dplyr::select(c(6,4)) %>% 
  group_by(miRNA_region) %>% 
  summarise(read_count = sum(reads_per_million_miRNA_mapped))

count_files = dir(add2,pattern = "*isoforms.quantification.txt$",recursive = T)

exp = list()
for(i in 1:length(count_files)){
  exp[[i]] <- read.table(paste0(add2,count_files[[i]]),sep = "\t",header = T) %>%
    dplyr::select(c(6,4)) %>% 
    group_by(miRNA_region) %>% 
    summarise(read_count = sum(reads_per_million_miRNA_mapped)) 
}
sapply(exp, nrow)
m = Reduce(function(x, y) merge(x, y, by= 'miRNA_region',all = T), exp)
m[is.na(m)]=0
exp <- column_to_rownames(m,var = "miRNA_region")
exp = exp[-((nrow(exp)-2):nrow(exp)),]
dim(exp)
#> [1] 1753   45
exp[1:4,1:4]

#> 行名转换
#> GDC数据库使用的mirBasev21版本的id，我们转换时需要使用相同的版本,使用miRBaseVersions.db包

table(str_detect(rownames(exp),"mature,"))
#> 
#> TRUE 
#> 1753
rownames(exp) = str_remove(rownames(exp),"mature,")
rowna = rownames(exp)
library(miRBaseVersions.db)
mh <- select(miRBaseVersions.db,
             keys = rownames(exp),
             keytype = "MIMAT",
             columns = c("ACCESSION","NAME","VERSION"))
head(mh)
mh = mh[mh$VERSION=="21",]

mh = mh[match(rownames(exp),mh$ACCESSION),]

identical(rownames(exp),mh$ACCESSION)

table(!duplicated(mh$NAME))

rownames(exp) = mh$NAME

exp[1:4,1:4]

#> 列名转换

x = read.table(fileadd,sep = "\t",header = T) %>%
  dplyr::select(c(6,1)) %>% 
  distinct()
x$miRNA_region = str_remove(x$miRNA_region,"mature,")
x = merge(x,mh,by.x = "miRNA_region",by.y = "ACCESSION")  
head(x)

meta <- jsonlite::fromJSON(metadatename)
ID = sapply(meta$associated_entities,
            function(x){x$entity_submitter_id})
file2id = data.frame(file_name = meta$file_name,
                     ID = ID)
count_files2 = stringr::str_split(count_files,"/",simplify = T)[,2]
file2id = file2id[match(count_files2,file2id$file_name),]
colnames(exp) = file2id$ID
exp[1:4,1:4]

write.table(exp,file=resultname,sep="\t",quote=F,row.names=T)


