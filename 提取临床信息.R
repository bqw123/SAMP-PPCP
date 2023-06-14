
install.packages("rjson")
library("rjson")

setwd("D:\\TCGA\\TCGA-BRCA mature")

#> 读入clinical.json文件

clinical_traits <- fromJSON(file = "clinical.cart.2021-09-07.json")

#> 计算文件长度n

n = length(clinical_traits)

#> 初始化变量

id = classfication_of_tumor = c(rep(0, n))
tumor_stage = gender = c(rep(0, n))
year_to_birth = year_to_death =  c(rep(0, n))
year_to_diagnosis = days_to_death = c(rep(0, n))
age = deadORlive = race = alcohol = smoked = c(rep(0, n))
ajcc_pathologic_t = ajcc_pathologic_n = ajcc_pathologic_m = c(rep(0, n))
ajcc_pathologic_stage = c(rep(0, n))
tumor_grade = c(rep(0, n))
#> 利用一个for循环由json文件中提取信息

for (i in 1:n) {
  id[i] = clinical_traits[[i]]$diagnoses[[1]]$submitter_id
  classfication_of_tumor[i]=clinical_traits[[i]]$diagnoses[[1]]$classification_of_tumor
  tumor_stage[i] = clinical_traits[[i]]$diagnoses[[1]]$tumor_stage
  gender[i] = clinical_traits[[i]]$demographic$gender
  year_to_birth[i] = ifelse(
    is.null(clinical_traits[[i]]$demographic$year_of_birth),
    "notReport",
    clinical_traits[[i]]$demographic$year_of_birth
  )
  year_to_death[i] = ifelse(
    is.null(clinical_traits[[i]]$demographic$year_of_death),
    "notReport",
    clinical_traits[[i]]$demographic$year_of_death
  )
  year_to_diagnosis[i] = ifelse(
    is.null(clinical_traits[[i]]$diagnoses[[1]]$year_of_diagnosis),
    "notReport",
    clinical_traits[[i]]$diagnoses[[1]]$year_of_diagnosis
  )
  days_to_death[i] = ifelse(
    is.null(clinical_traits[[i]]$demographic$days_to_death),
    ifelse(
      is.null(clinical_traits[[i]]$diagnoses[[1]]$days_to_last_follow_up),
      "notReport",               
      clinical_traits[[i]]$diagnoses[[1]]$days_to_last_follow_up),               
    clinical_traits[[i]]$demographic$days_to_death
  )
  age[i] = ifelse(
    is.null(clinical_traits[[i]]$demographic$age_at_index),
    "notReport",
    clinical_traits[[i]]$demographic$age_at_index
  )
  tumor_grade = ifelse(
    is.null(clinical_traits[[i]]$diagnoses[[1]]$tumor_grade),
    "notReport",
    clinical_traits[[i]]$diagnoses[[1]]$tumor_grade
  )
  ajcc_pathologic_t[i] = ifelse(
    is.null(clinical_traits[[i]]$diagnoses[[1]]$ajcc_pathologic_t),
    "notReport",
    clinical_traits[[i]]$diagnoses[[1]]$ajcc_pathologic_t
  )
  ajcc_pathologic_n[i] = ifelse(
    is.null(clinical_traits[[i]]$diagnoses[[1]]$ajcc_pathologic_n),
    "notReport",
    clinical_traits[[i]]$diagnoses[[1]]$ajcc_pathologic_n
  )
  ajcc_pathologic_m[i] = ifelse(
    is.null(clinical_traits[[i]]$diagnoses[[1]]$ajcc_pathologic_m),
    "notReport",
    clinical_traits[[i]]$diagnoses[[1]]$ajcc_pathologic_m
  )
  ajcc_pathologic_stage[i] = ifelse(
    is.null(clinical_traits[[i]]$diagnoses[[1]]$ajcc_pathologic_stage),
    "notReport",
    clinical_traits[[i]]$diagnoses[[1]]$ajcc_pathologic_stage
  )
  deadORlive[i] = ifelse(
    is.null(clinical_traits[[i]]$demographic$vital_status),
    "notReport",
    clinical_traits[[i]]$demographic$vital_status
  )
  race[i] = ifelse(
    is.null(clinical_traits[[i]]$demographic$race),
    "notReport",
    clinical_traits[[i]]$demographic$race
  )
  alcohol[i] = ifelse(
    is.null(clinical_traits[[i]]$exposures[[1]]$alcohol_history),
    "notReprot",
    clinical_traits[[i]]$exposures[[1]]$alcohol_history
  )
  smoked[i] = ifelse(
    is.null(clinical_traits[[i]]$exposures[[1]]$years_smoked),
    "notReport",
    clinical_traits[[i]]$exposures[[1]]$years_smoked
  )
}

#> 将提取的信息做成一个dataFrame
gastric_clinic <- data.frame(
  id,
  classfication_of_tumor,
  tumor_stage,
  gender,
  year_to_birth,
  year_to_death,
  year_to_diagnosis,
  days_to_death,
  tumor_grade,
  ajcc_pathologic_stage,
  ajcc_pathologic_t,
  ajcc_pathologic_n,
  ajcc_pathologic_m,
  age,
  deadORlive,
  race,
  alcohol,
  smoked
)

write.table(gastric_clinic,file="BRCA_clinical4.txt",sep="\t",quote=F,row.names=F)
