######### Venn diagram and contingency table: compare rejections aLKF, global KF and Fixed-LKF #############

# UKB GWAS
#
# Class: R script
#
# This script compares the rejections between the aLKF, global KF and Fixed-LKF

rm(list = ls())
library(tidyverse)

# set parameters
my_population = "british"
resolution = "4"
pheno.name = "whr"
covar_identifiyer = "sex"
suffix = ""
num.int = 1

stats.path = paste0("results_path/", my_population, 
                    "_res", resolution, "_", pheno.name, "_cov_", covar_identifiyer)

# rejections fixed-LKF
selections_jointly = read_delim(sprintf("%s_lasso_discoveries_male_female_joint_real.txt", stats.path), delim=" ")

# rejections aLKF
selections_sskf = read_csv(paste0("results_path/lasso_selections_interactions_weight0.25_" ,my_population, "_res",
                                  resolution, "_", pheno.name, "_cov_", covar_identifiyer, "_", suffix, "_", num.int ,".csv"))

####### venn diagram ###########

unique_chr_group_rejected_sskf = unique(selections_sskf$CHR_Group)
unique_chr_group_rejected_joint = unique(selections_jointly$CHR_Group)
unique_chr_group_rejected_vanilla = unique(selections_vanilla %>% mutate(CHR_Group = paste0(CHR, "_", Group)) %>% pull(CHR_Group))

# common elements 
common = Reduce(intersect, list(unique_chr_group_rejected_sskf, unique_chr_group_rejected_joint, unique_chr_group_rejected_vanilla))
length(common)

# unique 
unique_to_vanilla <- setdiff(unique_chr_group_rejected_vanilla, union(unique_chr_group_rejected_sskf, unique_chr_group_rejected_joint))
length(unique_to_vanilla)

unique_to_joint <- setdiff(unique_chr_group_rejected_sskf, union(unique_chr_group_rejected_vanilla, unique_chr_group_rejected_joint))
length(unique_to_joint)

unique_to_sskf <- setdiff(unique_chr_group_rejected_joint, union(unique_chr_group_rejected_sskf, unique_chr_group_rejected_vanilla))
length(unique_to_sskf)

# overlapping
length(intersect(unique_chr_group_rejected_sskf, unique_chr_group_rejected_joint ))
length(intersect(unique_chr_group_rejected_vanilla, unique_chr_group_rejected_joint ))
length(intersect(unique_chr_group_rejected_sskf, unique_chr_group_rejected_vanilla ))

sum(!(unique_chr_group_rejected_joint %in% unique_chr_group_rejected_sskf))
sum(!(unique_chr_group_rejected_sskf %in% unique_chr_group_rejected_joint))


############ CONTINGENCY TABLE ############### 

######## __aLKF ALL ###########

chr_group_selections_aLKF_all = selections_sskf %>% 
  filter(Interaction == "none") %>% 
  pull(CHR_Group)

# not detected 
sum(!(chr_group_selections_aLKF_all %in% selections_jointly$CHR_Group))

# detected
chr_group_sskf_all_in_naLKF = selections_jointly %>%
  filter(CHR_Group %in% chr_group_selections_aLKF_all) %>% 
  group_by(CHR_Group) %>% 
  mutate(rejected_by_group = n())

# rejected once
chr_group_sskf_all_in_naLKF_rej_once = chr_group_sskf_all_in_naLKF %>% filter(rejected_by_group == 1)
table(chr_group_sskf_all_in_naLKF_rej_once$sample)

# rejected twice
chr_group_sskf_all_in_naLKF_rej_twice = chr_group_sskf_all_in_naLKF %>%
  filter(rejected_by_group == 2)
dim(chr_group_sskf_all_in_naLKF_rej_twice)
length(unique(chr_group_sskf_all_in_naLKF_rej_twice$CHR_Group))

###### _____aLKF female #########

chr_group_selections_aLKF_female = selections_sskf %>% 
  group_by(CHR_Group) %>% 
  mutate(rejected_by_group = n()) %>% 
  filter(rejected_by_group == 1) %>%
  filter(Interaction == "sex" & environment == 1) %>% 
  pull(CHR_Group)

# not detected 
sum(!(chr_group_selections_aLKF_female %in% selections_jointly$CHR_Group))

# detected
chr_group_sskf_female_in_naLKF = selections_jointly %>%
  filter(CHR_Group %in% chr_group_selections_aLKF_female) %>% 
  group_by(CHR_Group) %>% 
  mutate(rejected_by_group = n())

# rejected once
chr_group_sskf_female_in_naLKF_rej_once = chr_group_sskf_female_in_naLKF %>% filter(rejected_by_group == 1)
table(chr_group_sskf_female_in_naLKF_rej_once$sample)

# rejected twice
chr_group_sskf_female_in_naLKF_rej_twice = chr_group_sskf_female_in_naLKF %>%
  filter(rejected_by_group == 2)
dim(chr_group_sskf_female_in_naLKF_rej_twice)
length(unique(chr_group_sskf_female_in_naLKF_rej_twice$CHR_Group))

###### _____aLKF male #########


chr_group_selections_aLKF_male = selections_sskf %>% 
  group_by(CHR_Group) %>% 
  mutate(rejected_by_group = n()) %>% 
  filter(rejected_by_group == 1) %>%
  filter(Interaction == "sex" & environment == 2) %>% 
  pull(CHR_Group)

# not detected 
sum(!(chr_group_selections_aLKF_male %in% selections_jointly$CHR_Group))

# detected
chr_group_sskf_male_in_naLKF = selections_jointly %>%
  filter(CHR_Group %in% chr_group_selections_aLKF_male) %>% 
  group_by(CHR_Group) %>% 
  mutate(rejected_by_group = n())

# rejected once
chr_group_sskf_male_in_naLKF_rej_once = chr_group_sskf_male_in_naLKF %>% filter(rejected_by_group == 1)
table(chr_group_sskf_male_in_naLKF_rej_once$sample)

# rejected twice
chr_group_sskf_male_in_naLKF_rej_twice = chr_group_sskf_male_in_naLKF %>%
  filter(rejected_by_group == 2)
dim(chr_group_sskf_male_in_naLKF_rej_twice)
length(unique(chr_group_sskf_male_in_naLKF_rej_twice$CHR_Group))

###### _____aLKF female and male #########


chr_group_selections_aLKF_fem_male = selections_sskf %>% 
  group_by(CHR_Group) %>% 
  mutate(rejected_by_group = n()) %>% 
  filter(rejected_by_group == 2) %>%
  pull(CHR_Group) %>% 
  unique()

# not detected 
sum(!(chr_group_selections_aLKF_fem_male %in% selections_jointly$CHR_Group))

# detected
chr_group_sskf_fem_male_in_naLKF = selections_jointly %>%
  filter(CHR_Group %in% chr_group_selections_aLKF_fem_male) %>% 
  group_by(CHR_Group) %>% 
  mutate(rejected_by_group = n())

# rejected once
chr_group_sskf_fem_male_in_naLKF_rej_once = chr_group_sskf_fem_male_in_naLKF %>% filter(rejected_by_group == 1)
table(chr_group_sskf_fem_male_in_naLKF_rej_once$sample)

# rejected twice
chr_group_sskf_fem_male_in_naLKF_rej_twice = chr_group_sskf_fem_male_in_naLKF %>%
  filter(rejected_by_group == 2)
dim(chr_group_sskf_fem_male_in_naLKF_rej_twice)
length(unique(chr_group_sskf_fem_male_in_naLKF_rej_twice$CHR_Group))


############# naLKF ######### 

selections_joint_once = selections_jointly %>% 
  group_by(CHR_Group) %>% 
  mutate(rej_by_group = n()) %>% 
  filter(rej_by_group == 1)

table(selections_joint_once$sample)

selections_joint_twice = selections_jointly %>% 
  group_by(CHR_Group) %>% 
  mutate(rej_by_group = n()) %>% 
  filter(rej_by_group == 2) 

nrow(selections_joint_twice)/2

######## __naLKF FEMALE ###########


chr_group_selections_naLKF_female = selections_joint_once %>% 
  filter(sample == "female") %>% 
  pull(CHR_Group)

# not detected 
sum(!(chr_group_selections_naLKF_female %in% selections_sskf$CHR_Group))

# detected
chr_group_aLKF_in_naLKF = selections_sskf %>%
  filter(CHR_Group %in% chr_group_selections_naLKF_female) %>% 
  group_by(CHR_Group) %>% 
  mutate(rejected_by_group = n()) %>% 
  mutate(sample = ifelse(Interaction == "none", "all", 
                         ifelse(environment == 1 & Interaction == "sex", "female", "male")))

# rejected once
chr_group_aLKF_in_naLKF_once = chr_group_aLKF_in_naLKF %>% 
  filter(rejected_by_group == 1) 
table(chr_group_aLKF_in_naLKF_once$sample)

# rejected twice
chr_group_aLKF_in_naLKF_twice = chr_group_aLKF_in_naLKF %>%
  filter(rejected_by_group == 2)
dim(chr_group_aLKF_in_naLKF_twice)
length(unique(chr_group_aLKF_in_naLKF_twice$CHR_Group))


######## __naLKF MALE ###########

chr_group_selections_naLKF_male = selections_joint_once %>% 
  filter(sample == "male") %>% 
  pull(CHR_Group)

# not detected 
sum(!(chr_group_selections_naLKF_male %in% selections_sskf$CHR_Group))

# detected
chr_group_aLKF_in_naLKF = selections_sskf %>%
  filter(CHR_Group %in% chr_group_selections_naLKF_male) %>% 
  group_by(CHR_Group) %>% 
  mutate(rejected_by_group = n()) %>% 
  mutate(sample = ifelse(Interaction == "none", "all", 
                         ifelse(environment == 1 & Interaction == "sex", "female", "male")))

# rejected once
chr_group_aLKF_in_naLKF_once = chr_group_aLKF_in_naLKF %>% 
  filter(rejected_by_group == 1) 
table(chr_group_aLKF_in_naLKF_once$sample)

# rejected twice
chr_group_aLKF_in_naLKF_twice = chr_group_aLKF_in_naLKF %>%
  filter(rejected_by_group == 2)
dim(chr_group_aLKF_in_naLKF_twice)
length(unique(chr_group_aLKF_in_naLKF_twice$CHR_Group))



######## __naLKF MALE + FEMALE ###########

chr_group_selections_naLKF_male_female = selections_joint_twice %>% 
  pull(CHR_Group) %>% 
  unique()

# not detected 
sum(!(chr_group_selections_naLKF_male_female %in% selections_sskf$CHR_Group))

# detected
chr_group_aLKF_in_naLKF = selections_sskf %>%
  filter(CHR_Group %in% chr_group_selections_naLKF_male_female) %>% 
  group_by(CHR_Group) %>% 
  mutate(rejected_by_group = n()) %>% 
  mutate(sample = ifelse(Interaction == "none", "all", 
                         ifelse(environment == 1 & Interaction == "sex", "female", "male")))

# rejected once
chr_group_aLKF_in_naLKF_once = chr_group_aLKF_in_naLKF %>% 
  filter(rejected_by_group == 1) 
table(chr_group_aLKF_in_naLKF_once$sample)

# rejected twice
chr_group_aLKF_in_naLKF_twice = chr_group_aLKF_in_naLKF %>%
  filter(rejected_by_group == 2)
dim(chr_group_aLKF_in_naLKF_twice)
length(unique(chr_group_aLKF_in_naLKF_twice$CHR_Group))

