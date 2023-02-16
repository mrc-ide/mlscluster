#rm(list=ls())
library(ape)
library(stringr)
library(tidyverse)
library(data.table)
#install.packages("splitstackshape")
library(splitstackshape)
library(dplyr)
library(lubridate)

options(scipen=999)

sc2_md_curated <- readRDS("rds/sc2_md_curated.rds")
sc2_tre_curated <- readRDS("rds/sc2_tre_curated.rds")
# cog_md has 1,729,244, sc2_md_curated has 1,275,670

cog_md <- fread("data/cog_metadata_pillar2_england_excl1stWave_09May.csv.gz", sep=",", header=F, select=c(1,5,6,7,14), nThread=6, tmpdir="/home/vinibfranc/")
colnames(cog_md) <- c("sequence_name","sample_date","epi_week","lineage","mutations")
cog_md <- include_major_lineage_column(cog_md) #1,727,840
cog_md <- cog_md %>% relocate(major_lineage, .after = lineage)
#cog_md = cog_md %>% add_row(sequence_name = "Wuhan/WH04/2020",  sample_date=as.Date("2019-12-30"), epi_week=1, lineage="B", mutations="synSNP:C8782T|orf8:L84S")
cog_md$sample_time <- decimal_date(cog_md$sample_date)
cog_md <- cog_md %>% relocate(sample_time, .after = sample_date)
#cog_md <- include_major_lineage_column(cog_md)

cog_md_muts <- cSplit(cog_md, 'mutations', '|', 'long')
cog_md_muts$mutations <- toupper(cog_md_muts$mutations)
cog_md_muts$type <- sub("\\:.*", "", cog_md_muts$mutations)
cog_md_muts$mut <- sub('.*:', "", cog_md_muts$mutations)
colnames(cog_md_muts) <- c("sequence_name","sample_date","sample_time","epi_week","major_lineage","lineage","mut_and_type","type","mut") #"major_lineage"

cog_md_muts$type  <- factor(cog_md_muts$type, levels = c("SYNSNP","ORF1AB","S","ORF3A","E","M","ORF6","ORF7A","ORF8","N","ORF10"))

# >50k unique mutations
seq_muts_df_counts <- cog_md_muts %>% group_by(mut_and_type)  %>% summarise(n = n()) %>% arrange(desc(n)) %>% mutate(percent = 100*(n / nrow(cog_md)))
seq_muts_df_counts$protein <- sub("\\:.*", "", seq_muts_df_counts$mut_and_type)

prot_lengths <- read.csv("config/aa_length_prots.tsv", sep="\t", header=TRUE)
prot_lengths$protein <- toupper(prot_lengths$protein)

seq_muts_df_counts <- seq_muts_df_counts %>% inner_join(prot_lengths, by="protein") # TODO keep SYNSNPs? I don't think it is needed here because all processing of positions and region attribution is done on `join_mut_sites_homoplasies.R`
seq_muts_df_counts$syn_non_syn <- ifelse(seq_muts_df_counts$protein == "SYNSNP", "syn", "non-syn")

saveRDS(seq_muts_df_counts, "rds/polymorphic_site_counts.rds")
