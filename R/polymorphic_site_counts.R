#rm(list=ls())
library(ape)
library(stringr)
library(tidyverse)
library(data.table)
#install.packages("splitstackshape")
library(splitstackshape)
library(dplyr)
library(lubridate)
library(readr)

#source("R/utils.R")

options(scipen=999)

sc2_md_curated <- readRDS("rds/sc2_md_curated.rds")
sc2_tre_curated <- readRDS("rds/sc2_tre_curated.rds")
# cog_md has 1,729,244, sc2_md_curated has 1,275,670

#cog_md <- fread("data/cog_metadata_pillar2_england_excl1stWave_09May.csv.gz", sep=",", header=F, select=c(1,5,6,7,14), nThread=6, tmpdir="/home/vinibfranc/")
#colnames(cog_md) <- c("sequence_name","sample_date","epi_week","lineage","mutations")
#cog_md <- cog_md %>% select(sequence_name,sample_date,epi_week,lineage,mutations)
#cog_md <- include_major_lineage_column(cog_md) #1,727,840
#cog_md <- cog_md %>% relocate(major_lineage, .after = lineage)
#cog_md = cog_md %>% add_row(sequence_name = "Wuhan/WH04/2020",  sample_date=as.Date("2019-12-30"), epi_week=1, lineage="B", mutations="synSNP:C8782T|orf8:L84S")
#cog_md$sample_time <- decimal_date(cog_md$sample_date)
#cog_md <- cog_md %>% relocate(sample_time, .after = sample_date)

cog_md_muts <- cSplit(sc2_md_curated, 'mutations', '|', 'long') #cog_md
#cog_md_muts$mutations <- toupper(cog_md_muts$mutations)
cog_md_muts$type <- sub("\\:.*", "", cog_md_muts$mutations)
cog_md_muts$mut <- sub('.*:', "", cog_md_muts$mutations)

#colnames(cog_md_muts) <- c("sequence_name","sample_date","sample_time","epi_week","major_lineage","lineage","mut_and_type","type","mut") #"major_lineage"

cog_md_muts$type  <- factor(cog_md_muts$type, levels = c("SYNSNP","ORF1AB","S","ORF3A","E","M","ORF6","ORF7A","ORF8","N","ORF10"))

# Start here to avoid computing everything again
cog_md_muts <- readRDS("rds/cog_md_muts_FDR.rds")

# >50k unique mutations
seq_muts_df_counts <- cog_md_muts[cog_md_muts$sample_date <= as.Date("2021-11-15"),] #"2020-12-31"
seq_muts_df_counts <- seq_muts_df_counts %>% group_by(mutations)  %>% summarise(n = n()) %>% arrange(desc(n)) %>% mutate(percent = 100*(n / nrow(seq_muts_df_counts)))
# seq_muts_df_counts$protein <- sub("\\:.*", "", seq_muts_df_counts$mut_and_type)
# 
# prot_lengths <- read.csv("config/aa_length_prots.tsv", sep="\t", header=TRUE)
# prot_lengths$protein <- toupper(prot_lengths$protein)
# 
# seq_muts_df_counts <- seq_muts_df_counts %>% inner_join(prot_lengths, by="protein") # TODO keep SYNSNPs? I don't think it is needed here because all processing of positions and region attribution is done on `join_mut_sites_homoplasies.R`
# seq_muts_df_counts$syn_non_syn <- ifelse(seq_muts_df_counts$protein == "SYNSNP", "syn", "non-syn")

#saveRDS(seq_muts_df_counts, "rds/polymorphic_site_counts.rds")

### SEQUENCING ARTIFACT INSPECTION

cog_md_inspection <- cog_md[(cog_md$sample_time > 2021.871 & cog_md$sample_time < 2022.326),] #& cog_md$major_lineage=="Omicron_BA.1.*"
cog_md_inspection <- cog_md_inspection[cog_md_inspection$sequence_name %in% sc2_tre_curated$tip.label, ] #875k to 594k

rm(cog_md, sc2_md_curated, sc2_tre_curated)
gc()

# # counts of potential artifacts without other variables
# artifacts <- cog_md_muts[cog_md_muts$mut_and_type=="S:K417N" | cog_md_muts$mut_and_type=="S:N440K" | cog_md_muts$mut_and_type=="S:G446S",]
# # TRIED X and "-"
# artifacts_t <- plyr::count(artifacts, "mut_and_type")

# counts of artifacts during Omicron period
artifacts2 <- cog_md_muts[(cog_md_muts$mut_and_type=="S:K417N" | cog_md_muts$mut_and_type=="S:N440K" | cog_md_muts$mut_and_type=="S:G446S") & (cog_md_muts$sample_time > 2021.871 & cog_md_muts$sample_time < 2022.326),] #Omicron period
#artifacts2_adj <- merge(artifacts2, cog_md_inspection, by="sequence_name")
artifacts2_417 <- cog_md_muts[(cog_md_muts$mut_and_type=="S:K417N") & (cog_md_muts$sample_time > 2021.871 & cog_md_muts$sample_time < 2022.326),] #Omicron period
artifacts2_417_t <- merge(artifacts2_417, cog_md_inspection, by="sequence_name")
artifacts2_417_a <- artifacts2_417[!artifacts2_417$sequence_name %in% artifacts2_417_t$sequence_name,]

artifacts2_440 <- cog_md_muts[(cog_md_muts$mut_and_type=="S:N440K") & (cog_md_muts$sample_time > 2021.871 & cog_md_muts$sample_time < 2022.326),] #Omicron period
artifacts2_440_t <- merge(artifacts2_440, cog_md_inspection, by="sequence_name")
artifacts2_446 <- cog_md_muts[(cog_md_muts$mut_and_type=="S:G446S") & (cog_md_muts$sample_time > 2021.871 & cog_md_muts$sample_time < 2022.326),] #Omicron period
artifacts2_446_t <- merge(artifacts2_446, cog_md_inspection, by="sequence_name")
#artifacts_t2 <- plyr::count(artifacts2, "mut_and_type")

# # counts of artifacts before Omicron period
# artifacts3 <- cog_md_muts[(cog_md_muts$mut_and_type=="S:K417N" | cog_md_muts$mut_and_type=="S:N440K" | cog_md_muts$mut_and_type=="S:G446S") & !(cog_md_muts$sample_time > 2021.871 & cog_md_muts$sample_time < 2022.326),] #Omicron period
# artifacts_t3 <- plyr::count(artifacts3, "mut_and_type")

### CHANGED PATHS ADDING new_try_updated_md below
write.csv(cog_md_inspection$sequence_name, file="data/omicron_artifact_inspection.txt", quote=F, sep="\n", row.names=F, col.names=F) # IMPORTANT: remove colname manually

#system("seqkit grep -r -f data/omicron_artifact_inspection.txt data/new/cog_2023-03-31_all_alignment.fa -o data/omicron_aln_inspection.fasta")

system("seqtk subseq data/new/cog_2023-03-31_all_alignment.fa data/omicron_artifact_inspection.txt > data/omicron_aln_inspection.fasta")

# https://www.ncbi.nlm.nih.gov/gene/43740568
# Spike positions nt: 21563 to 25384
# dropout of amplicon76 from 22677 to 23028 (this amplicon covers most of S-gene RBD: codons 372 to 489)
# region at 22786-22974 which contains K417N, N440K, G446S without calls
# https://codon2nucleotide.theo.io/ to find nt positions and subtract 1
# 417 = 22810 to 22812 --> AAT(N), AAG(K) and NNN
# 440 = 22879 to 22881 --> AAG(K), AAT(N), and NNN
# 446 = 22897 to 22899 --> AGT(S), GGT(G), and NNN

system("seqkit subseq -r 22811:22813 data/omicron_aln_inspection.fasta -o data/417_codon.fasta")
system("seqkit subseq -r 22880:22882 data/omicron_aln_inspection.fasta -o data/440_codon.fasta")
system("seqkit subseq -r 22898:22900 data/omicron_aln_inspection.fasta -o data/446_codon.fasta")

s417 <- read.FASTA("data/417_codon.fasta", type="DNA")
s440 <- read.FASTA("data/440_codon.fasta", type="DNA")
s446 <- read.FASTA("data/446_codon.fasta", type="DNA")

dnabin_to_df <- function(x) {
	df <- data.frame(ID=labels(x),
																		#sp=attr(x, "species"),
																		seq=sapply(as.character(x), paste, collapse=""))
	return(df)
}

s417_df <- dnabin_to_df(s417); s417_df$site = 417
s440_df <- dnabin_to_df(s440); s440_df$site = 440
s446_df <- dnabin_to_df(s446); s446_df$site = 446

art_sites <- rbind(s417_df,s440_df,s446_df)
patterns_count <- art_sites %>% group_by(site,seq) %>% summarise(n=n())

aln_sites <- art_sites
md_sites <- artifacts2
md_sites$site <- parse_number(md_sites$mut) 

# count of unique patterns (total: 593976)
# 417: aat -> N (347946), aag -> K (184735), nnn -> X (59932), act -> T (1234)
# 440: aag -> K (348375), aat -> N (185231), nnn -> X (60271)
# 446: agt -> S (180951), ggt -> G (352101), nnn -> X (60384), gtt -> V (336)

#md_sites #md
#aln_sites #aln
# confirm that 'nnn' from alignment is not present in the metadata
aln_sites_ns <- aln_sites[aln_sites$seq == "nnn",]
aln_sites_ns_join <- aln_sites_ns %>% inner_join(md_sites, by=c("ID"="sequence_name", "site"="site"))

aln_sites_mut_417 <- aln_sites[aln_sites$seq == "aat" & aln_sites$site==417,] #aag
aln_sites_mut_417_join <- aln_sites_mut_417 %>% inner_join(md_sites, by=c("ID"="sequence_name", "site"="site"))
aln_sites_mut_440 <- aln_sites[aln_sites$seq == "aag" & aln_sites$site==440,] #aat
aln_sites_mut_440_join <- aln_sites_mut_440 %>% inner_join(md_sites, by=c("ID"="sequence_name", "site"="site"))
aln_sites_mut_446 <- aln_sites[aln_sites$seq == "agt" & aln_sites$site==446,] #ggt
aln_sites_mut_446_join <- aln_sites_mut_446 %>% inner_join(md_sites, by=c("ID"="sequence_name", "site"="site"))
# No overlapping, so sequences with wt muts in aln not reported by metadata --> OK