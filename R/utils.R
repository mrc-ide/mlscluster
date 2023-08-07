libs_load <- c("dplyr","ggplot2", "glue","readr","data.table","ape","stringr")
invisible( lapply(libs_load, library, character.only=TRUE) )

### Remove >99% quantile threshold of Freq_homopl since these sites are clearly outliers
thr_pref <- "threshold_quantile"
path_thresholds <- c(glue("{thr_pref}_0.25"), glue("{thr_pref}_0.5"), glue("{thr_pref}_0.75"), glue("{thr_pref}_1"), glue("{thr_pref}_2"),
																				glue("{thr_pref}_3"), glue("{thr_pref}_4"), glue("{thr_pref}_5"), glue("{thr_pref}_10"), glue("{thr_pref}_25"))

remove_homopl_freq_outliers <- function(path_stats)  {
	clustered_dfs <- quant <- list()
	for(i in 1:length(path_thresholds)) {
		clustered_dfs[[i]] <- read.csv(glue("{path_stats}/{path_thresholds[i]}/clustered_all_df.csv"), header=T)
		clustered_dfs[[i]]$Freq_homopl <- as.integer(clustered_dfs[[i]]$Freq_homopl)
		quant[[i]] <- quantile(clustered_dfs[[i]]$Freq_homopl, probs=seq(0,1,1/100))
		#print(quant[[i]])
		#print(unname(quant[[i]][names(quant[[i]]) == "99%"]))
		clustered_dfs[[i]] <- clustered_dfs[[i]][ clustered_dfs[[i]]$Freq_homopl <= unname(quant[[i]][names(quant[[i]]) == "99%"]), ]
		print(max(clustered_dfs[[i]]$Freq_homopl))
	}
	clustered_dfs
}

#p2_rem_outl <- remove_homopl_freq_outliers("results/02_sc2_root_to_nov2021/")
#p3_rem_outl <- remove_homopl_freq_outliers("results/03_sc2_whole_period/")

include_major_lineage_column <- function(md_df) {
	md_df$major_lineage = md_df$lineage
	#md_df = md_df %>% relocate(major_lineage, .after = lineage) (REMOVED BECAUSE ERRORS oON `polymorphic_site_counts.R`)
	md_df$major_lineage <- ifelse(grepl("^A[A-S]\\.|^AZ\\.", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	#md_df$major_lineage <- ifelse(grepl("^AZ\\.", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^(AY\\.1$|AY\\.10|AY\\.11|AY\\.17|AY\\.19|AY\\.21|AY\\.120|AY\\.121|AY\\.122|AY\\.123|AY\\.124|AY\\.125|AY\\.126|AY\\.127|AY\\.128|AY\\.129|AY\\.13|AY\\.14|AY\\.16|AY\\.20|AY\\.22|AY\\.23|AY\\.24|AY\\.25|AY\\.26|AY\\.27|AY\\.28|AY\\.29|AY\\.3|AY\\.41|AY\\.42|AY\\.43|AY\\.44|AY\\.45|AY\\.46|AY\\.47|AY\\.5|AY\\.6|AY\\.7|AY\\.8|AY\\.9)", md_df$major_lineage), "Delta_other", as.character(md_df$major_lineage))
	#md_df$major_lineage <- ifelse(grepl("^AY\\.43", md_df$major_lineage), "Delta AY.43.*", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^AY\\.4$", md_df$major_lineage), "Delta_AY.4.*", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^AY\\.4.", md_df$major_lineage), "Delta_AY.4.*", as.character(md_df$major_lineage))
	#md_df$major_lineage <- ifelse(grepl("^AY\\.122", md_df$major_lineage), "Delta AY.122.*", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^A\\.", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^A$", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("B\\.1\\.617", md_df$major_lineage), "Delta_other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("B\\.1\\.351", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	#md_df$major_lineage <- ifelse(grepl("B\\.1\\.351", md_df$major_lineage), "Beta_B.1.351", as.character(md_df$major_lineage)) # OLD subdivision
	md_df$major_lineage <- ifelse(grepl("^B\\.1\\.1\\.7$", md_df$major_lineage), "Alpha_B.1.1.7", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B\\.1\\.177", md_df$major_lineage), "EU1_B.1.177", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B$", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B\\.1$", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B\\.1\\.1$", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	
	md_df$major_lineage <- ifelse(grepl("^BA\\.1", md_df$major_lineage), "Omicron_BA.1.*", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^BA\\.2", md_df$major_lineage), "Omicron_BA.2.*", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^BA\\.3|^BA\\.4|^BA\\.5", md_df$major_lineage), "Omicron_BA.*", as.character(md_df$major_lineage))
	#md_df$major_lineage <- ifelse(grepl("^BA\\.4", md_df$major_lineage), "Omicron_BA.4.*", as.character(md_df$major_lineage))
	#md_df$major_lineage <- ifelse(grepl("^BA\\.5", md_df$major_lineage), "Omicron_BA.5.*", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^C|^W|^P\\.2", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	#md_df$major_lineage <- ifelse(grepl("^W", md_df$major_lineage), "W.*", as.character(md_df$major_lineage))
	#md_df$major_lineage <- ifelse(grepl("^P\\.1", md_df$major_lineage), "Gamma_P.1", as.character(md_df$major_lineage)) # OLD subdivision
	md_df$major_lineage <- ifelse(grepl("^P\\.1", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	#md_df$major_lineage <- ifelse(grepl("^P\\.2", md_df$major_lineage), "Zeta P.2", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^(AV\\.|AT\\.|BB\\.|D\\.|G\\.|L\\.|N\\.|M\\.|P\\.|Q\\.|R\\.|S\\.|U\\.|V\\.|Unassigned|Y\\.|Z\\.)", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B\\.1\\.1\\.[1-9]{1,3}", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B\\.1\\.[2-9]{1,3}", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B\\.[2-9]{1,3}", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("B.18|B.1.104|B.1.105|B.1.110.2|B.1.110.3|B.1.111|B.1.115|B.1.117|B.1.12|B.1.128|B.1.13|B.1.131|B.1.142|B.1.146|B.1.147|B.1.149|B.1.153|B.1.157|B.1.160|B.1.160.|B.1.164|B.1.165|B.1.167|B.1.170|B.1.173|B.1.179|B.1.180|B.1.189|B.1.190|B.1.195|B.1.198|B.1.199", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^X", md_df$major_lineage), "Recombinant", as.character(md_df$major_lineage))
	#md_df$major_lineage <- ifelse(grepl("B.18", md_df$major_lineage), "B.*", as.character(md_df$major_lineage))
	#md_df$major_lineage <- ifelse(grepl("^X", md_df$major_lineage), "X.*", as.character(md_df$major_lineage))
	
	## UNCOMMENT IF WANTS PLOT
	# md_df_freq <- md_df %>% count(major_lineage, sort=T)
	# 
	# ggplot(data=md_df_freq, aes(x=n, y=major_lineage)) + geom_bar(stat="identity")
	# ggsave(glue("rds/major_lineage_freqs.png"), width=8, height=5, dpi=600, bg="white")
	# 
	md_df_res = subset(md_df, major_lineage != "Recombinant")
	md_df_res <- md_df_res[!is.na(md_df_res$major_lineage),]
	# 
	# md_df_res_freq <- md_df_res %>% count(major_lineage, sort=T)
	# ggplot(data=md_df_res_freq, aes(x=n, y=major_lineage)) + geom_bar(stat="identity")
	# ggsave(glue("rds/major_lineage_freqs_after_exclusions.png"), width=8, height=5, dpi=600, bg="white")
	
	return(md_df_res)
}

# codon2nucleotide_conversion to handle potential artifacts by checking alignment NNNs (based on https://github.com/theosanderson/Codon2Nucleotide/blob/main/src/App.js)
#nonsyn_p2 <- read.csv("stat_results/period2/non-syn_freq_g10.csv", header=T)
nonsyn_p3 <- read.csv("stat_results/period3/non-syn_freq_g10.csv", header=T)

#syn_p2 <- read.csv("stat_results/period2/syn_freq_g10.csv", header=T)
syn_p3 <- read.csv("stat_results/period3/syn_freq_g10.csv", header=T)

# Already done below, no need to do again (begin)
sc2_md_curated <- readRDS("rds/sc2_md_curated_WITHOUT_MISSING.rds")
system("mkdir -p data/insp/")
write.csv(sc2_md_curated$sequence_name, file="data/insp/incl_seqs.txt", quote=F, sep="\n", row.names=F, col.names=F)
# Important to remove header manually from file above
# Subset 2023 aln to include only sequences included in our study
system("seqtk subseq data/new/cog_2023-03-31_all_alignment.fa data/insp/incl_seqs.txt > data/insp/incl_seqs.fasta")
# (end)

codon2nucleotide_conversion <- function(df, non_syn=T) {
	
	# load coordinates
	coords_all <- read.csv("config/codon2nt_conv_all.csv", header=T)
	#coords_nsps <- read.csv("config/codon2nt_conv_nsps.csv", header=T)
	
	# Remove duplicates
	df <- df[!duplicated(df$defining_mut),]
	
	# Separate e.g protein=="S", ref_pos_alt=="E484K" and aa_pos==484
	df$protein <- sub("\\:.*", "", df$defining_mut)
	df$ref_pos_alt <- sub('.*:', "", df$defining_mut)
	df$aa_pos <- parse_number(df$ref_pos_alt)
	
	
	if(non_syn) {
		# join with coord df
		df_join <- df %>% inner_join(coords_all, by="protein")
		df_join$start_nt <- (df_join$start + (df_join$aa_pos - 1) * 3) # NOT DOING THE FOLLOWING ANYMORE: subtracting 1 in the end because COG-UK aln is 29902 nt (not 29903), but not sure this will work for every protein
		df_join$start_nt <- ifelse(df_join$protein == "ORF1AB" & (df_join$start_nt > 13468 & df_join$start_nt < 21557), df_join$start_nt-1, df_join$start_nt)
		#View(df_join)
		df_join$end_nt <- df_join$start_nt + 2
		system("mkdir -p data/insp/non-syn/")
		for(n in 1:nrow(df_join)) {
			print(glue("{df_join[n,7]}:{df_join[n,17]}"))
			system(glue("seqkit subseq -r {df_join[n,21]}:{df_join[n,22]} data/insp/incl_seqs.fasta -o data/insp/non-syn/{df_join[n,7]}:{df_join[n,17]}_nt{df_join[n,21]}-{df_join[n,22]}.fasta"))
		}
	} else { #SYN
		# identify from which region the mutation come from since only SYNSNP shown
		setDT(df); setDT(coords_all)
		df_join <- df[coords_all, on=c("protein==name_cog","aa_pos>=start","aa_pos<=end"), protein := protein]
		# in this case we will extract only the nt site and not the triplet
		system("mkdir -p data/insp/syn/")
		for(n in 1:nrow(df_join)) {
			system(glue("seqkit subseq -r {df_join[n,17]}:{df_join[n,17]} data/insp/incl_seqs.fasta -o data/insp/syn/{df_join[n,7]}:{df_join[n,17]}.fasta"))
		}
	}
	View(df_join)
}

#codon2nucleotide_conversion(nonsyn_p2, non_syn=T)
codon2nucleotide_conversion(nonsyn_p3, non_syn=T) #363

#codon2nucleotide_conversion(syn_p2, non_syn=F)
codon2nucleotide_conversion(syn_p3, non_syn=F) #329

dnabin2df <- function(x) {
	df <- data.frame(sequence_name=labels(x),
																		#sp=attr(x, "species"),
																		seq=sapply(as.character(x), paste, collapse=""))
	return(df)
}

# WORKS, BUT MEMORY ISSUES
load_aln_fill_missing_sites <- function(path_ns, path_s) {
	if(dir.exists(path_ns) & dir.exists(path_s)) {
		# extract path for all files
		files_ns <- list.files(path_ns, full.names=T, pattern="*.fasta")
		files_s <- list.files(path_s, full.names=T, pattern="*.fasta")
		#print(files_s)
		codon_seqs <- list()
		nt_seqs <- list()
		print("LOADING NON-SYN ALNS")
		for(i in 1:length(files_ns)) {
			# NON-SYN
			print(i)
			codon_seqs[[i]] <- read.FASTA(glue("{files_ns[i]}"), type="DNA")
			codon_seqs[[i]] <- dnabin2df(codon_seqs[[i]])
			rownames(codon_seqs[[i]]) <- NULL
			#print(head(codon_seqs[[i]]))
		}
		
		print("LOADING SYN ALNS")
		for(i in 1:length(files_s)) { #TODO change here when all working to do for all nts: length(files_s)
			# SYN
			print(i)
			nt_seqs[[i]] <- read.FASTA(glue("{files_s[i]}"), type="DNA")
			#print(files_s[[i]])
			nt_seqs[[i]] <- dnabin2df(nt_seqs[[i]])
			rownames(nt_seqs[[i]]) <- NULL
			#print(head(nt_seqs[[i]]))
		}
		
		# NON-SYN
		print("STARTING NON-SYN MATCHING")
		codon_seqs_df_nnn <- list()
		files_basen <- c()
		files_basen_df <- list()
		dfs_join <- intersect_seqnames <- adj_df1 <- list()
		for(j in 1:length(codon_seqs)) {
			print(j)
			# get only seqs with "nnn"
			codon_seqs_df_nnn[[j]] <- codon_seqs[[j]][codon_seqs[[j]]$seq == "nnn",]
			files_basen[j] <- basename(files_ns[j])
			files_basen[j] <- gsub("_.*","",files_basen[j])
			#print(files_basen[j])
			files_basen_df[[j]] <- data.frame(prot_site=files_basen[[j]])
			dfs_join[[j]] <- nonsyn_p3 %>% inner_join(files_basen_df[[j]], by="prot_site")
			dfs_join[[j]]$anc <- str_sub(dfs_join[[j]]$defining_mut, start=1, end=-2)
			intersect_seqnames[[j]] <- sc2_md_curated %>% inner_join(codon_seqs_df_nnn[[j]], by="sequence_name")
			adj_df1[[j]] <- sc2_md_curated %>% inner_join(intersect_seqnames[[j]], by="sequence_name") #left_join
			adj_df1[[j]] <- adj_df1[[j]] %>% mutate(nnn_mutations = if_else(seq == 'nnn', paste0("|",unique(dfs_join[[j]]$anc),"X"), mutations.x))
			adj_df1[[j]] <- adj_df1[[j]] %>% select(sequence_name, sample_date.x, lineage.x, major_lineage.x, region.x, mutations.x, seq, nnn_mutations)
			colnames(adj_df1[[j]]) <- c("sequence_name","sample_date","lineage","major_lineage","region","mutations","seq","n_mutations")
			#View(adj_df1[[j]])
		}
		
		# SYN
		nt_seqs_df_n <- list()
		files_basen2 <- c()
		files_basen2_df <- list()
		
		coords_all <- read.csv("config/codon2nt_conv_all.csv", header=T)
		setDT(syn_p3); setDT(coords_all)
		syn_p3 <- syn_p3[coords_all, on=c("protein==name_cog","site>=start","site<=end"), protein := protein]
		syn_p3$prot_site <- paste0(syn_p3$protein,":",syn_p3$site)
		#syn_p3$defining_mut2 <- gsub("SYNSNP:","synSNP:",syn_p3$defining_mut)
		syn_p3$defining_mut2 <- gsub("SYNSNP:","",syn_p3$defining_mut)
		syn_p3$defining_mut3 <- paste0(syn_p3$protein,":",syn_p3$defining_mut2)
		
		print("STARTING SYN MATCHING")
		
		dfs_join2 <- intersect_seqnames2 <- adj_df1s <- list()
		for(k in 1:length(nt_seqs)) {
			# get only seqs with "n"
			print(k)
			nt_seqs_df_n[[k]] <- nt_seqs[[k]][nt_seqs[[k]]$seq == "n",]
			files_basen2[k] <- basename(files_s[k])
			files_basen2[k] <- gsub("\\.fasta","",files_basen2[k])
			files_basen2_df[[k]] <- data.frame(prot_site=files_basen2[[k]])
			
			dfs_join2[[k]] <- syn_p3 %>% inner_join(files_basen2_df[[k]], by="prot_site")
			dfs_join2[[k]]$anc <- str_sub(dfs_join2[[k]]$defining_mut3, start=1, end=-2) #defining_mut3
			#print("BEFORE")
			#print(head(dfs_join2[[k]]))
			intersect_seqnames2[[k]] <- sc2_md_curated %>% inner_join(nt_seqs_df_n[[k]], by="sequence_name")
			adj_df1s[[k]] <- sc2_md_curated %>% inner_join(intersect_seqnames2[[k]], by="sequence_name") #left_join
			#dfs_join2[[k]]$anc <- gsub(".*:","synSNP:",dfs_join2[[k]]$anc)
			dfs_join2[[k]]$anc <- gsub(".*:","",dfs_join2[[k]]$anc)
			dfs_join2[[k]]$anc <- paste0("synSNP:",dfs_join2[[k]]$anc)
			#print("AFTER")
			#print(head(dfs_join2[[k]]))
			adj_df1s[[k]] <- adj_df1s[[k]] %>% mutate(n_mutations = if_else(seq == 'n', paste0("|",unique(dfs_join2[[k]]$anc),"N"), mutations.x))
			adj_df1s[[k]] <- adj_df1s[[k]] %>% select(sequence_name, sample_date.x, lineage.x, major_lineage.x, region.x, mutations.x, seq, n_mutations)
			colnames(adj_df1s[[k]]) <- c("sequence_name","sample_date","lineage","major_lineage","region","mutations","seq","n_mutations")
			#View(adj_df1s[[k]])
		}
		
		adj_df2ns <- rbindlist(adj_df1)
		adj_df2s <- rbindlist(adj_df1s)
		adj_df2 <- rbind(adj_df2ns, adj_df2s)
		#View(adj_df2)
		adj_df2 <- adj_df2 %>% group_by(sequence_name) %>% summarize(miss_mutations = paste0(na.omit(n_mutations), collapse="")) %>% ungroup()
		#View(adj_df2)
		adj_df3 <- sc2_md_curated %>% left_join(adj_df2, by="sequence_name")
		adj_df3$muts_joined <- paste0(adj_df3$mutations, adj_df3$miss_mutations) #na.omit
		#View(adj_df3)
		adj_df3 <- adj_df3 %>% select(sequence_name, sample_date, lineage, major_lineage, region, muts_joined)
		colnames(adj_df3) <- c("sequence_name", "sample_date", "lineage", "major_lineage", "region", "mutations")
		adj_df3$mutations <- gsub("NA","", as.character(adj_df3$mutations))
		# View(adj_df3)
	}else {
		print("One or both the specified paths do not exist!")
	}
	return(adj_df3)
}

incl_missing_muts <- load_aln_fill_missing_sites("data/insp/non-syn", "data/insp/syn")

#sc2_md_curated <- incl_missing_muts
saveRDS(incl_missing_muts, file="rds/sc2_md_curated.rds")

# incl_missing_muts[incl_missing_muts$sequence_name == "England/ALDP-13A2064/2021",]$mutations #both syn and non-syn
# incl_missing_muts[incl_missing_muts$sequence_name == "England/ALDP-14A91A5/2021",]$mutations #only syn
# incl_missing_muts[incl_missing_muts$sequence_name == "England/ALDP-14ADC37/2021",]$mutations #only non-syn
