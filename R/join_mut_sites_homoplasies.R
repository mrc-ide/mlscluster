library(glue)
library(stringr)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)

#load(file="rds/env.RData")
seq_muts_df_counts <- readRDS("rds/polymorphic_site_counts.rds")

pos_sel_sites <- read.csv("config/positive_selected_sites_20221006.tsv", sep="\t", header=TRUE)
pos_sel_sites$prot_site <- paste0(pos_sel_sites$protein,":",pos_sel_sites$site)

# Paths containing results for each time period
c_path <- "results"
path_periods <- c(glue("{c_path}/01_sc2_root_to_dec2020"), glue("{c_path}/02_sc2_jan2021_to_may2021"),glue("{c_path}/03_sc2_jun2021_to_dec2021"), 
																		glue("{c_path}/04_sc2_jan2022_to_apr2022"), glue("{c_path}/05_sc2_whole_period"))

# Paths containing results for each quantile threshold
thr_pref <- "threshold_quantile"
thr_path <- c(glue("{thr_pref}_0.25"), glue("{thr_pref}_0.5"), glue("{thr_pref}_0.75"), glue("{thr_pref}_1"), glue("{thr_pref}_2"),
														glue("{thr_pref}_3"), glue("{thr_pref}_4"), glue("{thr_pref}_5"), glue("{thr_pref}_10"), glue("{thr_pref}_25")) #, glue("{thr_pref}_50")

table_names_periods <- c(1,2,3,4,5)
table_names_thresholds <- c(0.25,0.5,0.75,1,2,3,4,5,10,25)
#table_combs <- paste(rep(table_names_periods, each = length(table_names_thresholds)), table_names_thresholds, sep = "_")
#table_combs <- paste(rep(table_names_thresholds, length.out = length(table_names_periods)), table_names_periods, sep = "_")
table_combs <- do.call(paste, c(expand.grid(table_names_periods, table_names_thresholds), sep="_"))
names(table_combs) <- seq(1:50)

annot_gen_range_positions <- read.csv("config/aa_ranges_stat_model.tsv",sep="\t")
annot_gen_length_positions <- read.csv("config/aa_length_stat_model.tsv",sep="\t")
sink(file = "rds/lm_output_POISSON_SPLITTED_ALL.txt")
clustered_dfs <- matrix(data=list(), nrow=length(path_periods), ncol=length(thr_path))
joined_mut_sites_clustered <- matrix(data=list(), nrow=length(path_periods), ncol=length(thr_path))
joined_mut_sites_clustered_syn <- joined_mut_sites_clustered_non_syn <- joined_mut_sites_clustered
var_all <- mean_all <- poisson_models_syn <- poisson_models_non_syn <- lr_models_syn <- lr_models_non_syn <- csv_cp <- csv_cp2 <- matrix(data=list(), nrow=length(path_periods), ncol=length(thr_path))
for(i in 1:length(path_periods)) { #length(path_periods)
	for(j in 1:length(thr_path)) {
		print("====================")
		print(glue("[{i},{j}]"))
		clustered_dfs[[i,j]] <- read.csv(glue("{path_periods[i]}/{thr_path[j]}/clustered_all_df.csv"), header=T)
		#View(clustered_dfs[[i,j]])
		
		joined_mut_sites_clustered[[i,j]] <- clustered_dfs[[i,j]] %>% full_join(seq_muts_df_counts, by=c("defining_mut"="mut_and_type")) #inner_join
		# Extract genomic index position
		rgx_genomic_idx <- regexpr(':[A-Z]{1}[0-9]{1,5}[A-Z*]{1}', joined_mut_sites_clustered[[i,j]]$defining_mut)
		comm_coord <- rep(NA,length(joined_mut_sites_clustered[[i,j]]$defining_mut))
		comm_coord[rgx_genomic_idx != -1] <- str_sub( regmatches(joined_mut_sites_clustered[[i,j]]$defining_mut, rgx_genomic_idx), 3, -2)
		comm_coord <- as.integer(comm_coord)
		comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("genomic_index")
		joined_mut_sites_clustered[[i,j]] <- cbind(joined_mut_sites_clustered[[i,j]], comm_coord_df)
		#joined_mut_sites_clustered[[i,j]]$genomic_index <- as.integer(joined_mut_sites_clustered[[i,j]]$genomic_index)
		joined_mut_sites_clustered[[i,j]] <- joined_mut_sites_clustered[[i,j]] %>% select(defining_mut, Freq_homopl, is_clustered, syn_non_syn.x, protein.x, aa_length.x, indep_found_pos_selection, genomic_index) #s_mut_region_interest
		colnames(joined_mut_sites_clustered[[i,j]]) <- c("defining_mut","Freq_homopl","is_clustered","syn","protein","aa_length","indep_found_pos_selection","genomic_index") #s_mut_region_interest
		
		joined_mut_sites_clustered[[i,j]]$protein <- sub("\\:.*", "", joined_mut_sites_clustered[[i,j]]$defining_mut)
		
		setDT(joined_mut_sites_clustered[[i,j]]); setDT(annot_gen_range_positions)
		joined_mut_sites_clustered[[i,j]] <- joined_mut_sites_clustered[[i,j]][annot_gen_range_positions, on=c("protein==name_cog_md","genomic_index>=start","genomic_index<=end"), spec_regions := new_prot_name]
		joined_mut_sites_clustered[[i,j]] <- joined_mut_sites_clustered[[i,j]][annot_gen_length_positions, on=c("spec_regions==new_prot_name"), aa_length := length]
		#View(joined_mut_sites_clustered[[i,j]])
		#homopl_df[regions_s, on=c("s_mut_coord>=start", "s_mut_coord<=end"), s_mut_region_interest := region]
		# joined_mut_sites_clustered[[i,j]][joined_mut_sites_clustered[[i,j]] == "syn"] <- 1
		# joined_mut_sites_clustered[[i,j]][joined_mut_sites_clustered[[i,j]] == "non-syn"] <- 0
		joined_mut_sites_clustered[[i,j]]$syn <- ifelse(joined_mut_sites_clustered[[i,j]]$syn == as.character("non-syn"), 0, 1)
		joined_mut_sites_clustered[[i,j]]$indep_found_pos_selection <- ifelse(joined_mut_sites_clustered[[i,j]]$indep_found_pos_selection == as.character("yes"), 1, 0)
		# Add zeros to Freq_homopl and is_clustered if polymorphic site==NA (not found on homoplasy tables)
		joined_mut_sites_clustered[[i,j]] <- joined_mut_sites_clustered[[i,j]] %>% mutate(Freq_homopl = ifelse(is.na(Freq_homopl), 0, Freq_homopl))
		joined_mut_sites_clustered[[i,j]] <- joined_mut_sites_clustered[[i,j]] %>% mutate(is_clustered = ifelse(is.na(is_clustered), 0, is_clustered))
		# Add syn and non-syn integer values
		joined_mut_sites_clustered[[i,j]]$syn <- ifelse(joined_mut_sites_clustered[[i,j]]$protein == "SYNSNP", 1, 0)
		# Add indep_found_under_selection flag
		
		rgx_scpss <- regexpr(':[A-Z]{1}[0-9]{1,5}[A-Z*]{1}', joined_mut_sites_clustered[[i,j]]$defining_mut)
		comm_coord <- rep(NA,length(joined_mut_sites_clustered[[i,j]]$defining_mut))
		# Mutation coordinate of homoplasies
		comm_coord[rgx_scpss != -1] <- str_sub( regmatches(joined_mut_sites_clustered[[i,j]]$defining_mut, rgx_scpss), 3, -2)
		comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("site")
		joined_mut_sites_clustered[[i,j]] <- cbind(joined_mut_sites_clustered[[i,j]], comm_coord_df)
		joined_mut_sites_clustered[[i,j]]$prot_site <- paste0(joined_mut_sites_clustered[[i,j]]$protein,":",joined_mut_sites_clustered[[i,j]]$site)
		joined_mut_sites_clustered[[i,j]]$indep_found_pos_selection <- ifelse(joined_mut_sites_clustered[[i,j]]$prot_site %in% pos_sel_sites$prot_site, 1, 0)
		#View(joined_mut_sites_clustered[[i,j]])
		
		# See if variance and mean approx equal: one of poisson regression assumptions
		# var_all[[i,j]] <- var(joined_mut_sites_clustered[[i,j]]$Freq, na.rm=T)
		# print(glue("Variance = {var_all[[i,j]]}"))
		# mean_all[[i,j]] <- mean(joined_mut_sites_clustered[[i,j]]$Freq, na.rm=T)
		# print(glue("Mean = {mean_all[[i,j]]}"))
		#hist(joined_mut_sites_clustered[[i,j]]$genomic_index)
		
		# joined_mut_sites_clustered[[i,j]]$is_clustered <- as.factor( joined_mut_sites_clustered[[i,j]]$is_clustered ) #, levels=c(0,1), labels=c('0','1')
		# joined_mut_sites_clustered[[i,j]]$indep_found_pos_selection <- as.factor( joined_mut_sites_clustered[[i,j]]$indep_found_pos_selection ) #, levels=c(0,1), labels=c('0','1')
		#joined_mut_sites_clustered[[i,j]]$protein <- as.factor(joined_mut_sites_clustered[[i,j]]$protein)
		#joined_mut_sites_clustered[[i,j]]$aa_length <- as.factor(joined_mut_sites_clustered[[i,j]]$aa_length)
		#joined_mut_sites_clustered[[i,j]]$syn <- as.factor(joined_mut_sites_clustered[[i,j]]$syn)
		#joined_mut_sites_clustered[[i,j]]$s_mut_region_interest <- as.factor(joined_mut_sites_clustered[[i,j]]$s_mut_region_interest)
		joined_mut_sites_clustered[[i,j]]$spec_regions <- as.factor(joined_mut_sites_clustered[[i,j]]$spec_regions)
		# Normalised counts not being used for regressions for now
		joined_mut_sites_clustered[[i,j]]$Freq_homopl_norm_size <- joined_mut_sites_clustered[[i,j]]$Freq_homopl / joined_mut_sites_clustered[[i,j]]$aa_length
		
		# Split into syn and non-syn
		system("mkdir -p rds/stat_ready_csvs_all_sites/")
		system("mkdir -p rds/stat_ready_csvs_only_clustered/")
		joined_mut_sites_clustered_syn[[i,j]] <- joined_mut_sites_clustered[[i,j]][joined_mut_sites_clustered[[i,j]]$syn == 1,]
		csv_cp[[i,j]] <- joined_mut_sites_clustered_syn[[i,j]][ order(joined_mut_sites_clustered_syn[[i,j]]$Freq_homopl, decreasing=T) ,]
		#csv_cp[[i,j]] <- joined_mut_sites_clustered_syn[[i,j]][ order(joined_mut_sites_clustered_syn[[i,j]]$Freq_homopl_norm_size, decreasing=T) ,]
		write.csv(csv_cp[[i,j]], glue("rds/stat_ready_csvs_all_sites/df_all_sites_P{i}_{thr_path[j]}_SYN.csv"), quote=F, row.names=F)
		csv_cp[[i,j]] <- csv_cp[[i,j]][(as.integer(csv_cp[[i,j]]$Freq_homopl) > 0) & (as.integer(csv_cp[[i,j]]$is_clustered)==1), ]
		#csv_cp[[i,j]] <- csv_cp[[i,j]][(as.integer(csv_cp[[i,j]]$Freq_homopl_norm_size) > 0) & (as.integer(csv_cp[[i,j]]$is_clustered)==1), ]
	 write.csv(csv_cp[[i,j]], glue("rds/stat_ready_csvs_only_clustered/df_only_clustered_P{i}_{thr_path[j]}_SYN.csv"), quote=F, row.names=F)
		joined_mut_sites_clustered_non_syn[[i,j]] <- joined_mut_sites_clustered[[i,j]][joined_mut_sites_clustered[[i,j]]$syn == 0,]
		csv_cp2[[i,j]] <- joined_mut_sites_clustered_non_syn[[i,j]][ order(joined_mut_sites_clustered_non_syn[[i,j]]$Freq_homopl, decreasing=T) ,]
		#csv_cp2[[i,j]] <- joined_mut_sites_clustered_non_syn[[i,j]][ order(joined_mut_sites_clustered_non_syn[[i,j]]$Freq_homopl_norm_size, decreasing=T) ,]
		write.csv(csv_cp2[[i,j]], glue("rds/stat_ready_csvs_all_sites/df_all_sites_P{i}_{thr_path[j]}_NON_SYN.csv"), quote=F, row.names=F)
		csv_cp2[[i,j]] <- csv_cp2[[i,j]][(as.integer(csv_cp2[[i,j]]$Freq_homopl) > 0) & (as.integer(csv_cp2[[i,j]]$is_clustered)==1), ]
		#csv_cp2[[i,j]] <- csv_cp2[[i,j]][(as.integer(csv_cp2[[i,j]]$Freq_homopl_norm_size) > 0) & (as.integer(csv_cp2[[i,j]]$is_clustered)==1), ]
		write.csv(csv_cp2[[i,j]], glue("rds/stat_ready_csvs_only_clustered/df_only_clustered_P{i}_{thr_path[j]}_NON_SYN.csv"), quote=F, row.names=F)
		
		print(glue("{path_periods[i]}/{thr_path[j]}/"))

		# Code for non-normalised freqs (poisson model applied for counts)
		poisson_models_syn[[i,j]] <- glm(Freq_homopl ~ is_clustered + spec_regions + aa_length + indep_found_pos_selection, data=joined_mut_sites_clustered_syn[[i,j]], na.action=na.omit, family=poisson(link = "log")) #defining_mut + genomic_index +s_mut_region_interest +syn
		print("Summary poisson (SYN): ")
		print(summary(poisson_models_syn[[i,j]]))

		poisson_models_non_syn[[i,j]] <- glm(Freq_homopl ~ is_clustered + spec_regions + aa_length + indep_found_pos_selection, data=joined_mut_sites_clustered_non_syn[[i,j]], na.action=na.omit, family=poisson(link = "log")) #defining_mut + genomic_index +s_mut_region_interest +syn
		print("Summary poisson (NON-SYN): ")
		print(summary(poisson_models_non_syn[[i,j]]))
		
		# gamma_models_syn[[i,j]] <- glm(Freq_homopl_norm_size ~ is_clustered + spec_regions + aa_length + indep_found_pos_selection, data=joined_mut_sites_clustered_syn[[i,j]], na.action=na.omit, family=gamma(link="log")) #defining_mut + genomic_index +s_mut_region_interest +syn
		# print("Summary gamma (SYN): ")
		# print(summary(gamma_models_syn[[i,j]]))
		# 
		# gamma_models_non_syn[[i,j]] <- glm(Freq_homopl_norm_size ~ is_clustered + spec_regions + aa_length + indep_found_pos_selection, data=joined_mut_sites_clustered_non_syn[[i,j]], na.action=na.omit, family=gamma(link="log")) #defining_mut + genomic_index +s_mut_region_interest +syn
		# print("Summary gamma (NON-SYN): ")
		# print(summary(gamm_models_non_syn[[i,j]]))

		# Not using anymore the code below for logistic regression
		# lr_models_syn[[i,j]] <- glm(is_clustered==1 ~ Freq_homopl:aa_length + spec_regions + indep_found_pos_selection, data=joined_mut_sites_clustered_syn[[i,j]], na.action=na.omit, family=binomial(link="logit"))
		# print("Summary logistic (SYN): ")
		# print(summary(lr_models_syn[[i,j]]))
		# 
		# lr_models_non_syn[[i,j]] <- glm(is_clustered==1 ~ Freq_homopl:aa_length + spec_regions + indep_found_pos_selection, data=joined_mut_sites_clustered_non_syn[[i,j]], na.action=na.omit, family=binomial(link="logit"))
		# print("Summary logistic (NON-SYN): ")
		# print(summary(lr_models_non_syn[[i,j]]))

		print("====================")
	}
	rownames(csv_cp) <- table_names_periods
	colnames(csv_cp) <- table_names_thresholds
	rownames(csv_cp2) <- table_names_periods
	colnames(csv_cp2) <- table_names_thresholds
}
sink(file = NULL)
#save.image("rds/env.RData")

df_syn_clustered_homopl <- rbindlist(csv_cp, use.names=T, idcol="period_threshold")
df_non_syn_clustered_homopl <- rbindlist(csv_cp2, use.names=T, idcol="period_threshold")
# indices 1 to 5 --> P1_t0.25, P2_t0.25 ...
# indices 6 to 10 --> P1_t0.5, P2_t0.5 ...

# Change period_threshold to match correct period (1 to 5) and threshold (0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 10, 25)
change_ids <- function(df) {
	for(i in 1:nrow(df)) {
		if(df[[i,1]] %in% names(table_combs)) {
			idx <- match(df[[i,1]],names(table_combs))
			df[[i,1]] <- table_combs[idx]
		}
	}
	return(df)
}

system("mkdir -p rds/only_clustered_ordered_homopl/")
# SYN
df_syn_clustered_homopl <- change_ids(df_syn_clustered_homopl)
df_syn_clustered_homopl <- df_syn_clustered_homopl[df_syn_clustered_homopl$Freq_homopl >= 5,] # Getting only homoplasies happening at least 5 times
df_syn_clustered_homopl <- df_syn_clustered_homopl[order(df_syn_clustered_homopl$defining_mut, decreasing=T),]
hist(df_syn_clustered_homopl$Freq_homopl)
write.csv(df_syn_clustered_homopl, glue("rds/only_clustered_ordered_homopl/syn_clustered_homopl_ordered.csv"), quote=F, row.names=F)
# NON-SYN
df_non_syn_clustered_homopl <- change_ids(df_non_syn_clustered_homopl)
df_non_syn_clustered_homopl <- df_non_syn_clustered_homopl[df_non_syn_clustered_homopl$Freq_homopl >= 5,]
df_non_syn_clustered_homopl <- df_non_syn_clustered_homopl[order(df_non_syn_clustered_homopl$defining_mut, decreasing=T),]
hist(df_non_syn_clustered_homopl$Freq_homopl)
write.csv(df_non_syn_clustered_homopl, glue("rds/only_clustered_ordered_homopl/non_syn_clustered_homopl_ordered.csv"), quote=F, row.names=F)

# Counts number of rows for each index
df_syn_clustered_homopl_counts <- df_syn_clustered_homopl %>% count(period_threshold)
df_non_syn_clustered_homopl_counts <- df_non_syn_clustered_homopl %>% count(period_threshold)

plot_distr_freqs_spec_regions <- function(df,aa_change,remove_if_x_rows_period_thr=10) {
	system(glue("mkdir -p rds/plot_mut_mean_freqs_{aa_change}"))
	# Split for each period+threshold
	df_splitted <- split(df, df$period_threshold)
	# Remove period+threshold with less than 10 rows
	for(i in 1:length(df_splitted)) {
		if(nrow(df_splitted[[i]]) < remove_if_x_rows_period_thr) {
			df_splitted[i] <- 0
		}
	}
	df_splitted <- Filter(is.data.frame, df_splitted)
	
	for(x in 1:length(df_splitted)) {
		# Normalise frequency of homoplasy by size of genomic region
		df_splitted[[x]]$Freq_homopl_norm_size <- df_splitted[[x]]$Freq_homopl / df_splitted[[x]]$aa_length
	}
	
	system(glue("mkdir -p rds/plot_mut_counts_{aa_change}"))
	df_splitted_mean <- df_splitted_mean_norm_size <- df_splitted
	# Mean of freqs per genomic region
	for(j in 1:length(df_splitted)) {
		df_splitted_mean[[j]] <- df_splitted_mean[[j]] %>% group_by(spec_regions) %>% summarise(mean_freq_homopl=mean(Freq_homopl))
		df_splitted_mean_norm_size[[j]] <- df_splitted_mean_norm_size[[j]] %>% group_by(spec_regions) %>% summarise(mean_freq_homopl_norm=mean(Freq_homopl_norm_size))
		if(aa_change=="syn") {
			df_splitted_mean[[j]]$spec_regions <- sub('.*:', "", df_splitted_mean[[j]]$spec_regions)
			df_splitted_mean_norm_size[[j]]$spec_regions <- sub('.*:', "", df_splitted_mean_norm_size[[j]]$spec_regions)
		}
		# ggplot(data=df_splitted_mean[[j]], aes(x=spec_regions, y=mean_freq_homopl, fill=spec_regions)) + geom_bar(stat="identity") + labs(x="Genomic region",y="Mean homoplasy Freq") + theme_classic() + theme(legend.position="none", axis.text.x=element_text(size=5))
		# ggsave(glue("rds/plot_mut_mean_freqs_{aa_change}/{names(df_splitted_mean)[j]}_mean_freq.png"), width=8, height=6, dpi=600, bg="white")
		# 
		# ggplot(data=df_splitted_mean_norm_size[[j]], aes(x=spec_regions, y=mean_freq_homopl_norm, fill=spec_regions)) + geom_bar(stat="identity") + labs(x="Genomic region",y="Mean homoplasy Freq normalised by genomic region") + theme_classic() + theme(legend.position="none", axis.text.x=element_text(size=5))
		# ggsave(glue("rds/plot_mut_mean_freqs_{aa_change}/{names(df_splitted_mean)[j]}_mean_freq_norm.png"), width=8, height=6, dpi=600, bg="white")
	}
	#View(df_splitted_mean[[1]])
	
	# Count of homoplasies per genomic region
	df_splitted_count_muts <- df_splitted
	for(k in 1:length(df_splitted_count_muts)) {
		df_splitted_count_muts[[k]] <- df_splitted_count_muts[[k]] %>% group_by(spec_regions) %>% summarise(count_muts=n())
		if(aa_change=="syn") {
			df_splitted_count_muts[[k]]$spec_regions <- sub('.*:', "", df_splitted_count_muts[[k]]$spec_regions)
		}
		# ggplot(data=df_splitted_count_muts[[k]], aes(x=spec_regions, y=count_muts, fill=spec_regions)) + geom_bar(stat="identity") + labs(x="Genomic region",y="# mutations") + theme_classic() + theme(legend.position="none", axis.text.x=element_text(size=5))
		# ggsave(glue("rds/plot_mut_counts_{aa_change}/{names(df_splitted_count_muts)[k]}_mut_counts.png"), width=8, height=6, dpi=600, bg="white")
	}
	
	#View(df_splitted_count_muts[[1]])
	return(df_splitted)
}

df_syn_clustered_homopl_plots <- plot_distr_freqs_spec_regions(df_syn_clustered_homopl,"syn",10)
df_non_syn_clustered_homopl_plots <- plot_distr_freqs_spec_regions(df_non_syn_clustered_homopl,"non-syn",10)

# For tables, combination period_threshold needs to have >1 row
df_syn_clustered_homopl_tables <- plot_distr_freqs_spec_regions(df_syn_clustered_homopl,"syn",2)
df_non_syn_clustered_homopl_tables <- plot_distr_freqs_spec_regions(df_non_syn_clustered_homopl,"syn",2)

get_most_common_homopls_period_interest <- function(df,period_interest) {
	df_homopl_binded <- rbindlist(df)
	# Split period_threshold into separated period and threshold columns
	df_homopl_binded <- df_homopl_binded %>% separate(period_threshold, c("period","threshold"), sep="_")
	# Get only period 5 (full period: end of 2019 to Apr 2022)
	df_homopl<- df_homopl_binded %>% filter(period==period_interest)
	# Make sure that columns are numeric
	df_homopl$period <- as.integer(df_homopl$period)
	df_homopl$threshold <- as.numeric(df_homopl$threshold)
	df_homopl$Freq_homopl <- as.numeric(df_homopl$Freq_homopl)

	filt_high_thresholds <- filt_low_thresholds <- df_homopl
	# Get TOP-50 frequent homoplasies when threshold >= 10 (n=2)
	filt_high_thresholds <- filt_high_thresholds %>% filter(threshold>=10) %>% slice_max(n=50, order_by=Freq_homopl) %>% arrange(period,threshold,desc(Freq_homopl))
	# Get TOP-50 frequent homoplasies when threshold <= 5 (n=8)
	filt_low_thresholds <- filt_low_thresholds %>% filter(threshold<=5) %>% slice_max(n=50, order_by=Freq_homopl) %>% arrange(period,threshold,desc(Freq_homopl))

	return(list(filt_high_thresholds, filt_low_thresholds))
}

# TOP50 for period 5 (end of 2019 to apr 2022)
# >= 10% quantile threshold
most_common_homopl_p5_syn <- get_most_common_homopls_period_interest(df_syn_clustered_homopl_tables, 5)
# <= 5% quantile threshold
most_common_homopl_p5_non_syn <- get_most_common_homopls_period_interest(df_non_syn_clustered_homopl_tables, 5)

# TOP50 for each threshold
get_most_common_homopls_each_threshold <- function(df) {
	df_homopl_binded <- rbindlist(df)
	df_homopl_binded <- df_homopl_binded %>% separate(period_threshold, c("period","threshold"), sep="_")
	df_homopl_binded$period <- as.integer(df_homopl_binded$period)
	df_homopl_binded$threshold <- as.numeric(df_homopl_binded$threshold)
	df_homopl_binded$Freq_homopl <- as.integer(df_homopl_binded$Freq_homopl)
	df_homopl_each_threshold <- df_homopl_binded %>% group_by(threshold) %>% slice_max(n=50, order_by=Freq_homopl) %>% arrange(period,threshold,desc(Freq_homopl))
	
	return(df_homopl_each_threshold)
}

most_common_homopl_each_thr_syn <- get_most_common_homopls_each_threshold(df_syn_clustered_homopl_tables)
most_common_homopl_each_thr_non_syn <- get_most_common_homopls_each_threshold(df_non_syn_clustered_homopl_tables)

# TOP50 for each period
get_most_common_homopls_each_period <- function(df) {
	df_homopl_binded <- rbindlist(df)
	df_homopl_binded <- df_homopl_binded %>% separate(period_threshold, c("period","threshold"), sep="_")
	df_homopl_binded$period <- as.integer(df_homopl_binded$period)
	df_homopl_binded$threshold <- as.numeric(df_homopl_binded$threshold)
	df_homopl_binded$Freq_homopl <- as.integer(df_homopl_binded$Freq_homopl)
	df_homopl_each_period <- df_homopl_binded %>% group_by(period) %>% slice_max(n=50, order_by=Freq_homopl) %>% arrange(period,threshold,desc(Freq_homopl))
	
	return(df_homopl_each_period)
}

most_common_homopl_each_period_syn <- get_most_common_homopls_each_period(df_syn_clustered_homopl_tables)
most_common_homopl_each_period_non_syn <- get_most_common_homopls_each_period(df_non_syn_clustered_homopl_tables)

get_most_common_homopls_abs_freq <- function(df) {
	df_homopl_binded <- rbindlist(df)
	df_homopl_binded <- df_homopl_binded %>% separate(period_threshold, c("period","threshold"), sep="_")
	df_homopl_binded$period <- as.integer(df_homopl_binded$period)
	df_homopl_binded$threshold <- as.numeric(df_homopl_binded$threshold)
	df_homopl_binded$Freq_homopl <- as.integer(df_homopl_binded$Freq_homopl)
	most_common_homopl_abs_freq <- df_homopl_binded %>% slice_max(n=50, order_by=Freq_homopl) %>% arrange(desc(Freq_homopl),period,threshold)
	return(most_common_homopl_abs_freq)
}

# TOP50 by absolute frequency
most_common_homopl_abs_freq_syn <- get_most_common_homopls_abs_freq(df_syn_clustered_homopl_tables)
most_common_homopl_abs_freq_non_syn <- get_most_common_homopls_abs_freq(df_non_syn_clustered_homopl_tables)

# TOP50 only considering S regions
get_most_common_homopls_s_regions <- function(df) {
	df_homopl_binded <- rbindlist(df)
	df_homopl_binded <- df_homopl_binded %>% separate(period_threshold, c("period","threshold"), sep="_")
	df_homopl_binded$period <- as.integer(df_homopl_binded$period)
	df_homopl_binded$threshold <- as.numeric(df_homopl_binded$threshold)
	df_homopl_binded$Freq_homopl <- as.integer(df_homopl_binded$Freq_homopl)
	most_common_homopl_s_regions <- df_homopl_binded %>% filter(protein=="S") %>% slice_max(n=50, order_by=Freq_homopl) %>% arrange(period,threshold,desc(Freq_homopl))
	return(most_common_homopl_s_regions)
}

most_common_homopl_spike_non_syn <- get_most_common_homopls_s_regions(df_non_syn_clustered_homopl_tables)

# TOP-10 for suitable combinations of period and threshold (IMPORTANT: using this now)
get_most_common_homopls <- function(df_list, out_suffix) {
	df_homopl_binded <- rbindlist(df_list)
	df_homopl_binded <- df_homopl_binded %>% separate(period_threshold, c("period","threshold"), sep="_")
	df_homopl_binded$id <- paste0(df_homopl_binded$period,"_",df_homopl_binded$threshold)
	df_splitted <- split(df_homopl_binded, df_homopl_binded$id)
	df_splitted_tops <- list()
	system(glue("mkdir -p rds/most_common_homopl_{out_suffix}/"))
	# Each combination will have its own file
	for(i in 1:length(df_splitted)) {
		df_splitted_tops[[i]] <- df_splitted[[i]] %>% slice_max(n=10, order_by=Freq_homopl) %>% arrange(period,threshold,desc(Freq_homopl))
		write.csv(df_splitted_tops[[i]], file=glue("rds/most_common_homopl_{out_suffix}/{unique(df_splitted_tops[[i]]$id)}.csv"), quote=F, row.names=F)
	}
	
	# File where all dfs above are concatenated
	df_splitted_tops_all <- rbindlist(df_splitted_tops)
	write.csv(df_splitted_tops_all, file=glue("rds/most_common_homopl_{out_suffix}/all_most_common_concatenated.csv"), quote=F, row.names=F)
	
	# Including all cols
	df_splitted_tops_all_cols <- df_splitted_tops_all %>% distinct(defining_mut, .keep_all = TRUE)
	
	# if `defining_mut` duplicated for diff periods or thresholds, append it to appropriate column (collapsed using |) and show only 1 row with each mut
	df_splitted_tops_all_no_rep <- df_splitted_tops_all %>% group_by(defining_mut) %>% summarize(period=paste0(na.omit(period), collapse="|"), threshold=paste0(na.omit(threshold), collapse="|"), Freq_homopl=paste0(na.omit(Freq_homopl), collapse="|")) %>% ungroup()
	df_splitted_tops_all_no_rep <- df_splitted_tops_all_no_rep %>% select(period, threshold, defining_mut, Freq_homopl)
	df_splitted_tops_all_no_rep <- df_splitted_tops_all_no_rep %>% inner_join(df_splitted_tops_all_cols, by="defining_mut")
	df_splitted_tops_all_no_rep <- df_splitted_tops_all_no_rep[,-c(5:7)]
	write.csv(df_splitted_tops_all_no_rep, file=glue("rds/most_common_homopl_{out_suffix}/all_most_common_concatenated_no_rep.csv"), quote=F, row.names=F)
	
	return(list(df_splitted_tops, df_splitted_tops_all, df_splitted_tops_all_no_rep))
}

most_common_homopl_syn <- get_most_common_homopls(df_syn_clustered_homopl_tables,"syn")
most_common_homopl_non_syn <- get_most_common_homopls(df_non_syn_clustered_homopl_tables, "non_syn")

# prepare dfs for bubble plot selection analysis from HypHy identified mutations
pos_sel_sites2 <- pos_sel_sites %>% select(protein, site, prot_site)
pos_sel_sites2$analysis_identified <- "HyPhy"
pos_sel_sites2$Freq_homopl <- NA
pos_sel_sites2$period <- NA
pos_sel_sites2$threshold <- NA
pos_sel_sites2 <- pos_sel_sites2 %>% select(period, threshold, protein, site, prot_site, analysis_identified, Freq_homopl)

most_common_homopl_non_syn_all <- most_common_homopl_non_syn[[2]]
most_common_homopl_non_syn_all$analysis_identified <- "TreeBasedClustering"
most_common_homopl_non_syn_all_adj <- most_common_homopl_non_syn_all %>% select(period, threshold, protein, site, prot_site, analysis_identified, Freq_homopl)

joined_most_common_indep <- rbind(pos_sel_sites2, most_common_homopl_non_syn_all_adj)

joined_most_common_indep1 <- joined_most_common_indep %>% group_by(prot_site) %>% summarize(analysis_identified=paste0(na.omit(analysis_identified), collapse=";")) %>% ungroup()
# Make sure there is only one TreeBasedClustering for each row (remove if more than one semicolon)
joined_most_common_indep1$analysis_identified <- sub(sprintf("^((.*?;){%d}.*?);.*", 1), "\\1", joined_most_common_indep1$analysis_identified)
#unique(joined_most_common_indep1$analysis_identified)
joined_most_common_indep1$analysis_identified <- ifelse(joined_most_common_indep1$analysis_identified=="TreeBasedClustering;TreeBasedClustering", "TreeBasedClustering", joined_most_common_indep1$analysis_identified)
#unique(joined_most_common_indep1$analysis_identified)
joined_most_common_indep_ok <- joined_most_common_indep %>% inner_join(joined_most_common_indep1, by="prot_site")

# Plot homoplasies identified (3 categories: (i) HyPhy, (ii) HyPhy;TreeBasedClustering, (iii) TreeBasedClustering)
ggplot(joined_most_common_indep_ok, aes(y=prot_site, x=analysis_identified.y, color=analysis_identified.y)) + #color=analysis_identified.y
	geom_point(alpha=0.8) + labs(y="Protein site", x="Analyses identified") + theme_minimal() +
	scale_color_manual(values=c("#7c9885","#e1ad01", "#033f63")) + theme(legend.position="none", axis.text.y=element_text(size=7.5))
ggsave(glue("rds/comparison_bubble/homopl_method_combinations.png"), width=12, height=10, dpi=600, bg="white")

joined_most_common_indep_ok_high_freq_period <- joined_most_common_indep_ok %>% group_by(period, prot_site) %>% slice_max(Freq_homopl)

#sort(unique(joined_most_common_indep_ok_high_freq_period$Freq_homopl))

system("mkdir -p rds/comparison_bubble/")
# Plot max frequency of homoplasy (bubble sizes) for each one of the 5 periods (coloured by analysis where identified)
ggplot(joined_most_common_indep_ok_high_freq_period, aes(x=period, y=prot_site, size=Freq_homopl, color=analysis_identified.y)) +
	geom_point(alpha=0.5) + labs(x="Period of pandemic", y="Protein site") + theme_minimal() +
	scale_size_area(name="Highest freq of homopl", breaks=c(20,50,150,500,1300)) + scale_color_manual(values=c("#7c9885","#e1ad01"), name="Analysis identified")
ggsave(glue("rds/comparison_bubble/max_freq_all_periods.png"), width=10, height=8, dpi=600, bg="white")

joined_most_common_indep_ok_each_period <- list()
breaks_period_thr <- list(c(5,10), c(5,10,15,30), c(10,25,50,100,230), c(10,20,40,80,150,250,400), c(10,25,50,100,200,300,500,1000,1300))
#scale_func <- function(x) sprintf("%.0f", x)

# For each one of the 5 periods, plot the thresholds met and freq of homoplasy in each case (coloured by analysis where identified)
for(i in 1:length(table_names_periods)) {
	print(paste("Index=",i))
	joined_most_common_indep_ok_each_period[[i]] <- joined_most_common_indep_ok %>% filter(period==i)
	joined_most_common_indep_ok_each_period[[i]]$threshold <- factor(joined_most_common_indep_ok_each_period[[i]]$threshold, levels=table_names_thresholds) 
	
	print(sort(unique(joined_most_common_indep_ok_each_period[[i]]$Freq_homopl)))
	
	p <- ggplot(joined_most_common_indep_ok_each_period[[i]], aes(x=as.factor(threshold), y=prot_site, size=Freq_homopl, color=analysis_identified.y)) +
		geom_point(alpha=0.5) + labs(x="Quantile threshold", y="Protein site", title=glue("Period {i}")) + theme_minimal() +
		scale_size_area(name="Homoplasy frequency", breaks=breaks_period_thr[[i]]) + scale_color_manual(values=c("#7c9885","#e1ad01"), name="Analysis identified") #+
		#scale_x_discrete(breaks=table_names_thresholds) #labels=scale_funcs
	ggsave(glue("rds/comparison_bubble/all_thr_period{i}.png"), width=10, height=8, dpi=600, bg="white")
}

