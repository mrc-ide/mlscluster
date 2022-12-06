library(glue)
library(stringr)
library(dplyr)
library(data.table)

#load(file="rds/env.RData")
seq_muts_df_counts <- readRDS("rds/polymorphic_site_counts.rds")

pos_sel_sites <- read.csv("config/positive_selected_sites_20221006.tsv", sep="\t", header=TRUE)
pos_sel_sites$prot_site <- paste0(pos_sel_sites$protein,":",pos_sel_sites$site)

c_path <- "results"
path_periods <- c(glue("{c_path}/01_sc2_root_to_dec2020"), glue("{c_path}/02_sc2_jan2021_to_may2021"),glue("{c_path}/03_sc2_jun2021_to_dec2021"), 
																		glue("{c_path}/04_sc2_jan2022_to_apr2022"), glue("{c_path}/05_sc2_whole_period"))

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
		
		# Split into syn and non-syn
		system("mkdir -p rds/stat_ready_csvs_all_sites/")
		system("mkdir -p rds/stat_ready_csvs_only_clustered/")
		joined_mut_sites_clustered_syn[[i,j]] <- joined_mut_sites_clustered[[i,j]][joined_mut_sites_clustered[[i,j]]$syn == 1,]
		csv_cp[[i,j]] <- joined_mut_sites_clustered_syn[[i,j]][ order(joined_mut_sites_clustered_syn[[i,j]]$Freq_homopl, decreasing=T) ,]
		write.csv(csv_cp[[i,j]], glue("rds/stat_ready_csvs_all_sites/df_all_sites_P{i}_{thr_path[j]}_SYN.csv"), quote=F, row.names=F)
		csv_cp[[i,j]] <- csv_cp[[i,j]][(as.integer(csv_cp[[i,j]]$Freq_homopl) > 0) & (as.integer(csv_cp[[i,j]]$is_clustered)==1), ]
	 write.csv(csv_cp[[i,j]], glue("rds/stat_ready_csvs_only_clustered/df_only_clustered_P{i}_{thr_path[j]}_SYN.csv"), quote=F, row.names=F)
		joined_mut_sites_clustered_non_syn[[i,j]] <- joined_mut_sites_clustered[[i,j]][joined_mut_sites_clustered[[i,j]]$syn == 0,]
		csv_cp2[[i,j]] <- joined_mut_sites_clustered_non_syn[[i,j]][ order(joined_mut_sites_clustered_non_syn[[i,j]]$Freq_homopl, decreasing=T) ,]
		write.csv(csv_cp2[[i,j]], glue("rds/stat_ready_csvs_all_sites/df_all_sites_P{i}_{thr_path[j]}_NON_SYN.csv"), quote=F, row.names=F)
		csv_cp2[[i,j]] <- csv_cp2[[i,j]][(as.integer(csv_cp2[[i,j]]$Freq_homopl) > 0) & (as.integer(csv_cp2[[i,j]]$is_clustered)==1), ]
		write.csv(csv_cp2[[i,j]], glue("rds/stat_ready_csvs_only_clustered/df_only_clustered_P{i}_{thr_path[j]}_NON_SYN.csv"), quote=F, row.names=F)
		
		print(glue("{path_periods[i]}/{thr_path[j]}/"))

		# genomic index
		poisson_models_syn[[i,j]] <- glm(Freq_homopl ~ is_clustered + spec_regions + aa_length + indep_found_pos_selection, data=joined_mut_sites_clustered_syn[[i,j]], na.action=na.omit, family=poisson(link = "log")) #defining_mut + genomic_index +s_mut_region_interest +syn
		print("Summary poisson (SYN): ")
		print(summary(poisson_models_syn[[i,j]]))

		poisson_models_non_syn[[i,j]] <- glm(Freq_homopl ~ is_clustered + spec_regions + aa_length + indep_found_pos_selection, data=joined_mut_sites_clustered_non_syn[[i,j]], na.action=na.omit, family=poisson(link = "log")) #defining_mut + genomic_index +s_mut_region_interest +syn
		print("Summary poisson (NON-SYN): ")
		print(summary(poisson_models_non_syn[[i,j]]))

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
df_syn_clustered_homopl <- change_ids(df_syn_clustered_homopl)
df_syn_clustered_homopl <- df_syn_clustered_homopl[df_syn_clustered_homopl$Freq_homopl >= 5,]
df_syn_clustered_homopl <- df_syn_clustered_homopl[order(df_syn_clustered_homopl$defining_mut, decreasing=T),]
hist(df_syn_clustered_homopl$Freq_homopl)
write.csv(df_syn_clustered_homopl, glue("rds/only_clustered_ordered_homopl/syn_clustered_homopl_ordered.csv"), quote=F, row.names=F)
df_non_syn_clustered_homopl <- change_ids(df_non_syn_clustered_homopl)
df_non_syn_clustered_homopl <- df_non_syn_clustered_homopl[df_non_syn_clustered_homopl$Freq_homopl >= 5,]
df_non_syn_clustered_homopl <- df_non_syn_clustered_homopl[order(df_non_syn_clustered_homopl$defining_mut, decreasing=T),]
hist(df_non_syn_clustered_homopl$Freq_homopl)
write.csv(df_non_syn_clustered_homopl, glue("rds/only_clustered_ordered_homopl/non_syn_clustered_homopl_ordered.csv"), quote=F, row.names=F)

# Counts number of rows for each index
df_syn_clustered_homopl_counts <- df_syn_clustered_homopl %>% count(period_threshold)
df_non_syn_clustered_homopl_counts <- df_non_syn_clustered_homopl %>% count(period_threshold)

plot_distr_freqs_spec_regions <- function(df,aa_change) {
	system(glue("mkdir -p rds/plot_mut_mean_freqs_{aa_change}"))
	# Split for each period+threshold
	df_splitted <- split(df, df$period_threshold)
	# Remove period+threshold with less than 10 rows
	for(i in 1:length(df_splitted)) {
		if(nrow(df_splitted[[i]]) < 10) {
			df_splitted[i] <- 0
		}
	}
	df_splitted <- Filter(is.data.frame, df_splitted)
	
	system(glue("mkdir -p rds/plot_mut_counts_{aa_change}"))
	df_splitted_mean <- df_splitted
	# Mean of freqs per genomic region
	for(j in 1:length(df_splitted)) {
		df_splitted_mean[[j]] <- df_splitted_mean[[j]] %>% group_by(spec_regions) %>% summarise(mean_freq_homopl=mean(Freq_homopl))
		if(aa_change=="syn") {
			df_splitted_mean[[j]]$spec_regions <- sub('.*:', "", df_splitted_mean[[j]]$spec_regions)
		}
		ggplot(data=df_splitted_mean[[j]], aes(x=spec_regions, y=mean_freq_homopl, fill=spec_regions)) + geom_bar(stat="identity") + labs(x="Genomic region",y="Mean homoplasy Freq") + theme_classic() + theme(legend.position="none", axis.text.x=element_text(size=5))
		ggsave(glue("rds/plot_mut_mean_freqs_{aa_change}/{names(df_splitted_mean)[j]}_mean_freq.png"), width=8, height=6, dpi=600, bg="white")
	}
	
	df_splitted_count_muts <- df_splitted
	for(k in 1:length(df_splitted_count_muts)) {
		df_splitted_count_muts[[k]] <- df_splitted_count_muts[[k]] %>% group_by(spec_regions) %>% summarise(count_muts=n())
		if(aa_change=="syn") {
			df_splitted_count_muts[[k]]$spec_regions <- sub('.*:', "", df_splitted_count_muts[[k]]$spec_regions)
		}
		ggplot(data=df_splitted_count_muts[[k]], aes(x=spec_regions, y=count_muts, fill=spec_regions)) + geom_bar(stat="identity") + labs(x="Genomic region",y="# mutations") + theme_classic() + theme(legend.position="none", axis.text.x=element_text(size=5))
		ggsave(glue("rds/plot_mut_counts_{aa_change}/{names(df_splitted_count_muts)[k]}_mut_counts.png"), width=8, height=6, dpi=600, bg="white")
	}
	return(df_splitted_count_muts)
}

df_syn_clustered_homopl_plots <- plot_distr_freqs_spec_regions(df_syn_clustered_homopl,"syn")
df_non_syn_clustered_homopl_plots <- plot_distr_freqs_spec_regions(df_non_syn_clustered_homopl,"non-syn")
