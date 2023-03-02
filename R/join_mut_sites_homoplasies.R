library(glue)
library(stringr)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggrepel)
library(viridis)
#library(emmeans) #install.packages("emmeans")
# library(ggpubr)
# library(rstatix)
# library(car)

#load(file="rds/env.RData")

# Palette for thresholds
pal_thresholds <- c("#1984c5", "#22a7f0", "#63bff0", "#a7d5ed", "#f0f8ff", "#ffcccb", "#e1a692", "#de6e56", "#e14b31", "#c23728") # removed e2e2e2

# Palette for major_lineages
pal_lineages <- c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#fdcce5") #removed "#beb9db", "#8bd3c7"

syn_p <- "SYNSNP:"
lvls_syn_proteins <- c("NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NPS12","NSP13","NSP14","NSP15",
																								"NSP16","S","ORF3A","E","M","ORF6","ORF7A","ORF7B","ORF8","N","ORF10")

lvls_nonsyn_proteins_annots <- c("NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NSP12","NSP13","NSP14","NSP15","NSP16",
																																	"S","S:NTD","S:RBD","S:FCS","ORF3A","E","M","ORF6","ORF7A","ORF7B","ORF8","N","N:Linker_domain","ORF10")

lvls_nonsyn_proteins <- c("NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NSP12","NSP13","NSP14","NSP15","NSP16",
																																	"S","ORF3A","E","M","ORF6","ORF7A","ORF7B","ORF8","N","ORF10")

# lvls_major_lineages <- 

# Load config csv's
seq_muts_df_counts <- readRDS("rds/polymorphic_site_counts.rds")
pos_sel_sites <- read.csv("config/positive_selected_sites_20221006.tsv", sep="\t", header=TRUE)
pos_sel_sites$prot_site <- paste0(pos_sel_sites$protein,":",pos_sel_sites$site)
annot_gen_range_positions <- read.csv("config/aa_ranges_stat_model.tsv",sep="\t")
annot_gen_length_positions <- read.csv("config/aa_length_stat_model.tsv",sep="\t")

#seq_muts_df_counts_adj <- seq_muts_df_counts %>% select(mut_and_type, protein, aa_length, syn_non_syn)

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


plot_distr_freqs_spec_regions <- function(df,aa_change,remove_if_x_rows_period_thr=2, all_clust=FALSE) {
	system(glue("mkdir -p stat_results/{out_folder}/plot_mut_mean_freqs_{aa_change}_regions"))
	if(all_clust) {
		df <- df
	} else {
		# if false will get only clustered == 1
		df <- df[df$is_clustered==1,]
	}
	# Split for each period+threshold
	df_splitted <- split(df, df$period_threshold)
	# Remove period+threshold with less than 2 rows
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
	
	system(glue("mkdir -p stat_results/{out_folder}/plot_mut_counts_{aa_change}_regions"))
	df_splitted_mean <- df_splitted_mean_norm_size <- df_splitted
	# Mean of freqs per genomic region
	for(j in 1:length(df_splitted)) {
		df_splitted_mean[[j]] <- df_splitted_mean[[j]] %>% group_by(spec_regions) %>% summarise(mean_freq_homopl=mean(Freq_homopl))
		df_splitted_mean_norm_size[[j]] <- df_splitted_mean_norm_size[[j]] %>% group_by(spec_regions) %>% summarise(mean_freq_homopl_norm=mean(Freq_homopl_norm_size))
		if(aa_change=="syn") {
			df_splitted_mean[[j]]$spec_regions <- sub('.*:', "", df_splitted_mean[[j]]$spec_regions)
			df_splitted_mean_norm_size[[j]]$spec_regions <- sub('.*:', "", df_splitted_mean_norm_size[[j]]$spec_regions)
		}
		ggplot(data=df_splitted_mean[[j]], aes(x=spec_regions, y=mean_freq_homopl, fill=spec_regions)) + geom_bar(stat="identity") + labs(x="Genomic region",y="Mean homoplasy Freq") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
		ggsave(glue("stat_results/{out_folder}/plot_mut_mean_freqs_{aa_change}_regions/{names(df_splitted_mean)[j]}_mean_freq.png"), width=6, height=5, dpi=300, bg="white")
		ggplot(data=df_splitted_mean_norm_size[[j]], aes(x=spec_regions, y=mean_freq_homopl_norm, fill=spec_regions)) + geom_bar(stat="identity") + labs(x="Genomic region",y="Mean homoplasy Freq normalised by genomic region") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
		ggsave(glue("stat_results/{out_folder}/plot_mut_mean_freqs_{aa_change}_regions/{names(df_splitted_mean)[j]}_mean_freq_norm.png"), width=6, height=5, dpi=300, bg="white")
	}
	#View(df_splitted_mean[[1]])
	
	# Count of homoplasies per genomic region
	df_splitted_count_muts <- df_splitted #df_splitted_count_muts_norm_size
	for(k in 1:length(df_splitted_count_muts)) {
		df_splitted_count_muts[[k]] <- df_splitted_count_muts[[k]] %>% group_by(spec_regions) %>% summarise(count_muts=n())
		if(aa_change=="syn") {
			df_splitted_count_muts[[k]]$spec_regions <- sub('.*:', "", df_splitted_count_muts[[k]]$spec_regions)
		}
		ggplot(data=df_splitted_count_muts[[k]], aes(x=spec_regions, y=count_muts, fill=spec_regions)) + geom_bar(stat="identity") + labs(x="Genomic region",y="# mutations") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
		ggsave(glue("stat_results/{out_folder}/plot_mut_counts_{aa_change}_regions/{names(df_splitted_count_muts)[k]}_mut_counts.png"), width=6, height=5, dpi=300, bg="white")
	}
	
	return(df_splitted)
}

plot_distr_freqs_major_lineage <- function(df,aa_change,remove_if_x_rows_period_thr=2,all_clust=FALSE) {
	system(glue("mkdir -p stat_results/{out_folder}/plot_mut_mean_freqs_{aa_change}_lineage"))
	if(all_clust) {
		df <- df
	} else {
		# if false will get only clustered == 1
		df <- df[df$is_clustered==1,]
	}
	# Split for each period+threshold
	df_splitted <- split(df, df$period_threshold)
	# Remove period+threshold with less than 2 rows
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
	
	system(glue("mkdir -p stat_results/{out_folder}/plot_mut_counts_{aa_change}_lineage"))
	df_splitted_mean <- df_splitted_mean_norm_size <- df_splitted
	# Mean of freqs per major_lineage
	for(j in 1:length(df_splitted)) {
		df_splitted_mean[[j]] <- df_splitted_mean[[j]] %>% group_by(major_lineage) %>% summarise(mean_freq_homopl=mean(Freq_homopl))
		df_splitted_mean_norm_size[[j]] <- df_splitted_mean_norm_size[[j]] %>% group_by(major_lineage) %>% summarise(mean_freq_homopl_norm=mean(Freq_homopl_norm_size))
		ggplot(data=df_splitted_mean[[j]], aes(x=major_lineage, y=mean_freq_homopl, fill=major_lineage)) + geom_bar(stat="identity") + labs(x="Major lineage",y="Mean homoplasy Freq") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
		ggsave(glue("stat_results/{out_folder}/plot_mut_mean_freqs_{aa_change}_lineage/{names(df_splitted_mean)[j]}_mean_freq.png"), width=6, height=5, dpi=300, bg="white")
		ggplot(data=df_splitted_mean_norm_size[[j]], aes(x=major_lineage, y=mean_freq_homopl_norm, fill=major_lineage)) + geom_bar(stat="identity") + labs(x="Major lineage",y="Mean homoplasy Freq normalised by genomic region") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
		ggsave(glue("stat_results/{out_folder}/plot_mut_mean_freqs_{aa_change}_lineage/{names(df_splitted_mean)[j]}_mean_freq_norm.png"), width=6, height=5, dpi=300, bg="white")
	}
	#View(df_splitted_mean[[1]])
	
	# Count of homoplasies per major_lineage
	df_splitted_count_muts <- df_splitted # df_splitted_count_muts_norm_size
	for(k in 1:length(df_splitted_count_muts)) {
		df_splitted_count_muts[[k]] <- df_splitted_count_muts[[k]] %>% group_by(major_lineage) %>% summarise(count_muts=n()) #%>% mutate(count_norm=count_muts/aa_length)
		ggplot(data=df_splitted_count_muts[[k]], aes(x=major_lineage, y=count_muts, fill=major_lineage)) + geom_bar(stat="identity") + labs(x="Major lineage",y="# mutations") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
		ggsave(glue("stat_results/{out_folder}/plot_mut_counts_{aa_change}_lineage/{names(df_splitted_count_muts)[k]}_mut_counts.png"), width=6, height=5, dpi=300, bg="white")
	}
	
	#View(df_splitted_count_muts[[1]])
	return(df_splitted)
}

get_most_common_homopls_lineage_interest <- function(df, period_interest, lineage_interest) {
	df_homopl_binded <- rbindlist(df)
	
	#df_homopl_binded <- df_homopl_binded[df_homopl_binded$is_clustered==1,]
	
	# Split period_threshold into separated period and threshold columns
	df_homopl_binded <- df_homopl_binded %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
	# # Get only period of interest
	df_homopl <- df_homopl_binded %>% filter(period==period_interest)
	# Make sure that columns are numeric
	df_homopl$period <- as.integer(df_homopl$period)
	df_homopl$threshold <- as.numeric(df_homopl$threshold)
	df_homopl$Freq_homopl <- as.numeric(df_homopl$Freq_homopl)
	
	filt_high_thresholds <- filt_low_thresholds <- df_homopl
	# Get TOP-50 frequent homoplasies when threshold >= 10 (n=2)
	filt_high_thresholds <- filt_high_thresholds %>% filter(major_lineage==lineage_interest & threshold>=10) %>% slice_max(n=50, order_by=Freq_homopl) %>% arrange(period,threshold,desc(as.numeric(Freq_homopl)))
	# Get TOP-50 frequent homoplasies when threshold <= 5 (n=8)
	filt_low_thresholds <- filt_low_thresholds %>% filter(major_lineage==lineage_interest & threshold<=5) %>% slice_max(n=50, order_by=Freq_homopl) %>% arrange(period,threshold,desc(as.numeric(Freq_homopl)))
	
	return(list(filt_high_thresholds, filt_low_thresholds))
}

# TOP50 for each threshold
get_most_common_homopls_each_threshold <- function(df) {
	df_homopl_binded <- rbindlist(df)
	
	#df_homopl_binded <- df_homopl_binded[df_homopl_binded$is_clustered==1,]
	
	df_homopl_binded <- df_homopl_binded %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
	df_homopl_binded$period <- as.integer(df_homopl_binded$period)
	df_homopl_binded$threshold <- as.numeric(df_homopl_binded$threshold)
	df_homopl_binded$Freq_homopl <- as.integer(df_homopl_binded$Freq_homopl)
	df_homopl_each_threshold <- df_homopl_binded %>% group_by(threshold) %>% slice_max(n=50, order_by=Freq_homopl) %>% arrange(period,threshold,desc(as.numeric(Freq_homopl)))
	
	return(df_homopl_each_threshold)
}

get_most_common_homopls_abs_freq <- function(df) {
	df_homopl_binded <- rbindlist(df)
	
	#df_homopl_binded <- df_homopl_binded[df_homopl_binded$is_clustered==1,]
	
	df_homopl_binded <- df_homopl_binded %>% separate(period_threshold, c("period","threshold"), sep="_")
	df_homopl_binded$period <- as.integer(df_homopl_binded$period)
	df_homopl_binded$threshold <- as.numeric(df_homopl_binded$threshold)
	df_homopl_binded$Freq_homopl <- as.integer(df_homopl_binded$Freq_homopl)
	most_common_homopl_abs_freq <- df_homopl_binded[!duplicated(df_homopl_binded$defining_mut),] # Removing duplicates (period and threshold)
	most_common_homopl_abs_freq <- most_common_homopl_abs_freq %>% slice_max(n=50, order_by=Freq_homopl) %>% arrange(desc(as.numeric(Freq_homopl)),period,threshold)
	return(most_common_homopl_abs_freq)
}

# TOP50 only considering S regions
get_most_common_homopls_s_regions <- function(df) {
	df_homopl_binded <- rbindlist(df)
	
	#df_homopl_binded <- df_homopl_binded[df_homopl_binded$is_clustered==1,]
	
	df_homopl_binded <- df_homopl_binded %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
	df_homopl_binded$period <- as.integer(df_homopl_binded$period)
	df_homopl_binded$threshold <- as.numeric(df_homopl_binded$threshold)
	df_homopl_binded$Freq_homopl <- as.integer(df_homopl_binded$Freq_homopl)
	most_common_homopl_s_regions <- df_homopl_binded[!duplicated(df_homopl_binded$defining_mut),] # Removing duplicates (period and threshold)
	most_common_homopl_s_regions <- most_common_homopl_s_regions %>% filter(protein=="S") %>% slice_max(n=50, order_by=Freq_homopl) %>% arrange(period,threshold,desc(as.numeric(Freq_homopl)))
	return(most_common_homopl_s_regions)
}

# Ordered list of homoplasies for suitable combinations of period and threshold (IMPORTANT: using this now) - previously TOP-10
get_most_common_homopls <- function(df_list, out_suffix, is_clustered_flag=TRUE) {
	df_homopl_binded <- rbindlist(df_list)
	if(is_clustered_flag) {
		df_homopl_binded <- df_homopl_binded[df_homopl_binded$is_clustered==1,]
	}else {
		df_homopl_binded <- df_homopl_binded[df_homopl_binded$is_clustered==0,]
		#View(df_homopl_binded)
	}
	df_homopl_binded <- df_homopl_binded %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
	df_homopl_binded$id <- paste0(df_homopl_binded$period,"_",df_homopl_binded$threshold)
	df_splitted <- split(df_homopl_binded, df_homopl_binded$id)
	df_splitted_tops <- list()
	system(glue("mkdir -p stat_results/{out_folder}/most_common_homopl_{out_suffix}/"))
	# Each combination will have its own file
	for(i in 1:length(df_splitted)) {
		df_splitted_tops[[i]] <- df_splitted[[i]] %>% slice_max(n=50, order_by=Freq_homopl) %>% arrange(period,threshold,desc(as.numeric(as.numeric(Freq_homopl))))
		write.csv(df_splitted_tops[[i]], file=glue("stat_results/{out_folder}/most_common_homopl_{out_suffix}/{unique(df_splitted_tops[[i]]$id)}.csv"), quote=F, row.names=F)
	}

	# File where all dfs above are concatenated
	df_splitted_tops_all <- rbindlist(df_splitted_tops)
	write.csv(df_splitted_tops_all, file=glue("stat_results/{out_folder}/most_common_homopl_{out_suffix}/all_most_common_concatenated.csv"), quote=F, row.names=F)

	# Including all cols
	df_splitted_tops_all_cols <- df_splitted_tops_all %>% distinct(defining_mut, .keep_all = TRUE)

	# if `defining_mut` duplicated for diff periods or thresholds, append it to appropriate column (collapsed using |) and show only 1 row with each mut
	df_splitted_tops_all_no_rep <- df_splitted_tops_all %>% group_by(defining_mut) %>% summarize(period=paste0(na.omit(period), collapse="|"), threshold=paste0(na.omit(threshold), collapse="|"), Freq_homopl=paste0(na.omit(Freq_homopl), collapse="|")) %>% ungroup()
	df_splitted_tops_all_no_rep <- df_splitted_tops_all_no_rep %>% select(period, threshold, defining_mut, Freq_homopl)
	df_splitted_tops_all_no_rep <- df_splitted_tops_all_no_rep %>% inner_join(df_splitted_tops_all_cols, by="defining_mut")
	df_splitted_tops_all_no_rep <- df_splitted_tops_all_no_rep[,-c(5:7)]
	write.csv(df_splitted_tops_all_no_rep, file=glue("stat_results/{out_folder}/most_common_homopl_{out_suffix}/all_most_common_concatenated_no_rep.csv"), quote=F, row.names=F)

	return(list(df_splitted_tops, df_splitted_tops_all, df_splitted_tops_all_no_rep))
}

bubble_plots_replacements_lineages <- function(df, show_replacements=TRUE, out_suffix) {
	
	most_common_lineages_rep <- most_common_lineages_norep <- list()
	most_common_lineages_rep_first <- most_common_lineages_norep_first <- list()
	# For each one of the major_lineage, plot the thresholds met and freq of homoplasy in each case (coloured by analysis where identified) + only first threshold met
	for(i in 1:length(major_lineages)) {
		most_common_lineages_rep[[i]] <- df %>% filter(major_lineage==major_lineages[i])
		most_common_lineages_rep[[i]]$threshold <- factor(most_common_lineages_rep[[i]]$threshold, levels=table_names_thresholds)
		most_common_lineages_norep[[i]] <- most_common_lineages_rep[[i]]
		
		if(show_replacements) {
			ggplot(most_common_lineages_rep[[i]], aes(x=as.factor(threshold), y=defining_mut, size=Freq_homopl)) +
				geom_point(alpha=0.5) + labs(x="Quantile threshold", y="Homoplasy", title=glue("Lineage {major_lineages[i]}")) + theme_minimal() + theme(axis.text.y=element_text(size=5)) +
				scale_size_area(name="Homoplasy frequency") #+ scale_color_manual(values=c("#7c9885","#e1ad01"), name="Analysis identified")
			ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/repl_thr_lineage_{major_lineages[i]}.png"), width=12, height=15, dpi=600, bg="white")
			
			most_common_lineages_rep_first[[i]] <- most_common_lineages_rep[[i]]
			most_common_lineages_rep_first[[i]] <- most_common_lineages_rep_first[[i]] %>% distinct() %>% group_by(defining_mut) %>% slice_min(n=1, as.numeric(threshold))
			
			ggplot(most_common_lineages_rep_first[[i]], aes(x=as.factor(threshold), y=defining_mut, size=Freq_homopl)) +
				geom_point(alpha=0.5) + labs(x="Quantile threshold (first detected only)", y="Homoplasy", title=glue("Lineage {major_lineages[i]}")) + theme_minimal() + theme(axis.text.y=element_text(size=5)) +
				scale_size_area(name="Homoplasy frequency") #+ scale_color_manual(values=c("#7c9885","#e1ad01"), name="Analysis identified")
			ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/repl_first_detected_thr_lineage_{major_lineages[i]}.png"), width=12, height=15, dpi=600, bg="white")
		} else {
			
			ggplot(most_common_lineages_norep[[i]], aes(x=as.factor(threshold), y=prot_site, size=Freq_homopl, color=analysis_identified.y)) +
				geom_point(alpha=0.5) + labs(x="Quantile threshold", y="Protein site", title=glue("Lineage {major_lineages[i]}")) + theme_minimal() + theme(axis.text.y=element_text(size=5)) +
				scale_size_area(name="Homoplasy frequency") + scale_color_manual(values=c("#7c9885","#e1ad01"), name="Analysis identified")
			ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/onlysite_thr_lineage_{major_lineages[i]}.png"), width=12, height=15, dpi=600, bg="white")
			
			most_common_lineages_norep_first[[i]] <- most_common_lineages_norep[[i]]
			most_common_lineages_norep_first[[i]] <- most_common_lineages_norep_first[[i]] %>% distinct() %>% group_by(prot_site) %>% slice_min(n=1, as.numeric(threshold))
			
			ggplot(most_common_lineages_norep_first[[i]], aes(x=as.factor(threshold), y=prot_site, size=Freq_homopl, color=analysis_identified.y)) +
				geom_point(alpha=0.5) + labs(x="Quantile threshold (first detected only)", y="Protein site", title=glue("Lineage {major_lineages[i]}")) + theme_minimal() + theme(axis.text.y=element_text(size=5)) +
				scale_size_area(name="Homoplasy frequency") + scale_color_manual(values=c("#7c9885","#e1ad01"), name="Analysis identified")
			ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/onlysite_first_detected_thr_lineage_{major_lineages[i]}.png"), width=12, height=15, dpi=600, bg="white")
		}
	
	}
		
}

bubble_plots <- function(most_common_homopl_df, out_suffix) {
	# prepare dfs for bubble plot selection analysis from HypHy identified mutations
	pos_sel_sites2 <- pos_sel_sites %>% select(protein, site, prot_site)
	pos_sel_sites2$analysis_identified <- "HyPhy"
	pos_sel_sites2$Freq_homopl <- NA
	pos_sel_sites2$period <- NA
	pos_sel_sites2$threshold <- NA
	pos_sel_sites2$major_lineage <- NA
	pos_sel_sites2 <- pos_sel_sites2 %>% select(period, threshold, protein, site, prot_site, major_lineage, analysis_identified, Freq_homopl)
	
	most_common_homopl_all <- most_common_homopl_df[[2]]
	most_common_homopl_all$analysis_identified <- "TreeBasedClustering"
	most_common_homopl_all_adj <- most_common_homopl_all %>% select(period, threshold, protein, site, prot_site, major_lineage, analysis_identified, Freq_homopl)
	
	joined_most_common_indep <- rbind(pos_sel_sites2, most_common_homopl_all_adj)
	
	joined_most_common_indep1 <- joined_most_common_indep %>% group_by(prot_site) %>% summarize(analysis_identified=paste0(na.omit(analysis_identified), collapse=";")) %>% ungroup()
	# Make sure there is only one TreeBasedClustering for each row (remove if more than one semicolon)
	joined_most_common_indep1$analysis_identified <- sub(sprintf("^((.*?;){%d}.*?);.*", 1), "\\1", joined_most_common_indep1$analysis_identified)
	#unique(joined_most_common_indep1$analysis_identified)
	joined_most_common_indep1$analysis_identified <- ifelse(joined_most_common_indep1$analysis_identified=="TreeBasedClustering;TreeBasedClustering", "TreeBasedClustering", joined_most_common_indep1$analysis_identified)
	#unique(joined_most_common_indep1$analysis_identified)
	joined_most_common_indep_ok <- joined_most_common_indep %>% inner_join(joined_most_common_indep1, by="prot_site")
	
	system(glue("mkdir -p stat_results/{out_folder}/comparison_bubble_{out_suffix}/"))
	
	# Plot homoplasic sites identified (3 categories: (i) HyPhy, (ii) HyPhy;TreeBasedClustering, (iii) TreeBasedClustering)
	ggplot(joined_most_common_indep_ok, aes(y=prot_site, x=analysis_identified.y, color=analysis_identified.y)) + #color=analysis_identified.y
		geom_point(alpha=0.8) + labs(y="Protein site", x="Analyses identified") + theme_minimal() +
		scale_color_manual(values=c("#7c9885","#e1ad01", "#033f63")) + theme(legend.position="none", axis.text.y=element_text(size=5), aspect.ratio=16/9)
	ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/homopl_method_combinations.png"), width=12, height=15, dpi=600, bg="white")
	
	joined_most_common_indep_ok_high_freq_lin <- joined_most_common_indep_ok %>% group_by(major_lineage, prot_site) %>% slice_max(Freq_homopl)
	
	# print( sort(unique(joined_most_common_indep_ok_high_freq_lin$Freq_homopl)) )
	
	# Plot max frequency of homoplasic sites (bubble sizes) for each one of the major_lineage (coloured by analysis where identified)
	#View(joined_most_common_indep_ok_high_freq_lin)
	joined_most_common_indep_ok_high_freq_lin <- joined_most_common_indep_ok_high_freq_lin[joined_most_common_indep_ok_high_freq_lin$analysis_identified.y != "HyPhy",]
	joined_most_common_indep_ok_high_freq_lin <- joined_most_common_indep_ok_high_freq_lin[!is.na(joined_most_common_indep_ok_high_freq_lin$major_lineage),]
	
	ggplot(joined_most_common_indep_ok_high_freq_lin, aes(x=major_lineage, y=prot_site, size=Freq_homopl, color=analysis_identified.y)) +
		geom_point(alpha=0.5) + labs(x="Major lineage", y="Protein site") + theme_minimal() + theme(axis.text.y = element_text(size=5), aspect.ratio=16/9) +
		scale_size_area(name="Highest freq of homopl", breaks=c(20,50,150,500,1300)) + scale_color_manual(values=c("#7c9885","#e1ad01"), name="Analysis identified") #"#033f63"
	ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_all_lineages.png"), width=12, height=15, dpi=600, bg="white")
	
	most_common_homopl_all_high_freq_lin <- most_common_homopl_all %>% group_by(major_lineage, defining_mut) %>% slice_max(Freq_homopl)
	most_common_homopl_all_high_freq_lin <- most_common_homopl_all_high_freq_lin[!is.na(most_common_homopl_all_high_freq_lin$major_lineage),]
	# Plot max frequency of each homoplasy across thresholds (including WT and replaced aa) identified (without comparison with HyPhy)
	ggplot(most_common_homopl_all_high_freq_lin, aes(x=major_lineage, y=defining_mut, size=Freq_homopl)) + #color=analysis_identified.y
		geom_point(alpha=0.5) + labs(x="Major lineage", y="Homoplasy") + theme_minimal() + theme(axis.text.y = element_text(size=5), aspect.ratio=16/9) +
		scale_size_area(name="Highest freq of homopl") # breaks=c(20,50,150,500,1300) + scale_color_manual(values=c("#7c9885","#e1ad01", "#033f63")) +
	ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_replacement_lineages.png"), width=12, height=15, dpi=600, bg="white")
	
	# Bubble plots per major lineage comparing with other results (HyPhy): just showing protein site and not actual replacements
	bubble_plots_replacements_lineages(joined_most_common_indep_ok, show_replacements=FALSE, out_suffix)
	# Bubble plots per major lineage  showing actual replacements and no comparison with other analyses
	bubble_plots_replacements_lineages(most_common_homopl_all, show_replacements=TRUE, out_suffix)

}

genomewise_plot <- function(df, mut_type) {
	# TODO? Annotate ORF1AB NSPs, S regions and N linker
	df <- df %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
	df$threshold <- factor(df$threshold, levels=table_names_thresholds)
	
	if(mut_type == "syn") {
		genomic_ranges_df <- read.csv("config/genomic_ranges_plot_syn.tsv", header=T, sep="\t")
		# No split for SYN
		df_freq_genomic_regions_syn <- df %>% inner_join(genomic_ranges_df, by=c("spec_regions"="new_prot_name"))
		df_freq_genomic_regions_syn_first <- df_freq_genomic_regions_syn %>% distinct() %>% group_by(defining_mut) %>% slice_min(n=1, as.numeric(threshold))
		#View(df_freq_genomic_regions_syn_first)
		
		system(glue("mkdir -p stat_results/{out_folder}/genomewise_plot_{mut_type}/"))
		
		ggplot(df_freq_genomic_regions_syn_first, aes(x=genomic_index, y=Freq_homopl, color=threshold)) +
			geom_point(alpha=0.9) + labs(x="Nucleotide position", y="Frequency in TFP clusters") + theme_minimal() + #theme(axis.text.y=element_text(size=5)) +
			scale_color_manual(values=pal_thresholds, name="Threshold where first detected") +
			coord_cartesian(xlim = c(1, 29903), ylim = c(2,10)) + scale_y_continuous(breaks=seq(2,10,by=1))
		ggsave(glue("stat_results/{out_folder}/genomewise_plot_{mut_type}/syn_freqs_1st_threshold.png"), width=15, height=5, dpi=600, bg="white")
		
		ggplot(df_freq_genomic_regions_syn_first, aes(x=genomic_index, y=Freq_homopl, color=major_lineage)) +
			geom_point(alpha=0.9) + labs(x="Nucleotide position", y="Frequency in TFP clusters") + theme_minimal() + #theme(axis.text.y=element_text(size=5)) +
			scale_color_manual(values=pal_lineages, name="Major lineage") +
			coord_cartesian(xlim = c(0, 29903), ylim = c(2,10)) + scale_y_continuous(breaks=seq(2,10,by=1))
		ggsave(glue("stat_results/{out_folder}/genomewise_plot_{mut_type}/syn_freqs_major_lineage.png"), width=15, height=5, dpi=600, bg="white")
		
	}else if(mut_type == "non-syn") {
		genomic_ranges_df <- read.csv("config/genomic_ranges_plot_non_syn.tsv", header=T, sep="\t")
		#print(unique(df_freq_genomic_regions$protein))
		df_freq_genomic_regions_nonsyn <- df %>% inner_join(genomic_ranges_df, by=c("spec_regions"="new_prot_name"))
		# Split because each protein (ORF1AB + others) has its own relative coordinates
		df_freq_genomic_regions_split_nonsyn <- split(df_freq_genomic_regions_nonsyn, df_freq_genomic_regions_nonsyn$protein)
		df_freq_genomic_regions_split_nonsyn_first <- list()
		for(i in 1:length(df_freq_genomic_regions_split_nonsyn)) {
			df_freq_genomic_regions_split_nonsyn_first[[i]] <- df_freq_genomic_regions_split_nonsyn[[i]] %>% distinct() %>% group_by(defining_mut) %>% slice_min(n=1, as.numeric(threshold))
			
			system(glue("mkdir -p stat_results/{out_folder}/genomewise_plot_{mut_type}/"))
			
			ggplot(df_freq_genomic_regions_split_nonsyn_first[[i]], aes(x=genomic_index, y=Freq_homopl, color=threshold)) +
				geom_point(alpha=0.9) + labs(x="Aminoacid position", y="Frequency in TFP clusters") + theme_minimal() +
				#facet_grid(~spec_regions, scales = "free", switch = "x") +
				scale_color_manual(values=pal_thresholds, name="Threshold where first detected") +
				geom_text_repel(aes(label = ifelse(Freq_homopl>5,as.character(defining_mut),'')), size=2, box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=15)+
				coord_cartesian(xlim = c(df_freq_genomic_regions_split_nonsyn_first[[i]]$start[1], df_freq_genomic_regions_split_nonsyn_first[[i]]$end[1])) #, ylim = c(2,10)) + scale_y_continuous(breaks=seq(2,10,by=1))
			ggsave(glue("stat_results/{out_folder}/genomewise_plot_{mut_type}/nonsyn_freqs_1st_threshold_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.png"), width=15, height=5, dpi=600, bg="white")
			
			ggplot(df_freq_genomic_regions_split_nonsyn_first[[i]], aes(x=genomic_index, y=Freq_homopl, color=major_lineage)) +
				geom_point(alpha=0.9) + labs(x="Aminoacid position", y="Frequency in TFP clusters") + theme_minimal() +
				#facet_grid(~spec_regions, scales = "free", switch = "x") +
				scale_color_manual(values=pal_lineages, name="Major lineage") +
				#geom_tile(aes(x=spec_regions, y=1, fill=)) +
				geom_text_repel(aes(label = ifelse(Freq_homopl>2,as.character(defining_mut),'')), size=2, box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=15)+
				coord_cartesian(xlim = c(df_freq_genomic_regions_split_nonsyn_first[[i]]$start[1], df_freq_genomic_regions_split_nonsyn_first[[i]]$end[1])) #, ylim = c(2,10)) + scale_y_continuous(breaks=seq(2,10,by=1))
			ggsave(glue("stat_results/{out_folder}/genomewise_plot_{mut_type}/nonsyn_freqs_major_lineage_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.png"), width=15, height=5, dpi=600, bg="white")
			
		}
		#View(df_freq_genomic_regions_split_nonsyn_first[[1]])
	} else {
		stop("Wrong choice of mut_type: options are `syn` and `non-syn`")
	}
	
}

violin_boxplot_freq_variations <- function(df, mut_type) {
	df <- df %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
	df$threshold <- factor(df$threshold, levels=table_names_thresholds)
	
	df_freq_genomic_regions_first <- df %>% distinct() %>% group_by(defining_mut) %>% slice_min(n=1, as.numeric(threshold))
	
	if(mut_type == "syn") {
		df_freq_genomic_regions_first$spec_regions <- sub("SYNSNP:","",df_freq_genomic_regions_first$spec_regions)
		#df_freq_genomic_regions_first[] <- lapply(df_freq_genomic_regions_first, function(x) str_remove(x, "SYNSNP:"))
		#View(df_freq_genomic_regions_first)
		df_freq_genomic_regions_first$spec_regions <- factor(df_freq_genomic_regions_first$spec_regions, levels=lvls_syn_proteins)
	} else {
		df_freq_genomic_regions_first$spec_regions <- factor(df_freq_genomic_regions_first$spec_regions, levels=lvls_nonsyn_proteins_annots)
	}
	
	df_freq_genomic_regions_first <- df_freq_genomic_regions_first[!is.na(df_freq_genomic_regions_first$spec_regions),]
	
	sample_size_regions <- df_freq_genomic_regions_first %>% group_by(spec_regions) %>% summarize(num=n())
	
	system(glue("mkdir -p stat_results/{out_folder}/violin_boxplot_freq_{mut_type}/"))
	
	#inner_join(sample_size_regions) %>% mutate(myaxis = paste0(spec_regions, "\n", "n=", num)) %>%
	df_freq_genomic_regions_first %>%
		ggplot( aes(x=spec_regions, y=Freq_homopl, fill=spec_regions)) + #myaxis
		geom_violin(width=1.4) + 
		#geom_boxplot(width=0.1, color="grey", alpha=0.8) + 
		coord_cartesian(ylim=c(2,10)) +
		scale_fill_viridis(discrete = TRUE) + theme_minimal() + theme(legend.position="none",axis.text.x = element_text(size=6)) +
		labs(x="Genomic region", y="Frequency in TFP clusters")
	ggsave(glue("stat_results/{out_folder}/violin_boxplot_freq_{mut_type}/variation_freq_regions.png"), width=15, height=10, dpi=600, bg="white")

	sample_size_lineages <- df_freq_genomic_regions_first %>% group_by(major_lineage) %>% summarize(num=n())
	
	#inner_join(sample_size_lineages) %>% mutate(myaxis = paste0(major_lineage, "\n", "n=", num)) %>%
	df_freq_genomic_regions_first %>%
		ggplot( aes(x=major_lineage, y=Freq_homopl, fill=major_lineage)) + #myaxis
		geom_violin(width=1.4) + 
		#geom_boxplot(width=0.1, color="grey", alpha=0.8) + 
		coord_cartesian(ylim=c(2,10)) +
		scale_fill_viridis(discrete = TRUE) + theme_minimal() + theme(legend.position="none",axis.text.x = element_text(size=10)) +
		labs(x="Genomic region", y="Frequency in TFP clusters")
	ggsave(glue("stat_results/{out_folder}/violin_boxplot_freq_{mut_type}/variation_freq_lineages.png"), width=15, height=10, dpi=600, bg="white")
	
}

stacked_bar_freqs <- function(df, mut_type) {
	df <- df %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
	df$threshold <- factor(df$threshold, levels=table_names_thresholds)
	
	df_freq_genomic_regions_first <- df %>% distinct() %>% group_by(defining_mut) %>% slice_min(n=1, as.numeric(threshold))
	df_freq_genomic_regions_first <- df_freq_genomic_regions_first[!is.na(df_freq_genomic_regions_first$spec_regions),]
	
	df_freq_genomic_regions_first$Freq_norm <- df_freq_genomic_regions_first$Freq_homopl / df_freq_genomic_regions_first$aa_length
	
	if(mut_type == "syn") {
		df_freq_genomic_regions_first$spec_regions <- sub(".*:","",df_freq_genomic_regions_first$spec_regions)
		#df_freq_genomic_regions_first[] <- lapply(df_freq_genomic_regions_first, function(x) str_remove(x, "SYNSNP:"))
		#View(df_freq_genomic_regions_first)
		df_freq_genomic_regions_first$spec_regions <- factor(df_freq_genomic_regions_first$spec_regions, levels=lvls_syn_proteins)
	} else {
		df_freq_genomic_regions_first$spec_regions <- factor(df_freq_genomic_regions_first$spec_regions, levels=lvls_nonsyn_proteins_annots)
	}
	
	df_freq_genomic_regions_first_counts <- df_freq_genomic_regions_first %>% group_by(spec_regions, major_lineage) %>% summarise(count_uniq_homopls=n())

	system(glue("mkdir -p stat_results/{out_folder}/stacked_bar_uniq_counts_{mut_type}/"))
	
	ggplot(df_freq_genomic_regions_first_counts, aes(fill=major_lineage, y=spec_regions, x=count_uniq_homopls)) + 
		geom_bar(position="stack", stat="identity") + scale_fill_manual(values=pal_lineages, name="Major_lineage") + theme(axis.text.y = element_text(size=10)) + theme_minimal() +
		labs(x="Count of unique homoplasies", y="Genomic region")
	ggsave(glue("stat_results/{out_folder}/stacked_bar_uniq_counts_{mut_type}/stacked_bar_uniq_counts.png"), width=8, height=6, dpi=600, bg="white")
	
	df_freq_genomic_regions_first_counts_norm <- df_freq_genomic_regions_first %>% group_by(spec_regions, major_lineage) %>% summarise(count_uniq_homopls_norm=sum(Freq_norm))
	
	ggplot(df_freq_genomic_regions_first_counts_norm, aes(fill=major_lineage, y=spec_regions, x=count_uniq_homopls_norm)) + 
		geom_bar(position="stack", stat="identity") + scale_fill_manual(values=pal_lineages, name="Major_lineage") + theme(axis.text.y = element_text(size=10)) + theme_minimal() +
		labs(x="Count of unique homoplasies (normalised by genomic size)", y="Genomic region")
	ggsave(glue("stat_results/{out_folder}/stacked_bar_uniq_counts_{mut_type}/stacked_bar_uniq_counts_norm.png"), width=8, height=6, dpi=600, bg="white")

}

# Function to facilitate comparison between periods including and excluding Omicron muts
clustering_model_fitting <- function(path_stats, out_folder) {
	
	if(dir.exists(path_stats))
		path <- path_stats
	else
		stop("Directory does not exist!")
	
	system(glue("mkdir -p stat_results/{out_folder}/"))
	
	clustered_dfs <- joined_mut_sites_clustered <- joined_mut_sites_clustered_syn <- joined_mut_sites_clustered_non_syn <- list()
	poisson_models_syn <- poisson_models_non_syn <- lr_models_syn <- lr_models_non_syn <- list()
	poisson_models_syn_ci <- poisson_models_non_syn_ci <- lr_models_syn_ci <- lr_models_non_syn_ci <- list()
	csv_cp <- csv_cp2 <- list()
	
	for(i in 1:length(path_thresholds)) {
		# print("====================")
		# print(path_thresholds[[i]])
		# print(glue("{path_stats}/{path_thresholds[i]}/clustered_all_df.csv"))
		# print("====================")
		# For each threshold read its csv into a df
		clustered_dfs[[i]] <- read.csv(glue("{path_stats}/{path_thresholds[i]}/clustered_all_df.csv"), header=T)
		
		# print(nrow(clustered_dfs[[i]][clustered_dfs[[i]]$is_clustered==1,]))
		# print(nrow(clustered_dfs[[i]][clustered_dfs[[i]]$is_clustered==0,]))
		
		clustered_dfs[[i]] <- clustered_dfs[[i]][clustered_dfs[[i]]$is_clustered==1,]
		# IMPORTANT: Make sure is_clustered==0 (not detected by stats) have zero Freq in tfp clusters (is_clustered==1)
		# clustered_dfs[[i]]$Freq_homopl <- ifelse(clustered_dfs[[i]]$is_clustered==0, yes=0, no=clustered_dfs[[i]]$Freq_homopl)
		
		clustered_dfs[[i]]$major_lineage <- ifelse(grepl("Beta_B.1.351", clustered_dfs[[i]]$major_lineage), "Other", as.character(clustered_dfs[[i]]$major_lineage))
		clustered_dfs[[i]]$major_lineage <- ifelse(grepl("Gamma_P.1", clustered_dfs[[i]]$major_lineage), "Other", as.character(clustered_dfs[[i]]$major_lineage))
		
		#print(nrow(clustered_dfs[[i]]))
		
		#View(clustered_dfs[[1]])
		
		# Join with seq_muts_df_counts from COG-UK (CAUTION: this adds all polymorphic sites, not only homoplasic ones)
		# joined_mut_sites_clustered[[i]] <- clustered_dfs[[i]] %>% full_join(seq_muts_df_counts, by=c("defining_mut"="mut_and_type")) #inner_join
		# print(nrow(joined_mut_sites_clustered[[i]]))
		joined_mut_sites_clustered[[i]] <- clustered_dfs[[i]]
		
		# Extract genome index position
		rgx_genomic_idx <- regexpr(':[A-Z]{1}[0-9]{1,5}[A-Z*]{1}', joined_mut_sites_clustered[[i]]$defining_mut)
		comm_coord <- rep(NA,length(joined_mut_sites_clustered[[i]]$defining_mut))
		comm_coord[rgx_genomic_idx != -1] <- str_sub( regmatches(joined_mut_sites_clustered[[i]]$defining_mut, rgx_genomic_idx), 3, -2)
		comm_coord <- as.integer(comm_coord)
		comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("genomic_index")
		joined_mut_sites_clustered[[i]] <- cbind(joined_mut_sites_clustered[[i]], comm_coord_df)
		
		# Adjust df columns
		joined_mut_sites_clustered[[i]] <- joined_mut_sites_clustered[[i]] %>% select(defining_mut, Freq_homopl, is_clustered, major_lineage, syn_non_syn, protein, aa_length, indep_found_pos_selection, genomic_index) #s_mut_region_interest
		colnames(joined_mut_sites_clustered[[i]]) <- c("defining_mut","Freq_homopl","is_clustered","major_lineage","syn","protein","aa_length","indep_found_pos_selection","genomic_index") #s_mut_region_interest
		
		joined_mut_sites_clustered[[i]]$protein <- sub("\\:.*", "", joined_mut_sites_clustered[[i]]$defining_mut)
		
		setDT(joined_mut_sites_clustered[[i]]); setDT(annot_gen_range_positions)
		joined_mut_sites_clustered[[i]] <- joined_mut_sites_clustered[[i]][annot_gen_range_positions, on=c("protein==name_cog_md","genomic_index>=start","genomic_index<=end"), spec_regions := new_prot_name]
		joined_mut_sites_clustered[[i]] <- joined_mut_sites_clustered[[i]][annot_gen_length_positions, on=c("spec_regions==new_prot_name"), aa_length := length]
		
		# syn flag
		joined_mut_sites_clustered[[i]]$syn <- ifelse(joined_mut_sites_clustered[[i]]$syn == as.character("non-syn"), 0, 1)
		#joined_mut_sites_clustered[[i]]$syn <- ifelse(joined_mut_sites_clustered[[i]]$protein == "SYNSNP", 1, 0)
		# independent found under selection flag
		joined_mut_sites_clustered[[i]]$indep_found_pos_selection <- ifelse(joined_mut_sites_clustered[[i]]$indep_found_pos_selection == as.character("yes"), 1, 0)
		# Add zeros to Freq_homopl and is_clustered if polymorphic site==NA (not found on homoplasy tables)
		joined_mut_sites_clustered[[i]] <- joined_mut_sites_clustered[[i]] %>% mutate(Freq_homopl = ifelse(is.na(Freq_homopl), 0, Freq_homopl))
		joined_mut_sites_clustered[[i]] <- joined_mut_sites_clustered[[i]] %>% mutate(is_clustered = ifelse(is.na(is_clustered), 0, is_clustered))
		
		# Not sure if should use this since flag is already set above
		rgx_scpss <- regexpr(':[A-Z]{1}[0-9]{1,5}[A-Z*]{1}', joined_mut_sites_clustered[[i]]$defining_mut)
		comm_coord <- rep(NA,length(joined_mut_sites_clustered[[i]]$defining_mut))
		# Mutation coordinate of homoplasies
		comm_coord[rgx_scpss != -1] <- str_sub( regmatches(joined_mut_sites_clustered[[i]]$defining_mut, rgx_scpss), 3, -2)
		comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("site")
		joined_mut_sites_clustered[[i]] <- cbind(joined_mut_sites_clustered[[i]], comm_coord_df)
		joined_mut_sites_clustered[[i]]$prot_site <- paste0(joined_mut_sites_clustered[[i]]$protein,":",joined_mut_sites_clustered[[i]]$site)
		#joined_mut_sites_clustered[[i]]$indep_found_pos_selection <- ifelse(joined_mut_sites_clustered[[i]]$prot_site %in% pos_sel_sites$prot_site, 1, 0)
		
		# For logistic regression, is_clustered_h3 will be used and it is true if Freq_homopl is > 3
		joined_mut_sites_clustered[[i]]$is_clustered_h3 <- ifelse( (as.integer(joined_mut_sites_clustered[[i]]$Freq_homopl) > 3) & (joined_mut_sites_clustered[[i]]$is_clustered==1), 1, 0)
		# Normalised counts (for inspection only)
		joined_mut_sites_clustered[[i]]$Freq_homopl_norm_size <- as.integer(joined_mut_sites_clustered[[i]]$Freq_homopl) / as.integer(joined_mut_sites_clustered[[i]]$aa_length)
		
		# Make sure to remove all rows with NAs in key columns
		# Freq_homopl, major_lineage, is_clustered, is_clustered_h3, spec_regions, aa_length, indep_found_pos_selection
		joined_mut_sites_clustered[[i]] <- joined_mut_sites_clustered[[i]][!is.na(joined_mut_sites_clustered[[i]]$Freq_homopl) | !is.na(joined_mut_sites_clustered[[i]]$major_lineage) | !is.na(joined_mut_sites_clustered[[i]]$is_clustered) | !is.na(joined_mut_sites_clustered[[i]]$is_clustered_h3) | !is.na(joined_mut_sites_clustered[[i]]$spec_regions) | !is.na(joined_mut_sites_clustered[[i]]$aa_length) | !is.na(joined_mut_sites_clustered[[i]]$indep_found_pos_selection), ]
		
		#joined_mut_sites_clustered[[i]]$major_lineage <- as.factor( joined_mut_sites_clustered[[i]]$major_lineage )
		joined_mut_sites_clustered[[i]]$site <- as.numeric(joined_mut_sites_clustered[[i]]$site)
		#joined_mut_sites_clustered[[i]]$spec_regions <- as.factor(joined_mut_sites_clustered[[i]]$spec_regions)
		#print(str(joined_mut_sites_clustered[[i]]))
		
		# Split into syn and non-syn
		system(glue("mkdir -p stat_results/{out_folder}/stat_ready_csvs_all_sites/"))
		system(glue("mkdir -p stat_results/{out_folder}/stat_ready_csvs_only_clustered/"))
		
		# SYN
		joined_mut_sites_clustered_syn[[i]] <- joined_mut_sites_clustered[[i]][joined_mut_sites_clustered[[i]]$syn == 1,]
		joined_mut_sites_clustered_syn[[i]]$major_lineage <- as.factor( joined_mut_sites_clustered_syn[[i]]$major_lineage )
		joined_mut_sites_clustered_syn[[i]]$spec_regions <- as.factor(joined_mut_sites_clustered_syn[[i]]$spec_regions)
		csv_cp[[i]] <- joined_mut_sites_clustered_syn[[i]][ order(joined_mut_sites_clustered_syn[[i]]$Freq_homopl, decreasing=T) ,]
		#csv_cp[[i]] <- joined_mut_sites_clustered_syn[[i]][ order(joined_mut_sites_clustered_syn[[i]]$Freq_homopl_norm_size, decreasing=T) ,]
		write.csv(csv_cp[[i]], glue("stat_results/{out_folder}/stat_ready_csvs_all_sites/df_all_sites_{path_thresholds[i]}_SYN.csv"), quote=F, row.names=F)
		csv_cp[[i]] <- csv_cp[[i]][(as.integer(csv_cp[[i]]$Freq_homopl) > 0) & (as.integer(csv_cp[[i]]$is_clustered)==1), ]
		#csv_cp[[i]] <- csv_cp[[i]][(as.integer(csv_cp[[i]]$Freq_homopl_norm_size) > 0) & (as.integer(csv_cp[[i]]$is_clustered)==1), ]
		write.csv(csv_cp[[i]], glue("stat_results/{out_folder}/stat_ready_csvs_only_clustered/df_only_clustered_{path_thresholds[i]}_SYN.csv"), quote=F, row.names=F)

		# NON-SYN
		joined_mut_sites_clustered_non_syn[[i]] <- joined_mut_sites_clustered[[i]][joined_mut_sites_clustered[[i]]$syn == 0,]
		joined_mut_sites_clustered_non_syn[[i]]$major_lineage <- as.factor( joined_mut_sites_clustered_non_syn[[i]]$major_lineage )
		joined_mut_sites_clustered_non_syn[[i]]$spec_regions <- as.factor(joined_mut_sites_clustered_non_syn[[i]]$spec_regions)
		csv_cp2[[i]] <- joined_mut_sites_clustered_non_syn[[i]][ order(joined_mut_sites_clustered_non_syn[[i]]$Freq_homopl, decreasing=T) ,]
		write.csv(csv_cp2[[i]], glue("stat_results/{out_folder}/stat_ready_csvs_all_sites/df_all_sites_{path_thresholds[i]}_NON_SYN.csv"), quote=F, row.names=F)
		csv_cp2[[i]] <- csv_cp2[[i]][(as.integer(csv_cp2[[i]]$Freq_homopl) > 0) & (as.integer(csv_cp2[[i]]$is_clustered)==1), ]
		write.csv(csv_cp2[[i]], glue("stat_results/{out_folder}/stat_ready_csvs_only_clustered/df_only_clustered_{path_thresholds[i]}_NON_SYN.csv"), quote=F, row.names=F)
	}
	
	#View(clustered_dfs[[1]])
	#View(joined_mut_sites_clustered[[1]])
	
	names(joined_mut_sites_clustered_syn) <- table_combs
	names(joined_mut_sites_clustered_non_syn) <- table_combs
	
	# Stats
	# sink(file = glue("stat_results/{out_folder}/lm_output_SPLITTED_ALL.txt"))
	# for(i in 1:length(path_thresholds)) {
	# 	print(glue("{path_stats}/{path_thresholds[i]}/"))
	# 
	# 	#print(str(joined_mut_sites_clustered_syn[[i]]))
	# 
	# 	# aa_length == NA (i.e., linearly related to the other variables) for all SYN and NON-SYN on poisson (then 'aa_length' dropped for poisson, kept for logistic)
	# 	# indep_found_pos_selection == NA for all SYN (then dropped 'indep_found_pos_selection' for poisson and logistic on SYN models)
	# 
	# 	# Summary stats
	# 	library(rstatix)
	# 	syn_lin_summ <- joined_mut_sites_clustered_syn[[i]] %>% group_by(major_lineage) %>% get_summary_stats(Freq_homopl, type = "mean_sd") #is_clustered
	# 	print("Summary major_lineage and is_clustered (SYN): ")
	# 	print(syn_lin_summ, n=20)
	# 	syn_reg_summ <- joined_mut_sites_clustered_syn[[i]] %>% group_by(spec_regions) %>% get_summary_stats(Freq_homopl, type = "mean_sd") #is_clustered
	# 	print("Summary spec_regions and is_clustered (SYN): ")
	# 	print(syn_reg_summ, n=60)
	# 	detach(package:rstatix,unload=TRUE)
	# 	
	# 	print("1111")
	# 	
	# 	# poisson regression (removed is_clustered)
	# 	poisson_models_syn[[i]] <- glm(Freq_homopl ~ major_lineage + spec_regions + -1, data=joined_mut_sites_clustered_syn[[i]], na.action=na.omit, family=poisson(link = "log")) #genomic_index +s_mut_region_interest +syn + aa_length + indep_found_pos_selection
	# 	poisson_models_syn_ci[[i]] <- confint(poisson_models_syn[[i]], level=0.95)
	# 	print("Summary poisson (SYN): ")
	# 	#print(nrow(joined_mut_sites_clustered_syn[[i]]))
	# 	print(summary(poisson_models_syn[[i]]))
	# 	print("Summary poisson (SYN) CIs for coeffs: ")
	# 	print(poisson_models_syn_ci[[i]])
	# 	
	# 	library(rstatix)
	# 	nonsyn_lin_summ <- joined_mut_sites_clustered_non_syn[[i]] %>% group_by(major_lineage) %>% get_summary_stats(Freq_homopl, type = "mean_sd") #is_clustered
	# 	print("Summary major_lineage and is_clustered (NON-SYN): ")
	# 	print(nonsyn_lin_summ, n=20)
	# 	nonsyn_reg_summ <- joined_mut_sites_clustered_non_syn[[i]] %>% group_by(spec_regions) %>% get_summary_stats(Freq_homopl, type = "mean_sd") #is_clustered
	# 	print("Summary spec_regions and is_clustered (NON-SYN): ")
	# 	print(nonsyn_reg_summ, n=60)
	# 	detach(package:rstatix,unload=TRUE)
	# 	
	# 	print("2222")
	# 
	# 	#removed is_clustered
	# 	poisson_models_non_syn[[i]] <- glm(Freq_homopl ~ major_lineage + spec_regions + indep_found_pos_selection + -1, data=joined_mut_sites_clustered_non_syn[[i]], na.action=na.omit, family=poisson(link = "log")) #genomic_index +s_mut_region_interest +syn + aa_length
	# 	poisson_models_non_syn_ci[[i]] <- confint(poisson_models_non_syn[[i]], level=0.95)
	# 	print("Summary poisson (NON-SYN): ")
	# 	#print(nrow(joined_mut_sites_clustered_non_syn[[i]]))
	# 	print(summary(poisson_models_non_syn[[i]]))
	# 	print("Summary poisson (NON-SYN) CIs for coeffs: ")
	# 	print(poisson_models_non_syn_ci[[i]])
	# 	# 
	# 	# # Contingency tables
	# 	# conting_syn_lin <- table(joined_mut_sites_clustered_syn[[i]]$major_lineage, joined_mut_sites_clustered_syn[[i]]$is_clustered)
	# 	# print("Contingency table of major_lineage and is_clustered (SYN): ")
	# 	# print(conting_syn_lin)
	# 	# conting_syn_reg <- table(joined_mut_sites_clustered_syn[[i]]$spec_regions, joined_mut_sites_clustered_syn[[i]]$is_clustered)
	# 	# print("Contingency table of spec_regions and is_clustered (SYN): ")
	# 	# print(conting_syn_reg)
	# 	# 
	# 	# # is_clustered_h3
	# 	# conting_syn_lin_h3 <- table(joined_mut_sites_clustered_syn[[i]]$major_lineage, joined_mut_sites_clustered_syn[[i]]$is_clustered_h3)
	# 	# print("Contingency table of major_lineage and is_clustered_h3 (SYN): ")
	# 	# print(conting_syn_lin_h3)
	# 	# conting_syn_reg_h3 <- table(joined_mut_sites_clustered_syn[[i]]$spec_regions, joined_mut_sites_clustered_syn[[i]]$is_clustered_h3)
	# 	# print("Contingency table of spec_regions and is_clustered_h3 (SYN): ")
	# 	# print(conting_syn_reg_h3)
	# 	# 
	# 	# print("3333")
	# 
	# 	# logistic regression (NOT USING NOW: is_clustered_h3 == 1 if Freq_homopl > 3, == 0 otherwise)
	# 	# lr_models_syn[[i]] <- glm(is_clustered==1 ~ major_lineage + aa_length + spec_regions + aa_length + -1, data=joined_mut_sites_clustered_syn[[i]], na.action=na.omit, family=binomial(link="logit")) # Freq_homopl:aa_length + indep_found_pos_selection
	# 	# #lr_models_syn_ci[[i]] <- confint(lr_models_syn[[i]], level=0.95)
	# 	# print("Summary logistic (SYN): ")
	# 	# #print(nrow(joined_mut_sites_clustered_syn[[i]]))
	# 	# print(summary(lr_models_syn[[i]]))
	# 	# # print("Summary logistic (SYN) CIs for coeffs: ")
	# 	# # print(lr_models_syn_ci[[i]])
	# 
	# 	# # Contingency tables
	# 	# conting_nonsyn_lin <- table(joined_mut_sites_clustered_non_syn[[i]]$major_lineage, joined_mut_sites_clustered_non_syn[[i]]$is_clustered)
	# 	# print("Contingency table of major_lineage and is_clustered (NON-SYN): ")
	# 	# print(conting_nonsyn_lin)
	# 	# conting_nonsyn_reg <- table(joined_mut_sites_clustered_non_syn[[i]]$spec_regions, joined_mut_sites_clustered_non_syn[[i]]$is_clustered)
	# 	# print("Contingency table of spec_regions and is_clustered (NON-SYN): ")
	# 	# print(conting_nonsyn_reg)
	# 	# 
	# 	# # is_clustered_h3
	# 	# conting_nonsyn_lin_h3 <- table(joined_mut_sites_clustered_non_syn[[i]]$major_lineage, joined_mut_sites_clustered_non_syn[[i]]$is_clustered_h3)
	# 	# print("Contingency table of major_lineage and is_clustered_h3 (NON-SYN): ")
	# 	# print(conting_nonsyn_lin_h3)
	# 	# conting_nonsyn_reg_h3 <- table(joined_mut_sites_clustered_non_syn[[i]]$spec_regions, joined_mut_sites_clustered_non_syn[[i]]$is_clustered_h3)
	# 	# print("Contingency table of spec_regions and is_clustered_h3 (NON-SYN): ")
	# 	# print(conting_nonsyn_reg_h3)
	# 	# 
	# 	# print("4444")
	# 	
	# 	# NOT USING NOW: is_clustered_h3
	# 	# lr_models_non_syn[[i]] <- glm(is_clustered==1 ~ major_lineage + spec_regions + aa_length + indep_found_pos_selection + -1, data=joined_mut_sites_clustered_non_syn[[i]], na.action=na.omit, family=binomial(link="logit"))
	# 	# #lr_models_non_syn_ci[[i]] <- confint(lr_models_non_syn[[i]], level=0.95)
	# 	# print("Summary logistic (NON-SYN): ")
	# 	# #print(nrow(joined_mut_sites_clustered_non_syn[[i]]))
	# 	# print(summary(lr_models_non_syn[[i]]))
	# 	# # print("Summary logistic (NON-SYN) CIs for coeffs: ")
	# 	# # print(lr_models_non_syn_ci[[i]])
	# 
	# }
	# sink(file = NULL)
	
	# # Plots
	# for(i in 1:length(path_thresholds)) {
	# 	ggplot(joined_mut_sites_clustered_non_syn[[i]], aes(is_clustered_h3)) + geom_bar(fill = "#0073C2FF")
	# 	ggsave(glue("stat_results/{out_folder}/freqs_clustered/logistic_nonsyn_{i}.png"), width=4, height=3, dpi=600, bg="white")
	# 
	# 	# Boxplots
	# 	joined_mut_sites_clustered_syn[[i]]$is_clustered <- as.factor(joined_mut_sites_clustered_syn[[i]]$is_clustered)
	# 	system(glue("mkdir -p stat_results/{out_folder}/boxplots/"))
	# 	ggboxplot(joined_mut_sites_clustered_syn[[i]], x = "major_lineage", y = "Freq_homopl", color="is_clustered", palette="jco") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	# 	ggsave(glue("stat_results/{out_folder}/boxplots/syn_lineage_freq_{i}.png"), width=6, height=4, dpi=600, bg="white")
	# 	ggboxplot(joined_mut_sites_clustered_syn[[i]], x = "spec_regions", y = "Freq_homopl", color="is_clustered", palette="jco") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	# 	ggsave(glue("stat_results/{out_folder}/boxplots/syn_region_freq_{i}.png"), width=8, height=6, dpi=600, bg="white")
	# 
	# 	joined_mut_sites_clustered_non_syn[[i]]$is_clustered <- as.factor(joined_mut_sites_clustered_non_syn[[i]]$is_clustered)
	# 	ggboxplot(joined_mut_sites_clustered_non_syn[[i]], x = "major_lineage", y = "Freq_homopl", color="is_clustered", palette="jco") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	# 	ggsave(glue("stat_results/{out_folder}/boxplots/nonsyn_lineage_freq_{i}.png"), width=6, height=4, dpi=600, bg="white")
	# 	ggboxplot(joined_mut_sites_clustered_non_syn[[i]], x = "spec_regions", y = "Freq_homopl", color="is_clustered", palette="jco") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	# 	ggsave(glue("stat_results/{out_folder}/boxplots/nonsyn_region_freq_{i}.png"), width=8, height=6, dpi=600, bg="white")
	# }
	
	# Inspection plots
	df_syn_clustered_homopl <- rbindlist(joined_mut_sites_clustered_syn, use.names=T, idcol="period_threshold")
	df_non_syn_clustered_homopl <- rbindlist(joined_mut_sites_clustered_non_syn, use.names=T, idcol="period_threshold")
	
	df_syn_clustered_homopl <- change_ids(df_syn_clustered_homopl)
	df_non_syn_clustered_homopl <- change_ids(df_non_syn_clustered_homopl)
	
	genomewise_plot(df_syn_clustered_homopl, mut_type="syn")
	genomewise_plot(df_non_syn_clustered_homopl, mut_type="non-syn")
	
	violin_boxplot_freq_variations(df_syn_clustered_homopl, mut_type="syn")
	violin_boxplot_freq_variations(df_non_syn_clustered_homopl, mut_type="non-syn")
	
	stacked_bar_freqs(df_syn_clustered_homopl, mut_type="syn")
	stacked_bar_freqs(df_non_syn_clustered_homopl, mut_type="non-syn")
	
	system(glue("mkdir -p stat_results/{out_folder}/only_clustered_ordered_homopl/"))
	# SYN
	# df_syn_clustered_homopl <- change_ids(df_syn_clustered_homopl)
	df_syn_clustered_homopl <- df_syn_clustered_homopl[df_syn_clustered_homopl$Freq_homopl >= 5,] # Getting only homoplasies happening at least 5 times
	df_syn_clustered_homopl <- df_syn_clustered_homopl[order(df_syn_clustered_homopl$defining_mut, decreasing=T),]
	#hist(df_syn_clustered_homopl$Freq_homopl)
	write.csv(df_syn_clustered_homopl, glue("stat_results/{out_folder}/only_clustered_ordered_homopl/syn_clustered_homopl_ordered.csv"), quote=F, row.names=F)
	# NON-SYN
	# df_non_syn_clustered_homopl <- change_ids(df_non_syn_clustered_homopl)
	df_non_syn_clustered_homopl <- df_non_syn_clustered_homopl[df_non_syn_clustered_homopl$Freq_homopl >= 5,]
	df_non_syn_clustered_homopl <- df_non_syn_clustered_homopl[order(df_non_syn_clustered_homopl$defining_mut, decreasing=T),]
	#hist(df_non_syn_clustered_homopl$Freq_homopl)
	write.csv(df_non_syn_clustered_homopl, glue("stat_results/{out_folder}/only_clustered_ordered_homopl/non_syn_clustered_homopl_ordered.csv"), quote=F, row.names=F)

	# Counts number of rows for each index
	df_syn_clustered_homopl_counts <- df_syn_clustered_homopl %>% count(period_threshold)
	df_non_syn_clustered_homopl_counts <- df_non_syn_clustered_homopl %>% count(period_threshold)

	# distribution of freqs per genomic region (JUST TO USE FOR BUBBLE PLOTS, FOLLOWING FUNCTIONS WILL BE ONLY FOR CLUSTERED)
	df_syn_clustered_homopl_reg_plots_all <- plot_distr_freqs_spec_regions(df_syn_clustered_homopl,"syn",2,all_clust=TRUE)
	df_non_syn_clustered_homopl_reg_plots_all <- plot_distr_freqs_spec_regions(df_non_syn_clustered_homopl,"non-syn",2,all_clust=TRUE)

	df_syn_clustered_homopl_reg_plots <- plot_distr_freqs_spec_regions(df_syn_clustered_homopl,"syn",2,all_clust=FALSE)
	df_non_syn_clustered_homopl_reg_plots <- plot_distr_freqs_spec_regions(df_non_syn_clustered_homopl,"non-syn",2,all_clust=FALSE)

	# distribution of freqs per major_lineage
	df_syn_clustered_homopl_lin_plots <- plot_distr_freqs_major_lineage(df_syn_clustered_homopl,"syn",2,all_clust=FALSE)
	df_non_syn_clustered_homopl_lin_plots <- plot_distr_freqs_major_lineage(df_non_syn_clustered_homopl,"non-syn",2,all_clust=FALSE)

	system(glue("mkdir -p stat_results/{out_folder}/higher_eq10perc_syn/"))
	system(glue("mkdir -p stat_results/{out_folder}/higher_eq10perc_non_syn/"))
	system(glue("mkdir -p stat_results/{out_folder}/lower_eq5perc_syn/"))
	system(glue("mkdir -p stat_results/{out_folder}/lower_eq5perc_non_syn/"))

	most_common_homopl_syn_period <- most_common_homopl_nonsyn_period <- list()
	for(k in 1:length(major_lineages)) {
		# SYN
		most_common_homopl_syn_period[[k]] <- get_most_common_homopls_lineage_interest(df_syn_clustered_homopl_reg_plots, period_interest=PERIOD_INTEREST, major_lineages[k])
		write.csv(most_common_homopl_syn_period[[k]][[1]], file=glue("stat_results/{out_folder}/higher_eq10perc_syn/{major_lineages[k]}_top_muts.csv"), quote=F, row.names=F)
		write.csv(most_common_homopl_syn_period[[k]][[2]], file=glue("stat_results/{out_folder}/lower_eq5perc_syn/{major_lineages[k]}_top_muts.csv"), quote=F, row.names=F)
		# NON-SYN
		most_common_homopl_nonsyn_period[[k]] <- get_most_common_homopls_lineage_interest(df_non_syn_clustered_homopl_reg_plots, period_interest=PERIOD_INTEREST, major_lineages[k])
		write.csv(most_common_homopl_nonsyn_period[[k]][[1]], file=glue("stat_results/{out_folder}/higher_eq10perc_non_syn/{major_lineages[k]}_top_muts.csv"), quote=F, row.names=F)
		write.csv(most_common_homopl_nonsyn_period[[k]][[2]], file=glue("stat_results/{out_folder}/lower_eq5perc_non_syn/{major_lineages[k]}_top_muts.csv"), quote=F, row.names=F)
		# TODO add Omicron_BA.* for full period?
	}

	# TOP50 by absolute frequency
	most_common_homopl_abs_freq_syn <- get_most_common_homopls_abs_freq(df_syn_clustered_homopl_reg_plots)
	system(glue("mkdir -p stat_results/{out_folder}/abs_freq/"))
	write.csv(most_common_homopl_abs_freq_syn, file=glue("stat_results/{out_folder}/abs_freq/syn_top_muts.csv"), quote=F, row.names=F)
	most_common_homopl_abs_freq_non_syn <- get_most_common_homopls_abs_freq(df_non_syn_clustered_homopl_reg_plots)
	write.csv(most_common_homopl_abs_freq_non_syn, file=glue("stat_results/{out_folder}/abs_freq/nonsyn_top_muts.csv"), quote=F, row.names=F)

	most_common_homopl_s_syn <- get_most_common_homopls_s_regions(df_syn_clustered_homopl_reg_plots)
	system(glue("mkdir -p stat_results/{out_folder}/abs_freq_spike/"))
	write.csv(most_common_homopl_s_syn, file=glue("stat_results/{out_folder}/abs_freq_spike/syn_top_spike_muts.csv"), quote=F, row.names=F)
	most_common_homopl_s_non_syn <- get_most_common_homopls_s_regions(df_non_syn_clustered_homopl_reg_plots)
	write.csv(most_common_homopl_s_non_syn, file=glue("stat_results/{out_folder}/abs_freq_spike/nonsyn_top_spike_muts.csv"), quote=F, row.names=F)

	# is_clustered = 1 (identified by quantile thresholds of stats)
	most_common_homopl_syn_clust1 <- get_most_common_homopls(df_syn_clustered_homopl_reg_plots_all,"syn_clust1",is_clustered_flag=TRUE)
	most_common_homopl_non_syn_clust1 <- get_most_common_homopls(df_non_syn_clustered_homopl_reg_plots_all, "non_syn_clust1",is_clustered_flag=TRUE)
	#View(most_common_homopl_non_syn_clust1[[2]])
	# 
	# # is_clustered = 0 (NOT identified by quantile thresholds of stats)
	# most_common_homopl_syn_clust0 <- get_most_common_homopls(df_syn_clustered_homopl_reg_plots_all,"syn_clust0",is_clustered_flag=FALSE)
	# most_common_homopl_non_syn_clust0 <- get_most_common_homopls(df_non_syn_clustered_homopl_reg_plots_all, "non_syn_clust0",is_clustered_flag=FALSE)
	# #View(most_common_homopl_non_syn_clust0[[2]])
	#
	# Bubbles for is_clustered = 1
	bubble_syn_clust1 <- bubble_plots(most_common_homopl_syn_clust1,"syn_clust1")
	bubble_non_syn_clust1 <- bubble_plots(most_common_homopl_non_syn_clust1,"non-syn_clust1")

	# # Bubbles for is_clustered = 0
	# bubble_syn_clust0 <- bubble_plots(most_common_homopl_syn_clust0,"syn_clust0")
	# bubble_non_syn_clust0 <- bubble_plots(most_common_homopl_non_syn_clust0,"non-syn_clust0")
	
}

# FOR PERIOD 2
# Expected paths containing results for each quantile threshold
thr_pref <- "threshold_quantile"
path_thresholds <- c(glue("{thr_pref}_0.25"), glue("{thr_pref}_0.5"), glue("{thr_pref}_0.75"), glue("{thr_pref}_1"), glue("{thr_pref}_2"),
																					glue("{thr_pref}_3"), glue("{thr_pref}_4"), glue("{thr_pref}_5"), glue("{thr_pref}_10"), glue("{thr_pref}_25")) #, glue("{thr_pref}_50")
table_names_periods <- c(2)
table_names_thresholds <- c(0.25,0.5,0.75,1,2,3,4,5,10,25)
table_combs <- do.call(paste, c(expand.grid(table_names_periods, table_names_thresholds), sep="_"))
out_folder <- "period2"
PERIOD_INTEREST <- 2
major_lineages <- c("EU1_B.1.177","Alpha_B.1.1.7","Delta_AY.4.*", "Delta_other","Other") #"Beta_B.1.351","Gamma_P.1"

clustering_model_fitting("results/02_sc2_root_to_nov2021","period2")

# FOR PERIOD 3
thr_pref <- "threshold_quantile"
path_thresholds <- c(glue("{thr_pref}_0.25"), glue("{thr_pref}_0.5"), glue("{thr_pref}_0.75"), glue("{thr_pref}_1"), glue("{thr_pref}_2"), glue("{thr_pref}_3"), glue("{thr_pref}_4")) #, glue("{thr_pref}_4"), glue("{thr_pref}_5"), glue("{thr_pref}_10"), glue("{thr_pref}_25"))
table_names_periods <- c(3)
table_names_thresholds <- c(0.25,0.5,0.75,1,2,3,4)#,5,10,25) # TODO update based on thresholds
table_combs <- do.call(paste, c(expand.grid(table_names_periods, table_names_thresholds), sep="_"))
out_folder <- "period3"
PERIOD_INTEREST <- 3
major_lineages <- c("EU1_B.1.177","Alpha_B.1.1.7","Delta_AY.4.*", "Delta_other","Omicron_BA.1.*","Omicron_BA.2.*","Other") #"Beta_B.1.351","Gamma_P.1","Omicron_BA.*",

clustering_model_fitting("results/03_sc2_whole_period","period3")

