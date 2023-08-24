libs_load <- c("glue","stringr", "dplyr","data.table","ggplot2","tidyr","tidyverse","ggrepel","ggforce", "viridis", "ggpubr", "htmlwidgets", "plotly", "olsrr")
invisible( lapply(libs_load, library, character.only=TRUE) )

# Function to facilitate comparison between periods and across thresholds
#' Title
#'
#' @param path_stats 
#' @param rm_freq_outliers 
#' @param out_folder 
#' @param out_poisson_file 
#'
#' @return
#' @export
#'
#' @examples
stats_multiple_thresholds <- function(path_stats, rm_freq_outliers=TRUE, out_folder, out_poisson_file) {
	
	if(dir.exists(path_stats))
		path <- path_stats
	else
		stop("Directory does not exist!")
	
	system(glue("mkdir -p stat_results/{out_folder}/"))

	# Change period_threshold to match correct period and threshold (0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 10, 25)
	.change_ids <- function(df) {
		for(i in 1:nrow(df)) {
			if(df[[i,1]] %in% names(table_combs)) {
				idx <- match(df[[i,1]],names(table_combs))
				df[[i,1]] <- table_combs[idx]
			}
		}
		return(df)
	}
	
	.plot_distr_freqs_spec_regions <- function(df,aa_change,remove_if_x_rows_period_thr=2, all_clust=FALSE) {
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
			# Normalize frequency of homoplasy by size of genomic region
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
			ggplot(data=df_splitted_mean_norm_size[[j]], aes(x=spec_regions, y=mean_freq_homopl_norm, fill=spec_regions)) + geom_bar(stat="identity") + labs(x="Genomic region",y="Mean homoplasy Freq normalized by genomic region") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
			ggsave(glue("stat_results/{out_folder}/plot_mut_mean_freqs_{aa_change}_regions/{names(df_splitted_mean)[j]}_mean_freq_norm.png"), width=6, height=5, dpi=300, bg="white")
		}
		
		# Count of homoplasies per genomic region
		df_splitted_count_muts <- df_splitted
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
	
	.plot_distr_freqs_major_lineage <- function(df,aa_change,remove_if_x_rows_period_thr=2,all_clust=FALSE) {
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
			# Normalize frequency of homoplasy by size of genomic region
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
			ggplot(data=df_splitted_mean_norm_size[[j]], aes(x=major_lineage, y=mean_freq_homopl_norm, fill=major_lineage)) + geom_bar(stat="identity") + labs(x="Major lineage",y="Mean homoplasy Freq normalized by genomic region") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
			ggsave(glue("stat_results/{out_folder}/plot_mut_mean_freqs_{aa_change}_lineage/{names(df_splitted_mean)[j]}_mean_freq_norm.png"), width=6, height=5, dpi=300, bg="white")
		}
		#View(df_splitted_mean[[1]])
		
		# Count of homoplasies per major_lineage
		df_splitted_count_muts <- df_splitted
		for(k in 1:length(df_splitted_count_muts)) {
			df_splitted_count_muts[[k]] <- df_splitted_count_muts[[k]] %>% group_by(major_lineage) %>% summarise(count_muts=n()) #%>% mutate(count_norm=count_muts/aa_length)
			ggplot(data=df_splitted_count_muts[[k]], aes(x=major_lineage, y=count_muts, fill=major_lineage)) + geom_bar(stat="identity") + labs(x="Major lineage",y="# mutations") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
			ggsave(glue("stat_results/{out_folder}/plot_mut_counts_{aa_change}_lineage/{names(df_splitted_count_muts)[k]}_mut_counts.png"), width=6, height=5, dpi=300, bg="white")
		}
		
		return(df_splitted)
	}
	
	.get_most_common_homopls_lineage_interest <- function(df, period_interest, lineage_interest) {
		df_homopl_binded <- rbindlist(df)
		
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
		filt_high_thresholds <- filt_high_thresholds %>% filter(major_lineage==lineage_interest & threshold>=10) %>% slice_max(n=30, order_by=Freq_homopl) %>% arrange(period,threshold,desc(as.numeric(Freq_homopl)))
		# Get TOP-50 frequent homoplasies when threshold <= 5 (n=8)
		filt_low_thresholds <- filt_low_thresholds %>% filter(major_lineage==lineage_interest & threshold<=5) %>% slice_max(n=30, order_by=Freq_homopl) %>% arrange(period,threshold,desc(as.numeric(Freq_homopl)))
		
		return(list(filt_high_thresholds, filt_low_thresholds))
	}
	
	.get_most_common_homopls_abs_freq <- function(df) {
		df_homopl_binded <- rbindlist(df)
		
		df_homopl_binded <- df_homopl_binded %>% separate(period_threshold, c("period","threshold"), sep="_")
		df_homopl_binded$period <- as.integer(df_homopl_binded$period)
		df_homopl_binded$threshold <- as.numeric(df_homopl_binded$threshold)
		df_homopl_binded$Freq_homopl <- as.integer(df_homopl_binded$Freq_homopl)
		most_common_homopl_abs_freq <- df_homopl_binded[!duplicated(df_homopl_binded$defining_mut),] # Removing duplicates (period and threshold)
		most_common_homopl_abs_freq <- most_common_homopl_abs_freq %>% slice_max(n=30, order_by=Freq_homopl) %>% arrange(desc(as.numeric(Freq_homopl)),period,threshold)
		return(most_common_homopl_abs_freq)
	}
	
	# TOP50 only considering S regions
	.get_most_common_homopls_s_regions <- function(df) {
		df_homopl_binded <- rbindlist(df)
		
		df_homopl_binded <- df_homopl_binded %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
		df_homopl_binded$period <- as.integer(df_homopl_binded$period)
		df_homopl_binded$threshold <- as.numeric(df_homopl_binded$threshold)
		df_homopl_binded$Freq_homopl <- as.integer(df_homopl_binded$Freq_homopl)
		most_common_homopl_s_regions <- df_homopl_binded[!duplicated(df_homopl_binded$defining_mut),] # Removing duplicates (period and threshold)
		most_common_homopl_s_regions <- most_common_homopl_s_regions %>% filter(protein=="S") %>% slice_max(n=30, order_by=Freq_homopl) %>% arrange(period,threshold,desc(as.numeric(Freq_homopl)))
		return(most_common_homopl_s_regions)
	}
	
	# Ordered list of homoplasies for suitable combinations of period and threshold (IMPORTANT: using this now) - previously TOP-10
	.get_most_common_homopls <- function(df_list, out_suffix, is_clustered_flag=TRUE) {
		df_homopl_binded <- rbindlist(df_list)
		if(is_clustered_flag) {
			df_homopl_binded <- df_homopl_binded[df_homopl_binded$is_clustered==1,]
		}else {
			df_homopl_binded <- df_homopl_binded[df_homopl_binded$is_clustered==0,]
		}
		df_homopl_binded <- df_homopl_binded %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
		df_homopl_binded$id <- paste0(df_homopl_binded$period,"_",df_homopl_binded$threshold)
		df_splitted <- split(df_homopl_binded, df_homopl_binded$id)
		df_splitted_tops <- list()
		system(glue("mkdir -p stat_results/{out_folder}/most_common_homopl_{out_suffix}/"))
		# Each combination will have its own file
		for(i in 1:length(df_splitted)) {
			df_splitted_tops[[i]] <- df_splitted[[i]] %>% slice_max(n=30, order_by=Freq_homopl) %>% arrange(period,threshold,desc(as.numeric(as.numeric(Freq_homopl))))
			write.csv(df_splitted_tops[[i]], file=glue("stat_results/{out_folder}/most_common_homopl_{out_suffix}/{unique(df_splitted_tops[[i]]$id)}.csv"), quote=F, row.names=F)
		}
		
		# File where all dfs above are concatenated
		df_splitted_tops_all <- rbindlist(df_splitted_tops)
		write.csv(df_splitted_tops_all, file=glue("stat_results/{out_folder}/most_common_homopl_{out_suffix}/all_most_common_concatenated.csv"), quote=F, row.names=F)
		
		# Including all cols
		df_splitted_tops_all_cols <- df_splitted_tops_all %>% distinct(defining_mut, .keep_all = TRUE)
		
		# if `defining_mut` duplicated for diff periods or thresholds, append it to appropriate column (collapsed using |) and show only 1 row with each mut
		df_splitted_tops_all_no_rep <- df_splitted_tops_all %>% group_by(defining_mut) %>% summarize(period=paste0(na.omit(period), collapse="|"), threshold=paste0(na.omit(threshold), collapse="|"), Freq_homopl=paste0(na.omit(Freq_homopl), collapse="|")) %>% ungroup()
		df_splitted_tops_all_no_rep <- df_splitted_tops_all_no_rep %>% dplyr::select(period, threshold, defining_mut, Freq_homopl)
		df_splitted_tops_all_no_rep <- df_splitted_tops_all_no_rep %>% inner_join(df_splitted_tops_all_cols, by="defining_mut")
		df_splitted_tops_all_no_rep <- df_splitted_tops_all_no_rep[,-c(5:7)]
		write.csv(df_splitted_tops_all_no_rep, file=glue("stat_results/{out_folder}/most_common_homopl_{out_suffix}/all_most_common_concatenated_no_rep.csv"), quote=F, row.names=F)
		
		return(list(df_splitted_tops, df_splitted_tops_all, df_splitted_tops_all_no_rep))
	}
	
	.bubble_plots_replacements_lineages <- function(df, show_replacements=TRUE, out_suffix) {
		
		most_common_lineages_rep <- most_common_lineages_norep <- list()
		most_common_lineages_rep_first <- most_common_lineages_norep_first <- list()
		p4 <- list()
		# For each one of the major_lineage, plot the thresholds met and freq of homoplasy in each case (coloured by analysis where identified) + only first threshold met
		for(i in 1:length(major_lineages)) {
			most_common_lineages_rep[[i]] <- df %>% filter(major_lineage==major_lineages[i])
			most_common_lineages_rep[[i]]$threshold <- factor(most_common_lineages_rep[[i]]$threshold, levels=table_names_thresholds)
			most_common_lineages_norep[[i]] <- most_common_lineages_rep[[i]]
			
			if(show_replacements) {
				p <- ggplot(most_common_lineages_rep[[i]], aes(x=as.factor(threshold), y=defining_mut, size=Freq_homopl)) +
					geom_point(alpha=0.5) + labs(x="Quantile threshold", y="Homoplasy", title=glue("Lineage {major_lineages[i]}")) + theme_minimal() + theme(axis.text.y=element_text(size=5)) +
					scale_size_area(name="Homoplasy frequency")
				ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/repl_thr_lineage_{major_lineages[i]}.png"), plot=p, width=12, height=15, dpi=600, bg="white")
				
				#p_ia <- ggplotly(p)
				#htmlwidgets::saveWidget(p_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/repl_thr_lineage_{major_lineages[i]}.html"))
				
				most_common_lineages_rep_first[[i]] <- most_common_lineages_rep[[i]]
				most_common_lineages_rep_first[[i]] <- most_common_lineages_rep_first[[i]] %>% distinct() %>% group_by(defining_mut) %>% slice_min(n=1, as.numeric(threshold))
				
				p2 <- ggplot(most_common_lineages_rep_first[[i]], aes(x=as.factor(threshold), y=defining_mut, size=Freq_homopl)) +
					geom_point(alpha=0.5) + labs(x="Quantile threshold (first detected only)", y="Homoplasy", title=glue("Lineage {major_lineages[i]}")) + theme_minimal() + theme(axis.text.y=element_text(size=5)) +
					scale_size_area(name="Homoplasy frequency") #+ scale_color_manual(values=c("#7c9885","#e1ad01"), name="Analysis identified")
				ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/repl_first_detected_thr_lineage_{major_lineages[i]}.png"), plot=p2, width=12, height=15, dpi=600, bg="white")
				
				#p2_ia <- ggplotly(p2)
				#htmlwidgets::saveWidget(p2_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/repl_first_detected_thr_lineage_{major_lineages[i]}.html"))
			} else {
				
				p3 <- ggplot(most_common_lineages_norep[[i]], aes(x=as.factor(threshold), y=prot_site, size=Freq_homopl, color=analysis_identified.y)) +
					geom_point(alpha=0.5) + labs(x="Quantile threshold", y="Protein site") + theme_minimal() + #title=glue("Lineage {major_lineages[i]}"))
					theme(axis.text.y = element_text(size=7, color="black"), axis.text.x = element_text(size=8, color="black")) + #theme(axis.text.y=element_text(size=5)) +
					scale_size_area(name="Homoplasy frequency") + scale_color_manual(values=pal_id, name="Analysis identified")
				ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/onlysite_thr_lineage_{major_lineages[i]}.png"), plot=p3, width=12, height=15, dpi=600, bg="white")
				#p3_ia <- ggplotly(p3)
				#htmlwidgets::saveWidget(p3_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/onlysite_thr_lineage_{major_lineages[i]}.html"))
				
				most_common_lineages_norep_first[[i]] <- most_common_lineages_norep[[i]]
				most_common_lineages_norep_first[[i]] <- most_common_lineages_norep_first[[i]] %>% distinct() %>% group_by(prot_site) %>% slice_min(n=1, as.numeric(threshold))
				
				p4[[i]] <- ggplot(most_common_lineages_norep_first[[i]], aes(x=as.factor(threshold), y=prot_site, size=Freq_homopl, color=analysis_identified.y)) +
					geom_point(alpha=0.5) + labs(x="Quantile threshold (first detected only)", y="Protein site") + theme_minimal() + #, title=glue("Lineage {major_lineages[i]}") 
					theme(axis.text.y = element_text(size=7, color="black"), axis.text.x = element_text(size=8, color="black")) + #theme(axis.text.y=element_text(size=5)) +
					scale_size_area(name="Homoplasy frequency") + scale_color_manual(values=pal_id, name="Analysis identified")
				ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/onlysite_first_detected_thr_lineage_{major_lineages[i]}.png"), plot=p4[[i]], width=10, height=12, dpi=600, bg="white") #12,15
				#p4_ia <- ggplotly(p4)
				#htmlwidgets::saveWidget(p4_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/onlysite_first_detected_thr_lineage_{major_lineages[i]}.html"))
			}
		}
		return(p4)
	}
	
	.bubble_plots <- function(most_common_homopl_df, mut_type, out_suffix) {
		# prepare dfs for bubble plot selection analysis from HypHy identified mutations
		pos_sel_sites2 <- pos_sel_sites %>% dplyr::select(protein, site, prot_site)
		pos_sel_sites2$analysis_identified <- "HyPhy"
		pos_sel_sites2$Freq_homopl <- NA
		pos_sel_sites2$period <- NA
		pos_sel_sites2$threshold <- NA
		pos_sel_sites2$major_lineage <- NA
		pos_sel_sites2 <- pos_sel_sites2 %>% dplyr::select(period, threshold, protein, site, prot_site, major_lineage, analysis_identified, Freq_homopl)
		
		most_common_homopl_all <- most_common_homopl_df[[2]]
		most_common_homopl_all$analysis_identified <- "mlscluster"
		most_common_homopl_all_adj <- most_common_homopl_all %>% dplyr::select(period, threshold, protein, site, prot_site, major_lineage, analysis_identified, Freq_homopl)
		
		joined_most_common_indep <- rbind(pos_sel_sites2, most_common_homopl_all_adj)
		
		joined_most_common_indep1 <- joined_most_common_indep %>% group_by(prot_site) %>% summarize(analysis_identified=paste0(na.omit(analysis_identified), collapse=";")) %>% ungroup()
		# Make sure there is only one mlscluster for each row (remove if more than one semicolon)
		joined_most_common_indep1$analysis_identified <- sub(sprintf("^((.*?;){%d}.*?);.*", 1), "\\1", joined_most_common_indep1$analysis_identified)
		joined_most_common_indep1$analysis_identified <- ifelse(joined_most_common_indep1$analysis_identified=="mlscluster;mlscluster", "mlscluster", joined_most_common_indep1$analysis_identified)
		joined_most_common_indep_ok <- joined_most_common_indep %>% inner_join(joined_most_common_indep1, by="prot_site")
		joined_most_common_indep_ok[joined_most_common_indep_ok == "HyPhy;mlscluster"] <- "Intersection"
		
		if(mut_type=="non-syn") {
			joined_most_common_indep_ok$protein <- factor(joined_most_common_indep_ok$protein, levels=lvls_bubble)
			joined_most_common_indep_ok <- joined_most_common_indep_ok[order(joined_most_common_indep_ok$protein, as.numeric(joined_most_common_indep_ok$site) ),]
			joined_most_common_indep_ok$prot_site <- factor(joined_most_common_indep_ok$prot_site, levels=unique( joined_most_common_indep_ok$prot_site[order(joined_most_common_indep_ok$protein, as.numeric(joined_most_common_indep_ok$site) )] ))
		}
		
		system(glue("mkdir -p stat_results/{out_folder}/comparison_bubble_{out_suffix}/"))
		
		p <- ggplot(joined_most_common_indep_ok, aes(y=analysis_identified.y, x=prot_site, color=analysis_identified.y)) + #color=analysis_identified.y
			geom_point(alpha=0.8) + labs(y="Analyses identified",x="Protein site") + theme_minimal() +
			scale_color_manual(values=pal_id, name="Analysis identified") + theme(legend.position="none", axis.text.x=element_text(size=5, angle=90, vjust = 0.5, hjust=1)) #, aspect.ratio=16/9
		ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/homopl_method_combinations.png"), plot=p, width=15, height=5, dpi=900, bg="white")
		#p_ia <- ggplotly(p)
		#htmlwidgets::saveWidget(p_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/homopl_method_combinations.html"))
		
		joined_most_common_indep_ok_high_freq_lin <- joined_most_common_indep_ok %>% group_by(major_lineage, prot_site) %>% slice_max(Freq_homopl)
		
		# Plot max frequency of homoplasic sites (bubble sizes) for each one of the major_lineage (coloured by analysis where identified)
		joined_most_common_indep_ok_high_freq_lin <- joined_most_common_indep_ok_high_freq_lin[joined_most_common_indep_ok_high_freq_lin$analysis_identified.y != "HyPhy",]
		joined_most_common_indep_ok_high_freq_lin <- joined_most_common_indep_ok_high_freq_lin[!is.na(joined_most_common_indep_ok_high_freq_lin$major_lineage),]
		
		p2 <- ggplot(joined_most_common_indep_ok_high_freq_lin, aes(x=major_lineage, y=prot_site, size=Freq_homopl, color=analysis_identified.y)) +
			geom_point(alpha=0.5) + labs(x="Major lineage", y="Protein site") + theme_minimal() + theme(axis.text.y = element_text(size=5), aspect.ratio=16/9) +
			scale_size_area(name="Highest freq of homopl", breaks=c(20,50,150,500,1300)) + scale_color_manual(values=c("#7c9885","#e1ad01"), name="Analysis identified") #"#033f63"
		ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_all_lineages.png"), plot=p2, width=12, height=15, dpi=600, bg="white")
		#p2_ia <- ggplotly(p2)
		#htmlwidgets::saveWidget(p2_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_all_lineages.html"))
		
		most_common_homopl_all_high_freq_lin <- most_common_homopl_all %>% group_by(major_lineage, defining_mut) %>% slice_max(Freq_homopl)
		most_common_homopl_all_high_freq_lin <- most_common_homopl_all_high_freq_lin[!is.na(most_common_homopl_all_high_freq_lin$major_lineage),]
		# Plot max frequency of each homoplasy across thresholds (including WT and replaced aa) identified (without comparison with HyPhy)
		p3 <- ggplot(most_common_homopl_all_high_freq_lin, aes(x=major_lineage, y=defining_mut, size=Freq_homopl)) + #color=analysis_identified.y
			geom_point(alpha=0.5) + labs(x="Major lineage", y="Homoplasy") + theme_minimal() + theme(axis.text.y = element_text(size=5), aspect.ratio=16/9) +
			scale_size_area(name="Highest freq of homopl") # breaks=c(20,50,150,500,1300) + scale_color_manual(values=c("#7c9885","#e1ad01", "#033f63")) +
		ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_replacement_lineages.png"), plot=p3, width=12, height=15, dpi=600, bg="white")
		#p3_ia <- ggplotly(p3)
		#htmlwidgets::saveWidget(p3_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_replacement_lineages.html"))
		
		# Bubble plots per major lineage comparing with other results (HyPhy): just showing protein site and not actual replacements
		s3_s8 <- .bubble_plots_replacements_lineages(joined_most_common_indep_ok, show_replacements=FALSE, out_suffix)
		# Bubble plots per major lineage  showing actual replacements and no comparison with other analyses
		.bubble_plots_replacements_lineages(most_common_homopl_all, show_replacements=TRUE, out_suffix)
		
		s3_s8
	}
	
	.genomewide_plot <- function(df, mut_type) {
		df <- df %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
		df$threshold <- factor(df$threshold, levels=table_names_thresholds)
		
		if(mut_type == "syn") {
			genomic_ranges_df <- read.csv("config/genomic_ranges_plot_syn.tsv", header=T, sep="\t")
			# No split for SYN
			df_freq_genomic_regions_syn <- df %>% inner_join(genomic_ranges_df, by=c("spec_regions"="new_prot_name"))
			df_freq_genomic_regions_syn_first <- df_freq_genomic_regions_syn %>% distinct() %>% group_by(defining_mut) %>% slice_min(n=1, as.numeric(threshold))
			df_freq_genomic_regions_syn_first <- df_freq_genomic_regions_syn_first[!is.na(df_freq_genomic_regions_syn_first$threshold),]
			
			system(glue("mkdir -p stat_results/{out_folder}/genomewide_plot_{mut_type}/"))
			
			# Get table with all TFP-homoplasies
			df_freq_genomic_regions_syn_first_t <- df_freq_genomic_regions_syn_first[order(as.numeric(df_freq_genomic_regions_syn_first$Freq_homopl), decreasing=TRUE),]
			# Get table grouped by spec_regions and threshold
			df_freq_genomic_regions_syn_first_t2 <- df_freq_genomic_regions_syn_first_t %>% group_by(spec_regions, threshold) %>% summarise(n=n())
			
			write.csv(df_freq_genomic_regions_syn_first_t, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_1st_threshold.csv"), quote=F, row.names=F)
			write.csv(df_freq_genomic_regions_syn_first_t2, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_1st_threshold_REG_THR.csv"), quote=F, row.names=F)
			
			p <- ggplot(df_freq_genomic_regions_syn_first, aes(x=genomic_index, y=Freq_homopl, color=threshold)) +
				geom_point(alpha=0.9) + labs(x="Nucleotide position", y="Frequency in TFP clusters") + theme_minimal() + #theme(axis.text.y=element_text(size=5)) +
				scale_color_manual(values=pal_thresholds, name="Threshold first\ndetected") +
				coord_cartesian(xlim = c(1, 29903), ylim = c(2,10)) + scale_y_continuous(breaks=seq(2,10,by=1))
			ggsave(glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_1st_threshold.png"), plot=p, width=15, height=5, dpi=600, bg="white")
			p_ia <- ggplotly(p)
			htmlwidgets::saveWidget(p_ia, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_1st_threshold.html"))
			
			p2 <- ggplot(df_freq_genomic_regions_syn_first, aes(x=genomic_index, y=Freq_homopl, color=major_lineage)) +
				geom_point(alpha=0.9) + labs(x="Nucleotide position", y="Frequency in TFP clusters") + theme_minimal() + #theme(axis.text.y=element_text(size=5)) +
				scale_color_manual(values=pal_lineages, name="Major lineage") +
				coord_cartesian(xlim = c(0, 29903), ylim = c(2,10)) + scale_y_continuous(breaks=seq(2,10,by=1))
			ggsave(glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_major_lineage.png"), plot=p2, width=15, height=5, dpi=600, bg="white")
			p2_ia <- ggplotly(p2)
			htmlwidgets::saveWidget(p2_ia, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_major_lineage.html"))
			
			return(list(p, p2))
		}else if(mut_type == "non-syn") {
			genomic_ranges_df <- read.csv("config/genomic_ranges_plot_non_syn.tsv", header=T, sep="\t")
			df_freq_genomic_regions_nonsyn <- df %>% inner_join(genomic_ranges_df, by=c("spec_regions"="new_prot_name"))
			# Split because each protein (ORF1AB + others) has its own relative coordinates
			df_freq_genomic_regions_split_nonsyn <- split(df_freq_genomic_regions_nonsyn, df_freq_genomic_regions_nonsyn$protein)
			df_freq_genomic_regions_split_nonsyn_first <- df_freq_genomic_regions_split_nonsyn_first_t <- df_freq_genomic_regions_split_nonsyn_first_t2 <- list()
			p3 <- p4 <- list()
			for(i in 1:length(df_freq_genomic_regions_split_nonsyn)) {
				df_freq_genomic_regions_split_nonsyn_first[[i]] <- df_freq_genomic_regions_split_nonsyn[[i]] %>% distinct() %>% group_by(defining_mut) %>% slice_min(n=1, as.numeric(threshold))
				df_freq_genomic_regions_split_nonsyn_first[[i]] <- df_freq_genomic_regions_split_nonsyn_first[[i]][!is.na(df_freq_genomic_regions_split_nonsyn_first[[i]]$threshold),]
				
				# Get tables with all TFP-homoplasies
				df_freq_genomic_regions_split_nonsyn_first_t[[i]] <- df_freq_genomic_regions_split_nonsyn_first[[i]][order(as.numeric(df_freq_genomic_regions_split_nonsyn_first[[i]]$Freq_homopl), decreasing=TRUE),]
				# Get tables grouped by spec_regions and threshold
				df_freq_genomic_regions_split_nonsyn_first_t2[[i]] <- df_freq_genomic_regions_split_nonsyn_first_t[[i]] %>% group_by(spec_regions, threshold) %>% summarise(n=n())
				
				system(glue("mkdir -p stat_results/{out_folder}/genomewide_plot_{mut_type}/"))
				
				write.csv(df_freq_genomic_regions_split_nonsyn_first_t[[i]], glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_1st_threshold_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.csv"), quote=F, row.names=F)
				write.csv(df_freq_genomic_regions_split_nonsyn_first_t2[[i]], glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_1st_threshold_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}_REG_THR.csv"), quote=F, row.names=F)
				
				p3[[i]] <- ggplot(df_freq_genomic_regions_split_nonsyn_first[[i]], aes(x=genomic_index, y=Freq_homopl, color=threshold)) +
					geom_point(alpha=0.9) + labs(x="Aminoacid position", y="Frequency in TFP clusters") + theme_minimal() +
					scale_color_manual(values=pal_thresholds, name="Threshold first\ndetected") +
					geom_text_repel(aes(label = ifelse(Freq_homopl>5,as.character(defining_mut),'')), size=2, box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=15)+
					coord_cartesian(xlim = c(df_freq_genomic_regions_split_nonsyn_first[[i]]$start[1], df_freq_genomic_regions_split_nonsyn_first[[i]]$end[1])) #, ylim = c(2,10)) + scale_y_continuous(breaks=seq(2,10,by=1))
				ggsave(glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_1st_threshold_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.png"), plot=p3[[i]], width=15, height=5, dpi=600, bg="white")
				p3_ia <- ggplotly(p3[[i]])
				htmlwidgets::saveWidget(p3_ia, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_1st_threshold_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.html"))
				
				p4[[i]] <- ggplot(df_freq_genomic_regions_split_nonsyn_first[[i]], aes(x=genomic_index, y=Freq_homopl, color=major_lineage)) +
					geom_point(alpha=0.9) + labs(x="Aminoacid position", y="Frequency in TFP clusters") + theme_minimal() +
					scale_color_manual(values=pal_lineages, name="Major lineage") +
					geom_text_repel(aes(label = ifelse(Freq_homopl>2,as.character(defining_mut),'')), size=2, box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=15)+
					coord_cartesian(xlim = c(df_freq_genomic_regions_split_nonsyn_first[[i]]$start[1], df_freq_genomic_regions_split_nonsyn_first[[i]]$end[1])) #, ylim = c(2,10)) + scale_y_continuous(breaks=seq(2,10,by=1))
				ggsave(glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_major_lineage_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.png"), plot=p4[[i]], width=15, height=5, dpi=600, bg="white")
				p4_ia <- ggplotly(p4[[i]])
				htmlwidgets::saveWidget(p4_ia, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_major_lineage_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.html"))
				
			}
			return(list(p3, p4))
			#View(df_freq_genomic_regions_split_nonsyn_first[[1]])
		} else {
			stop("Wrong choice of mut_type: options are `syn` and `non-syn`")
		}
		
	}
	
	.violin_boxplot_freq_variations <- function(df, mut_type) {
		df <- df %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
		df$threshold <- factor(df$threshold, levels=table_names_thresholds)
		
		df_freq_genomic_regions_first <- df %>% distinct() %>% group_by(defining_mut) %>% slice_min(n=1, as.numeric(threshold))
		
		if(mut_type == "syn") {
			df_freq_genomic_regions_first$spec_regions <- sub("SYNSNP:","",df_freq_genomic_regions_first$spec_regions)
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
			geom_violin(width=1.4) + coord_cartesian(ylim=c(2,10)) +
			scale_fill_viridis(discrete = TRUE) + theme_minimal() + theme(legend.position="none",axis.text.x = element_text(size=6)) +
			labs(x="Genomic region", y="Frequency in TFP clusters")
		ggsave(glue("stat_results/{out_folder}/violin_boxplot_freq_{mut_type}/variation_freq_regions.png"), width=15, height=10, dpi=600, bg="white")
		
		sample_size_lineages <- df_freq_genomic_regions_first %>% group_by(major_lineage) %>% summarize(num=n())
		
		#inner_join(sample_size_lineages) %>% mutate(myaxis = paste0(major_lineage, "\n", "n=", num)) %>%
		df_freq_genomic_regions_first %>%
			ggplot( aes(x=major_lineage, y=Freq_homopl, fill=major_lineage)) + #myaxis
			geom_violin(width=1.4) + coord_cartesian(ylim=c(2,10)) +
			scale_fill_viridis(discrete = TRUE) + theme_minimal() + theme(legend.position="none",axis.text.x = element_text(size=10)) +
			labs(x="Genomic region", y="Frequency in TFP clusters")
		ggsave(glue("stat_results/{out_folder}/violin_boxplot_freq_{mut_type}/variation_freq_lineages.png"), width=15, height=10, dpi=600, bg="white")
		
	}
	
	.stacked_bar_freqs <- function(df, mut_type) {
		df <- df %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
		df$threshold <- factor(df$threshold, levels=table_names_thresholds)
		
		df_freq_genomic_regions_first <- df %>% distinct() %>% group_by(defining_mut) %>% slice_min(n=1, as.numeric(threshold))
		df_freq_genomic_regions_first <- df_freq_genomic_regions_first[!is.na(df_freq_genomic_regions_first$spec_regions),]
		
		df_freq_genomic_regions_first$Freq_norm <- df_freq_genomic_regions_first$Freq_homopl / df_freq_genomic_regions_first$aa_length
		
		if(mut_type == "syn") {
			df_freq_genomic_regions_first$spec_regions <- sub(".*:","",df_freq_genomic_regions_first$spec_regions)
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
			labs(x="Count of unique homoplasies (normalized by genomic size)", y="Genomic region")
		ggsave(glue("stat_results/{out_folder}/stacked_bar_uniq_counts_{mut_type}/stacked_bar_uniq_counts_norm.png"), width=8, height=6, dpi=600, bg="white")
		
	}
	
	.total_amount_unique_homoplasies <- function(df, mut_type) {
		df <- df %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
		df$threshold <- as.numeric(df$threshold)
		
		df_freq_genomic_regions_first <- df %>% distinct() %>% group_by(defining_mut) %>% slice_min(n=1, as.numeric(threshold))
		
		if(mut_type == "syn") {
			df_freq_genomic_regions_first$spec_regions <- sub(".*:","",df_freq_genomic_regions_first$spec_regions)
		}
		
		system(glue("mkdir -p stat_results/{out_folder}/total_homopls_{mut_type}/"))
		
		total_homopls <- nrow(df_freq_genomic_regions_first)
		
		total_homopls_region <- df_freq_genomic_regions_first %>% group_by(spec_regions) %>% summarise(total=n())
		write.csv(total_homopls_region, glue("stat_results/{out_folder}/total_homopls_{mut_type}/total_homopls_region.csv"), quote=F, row.names=F)
		
		total_homopls_region_lineage <- df_freq_genomic_regions_first %>% group_by(spec_regions, major_lineage) %>% summarise(total=n())
		write.csv(total_homopls_region_lineage, glue("stat_results/{out_folder}/total_homopls_{mut_type}/total_homopls_region_lineage.csv"), quote=F, row.names=F)
		
		total_homopls_region_lower_thr <- df_freq_genomic_regions_first %>% filter(threshold<=1) %>% group_by(spec_regions) %>% summarise(total=n())
		write.csv(total_homopls_region_lower_thr, glue("stat_results/{out_folder}/total_homopls_{mut_type}/total_homopls_region_lower_thr.csv"), quote=F, row.names=F)
		
		total_homopls_region_higher_thr <- df_freq_genomic_regions_first %>% filter(threshold>=2) %>% group_by(spec_regions) %>% summarise(total=n())
		write.csv(total_homopls_region_higher_thr, glue("stat_results/{out_folder}/total_homopls_{mut_type}/total_homopls_region_higher_thr.csv"), quote=F, row.names=F)
	}
	
	
	clustered_dfs <- joined_mut_sites_clustered <- joined_mut_sites_clustered_syn <- joined_mut_sites_clustered_non_syn <- list()
	poisson_models_syn <- poisson_models_non_syn <- lr_models_syn <- lr_models_non_syn <- list()
	poisson_models_syn_ci <- poisson_models_non_syn_ci <- lr_models_syn_ci <- lr_models_non_syn_ci <- list()
	csv_cp <- csv_cp2 <- list()
	
	for(i in 1:length(path_thresholds)) {
		clustered_dfs[[i]] <- read.csv(glue("{path_stats}/{path_thresholds[i]}/clustered_all_df.csv"), header=T)
	}
	
	# If flag to remove 99% quantile frequency outliers is TRUE, then remove
	if(rm_freq_outliers) {
		clustered_dfs <- remove_homopl_freq_outliers(path_stats)
	}
	
	for(i in 1:length(path_thresholds)) {
		
		clustered_dfs[[i]] <- clustered_dfs[[i]][clustered_dfs[[i]]$is_clustered==1,]
		
		clustered_dfs[[i]]$major_lineage <- ifelse(grepl("Beta_B.1.351", clustered_dfs[[i]]$major_lineage), "Other", as.character(clustered_dfs[[i]]$major_lineage))
		clustered_dfs[[i]]$major_lineage <- ifelse(grepl("Gamma_P.1", clustered_dfs[[i]]$major_lineage), "Other", as.character(clustered_dfs[[i]]$major_lineage))
		
		joined_mut_sites_clustered[[i]] <- clustered_dfs[[i]]
		
		# Extract genome index position
		rgx_genomic_idx <- regexpr(':[A-Z]{1}[0-9]{1,5}[A-Z*]{1}', joined_mut_sites_clustered[[i]]$defining_mut)
		comm_coord <- rep(NA,length(joined_mut_sites_clustered[[i]]$defining_mut))
		comm_coord[rgx_genomic_idx != -1] <- str_sub( regmatches(joined_mut_sites_clustered[[i]]$defining_mut, rgx_genomic_idx), 3, -2)
		comm_coord <- as.integer(comm_coord)
		comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("genomic_index")
		joined_mut_sites_clustered[[i]] <- cbind(joined_mut_sites_clustered[[i]], comm_coord_df)
		
		# Adjust df columns
		joined_mut_sites_clustered[[i]] <- joined_mut_sites_clustered[[i]] %>% dplyr::select(defining_mut, Freq_homopl, is_clustered, major_lineage, syn_non_syn, protein, aa_length, indep_found_pos_selection, genomic_index) #s_mut_region_interest
		colnames(joined_mut_sites_clustered[[i]]) <- c("defining_mut","Freq_homopl","is_clustered","major_lineage","syn","protein","aa_length","indep_found_pos_selection","genomic_index") #s_mut_region_interest
		
		joined_mut_sites_clustered[[i]]$protein <- sub("\\:.*", "", joined_mut_sites_clustered[[i]]$defining_mut)
		
		setDT(joined_mut_sites_clustered[[i]]); setDT(annot_gen_range_positions)
		joined_mut_sites_clustered[[i]] <- joined_mut_sites_clustered[[i]][annot_gen_range_positions, on=c("protein==name_cog_md","genomic_index>=start","genomic_index<=end"), spec_regions := new_prot_name]
		joined_mut_sites_clustered[[i]] <- joined_mut_sites_clustered[[i]][annot_gen_length_positions, on=c("spec_regions==new_prot_name"), aa_length := length]
		
		# syn flag
		joined_mut_sites_clustered[[i]]$syn <- ifelse(joined_mut_sites_clustered[[i]]$syn == as.character("non-syn"), 0, 1)
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
		
		# For logistic regression, is_clustered_h3 will be used and it is true if Freq_homopl is > 3
		joined_mut_sites_clustered[[i]]$is_clustered_h3 <- ifelse( (as.integer(joined_mut_sites_clustered[[i]]$Freq_homopl) > 3) & (joined_mut_sites_clustered[[i]]$is_clustered==1), 1, 0)
		# Normalized counts (for inspection only)
		joined_mut_sites_clustered[[i]]$Freq_homopl_norm_size <- as.integer(joined_mut_sites_clustered[[i]]$Freq_homopl) / as.integer(joined_mut_sites_clustered[[i]]$aa_length)
		
		# Make sure to remove all rows with NAs in key columns
		# Freq_homopl, major_lineage, is_clustered, is_clustered_h3, spec_regions, aa_length, indep_found_pos_selection
		joined_mut_sites_clustered[[i]] <- joined_mut_sites_clustered[[i]][!is.na(joined_mut_sites_clustered[[i]]$Freq_homopl) | !is.na(joined_mut_sites_clustered[[i]]$major_lineage) | !is.na(joined_mut_sites_clustered[[i]]$is_clustered) | !is.na(joined_mut_sites_clustered[[i]]$is_clustered_h3) | !is.na(joined_mut_sites_clustered[[i]]$spec_regions) | !is.na(joined_mut_sites_clustered[[i]]$aa_length) | !is.na(joined_mut_sites_clustered[[i]]$indep_found_pos_selection), ]
		
		joined_mut_sites_clustered[[i]]$site <- as.numeric(joined_mut_sites_clustered[[i]]$site)
		
		# Split into syn and non-syn
		system(glue("mkdir -p stat_results/{out_folder}/stat_ready_csvs_all_sites/"))
		system(glue("mkdir -p stat_results/{out_folder}/stat_ready_csvs_only_clustered/"))
		
		# SYN
		joined_mut_sites_clustered_syn[[i]] <- joined_mut_sites_clustered[[i]][joined_mut_sites_clustered[[i]]$syn == 1,]
		joined_mut_sites_clustered_syn[[i]]$major_lineage <- as.factor( joined_mut_sites_clustered_syn[[i]]$major_lineage )
		joined_mut_sites_clustered_syn[[i]]$spec_regions <- as.factor(joined_mut_sites_clustered_syn[[i]]$spec_regions)
		csv_cp[[i]] <- joined_mut_sites_clustered_syn[[i]][ order(joined_mut_sites_clustered_syn[[i]]$Freq_homopl, decreasing=T) ,]
		write.csv(csv_cp[[i]], glue("stat_results/{out_folder}/stat_ready_csvs_all_sites/df_all_sites_{path_thresholds[i]}_SYN.csv"), quote=F, row.names=F)
		csv_cp[[i]] <- csv_cp[[i]][(as.integer(csv_cp[[i]]$Freq_homopl) > 0) & (as.integer(csv_cp[[i]]$is_clustered)==1), ]
		write.csv(csv_cp[[i]], glue("stat_results/{out_folder}/stat_ready_csvs_only_clustered/df_only_clustered_{path_thresholds[i]}_SYN.csv"), quote=F, row.names=F)

		# NON-SYN
		joined_mut_sites_clustered_non_syn[[i]] <- joined_mut_sites_clustered[[i]][joined_mut_sites_clustered[[i]]$syn == 0,]
		joined_mut_sites_clustered_non_syn[[i]]$major_lineage <- as.factor( joined_mut_sites_clustered_non_syn[[i]]$major_lineage )
		joined_mut_sites_clustered_non_syn[[i]]$spec_regions <- as.factor(joined_mut_sites_clustered_non_syn[[i]]$spec_regions)
		
	}
	
	names(joined_mut_sites_clustered_syn) <- table_combs
	names(joined_mut_sites_clustered_non_syn) <- table_combs
	
	print("Performing poisson regressions to see if lineages and genomic regions enriched for TFPs:")
	
	# Stats
	sink(file = glue("stat_results/plots_paper/{out_poisson_file}"))
	for(i in 1:length(path_thresholds)) {

		poisson_models_syn[[i]] <- glm(Freq_homopl ~ major_lineage + spec_regions + -1, data=joined_mut_sites_clustered_syn[[i]], na.action=na.omit, family=poisson(link = "log")) #genomic_index +s_mut_region_interest +syn + aa_length + indep_found_pos_selection
		poisson_models_syn_ci[[i]] <- confint(poisson_models_syn[[i]], level=0.95)
		#print("Summary poisson (SYN): ")
		#print(nrow(joined_mut_sites_clustered_syn[[i]]))
		#print(summary(poisson_models_syn[[i]]))
		#print("Summary poisson (SYN) CIs for coeffs: ")
		#print(poisson_models_syn_ci[[i]])

		poisson_models_non_syn[[i]] <- glm(Freq_homopl ~ major_lineage + spec_regions + indep_found_pos_selection + -1, data=joined_mut_sites_clustered_non_syn[[i]], na.action=na.omit, family=poisson(link = "log")) #genomic_index +s_mut_region_interest +syn + aa_length
		poisson_models_non_syn_ci[[i]] <- confint(poisson_models_non_syn[[i]], level=0.95)
		print(glue("Summary poisson regression (NON-SYN) for {path_thresholds[i]}%: "))
		#print(nrow(joined_mut_sites_clustered_non_syn[[i]]))
		print(summary(poisson_models_non_syn[[i]]))
		# print("Summary poisson (NON-SYN) CIs for coeffs: ")
		# print(poisson_models_non_syn_ci[[i]])

	}
	sink(file = NULL)

	# Inspection plots
	df_syn_clustered_homopl <- rbindlist(joined_mut_sites_clustered_syn, use.names=T, idcol="period_threshold")
	df_non_syn_clustered_homopl <- rbindlist(joined_mut_sites_clustered_non_syn, use.names=T, idcol="period_threshold")
	
	df_syn_clustered_homopl <- .change_ids(df_syn_clustered_homopl)
	df_non_syn_clustered_homopl <- .change_ids(df_non_syn_clustered_homopl)

	print("Genome-wide plots per genomic region:")
	s13 <- .genomewide_plot(df_syn_clustered_homopl, mut_type="syn")
	s9 <- .genomewide_plot(df_non_syn_clustered_homopl, mut_type="non-syn")

	print("Violin plots of frequencies by lineage and genomic region:")
	.violin_boxplot_freq_variations(df_syn_clustered_homopl, mut_type="syn")
	.violin_boxplot_freq_variations(df_non_syn_clustered_homopl, mut_type="non-syn")

	print("Stacked bar of unique sites under selection by lineage and genomic region:")
	.stacked_bar_freqs(df_syn_clustered_homopl, mut_type="syn")
	.stacked_bar_freqs(df_non_syn_clustered_homopl, mut_type="non-syn")

	print("Tables of unique homoplasies:")
	.total_amount_unique_homoplasies(df_syn_clustered_homopl, mut_type="syn")
	.total_amount_unique_homoplasies(df_non_syn_clustered_homopl, mut_type="non-syn")

	system(glue("mkdir -p stat_results/{out_folder}/only_clustered_ordered_homopl/"))
	# SYN
	df_syn_clustered_homopl <- .change_ids(df_syn_clustered_homopl)
	df_syn_clustered_homopl <- df_syn_clustered_homopl[df_syn_clustered_homopl$Freq_homopl >= 5,] # Getting only homoplasies happening at least 5 times
	df_syn_clustered_homopl <- df_syn_clustered_homopl[order(df_syn_clustered_homopl$defining_mut, decreasing=T),]
	write.csv(df_syn_clustered_homopl, glue("stat_results/{out_folder}/only_clustered_ordered_homopl/syn_clustered_homopl_ordered.csv"), quote=F, row.names=F)
	# NON-SYN
	df_non_syn_clustered_homopl <- .change_ids(df_non_syn_clustered_homopl)
	df_non_syn_clustered_homopl <- df_non_syn_clustered_homopl[df_non_syn_clustered_homopl$Freq_homopl >= 5,]
	df_non_syn_clustered_homopl <- df_non_syn_clustered_homopl[order(df_non_syn_clustered_homopl$defining_mut, decreasing=T),]
	write.csv(df_non_syn_clustered_homopl, glue("stat_results/{out_folder}/only_clustered_ordered_homopl/non_syn_clustered_homopl_ordered.csv"), quote=F, row.names=F)

	# Counts number of rows for each index
	df_syn_clustered_homopl_counts <- df_syn_clustered_homopl %>% count(period_threshold)
	df_non_syn_clustered_homopl_counts <- df_non_syn_clustered_homopl %>% count(period_threshold)

	#distribution of freqs per genomic region
	print("Simple distributions per genomic region:")
	df_syn_clustered_homopl_reg_plots_all <- .plot_distr_freqs_spec_regions(df_syn_clustered_homopl,"syn",2,all_clust=TRUE)
	df_non_syn_clustered_homopl_reg_plots_all <- .plot_distr_freqs_spec_regions(df_non_syn_clustered_homopl,"non-syn",2,all_clust=TRUE)

	df_syn_clustered_homopl_reg_plots <- .plot_distr_freqs_spec_regions(df_syn_clustered_homopl,"syn",2,all_clust=FALSE)
	df_non_syn_clustered_homopl_reg_plots <- .plot_distr_freqs_spec_regions(df_non_syn_clustered_homopl,"non-syn",2,all_clust=FALSE)

	# distribution of freqs per major_lineage
	print("Simple distributions per major lineage:")
	df_syn_clustered_homopl_lin_plots <- .plot_distr_freqs_major_lineage(df_syn_clustered_homopl,"syn",2,all_clust=FALSE)
	df_non_syn_clustered_homopl_lin_plots <- .plot_distr_freqs_major_lineage(df_non_syn_clustered_homopl,"non-syn",2,all_clust=FALSE)

	system(glue("mkdir -p stat_results/{out_folder}/higher_eq10perc_syn/"))
	system(glue("mkdir -p stat_results/{out_folder}/higher_eq10perc_non_syn/"))
	system(glue("mkdir -p stat_results/{out_folder}/lower_eq5perc_syn/"))
	system(glue("mkdir -p stat_results/{out_folder}/lower_eq5perc_non_syn/"))

	print("Tables splitting sites resulting from >= 10% and <=5% thresholds:")
	most_common_homopl_syn_period <- most_common_homopl_nonsyn_period <- list()
	for(k in 1:length(major_lineages)) {
		# SYN
		most_common_homopl_syn_period[[k]] <- .get_most_common_homopls_lineage_interest(df_syn_clustered_homopl_reg_plots, period_interest=PERIOD_INTEREST, major_lineages[k])
		write.csv(most_common_homopl_syn_period[[k]][[1]], file=glue("stat_results/{out_folder}/higher_eq10perc_syn/{major_lineages[k]}_top_muts.csv"), quote=F, row.names=F)
		write.csv(most_common_homopl_syn_period[[k]][[2]], file=glue("stat_results/{out_folder}/lower_eq5perc_syn/{major_lineages[k]}_top_muts.csv"), quote=F, row.names=F)
		# NON-SYN
		most_common_homopl_nonsyn_period[[k]] <- .get_most_common_homopls_lineage_interest(df_non_syn_clustered_homopl_reg_plots, period_interest=PERIOD_INTEREST, major_lineages[k])
		write.csv(most_common_homopl_nonsyn_period[[k]][[1]], file=glue("stat_results/{out_folder}/higher_eq10perc_non_syn/{major_lineages[k]}_top_muts.csv"), quote=F, row.names=F)
		write.csv(most_common_homopl_nonsyn_period[[k]][[2]], file=glue("stat_results/{out_folder}/lower_eq5perc_non_syn/{major_lineages[k]}_top_muts.csv"), quote=F, row.names=F)
	}

	# TOP30 by absolute frequency
	print("Top 30 absolute frequency:")
	most_common_homopl_abs_freq_syn <- .get_most_common_homopls_abs_freq(df_syn_clustered_homopl_reg_plots)
	system(glue("mkdir -p stat_results/{out_folder}/abs_freq/"))
	write.csv(most_common_homopl_abs_freq_syn, file=glue("stat_results/{out_folder}/abs_freq/syn_top_muts.csv"), quote=F, row.names=F)
	most_common_homopl_abs_freq_non_syn <- .get_most_common_homopls_abs_freq(df_non_syn_clustered_homopl_reg_plots)
	write.csv(most_common_homopl_abs_freq_non_syn, file=glue("stat_results/{out_folder}/abs_freq/nonsyn_top_muts.csv"), quote=F, row.names=F)
	
	print("Top 30 absolute frequency within Spike:")
	most_common_homopl_s_syn <- .get_most_common_homopls_s_regions(df_syn_clustered_homopl_reg_plots)
	system(glue("mkdir -p stat_results/{out_folder}/abs_freq_spike/"))
	write.csv(most_common_homopl_s_syn, file=glue("stat_results/{out_folder}/abs_freq_spike/syn_top_spike_muts.csv"), quote=F, row.names=F)
	most_common_homopl_s_non_syn <- .get_most_common_homopls_s_regions(df_non_syn_clustered_homopl_reg_plots)
	write.csv(most_common_homopl_s_non_syn, file=glue("stat_results/{out_folder}/abs_freq_spike/nonsyn_top_spike_muts.csv"), quote=F, row.names=F)
	
	# is_clustered = 1 (identified by quantile thresholds of stats)
	most_common_homopl_syn_clust1 <- .get_most_common_homopls(df_syn_clustered_homopl_reg_plots_all,"syn_clust1",is_clustered_flag=TRUE)
	most_common_homopl_non_syn_clust1 <- .get_most_common_homopls(df_non_syn_clustered_homopl_reg_plots_all, "non_syn_clust1",is_clustered_flag=TRUE)
	
	# Bubbles for is_clustered = 1
	print("Bubble plots:")
	bubble_syn_clust1 <- .bubble_plots(most_common_homopl_syn_clust1,"syn","syn_clust1")
	bubble_non_syn_clust1 <- .bubble_plots(most_common_homopl_non_syn_clust1, "non-syn", "non-syn_clust1")
	
	list(bubble_non_syn_clust1, s9, s13)
}

# for union (at least one statistic) of tree-based clustering stats, plot them
plot_tbc_stats <- function(path_stats, out_folder, ymin_vec=c(0,0,-50), ymax_vec=c(10,10,100)) {
	clustered_dfs <- list()
	for(i in 1:length(path_thresholds)) {
		clustered_dfs[[i]] <- read.csv(glue("{path_stats}/{path_thresholds[i]}/stats_union.csv"), header=T)
		clustered_dfs[[i]]$threshold <- table_names_thresholds[i]
		clustered_dfs[[i]]$threshold <- factor(clustered_dfs[[i]]$threshold, levels=table_names_thresholds)
		
	}
	
	clustered_dfs_all_thr <- rbindlist(clustered_dfs)
	clustered_dfs_all_thr <- clustered_dfs_all_thr[clustered_dfs_all_thr$homoplasies == "yes", ]
	
	system(glue("mkdir -p stat_results/{out_folder}/tbc_stats/"))
	
	# RATIO OF SIZES
	rs <- ggplot(clustered_dfs_all_thr, aes(x=threshold, y=ratio_sizes, color=threshold)) + 
		#geom_violin() + geom_sina(size=0.3) +
		geom_boxplot(outlier.size=0.25) +
		scale_color_manual(name = "Threshold", values = pal_thresholds) +
		labs(x="Threshold", y="Ratio of sizes") + coord_cartesian(ylim=c(ymin_vec[1],ymax_vec[1])) +
		theme_bw() + theme(plot.title = element_text(hjust=0.5), axis.text=element_text(size=7), axis.title=element_text(size=7), legend.position="none")
	ggsave(glue("stat_results/{out_folder}/tbc_stats/ratio_sizes.png"), plot=rs, units="px", width=1500, height=1000, dpi=300, bg="white")
	
	# RATIO OF PERSISTENCE TIME
	rpt <- ggplot(clustered_dfs_all_thr, aes(x=threshold, y=ratio_persist_time, color=threshold)) + 
		#geom_violin() + geom_sina(size=0.3) +
		geom_boxplot(outlier.size=0.25) +
		scale_color_manual(name = "Threshold", values = pal_thresholds) +
		labs(x="Threshold", y="Ratio of persistence time") + coord_cartesian(ylim=c(ymin_vec[2],ymax_vec[2])) +
		theme_bw() + theme(plot.title = element_text(hjust=0.5), axis.text=element_text(size=7), axis.title=element_text(size=7), legend.position="none")
	ggsave(glue("stat_results/{out_folder}/tbc_stats/ratio_persist_time.png"), plot=rpt, units="px", width=1500, height=1000, dpi=300, bg="white")
	
	# LOGISTIC GROWTH RATE
	lgr <- ggplot(clustered_dfs_all_thr, aes(x=threshold, y=logistic_growth, color=threshold)) + 
		#geom_violin() + geom_sina(size=0.3) +
		geom_boxplot(outlier.size=0.25) +
		scale_color_manual(name = "Threshold", values = pal_thresholds) +
		labs(x="Threshold", y="Logistic growth rate") + coord_cartesian(ylim=c(ymin_vec[3],ymax_vec[3])) +
		theme_bw() + theme(plot.title = element_text(hjust=0.5), axis.text=element_text(size=7), axis.title=element_text(size=7), legend.position="none")
	ggsave(glue("stat_results/{out_folder}/tbc_stats/logit_growth.png"), plot=lgr, units="px", width=1500, height=1000, dpi=300, bg="white")
	
	rm_x <- theme(axis.title.x=element_blank(), axis.text.x=element_blank())
	
	p_tbc_stats <- ggarrange(rs+rm_x, rpt+rm_x, lgr, nrow=3, ncol=1, labels=c("A","B","C"), common.legend=TRUE, legend="right")
	p_tbc_stats
}

stacked_nsites_genomic_region_mult_thresholds <- function(path) {
	reg_thr_files <- list.files(path = path, pattern = "*REG_THR.csv", full.names = TRUE)
	reg_thr_files_list <- list()
	for(i in 1:length(reg_thr_files)) {
		reg_thr_files_list[[i]] <- read.csv(reg_thr_files[i], header=T)
	}
	
	reg_thr_files_df <- rbindlist(reg_thr_files_list)
	reg_thr_files_df$threshold <- factor(reg_thr_files_df$threshold, levels=table_names_thresholds)
	reg_thr_files_df$spec_regions <- factor(reg_thr_files_df$spec_regions, levels=lvls_nonsyn_proteins_annots)
	#View(reg_thr_files_df)
	View(reg_thr_files_df)
	
	system(glue("mkdir -p {path}/stacked_nsites_threshold/"))
	
	p <- ggplot(reg_thr_files_df, aes(fill=threshold, x=n, y=spec_regions)) + 
		geom_bar(position="stack", stat="identity") + scale_fill_manual(values=pal_thresholds, name="Threshold first\ndetected") + theme(axis.text.y = element_text(size=8,color="black"), axis.text.x=element_text(color="black")) + theme_minimal() +
		labs(y="Genomic region", x="Number of identified sites")
	ggsave(glue("{path}/stacked_nsites_threshold/stacked_nsites_region_thr.pdf"), plot=p, width=8, height=6, dpi=600, bg="white")
	
	p_ia <- ggplotly(p)
	htmlwidgets::saveWidget(p_ia,glue("{path}/stacked_nsites_threshold/stacked_nsites_region_thr.html"))
	
	p
}