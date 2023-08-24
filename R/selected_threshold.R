libs_load <- c("glue","stringr", "dplyr","data.table","ggplot2","tidyr","tidyverse","ggrepel","ggforce", "viridis", "ggpubr", "htmlwidgets", "plotly")
invisible( lapply(libs_load, library, character.only=TRUE) )

# Function to facilitate comparison between periods including and excluding Omicron muts
#' Title
#'
#' @param path_stats 
#' @param rm_freq_outliers 
#' @param thr_index 
#' @param out_folder 
#'
#' @return
#' @export
#'
#' @examples
stats_selected_threshold <- function(path_stats, rm_freq_outliers=TRUE, thr_index, out_folder) {
	
	if(dir.exists(path_stats))
		path <- path_stats
	else
		stop("Directory does not exist!")
	
	.genomewide_plot <- function(df, mut_type) {
		
		if(mut_type == "syn") {
			genomic_ranges_df <- read.csv("config/genomic_ranges_plot_syn.tsv", header=T, sep="\t")
			# No split for SYN
			df_freq_genomic_regions_syn <- df %>% inner_join(genomic_ranges_df, by=c("spec_regions"="new_prot_name"))
			df_freq_genomic_regions_syn_first <- df_freq_genomic_regions_syn %>% distinct(defining_mut, major_lineage, .keep_all=TRUE) %>% group_by(defining_mut) %>% slice_min(n=1, defining_mut) #as.numeric(threshold)
			
			system(glue("mkdir -p stat_results/{out_folder}/genomewide_plot_{mut_type}/"))
			
			# Get table with all TFP-homoplasies
			df_freq_genomic_regions_syn_first_t <- df_freq_genomic_regions_syn_first[order(as.integer(df_freq_genomic_regions_syn_first$Freq_homopl), decreasing=TRUE),]
			
			write.csv(df_freq_genomic_regions_syn_first_t, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_1st_threshold.csv"), quote=F, row.names=F)
			
			# p <- ggplot(df_freq_genomic_regions_syn_first, aes(x=genomic_index, y=Freq_homopl, color=threshold)) +
			# 	geom_point(alpha=0.9) + labs(x="Nucleotide position", y="Frequency in TFP clusters") + theme_minimal() + #theme(axis.text.y=element_text(size=5)) +
			# 	scale_color_manual(values=pal_thresholds, name="Threshold where first detected") +
			# 	coord_cartesian(xlim = c(1, 29903), ylim = c(2,10)) + scale_y_discrete(breaks=seq(2,10,by=1))
			# ggsave(glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_1st_threshold.png"), plot=p, width=15, height=5, dpi=600, bg="white")
			# #p_ia <- ggplotly(p)
			# #htmlwidgets::saveWidget(p_ia, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_1st_threshold.html"))
			
			p2 <- ggplot(df_freq_genomic_regions_syn_first, aes(x=genomic_index, y=Freq_homopl, color=major_lineage)) +
				geom_point(alpha=0.9) + labs(x="Nucleotide position", y="Frequency in TFP clusters") + theme_minimal() + #theme(axis.text.y=element_text(size=5)) +
				scale_color_manual(values=pal_lineages, name="Major lineage") +
				coord_cartesian(xlim = c(0, 29903), ylim = c(2,10)) + scale_y_discrete(breaks=seq(2,10,by=1))
			ggsave(glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_major_lineage.png"), plot=p2, width=15, height=5, dpi=600, bg="white")
			#p2_ia <- ggplotly(p2)
			#htmlwidgets::saveWidget(p2_ia, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_major_lineage.html"))
			
		}else if(mut_type == "non-syn") {
			
			genomic_ranges_df <- read.csv("config/genomic_ranges_plot_non_syn.tsv", header=T, sep="\t")
			df_freq_genomic_regions_nonsyn <- df %>% inner_join(genomic_ranges_df, by=c("spec_regions"="new_prot_name"))
			# Split because each protein (ORF1AB + others) has its own relative coordinates 
			df_freq_genomic_regions_split_nonsyn <- split(df_freq_genomic_regions_nonsyn, df_freq_genomic_regions_nonsyn$protein)
			df_freq_genomic_regions_split_nonsyn_first <- df_freq_genomic_regions_split_nonsyn_first_t <- df_freq_genomic_regions_split_nonsyn_first_t2 <- list()
			for(i in 1:length(df_freq_genomic_regions_split_nonsyn)) {
				df_freq_genomic_regions_split_nonsyn_first[[i]] <- df_freq_genomic_regions_split_nonsyn[[i]] %>% distinct(defining_mut, major_lineage, .keep_all=TRUE) %>% group_by(defining_mut) %>% slice_min(n=1, defining_mut) #, as.numeric(threshold)
				
				# Get tables with all TFP-homoplasies
				df_freq_genomic_regions_split_nonsyn_first_t[[i]] <- df_freq_genomic_regions_split_nonsyn_first[[i]][order(as.integer(df_freq_genomic_regions_split_nonsyn_first[[i]]$Freq_homopl), decreasing=TRUE),]
				
				system(glue("mkdir -p stat_results/{out_folder}/genomewide_plot_{mut_type}/"))
				
				write.csv(df_freq_genomic_regions_split_nonsyn_first_t[[i]], glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_1st_threshold_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.csv"), quote=F, row.names=F)
				
				# p3 <- ggplot(df_freq_genomic_regions_split_nonsyn_first[[i]], aes(x=genomic_index, y=Freq_homopl, color=threshold)) +
				# 	geom_point(alpha=0.9) + labs(x="Aminoacid position", y="Frequency in TFP clusters") + theme_minimal() +
				# 	#facet_grid(~spec_regions, scales = "free", switch = "x") +
				# 	scale_color_manual(values=pal_thresholds, name="Threshold where first detected") +
				# 	geom_text_repel(aes(label = ifelse(Freq_homopl>5,as.character(defining_mut),'')), size=2, box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=15)+
				# 	coord_cartesian(xlim = c(df_freq_genomic_regions_split_nonsyn_first[[i]]$start[1], df_freq_genomic_regions_split_nonsyn_first[[i]]$end[1]))
				# ggsave(glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_1st_threshold_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.png"), plot=p3, width=15, height=5, dpi=600, bg="white")
				# # p3_ia <- ggplotly(p3)
				# # htmlwidgets::saveWidget(p3_ia, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_1st_threshold_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.html"))
				
				p4 <- ggplot(df_freq_genomic_regions_split_nonsyn_first[[i]], aes(x=genomic_index, y=Freq_homopl, color=major_lineage)) +
					geom_point(alpha=0.9) + labs(x="Aminoacid position", y="Frequency in TFP clusters") + theme_minimal() +
					scale_color_manual(values=pal_lineages, name="Major lineage") +
					geom_text_repel(aes(label = ifelse(Freq_homopl>2,as.character(defining_mut),'')), size=2, box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=15)+
					coord_cartesian(xlim = c(df_freq_genomic_regions_split_nonsyn_first[[i]]$start[1], df_freq_genomic_regions_split_nonsyn_first[[i]]$end[1]))
				ggsave(glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_major_lineage_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.png"), plot=p4, width=15, height=5, dpi=600, bg="white")
				# p4_ia <- ggplotly(p4)
				# htmlwidgets::saveWidget(p4_ia, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_major_lineage_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.html"))
				
			}
		} else {
			stop("Wrong choice of mut_type: options are `syn` and `non-syn`")
		}
		
		if(exists("p2")){
			return(p2)
		}
		if(exists("p4")) {
			return(p4)
		}
		
	}
	
	.stacked_bar_freqs <- function(df, mut_type) {
		
		df_freq_genomic_regions_first <- df %>% distinct(defining_mut, major_lineage, .keep_all=TRUE) %>% group_by(defining_mut) %>% slice_min(n=1, defining_mut)
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
		
		p1 <- ggplot(df_freq_genomic_regions_first_counts, aes(fill=major_lineage, y=spec_regions, x=count_uniq_homopls)) + 
			geom_bar(position="stack", stat="identity") + scale_fill_manual(values=pal_lineages, name="Major lineage") + theme_minimal() + theme(axis.text = element_text(size=6,color="black"), axis.title=element_text(size=7,color="black")) +
			labs(x="Count of unique homoplasies", y="Genomic region")
		ggsave(glue("stat_results/{out_folder}/stacked_bar_uniq_counts_{mut_type}/stacked_bar_uniq_counts.png"), width=8, height=6, dpi=600, bg="white")
		
		df_freq_genomic_regions_first_counts_norm <- df_freq_genomic_regions_first %>% group_by(spec_regions, major_lineage) %>% summarise(count_uniq_homopls_norm=sum(Freq_norm))
		
		p2 <- ggplot(df_freq_genomic_regions_first_counts_norm, aes(fill=major_lineage, y=spec_regions, x=count_uniq_homopls_norm)) + 
			geom_bar(position="stack", stat="identity") + scale_fill_manual(values=pal_lineages, name="Major lineage") + theme_minimal() + theme(axis.text = element_text(size=6,color="black"), axis.title=element_text(size=7,color="black")) +
			labs(x="Normalized count of unique homoplasies", y="Genomic region")
		ggsave(glue("stat_results/{out_folder}/stacked_bar_uniq_counts_{mut_type}/stacked_bar_uniq_counts_norm.png"), width=8, height=6, dpi=600, bg="white")
		return(list(p1, p2))
	}
	
	# Top30 and 100
	.get_most_common_homopls_abs_freq <- function(df, thr=THRESHOLD_INTEREST, per=PERIOD_INTEREST) {
		
		df_homopl_binded <- df
		df_homopl_binded$period <- as.integer(df_homopl_binded$period)
		df_homopl_binded$threshold <- as.numeric(df_homopl_binded$threshold)
		df_homopl_binded$Freq_homopl <- as.integer(df_homopl_binded$Freq_homopl)
		
		df_homopl_binded <- df_homopl_binded %>% distinct(defining_mut, major_lineage, .keep_all=TRUE)
		# Top30 for plot
		most_common_homopl_abs_freq <- df_homopl_binded %>% slice_max(n=30, order_by=Freq_homopl) %>% arrange(desc(as.integer(Freq_homopl)))
		# Top100 for table
		most_common_homopl_abs_freq100 <- df_homopl_binded %>% slice_max(n=100, order_by=Freq_homopl) %>% arrange(desc(as.integer(Freq_homopl)))
		
		return(list(most_common_homopl_abs_freq,most_common_homopl_abs_freq100))
	}
	
	# TOP30 only considering S regions
	.get_most_common_homopls_s_regions <- function(df, thr=THRESHOLD_INTEREST, per=PERIOD_INTEREST) {
		
		df_homopl_binded <- df
		df_homopl_binded$period <- as.integer(df_homopl_binded$period)
		df_homopl_binded$threshold <- as.numeric(df_homopl_binded$threshold)
		df_homopl_binded$Freq_homopl <- as.integer(df_homopl_binded$Freq_homopl)
		
		most_common_homopl_s_regions <- df_homopl_binded %>% distinct(defining_mut, major_lineage, .keep_all=TRUE)
		most_common_homopl_s_regions <- most_common_homopl_s_regions %>% filter(protein=="S") %>% slice_max(n=30, order_by=Freq_homopl)
		
		return(most_common_homopl_s_regions)
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
		
		most_common_homopl_all <- most_common_homopl_df #most_common_homopl_df[[2]]
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
			#joined_most_common_indep_ok$protein <- factor(joined_most_common_indep_ok$protein, levels=joined_most_common_indep_ok$protein[ order( unique(as.numeric(joined_most_common_indep_ok$site)) ) ])
			joined_most_common_indep_ok <- joined_most_common_indep_ok[order(joined_most_common_indep_ok$protein, as.numeric(joined_most_common_indep_ok$site) ),]
			joined_most_common_indep_ok$prot_site <- factor(joined_most_common_indep_ok$prot_site, levels=unique( joined_most_common_indep_ok$prot_site[order(joined_most_common_indep_ok$protein, as.numeric(joined_most_common_indep_ok$site) )] ))
			#View(joined_most_common_indep_ok)
		}
		
		system(glue("mkdir -p stat_results/{out_folder}/comparison_bubble_{out_suffix}/"))
		
		pal_id <- c("mlscluster"="#7c9885","Intersection"="#e1ad01", "HyPhy"="#033f63")
		
		p <- ggplot(joined_most_common_indep_ok, aes(y=analysis_identified.y, x=prot_site, color=analysis_identified.y)) + #color=analysis_identified.y
			geom_point(alpha=0.8) + labs(y="Analyses identified",x="Protein site") + theme_minimal() +
			scale_color_manual(values=pal_id, name="Analysis identified") + theme(legend.position="none", axis.text.x=element_text(size=5, angle=90, vjust = 0.5, hjust=1,color="black"), axis.text.y=element_text(size=7,color="black"), axis.title=element_text(size=9)) + marg
		ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/homopl_method_combinations.png"), plot=p, width=12, height=5, dpi=600, bg="white")
		#p_ia <- ggplotly(p)
		#htmlwidgets::saveWidget(p_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/homopl_method_combinations.html"))
		
		joined_most_common_indep_ok_high_freq_lin <- joined_most_common_indep_ok %>% group_by(major_lineage, prot_site) %>% slice_max(Freq_homopl)
		
		# Plot max frequency of homoplasic sites (bubble sizes) for each one of the major_lineage (coloured by analysis where identified)
		joined_most_common_indep_ok_high_freq_lin <- joined_most_common_indep_ok_high_freq_lin[joined_most_common_indep_ok_high_freq_lin$analysis_identified.y != "HyPhy",]
		joined_most_common_indep_ok_high_freq_lin <- joined_most_common_indep_ok_high_freq_lin[!is.na(joined_most_common_indep_ok_high_freq_lin$major_lineage),]
		
		p2 <- ggplot(joined_most_common_indep_ok_high_freq_lin, aes(x=major_lineage, y=prot_site, size=Freq_homopl, color=analysis_identified.y)) +
			geom_point(alpha=0.5) + labs(x="Major lineage", y="Protein site") + theme_minimal() + theme(axis.text.y = element_text(size=5,color="black"), axis.text.x = element_text(size=7, color="black", angle = 20,hjust = 1), axis.title=element_text(size=9)) + marg +
			scale_size_area(name="TFP-homoplasy\nfrequency") + scale_color_manual(values=pal_id, name="Analysis identified")
		ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_all_lineages.png"), plot=p2, width=12, height=15, dpi=600, bg="white")
		#p2_ia <- ggplotly(p2)
		#htmlwidgets::saveWidget(p2_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_all_lineages.html"))
		
		f4 <- ggarrange(p, p2, nrow=2, ncol=1, labels=c('A', 'B'), heights=c(1,2))
		ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/compHyPhy_lineages.pdf"), width=9, height=8, dpi=600, bg="white")
		
		most_common_homopl_all_high_freq_lin <- most_common_homopl_all %>% group_by(major_lineage, defining_mut) %>% slice_max(Freq_homopl)
		most_common_homopl_all_high_freq_lin <- most_common_homopl_all_high_freq_lin[!is.na(most_common_homopl_all_high_freq_lin$major_lineage),]
		# Plot max frequency of each homoplasy across thresholds (including WT and replaced aa) identified (without comparison with HyPhy)
		p3 <- ggplot(most_common_homopl_all_high_freq_lin, aes(x=major_lineage, y=defining_mut, size=Freq_homopl)) + #color=analysis_identified.y
			geom_point(alpha=0.5) + labs(x="Major lineage", y="Homoplasy") + theme_minimal() + theme(axis.text.y = element_text(size=5), aspect.ratio=16/9) +
			scale_size_area(name="Highest freq of homopl") # breaks=c(20,50,150,500,1300) + scale_color_manual(values=c("#7c9885","#e1ad01", "#033f63")) +
		ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_replacement_lineages.png"), plot=p3, width=12, height=15, dpi=600, bg="white")
		#p3_ia <- ggplotly(p3)
		#htmlwidgets::saveWidget(p3_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_replacement_lineages.html"))
		
		return(f4)
	}
	
	system(glue("mkdir -p stat_results/{out_folder}/"))
	
	joined_mut_sites_clustered <- joined_mut_sites_clustered_syn <- joined_mut_sites_clustered_non_syn <- list()
	csv_cp <- thr2 <- list()
	
	clustered_df <- read.csv(glue("{path_stats}/{path_thresholds[1]}/clustered_all_df.csv"), header=T)
	# If flag to remove 99% quantile frequency outliers is TRUE, then remove
	if(rm_freq_outliers) {
		clustered_df <- remove_homopl_freq_outliers(path_stats)
	}
	clustered_df <- clustered_df[[thr_index]]
	
	clustered_df <- clustered_df[clustered_df$is_clustered==1,]
	
	clustered_df$major_lineage <- ifelse(grepl("Beta_B.1.351", clustered_df$major_lineage), "Other", as.character(clustered_df$major_lineage))
	clustered_df$major_lineage <- ifelse(grepl("Gamma_P.1", clustered_df$major_lineage), "Other", as.character(clustered_df$major_lineage))
	
	joined_mut_sites_clustered <- clustered_df
	
	# Extract genome index position
	rgx_genomic_idx <- regexpr(':[A-Z]{1}[0-9]{1,5}[A-Z*]{1}', joined_mut_sites_clustered$defining_mut)
	comm_coord <- rep(NA,length(joined_mut_sites_clustered$defining_mut))
	comm_coord[rgx_genomic_idx != -1] <- str_sub( regmatches(joined_mut_sites_clustered$defining_mut, rgx_genomic_idx), 3, -2)
	comm_coord <- as.integer(comm_coord)
	comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("genomic_index")
	joined_mut_sites_clustered <- cbind(joined_mut_sites_clustered, comm_coord_df)
	
	# Adjust df columns
	joined_mut_sites_clustered <- joined_mut_sites_clustered %>% dplyr::select(defining_mut, Freq_homopl, is_clustered, major_lineage, syn_non_syn, protein, aa_length, indep_found_pos_selection, genomic_index) #s_mut_region_interest
	colnames(joined_mut_sites_clustered) <- c("defining_mut","Freq_homopl","is_clustered","major_lineage","syn","protein","aa_length","indep_found_pos_selection","genomic_index") #s_mut_region_interest
	
	joined_mut_sites_clustered$protein <- sub("\\:.*", "", joined_mut_sites_clustered$defining_mut)
	
	setDT(joined_mut_sites_clustered); setDT(annot_gen_range_positions)
	joined_mut_sites_clustered <- joined_mut_sites_clustered[annot_gen_range_positions, on=c("protein==name_cog_md","genomic_index>=start","genomic_index<=end"), spec_regions := new_prot_name]
	joined_mut_sites_clustered <- joined_mut_sites_clustered[annot_gen_length_positions, on=c("spec_regions==new_prot_name"), aa_length := length]
	
	# syn flag
	joined_mut_sites_clustered$syn <- ifelse(joined_mut_sites_clustered$syn == as.character("non-syn"), 0, 1)
	# independent found under selection flag
	joined_mut_sites_clustered$indep_found_pos_selection <- ifelse(joined_mut_sites_clustered$indep_found_pos_selection == as.character("yes"), 1, 0)
	# Add zeros to Freq_homopl and is_clustered if polymorphic site==NA (not found on homoplasy tables)
	joined_mut_sites_clustered <- joined_mut_sites_clustered %>% mutate(Freq_homopl = ifelse(is.na(Freq_homopl), 0, Freq_homopl))
	joined_mut_sites_clustered <- joined_mut_sites_clustered %>% mutate(is_clustered = ifelse(is.na(is_clustered), 0, is_clustered))
	
	# Not sure if should use this since flag is already set above
	rgx_scpss <- regexpr(':[A-Z]{1}[0-9]{1,5}[A-Z*]{1}', joined_mut_sites_clustered$defining_mut)
	comm_coord <- rep(NA,length(joined_mut_sites_clustered$defining_mut))
	# Mutation coordinate of homoplasies
	comm_coord[rgx_scpss != -1] <- str_sub( regmatches(joined_mut_sites_clustered$defining_mut, rgx_scpss), 3, -2)
	comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("site")
	joined_mut_sites_clustered <- cbind(joined_mut_sites_clustered, comm_coord_df)
	joined_mut_sites_clustered$prot_site <- paste0(joined_mut_sites_clustered$protein,":",joined_mut_sites_clustered$site)
	
	joined_mut_sites_clustered$Freq_homopl <- as.integer(joined_mut_sites_clustered$Freq_homopl)
	
	joined_mut_sites_clustered$is_clustered_h3 <- ifelse( (as.integer(joined_mut_sites_clustered$Freq_homopl) > 3) & (joined_mut_sites_clustered$is_clustered==1), 1, 0)
	# Normalized counts (for inspection only)
	joined_mut_sites_clustered$Freq_homopl_norm_size <- as.integer(joined_mut_sites_clustered$Freq_homopl) / as.integer(joined_mut_sites_clustered$aa_length)
	
	# Make sure to remove all rows with NAs in key columns
	joined_mut_sites_clustered <- joined_mut_sites_clustered[!is.na(joined_mut_sites_clustered$Freq_homopl) | !is.na(joined_mut_sites_clustered$major_lineage) | !is.na(joined_mut_sites_clustered$is_clustered) | !is.na(joined_mut_sites_clustered$is_clustered_h3) | !is.na(joined_mut_sites_clustered$spec_regions) | !is.na(joined_mut_sites_clustered$aa_length) | !is.na(joined_mut_sites_clustered$indep_found_pos_selection), ]
	
	joined_mut_sites_clustered$site <- as.numeric(joined_mut_sites_clustered$site)
	
	# Split into syn and non-syn
	system(glue("mkdir -p stat_results/{out_folder}/stat_ready_csvs_all_sites/"))
	system(glue("mkdir -p stat_results/{out_folder}/stat_ready_csvs_only_clustered/"))
	
	# SYN
	joined_mut_sites_clustered_syn <- joined_mut_sites_clustered[joined_mut_sites_clustered$syn == 1,]
	joined_mut_sites_clustered_syn$major_lineage <- as.factor( joined_mut_sites_clustered_syn$major_lineage )
	joined_mut_sites_clustered_syn$spec_regions <- as.factor(joined_mut_sites_clustered_syn$spec_regions)
	csv_cp <- joined_mut_sites_clustered_syn[ order(joined_mut_sites_clustered_syn$Freq_homopl, decreasing=T) ,]
	write.csv(csv_cp, glue("stat_results/{out_folder}/stat_ready_csvs_all_sites/df_all_sites_{path_thresholds[thr_index]}_SYN.csv"), quote=F, row.names=F)
	csv_cp <- csv_cp[(as.integer(csv_cp$Freq_homopl) > 0) & (as.integer(csv_cp$is_clustered)==1), ]
	write.csv(csv_cp, glue("stat_results/{out_folder}/stat_ready_csvs_only_clustered/df_only_clustered_{path_thresholds[thr_index]}_SYN.csv"), quote=F, row.names=F)
	
	# NON-SYN
	joined_mut_sites_clustered_non_syn <- joined_mut_sites_clustered[joined_mut_sites_clustered$syn == 0,]
	joined_mut_sites_clustered_non_syn$major_lineage <- as.factor( joined_mut_sites_clustered_non_syn$major_lineage )
	joined_mut_sites_clustered_non_syn$spec_regions <- as.factor(joined_mut_sites_clustered_non_syn$spec_regions)
	csv_cp2 <- joined_mut_sites_clustered_non_syn[ order(joined_mut_sites_clustered_non_syn$Freq_homopl, decreasing=T) ,]
	write.csv(csv_cp2, glue("stat_results/{out_folder}/stat_ready_csvs_all_sites/df_all_sites_{path_thresholds[thr_index]}_NON_SYN.csv"), quote=F, row.names=F)
	csv_cp2 <- csv_cp2[(as.integer(csv_cp2$Freq_homopl) > 0) & (as.integer(csv_cp2$is_clustered)==1), ]
	write.csv(csv_cp2, glue("stat_results/{out_folder}/stat_ready_csvs_only_clustered/df_only_clustered_{path_thresholds[thr_index]}_NON_SYN.csv"), quote=F, row.names=F)
	
	df_syn_clustered_homopl <- joined_mut_sites_clustered_syn
	df_non_syn_clustered_homopl <- joined_mut_sites_clustered_non_syn
	
	df_syn_clustered_homopl$period <- PERIOD_INTEREST
	df_non_syn_clustered_homopl$period <- PERIOD_INTEREST
	df_syn_clustered_homopl$threshold <- THRESHOLD_INTEREST
	df_non_syn_clustered_homopl$threshold <- THRESHOLD_INTEREST
	
	print("Genome-wide plots per genomic region:")
	.genomewide_plot(df_syn_clustered_homopl, mut_type="syn")
	gp_ns <- .genomewide_plot(df_non_syn_clustered_homopl, mut_type="non-syn")

	print("Stacked bar of unique sites under selection by lineage and genomic region:")
	sbf_s <- .stacked_bar_freqs(df_syn_clustered_homopl, mut_type="syn")
	sbf_ns <- .stacked_bar_freqs(df_non_syn_clustered_homopl, mut_type="non-syn")

	# Top30 and 100 by absolute frequency
	print("Top30 and 100 regardless of genomic region:")
	most_common_homopl_abs_freq_syn <- .get_most_common_homopls_abs_freq(df_syn_clustered_homopl)
	system(glue("mkdir -p stat_results/{out_folder}/abs_freq/"))
	write.csv(most_common_homopl_abs_freq_syn[[1]], file=glue("stat_results/{out_folder}/abs_freq/syn_top30_muts.csv"), quote=F, row.names=F)
	write.csv(most_common_homopl_abs_freq_syn[[2]], file=glue("stat_results/{out_folder}/abs_freq/syn_top100_muts.csv"), quote=F, row.names=F)
	most_common_homopl_abs_freq_non_syn <- .get_most_common_homopls_abs_freq(df_non_syn_clustered_homopl)
	write.csv(most_common_homopl_abs_freq_non_syn[[1]], file=glue("stat_results/{out_folder}/abs_freq/nonsyn_top30_muts.csv"), quote=F, row.names=F)
	write.csv(most_common_homopl_abs_freq_non_syn[[2]], file=glue("stat_results/{out_folder}/abs_freq/nonsyn_top100_muts.csv"), quote=F, row.names=F)

	# Top30 spike
	print("Top30 on Spike:")
	most_common_homopl_s_non_syn <- .get_most_common_homopls_s_regions(df_non_syn_clustered_homopl)
	system(glue("mkdir -p stat_results/{out_folder}/abs_freq_s/"))
	write.csv(most_common_homopl_s_non_syn, file=glue("stat_results/{out_folder}/abs_freq_s/nonsyn_top_muts_spike.csv"), quote=F, row.names=F)

	# Bubbles for is_clustered = 1
	print("Bubble plots:")
	bubble_syn_clust1 <- .bubble_plots(most_common_homopl_abs_freq_syn[[1]],"syn","syn_clust1")
	bubble_non_syn_clust1 <- .bubble_plots(most_common_homopl_abs_freq_non_syn[[1]], "non-syn", "non-syn_clust1")
	
	return(list(sbf_s, sbf_ns, gp_ns, bubble_non_syn_clust1))
}

