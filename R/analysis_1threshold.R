library(glue)
library(stringr)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggrepel)
library(ggforce)
library(viridis)
library(ggpubr)
library(htmlwidgets) #install.packages("htmlwidgets")
library(plotly)

#load(file="rds/env.RData")

# Palette for thresholds (OLD: f0f8ff)
pal_thresholds <- c("#1984c5", "#22a7f0", "#63bff0", "#a7d5ed", "#cbfeff", "#ffcccb", "#e1a692", "#de6e56", "#e14b31", "#c23728") # removed e2e2e2

# Palette for major_lineages
pal_lineages <- c("Alpha_B.1.1.7"="#fd7f6f", "Delta_AY.4.*"="#7eb0d5", "Delta_other"="#b2e061", "EU1_B.1.177"="#bd7ebe", "Omicron_BA.1.*"="#ffb55a", "Omicron_BA.2.*"="#ffee65", "Other"="#fdcce5") #removed "#beb9db", "#8bd3c7"

syn_p <- "SYNSNP:"
lvls_syn_proteins <- c("NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NPS12","NSP13","NSP14","NSP15",
																							"NSP16","S","ORF3A","E","M","ORF6","ORF7A","ORF7B","ORF8","N","ORF10")

lvls_nonsyn_proteins_annots <- c("NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NSP12","NSP13","NSP14","NSP15","NSP16",
																																	"S","S:NTD","S:RBD","S:FCS","ORF3A","E","M","ORF6","ORF7A","ORF7B","ORF8","N","N:Linker_domain","ORF10")

lvls_nonsyn_proteins <- c("NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NSP12","NSP13","NSP14","NSP15","NSP16",
																										"S","ORF3A","E","M","ORF6","ORF7A","ORF7B","ORF8","N","ORF10")

lvls_bubble <- c("ORF1AB","S","ORF3A","E","M","ORF6","ORF7A","ORF7B","ORF8","N","ORF10")

# Load config csv's
seq_muts_df_counts <- readRDS("rds/polymorphic_site_counts.rds")
pos_sel_sites <- read.csv("config/positive_selected_sites_20221006.tsv", sep="\t", header=TRUE)
pos_sel_sites$prot_site <- paste0(pos_sel_sites$protein,":",pos_sel_sites$site)
annot_gen_range_positions <- read.csv("config/aa_ranges_stat_model.tsv",sep="\t")
annot_gen_length_positions <- read.csv("config/aa_length_stat_model.tsv",sep="\t")

# # Change period_threshold to match correct period (1 to 5) and threshold (0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 10, 25)
# change_ids <- function(df) {
# 	for(i in 1:nrow(df)) {
# 		if(df[[i,1]] %in% names(table_combs)) {
# 			idx <- match(df[[i,1]],names(table_combs))
# 			df[[i,1]] <- table_combs[idx]
# 		}
# 	}
# 	return(df)
# }

genomewide_plot <- function(df, mut_type) {
	# TODO? Annotate ORF1AB NSPs, S regions and N linker
	#df <- df %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
	#df$threshold <- factor(df$threshold, levels=table_names_thresholds)
	
	if(mut_type == "syn") {
		genomic_ranges_df <- read.csv("config/genomic_ranges_plot_syn.tsv", header=T, sep="\t")
		# No split for SYN
		df_freq_genomic_regions_syn <- df %>% inner_join(genomic_ranges_df, by=c("spec_regions"="new_prot_name"))
		#df_freq_genomic_regions_syn <- factor(df_freq_genomic_regions_syn$spec_regions, levels=lvls_syn_proteins)
		df_freq_genomic_regions_syn_first <- df_freq_genomic_regions_syn %>% distinct(defining_mut, major_lineage, .keep_all=TRUE) %>% group_by(defining_mut) %>% slice_min(n=1, defining_mut) #as.numeric(threshold)
		#df_freq_genomic_regions_syn_first <- df_freq_genomic_regions_syn_first[!is.na(df_freq_genomic_regions_syn_first$threshold),]

		system(glue("mkdir -p stat_results/{out_folder}/genomewide_plot_{mut_type}/"))
		
		# Get table with all TFP-homoplasies
		df_freq_genomic_regions_syn_first_t <- df_freq_genomic_regions_syn_first[order(as.numeric(df_freq_genomic_regions_syn_first$Freq_homopl), decreasing=TRUE),]
		# Get table grouped by spec_regions and threshold
		#df_freq_genomic_regions_syn_first_t2 <- df_freq_genomic_regions_syn_first_t %>% group_by(spec_regions, threshold) %>% summarise(n=n())
		
		write.csv(df_freq_genomic_regions_syn_first_t, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_1st_threshold.csv"), quote=F, row.names=F)
		#write.csv(df_freq_genomic_regions_syn_first_t2, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_1st_threshold_REG_THR.csv"), quote=F, row.names=F)
		
		# p <- ggplot(df_freq_genomic_regions_syn_first, aes(x=genomic_index, y=Freq_homopl, color=threshold)) +
		# 	geom_point(alpha=0.9) + labs(x="Nucleotide position", y="Frequency in TFP clusters") + theme_minimal() + #theme(axis.text.y=element_text(size=5)) +
		# 	scale_color_manual(values=pal_thresholds, name="Threshold where first detected") +
		# 	coord_cartesian(xlim = c(1, 29903), ylim = c(2,10)) + scale_y_continuous(breaks=seq(2,10,by=1))
		# ggsave(glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_1st_threshold.png"), plot=p, width=15, height=5, dpi=600, bg="white")
		# p_ia <- ggplotly(p)
		# htmlwidgets::saveWidget(p_ia, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_1st_threshold.html"))
		
		p2 <- ggplot(df_freq_genomic_regions_syn_first, aes(x=genomic_index, y=Freq_homopl, color=major_lineage)) +
			geom_point(alpha=0.9) + labs(x="Nucleotide position", y="Frequency in TFP clusters") + theme_minimal() + #theme(axis.text.y=element_text(size=5)) +
			scale_color_manual(values=pal_lineages, name="Major lineage") +
			coord_cartesian(xlim = c(0, 29903), ylim = c(2,10)) + scale_y_continuous(breaks=seq(2,10,by=1))
		ggsave(glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_major_lineage.png"), plot=p2, width=15, height=5, dpi=600, bg="white")
		p2_ia <- ggplotly(p2)
		htmlwidgets::saveWidget(p2_ia, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/syn_freqs_major_lineage.html"))
		
	}else if(mut_type == "non-syn") {
		genomic_ranges_df <- read.csv("config/genomic_ranges_plot_non_syn.tsv", header=T, sep="\t")
		#print(unique(df_freq_genomic_regions$protein))
		df_freq_genomic_regions_nonsyn <- df %>% inner_join(genomic_ranges_df, by=c("spec_regions"="new_prot_name"))
		#df_freq_genomic_regions_nonsyn <- factor(df_freq_genomic_regions_nonsyn$spec_regions, levels=lvls_bubble) #lvls_nonsyn_proteins_annots
		# Split because each protein (ORF1AB + others) has its own relative coordinates 
		df_freq_genomic_regions_split_nonsyn <- split(df_freq_genomic_regions_nonsyn, df_freq_genomic_regions_nonsyn$protein)
		df_freq_genomic_regions_split_nonsyn_first <- df_freq_genomic_regions_split_nonsyn_first_t <- df_freq_genomic_regions_split_nonsyn_first_t2 <- list()
		for(i in 1:length(df_freq_genomic_regions_split_nonsyn)) {
			df_freq_genomic_regions_split_nonsyn_first[[i]] <- df_freq_genomic_regions_split_nonsyn[[i]] %>% distinct(defining_mut, major_lineage, .keep_all=TRUE) %>% group_by(defining_mut) %>% slice_min(n=1, defining_mut) #, as.numeric(threshold)
			#df_freq_genomic_regions_split_nonsyn_first[[i]] <- df_freq_genomic_regions_split_nonsyn_first[[i]][!is.na(df_freq_genomic_regions_split_nonsyn_first[[i]]$threshold),]
			
			# DROP OUTLIERS WITHOUT MANUAL CURATION DONE YET
			#df_freq_genomic_regions_split_nonsyn_first[[i]] <- df_freq_genomic_regions_split_nonsyn_first[[i]] %>% filter(!(Freq_homopl > 20 & defining_mut %in% OUTLIERS_DROP))
			
			# Get tables with all TFP-homoplasies
			df_freq_genomic_regions_split_nonsyn_first_t[[i]] <- df_freq_genomic_regions_split_nonsyn_first[[i]][order(as.numeric(df_freq_genomic_regions_split_nonsyn_first[[i]]$Freq_homopl), decreasing=TRUE),]
			# Get tables grouped by spec_regions and threshold
			#df_freq_genomic_regions_split_nonsyn_first_t2[[i]] <- df_freq_genomic_regions_split_nonsyn_first_t[[i]] %>% group_by(spec_regions, threshold) %>% summarise(n=n())
			
			system(glue("mkdir -p stat_results/{out_folder}/genomewide_plot_{mut_type}/"))
			
			write.csv(df_freq_genomic_regions_split_nonsyn_first_t[[i]], glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_1st_threshold_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.csv"), quote=F, row.names=F)
			#write.csv(df_freq_genomic_regions_split_nonsyn_first_t2[[i]], glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_1st_threshold_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}_REG_THR.csv"), quote=F, row.names=F)
			
			# p3 <- ggplot(df_freq_genomic_regions_split_nonsyn_first[[i]], aes(x=genomic_index, y=Freq_homopl, color=threshold)) +
			# 	geom_point(alpha=0.9) + labs(x="Aminoacid position", y="Frequency in TFP clusters") + theme_minimal() +
			# 	#facet_grid(~spec_regions, scales = "free", switch = "x") +
			# 	scale_color_manual(values=pal_thresholds, name="Threshold where first detected") +
			# 	geom_text_repel(aes(label = ifelse(Freq_homopl>5,as.character(defining_mut),'')), size=2, box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=15)+
			# 	coord_cartesian(xlim = c(df_freq_genomic_regions_split_nonsyn_first[[i]]$start[1], df_freq_genomic_regions_split_nonsyn_first[[i]]$end[1])) #, ylim = c(2,10)) + scale_y_continuous(breaks=seq(2,10,by=1))
			# ggsave(glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_1st_threshold_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.png"), plot=p3, width=15, height=5, dpi=600, bg="white")
			# p3_ia <- ggplotly(p3)
			# htmlwidgets::saveWidget(p3_ia, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_1st_threshold_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.html"))
			
			p4 <- ggplot(df_freq_genomic_regions_split_nonsyn_first[[i]], aes(x=genomic_index, y=Freq_homopl, color=major_lineage)) +
				geom_point(alpha=0.9) + labs(x="Aminoacid position", y="Frequency in TFP clusters") + theme_minimal() +
				#facet_grid(~spec_regions, scales = "free", switch = "x") +
				scale_color_manual(values=pal_lineages, name="Major lineage") +
				#geom_tile(aes(x=spec_regions, y=1, fill=)) +
				geom_text_repel(aes(label = ifelse(Freq_homopl>2,as.character(defining_mut),'')), size=2, box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps=15)+
				coord_cartesian(xlim = c(df_freq_genomic_regions_split_nonsyn_first[[i]]$start[1], df_freq_genomic_regions_split_nonsyn_first[[i]]$end[1])) #, ylim = c(2,10)) + scale_y_continuous(breaks=seq(2,10,by=1))
			ggsave(glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_major_lineage_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.png"), plot=p4, width=15, height=5, dpi=600, bg="white")
			p4_ia <- ggplotly(p4)
			htmlwidgets::saveWidget(p4_ia, glue("stat_results/{out_folder}/genomewide_plot_{mut_type}/nonsyn_freqs_major_lineage_{unique(df_freq_genomic_regions_split_nonsyn_first[[i]]$protein)}.html"))
			
		}
		#View(df_freq_genomic_regions_split_nonsyn_first[[1]])
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

violin_boxplot_freq_variations <- function(df, mut_type) {
	#df <- df %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
	#df$threshold <- factor(df$threshold, levels=table_names_thresholds)
	
	df_freq_genomic_regions_first <- df %>% distinct(defining_mut, major_lineage, .keep_all=TRUE) %>% group_by(defining_mut) %>% slice_min(n=1, defining_mut)
	
	if(mut_type == "syn") {
		df_freq_genomic_regions_first$spec_regions <- sub("SYNSNP:","",df_freq_genomic_regions_first$spec_regions)
		#df_freq_genomic_regions_first[] <- lapply(df_freq_genomic_regions_first, function(x) str_remove(x, "SYNSNP:"))
		#View(df_freq_genomic_regions_first)
		df_freq_genomic_regions_first$spec_regions <- factor(df_freq_genomic_regions_first$spec_regions, levels=lvls_syn_proteins)
	} else {
		df_freq_genomic_regions_first$spec_regions <- factor(df_freq_genomic_regions_first$spec_regions, levels=lvls_nonsyn_proteins_annots)
		#df_freq_genomic_regions_first <- df_freq_genomic_regions_first %>% filter(!(Freq_homopl > 20 & defining_mut %in% OUTLIERS_DROP))
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
	#df <- df %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
	#df$threshold <- factor(df$threshold, levels=table_names_thresholds)
	
	df_freq_genomic_regions_first <- df %>% distinct(defining_mut, major_lineage, .keep_all=TRUE) %>% group_by(defining_mut) %>% slice_min(n=1, defining_mut)
	df_freq_genomic_regions_first <- df_freq_genomic_regions_first[!is.na(df_freq_genomic_regions_first$spec_regions),]
	
	df_freq_genomic_regions_first$Freq_norm <- df_freq_genomic_regions_first$Freq_homopl / df_freq_genomic_regions_first$aa_length
	
	if(mut_type == "syn") {
		df_freq_genomic_regions_first$spec_regions <- sub(".*:","",df_freq_genomic_regions_first$spec_regions)
		#df_freq_genomic_regions_first[] <- lapply(df_freq_genomic_regions_first, function(x) str_remove(x, "SYNSNP:"))
		#View(df_freq_genomic_regions_first)
		df_freq_genomic_regions_first$spec_regions <- factor(df_freq_genomic_regions_first$spec_regions, levels=lvls_syn_proteins)
	} else {
		df_freq_genomic_regions_first$spec_regions <- factor(df_freq_genomic_regions_first$spec_regions, levels=lvls_nonsyn_proteins_annots)
		#df_freq_genomic_regions_first <- df_freq_genomic_regions_first %>% filter(!(Freq_homopl > 20 & defining_mut %in% OUTLIERS_DROP))
	}
	
	df_freq_genomic_regions_first_counts <- df_freq_genomic_regions_first %>% group_by(spec_regions, major_lineage) %>% summarise(count_uniq_homopls=n())
	
	system(glue("mkdir -p stat_results/{out_folder}/stacked_bar_uniq_counts_{mut_type}/"))
	
	p1 <- ggplot(df_freq_genomic_regions_first_counts, aes(fill=major_lineage, y=spec_regions, x=count_uniq_homopls)) + 
		geom_bar(position="stack", stat="identity") + scale_fill_manual(values=pal_lineages, name="Major lineage") + theme_minimal() + theme(axis.text = element_text(size=6,color="black"), axis.title=element_text(size=7,color="black")) +
		labs(x="Count of unique homoplasies", y="Genomic region")
	#ggsave(glue("stat_results/{out_folder}/stacked_bar_uniq_counts_{mut_type}/stacked_bar_uniq_counts.png"), width=8, height=6, dpi=600, bg="white")
	
	df_freq_genomic_regions_first_counts_norm <- df_freq_genomic_regions_first %>% group_by(spec_regions, major_lineage) %>% summarise(count_uniq_homopls_norm=sum(Freq_norm))
	
	p2 <- ggplot(df_freq_genomic_regions_first_counts_norm, aes(fill=major_lineage, y=spec_regions, x=count_uniq_homopls_norm)) + 
		geom_bar(position="stack", stat="identity") + scale_fill_manual(values=pal_lineages, name="Major lineage") + theme_minimal() + theme(axis.text = element_text(size=6,color="black"), axis.title=element_text(size=7,color="black")) +
		labs(x="Normalised count of unique homoplasies", y="Genomic region")
	#ggsave(glue("stat_results/{out_folder}/stacked_bar_uniq_counts_{mut_type}/stacked_bar_uniq_counts_norm.png"), width=8, height=6, dpi=600, bg="white")
	return(list(p1, p2))
}

total_amount_unique_homoplasies <- function(df, mut_type) {
	#df <- df %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
	#df$threshold <- as.numeric(df$threshold)
	#df$threshold <- factor(df$threshold, levels=table_names_thresholds)
	
	df_freq_genomic_regions_first <- df %>% distinct(defining_mut, major_lineage, .keep_all=TRUE) %>% group_by(defining_mut) %>% slice_min(n=1, defining_mut)
	
	if(mut_type == "syn") {
		df_freq_genomic_regions_first$spec_regions <- sub(".*:","",df_freq_genomic_regions_first$spec_regions)
	}
	
	#df_freq_genomic_regions_first <- df_freq_genomic_regions_first %>% filter(!(Freq_homopl > 20 & defining_mut %in% OUTLIERS_DROP))
	
	system(glue("mkdir -p stat_results/{out_folder}/total_homopls_{mut_type}/"))
	
	total_homopls <- nrow(df_freq_genomic_regions_first)
	print(glue("Total number of unique homoplasies for mut_type={mut_type} = {total_homopls}"))
	
	total_homopls_region <- df_freq_genomic_regions_first %>% group_by(spec_regions) %>% summarise(total=n())
	write.csv(total_homopls_region, glue("stat_results/{out_folder}/total_homopls_{mut_type}/total_homopls_region.csv"), quote=F, row.names=F)
	
	total_homopls_region_lineage <- df_freq_genomic_regions_first %>% group_by(spec_regions, major_lineage) %>% summarise(total=n())
	write.csv(total_homopls_region_lineage, glue("stat_results/{out_folder}/total_homopls_{mut_type}/total_homopls_region_lineage.csv"), quote=F, row.names=F)
	
	# total_homopls_region_lower_thr <- df_freq_genomic_regions_first %>% filter(threshold<=1) %>% group_by(spec_regions) %>% summarise(total=n())
	# write.csv(total_homopls_region_lower_thr, glue("stat_results/{out_folder}/total_homopls_{mut_type}/total_homopls_region_lower_thr.csv"), quote=F, row.names=F)
	# 
	# total_homopls_region_higher_thr <- df_freq_genomic_regions_first %>% filter(threshold>=2) %>% group_by(spec_regions) %>% summarise(total=n())
	# write.csv(total_homopls_region_higher_thr, glue("stat_results/{out_folder}/total_homopls_{mut_type}/total_homopls_region_higher_thr.csv"), quote=F, row.names=F)
}


# TOP30
get_most_common_homopls_abs_freq <- function(df, thr=THRESHOLD_INTEREST, per=PERIOD_INTEREST) {
	#df_homopl_binded <- rbindlist(df)
	
	#df_homopl_binded <- df_homopl_binded[df_homopl_binded$is_clustered==1,]
	
	#df_homopl_binded <- df_homopl_binded %>% separate(period_threshold, c("period","threshold"), sep="_")
	df_homopl_binded <- df
	df_homopl_binded$period <- as.integer(df_homopl_binded$period)
	df_homopl_binded$threshold <- as.numeric(df_homopl_binded$threshold)
	df_homopl_binded$Freq_homopl <- as.integer(df_homopl_binded$Freq_homopl)
	
	#most_common_homopl_abs_freq <- df_homopl_binded[(df_homopl_binded$threshold==thr) & (df_homopl_binded$period==per),]
	
	#most_common_homopl_abs_freq <- most_common_homopl_abs_freq[!duplicated(df_homopl_binded$defining_mut),]
	most_common_homopl_abs_freq <- df_homopl_binded %>% distinct(defining_mut, major_lineage, .keep_all=TRUE)
	most_common_homopl_abs_freq <- most_common_homopl_abs_freq %>% slice_max(n=30, order_by=Freq_homopl) %>% arrange(desc(as.numeric(Freq_homopl)))
	return(most_common_homopl_abs_freq)
}

# TOP30 only considering S regions
get_most_common_homopls_s_regions <- function(df, thr=THRESHOLD_INTEREST, per=PERIOD_INTEREST) {
	#df_homopl_binded <- rbindlist(df)
	
	#df_homopl_binded <- df_homopl_binded[df_homopl_binded$is_clustered==1,]
	
	#df_homopl_binded <- df_homopl_binded %>% tidyr::separate(period_threshold, c("period","threshold"), sep="_")
	df_homopl_binded <- df
	df_homopl_binded$period <- as.integer(df_homopl_binded$period)
	df_homopl_binded$threshold <- as.numeric(df_homopl_binded$threshold)
	df_homopl_binded$Freq_homopl <- as.integer(df_homopl_binded$Freq_homopl)
	
	#most_common_homopl_s_regions <- df_homopl_binded[(df_homopl_binded$threshold==thr) & (df_homopl_binded$period==per),]
	
	#most_common_homopl_s_regions <- df_homopl_binded[!duplicated(df_homopl_binded$defining_mut),] # Removing duplicates (period and threshold)
	most_common_homopl_s_regions <- df_homopl_binded %>% distinct(defining_mut, major_lineage, .keep_all=TRUE)
	
	most_common_homopl_s_regions <- most_common_homopl_s_regions %>% filter(protein=="S") %>% slice_max(n=30, order_by=Freq_homopl) %>% arrange(desc(as.numeric(Freq_homopl)))
	return(most_common_homopl_s_regions)
}

bubble_plots_replacements_lineages <- function(df, show_replacements=TRUE, out_suffix) {
	
	most_common_lineages_rep <- most_common_lineages_norep <- list()
	most_common_lineages_rep_first <- most_common_lineages_norep_first <- list()
	# For each one of the major_lineage, plot the thresholds met and freq of homoplasy in each case (coloured by analysis where identified) + only first threshold met
	for(i in 1:length(major_lineages)) {
		most_common_lineages_rep[[i]] <- df %>% filter(major_lineage==major_lineages[i])
		#most_common_lineages_rep[[i]]$threshold <- factor(most_common_lineages_rep[[i]]$threshold, levels=table_names_thresholds)
		most_common_lineages_norep[[i]] <- most_common_lineages_rep[[i]]
		
		# if(show_replacements) {
		# 	p <- ggplot(most_common_lineages_rep[[i]], aes(x=as.factor(threshold), y=defining_mut, size=Freq_homopl)) +
		# 		geom_point(alpha=0.5) + labs(x="Quantile threshold", y="Homoplasy", title=glue("Lineage {major_lineages[i]}")) + theme_minimal() + theme(axis.text.y=element_text(size=5)) +
		# 		scale_size_area(name="Homoplasy frequency") #+ scale_color_manual(values=c("#7c9885","#e1ad01"), name="Analysis identified")
		# 	ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/repl_thr_lineage_{major_lineages[i]}.png"), plot=p, width=12, height=15, dpi=600, bg="white")
		# 	
		# 	p_ia <- ggplotly(p)
		# 	htmlwidgets::saveWidget(p_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/repl_thr_lineage_{major_lineages[i]}.html"))
		# 	
		# 	most_common_lineages_rep_first[[i]] <- most_common_lineages_rep[[i]]
		# 	most_common_lineages_rep_first[[i]] <- most_common_lineages_rep_first[[i]] %>% distinct() %>% group_by(defining_mut) %>% slice_min(n=1, as.numeric(threshold))
		# 	
		# 	# p2 <- ggplot(most_common_lineages_rep_first[[i]], aes(x=as.factor(threshold), y=defining_mut, size=Freq_homopl)) +
		# 	# 	geom_point(alpha=0.5) + labs(x="Quantile threshold (first detected only)", y="Homoplasy", title=glue("Lineage {major_lineages[i]}")) + theme_minimal() + theme(axis.text.y=element_text(size=5)) +
		# 	# 	scale_size_area(name="Homoplasy frequency") #+ scale_color_manual(values=c("#7c9885","#e1ad01"), name="Analysis identified")
		# 	# ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/repl_first_detected_thr_lineage_{major_lineages[i]}.png"), plot=p2, width=12, height=15, dpi=600, bg="white")
		# 	# 
		# 	# p2_ia <- ggplotly(p2)
		# 	# htmlwidgets::saveWidget(p2_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/repl_first_detected_thr_lineage_{major_lineages[i]}.html"))
		# } else {
		# 	
		# 	p3 <- ggplot(most_common_lineages_norep[[i]], aes(x=as.factor(threshold), y=prot_site, size=Freq_homopl, color=analysis_identified.y)) +
		# 		geom_point(alpha=0.5) + labs(x="Quantile threshold", y="Protein site", title=glue("Lineage {major_lineages[i]}")) + theme_minimal() + theme(axis.text.y=element_text(size=5)) +
		# 		scale_size_area(name="Homoplasy frequency") + scale_color_manual(values=c("#7c9885","#e1ad01"), name="Analysis identified")
		# 	ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/onlysite_thr_lineage_{major_lineages[i]}.png"), plot=p3, width=12, height=15, dpi=600, bg="white")
		# 	p3_ia <- ggplotly(p3)
		# 	htmlwidgets::saveWidget(p3_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/onlysite_thr_lineage_{major_lineages[i]}.html"))
		# 	
		# 	most_common_lineages_norep_first[[i]] <- most_common_lineages_norep[[i]]
		# 	most_common_lineages_norep_first[[i]] <- most_common_lineages_norep_first[[i]] %>% distinct() %>% group_by(prot_site) %>% slice_min(n=1, as.numeric(threshold))
		# 	
		# 	# p4 <- ggplot(most_common_lineages_norep_first[[i]], aes(x=as.factor(threshold), y=prot_site, size=Freq_homopl, color=analysis_identified.y)) +
		# 	# 	geom_point(alpha=0.5) + labs(x="Quantile threshold (first detected only)", y="Protein site", title=glue("Lineage {major_lineages[i]}")) + theme_minimal() + theme(axis.text.y=element_text(size=5)) +
		# 	# 	scale_size_area(name="Homoplasy frequency") + scale_color_manual(values=c("#7c9885","#e1ad01"), name="Analysis identified")
		# 	# ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/onlysite_first_detected_thr_lineage_{major_lineages[i]}.png"), plot=p4, width=12, height=15, dpi=600, bg="white")
		# 	# p4_ia <- ggplotly(p4)
		# 	# htmlwidgets::saveWidget(p4_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/onlysite_first_detected_thr_lineage_{major_lineages[i]}.html"))
		# }
		
	}
	
}

bubble_plots <- function(most_common_homopl_df, mut_type, out_suffix) {
	# prepare dfs for bubble plot selection analysis from HypHy identified mutations
	pos_sel_sites2 <- pos_sel_sites %>% select(protein, site, prot_site)
	pos_sel_sites2$analysis_identified <- "HyPhy"
	pos_sel_sites2$Freq_homopl <- NA
	pos_sel_sites2$period <- NA
	pos_sel_sites2$threshold <- NA
	pos_sel_sites2$major_lineage <- NA
	pos_sel_sites2 <- pos_sel_sites2 %>% select(period, threshold, protein, site, prot_site, major_lineage, analysis_identified, Freq_homopl)
	
	most_common_homopl_all <- most_common_homopl_df #most_common_homopl_df[[2]]
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
	joined_most_common_indep_ok[joined_most_common_indep_ok == "HyPhy;TreeBasedClustering"] <- "Intersection"
	
	if(mut_type=="non-syn") {
		joined_most_common_indep_ok$protein <- factor(joined_most_common_indep_ok$protein, levels=lvls_bubble)
		#joined_most_common_indep_ok$protein <- factor(joined_most_common_indep_ok$protein, levels=joined_most_common_indep_ok$protein[ order( unique(as.numeric(joined_most_common_indep_ok$site)) ) ])
		joined_most_common_indep_ok <- joined_most_common_indep_ok[order(joined_most_common_indep_ok$protein, as.numeric(joined_most_common_indep_ok$site) ),]
		joined_most_common_indep_ok$prot_site <- factor(joined_most_common_indep_ok$prot_site, levels=unique( joined_most_common_indep_ok$prot_site[order(joined_most_common_indep_ok$protein, as.numeric(joined_most_common_indep_ok$site) )] ))
		#View(joined_most_common_indep_ok)
	}
	
	
	system(glue("mkdir -p stat_results/{out_folder}/comparison_bubble_{out_suffix}/"))
	
	pal_id <- c("TreeBasedClustering"="#7c9885","Intersection"="#e1ad01", "HyPhy"="#033f63")
	
	p <- ggplot(joined_most_common_indep_ok, aes(y=analysis_identified.y, x=prot_site, color=analysis_identified.y)) + #color=analysis_identified.y
		geom_point(alpha=0.8) + labs(y="Analyses identified",x="Protein site") + theme_minimal() +
		scale_color_manual(values=pal_id, name="Analysis identified") + theme(legend.position="none", axis.text.x=element_text(size=6, angle=90, vjust = 0.5, hjust=1,color="black"), axis.text.y=element_text(size=7,color="black"), axis.title=element_text(size=9)) #, aspect.ratio=16/9
	#ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/homopl_method_combinations.png"), plot=p, width=12, height=5, dpi=600, bg="white")
	p_ia <- ggplotly(p)
	htmlwidgets::saveWidget(p_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/homopl_method_combinations.html"))
	
	joined_most_common_indep_ok_high_freq_lin <- joined_most_common_indep_ok %>% group_by(major_lineage, prot_site) %>% slice_max(Freq_homopl)
	
	# print( sort(unique(joined_most_common_indep_ok_high_freq_lin$Freq_homopl)) )
	
	# Plot max frequency of homoplasic sites (bubble sizes) for each one of the major_lineage (coloured by analysis where identified)
	#View(joined_most_common_indep_ok_high_freq_lin)
	joined_most_common_indep_ok_high_freq_lin <- joined_most_common_indep_ok_high_freq_lin[joined_most_common_indep_ok_high_freq_lin$analysis_identified.y != "HyPhy",]
	joined_most_common_indep_ok_high_freq_lin <- joined_most_common_indep_ok_high_freq_lin[!is.na(joined_most_common_indep_ok_high_freq_lin$major_lineage),]
	
	p2 <- ggplot(joined_most_common_indep_ok_high_freq_lin, aes(x=major_lineage, y=prot_site, size=Freq_homopl, color=analysis_identified.y)) +
		geom_point(alpha=0.5) + labs(x="Major lineage", y="Protein site") + theme_minimal() + theme(axis.text.y = element_text(size=6,color="black"), axis.text.x = element_text(size=7, color="black"), axis.title=element_text(size=9)) + #aspect.ratio=16/9
		scale_size_area(name="TFP-homoplasy frequency", breaks=c(2,5,10,15,20)) + scale_color_manual(values=pal_id, name="Analysis identified") #c(1,50,150,500,1300)
	#ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_all_lineages.png"), plot=p2, width=12, height=15, dpi=600, bg="white")
	p2_ia <- ggplotly(p2)
	htmlwidgets::saveWidget(p2_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_all_lineages.html"))
	
	ggarrange(p, p2, nrow=2, ncol=1, labels=c('A', 'B'), heights=c(1,2))
	ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/compHyPhy_lineages.pdf"), width=9, height=8, dpi=600, bg="white")
	
	most_common_homopl_all_high_freq_lin <- most_common_homopl_all %>% group_by(major_lineage, defining_mut) %>% slice_max(Freq_homopl)
	most_common_homopl_all_high_freq_lin <- most_common_homopl_all_high_freq_lin[!is.na(most_common_homopl_all_high_freq_lin$major_lineage),]
	# Plot max frequency of each homoplasy across thresholds (including WT and replaced aa) identified (without comparison with HyPhy)
	p3 <- ggplot(most_common_homopl_all_high_freq_lin, aes(x=major_lineage, y=defining_mut, size=Freq_homopl)) + #color=analysis_identified.y
		geom_point(alpha=0.5) + labs(x="Major lineage", y="Homoplasy") + theme_minimal() + theme(axis.text.y = element_text(size=5), aspect.ratio=16/9) +
		scale_size_area(name="Highest freq of homopl") # breaks=c(20,50,150,500,1300) + scale_color_manual(values=c("#7c9885","#e1ad01", "#033f63")) +
	ggsave(glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_replacement_lineages.png"), plot=p3, width=12, height=15, dpi=600, bg="white")
	p3_ia <- ggplotly(p3)
	htmlwidgets::saveWidget(p3_ia, glue("stat_results/{out_folder}/comparison_bubble_{out_suffix}/max_freq_replacement_lineages.html"))
	
	# # Bubble plots per major lineage comparing with other results (HyPhy): just showing protein site and not actual replacements
	# bubble_plots_replacements_lineages(joined_most_common_indep_ok, show_replacements=FALSE, out_suffix)
	# # Bubble plots per major lineage  showing actual replacements and no comparison with other analyses
	# bubble_plots_replacements_lineages(most_common_homopl_all, show_replacements=TRUE, out_suffix)
	
}

# Function to facilitate comparison between periods including and excluding Omicron muts
clustering_model_fitting <- function(path_stats, out_folder) {
	
	if(dir.exists(path_stats))
		path <- path_stats
	else
		stop("Directory does not exist!")
	
	system(glue("mkdir -p stat_results/{out_folder}/"))
	
	clustered_dfs <- joined_mut_sites_clustered <- joined_mut_sites_clustered_syn <- joined_mut_sites_clustered_non_syn <- list()
	csv_cp <- thr2 <- list()
	
	clustered_df <- read.csv(glue("{path_stats}/{path_thresholds[1]}/clustered_all_df.csv"), header=T)
	
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
	joined_mut_sites_clustered <- joined_mut_sites_clustered %>% select(defining_mut, Freq_homopl, is_clustered, major_lineage, syn_non_syn, protein, aa_length, indep_found_pos_selection, genomic_index) #s_mut_region_interest
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
	
	joined_mut_sites_clustered$is_clustered_h3 <- ifelse( (as.integer(joined_mut_sites_clustered$Freq_homopl) > 3) & (joined_mut_sites_clustered$is_clustered==1), 1, 0)
	# Normalised counts (for inspection only)
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
	write.csv(csv_cp, glue("stat_results/{out_folder}/stat_ready_csvs_all_sites/df_all_sites_{path_thresholds[1]}_SYN.csv"), quote=F, row.names=F)
	csv_cp <- csv_cp[(as.integer(csv_cp$Freq_homopl) > 0) & (as.integer(csv_cp$is_clustered)==1), ]
	write.csv(csv_cp, glue("stat_results/{out_folder}/stat_ready_csvs_only_clustered/df_only_clustered_{path_thresholds[1]}_SYN.csv"), quote=F, row.names=F)
	
	# NON-SYN
	joined_mut_sites_clustered_non_syn <- joined_mut_sites_clustered[joined_mut_sites_clustered$syn == 0,]
	joined_mut_sites_clustered_non_syn$major_lineage <- as.factor( joined_mut_sites_clustered_non_syn$major_lineage )
	joined_mut_sites_clustered_non_syn$spec_regions <- as.factor(joined_mut_sites_clustered_non_syn$spec_regions)
	csv_cp2 <- joined_mut_sites_clustered_non_syn[ order(joined_mut_sites_clustered_non_syn$Freq_homopl, decreasing=T) ,]
	write.csv(csv_cp2, glue("stat_results/{out_folder}/stat_ready_csvs_all_sites/df_all_sites_{path_thresholds[1]}_NON_SYN.csv"), quote=F, row.names=F)
	csv_cp2 <- csv_cp2[(as.integer(csv_cp2$Freq_homopl) > 0) & (as.integer(csv_cp2$is_clustered)==1), ]
	write.csv(csv_cp2, glue("stat_results/{out_folder}/stat_ready_csvs_only_clustered/df_only_clustered_{path_thresholds[1]}_NON_SYN.csv"), quote=F, row.names=F)
	
	df_syn_clustered_homopl <- joined_mut_sites_clustered_syn
	df_non_syn_clustered_homopl <- joined_mut_sites_clustered_non_syn
	
	df_syn_clustered_homopl$period <- PERIOD_INTEREST
	df_non_syn_clustered_homopl$period <- PERIOD_INTEREST
	df_syn_clustered_homopl$threshold <- THRESHOLD_INTEREST
	df_non_syn_clustered_homopl$threshold <- THRESHOLD_INTEREST
	
	# DROP OUTLIERS WITHOUT MANUAL CURATION DONE YET
	df_syn_clustered_homopl <- df_syn_clustered_homopl %>% filter(!(Freq_homopl > 20 & defining_mut %in% OUTLIERS_DROP))
	df_non_syn_clustered_homopl <- df_non_syn_clustered_homopl %>% filter(!(Freq_homopl > 20 & defining_mut %in% OUTLIERS_DROP))
	
	genomewide_plot(df_syn_clustered_homopl, mut_type="syn")
	gp_ns <- genomewide_plot(df_non_syn_clustered_homopl, mut_type="non-syn")
	
	# violin_boxplot_freq_variations(df_syn_clustered_homopl, mut_type="syn")
	# violin_boxplot_freq_variations(df_non_syn_clustered_homopl, mut_type="non-syn")

	sbf_s <- stacked_bar_freqs(df_syn_clustered_homopl, mut_type="syn")
	sbf_ns <- stacked_bar_freqs(df_non_syn_clustered_homopl, mut_type="non-syn")

	# total_amount_unique_homoplasies(df_syn_clustered_homopl, mut_type="syn")
	# total_amount_unique_homoplasies(df_non_syn_clustered_homopl, mut_type="non-syn")
	# 
	# system(glue("mkdir -p stat_results/{out_folder}/only_clustered_ordered_homopl/"))
	# # SYN
	# # df_syn_clustered_homopl <- df_syn_clustered_homopl[df_syn_clustered_homopl$Freq_homopl >= 2,] # Getting only homoplasies happening at least 2 times
	# df_syn_clustered_homopl <- df_syn_clustered_homopl[order(df_syn_clustered_homopl$Freq_homopl, decreasing=T),]
	# #hist(df_syn_clustered_homopl$Freq_homopl)
	# write.csv(df_syn_clustered_homopl, glue("stat_results/{out_folder}/only_clustered_ordered_homopl/syn_clustered_homopl_ordered.csv"), quote=F, row.names=F)
	# # NON-SYN
	# # df_non_syn_clustered_homopl <- df_non_syn_clustered_homopl[df_non_syn_clustered_homopl$Freq_homopl >= 2,]
	# df_non_syn_clustered_homopl <- df_non_syn_clustered_homopl[order(df_non_syn_clustered_homopl$Freq_homopl, decreasing=T),]
	# # #hist(df_non_syn_clustered_homopl$Freq_homopl)
	# write.csv(df_non_syn_clustered_homopl, glue("stat_results/{out_folder}/only_clustered_ordered_homopl/non_syn_clustered_homopl_ordered.csv"), quote=F, row.names=F)
	#  
	# TOP30 by absolute frequency
	most_common_homopl_abs_freq_syn <- get_most_common_homopls_abs_freq(df_syn_clustered_homopl)
	system(glue("mkdir -p stat_results/{out_folder}/abs_freq/"))
	write.csv(most_common_homopl_abs_freq_syn, file=glue("stat_results/{out_folder}/abs_freq/syn_top_muts.csv"), quote=F, row.names=F)
	most_common_homopl_abs_freq_non_syn <- get_most_common_homopls_abs_freq(df_non_syn_clustered_homopl)
	write.csv(most_common_homopl_abs_freq_non_syn, file=glue("stat_results/{out_folder}/abs_freq/nonsyn_top_muts.csv"), quote=F, row.names=F)

	# TOP30 spike
	most_common_homopl_s_non_syn <- get_most_common_homopls_s_regions(df_non_syn_clustered_homopl)
	system(glue("mkdir -p stat_results/{out_folder}/abs_freq_s/"))
	write.csv(most_common_homopl_s_non_syn, file=glue("stat_results/{out_folder}/abs_freq_s/nonsyn_top_muts_spike.csv"), quote=F, row.names=F)

	# Bubbles for is_clustered = 1
	bubble_syn_clust1 <- bubble_plots(most_common_homopl_abs_freq_syn,"syn","syn_clust1")
	bubble_non_syn_clust1 <- bubble_plots(most_common_homopl_abs_freq_non_syn, "non-syn", "non-syn_clust1")
	
	return(list(sbf_s, sbf_ns, gp_ns))
}

OUTLIERS_DROP <- c("S:T95I","S:K417N","S:N440K","S:G446S")

print("=========")
print("PERIOD 2")
print("=========")
# FOR PERIOD 2
# Expected paths containing results for each quantile threshold
thr_pref <- "threshold_quantile"
#path_thresholds <- c(glue("{thr_pref}_0.25"), glue("{thr_pref}_0.5"), glue("{thr_pref}_0.75"), glue("{thr_pref}_1"), glue("{thr_pref}_2"),
#																					glue("{thr_pref}_3"), glue("{thr_pref}_4"), glue("{thr_pref}_5"), glue("{thr_pref}_10"), glue("{thr_pref}_25"))
path_thresholds <- c(glue("{thr_pref}_2"))
table_names_periods <- c(2)
#table_names_thresholds <- c(0.25,0.5,0.75,1,2,3,4,5,10,25)
#table_names_thresholds <- c(2)
#table_combs <- do.call(paste, c(expand.grid(table_names_periods, table_names_thresholds), sep="_"))
#table_combs <- c("2")
out_folder <- "period2_thr2"
PERIOD_INTEREST <- 2
THRESHOLD_INTEREST <- 2
major_lineages <- c("EU1_B.1.177","Alpha_B.1.1.7","Delta_AY.4.*", "Delta_other","Other") #"Beta_B.1.351","Gamma_P.1"

r2 <- clustering_model_fitting("results/02_sc2_root_to_nov2021",out_folder)
# false_discovery_rate_syn("results/02_sc2_root_to_nov2021",out_folder)
# plot_tbc_stats("results/02_sc2_root_to_nov2021",out_folder)
# stacked_nsites_genomic_region_threshold("stat_results/period2_thr2/genomewide_plot_non-syn/")

print("=========")
print("PERIOD 3")
print("=========")
# FOR PERIOD 3
thr_pref <- "threshold_quantile"
#path_thresholds <- c(glue("{thr_pref}_0.25"), glue("{thr_pref}_0.5"), glue("{thr_pref}_0.75"), glue("{thr_pref}_1"), glue("{thr_pref}_2"), 
#																					glue("{thr_pref}_3"), glue("{thr_pref}_4"), glue("{thr_pref}_5"), glue("{thr_pref}_10"), glue("{thr_pref}_25"))
path_thresholds <- c(glue("{thr_pref}_2"))
#table_names_periods <- c(3)
#table_names_thresholds <- c(0.25,0.5,0.75,1,2,3,4,5,10,25)
#table_names_thresholds <- c(2)
#table_combs <- do.call(paste, c(expand.grid(table_names_periods, table_names_thresholds), sep="_"))
#table_combs <- c("2")
out_folder <- "period3_thr2"
PERIOD_INTEREST <- 3
THRESHOLD_INTEREST <- 2
major_lineages <- c("EU1_B.1.177","Alpha_B.1.1.7","Delta_AY.4.*", "Delta_other","Omicron_BA.1.*","Omicron_BA.2.*","Other") #"Beta_B.1.351","Gamma_P.1","Omicron_BA.*",

r3 <- clustering_model_fitting("results/03_sc2_whole_period",out_folder)
# false_discovery_rate_syn("results/03_sc2_whole_period",out_folder)
# plot_tbc_stats("results/03_sc2_whole_period",out_folder)
# stacked_nsites_genomic_region_threshold("stat_results/period3_thr2/genomewide_plot_non-syn/")

# Figure 3: stacked bar counts of unique homoplasies total & norm
r2[[2]][[1]] <- r2[[2]][[1]] + scale_x_continuous(limits=c(0, 110))
r3[[2]][[1]] <- r3[[2]][[1]] + scale_x_continuous(limits=c(0, 110))
my_legend <- get_legend(r3[[2]][[1]])
#ggarrange(legend_1, legend_2, legend_3, nrow=3)
ggarrange(r2[[2]][[1]], r2[[2]][[2]], r3[[2]][[1]], r3[[2]][[2]], align='h', nrow=2, ncol=2, labels=c('A', 'B','C','D'), legend.grob=my_legend, legend="right") #common.legend = T
ggsave("stat_results/plots_paper/Fig3_uniqueHomoplasiesRegions.pdf", width=9, height=8, dpi=600, bg="white")

# Figure 5: spike-wide frequency of TFP-homoplasies & number of identified sites for each genomic position stacked by threshold

# sngrt_r3 variable from analysis_mult_thresholds.R
sngrt_r3 <- readRDS("rds/sngrt_r3.rds")
ggarrange(r3[[3]], sngrt_r3, nrow=2, ncol=1, labels=c('A', 'B'), legend="right") #common.legend = T
ggsave("stat_results/plots_paper/Fig5_freqsSpike_sitesThr.pdf", width=9, height=8, dpi=600, bg="white")
