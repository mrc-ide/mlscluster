
#' Extract all polymorphisms (mutations) in a database for a determined period
#' 
#' @param sc2_md_curated_rds Same as the `amd` parameter from [mlsclust()]
#' @param cut_date Only extract mutations for the samples before (and including) this date (Date object, default: NULL).
#'    E.g. `cut_date=as.Date("2021-11-15")`
#'
#' @return A list of 2 elements:
#'    * `cog_md_muts` is the data.frame with the protein/type and mutation
#'    * `nseqs` is the number of rows (sequences) in the `sc2_md_curated_rds` loaded data.frame
#' @export
extract_muts_period <- function(sc2_md_curated_rds, cut_date) {
	sc2_md_curated <- readRDS(sc2_md_curated_rds)
	sc2_md_curated$mutations <- toupper(sc2_md_curated$mutations)
	sc2_md_curated <- sc2_md_curated %>% dplyr::select(sequence_name, sample_date, major_lineage, mutations)
	sc2_md_curated <- sc2_md_curated[sc2_md_curated$sample_date <= cut_date,]
	nseqs <- nrow(sc2_md_curated)
	
	cog_md_muts <- splitstackshape::cSplit(sc2_md_curated, 'mutations', '|', 'long')
	cog_md_muts$type <- sub("\\:.*", "", cog_md_muts$mutations)
	cog_md_muts$mut <- sub('.*:', "", cog_md_muts$mutations)
	cog_md_muts$mut_and_type <- paste0(cog_md_muts$type,":",cog_md_muts$mut)
	
	cog_md_muts$type  <- factor(cog_md_muts$type, levels = c("SYNSNP","ORF1AB","S","ORF3A","E","M","ORF6","ORF7A","ORF8","N","ORF10"))
	return(list(cog_md_muts, nseqs))
}

#' Summarise mutation counts/percents, and extract codon positions of mutations
#'
#' @param cog_md_muts Data.frame of mutations returned by [extract_muts_period()] (list position 1).
#' @param nseqs Number of sequences in the database returned by [extract_muts_period()] (list position 2).
#'
#' @return A list of 3 elements: 
#'    * `cog_md_muts_syn_ranges` generates a data.frame with the respective 
#'    genomic regions / NSPs that mutations occur with coordinates of start and end
#'    * `codon_3rd_pos` gives 3rd codon mutations found in >100 sequences of the provided database
#'    * `codon_1st_2nd_pos` returns the same as `codon_3rd_pos` but for both 1st and 2nd codon mutations 
#'    instead of only 3rd position
#' @export
compute_muts_match_codons <- function(cog_md_muts, nseqs) {
	# Get counts and percents
	cog_md_muts <- cog_md_muts %>% dplyr::group_by(mut_and_type) %>% dplyr::summarise(n=dplyr::n(), percent = 100*(n / nseqs)) %>% dplyr::arrange(dplyr::desc(percent))
	#cog_md_muts <- cog_md_muts[order(cog_md_muts$percent,decreasing=TRUE),]
	cog_md_muts <- cog_md_muts[cog_md_muts$n > 100,] #when filter 1000 -> 48394 to 2073; when filter 100 -> 48394 to 9456
	cog_md_muts$protein <- sub("\\:.*", "", cog_md_muts$mut_and_type)
	cog_md_muts$mut <- sub('.*:', "", cog_md_muts$mut_and_type)
	cog_md_muts$mutsite <- stringr::str_sub(cog_md_muts$mut_and_type, -1)
	cog_md_muts$genomic_index <- readr::parse_number(cog_md_muts$mut)
	# Remove X/N missing sites from NON-SYN and SYN respectively
	cog_md_muts_non_xn <- cog_md_muts[(cog_md_muts$protein != "SYNSNP" & cog_md_muts$mutsite != "X") | (cog_md_muts$protein == "SYNSNP" & cog_md_muts$mutsite != "N"),]
	
	# Filter only SYN
	cog_md_muts_non_xn_syn <- cog_md_muts_non_xn[cog_md_muts_non_xn$protein == "SYNSNP",]
	# Get exact genomic region (incl NSPs) from mutation
	
	syn_ranges_df <- utils::read.csv(system.file("extdata", "genomic_ranges_plot_syn.tsv", package="mlscluster"),sep="\t", header=T)
	data.table::setDT(cog_md_muts_non_xn_syn); data.table::setDT(syn_ranges_df)
	cog_md_muts_syn_ranges <- cog_md_muts_non_xn_syn[syn_ranges_df, on=c("protein==name_cog_md","genomic_index>=start","genomic_index<=end"), genomic_region := new_prot_name]
	
	# dfs of codon position (1 to 3) and respective nt coordinates
	codon_1st_pos <- syn_ranges_df %>% dplyr::rowwise() %>% dplyr::transmute(new_prot_name, codon_pos=1, coord = list(seq(start, end, by=3))) %>% tidyr::unnest_longer(coord)
	codon_2nd_pos <- syn_ranges_df %>% dplyr::rowwise() %>% dplyr::transmute(new_prot_name, codon_pos=2, coord = list(seq(start+1, end, by=3))) %>% tidyr::unnest_longer(coord)
	codon_3rd_pos <- syn_ranges_df %>% dplyr::rowwise() %>% dplyr::transmute(new_prot_name, codon_pos=3, coord = list(seq(start+2, end, by=3))) %>% tidyr::unnest_longer(coord)
	all_codon_pos <- rbind(codon_1st_pos, codon_2nd_pos, codon_3rd_pos)
	
	# paste 2 columns to make it easier to merge
	cog_md_muts_syn_ranges$prot_site <- paste0(cog_md_muts_syn_ranges$genomic_region,":",cog_md_muts_syn_ranges$genomic_index) 
	all_codon_pos$prot_site <- paste0(all_codon_pos$new_prot_name,":",all_codon_pos$coord)
	
	# Get polymorphic sites only with respective codon positions (all codons)
	muts_codons <- base::merge(cog_md_muts_syn_ranges, all_codon_pos, by="prot_site")
	nrow(muts_codons[muts_codons$codon_pos == 1,])
	nrow(muts_codons[muts_codons$codon_pos == 2,])
	nrow(muts_codons[muts_codons$codon_pos == 3,])
	
	# Get 3rd codon polymorphic sites only with respective codon positions
	codon_3rd_pos$prot_site <- paste0(codon_3rd_pos$new_prot_name,":",codon_3rd_pos$coord)
	codon_3rd_pos <- base::merge(cog_md_muts_syn_ranges, codon_3rd_pos, by="prot_site")
	
	# Get 2nd codon polymorphic sites only with respective codon positions
	codon_2nd_pos$prot_site <- paste0(codon_2nd_pos$new_prot_name,":",codon_2nd_pos$coord)
	codon_2nd_pos <- base::merge(cog_md_muts_syn_ranges, codon_2nd_pos, by="prot_site")
	
	# Get 1st codon polymorphic sites only with respective codon positions
	codon_1st_pos$prot_site <- paste0(codon_1st_pos$new_prot_name,":",codon_1st_pos$coord)
	codon_1st_pos <- base::merge(cog_md_muts_syn_ranges, codon_1st_pos, by="prot_site")
	
	# Get 1st & 2nd codon polymorphic sites only with respective codon positions
	codon_1st_2nd_pos <- rbind(codon_1st_pos, codon_2nd_pos)
	return(list(cog_md_muts_syn_ranges, codon_3rd_pos, codon_1st_2nd_pos))
}


#' Compute the false discovery rate (FDR) of \emph{mlscluster}.
#' 
#' @description A key assumption is that synonymous mutations do not provide a fitness advantage
#'     and would then represent an erroneous TFP call. Since mutations occur much more frequently 
#'     in 3rd codon positions, FDR is calculated only considering the proportion of erroneous calls 
#'     across polymorphic 3rd codon mutated sites. A more conservative approach (epsilon) considers 
#'     both the FDR and the erroneous proportion of calls on sites at first and second codon positions.
#'     
#' @param path_stats Path with the output files generated by [run_diff_thresholds()]. Identical to the supplied 
#'     `output_dir` in the [stats_multiple_thresholds()] function (default: NULL)
#' @param cog_md_muts_syn_ranges First list element returned by [compute_muts_match_codons()]
#' @param codon_3rd_pos Second list element returned by [compute_muts_match_codons()]
#' @param codon_1st_2nd_pos Third list element returned by [compute_muts_match_codons()]
#' @param out_folder Simular to parameter with same name from [stats_multiple_thresholds()]. Outputs will be generated
#'     inside `stat_results/{out_folder}/fdr_codons/`
#'
#' @return Arranged plots of FDR and epsilon overall and across different genomic regions
#' @export
fdr_mlsclust <- function(path_stats, cog_md_muts_syn_ranges, codon_3rd_pos, codon_1st_2nd_pos, out_folder) {
	
	clustered_dfs <- clustered_dfs_syn <- first_codon_threshold <- sec_codon_threshold <- third_codon_threshold <- first_sec_codon_threshold <- fdr <- epsilon <- list()
	fdr_regions <- epsilon_regions <- matrix(list(),nrow=length(thr), ncol=length(unique(cog_md_muts_syn_ranges$genomic_region)))
	#print(dim(fdr_regions))
	# clustered_dfs <- remove_homopl_freq_outliers(path_stats)
	
	print("Calculating overall and gene-specific FDR and epsilon:")
	for(i in 1:length(path_thresholds)) {
		clustered_dfs[[i]] <- utils::read.csv(glue::glue("{path_stats}/{path_thresholds[i]}/clustered_all_df.csv"), header=T)
		print("threshold")
		print(thr[i])
		clustered_dfs[[i]]$syn <- ifelse(clustered_dfs[[i]]$syn == as.character("syn"), 1, 0)
		clustered_dfs_syn[[i]] <- clustered_dfs[[i]][(clustered_dfs[[i]]$syn == 1) & (clustered_dfs[[i]]$is_clustered==1),]
		# print("nrow")
		# print(nrow(clustered_dfs_syn[[i]]))
		
		# Extract genome index position
		rgx_genomic_idx <- regexpr(':[A-Z]{1}[0-9]{1,5}[A-Z*]{1}', clustered_dfs_syn[[i]]$defining_mut)
		comm_coord <- rep(NA,length(clustered_dfs_syn[[i]]$defining_mut))
		comm_coord[rgx_genomic_idx != -1] <- stringr::str_sub( regmatches(clustered_dfs_syn[[i]]$defining_mut, rgx_genomic_idx), 3, -2)
		comm_coord <- as.integer(comm_coord)
		comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("genomic_index")
		clustered_dfs_syn[[i]] <- cbind(clustered_dfs_syn[[i]], comm_coord_df)
		
		clustered_dfs_syn[[i]]$protein <- sub("\\:.*", "", clustered_dfs_syn[[i]]$defining_mut)
		
		#View(clustered_dfs_syn[[i]])
		
		annot_gen_range_positions <- utils::read.csv(system.file("extdata", "genomic_ranges_plot_syn.tsv", package="mlscluster"), sep="\t", header=T)
		data.table::setDT(clustered_dfs_syn[[i]]); data.table::setDT(annot_gen_range_positions)
		clustered_dfs_syn[[i]] <- clustered_dfs_syn[[i]][annot_gen_range_positions, on=c("protein==name_cog_md","genomic_index>=start","genomic_index<=end"), genomic_region := new_prot_name]
		# clustered_dfs_syn[[i]] <- clustered_dfs_syn[[i]][annot_gen_length_positions, on=c("spec_regions==new_prot_name"), aa_length := length]
		
		clustered_dfs_syn[[i]]$prot_site <- paste0(clustered_dfs_syn[[i]]$genomic_region,":",clustered_dfs_syn[[i]]$genomic_index)
		
		# Merge with 3rd codon positions
		third_codon_threshold[[i]] <- base::merge(codon_3rd_pos, clustered_dfs_syn[[i]], by="prot_site", all.x=T)
		third_codon_threshold[[i]]$threshold <- thr[i]
		
		# success = call = 1
		# failure = no call = 0
		# if NA, then no call
		third_codon_threshold[[i]]$success <- ifelse(!is.na(third_codon_threshold[[i]]$defining_mut), yes=1, no=0)
		third_codon_threshold[[i]]$new_prot_name <- sub('.*:', "", third_codon_threshold[[i]]$new_prot_name)
		third_codon_threshold[[i]] <- third_codon_threshold[[i]] %>% dplyr::select(new_prot_name, threshold, success)
		third_codon_threshold[[i]]$threshold <- factor(third_codon_threshold[[i]]$threshold, levels=thr)
		
		# Model for the probability (bernoulli outcome) y ~ threshold , where y is the outcome of being classified a TFP specifically among all 
		# polymorphic 3rd codon position sites in the genome.
		
		# Pick a site at random. The probability that it is called a TFP (success) is n_sites_called (success==1) / n_polymorphic_sites_3rd_codon
		
		# FDR = E[y(threshold)] (expected value equals to probability of success, erroneous call)
		# FDR is the estimated probability that a given site with at least one mutation (ie  a polymorphic site) in the database is erroneously called a TFP
		fdr[[i]] <- nrow(third_codon_threshold[[i]][third_codon_threshold[[i]]$success == 1,]) / nrow(third_codon_threshold[[i]])
		fdr[[i]] <- fdr[[i]]*100
		fdr[[i]] <- as.data.frame(fdr[[i]]); fdr[[i]]$threshold <- thr[i]; colnames(fdr[[i]]) <- c("FDR","threshold")
		#print(fdr[[i]])
		
		gen_regions <- unique(third_codon_threshold[[i]]$new_prot_name)
		
		# FDR for each genomic region
		for(j in 1:length(gen_regions)) {
			print(gen_regions[j])
			fdr_regions[[i,j]] <- third_codon_threshold[[i]][ third_codon_threshold[[i]]$new_prot_name == gen_regions[j], ]
			fdr_regions[[i,j]] <- nrow(fdr_regions[[i,j]][fdr_regions[[i,j]]$success == 1,]) / nrow(fdr_regions[[i,j]])
			fdr_regions[[i,j]] <- fdr_regions[[i,j]]*100
			fdr_regions[[i,j]] <- as.data.frame(fdr_regions[[i,j]]); fdr_regions[[i,j]]$threshold <- thr[i]; fdr_regions[[i,j]]$gen_region <- gen_regions[j]; colnames(fdr_regions[[i,j]]) <- c("FDR","threshold","gen_region")
			#print(fdr_regions[[i,j]])
		}
		
		# Separate error rate can be computed relative to the sites at codon positions 1-2
		# epsilon = FDR* <number TFP calls at positions 1-2> / <number polymorphic sites at position 1-2>
		
		first_sec_codon_threshold[[i]] <- base::merge(codon_1st_2nd_pos, clustered_dfs_syn[[i]], by="prot_site", all.x=T)
		first_sec_codon_threshold[[i]]$threshold <- thr[i]
		first_sec_codon_threshold[[i]]$success <- ifelse(!is.na(first_sec_codon_threshold[[i]]$defining_mut), yes=1, no=0)
		first_sec_codon_threshold[[i]]$new_prot_name <- sub('.*:', "", first_sec_codon_threshold[[i]]$new_prot_name)
		first_sec_codon_threshold[[i]] <- first_sec_codon_threshold[[i]] %>% dplyr::select(new_prot_name, threshold, success)
		
		epsilon[[i]] <- fdr[[i]] * nrow(first_sec_codon_threshold[[i]][first_sec_codon_threshold[[i]]$success == 1,]) / nrow(first_sec_codon_threshold[[i]])
		epsilon[[i]] <- as.data.frame(epsilon[[i]]); epsilon[[i]]$threshold <- thr[i]; colnames(epsilon[[i]]) <- c("epsilon","threshold")
		
		# epsilon for each genomic region
		for(j in 1:length(gen_regions)) {
			print(gen_regions[j])
			epsilon_regions[[i,j]] <- first_sec_codon_threshold[[i]][ first_sec_codon_threshold[[i]]$new_prot_name == gen_regions[j], ]
			epsilon_regions[[i,j]] <- fdr_regions[[i,j]]$FDR[1] * nrow(epsilon_regions[[i,j]][ epsilon_regions[[i,j]]$success == 1, ]) / nrow(epsilon_regions[[i,j]])
			epsilon_regions[[i,j]] <- as.data.frame(epsilon_regions[[i,j]]); epsilon_regions[[i,j]]$threshold <- thr[i]; epsilon_regions[[i,j]]$gen_region <- gen_regions[j]; colnames(epsilon_regions[[i,j]]) <- c("epsilon","threshold","gen_region")
		}
	}
	
	# Overall FDR and epsilon
	fdr_join <- data.table::rbindlist(fdr)
	epsilon_join <- data.table::rbindlist(epsilon)
	
	fdr_epsilon_join <- fdr_join %>% dplyr::inner_join(epsilon_join, by="threshold")
	fdr_epsilon_join <- reshape2::melt(fdr_epsilon_join, id="threshold")  
	
	system(glue::glue("mkdir -p stat_results/{out_folder}/fdr_codons/"))
	p_overall <- ggplot2::ggplot(fdr_epsilon_join, ggplot2::aes(x=as.factor(threshold), y=value, colour=variable, group=variable)) +
		ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::ggtitle("Overall") + ggplot2::theme_minimal() +
		ggplot2::labs(x="Quantile thresholds", y="FDR estimates (%)") + 
		ggplot2::coord_cartesian(ylim=c(0,60)) +
		ggplot2::scale_y_continuous(breaks=c(5,10,15,20,40,60)) +
		ggplot2::theme(axis.text = ggplot2::element_text(size=5,color="black"), plot.title=ggplot2::element_text(size=6), axis.title=ggplot2::element_text(size=7,color="black"))
	ggplot2::ggsave(glue::glue("stat_results/{out_folder}/fdr_codons/fdr.png"), plot=p_overall, width=5, height=4, dpi=600, bg="white")
	
	# FDR and epsilon for each region
	fdr_regions_join <- do.call(rbind, fdr_regions)
	
	epsilon_regions_join <- do.call(rbind, epsilon_regions)
	
	fdr_epsilon_regions_join <- fdr_regions_join %>% dplyr::inner_join(epsilon_regions_join, by=c("threshold","gen_region"))
	fdr_epsilon_regions_join <- reshape2::melt(fdr_epsilon_regions_join, id=c("threshold","gen_region")) 
	
	gen_regions <- unique(fdr_epsilon_regions_join$gen_region)
	gen_regions <- gen_regions[base::order(match(gen_regions,lvls_nonsyn_proteins))]
	gen_regions <- gen_regions[ !gen_regions == 'ORF7B'] #remove ORF7b (missing data)
	
	# NSP5, NSP16 & ORF8 do not have TFP calls overlapping with 1st and 2nd codon sites (epsilon = 0)
	# ORF10 does not have 1st+2nd polymorphic sites with Freq>100 (epsilon = NA)
	p_fdr <- list()
	for(k in 1:length(gen_regions)) {
		p_fdr[[k]] <-	fdr_epsilon_regions_join %>% dplyr::filter(gen_region == gen_regions[k]) %>%
			ggplot2::ggplot(ggplot2::aes(x=as.factor(threshold),y=value, colour=variable, group=variable)) +
			ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::ggtitle(gen_regions[k]) + ggplot2::theme_minimal() +
			ggplot2::coord_cartesian(ylim=c(0,60)) +
			ggplot2::scale_y_continuous(breaks=c(5,10,15,20,40,60)) +
			ggplot2::labs(x="Quantile thresholds", y="FDR estimates (%)") + ggplot2::theme(axis.text = ggplot2::element_text(size=5,color="black"), plot.title=ggplot2::element_text(size=6), axis.title=ggplot2::element_text(size=7,color="black"))
	}
	p_all <- c(list(p_overall),p_fdr)
	p_all
}
