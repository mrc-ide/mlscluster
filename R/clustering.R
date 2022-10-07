load("rds/chron_timetree_with_regional_and_muts_md.RData")

cladeScore <- function(tre, amd, min_descendants=10, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=NULL, max_date=NULL, branch_length_unit="years", output_dir=paste0("clade_score-",Sys.Date()), quantile_choice=1/100, quantile_threshold_ratio_sizes="1%", quantile_threshold_ratio_persist_time="1%", threshold_keep_lower=TRUE, defining_mut_threshold=0.75, root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995, compute_tree_annots=FALSE, plot_global_mut_freq=FALSE) { #ncpu=1
	
	library(ape)
	library(lubridate)
	library(glue)
	library(ggtree)
	library(ggplot2)
	library(data.table)
	library(dplyr)
	library(tidyr)
	library(stringr)
	library(outbreakinfo) #needs libudunits2-dev and libgdal-dev (Ubuntu: sudo apt install)
	
	# if(ncpu > 1) {
	# 	library(foreach)
	# 	library(doParallel)
	# }
	if(plot_global_mut_freq)
		outbreakinfo::authenticateUser()
	
	class(tre) <- "phylo"
	
	if (!dir.exists(output_dir))
		suppressWarnings( dir.create(output_dir) )
	
	max_time <- Inf
	if(!is.null(max_date))
		max_time <- decimal_date(max_date)
	else
		max_time <- Sys.Date()
	
	min_time <- -Inf
	if(!is.null(min_date))
		min_time <- decimal_date(min_date)
	
	max_time_full <- max(decimal_date(as.Date(amd$sample_date)))
	min_time_full <- min(decimal_date(as.Date(amd$sample_date)))
	diff_stimes <- max_time_full - min_time_full
	
	# load tree and md data
	message("Loading tree and metadata and filtering out...")
	amd <- amd[!is.na(amd$sequence_name),]
	amd$sample_date <- as.Date(amd$sample_date)
	stopifnot(all(tre$tip.label %in% amd$sequence_name)) # enforce that all tips are listed on metadata sequence_name column
	if(!('sample_time') %in% colnames(amd)) {
		amd$sample_time <- decimal_date(amd$sample_date)
	}
	if(!('lineage' %in% colnames(amd))) {
		amd$lineage <- 'lineage_not_provided'
	}
	if(!('region' %in% colnames(amd))) {
		amd$region <- 'region_not_provided'
	}
	
	#amd$sts <- amd$sample_time
	# remove missing dates and filter by sample times
	amd <- amd[!is.na(amd$sample_time),]
	amd <- amd[(amd$sample_time >= min_time) & (amd$sample_time <= max_time),]
	
	#if(!is.rooted(tre)) {
	if(!(root_on_tip %in% amd$sequence_name)) {
		stopifnot(root_on_tip %in% tre$tip.label)
		# amd <- rbind(data.frame( sequence_name=root_on_tip, sample_date=date_decimal(root_on_tip_sample_time), lineage=NA, region=NA, sample_time=root_on_tip_sample_time ))
		amd[nrow(amd) + 1,] <- data.frame( sequence_name=root_on_tip, sample_date=date_decimal(root_on_tip_sample_time), lineage=NA, region=NA, sample_time=root_on_tip_sample_time )	
	}
	#}
	message(paste0("Number of sequences included: ", nrow(amd)))
	
	# named vector of sequence_names and sample_times
	sts <- setNames(amd$sample_time , amd$sequence_name)
	amd <- amd[amd$sequence_name %in% tre$tip.label,]
	
	dropped_tips <- setdiff(tre$tip.label, amd$sequence_name)
	message(paste0("Dropped tips: ", length(dropped_tips)))
	#message("Root sequence Wuhan present?")
	#print("Wuhan/WH04/2020" %in% intersect(tre$tip.label, amd$sequence_name))
	
	tre <- keep.tip(tre, intersect(tre$tip.label, amd$sequence_name))
	tre2 <- tre
	if(!is.rooted(tre)) {
		stopifnot(root_on_tip %in% tre$tip.label)
		if(!(root_on_tip %in% tre$tip.label)) {
			stop("Outgroup sequence missing in the input tree")
		}
		tre2 <- root(tre, outgroup=root_on_tip, resolve.root=TRUE)
		tre <- tre2
	}
	message(paste0("Tips in the tree: ", length(tre$tip.label)))
	
	# variables and data structures to quickly look up tree data
	# based on treestructure/tfpscanner code
	message("Computing node descendants and times...")
	ntip <- Ntip(tre)
	nnode <- Nnode(tre)
	
	# getting 'coalescent' times
	if(branch_length_unit == "years") {
		nde <- node.depth.edgelength(tre) # depth of a node using branch lengths
	}else if(branch_length_unit == "days") {
		tre$edge.length <- tre$edge.length / 365
		nde <- node.depth.edgelength(tre)
	}else {
		stop("Choices for 'branch_length_unit' are: 'days' and 'years'")
	}
	rh <- max(nde[1:ntip]) # root heights
	stimes <- nde[1:ntip]
	shs <- rh - stimes # time to most recent sample
	nhs = nodeheights <- rh - nde # node heights
	mrst <- max(sts) # most recent sample time
	tmrca <- mrst - nhs # decimal times of coalescent events
	hist(tmrca)
	
	poedges <- tre$edge[postorder(tre),] #postorder: traverse from the left subtree to the right subtree then to the root
	preedges <- tre$edge[rev(postorder(tre)),] #preorder: traverse from the root to the left subtree then to the right subtree
	
	# number of descendants
	ndesc <- rep(0, ntip+nnode) # tips descending
	ndesc[1:ntip] <- 1 # tip indices will have only 1 descendant
	Ndesc <- rep(1, ntip+nnode) # including internal nodes
	
	# compute times (max and min) of descendants
	tlsts <- sts[ tre$tip.label ]
	max_desc_time <- rep(-Inf, ntip+nnode)
	max_desc_time[1:ntip] <- tlsts 
	min_desc_time <- rep(Inf, ntip+nnode)
	min_desc_time[1:ntip] <- tlsts
	persistence_time <- rep(-Inf, ntip+nnode)
	persistence_time[1:ntip] <- 0
	
	for (e in 1:nrow(poedges)){
		n <- poedges[e,1] # postorder traversal of nodes column
		tp <- poedges[e,2] # postorder traversal of tips column
		ndesc[n] <- ndesc[n] + ndesc[tp]
		Ndesc[n] <- Ndesc[n] + Ndesc[tp]
		max_desc_time[n] <- max(max_desc_time[n], max_desc_time[tp])
		min_desc_time[n] <- min(min_desc_time[n], min_desc_time[tp])
		persistence_time[n] <- max_desc_time[n] - tmrca[n]
	}
	clade_age <- max_desc_time - min_desc_time # time span of descendant tips
	
	# ancestors
	ancestors <- lapply(1:(ntip+nnode), function(tp) integer() )
	for(e in 1:nrow(preedges)) {
		n <- preedges[e,1]
		tp <- preedges[e,2]
		ancestors[[tp]] <- c(ancestors[[n]], n)
	}
	
	# descendant nodes including internal ones
	descendants <- lapply(1:(ntip+nnode), function(tp) integer(Ndesc[tp]))  # pre-allocate list to fill with node descendant IDs
	for(tp in 1:(ntip+nnode)) {
		descendants[[tp]][1] <- tp
	}
	
	Ndesc_idx <- rep(2, ntip+nnode)
	
	for(e in 1:nrow(poedges)) {
		n <- poedges[e,1]
		tp <- poedges[e,2]
		i0 <- Ndesc_idx[n]
		i1 <- Ndesc[tp] + i0 - 1
		descendants[[n]][i0:i1] <- descendants[[tp]]
		Ndesc_idx[n] <- i1 + 1
	}
	
	# descendant tip indices
	descendant_tips <- descendants
	# daughter nodes: direct descendants of a given node
	daughter_nodes <- descendants
	
	for(n in (ntip+1):(ntip+nnode)) {
		descendant_tips[[n]] <- descendant_tips[[n]][ descendant_tips[[n]] <= ntip ]
		daughter_nodes[[n]] <- daughter_nodes[[n]][ daughter_nodes[[n]] > ntip ]
		# sisters[[n]] <- tre$edge[ tre$edge[,1] == n, 2 ] # not scaling well with larger tree
	}
	
	# tip labels
	descendant_ids <- lapply(1:(ntip+nnode), function(tp) {
		na.omit(tre$tip.label[descendant_tips[[tp]]])
	})
	
	.get_comparison_sister_node <- function(node) {
		imed_ancestor <- tail(ancestors[[node]], 1)
		sisters <- tre$edge[ tre$edge[,1] == imed_ancestor, 2 ]
		comparison_node <- setdiff(sisters, node)
		desc_tips_sisters <- lapply(1:length(sisters), function(tp) descendant_ids[[ sisters[[tp]] ]] )
		names(desc_tips_sisters) <- lapply(1:length(desc_tips_sisters), function(tp) sisters[[tp]] )
		
		sts_tips_sisters <- lapply(1:length(desc_tips_sisters), function(tp) sts[ desc_tips_sisters[[tp]] ] )
		names(sts_tips_sisters) <- lapply(1:length(sts_tips_sisters), function(tp) sisters[[ tp ]] )
		
		desc_tips_sisters_node <- unlist(unname(desc_tips_sisters[names(desc_tips_sisters) == node]))
		desc_tips_sisters_control <- unlist(unname(desc_tips_sisters[names(desc_tips_sisters) == comparison_node]))
		sts_tips_sisters_node <- unlist(unname(sts_tips_sisters[names(sts_tips_sisters) == node]))
		sts_tips_sisters_control <- unlist(unname(sts_tips_sisters[names(sts_tips_sisters) == comparison_node]))
		
		res <- list()
		if((length(desc_tips_sisters_node) >= min_descendants) && (length(desc_tips_sisters_control) >= min_descendants)) {
			res <- list(sisters, comparison_node, desc_tips_sisters, sts_tips_sisters, desc_tips_sisters_node,	desc_tips_sisters_control, sts_tips_sisters_node, sts_tips_sisters_control)
		}else {
			res <- NULL
		}
		return(res)
	}
	
	# ratio sizes function
	.ratio_sizes_stat <- function(node) {
		comp_res <- .get_comparison_sister_node(node)
		length_sisters <- rep(1,2)
		if(!is.null(comp_res)) {
			length_sisters[1] <- length(comp_res[[5]])
			length_sisters[2] <- length(comp_res[[6]])
		}else {
			length_sisters <- c(-1, -1)
		}
		
		ratio_sizes <- round( length_sisters[1] / length_sisters[2], digits=2)
		res_ratio_sizes <- c(length_sisters, ratio_sizes)
		return(res_ratio_sizes)
	}
	
	# ratio persistence time function
	.ratio_persist_time_stat <- function(node) {
		comp_res <- .get_comparison_sister_node(node)
		persist_time <- rep(0,2)
		max_persist_time <- rep(0,2)
		
		if(!is.null(comp_res)) {
			persist_time[1] <- persistence_time[ node ]
			persist_time[2] <- max( persistence_time[ comp_res[[2]] ])
			max_persist_time[1] <- max_desc_time[ node ]
			max_persist_time[2] <- max( max_desc_time[ comp_res[[2]] ])
		}else {
			persist_time <- c(-1, -1)
			max_persist_time <- c(-1, -1)
		}
		#print(persist_time)
		
		ratio_persist_time <- round( persist_time[1] / persist_time[2], digits=2)
		res_ratio_persist_time <- c(tmrca[node], max_persist_time, ratio_persist_time)
		
		return(res_ratio_persist_time)
	}
	
	# logistic growth of sister clades function
	.logistic_growth_stat <- function(node) {
		lg_res_list <- list()
		comp_res <- .get_comparison_sister_node(node)
		
		if(!is.null(comp_res)) {
			X <- data.frame(sequence_name=c(comp_res[[5]], comp_res[[6]]), time=c(comp_res[[7]], comp_res[[8]]), type=c(rep("node",length(comp_res[[5]])), rep("control",length(comp_res[[6]]))))
			X <- na.omit(X)
			model <- suppressWarnings( glm(type=="node" ~ time, data=X, family=binomial(link="logit")) )
			s <- summary(model)
			rel_growth <- unname(coef(model)[2]) # * gen_time
			p <- s$coefficients[2,4]
			lg_res_list <- cbind(rel_growth, p, comp_res[[2]]) #node, sisters[[node]][1], sisters[[node]][2],
		}else {
			lg_res_list <- cbind(-1, -1, -1)
		}
		return(lg_res_list)
	}
	
	# Extract node defining muts
	.node_muts <- function(node, mut_var="mutations") {
		comp_res <- .get_comparison_sister_node(node)
		md_itn <- amd[ amd$sequence_name %in% comp_res[[5]], ]
		md_its <- amd[ amd$sequence_name %in% comp_res[[6]], ]
		
		if( nrow(md_itn) == 0 | nrow(md_its) == 0 )
			return(list(defining=NA, all=NA))
		
		vtab_node = sort( table( do.call( c, strsplit( md_itn[[mut_var]], split='\\|' )  ) ) / nrow( md_itn ) )
		#print(vtab_node)
		vtab_sister = sort( table( do.call( c, strsplit( md_its[[mut_var]], split='\\|' )  ) ) / nrow( md_its ) )
		#print(vtab_sister)
		
		all_muts <- names(vtab_node[vtab_node > defining_mut_threshold])
		defining_muts <- setdiff( names(vtab_node[vtab_node > defining_mut_threshold]), names(vtab_sister[vtab_sister > defining_mut_threshold]) )
		
		list(defining=defining_muts, all=all_muts)
	}
	
	#quantile_options <- c(1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10)
	.extract_value_below_quantile_threshold <- function(df, var, quantile_choice, quantile_threshold) {
		if(quantile_choice >= 1/100) {
			quant <- quantile(as.numeric(var), probs=seq(0,1,quantile_choice)) #na.rm=TRUE
			message("Quantiles")
			print(quant)
			value_chosen_quant <- unname(quant[names(quant) == quantile_threshold])
			value_chosen_quant <- round(value_chosen_quant, digits=2)
			message(glue("Threshold value at quantile {quantile_threshold}"))
			message(value_chosen_quant)
			if(threshold_keep_lower)
				res_df <- df[(as.numeric(var) < value_chosen_quant),]
			else {
				res_df <- df[(as.numeric(var) > value_chosen_quant),]
			}
		}else {
			stop("Choices for quantile are: {1/2, ... , 1/100}")
		}
		return(res_df)
	}
	
	# Detect independent occurences of defining mutations (homoplasies) across different nodes
	.find_homoplasies <- function(node_list) {
		def_muts_nodes <- lapply( 1:length(node_list), function(tp) .node_muts(node_list[[tp]]) )
		names(def_muts_nodes) <- lapply( 1:length(node_list), function(tp) node_list[[tp]] )
		
		def_muts_nodes <- sapply(def_muts_nodes, "[[",1)
		def_muts_nodes <- lapply(def_muts_nodes, function(x) if(identical(x, character(0))) NA_character_ else x)
		
		def_muts_nodes_df <- rbindlist( lapply(def_muts_nodes, function(x) data.table(x)), idcol="node") 
		def_muts_nodes_df <- na.omit(def_muts_nodes_df)
		names(def_muts_nodes_df) <- c("node","defining_mut")
		def_muts_nodes_df$defining_mut <- toupper( def_muts_nodes_df$defining_mut )
		tab_def_muts_df <- setDT(def_muts_nodes_df)[, .(Freq_homopl = .N), by = .(defining_mut)]
		homoplasy_count_df <- merge(def_muts_nodes_df, tab_def_muts_df, by="defining_mut")
		homoplasy_count_df <- homoplasy_count_df[as.numeric(homoplasy_count_df$Freq_homopl) > 1,]
		homoplasy_count_df$node <- as.integer(homoplasy_count_df$node)
		
		homoplasy_count_df
	}
	
	# Annotate homoplasies in: (i) S regions of potential significance
	.annotate_s_homoplasies <- function(homopl_df) {
		regions_s <- read.csv("config/spike_regions.tsv", sep="\t", header=TRUE)
		# If has homoplasy in S
		homopl_df$s_mut <- ifelse( grepl(homopl_df$defining_mut, pattern ='^[S]:') , "Yes", "No")
		# Get S mut coordinate only
		rgx_s <- regexpr('^[S]:[A-Z]{1}[0-9]{1,4}[A-Z]{1}', homopl_df$defining_mut)
		#print(rgx_s)
		comm_coord <- rep(NA,length(homopl_df$defining_mut))
		comm_coord[rgx_s != -1] <- str_sub( regmatches(homopl_df$defining_mut, rgx_s) , 4, -2)
		comm_coord <- as.integer(comm_coord)
		comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("s_mut_coord")
		homopl_df <- cbind(homopl_df, comm_coord_df)
		setDT(homopl_df); setDT(regions_s)
		# Check if coordinate within boundaries of regions of interest (NTD, RBD, FCS)
		homopl_df[regions_s, on=c("s_mut_coord>=start", "s_mut_coord<=end"), s_mut_region_interest := region]
		
		homopl_df
	}
	
	# Annotate homoplasies in: (ii) different substitutions at the same site
	.annotate_diff_aa_mut_same_site <- function(homopl_df) {
		rgx_aadns <- regexpr('^[A-Za-z]:[A-Z]{1}[0-9]{1,4}[A-Z]{1}', homopl_df$defining_mut)
		#print(rgx_aadns)
		comm_coord_aadns1 <- rep(NA,length(homopl_df$defining_mut))
		# All genomic regions + wt aa + position (without mutated aa)
		comm_coord_aadns1[rgx_aadns != -1] <- str_sub( regmatches(homopl_df$defining_mut, rgx_aadns) , 1, -2)
		#print(comm_coord_aadns1)
		comm_coord_aadns1_df <- as.data.frame(comm_coord_aadns1); colnames(comm_coord_aadns1_df) <- c("prot_coord_without_mut_site")
		homopl_df <- cbind(homopl_df, comm_coord_aadns1_df)
		homopl_df$prot_coord_without_mut_site <- str_sub( homopl_df$defining_mut, 1, -2)
		# Assign nodes with each prot_coord_without_mut_site
		homopl_aadns_nodes <- homopl_df %>% group_by(prot_coord_without_mut_site) %>% summarize(nodes_homopl = paste0(na.omit(node), collapse="|")) %>% ungroup()
		homopl_aadns <- merge(homopl_df, homopl_aadns_nodes, by="prot_coord_without_mut_site")
		homopl_aadns$mut_site <- str_sub( homopl_aadns$defining_mut, -1, -1)
		# Filter repeated rows
		homopl_aadns <- homopl_aadns[!duplicated(homopl_aadns$defining_mut),]
		homopl_aadns_freq <- homopl_aadns %>% group_by(prot_coord_without_mut_site) %>% summarise(Freq=n(), diff_mutated_aa=paste0(na.omit(mut_site), collapse="|"))
		homopl_aadns_final <- merge(homopl_aadns, homopl_aadns_freq, by="prot_coord_without_mut_site")
		homopl_aadns_final <- homopl_aadns_final %>% select(nodes_homopl, prot_coord_without_mut_site, diff_mutated_aa, Freq)
		homopl_aadns_final <- homopl_aadns_final[!duplicated(homopl_aadns_final$prot_coord_without_mut_site),]
		homopl_aadns_final <- homopl_aadns_final[as.numeric(homopl_aadns_final$Freq) > 1,]

		homopl_aadns_final
	}
	
	# Annotate homoplasies in: (iii) moving window of neighbouring 3 residues
	.annotate_adjacent_muts_window_s3 <- function(homopl_df) {
		window_size <- 3
		
		prot_lengths <- read.csv("config/aa_length_prots.tsv", sep="\t", header=TRUE)
		prot_lengths$protein <- toupper(prot_lengths$protein)
		
		# Create moving windows for each protein size
		prot_len_ids <- list()
		prot_windows <- list()
		for(n in 1:nrow(prot_lengths)) {
			prot_len_ids[[n]] <- rep.int(prot_lengths[n,1], times=(prot_lengths[n,2]-2))
			prot_len_ids[[n]] <- as.data.frame(prot_len_ids[[n]]); colnames(prot_len_ids[[n]]) <- "protein"
			prot_windows[[n]] <- prot_len_ids[[n]]
			for(i in 1:( nrow(prot_len_ids[[n]]) -2)) {
				prot_windows[[n]][i:(i+window_size-1),] <- paste0(i,"-",(i+window_size-2),"-",(i+window_size-1))
				prot_windows[[n]] <- as.data.frame(prot_windows[[n]]); colnames(prot_windows[[n]]) <- "window"
			}
		}
		prot_len_ids_df <- rbindlist(prot_len_ids)
		prot_windows_df <- rbindlist(prot_windows)
		prot_len_window_combined <- cbind(prot_len_ids_df, prot_windows_df)
		prot_len_window_combined <- data.frame(prot_len_window_combined, do.call(rbind,str_split(prot_len_window_combined$window,"-")))
		colnames(prot_len_window_combined) <- c("protein","window","window_pos1","window_pos2","window_pos3")
		prot_len_window_combined$window_pos1 <- as.integer(prot_len_window_combined$window_pos1)
		prot_len_window_combined$window_pos2 <- as.integer(prot_len_window_combined$window_pos2)
		prot_len_window_combined$window_pos3 <- as.integer(prot_len_window_combined$window_pos3)
		prot_len_window_combined <- prot_len_window_combined %>% group_by(protein) %>% mutate(window_id = rep(seq(n()), each = window_size, length = n())) %>% ungroup()
		
		# Get protein name of homoplasies
		homopl_df$protein <- sub("\\:.*", "", homopl_df$defining_mut)
		homopl_df$protein <- toupper(homopl_df$prot)
		
		rgx_amw3 <- regexpr('^[A-Za-z]:[A-Z]{1}[0-9]{1,4}[A-Z]{1}', homopl_df$defining_mut)
		comm_coord <- rep(NA,length(homopl_df$defining_mut))
		# Mutation coordinate of homoplasies
		comm_coord[rgx_amw3 != -1] <- str_sub( regmatches(homopl_df$defining_mut, rgx_amw3) , 4, -2)
		comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("mut_coord")
		homopl_df <- cbind(homopl_df, comm_coord_df)
		homopl_df <- homopl_df[(homopl_df$protein != "SYNSNP"),]
		homopl_df$mut_coord <- as.integer(homopl_df$mut_coord)
		homopl_df <- na.omit(homopl_df)
		
		homopl_df_merge <- homopl_df %>% left_join(prot_len_window_combined, by="protein")
		
		homopl_df_merge$window_pos1_match <- ifelse(homopl_df_merge$mut_coord == homopl_df_merge$window_pos1, TRUE, FALSE)
		homopl_df_merge$window_pos2_match <- ifelse(homopl_df_merge$mut_coord == homopl_df_merge$window_pos2, TRUE, FALSE)
		homopl_df_merge$window_pos3_match <- ifelse(homopl_df_merge$mut_coord == homopl_df_merge$window_pos3, TRUE, FALSE)
		
		homopl_df_matches <- homopl_df_merge[( (homopl_df_merge$window_pos1_match == TRUE) | (homopl_df_merge$window_pos2_match == TRUE) | (homopl_df_merge$window_pos3_match == TRUE)),]
		
		homopl_df_matches_nodes <- homopl_df_matches %>% group_by(defining_mut) %>% summarize(nodes_homopl = paste0(unique(node), collapse="|"), Freq=length(unique(node))) %>% ungroup()
		homopl_df_matches_nodes <- homopl_df_matches_nodes[as.numeric(homopl_df_matches_nodes$Freq) > 1,]
		
		homopl_df_matches_neigh_mut <- homopl_df_matches %>% inner_join(homopl_df_matches_nodes, by="defining_mut")
		homopl_df_matches_neigh_mut <- homopl_df_matches_neigh_mut %>% select(nodes_homopl, protein, mut_coord, defining_mut, Freq, window, window_id)
		homopl_df_matches_neigh_mut <- homopl_df_matches_neigh_mut[order(homopl_df_matches_neigh_mut$protein, as.numeric(homopl_df_matches_neigh_mut$mut_coord)),]
		homopl_df_matches_neigh_mut$window_id2 <- paste0(homopl_df_matches_neigh_mut$protein,":",homopl_df_matches_neigh_mut$window_id)
		
		homopl_df_matches_neigh_mut_final <- homopl_df_matches_neigh_mut[!duplicated(homopl_df_matches_neigh_mut[c("defining_mut","window_id2")]),] #4,8
		homopl_df_matches_neigh_mut_final <- homopl_df_matches_neigh_mut_final %>% group_by(window_id2) %>% filter(n() != 1)
		homopl_df_matches_neigh_mut_final <- homopl_df_matches_neigh_mut_final[!duplicated(homopl_df_matches_neigh_mut_final$defining_mut),]
		homopl_df_matches_neigh_mut_final <- homopl_df_matches_neigh_mut_final %>% select(nodes_homopl, defining_mut, Freq, window, window_id, window_id2)
		
		homopl_df_matches_neigh_mut_final
	}
	
	# Flag know problematic sites listed here: https://github.com/W-L/ProblematicSites_SARS-CoV2
	.remove_known_problematic_sites <- function(homopl_df) {
		# Process problematic sites file
		probl_sites <- read.csv("config/problematic_sites_20221006.tsv", sep="\t", header=TRUE)
		# Split multiple potential mutated sites (separated by commas) for each site into unique row
		probl_sites <- probl_sites %>% mutate(aa_alt = strsplit(as.character(aa_alt), ",")) %>% unnest(aa_alt)
		probl_sites$aa_prot_site <- toupper( paste0(probl_sites$gene,":",probl_sites$aa_ref,probl_sites$aa_pos,probl_sites$aa_alt) )
		probl_sites <- probl_sites[probl_sites$filter == "mask",] #only removing sites flagged as 'mask'
		
		# Remove sites listed as mask
		homopl_df <- dplyr::anti_join(homopl_df, probl_sites, by=c("defining_mut"="aa_prot_site"))
		homopl_df
	}
	
	# Plot tree with defining mutations
	.cluster_tree <- function(tips) {
		library(phangorn)
		
		tr <- keep.tip(tre, tips)
		ggtr <- ggtree(tr)
		mut_list <- strsplit( amd[match(tips, amd$sequence_name),]$mutations, split="\\|" )
		tab_mut <- table(do.call(c, mut_list))
		tab_mut <- tab_mut / Ntip(tr)
		shared_muts <- names(tab_mut)[tab_mut >= defining_mut_threshold]
		segregating_muts <- lapply(mut_list, function(x) setdiff(x, shared_muts) )
		all_segregating <- Reduce(union, segregating_muts)
		
		# remove stops
		all_segregating <- all_segregating[!grepl(all_segregating, pattern="[*]$")]
		all_segregating <- all_segregating[!grepl(all_segregating, pattern=":[*]")]
		annots <- rep("", Ntip(tr) + Nnode(tr))
		
		if(length(all_segregating) > 0) {
			allseg1 <- substr( regmatches( all_segregating, regexpr(all_segregating, pattern=':[A-Z]' )), 2, 2)
			allseg2 <- regmatches( all_segregating, regexpr(all_segregating, pattern='[A-Z*]$' ))
			sites_post <- regmatches( all_segregating, regexpr(all_segregating, pattern=':.*$' ))
			sites_post <- substr( sites_post, 3, nchar( sites_post)-1 )
			sites_pre <- regmatches( all_segregating, regexpr(all_segregating, pattern='^.*:' ))
			sites <- paste0(sites_pre, sites_post)
			
			aas <- c()
			aas_list = lapply( seq_along( all_segregating ), function(i){
				do.call( rbind, lapply( mut_list,	function(x) ifelse(all_segregating[i] %in% x, allseg1[i], allseg2[i])))
			})
			
			aas <- do.call(cbind, aas_list)
			colnames(aas) <- all_segregating
			rownames(aas) <- tr$tip.label
			aas <- aas[,order(sites)]
			aas[is.na(aas)] <- "X"
			sites <- sort(sites)
			# aas1 <- as.AAbin(aas)
			# print(aas1)
			# print(class(aas1))
			aas2 <- as.phyDat(aas)
			#print(aas2)
			ap <- ancestral.pars(tr, aas2, return="phyDat")
			ap1 <- as.character(ap)
			for(e in postorder(tr)) {
				n <- tr$edge[e, 1]
				tp <- tr$edge[e, 2]
				j <- which(ap1[n,] != ap1[tp,])
				# keep only S and N annots
				j <- j [ grepl(sites[j], pattern ='^[SN]:') ]
				annots[[tp]] <- paste( paste0(sites[j], ap1[tp,j]), collapse="," )
				if(nchar(annots[tp]) == 0) annots[tp] <- NA
			}
		}
		
		node_df <- data.frame(node=1:(Ntip(tr) + Nnode(tr)), annot=annots, stringsAsFactors=FALSE)
		ggtr1 <- ggtr %<+% node_df
		ggtr1 <- ggtr1 + geom_label( aes(x = branch, label = annot, size = 5)) #+ geom_tiplab(align=TRUE)
		ggtr2 <- ggtr1
		if(length(all_segregating) < 100) {
			suppressMessages( ggtr2 <- gheatmap(ggtr1, as.data.frame(aas), width=0.66, offset=0.0005, colnames=FALSE, colnames_angle=-90, colnames_position="top", colnames_offset_y=-2) + theme(legend.position="none") )
		}
		ggtr2
	}
	
	# summarise lineage frequencies inside the node of interest
	.lineage_summary <- function(tips, max_rows=5) {
		if(is.null(tips)) {
			return("")
		}
		lineages <- amd$lineage[ match(tips, amd$sequence_name) ]
		lins_table <- sort(table(lineages), decreasing=TRUE) / length(lineages)
		if(length(lins_table) > 1) {
			perc_df <- as.data.frame(lins_table)
			colnames(c("Lineage","Freq"))
		} else {
			perc_df <- data.frame(Lineage=lineages[1], Freq=1)
		}
		perc_df <- perc_df[ 1:min(nrow(perc_df),max_rows), ]
		perc_df$Freq <- paste0(round(perc_df$Freq*100), "%")
		perc_df
	}
	
	.region_summary <- function(tips, max_rows=5) {
		if(is.null(tips)) {
			return("")
		}
		regions <- amd$region[ match(tips, amd$sequence_name) ]
		regs_table <- sort(table(regions), decreasing=TRUE) / length(regions)
		if(length(regs_table) > 1) {
			perc_df <- as.data.frame(regs_table)
			colnames(perc_df) <- c("Region", "Freq")
		}else {
			perc_df <- data.frame(Region=regions[1], Freq=1)
		}
		perc_df <- perc_df[ 1:min(nrow(perc_df), max_rows), ]
		perc_df$Freq <- paste0(round(perc_df$Freq*100), "%")
		perc_df
	}
	
	.plot_global_mut_freqs_worldwide <- function(homopl_df) {
		
		homopl_df$defining_mut <- toupper(homopl_df$defining_mut)
		homopl_df$protein <- sub("\\:.*", "", homopl_df$defining_mut)
		homopl_df <- homopl_df[(homopl_df$protein != "SYNSNP"),]
		homopl_df <- homopl_df[(homopl_df$protein != "ORF1AB"),]
		#View(homopl_df)
		
		for(i in 1:nrow(homopl_df)) {
			mut_freq <- getPrevalence(mutations = homopl_df[i,2], logInfo=FALSE)
			pl <- plotPrevalenceOverTime(mut_freq, title = glue("Prevalence of {homopl_df[i,2]} worldwide"))
			ggsave(pl, file=glue("{output_dir}/global_mut_freqs/{homopl_df[i,2]}_prevalence.jpg"), width=6, height=5)
		}
	}
	
	# compute node stats based on conditions
	tgt_nodes <- which((ndesc >= min_descendants) & (ndesc <= max_descendants) & (clade_age >= min_cluster_age_yrs))
	message(glue("Number of target nodes: {length(tgt_nodes)}"))
	message("Target nodes are:")
	print(tgt_nodes)
	homoplasies <- .find_homoplasies(tgt_nodes) # will consider homoplasies only for nodes passing number of descendants and clade age filters
	st0 <- st1 <- st2 <- st3 <- st4 <- comb_stats <- list()
	for(n in 1:(ntip+nnode)) {
		if(n <= (ntip+1)) {
			st0[[n]] <- st1[[n]] <- st2[[n]] <- st3[[n]] <- st4[[n]] <- -1
		}else {
			if(n %% 1000 == 0) message(glue("Progress: {n} / {(ntip+nnode)}"))
			sisters <- .get_comparison_sister_node(n)[[1]]
			if(!is.null(sisters)) {
				if(n %in% tgt_nodes) {
					st0[[n]] <- .node_muts(n)
					st1[[n]] <- .ratio_sizes_stat(n)
					st2[[n]] <- .ratio_persist_time_stat(n)
					st3[[n]] <- .logistic_growth_stat(n)
					st4[[n]] <- ""
					if(n %in% homoplasies$node) {
						st4[[n]] <- "Yes"
					}else{
						st4[[n]] <- "No"
					}
					comb_stats[[n]] <- cbind(n, paste(st0[[n]]$defining, collapse="|"), sisters[1], sisters[2], st1[[n]][1], st1[[n]][2], st1[[n]][3], st2[[n]][1], st2[[n]][2], st2[[n]][3], st2[[n]][4], st3[[n]][,1], st3[[n]][,2], st4[[n]]) #st3[[n]][,3]
				}else {
					comb_stats[[n]] <- cbind(n, "", sisters[1], sisters[2], -1,-1,-1,-1,-1,-1,-1,-1,-1, "No")
					comb_stats[[n]] <- as.data.frame(comb_stats[[n]])
				}
				rownames(comb_stats[[n]]) <- NULL
			}
		}
	}
	
	stats_df <- as.data.frame(do.call(rbind, comb_stats))
	colnames(stats_df) <- c("node","defining_muts","comp_sister1","comp_sister2","size1","size2","ratio_sizes","min_time_node","max_time1","max_time2","ratio_persist_time","logistic_growth", "logistic_growth_p", "homoplasies") #"diff_time_sister1","diff_time_sister2", "control_clade(s)"
	
	# Filter and write to CSV files
	# unfiltered output (just removing filtered based on overall min_descendants => -1 on all values)
	stats_df_unfilt <- stats_df[(stats_df$ratio_sizes > 1),]
	# stats_df_unfilt <- as.data.table(stats_df_unfilt)
	stats_df_unfilt <- stats_df_unfilt[order(as.numeric(stats_df_unfilt$ratio_sizes), as.numeric(stats_df_unfilt$ratio_persist_time), as.numeric(stats_df_unfilt$logistic_growth), decreasing=TRUE),]
	write.csv(stats_df_unfilt, file=glue('{output_dir}/stats_unfiltered.csv'), quote=FALSE, row.names=FALSE)
	
	# ratio sizes threshold
	stats_df_rs_all <- .extract_value_below_quantile_threshold(stats_df_unfilt, stats_df_unfilt$ratio_sizes, quantile_choice, quantile_threshold_ratio_sizes)
	stats_df_rs <- stats_df_rs_all[, c(1:7,14)]
	write.csv(stats_df_rs, file=glue('{output_dir}/stats_ratio_size_threshold.csv'), quote=FALSE, row.names=FALSE)
	
	# ratio persistence time threshold
	stats_df_rpt_all <- .extract_value_below_quantile_threshold(stats_df_unfilt, stats_df_unfilt$ratio_persist_time, quantile_choice, quantile_threshold_ratio_persist_time)
	stats_df_rpt_all <- stats_df_rpt_all[order(as.numeric(stats_df_rpt_all$ratio_persist_time), decreasing=TRUE),]
	stats_df_rpt <- stats_df_rpt_all[, c(1:4,8:11,14)]
	write.csv(stats_df_rpt, file=glue('{output_dir}/stats_ratio_persist_time_threshold.csv'), quote=FALSE, row.names=FALSE)
	
	# logistic growth p-filtered
	stats_df_lg_p <- stats_df_unfilt[(stats_df_unfilt$logistic_growth_p > -1) & (stats_df_unfilt$logistic_growth_p <= 0.05),]
	stats_df_lg_p_all <- stats_df_lg_p[order(as.numeric(stats_df_lg_p$logistic_growth), decreasing=TRUE),]
	stats_df_lg_p <- stats_df_lg_p_all[, c(1,12:13,14)] #1:3,13:14
	write.csv(stats_df_lg_p, file=glue('{output_dir}/stats_logit_growth_p_threshold.csv'), quote=FALSE, row.names=FALSE)
	
	# Get intersection between applied filters
	stats_df_intersect <- merge(stats_df_rs, stats_df_rpt) # by="node"
	stats_df_intersect <- merge(stats_df_intersect, stats_df_lg_p)
	stats_df_intersect <- stats_df_intersect[order(as.numeric(stats_df_intersect$ratio_sizes), as.numeric(stats_df_intersect$ratio_persist_time), as.numeric(stats_df_intersect$logistic_growth), decreasing=TRUE),]
	write.csv(stats_df_intersect, file=glue('{output_dir}/stats_intersection.csv'), quote=FALSE, row.names=FALSE)
	
	# Get union between applied filters to plot trees starting on node
	stats_df_union <- rbind(stats_df_rs_all,stats_df_rpt_all,stats_df_lg_p_all)
	stats_df_union <- unique(stats_df_union); rownames(stats_df_union) <- NULL
	stats_df_union <- stats_df_union[order(as.numeric(stats_df_union$ratio_sizes), as.numeric(stats_df_union$ratio_persist_time), as.numeric(stats_df_union$logistic_growth), decreasing=TRUE),]
	write.csv(stats_df_union, file=glue('{output_dir}/stats_union.csv'), quote=FALSE, row.names=FALSE)
	
	# Homoplasies (all target nodes considered)
	homoplasies1_all_tgt_nodes <- data.frame(tgt_nodes); colnames(homoplasies1_all_tgt_nodes) <- c("node")
	homoplasies1_all_tgt_nodes_df <- merge(homoplasies1_all_tgt_nodes, homoplasies, by="node", all.x=TRUE, all.y=FALSE)
	homoplasies1_all_tgt_nodes_df <- na.omit(homoplasies1_all_tgt_nodes_df)
	homoplasies1_all_tgt_nodes_df <- .remove_known_problematic_sites(homoplasies1_all_tgt_nodes_df)
	stats_df_homopl_nodes1 <- homoplasies1_all_tgt_nodes_df %>% group_by(defining_mut) %>% summarize(nodes_homopl = paste0(na.omit(node), collapse="|")) %>% ungroup()
	stats_df_homopl1 <- merge(homoplasies1_all_tgt_nodes_df, stats_df_homopl_nodes1, by="defining_mut")
	stats_df_homopl1 <- stats_df_homopl1[!duplicated(stats_df_homopl1$defining_mut),]
	stats_df_homopl1 <- .annotate_s_homoplasies(stats_df_homopl1)
	stats_df_homopl1 <- stats_df_homopl1 %>% select(nodes_homopl, defining_mut, Freq_homopl, s_mut,	s_mut_region_interest)
	stats_df_homopl1[is.na(stats_df_homopl1)] = ""
	stats_df_homopl1 <- stats_df_homopl1[order(stats_df_homopl1$s_mut, stats_df_homopl1$s_mut_region_interest, decreasing=TRUE),]
	#View(stats_df_homopl1)
	stats_df_aadns1 <- .annotate_diff_aa_mut_same_site(homoplasies1_all_tgt_nodes_df)
	stats_df_aamw1 <- .annotate_adjacent_muts_window_s3(homoplasies1_all_tgt_nodes_df)
	#View(stats_df_aamw1)
	
	# Plot global prevalence of mutation for the homoplasies detected
	if(plot_global_mut_freq) {
		
		if(!dir.exists(glue("{output_dir}/global_mut_freqs/")))
			suppressWarnings( dir.create(glue("{output_dir}/global_mut_freqs/")) )
		
		.plot_global_mut_freqs_worldwide(stats_df_homopl1)
		#.plot_global_mut_freqs_worldwide(stats_df_aadns1)
		.plot_global_mut_freqs_worldwide(stats_df_aamw1)
	}
	
	# Homoplasies (node detected by stat)
	stats_df_union$node <- as.integer(stats_df_union$node)
	homoplasies2_detect_df <- stats_df_union	%>% left_join(homoplasies, by="node")
	homoplasies2_detect_df <- na.omit(homoplasies2_detect_df)
	homoplasies2_detect_df <- .remove_known_problematic_sites(homoplasies2_detect_df)
	stats_df_homopl_nodes2 <- homoplasies2_detect_df %>% group_by(defining_mut) %>% summarize(nodes_homopl = paste0(na.omit(node), collapse="|")) %>% ungroup()
	stats_df_homopl2 <- merge(homoplasies2_detect_df, stats_df_homopl_nodes2, by="defining_mut")
	stats_df_homopl2$Freq_homopl <- NULL
	stats_df_homopl2_freq <- setDT(stats_df_homopl2)[, .(Freq_homopl = .N), by = .(defining_mut)]
	stats_df_homopl2_freq_df <- merge(stats_df_homopl2, stats_df_homopl2_freq, by="defining_mut")
	stats_df_homopl2_freq_df <- stats_df_homopl2_freq_df[as.numeric(stats_df_homopl2_freq_df$Freq_homopl) > 1,]
	stats_df_homopl2_freq_df <- .annotate_s_homoplasies(stats_df_homopl2_freq_df)
	stats_df_homopl2_freq_df <- stats_df_homopl2_freq_df[, c(2:3,1,15:17,18,20,4:14)]
	stats_df_homopl2_freq_df[is.na(stats_df_homopl2_freq_df)] = ""
	stats_df_homopl2_freq_df <- stats_df_homopl2_freq_df[order(stats_df_homopl2_freq_df$s_mut, stats_df_homopl2_freq_df$s_mut_region_interest, decreasing=TRUE),]
	stats_df_aadns2 <- .annotate_diff_aa_mut_same_site(homoplasies2_detect_df)
	stats_df_aamw2 <- .annotate_adjacent_muts_window_s3(homoplasies2_detect_df)
	#View(stats_df_aamw2)
	
	# Homoplasies (node NOT detected by stat)
	homoplasies3_not_detect <- setdiff(tgt_nodes, stats_df_union$node)
	homoplasies3_not_detect <- data.frame(homoplasies3_not_detect); colnames(homoplasies3_not_detect) <- c("node")
	homoplasies3_not_detect_df <- homoplasies3_not_detect	%>% left_join(homoplasies, by="node")
	homoplasies3_not_detect_df <- na.omit(homoplasies3_not_detect_df)
	homoplasies3_not_detect_df <- .remove_known_problematic_sites(homoplasies3_not_detect_df)
	stats_df_homopl_nodes3 <- homoplasies3_not_detect_df %>% group_by(defining_mut) %>% summarize(nodes_homopl = paste0(na.omit(node), collapse="|")) %>% ungroup()
	stats_df_homopl3 <- merge(homoplasies3_not_detect_df, stats_df_homopl_nodes3, by="defining_mut")
	stats_df_homopl3$Freq_homopl <- NULL
	stats_df_homopl3_freq <- setDT(stats_df_homopl3)[, .(Freq_homopl = .N), by = .(defining_mut)]
	stats_df_homopl3_freq_df <- merge(stats_df_homopl3, stats_df_homopl3_freq, by="defining_mut")
	stats_df_homopl3_freq_df <- stats_df_homopl3_freq_df[as.numeric(stats_df_homopl3_freq_df$Freq_homopl) > 1,]
	stats_df_homopl3_freq_df <- stats_df_homopl3_freq_df[!duplicated(stats_df_homopl3_freq_df$defining_mut),]
	stats_df_homopl3_freq_df <- .annotate_s_homoplasies(stats_df_homopl3_freq_df)
	stats_df_homopl3_freq_df <- stats_df_homopl3_freq_df %>% select(nodes_homopl, defining_mut, Freq_homopl, s_mut, s_mut_region_interest)
	stats_df_homopl3_freq_df[is.na(stats_df_homopl3_freq_df)] = ""
	stats_df_homopl3_freq_df <- stats_df_homopl3_freq_df[order(stats_df_homopl3_freq_df$s_mut, stats_df_homopl3_freq_df$s_mut_region_interest, decreasing=TRUE),]
	stats_df_aadns3 <- .annotate_diff_aa_mut_same_site(homoplasies3_not_detect_df)
	stats_df_aamw3 <- .annotate_adjacent_muts_window_s3(homoplasies3_not_detect_df)
	#View(stats_df_aamw3)
	
	write.csv(stats_df_homopl1, file=glue('{output_dir}/homoplasy_DETAILS_all_tgt_nodes.csv'), quote=FALSE, row.names=FALSE)
	write.csv(stats_df_aadns1, file=glue('{output_dir}/homoplasy_SAME_SITE_muts_all_tgt_nodes.csv'), quote=FALSE, row.names=FALSE)
	write.csv(stats_df_aamw1, file=glue('{output_dir}/homoplasy_NEIGHBOUR_WINDOW_muts_all_tgt_nodes.csv'), quote=FALSE, row.names=FALSE)
	
	write.csv(stats_df_homopl2_freq_df, file=glue('{output_dir}/homoplasy_DETAILS_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	write.csv(stats_df_aadns2, file=glue('{output_dir}/homoplasy_SAME_SITE_muts_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	write.csv(stats_df_aamw2, file=glue('{output_dir}/homoplasy_NEIGHBOUR_WINDOW_muts_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	
	write.csv(stats_df_homopl3_freq_df, file=glue('{output_dir}/homoplasy_DETAILS_NOT_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	write.csv(stats_df_aadns3, file=glue('{output_dir}/homoplasy_SAME_SITE_muts_NOT_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	write.csv(stats_df_aamw3, file=glue('{output_dir}/homoplasy_NEIGHBOUR_WINDOW_muts_NOT_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	
	message(glue("Nodes matching threshold {ifelse(threshold_keep_lower, '<', '>')} {quantile_threshold_ratio_sizes} quantile of ratio sizes: {nrow(stats_df_rs)}"))
	message(glue("Nodes matching threshold {ifelse(threshold_keep_lower, '<', '>')} {quantile_threshold_ratio_persist_time} quantile of ratio persistence time: {nrow(stats_df_rpt)}"))
	message(glue("Nodes matching threshold <= 0.05 of logistic regression p-value: {nrow(stats_df_lg_p)}"))
	message(glue("Nodes matching all thresholds: {nrow(stats_df_intersect)}"))
	message(glue("Nodes matching at least one threshold: {nrow(stats_df_union)}"))
	
	if(!dir.exists(glue("{output_dir}/node_specific/")))
		suppressWarnings( dir.create(glue("{output_dir}/node_specific/")) )
	
	message(glue("Plotting trees + sequences + lineage + region summaries for nodes matching at least one of the thresholds..."))
	# Only plotting annotated trees, lineage and region summaries for nodes passing at least 1 of 3 criteria
	for(i in 1:nrow(stats_df_union)) { #length(tgt_nodes)
		comp_res <- .get_comparison_sister_node(as.numeric(stats_df_union[i,1]))
		tips <- c(comp_res[[5]], comp_res[[6]])
		if(compute_tree_annots) {
			ggtr <- .cluster_tree(tips)
			suppressMessages(ggsave(ggtr, file=glue("{output_dir}/node_specific/{stats_df_union[i,1]}_clust_tree.pdf"),height=max(6, floor(length(tips) / 5 )), width = min(44, max(24, sqrt(length(tips)))), limitsize = FALSE  ))
		}
		# else {
		# 	ggtr <- ggtree(tre, mrsd=max_date) + scale_x_continuous(expand = c(0, 0)) + theme_tree2(axis.text.x=element_text(size=5, angle=45, vjust=1, hjust=1))
		# 	ggtr_p1 <- zoomClade(ggtr, node=stats_df_union[i,1])
		# 	suppressMessages( ggsave(file=glue('{output_dir}/node_specific/{stats_df_union[i,1]}_zoomed_tree.pdf'), plot=ggtr_p1, dpi=600, limitsize=FALSE) )
		# }
		tips_df <- amd[ amd$sequence_name %in% tips, ]
		write.csv(tips_df, file=glue('{output_dir}/node_specific/{stats_df_union[i,1]}_sequences.csv'), quote=FALSE, row.names=FALSE)
		lineage_summary <- .lineage_summary(tips=tips)
		write.csv(lineage_summary, file=glue('{output_dir}/node_specific/{stats_df_union[i,1]}_lineage_summary.csv'), quote=FALSE, row.names=FALSE)
		region_summary <- .region_summary(tips=tips)
		write.csv(region_summary, file=glue('{output_dir}/node_specific/{stats_df_union[i,1]}_region_summary.csv'), quote=FALSE, row.names=FALSE)
	}
	
	message(paste0("Number of descendants vector length: ", length(ndesc)))
	
	return(stats_df)
}

# 2019-12-30 to 2020-06-30 (first lineages) + 2020-07-01 to 2020-12-31 (B.1.177 + start Alpha)
start <- Sys.time()
cladeScore1 <- cladeScore(sc2_tre, sc2_md, min_descendants=10, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=as.Date("2019-12-30"), max_date=as.Date("2020-12-31"), branch_length_unit="days", output_dir="results/01_sc2_root_to_dec2020_removing_artifacts", quantile_choice=1/100, quantile_threshold_ratio_sizes="1%", quantile_threshold_ratio_persist_time="1%", threshold_keep_lower=TRUE, defining_mut_threshold=0.75, root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995, compute_tree_annots=FALSE, plot_global_mut_freq=FALSE) #, ncpu=6
end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time,"mins")) # 12 mins (min_desc=10, all 2020) 

# 2021-01-01 to 2021-05-31 (Alpha + Delta rapidly replacing)
start <- Sys.time()
cladeScore2 <- cladeScore(sc2_tre, sc2_md, min_descendants=100, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=as.Date("2021-01-01"),max_date=as.Date("2021-05-31"),branch_length_unit="days", output_dir="results/02_sc2_jan2021_to_may2021_removing_artifacts", quantile_choice=1/100, quantile_threshold_ratio_sizes="1%", quantile_threshold_ratio_persist_time="1%", threshold_keep_lower=TRUE,  defining_mut_threshold=0.75, root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995,compute_tree_annots=FALSE, plot_global_mut_freq=FALSE)
end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time,"mins")) # 80 mins

# 2021-06-01 to 2021-12-31 (Delta + Omicron BA.1 rapidly replacing)
start <- Sys.time()
cladeScore3 <- cladeScore(sc2_tre, sc2_md, min_descendants=100, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=as.Date("2021-06-01"),max_date=as.Date("2021-12-31"),branch_length_unit="days", output_dir="results/03_sc2_jun2021_to_dec2021", quantile_choice=1/100, quantile_threshold_ratio_sizes="1%", quantile_threshold_ratio_persist_time="1%", threshold_keep_lower=TRUE,  defining_mut_threshold=0.75, root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995,compute_tree_annots=FALSE, plot_global_mut_freq=TRUE)
end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time,"mins")) # 1339 mins

# 2022-01-01 to 2022-04-30 (Omicron BA.1 + BA.2 rapidly replacing)
start <- Sys.time()
cladeScore4 <- cladeScore(sc2_tre, sc2_md, min_descendants=100, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=as.Date("2022-01-01"),max_date=as.Date("2022-04-30"),branch_length_unit="days", output_dir="results/04_sc2_jan2022_to_apr2022", quantile_choice=1/100, quantile_threshold_ratio_sizes="1%", quantile_threshold_ratio_persist_time="1%", threshold_keep_lower=TRUE,  defining_mut_threshold=0.75, root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995,compute_tree_annots=FALSE, plot_global_mut_freq=TRUE)
end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time,"mins"))

# Whole tree period (2019-12-30 to 2022-04-30)
start <- Sys.time()
cladeScore_all <- cladeScore(sc2_tre, sc2_md, min_descendants=10, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=as.Date("2019-12-30"), max_date=as.Date("2022-04-30"),branch_length_unit="days", output_dir="results/05_sc2_whole_period", quantile_choice=1/100, quantile_threshold_ratio_sizes="1%", quantile_threshold_ratio_persist_time="1%", threshold_keep_lower=TRUE, defining_mut_threshold=0.75, root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995,compute_tree_annots=FALSE, plot_global_mut_freq=FALSE)
end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time,"mins"))

saveRDS(cladeScore1, "rds/results/period2019-2020.rds")
saveRDS(cladeScore2, "rds/results/periodjan2021-may2021.rds")
saveRDS(cladeScore3, "rds/results/periodjun2021-dec2021.rds")
saveRDS(cladeScore4, "rds/results/periodjan2022-apr2022.rds")
saveRDS(cladeScore_all, "rds/results/whole_period.rds")
