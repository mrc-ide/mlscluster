library(ape)
library(lubridate)
library(tidyr)
library(dplyr)
library(glue)
library(ggtree)
library(ggplot2)

tre2 = readRDS("rds/timetree_ebola.rds")
amd2 = readRDS("rds/md_ebola.rds")

# load("rds/tree_cleaned.RData")
# tre <- eng_tre_bin # not ideal because does not include root from China and chronumental estimate of root goes to end of 2017
# md_eng_seqnames_match$lineage <- NULL
# amd <- md_eng_seqnames_match
# rm(eng_tre, eng_tre_bin, md_eng_seqnames_match, min_blen)

chronum_tre <- readRDS("rds/tree_cog_chronumental.rds") # OK, using Wuhan sequence and getting min coal time on 2019.936
chronum_amd <- readRDS("rds/amd_cog_chronumental.rds")

cladeScore <- function(tre, amd, min_descendants=100, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=NULL, max_date=NULL,
																							output_dir=paste0("clade_score-",Sys.Date()), threshold_ratio_sizes=3,threshold_ratio_persist_time=3, gen_time=7/365) { 
																																																																							# upper_threshold_sizes_logit_regr=5, ncpu=1
	
	stopifnot(is.rooted(tre))
	
	class(tre) <- "phylo"
	
	if (!dir.exists(output_dir))
		dir.create(output_dir)
	
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
	
	#amd$sts <- amd$sample_time
	# remove missing dates and filter by sample times
	amd <- amd[!is.na(amd$sample_time),]
	amd <- amd[(amd$sample_time >= min_time) & (amd$sample_time <= max_time),] 
	
	# named vector of sequence_names and sample_times
	sts <- setNames(amd$sample_time , amd$sequence_name)
	amd <- amd[amd$sequence_name %in% tre$tip.label,]
	message(paste0("Number of sequences included: ", nrow(amd)))
	
	dropped_tips <- setdiff(tre$tip.label, amd$sequence_name)
	message(paste0("Dropped tips: ", length(dropped_tips)))
	
	tre <- keep.tip(tre, intersect(tre$tip.label, amd$sequence_name))
	message(paste0("Tips in the tree: ", length(tre$tip.label)))
	
	# variables and data structures to quickly look up tree data
	# based on treestructure/tfpscanner code
	message("Computing node descendants and times...")
	ntip <- Ntip(tre)
	nnode <- Nnode(tre)
	
	# getting 'coalescent' times
	nde <- node.depth.edgelength(tre) # depth of a node using branch lengths
	message("Difference between max and min sample times before applying filters")
	message(diff_stimes)
	message("node.depth.edgelength max value")
	message(max(nde))
	# deal with chronumental trees in day scale
	if(max(nde) > (10 * diff_stimes)) { # since root time can be far in the past multiply by 10
		nde <- node.depth.edgelength(tre) / (29903/365)
		message("node.depth.edgelength max value after adjusting")
		message(max(nde))
	}
	rh <- max(nde[1:ntip]) # root heights
	stimes <- nde[1:ntip]
	# max_height <- Inf; max_height <- min(rh, max_height)
	shs <- rh - stimes # time to most recent sample
	#inhs <- sort(rh - nde[ (ntip+1):(ntip+nnode) ])
	nhs = nodeheights <- rh - nde # node heights
	mrst <- max(sts) # most recent sample time
	tmrca <- mrst - nhs # decimal times of coalescent events
	
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
	
	for (e in 1:nrow(poedges)){
		n <- poedges[e,1] # postorder traversal of nodes column
		tp <- poedges[e,2] # postorder traversal of tips column
		ndesc[n] <- ndesc[n] + ndesc[tp]
		Ndesc[n] <- Ndesc[n] + Ndesc[tp]
		max_desc_time[n] <- max(max_desc_time[n], max_desc_time[tp])
		min_desc_time[n] <- min(min_desc_time[n], min_desc_time[tp])
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
	
	# compute node stats based on conditions
	# TODO: consider these filters and wrap inside function(s)
	tgt_nodes <- which((ndesc >= min_descendants) & (ndesc <= max_descendants) & (clade_age >= min_cluster_age_yrs))
	message(paste0("Number of target nodes: ", length(tgt_nodes)))
	#tgt_node_ancestors <- do.call(c, lapply(tgt_nodes, function(tp) ancestors[[tp]]))
	
	# descendant tip indices
	descendant_tips <- descendants
	# daughter nodes: direct descendants of a given node
	daughter_nodes <- descendants
	sisters <- daughter_nodes
	length_sisters <- lapply(1:(ntip+nnode), function(tp) rep(1,2))
	ratio_sizes <- lapply(1:(ntip+nnode), function(tp) rep(1,1))
	min_persist_time <- lapply(1:(ntip+nnode), function(tp) rep(min_time_full,2))
	max_persist_time <- lapply(1:(ntip+nnode), function(tp) rep(max_time_full,2))
	diff_persist_time <- lapply(1:(ntip+nnode), function(tp) rep(0,2))
	ratio_persist_time <- lapply(1:(ntip+nnode), function(tp) rep(1,1))
	
	for(n in (ntip+1):(ntip+nnode)) {
		descendant_tips[[n]] <- descendant_tips[[n]][ descendant_tips[[n]] <= ntip ]
		daughter_nodes[[n]] <- daughter_nodes[[n]][ daughter_nodes[[n]] > ntip ]
		sisters[[n]] <- tre$edge[ tre$edge[,1] == n, 2 ]
		
		for(ll in 1:length(sisters[[n]])) {
			# size of sister clades
			length_sisters[[n]][ll] <- ndesc[ sisters[[n]][ll] ]
			# estimated tmrca of comparison node
			min_persist_time[[n]][ll] <- tmrca[n]
			desc_same_sisters <- descendants[[ sisters[[n]][ll] ]]
			# max estimated time of descending tips
			max_persist_time[[n]][ll] <- max(tmrca[ desc_same_sisters ])
			diff_persist_time[[n]][ll] <- max_persist_time[[n]][ll] - min_persist_time[[n]][ll]
		}
		
		# ratios are inflated for some nodes because sisters can have only one tip descending while the other has lots of internal nodes and resulting tips
		length_sisters[[n]] <- sort(length_sisters[[n]], decreasing=TRUE)
		ratio_sizes[[n]][1] <- round( length_sisters[[n]][1] / length_sisters[[n]][2], digits=2)
		diff_persist_time[[n]] <- sort(diff_persist_time[[n]], decreasing=TRUE)
		ratio_persist_time[[n]] <- round( diff_persist_time[[n]][1] / diff_persist_time[[n]][2], digits=2)
	}
	
	# Filter and write ratio of sizes to file
	ratio_sizes_df <- as.data.frame(do.call(rbind, ratio_sizes))
	ratio_sizes_df$node <- rownames(ratio_sizes_df); rownames(ratio_sizes_df) <- NULL
	colnames(ratio_sizes_df) <- c("ratio_sizes","node")
	ratio_sizes_df <- ratio_sizes_df[(ratio_sizes_df$ratio_sizes >= threshold_ratio_sizes),]
	ratio_sizes_df <- ratio_sizes_df[order(-ratio_sizes_df$ratio_sizes),]
	ratio_sizes_df <- subset(ratio_sizes_df, select=c(2,1))
	write.csv(ratio_sizes_df, file=glue('{output_dir}/ratio_sizes.csv'), quote=FALSE, row.names=FALSE)
	
	# Filter and write ratio of persistence time to file
	ratio_persist_time_df <- as.data.frame(do.call(rbind, ratio_persist_time))
	ratio_persist_time_df$node <- rownames(ratio_persist_time_df); rownames(ratio_persist_time_df) <- NULL
	colnames(ratio_persist_time_df) <- c("ratio_persist_time","node")
	ratio_persist_time_df <- ratio_persist_time_df[(ratio_persist_time_df$ratio_persist_time >= threshold_ratio_persist_time),]
	ratio_persist_time_df <- ratio_persist_time_df[order(-ratio_persist_time_df$ratio_persist_time),]
	ratio_persist_time_df <- subset(ratio_persist_time_df, select=c(2,1))
	write.csv(ratio_persist_time_df, file=glue('{output_dir}/ratio_persistence_time.csv'), quote=FALSE, row.names=FALSE)
	
	# Get intersection between ratio sizes and persistence time to plot trees starting on node
	nodes_intersect <- merge(ratio_sizes_df, ratio_persist_time_df, by="node")
	message(glue("Nodes matching threshold = {threshold_ratio_sizes} of ratio sizes: {nrow(ratio_sizes_df)}"))
	message(glue("Nodes matching threshold = {threshold_ratio_persist_time} of ratio persistence time: {nrow(ratio_persist_time_df)}"))
	message(glue("Nodes matching both thresholds: {nrow(nodes_intersect)}"))
	write.csv(nodes_intersect, file=glue('{output_dir}/intersect_ratios_threshold.csv'), quote=FALSE, row.names=FALSE)
	
	# tip labels
	descendant_ids <- lapply(1:(ntip+nnode), function(tp) {
		na.omit(tre$tip.label[descendant_tips[[tp]]])
	})
	
	desc_tips_same_sisters1 <- list()
	desc_tips_same_sisters2 <- list()
	lg_res_list <- list()
	for(n in (ntip+1):(ntip+nnode)) {
		desc_tips_same_sisters1[[n]] <- rep( 1, ndesc[ sisters[[n]][1] ]) #descendant_tips[[ sisters[[n]][ll] ]]
		desc_tips_same_sisters1[[n]] <- descendant_ids[[ sisters[[n]][1] ]]
		desc_tips_same_sisters2[[n]] <- rep( 1, ndesc[ sisters[[n]][2] ] )
		desc_tips_same_sisters2[[n]] <- descendant_ids[[ sisters[[n]][2] ]]
		# TODO remove this when have incorporated ndesc filters?
		if((length(desc_tips_same_sisters1[[n]]) >= min_descendants) && (length(desc_tips_same_sisters2[[n]]) >= min_descendants)) {
			sts_tips_same_sisters1 <- sts[ desc_tips_same_sisters1[[n]] ]
			sts_tips_same_sisters2 <- sts[ desc_tips_same_sisters2[[n]] ]
			X <- data.frame(sequence_name=c(names(sts_tips_same_sisters1),names(sts_tips_same_sisters2)),time=c(sts_tips_same_sisters1,sts_tips_same_sisters2), type=c(rep("sister1", length(sts_tips_same_sisters1)), rep("sister2",length(sts_tips_same_sisters2))))
			X <- na.omit(X)
			model <- glm(type=="sister1" ~ time, data=X, family=binomial(link="logit"))
			s <- summary(model)
			rel_growth_gen_time <- unname(coef(model)[2] * gen_time)
			p <- s$coefficients[2,4]
			lg_res_list[[n]] <- cbind(n, sisters[[n]][1], sisters[[n]][2], rel_growth_gen_time, p)
		}
	}
	lg_res_df <- do.call(rbind, lg_res_list)
	lg_res_df <- as.data.frame(lg_res_df)
	lg_res_df <- lg_res_df[order(-lg_res_df$rel_growth_gen_time),]
	colnames(lg_res_df) <- c("node","sister1","sister2","rel_growth_gen_time","p")
	write.csv(lg_res_df, file=glue('{output_dir}/logistic_growth_rates.csv'), quote=FALSE, row.names=FALSE)
	
	if ( length(tre$tip.label) < 1e3 ){
		message(glue("Plotting trees zoomed on node matching thresholds..."))
		p <- ggtree(tre, mrsd=max_date, as.Date=TRUE) + scale_x_date(date_labels="%b\n%Y", date_breaks="4 months") + theme_tree2(axis.text.x=element_text(size=6))
		for(i in 1:nrow(nodes_intersect)) {
			p1 <- zoomClade(p, node=nodes_intersect[i,1])
			suppressMessages( ggsave(file=glue('{output_dir}/{nodes_intersect[i,1]}.pdf'), plot=p1, dpi=600) )
		}
	}
	
	# max_dist_to_tip <- rep(0, ntip+nnode)
	# for(e in 1:nrow(poedges)) {
	# 	n <- poedges[e,1]
	# 	tp <- poedges[e,2]
	# 	max_dist_to_tip[n] <- max(max_dist_to_tip[n], max_dist_to_tip[tp] + 1)
	# }
	
	message(paste0("Number of descendants vector length: ", length(ndesc)))
	
	#return(list(ratio_sizes_df, ratio_persist_time_df, nodes_intersect, descendants, descendant_tips, daughter_nodes, sisters, length_sisters, tmrca, ndesc))
	return(list(ratio_sizes_df, ratio_persist_time_df, nodes_intersect, lg_res_df, sisters, length_sisters, tmrca, max_persist_time))
}

test_ebola <- cladeScore(tre2, amd2, min_descendants=5, max_descendants=75, min_cluster_age_yrs=0.2/12, min_date=as.Date("2014-07-01"),
																									max_date=as.Date("2015-10-24"),output_dir="test_ebola", threshold_ratio_sizes=2, threshold_ratio_persist_time=1.5, gen_time=16.6)

test_chron <- cladeScore(chronum_tre, chronum_amd, min_descendants=25, max_descendants=500, min_cluster_age_yrs=1/12, min_date=as.Date("2019-12-01"),
																			max_date=as.Date("2020-12-31"), output_dir="test_new_ctre", threshold_ratio_sizes=3, threshold_ratio_persist_time=2, gen_time=7/365)

# TODO: run for these periods
# DONE: 2019-12-01 to 2020-06-30 (first lineages)
# 2020-07-01 to 2020-12-31 (B.1.177 + start Alpha); 
# 2021-01-01 to 2021-05-31 (Alpha + Delta rapidly replacing);
# 2021-06-01 to 2021-12-31 (Delta + Omicron BA.1 rapidly replacing)
# 2022-01-01 to 2022-04-30 (Omicron BA.1 + BA.2 rapidly replacing)
# Whole tree period (not feasible now due to mem alloc issues on matrix)
