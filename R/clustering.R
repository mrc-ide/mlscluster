ebov_tre = readRDS("rds/timetree_ebola.rds")
ebov_md = readRDS("rds/md_ebola.rds")

load("rds/chron_timetree.RData")
sc2_tre <- eng_tre_bin # not ideal because does not include root from China and chronumental estimate of root goes to end of 2017
md_eng_seqnames_match$lineage <- NULL
md_eng_seqnames_match$sample_time <- NULL
sc2_md <- md_eng_seqnames_match
rm(eng_tre, eng_tre_bin, md_eng_seqnames_match, min_blen)

cladeScore <- function(tre, amd, min_descendants=100, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=NULL, max_date=NULL
																							, output_dir=paste0("clade_score-",Sys.Date()), threshold_ratio_sizes=3,threshold_ratio_persist_time=3
																							, gen_time=7/365, root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995) {
	
	library(ape)
	library(lubridate)
	library(tidyr)
	library(dplyr)
	library(glue)
	library(ggtree)
	library(ggplot2)
	
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
	
	#amd$sts <- amd$sample_time
	# remove missing dates and filter by sample times
	amd <- amd[!is.na(amd$sample_time),]
	amd <- amd[(amd$sample_time >= min_time) & (amd$sample_time <= max_time),] 
	
	# named vector of sequence_names and sample_times
	sts <- setNames(amd$sample_time , amd$sequence_name)
	amd <- amd[amd$sequence_name %in% tre$tip.label,]
	
	if(!is.rooted(tre)) {
		if(!(root_on_tip %in% amd$sequence_name)) {
			stopifnot(root_on_tip %in% tre$tip.label)
			amd <- rbind(data.frame( sequence_name=root_on_tip, sample_time=root_on_tip_sample_time, sample_date=date_decimal(root_on_tip_sample_time) ))
		}
	}
	message(paste0("Number of sequences included: ", nrow(amd)))
	
	dropped_tips <- setdiff(tre$tip.label, amd$sequence_name)
	message(paste0("Dropped tips: ", length(dropped_tips)))
	message("Root sequence Wuhan present?")
	print("Wuhan/WH04/2020" %in% intersect(tre$tip.label, amd$sequence_name))
	
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
	nde <- node.depth.edgelength(tre) # depth of a node using branch lengths
	message("Difference between max and min sample times before applying filters")
	message(diff_stimes)
	message("node.depth.edgelength max value")
	message(max(nde))
	# deal with chronumental trees in day scale (not really working as it should)
	if(max(nde) > (10 * diff_stimes)) { # since root time can be far in the past multiply by 10
		nde <- node.depth.edgelength(tre) / 365 # (29903/(diff_stimes * 365)) hardcoded for SARS-CoV-2 genomic size
		message("node.depth.edgelength max value after adjusting")
		message(max(nde))
	}
	rh <- max(nde[1:ntip]) # root heights
	stimes <- nde[1:ntip]
	shs <- rh - stimes # time to most recent sample
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
	
	# descendant tip indices
	descendant_tips <- descendants
	# daughter nodes: direct descendants of a given node
	daughter_nodes <- descendants
	sisters <- daughter_nodes
	
	for(n in (ntip+1):(ntip+nnode)) {
		descendant_tips[[n]] <- descendant_tips[[n]][ descendant_tips[[n]] <= ntip ]
		daughter_nodes[[n]] <- daughter_nodes[[n]][ daughter_nodes[[n]] > ntip ]
		# sisters[[n]] <- tre$edge[ tre$edge[,1] == n, 2 ] # not scaling well with larger tree
		#message(n)
	}
	
	# tip labels
	descendant_ids <- lapply(1:(ntip+nnode), function(tp) {
		na.omit(tre$tip.label[descendant_tips[[tp]]])
	})
	
	.get_sisters_node <- function(node) {
		sisters <- tre$edge[ tre$edge[,1] == node, 2 ]
		sisters
	}
	
	# ratio sizes function
	.ratio_sizes_stat <- function(node, sisters) {
		
		length_sisters <- rep(1,2)
		#ratio_sizes <- rep(1,1)
		for(ll in 1:length(sisters)) {
			if(ndesc[ sisters[ll] ] >= min_descendants) {
				# size of sister clades
				length_sisters[ll] <- ndesc[ sisters[ll] ]
			}else {
				length_sisters[ll] <- -1
			}
		}
		# TODO: way to calculate ratios considering polytomies
		length_sisters <- sort(length_sisters, decreasing=TRUE)
		ratio_sizes <- round( length_sisters[1] / length_sisters[2], digits=2)
		#print(ratio_sizes)
		
		return(ratio_sizes)
	}
	
	# ratio persistence time function
	.ratio_persist_time_stat <- function(node, sisters) {
		
		min_persist_time <- rep(min_time_full,2)
		max_persist_time <- rep(max_time_full,2)
		diff_persist_time <- rep(0,2)
		for(ll in 1:length(sisters)) {
			if(ndesc[ sisters[ll] ] >= min_descendants) {
				desc_same_sisters <- descendants[[ sisters[ll] ]]
				# estimated tmrca of comparison node
				min_persist_time[ll] <- tmrca[node]
				# max estimated time of descending tips
				max_persist_time[ll] <- max(tmrca[ desc_same_sisters ])
				diff_persist_time[ll] <- max_persist_time[ll] - min_persist_time[ll]
			}else {
				min_persist_time[ll] <- max_persist_time[ll] <- diff_persist_time[ll] <- -1
			}
		}
		# TODO: way to calculate ratios considering polytomies
		diff_persist_time <- sort(diff_persist_time, decreasing=TRUE)
		ratio_persist_time <- round( diff_persist_time[1] / diff_persist_time[2], digits=2)
		
		return(ratio_persist_time)
	}
	
	# logistic growth of sister clades function
	.logistic_growth_stat <- function(node, sisters) {
		# TODO: way to calculate logistic regression considering polytomies
		lg_res_list <- list()
		desc_tips_same_sisters1 <- rep( 1, ndesc[ sisters[1] ])
		desc_tips_same_sisters1 <- descendant_ids[[ sisters[1] ]]
		desc_tips_same_sisters2 <- rep( 1, ndesc[ sisters[2]] )
		desc_tips_same_sisters2 <- descendant_ids[[ sisters[2] ]]
		lg_res_list <- list()

		if((length(desc_tips_same_sisters1) >= min_descendants) && (length(desc_tips_same_sisters2) >= min_descendants)) {
			sts_tips_same_sisters1 <- sts[ desc_tips_same_sisters1 ]
			sts_tips_same_sisters2 <- sts[ desc_tips_same_sisters2 ]
			X <- data.frame(sequence_name=c(names(sts_tips_same_sisters1),names(sts_tips_same_sisters2)),time=c(sts_tips_same_sisters1,sts_tips_same_sisters2), type=c(rep("sister1", length(sts_tips_same_sisters1)), rep("sister2",length(sts_tips_same_sisters2))))
			X <- na.omit(X)
			model <- suppressWarnings ( glm(type=="sister1" ~ time, data=X, family=binomial(link="logit")) )
			s <- summary(model)
			rel_growth_gen_time <- unname(coef(model)[2] * gen_time)
			p <- s$coefficients[2,4]
			lg_res_list <- cbind(rel_growth_gen_time, p) #node, sisters[[node]][1], sisters[[node]][2], 
		}else {
			lg_res_list <- cbind(-1, -1)
		}
		return(lg_res_list)
	}
	
	# compute node stats based on conditions
	tgt_nodes <- which((ndesc >= min_descendants) & (ndesc <= max_descendants) & (clade_age >= min_cluster_age_yrs))
	message(glue("Number of target nodes: {length(tgt_nodes)}"))
	message("Target nodes are:")
	print(tgt_nodes)
	#tgt_node_ancestors <- do.call(c, lapply(tgt_nodes, function(tp) ancestors[[tp]]))
	st1 <- st2 <- st3 <- comb_stats <- list()
	for(n in 1:(ntip+nnode)) {
		if(n <= (ntip+1)) {
			st1[[n]] <- st2[[n]] <- st3[[n]] <- -1 
		}else {
			if(n %% 1000 == 0) message(glue("Progress: {n} / {(ntip+nnode)}"))
			sisters <- .get_sisters_node(n)
			if(n %in% tgt_nodes) {
				st1[[n]] <- .ratio_sizes_stat(n, sisters)
				st2[[n]] <- .ratio_persist_time_stat(n, sisters)
				st3[[n]] <- .logistic_growth_stat(n, sisters)
				comb_stats[[n]] <- cbind(n, sisters[1], sisters[2], st1[[n]], st2[[n]], st3[[n]][,1], st3[[n]][,2])
			}else {
				st1[[n]] <- st2[[n]] <- st3[[n]] <- -1
				comb_stats[[n]] <- cbind(n, sisters[1], sisters[2], -1, -1, -1, -1)
				comb_stats[[n]] <- as.data.frame(comb_stats[[n]])
			}
			rownames(comb_stats[[n]]) <- NULL
		}
	}
	stats_df <- as.data.frame(do.call(rbind, comb_stats))
	colnames(stats_df) <- c("node","sister1","sister2","ratio_sizes","ratio_persist_time","logistic_growth", "logistic_growth_p") 
	#View(stats_df)
	
	# Filter and write to CSV files
	# unfiltered output (just removing filtered based on overall min_descendants => -1 on all values)
	stats_df_unfilt <- stats_df[(stats_df$ratio_sizes > -1),]
	stats_df_unfilt <- stats_df_unfilt[order(-stats_df_unfilt$ratio_sizes, -stats_df_unfilt$ratio_persist_time, -stats_df_unfilt$logistic_growth, stats_df_unfilt$logistic_growth_p),]
	write.csv(stats_df_unfilt, file=glue('{output_dir}/stats_unfiltered.csv'), quote=FALSE, row.names=FALSE)
	
	# ratio sizes threshold
	stats_df_rs <- stats_df_unfilt[(stats_df_unfilt$ratio_sizes >= threshold_ratio_sizes),]
	write.csv(stats_df_rs, file=glue('{output_dir}/stats_ratio_size_threshold.csv'), quote=FALSE, row.names=FALSE)
	
	# ratio persistence time threshold
	stats_df_rpt <- stats_df_unfilt[(stats_df_unfilt$ratio_persist_time >= threshold_ratio_persist_time),]
	stats_df_rpt <- stats_df_rpt[order(-stats_df_rpt$ratio_persist_time, -stats_df_rpt$ratio_sizes, -stats_df_rpt$logistic_growth, stats_df_rpt$logistic_growth_p),]
	write.csv(stats_df_rpt, file=glue('{output_dir}/stats_ratio_persist_time_threshold.csv'), quote=FALSE, row.names=FALSE)
	
	# logistic growth p-filtered
	stats_df_lg_p <- stats_df_unfilt[(stats_df_unfilt$logistic_growth_p > -1) & (stats_df_unfilt$logistic_growth_p <= 0.05),]
	stats_df_lg_p <- stats_df_lg_p[order(-stats_df_lg_p$logistic_growth, stats_df_lg_p$logistic_growth_p, -stats_df_lg_p$ratio_persist_time, -stats_df_lg_p$ratio_sizes),]
	write.csv(stats_df_lg_p, file=glue('{output_dir}/stats_logit_growth_p_threshold.csv'), quote=FALSE, row.names=FALSE)
	
	# Get intersection between applied filters
	stats_df_intersect <- merge(stats_df_rs, stats_df_rpt) # by="node" 
	stats_df_intersect <- merge(stats_df_intersect, stats_df_lg_p)
	write.csv(stats_df_intersect, file=glue('{output_dir}/stats_intersection.csv'), quote=FALSE, row.names=FALSE)
	
	# Get union between applied filters to plot trees starting on node
	stats_df_union <- rbind(stats_df_rs,stats_df_rpt,stats_df_lg_p) 
	stats_df_union <- unique(stats_df_union); rownames(stats_df_union) <- NULL
	stats_df_union <- stats_df_union[order(-stats_df_union$ratio_sizes, -stats_df_union$ratio_persist_time, -stats_df_union$logistic_growth, stats_df_union$logistic_growth_p),]
	write.csv(stats_df_union, file=glue('{output_dir}/stats_union.csv'), quote=FALSE, row.names=FALSE)
	
	message(glue("Nodes matching threshold >= {threshold_ratio_sizes} of ratio sizes: {nrow(stats_df_rs)}"))
	message(glue("Nodes matching threshold >= {threshold_ratio_persist_time} of ratio persistence time: {nrow(stats_df_rpt)}"))
	message(glue("Nodes matching threshold <= 0.05 of logistic regression p-value: {nrow(stats_df_lg_p)}"))
	message(glue("Nodes matching all thresholds: {nrow(stats_df_intersect)}"))
	message(glue("Nodes matching at least one threshold: {nrow(stats_df_union)}"))
	
	if(length(tre$tip.label) < 5000){
		if (!dir.exists(glue("{output_dir}/trees/")))
			suppressWarnings( dir.create(glue("{output_dir}/trees/")) )
		message(glue("Plotting trees zoomed on node matching at least one of the thresholds..."))
		p <- ggtree(tre, mrsd=max_date, as.Date=TRUE) + scale_x_date(date_labels="%b\n%Y", date_breaks="4 months") + theme_tree2(axis.text.x=element_text(size=6))
		for(i in 1:nrow(stats_df_union)) { #length(tgt_nodes)
			p1 <- zoomClade(p, node=stats_df_union[i,1])
			suppressMessages( ggsave(file=glue('{output_dir}/trees/{stats_df_union[i,1]}.pdf'), plot=p1, dpi=600) )
		}
	}
	
	message(paste0("Number of descendants vector length: ", length(ndesc)))

	#return(list(sisters, ndesc, tmrca))
	return(list(stats_df, sisters, ndesc, tmrca))
}

# test_ebola <- cladeScore(ebov_tre, ebov_md, min_descendants=5, max_descendants=75, min_cluster_age_yrs=0.2/12, min_date=as.Date("2014-07-01"),
# 																									max_date=as.Date("2015-10-24"),output_dir="test_ebola", threshold_ratio_sizes=2, threshold_ratio_persist_time=1.5, gen_time=16.6,
# 																									root_on_tip="LIBR10245_2014-07-01", root_on_tip_sample_time=2014.496)

# Problems in branch lengths here (times are not reliable since root estimate is fluctuating from 2012 to 2019 depending on the inclusion criteria)
# 2019-12-30 to 2020-06-30 (first lineages)
start <- Sys.time()
cladeScore_1st <- cladeScore(sc2_tre, sc2_md, min_descendants=10, max_descendants=1e3, min_cluster_age_yrs=0.5/12, min_date=as.Date("2019-12-30"),
																			max_date=as.Date("2020-06-30"), output_dir="sc2_root_to_jun2020", threshold_ratio_sizes=1.5, threshold_ratio_persist_time=1.25, gen_time=7/365,
																			root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995)
end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed: ",total_time,"mins"))

# 2020-07-01 to 2020-12-31 (B.1.177 + start Alpha)
start <- Sys.time()
cladeScore_2nd <- cladeScore(sc2_tre, sc2_md, min_descendants=25, max_descendants=5e3, min_cluster_age_yrs=1/12, min_date=as.Date("2020-07-01"),
																													max_date=as.Date("2020-12-31"), output_dir="sc2_jul2020_to_dec2020", threshold_ratio_sizes=2, threshold_ratio_persist_time=1.5, gen_time=7/365,
																													root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995)
end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed: ",total_time,"mins"))

# 2021-01-01 to 2021-05-31 (Alpha + Delta rapidly replacing)
start <- Sys.time()
cladeScore_3rd <- cladeScore(sc2_tre, sc2_md, min_descendants=50, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=as.Date("2021-01-01"),
																													max_date=as.Date("2021-05-31"), output_dir="sc2_jan2021_to_may2021", threshold_ratio_sizes=2, threshold_ratio_persist_time=1.5, gen_time=7/365,
																													root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995)
end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed: ",total_time,"mins"))

# 2021-06-01 to 2021-12-31 (Delta + Omicron BA.1 rapidly replacing)
start <- Sys.time()
cladeScore_4th <- cladeScore(sc2_tre, sc2_md, min_descendants=50, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=as.Date("2021-06-01"),
																													max_date=as.Date("2021-12-31"), output_dir="sc2_jun2021_to_dec2021", threshold_ratio_sizes=2, threshold_ratio_persist_time=1.5, gen_time=7/365,
																													root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995)
end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed: ",total_time,"mins"))

# 2022-01-01 to 2022-04-30 (Omicron BA.1 + BA.2 rapidly replacing)
start <- Sys.time()
cladeScore_5th <- cladeScore(sc2_tre, sc2_md, min_descendants=50, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=as.Date("2022-01-01"),
																													max_date=as.Date("2022-04-30"), output_dir="sc2_jan2022_to_apr2022", threshold_ratio_sizes=2, threshold_ratio_persist_time=1.5, gen_time=7/365,
																													root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995)
end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed: ",total_time,"mins"))

saveRDS(cladeScore_1st, "rds/results/period1.rds")
saveRDS(cladeScore_2nd, "rds/results/period2.rds")
saveRDS(cladeScore_3rd, "rds/results/period3.rds")
saveRDS(cladeScore_4th, "rds/results/period4.rds")
saveRDS(cladeScore_5th, "rds/results/period5.rds")

# Whole tree period (2019-12-30 to 2022-04-30)
start <- Sys.time()
cladeScore_all <- cladeScore(sc2_tre, sc2_md, min_descendants=50, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=as.Date("2019-12-30"),
																													max_date=as.Date("2022-04-30"), output_dir="sc2_whole_period", threshold_ratio_sizes=2, threshold_ratio_persist_time=1.5, gen_time=7/365,
																													root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995) # around 10 minutes to run
end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed: ",total_time,"mins"))

saveRDS(cladeScore_all, "rds/results/whole_period.rds")
