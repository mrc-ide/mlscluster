library(ape)
library(lubridate)
library(tidyr)
library(dplyr)
library(glue)

# tre = readRDS("rds/timetree_ebola.rds")
# amd = readRDS("rds/md_ebola.rds")

# TODO: improve memory allocation (avoid matrix?) since it failed to create matrix for entire tree composed by 1.2M tips
# OR set maximum number of nodes for comparison
load("rds/tree_cleaned.RData")
tre <- eng_tre_bin
md_eng_seqnames_match$lineage <- NULL
amd <- md_eng_seqnames_match
rm(eng_tre, eng_tre_bin, md_eng_seqnames_match, min_blen)

cladeScore <- function(tre, amd, min_descendants=100, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=NULL, max_date=NULL,
																							output_dir=paste0("clade_score-",Sys.Date()), gen_time=7/365, threshold_ratio_sizes=3, threshold_ratio_persist_time=3,
																							upper_threshold_sizes_logit_regr=5) { #ncpu=1
	
	stopifnot(is.rooted(tre))
	
	if (!dir.exists(output_dir))
		dir.create(output_dir)
	
	max_time <- Inf
	if(!is.null(max_date)) {
		max_time <- decimal_date(max_date)
	}else{
		max_time <- Sys.Date()
	}
	
	min_time <- -Inf
	if(!is.null(min_date)) {
		min_time <- decimal_date(min_date)
	}
	
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
	nde <- node.depth.edgelength(tre)
	ntip <- Ntip(tre)
	nnode <- Nnode(tre)
	
	poedges <- tre$edge[postorder(tre),] #postorder: traverse from the left subtree to the right subtree then to the root
	preedges <- tre$edge[rev(postorder(tre)),] #preorder: traverse from the root to the left subtree then to the right subtree
	
	# number of descendants
	ndesc <- rep(0, ntip+nnode) # tips descending
	ndesc[1:ntip] <- 1 # tips index will have only 1 descendant
	Ndesc <- rep(1, ntip + nnode) # including internal nodes
	
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
	
	# descendant nodes
	descendants <- lapply(1:(ntip+nnode), function(t) integer(Ndesc[tp]))  # pre-allocate list to fill with node descendant IDs
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
	for(n in (ntip+1):(ntip+nnode)) {
		descendant_tips[[n]] <- descendant_tips[[n]][ descendant_tips[[n]] <= ntip ]
	}
	
	# tip labels
	descendant_ids <- lapply(1:(ntip+nnode), function(tp) {
		na.omit(tre$tip.label[descendant_tips[[tp]]])
	})
	
	# ancestors
	ancestors <- lapply(1:(ntip+nnode), function(tp) integer() )
	for(e in 1:nrow(preedges)) {
		n <- preedges[e,1]
		tp <- preedges[e,2]
		ancestors[[tp]] <- c(ancestors[[n]], n)
	}
	
	message(paste0("Number of descendants vector length: ", length(ndesc)))
	message(paste0("Clade ages vector length: ", length(clade_age)))
	
	message(min_descendants)
	message(max_descendants)
	message(min_cluster_age_yrs)
	
	# compute node stats based on conditions
	# message(which((ndesc >= min_descendants)))
	# message(which((ndesc >= max_descendants)))
	# message(which((clade_age >= min_cluster_age_yrs)))
	tgt_nodes <- which((ndesc >= min_descendants) & (ndesc <= max_descendants) & (clade_age >= min_cluster_age_yrs))
	message(paste0("Number of target nodes: ", length(tgt_nodes)))
	#tgt_nodes <- tgt_nodes[-1] # exclude root since ndesc = ntip
	tgt_node_ancestors <- do.call(c, lapply(tgt_nodes, function(tp) ancestors[[tp]]))
	
	# compute ratios general function
	.compute_ratios <- function(target_nodes, var, var_name, threshold=3) { #threshold
		var_tgt_nodes <- rep(0, length(target_nodes))
		i <- 1
		for(n in target_nodes) {
			var_tgt_nodes[i] <- var[n]
			i <- i+1
		}
		names(var_tgt_nodes) <- target_nodes
		
		# order vector to make sure division >= 1
		var_tgt_nodes <- sort(var_tgt_nodes, decreasing=TRUE) #rev(sort(var_tgt_nodes))
		
		# compute ratio of var matrix 
		ratio_var <- matrix(0, nrow = length(var_tgt_nodes), ncol = length(var_tgt_nodes))
		print(object.size(ratio_var),units="Mb")
		d <- length(var_tgt_nodes)-1
		for(i in 1:(length(var_tgt_nodes)-1)) {
			if(i %% 1000 == 0) print(i)
			for(j in i+1:d) {
				#print(paste0("[i=",i,", j=",j,"]"))
				ratio_var[[i,j]] <- var_tgt_nodes[i]/var_tgt_nodes[j]
			}
			d <- d - 1
			#print(paste0("d = ",d))
		}
		
		ratio_var <- ratio_var[,-1]
		rownames(ratio_var) <- names(var_tgt_nodes)
		colnames(ratio_var) <- names(var_tgt_nodes[2:length(var_tgt_nodes)])
		
		ratio_var_df <- as.data.frame(ratio_var)
		ratio_var_df$node_id <- rownames(ratio_var)
		rownames(ratio_var_df) <- NULL
		
		ratio_var_df <- ratio_var_df %>% pivot_longer(!node_id, names_to = "node_comp", values_to = var_name)
		ratio_var_df <- ratio_var_df[(ratio_var_df[[var_name]] > 0) & (ratio_var_df[[var_name]] >= threshold),]
		
		ratio_var_df <- ratio_var_df[order(-ratio_var_df[[var_name]]),]
		
		write.csv(ratio_var_df, file=glue('{output_dir}/{var_name}.csv'), quote=FALSE, row.names=FALSE)
		
		ratio_var_df
	}
	
	message("Computing ratio of sizes...")
	start <- Sys.time()
	ratio_sizes <- .compute_ratios(tgt_nodes, ndesc, "ratio_sizes", threshold=threshold_ratio_sizes)
	end <- Sys.time()
	total_time <- as.numeric (end - start, units = "mins")
	print(paste("Total time elapsed (ratio of sizes): ",total_time,"mins"))
	
	message("Computing ratio of persistence time of sister clades...")
	start <- Sys.time()
	ratio_time_persist <- .compute_ratios(tgt_nodes, clade_age, "ratio_persistence_time", threshold=threshold_ratio_persist_time) #threshold=3
	end <- Sys.time()
	total_time <- as.numeric (end - start, units = "mins")
	print(paste("Total time elapsed (ratio of persistence time): ",total_time,"mins"))
	
	# logistic growth of sister clades
	# using ratio sizes because it already has all combination of nodes to iterate over 
	# and can filter out comparisons of clades with e. g. size ratio > 2 (very unbalanced sizes)
	# ratio_size_upper_threshold: values between 1 and this value (including it) will be kept
	# TODO adjust by time: get only clades with compatible range of min and max sample times?
	.logistic_growth_sister_clades <- function(ratio_sizes, ratio_size_upper_threshold=upper_threshold_sizes_logit_regr, gen_time=6.5) {
		ratio_sizes <- ratio_sizes[(ratio_sizes$ratio_sizes >= 1) & (ratio_sizes$ratio_sizes <= ratio_size_upper_threshold),]
		#X <- model <- s <- list()
		#rel_growth_gen_time <- p <- c()
		res_list <- list()
		for(i in 1:nrow(ratio_sizes)) {
			if(i %% 1000 == 0) print(i)
			node <- unlist(as.integer(ratio_sizes[i,1]))
			node_comp <- unlist(as.integer(ratio_sizes[i,2]))
			lg_tgt_node <- descendant_ids[[node]]
			lg_sister_node  <- descendant_ids[[node_comp]]
			sts_tgt_node <- sts[lg_tgt_node]
			sts_sist_node <- sts[lg_sister_node]
			X <- data.frame(sequence_name=c(names(sts_tgt_node),names(sts_sist_node)),time=c(sts_tgt_node,sts_sist_node), type=c(rep("clade", length(sts_tgt_node)), rep("sister",length(sts_sist_node))))
			X <- na.omit(X)
			model <- glm(type=="clade" ~ time, data=X, family=binomial(link="logit"))
			s <- summary(model)
			rel_growth_gen_time <- unname(coef(model)[2] * gen_time)
			p <- s$coefficients[2,4]
			#print(paste0("Relative growth per generation = ",rel_growth_gen_time,"; p-value = ",p))
			res_list[[i]] <- cbind(node, node_comp, rel_growth_gen_time, p)
		}
		#list(rel_growth_gen_time, p)
		res_df <- do.call(rbind, res_list)
		res_df <- as.data.frame(res_df)
		res_df <- res_df[order(-res_df$rel_growth_gen_time),]
		write.csv(res_df, file=glue('{output_dir}/logistic_regression.csv'), quote=FALSE, row.names=FALSE)
		res_df
	}
	
	message("Computing logistic growth of sister clades...")
	start <- Sys.time()
	logit_growth <- .logistic_growth_sister_clades(ratio_sizes, ratio_size_upper_threshold=upper_threshold_sizes_logit_regr, gen_time=6.5)
	end <- Sys.time()
	total_time <- as.numeric (end - start, units = "mins")
	print(paste("Total time elapsed (logistic growth rate of sister clades): ",total_time,"mins"))
}

cladeScore(tre, amd, min_descendants=25, max_descendants=500, min_cluster_age_yrs=0.2/12, min_date=as.Date("2019-12-01"), 
											max_date=as.Date("2020-06-30"),output_dir="clade_score-end2019-june2020", gen_time=7/365, threshold_ratio_sizes=2, 
											threshold_ratio_persist_time=2,upper_threshold_sizes_logit_regr=3)

# TODO: run for these periods (logistic regression too slow, try to parallelise maybe): 
# DONE: 2019-12-01 to 2020-06-30 (first lineages)
# 2020-07-01 to 2020-12-31 (B.1.177 + start Alpha); 
# 2021-01-01 to 2021-05-31 (Alpha + Delta rapidly replacing);
# 2021-06-01 to 2021-12-31 (Delta + Omicron BA.1 rapidly replacing)
# 2022-01-01 to 2022-04-30 (Omicron BA.1 + BA.2 rapidly replacing)
# Whole tree period (not feasible now due to mem alloc issues on matrix)

# saveRDS(ratio_sizes, "rds/results/ratio_sizes.rds")
# saveRDS(ratio_time_persist, "rds/results/ratio_time_persist.rds")
# saveRDS(logit_growth, "rds/results/logit_growth.rds")
