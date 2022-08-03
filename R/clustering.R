library(ape)
library(lubridate)
library(tidyr)
library(dplyr)

tre = readRDS("rds/timetree_ebola.rds")
amd = readRDS("rds/md_ebola.rds")

# TODO: improve memory allocation (avoid matrix?) since it failed for tree composed by 1.2M tips
# load("rds/tree_cleaned.RData")
# tre <- eng_tre_bin
# md_eng_seqnames_match$lineage <- NULL
# amd <- md_eng_seqnames_match
# rm(eng_tre, eng_tre_bin, md_eng_seqnames_match, min_blen)

# TODO: include following parameters in function? max descendants, min_date, max_date, min_blen
# TODO: unique function encapsulating all parameters
min_descentants = 25
min_cluster_age_yrs = 1/12
# ncpu <- 1
# output_dir <- paste0('treeclust-', Sys.Date())

# if (!dir.exists(output_dir))
# 	dir.create(output_dir)

# load tree and md data
amd <- amd[!is.na(amd$sequence_name),]
amd$sample_date <- as.Date(amd$sample_date)
stopifnot(all(tre$tip.label %in% amd$sequence_name)) # enforce that all tips are listed on metadata sequence_name column
if(!('sample_time') %in% colnames(amd)) {
	amd$sample_time <- decimal_date(amd$sample_date)
}

amd$sts <- amd$sample_time
# remove missing dates
amd <- amd[!is.na(amd$sample_time),]

# named vector of sequence_names and sample_times
sts <- setNames(amd$sample_time , amd$sequence_name)
amd <- amd [amd$sequence_name %in% tre$tip.label,] 

stopifnot(is.rooted(tre))

tre <- keep.tip(tre, intersect(tre$tip.label, amd$sequence_name))

# variables and data structures to quickly look up tree data
# based on treestructure/tfpscanner code
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
	t <- poedges[e,2] # postorder traversal of tips column
	ndesc[n] <- ndesc[n] + ndesc[t]
	Ndesc[n] <- Ndesc[n] + Ndesc[t]
	max_desc_time[n] <- max(max_desc_time[n], max_desc_time[t])
	min_desc_time[n] <- min(min_desc_time[n], min_desc_time[t])
}
clade_age <- max_desc_time - min_desc_time # time span of descendant tips

# descendant nodes
descendants <- lapply(1:(ntip+nnode), function(t) integer(Ndesc[t]))  # pre-allocate list to fill with node descendant IDs
for(t in 1:(ntip+nnode)) {
	descendants[[t]][1] <- t
}

Ndesc_idx <- rep(2, ntip+nnode)

for(e in 1:nrow(poedges)) {
	n <- poedges[e,1]
	t <- poedges[e,2]
	i0 <- Ndesc_idx[n]
	i1 <- Ndesc[t] + i0 - 1
	descendants[[n]][i0:i1] <- descendants[[t]]
	Ndesc_idx[n] <- i1 + 1
}

# descendant tips indexes
descendant_tips <- descendants
for(n in (ntip+1):(ntip+nnode)) {
	descendant_tips[[n]] <- descendant_tips[[n]][ descendant_tips[[n]] <= ntip ]
}

# tip labels (20 secs)
descendant_ids <- lapply(1:(ntip+nnode), function(t) {
	na.omit(tre$tip.label[descendant_tips[[t]]])
})

# ancestors (7 secs)
ancestors <- lapply(1:(ntip+nnode), function(t) integer() )
for(e in 1:nrow(preedges)) {
	n <- preedges[e,1]
	t <- preedges[e,2]
	ancestors[[t]] <- c(ancestors[[n]], n)
}

# compute node stats based on conditions
tgt_nodes <- which((ndesc >= min_descentants) & (clade_age >= min_cluster_age_yrs))
tgt_nodes <- tgt_nodes[-1] # exclude root since ndesc = ntip
tgt_node_ancestors <- do.call(c, lapply(tgt_nodes, function(t) ancestors[[t]]))

.compute_ratios <- function(target_nodes, var, var_name) { #threshold
	var_tgt_nodes <- rep(0, length(target_nodes))
	i <- 1
	for(n in target_nodes) {
		var_tgt_nodes[i] <- var[n]
		i <- i+1
	}
	names(var_tgt_nodes) <- target_nodes
	
	# order vector to make sure division >= 1
	var_tgt_nodes <- rev(sort(var_tgt_nodes))
	print(var_tgt_nodes)
	
	# compute ratio of var matrix 
	ratio_var <- matrix(0, nrow = length(var_tgt_nodes), ncol = length(var_tgt_nodes))
	d <- length(var_tgt_nodes)-1
	for(i in 1:(length(var_tgt_nodes)-1)) {
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
	ratio_var_df <- ratio_var_df[ratio_var_df[[var_name]] > 0,] #threshold
	
	ratio_var_df
}

ratio_sizes <- .compute_ratios(tgt_nodes, ndesc, "ratio_sizes") #threshold=2
ratio_time_persist <- .compute_ratios(tgt_nodes, clade_age, "ratio_persistence_time") #threshold=3

# logistic growth of sister clades
# using ratio sizes because it already has all combination of nodes to iterate over 
# and can filter out comparisons of clades with e. g. size ratio > 2 (very unbalanced sizes)
# ratio_size_upper_threshold: values between 1 and this value (including it) will be kept
# TODO adjust by time: get only clades with compatible range of min and max sample times?
.logistic_growth_sister_clades <- function(ratio_sizes, ratio_size_upper_threshold=2, gen_time=6.5) {
	ratio_sizes <- ratio_sizes[ratio_sizes$ratio_sizes <= ratio_size_upper_threshold,]
	#X <- model <- s <- list()
	#rel_growth_gen_time <- p <- c()
	res_list <- list()
	for(i in 1:nrow(ratio_sizes)) {
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
		gen_time <- gen_time
		rel_growth_gen_time <- unname(coef(model)[2] * gen_time)
		p <- s$coefficients[2,4]
		print(paste0("Relative growth per generation = ",rel_growth_gen_time,"; p-value = ",p))
		res_list[[i]] <- cbind(node, node_comp, rel_growth_gen_time, p)
	}
	#list(rel_growth_gen_time, p)
	res_df <- do.call(rbind, res_list)
	res_df
}

logit_growth <- .logistic_growth_sister_clades(ratio_sizes, ratio_size_upper_threshold=2, gen_time=6.5)

# start <- Sys.time()
# end <- Sys.time()
# total_time <- as.numeric (end - start, units = "secs")
# print(paste("Total time elapsed: ",total_time,"secs"))