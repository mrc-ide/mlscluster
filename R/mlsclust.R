#' Compute statistics for sister clades in a time-scaled tree.
#' 
#' @description Calculate three statistics for sister clades (ratio of descendant clade sizes, ratio of persistence times, 
#'     and logistic growth rates) in huge time-scaled phylogenies. Clade-defining homoplasies are matched with these statistics 
#'     to represent a more reliable set of potential mutations that have an advantage for replication but not for transmission.
#'
#' @param tre A time-scaled phylogeny in `ape::phylo` form, which can be estimated using treedater, 
#'     treetime, chronumental, etc. If not rooted, an outgroup contained in the tree and its sample time
#'     must be provided
#' @param amd A data frame containing metadata for each tip of the tree. Mandatory columns are: 
#'     sequence_name, sample_date, lineage, major_lineage, mutations. Mutations must be amino acid replacements separated by the 
#'     pipe "|" character. E.g. "orf1ab:Q676H|synSNP:C3037T|orf1ab:T3750I|synSNP:C13168T|orf1ab:P4715L|synSNP:C18877T|
#'     synSNP:C22444T|S:D614G|S:V1068F|ORF3a:Q57H|synSNP:C26735T|synSNP:C27059T|N:S194L"
#'     An optional metadata column is: region
#' @param min_descendants Clade must have at least this quantity of tips (integer, default: 10)
#' @param max_descendants Clade must have at most this quantity of tips (integer, default: 20e3)
#' @param min_cluster_age_yrs Minimum time span of clade to be included (numeric, default: 1/12)
#' @param min_date Only include samples after (and including) this date (Date object, default: NULL)
#' @param max_date Only include samples before (and including) this date (Date object, default: NULL)
#' @param branch_length_unit How branch lengths were estimated (character) \cr
#'     Choices are 'years' (default) or 'days' \cr
#'     'years' should be used for timetrees estimated using tools such as \emph{treedater} and \emph{treetime} \cr
#'     'days' is appropriate for timetrees estimated using \emph{chronumental}
#' @param rm_seq_artifacts Remove sequencing artifacts based on X (non-syn) and N (SYN) annotations of missing sites from the 
#'     'mutations' metadata column. Default is TRUE and we highly recommend keeping this option. Only change to FALSE if you
#'     have extremely high-quality sequencing (e.g. without ambiguous and missing calls in the alignment used to build the time-
#'     scaled tree) 
#' @param defining_mut_threshold Frequency threshold to for a mutation within a clade to be
#'     considered a defining mutation (numeric between 0 and 1, default: 0.75)
#' @param root_on_tip Will root on this tip if the input tree is still not rooted (character, default: "Wuhan/WH04/2020")
#' @param root_on_tip_sample_time Numeric time that root tip was sampled (numeric, default: 2019.995)
#' @param detailed_output Boolean. If TRUE, generates objects for clade-, lineage-, and region-specific 
#'     tables and figures (default: FALSE). These objects are posteriorly used in the [run_diff_thresholds()] to actually
#'     generate the outputs
#' @param ncpu Number of CPUs for multicore processing (integer, default: 1)
#' 
#' @importFrom ggtree %<+%
#' @importFrom magrittr %>%
#' @importFrom data.table := .N
#'
#' @return Invisibly returns a list with 3 elements: 
#'    * the clustering statistics (data.frame), 
#'    * the target nodes that passed filtering (vector), and
#'    * homoplasy frequency table (data.frame). If `detailed_output=TRUE` returns 3 additional objects useful for the
#'    the [run_diff_thresholds()] function.
#' @export
mlsclust <- function(tre, amd, min_descendants=10, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=NULL, max_date=NULL, branch_length_unit="years", rm_seq_artifacts=TRUE, defining_mut_threshold=0.75, root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995, detailed_output=FALSE, ncpu=1) {
	
	class(tre) <- "phylo"
	
	max_time <- Inf
	if(!is.null(max_date))
		max_time <- lubridate::decimal_date(max_date)
	else
		max_time <- Sys.Date()
	
	min_time <- -Inf
	if(!is.null(min_date))
		min_time <- lubridate::decimal_date(min_date)
	
	#max_time_full <- max(lubridate::decimal_date(as.Date(amd$sample_date)))
	#min_time_full <- min(lubridate::decimal_date(as.Date(amd$sample_date)))
	#diff_stimes <- max_time_full - min_time_full
	
	# load tree and md data
	message("Loading tree and metadata and filtering out...")
	amd <- amd[!is.na(amd$sequence_name),]
	amd$sample_date <- as.Date(amd$sample_date)
	stopifnot(all(tre$tip.label %in% amd$sequence_name)) # enforce that all tips are listed on metadata sequence_name column
	if(!('sample_time') %in% colnames(amd)) {
		amd$sample_time <- lubridate::decimal_date(amd$sample_date)
	}
	if(!('lineage' %in% colnames(amd))) {
		amd$lineage <- 'lineage_not_provided'
	}
	if(!('major_lineage' %in% colnames(amd))) {
		amd$major_lineage <- 'major_lineage_not_provided'
	}
	if(!('region' %in% colnames(amd))) {
		amd$region <- 'region_not_provided'
	}
	
	# remove missing dates and filter by sample times
	amd <- amd[!is.na(amd$sample_time),]
	amd <- amd[(amd$sample_time >= min_time) & (amd$sample_time <= max_time),]
	
	if(!(root_on_tip %in% amd$sequence_name)) {
		stopifnot(root_on_tip %in% tre$tip.label)
		amd[nrow(amd) + 1,] <- data.frame( sequence_name=root_on_tip, sample_date=lubridate::date_decimal(root_on_tip_sample_time), lineage=NA, major_lineage=NA, region=NA, sample_time=root_on_tip_sample_time )	
	}
	message(paste0("Number of sequences included: ", nrow(amd)))
	
	# named vector of sequence_names and sample_times
	sts <- stats::setNames(amd$sample_time , amd$sequence_name)
	amd <- amd[amd$sequence_name %in% tre$tip.label,]
	
	#dropped_tips <- base::setdiff(tre$tip.label, amd$sequence_name)
	#message(paste0("Dropped tips: ", length(dropped_tips)))
	
	#message("Root sequence Wuhan present?")
	#print("Wuhan/WH04/2020" %in% base::intersect(tre$tip.label, amd$sequence_name))
	
	tre <- ape::keep.tip(tre, base::intersect(tre$tip.label, amd$sequence_name))
	tre2 <- tre
	if(!ape::is.rooted(tre)) {
		stopifnot(root_on_tip %in% tre$tip.label)
		if(!(root_on_tip %in% tre$tip.label)) {
			stop("Outgroup sequence missing in the input tree")
		}
		tre2 <- ape::root(tre, outgroup=root_on_tip, resolve.root=TRUE)
		tre <- tre2
	}
	#message(paste0("Tips in the tree: ", length(tre$tip.label)))
	rm(tre2); gc()
	
	# variables and data structures to quickly look up tree data
	# based on treestructure/tfpscanner code
	message("Computing node descendants and times...")
	ntip <- ape::Ntip(tre)
	nnode <- ape::Nnode(tre)
	
	# getting 'coalescent' times
	if(branch_length_unit == "years") {
		nde <- ape::node.depth.edgelength(tre) # depth of a node using branch lengths
	}else if(branch_length_unit == "days") {
		tre$edge.length <- tre$edge.length / 365
		nde <- ape::node.depth.edgelength(tre)
	}else {
		stop("Choices for 'branch_length_unit' are: 'days' and 'years'")
	}
	rh <- max(nde[1:ntip]) # root heights
	stimes <- nde[1:ntip]
	shs <- rh - stimes # time to most recent sample
	nhs = nodeheights <- rh - nde # node heights
	mrst <- max(sts) # most recent sample time
	tmrca <- mrst - nhs # decimal times of coalescent events
	
	poedges <- tre$edge[ape::postorder(tre),] #postorder: traverse from the left subtree to the right subtree then to the root
	preedges <- tre$edge[rev(ape::postorder(tre)),] #preorder: traverse from the root to the left subtree then to the right subtree
	
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
	descendants <- pbmcapply::pbmclapply(1:(ntip+nnode), function(tp) { integer(Ndesc[tp]) }, mc.cores = ncpu)
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
	#daughter_nodes <- descendants
	
	start <- Sys.time() ######
	for(n in (ntip+1):(ntip+nnode)) {
		descendant_tips[[n]] <- descendant_tips[[n]][ descendant_tips[[n]] <= ntip ]
		#daughter_nodes[[n]] <- daughter_nodes[[n]][ daughter_nodes[[n]] > ntip ]
		# sisters[[n]] <- tre$edge[ tre$edge[,1] == n, 2 ] # not scaling well with larger tree
	}
	
	# tip labels
	descendant_ids <- pbmcapply::pbmclapply(1:(ntip+nnode), function(tp) {
		stats::na.omit(tre$tip.label[descendant_tips[[tp]]])
	}, mc.cores = ncpu)
	
	# Get descendant tips for node 
	.get_tips_node <- function(node) {
		tips_node <- descendant_ids[[node]]
		return(tips_node)
	}
	
	.get_comparison_sister_node <- function(node) {
		imed_ancestor <- utils::tail(ancestors[[node]], 1)
		sisters <- tre$edge[ tre$edge[,1] == imed_ancestor, 2 ]
		comparison_node <- base::setdiff(sisters, node)
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
		
		ratio_sizes <- round( length_sisters[1] / length_sisters[2], digits=5)
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
		
		ratio_persist_time <- round( persist_time[1] / persist_time[2], digits=5)
		res_ratio_persist_time <- c(tmrca[node], max_persist_time, ratio_persist_time)
		
		return(res_ratio_persist_time)
	}
	
	# logistic growth of sister clades function
	.logistic_growth_stat <- function(node) {
		lg_res_list <- list()
		comp_res <- .get_comparison_sister_node(node)
		
		if(!is.null(comp_res)) {
			X <- data.frame(sequence_name=c(comp_res[[5]], comp_res[[6]]), time=c(comp_res[[7]], comp_res[[8]]), type=c(rep("node",length(comp_res[[5]])), rep("control",length(comp_res[[6]]))))
			X <- stats::na.omit(X)
			model <- suppressWarnings( stats::glm(type=="node" ~ time, data=X, family=stats::binomial(link="logit")) )
			s <- summary(model)
			rel_growth <- unname(stats::coef(model)[2]) # * gen_time
			p <- s$coefficients[2,4]
			lg_res_list <- cbind(rel_growth, p, comp_res[[2]])
		}else {
			lg_res_list <- cbind(-1, -1, -1)
		}
		return(lg_res_list)
	}
	
	# Extract node defining muts
	.node_muts <- function(node, mut_var="mutations") {
		comp_res <- .get_comparison_sister_node(node)
		
		# metadata from node/clade sequences
		md_itn <- amd[ amd$sequence_name %in% comp_res[[5]], ]
		# metadata from sister sequences
		md_its <- amd[ amd$sequence_name %in% comp_res[[6]], ]
		
		if( nrow(md_itn) == 0 | nrow(md_its) == 0 )
			return(NA)
		
		if(rm_seq_artifacts) {
			# Table and df of clade mutations and their freqs
			vtab_node_freqs <- sort( table( do.call( c, strsplit( md_itn[[mut_var]], split='\\|' )  ) ) / nrow( md_itn ) )
			vtab_node_freqs <- as.data.frame(vtab_node_freqs)
			colnames(vtab_node_freqs) <- c("defining_mut","Freq")
			vtab_node_freqs$prot <- sub("\\:.*", "", vtab_node_freqs$defining_mut)
			vtab_node_freqs$site <- sub('.*:', "", vtab_node_freqs$defining_mut)
			vtab_node_freqs$site <- parse_number(vtab_node_freqs$site)
			vtab_node_freqs$prot_site <- paste0(vtab_node_freqs$prot,":",vtab_node_freqs$site)
			vtab_node_freqs$prot_anc_site <- str_sub(vtab_node_freqs$defining_mut, start=1, end=-2)
			vtab_node_freqs$mutsite <- str_sub(vtab_node_freqs$defining_mut, -1)
			
			# Table and df of sister mutations and their freqs
			vtab_sister_freqs <- sort( table( do.call( c, strsplit( md_its[[mut_var]], split='\\|' )  ) ) / nrow( md_its ) )
			vtab_sister_freqs <- as.data.frame(vtab_sister_freqs)
			colnames(vtab_sister_freqs) <- c("defining_mut","Freq")
			vtab_sister_freqs$prot <- sub("\\:.*", "", vtab_sister_freqs$defining_mut)
			vtab_sister_freqs$site <- sub('.*:', "", vtab_sister_freqs$defining_mut)
			vtab_sister_freqs$site <- parse_number(vtab_sister_freqs$site)
			vtab_sister_freqs$prot_site <- paste0(vtab_sister_freqs$prot,":",vtab_sister_freqs$site)
			vtab_sister_freqs$prot_anc_site <- str_sub(vtab_sister_freqs$defining_mut, start=1, end=-2)
			vtab_sister_freqs$mutsite <- str_sub(vtab_sister_freqs$defining_mut, -1)
			
			# Get common sites between node and sister that can be changed if X/N missing sites
			df_common_sites <- vtab_node_freqs %>% dplyr::left_join(vtab_sister_freqs, by="prot_site", multiple="all") #left_join
			df_common_sites <- df_common_sites[!is.na(df_common_sites$mutsite.y),]
			#print("df_common_sites")
			#print(df_common_sites) # IMPORTANT: this only get 4 first possibilities
			
			# Get SYN sites with N in node and other char in sister
			df_common_sites_syn <- df_common_sites[(df_common_sites$prot.x == "SYNSNP" & df_common_sites$mutsite.x == "N")|(df_common_sites$prot.x == "SYNSNP" & df_common_sites$mutsite.y == "N"),]
			df_common_sites_syn <- df_common_sites_syn[(df_common_sites_syn$mutsite.x != df_common_sites_syn$mutsite.y)|(df_common_sites_syn$mutsite.y != df_common_sites_syn$mutsite.x),]
			# Get NON-SYN sites with X in node and other char in sister
			df_common_sites_nonsyn <- df_common_sites[(df_common_sites$prot.x != "SYNSNP" & df_common_sites$mutsite.x == "X")|(df_common_sites$prot.x != "SYNSNP" & df_common_sites$mutsite.y == "X"),]
			df_common_sites_nonsyn <- df_common_sites_nonsyn[(df_common_sites_nonsyn$mutsite.x != df_common_sites_nonsyn$mutsite.y)|(df_common_sites_nonsyn$mutsite.y != df_common_sites_nonsyn$mutsite.x),]
			#print(df_common_sites_nonsyn) # IMPORTANT: this only gets 1 possibility (no 1, X in node and N in sister)
			# Non-N / non-X adjusted based on node/sister without missing mutation
			df_common_sites_syn$final_mut <- ifelse(df_common_sites_syn$mutsite.x == "N", yes=as.character(df_common_sites_syn$mutsite.y), no=as.character(df_common_sites_syn$mutsite.x))
			df_common_sites_syn$defining_mut_final <- paste0(df_common_sites_syn$prot_anc_site.x, df_common_sites_syn$final_mut)
			
			df_common_sites_nonsyn$final_mut <- ifelse(df_common_sites_nonsyn$mutsite.x == "X", yes=as.character(df_common_sites_nonsyn$mutsite.y), no=as.character(df_common_sites_nonsyn$mutsite.x))
			df_common_sites_nonsyn$defining_mut_final <- paste0(df_common_sites_nonsyn$prot_anc_site.x, df_common_sites_nonsyn$final_mut)
			
			# Fixed SYN and NON-SYN sites
			fixed_sites_node <- rbind(df_common_sites_syn, df_common_sites_nonsyn) #bind_rows
			# Fixed sites matching with node mutations
			match_fixed_node <- fixed_sites_node %>% left_join(vtab_node_freqs, by=c("defining_mut_final"="defining_mut"))
			match_fixed_node <- match_fixed_node[!is.na(match_fixed_node$mutsite),]
			match_fixed_node$prot_anc_site <- str_sub(match_fixed_node$defining_mut_final, start=1, end=-2)
			match_fixed_node$mutsite <- str_sub(match_fixed_node$defining_mut_final, -1)
			# Fix proportion of adjusted muts on node
			match_fixed_node_prop <- match_fixed_node %>% group_by(prot_anc_site) %>% mutate(Freq_adj=sum(Freq.x, Freq.y)) %>% ungroup()
			# If >100%, maximum is 1
			match_fixed_node_prop$Freq_adj <- ifelse(match_fixed_node_prop$Freq_adj>1, 1, as.numeric(match_fixed_node_prop$Freq_adj))
			match_fixed_node_prop <- match_fixed_node_prop %>% select(defining_mut_final, Freq_adj)
			colnames(match_fixed_node_prop) <- c("defining_mut","Freq")
			vtab_node_freqs <- vtab_node_freqs %>% select(defining_mut, Freq)

			all_node_muts <- dplyr::setdiff(vtab_node_freqs, match_fixed_node_prop)
			all_node_muts <- all_node_muts[!duplicated(all_node_muts$defining_mut), ]
			#print("all_node_muts AFTER removing dups") # IMPORTANT: this still includes X (probably because not considering no 3 from notebook)
			#print(all_node_muts)
			all_node_muts <- all_node_muts[all_node_muts$Freq > defining_mut_threshold,]
			
			defining_mut_df <- setDT(all_node_muts)[!vtab_sister_freqs, on = "defining_mut"] #all_sister_muts
			defining_mut_df$mutsite <- str_sub(defining_mut_df$defining_mut, -1)
			defining_mut_df$prot <- sub("\\:.*", "", defining_mut_df$defining_mut)
			defining_mut_df <- as.data.frame(defining_mut_df)
			defining_mut_df$defining_mut <- as.character(defining_mut_df$defining_mut)
			
			# Removes X/N muts not found in both (opts 6 & 8)
			defining_mut_df1 <- defining_mut_df[(defining_mut_df$prot != "SYNSNP" & defining_mut_df$mutsite != "X"),]
			defining_mut_df2 <- defining_mut_df[(defining_mut_df$prot == "SYNSNP" & defining_mut_df$mutsite != "N"),]
			defining_mut_df <- rbind(defining_mut_df1, defining_mut_df2)
			defining_muts <- unname(defining_mut_df$defining_mut)
			
			return(defining_muts)
		} else {
			vtab_node = sort( table( do.call( c, strsplit( md_itn[[mut_var]], split='\\|' )  ) ) / nrow( md_itn ) )
			vtab_sister = sort( table( do.call( c, strsplit( md_its[[mut_var]], split='\\|' )  ) ) / nrow( md_its ) )
			defining_muts <- base::setdiff( names(vtab_node[vtab_node > defining_mut_threshold]), names(vtab_sister[vtab_sister > defining_mut_threshold]) )
			return(defining_muts)
		}
		
	}
	
	# Detect independent occurences of defining mutations (homoplasies) across different nodes
	.find_homoplasies <- function(node_list, major_lineage_nodes) {
		def_muts_nodes <- pbmcapply::pbmclapply(1:length(node_list), function(tp) { .node_muts(node_list[[tp]]) }, mc.cores = ncpu)
		
		names(def_muts_nodes) <- pbmcapply::pbmclapply( 1:length(node_list), function(tp) { node_list[[tp]] }, mc.cores=ncpu )
		
		def_muts_nodes_df <- data.table::rbindlist( lapply(def_muts_nodes, function(x) { data.table::data.table(x) }), idcol="node")
		
		def_muts_nodes_df <- stats::na.omit(def_muts_nodes_df)
		names(def_muts_nodes_df) <- c("node","defining_mut")
		def_muts_nodes_df$defining_mut <- toupper( def_muts_nodes_df$defining_mut )
		
		def_muts_nodes_df <- base::merge(def_muts_nodes_df, major_lineage_nodes, by="node")
		
		tab_def_muts_df <- data.table::setDT(def_muts_nodes_df)[, .(Freq_homopl = .N), by = .(defining_mut)]
		
		homoplasy_count_df <- base::merge(def_muts_nodes_df, tab_def_muts_df, by="defining_mut")
		homoplasy_count_df <- homoplasy_count_df[as.numeric(homoplasy_count_df$Freq_homopl) > 1,]
		homoplasy_count_df <- homoplasy_count_df %>% dplyr::select(defining_mut, node, major_lineage, Freq_homopl)
		homoplasy_count_df$node <- as.integer(homoplasy_count_df$node)
		
		homoplasy_count_df <- as.data.frame(homoplasy_count_df)
		
		homoplasy_count_df
	}
	
	# Function to extract major_lineage for target nodes
	.major_lineage_all_nodes <- function(tips_nodes_list, max_rows=1) {
		
		lineages_node <- pbmcapply::pbmclapply(1:length(tips_nodes_list), function(tp) {
			amd$major_lineage[ match(tips_nodes_list[[tp]], amd$sequence_name) ]
		}, mc.cores = ncpu)
		
		lineage_node_table <- pbmcapply::pbmclapply(1:length(tips_nodes_list), function(tp) {
			sort(table(lineages_node[[tp]]), decreasing=TRUE) / length(lineages_node[[tp]])
		}, mc.cores = ncpu)
		
		perc_df <- list()
		for(i in 1:length(tips_nodes_list)) {
			if(length(lineage_node_table[[i]]) > 1) {
				perc_df[[i]] <- as.data.frame(lineage_node_table[[i]])
				colnames(perc_df[[i]]) <- c("major_lineage","Freq")
			} else {
				perc_df[[i]] <- data.frame(major_lineage=lineages_node[[i]][1], Freq=1)
			}
			perc_df[[i]] <- perc_df[[i]][ 1:min(nrow(perc_df[[i]]),max_rows), ] # get only most prevalent
			perc_df[[i]]$Freq <- paste0(round(perc_df[[i]]$Freq*100), "%")
			perc_df[[i]]$node <- names(tips_nodes_list)[[i]] # node column receive node id
		}
		perc_df_res <- data.table::rbindlist(perc_df)
		return(perc_df_res)
	}
	
	# compute node stats based on conditions
	tgt_nodes <- which((ndesc >= min_descendants) & (ndesc <= max_descendants) & (clade_age >= min_cluster_age_yrs))
	message(glue::glue("Number of target nodes: {length(tgt_nodes)}"))
	#message("Target nodes are:")
	#print(tgt_nodes)
	
	# Iterate over tgt_nodes (preserving node ids) to get respective descendant tips
	tips_tgt_nodes <- list()
	idx <- 1
	for(x in tgt_nodes) {
		tips_tgt_nodes[[idx]] <- .get_tips_node(x)
		idx <- idx + 1
	}
	names(tips_tgt_nodes) <- tgt_nodes
	
	message("Getting major lineages associated with nodes: ")
	maj_lin_nodes <- .major_lineage_all_nodes(tips_tgt_nodes)
	rm(tips_tgt_nodes); gc()
	
	message("Getting homoplasies: ")
	homoplasies <- .find_homoplasies(tgt_nodes, maj_lin_nodes) # will consider homoplasies only for nodes passing number of descendants and clade age filters
	rm(maj_lin_nodes); gc()
	
	message("Time to populate stats:") ######
	start <- Sys.time() ######
	message("Calculating statistics: ")
	
	report_nodes <- tgt_nodes[seq(1, length(tgt_nodes), by = 5000)]
	
	st0 <- st1 <- st2 <- st3 <- st4 <- comb_stats <- list()
	for(n in seq_along(tgt_nodes)) { #(ntip+2):(ntip+nnode)
		if(tgt_nodes[n] %in% report_nodes) {
			message(glue::glue("Progress: {round(100 * n / length(tgt_nodes))}%"))
		}
		sisters <- .get_comparison_sister_node(tgt_nodes[n])[[1]]
		if(!is.null(sisters)) {
			st0 <- .node_muts(tgt_nodes[n])
			st1 <- .ratio_sizes_stat(tgt_nodes[n])
			st2 <- .ratio_persist_time_stat(tgt_nodes[n])
			st3 <- .logistic_growth_stat(tgt_nodes[n])
			st4 <- ""
			if(tgt_nodes[n] %in% homoplasies$node) {
				st4 <- "yes"
			}else{
				st4 <- "no"
			}
			comb_stats[[n]] <- data.frame(tgt_nodes[n], paste(st0, collapse="|"), sisters[1], sisters[2], st1[1], st1[2], st1[3], st2[1], st2[2], st2[3], st2[4], st3[,1], st3[,2], st4)
			colnames(comb_stats[[n]]) <- NULL
			rownames(comb_stats[[n]]) <- NULL
		}
	}
	
	end <- Sys.time(); total_time <- as.numeric (end - start, units = "mins"); message(paste("Total time elapsed:",total_time,"mins")) ######
	
	stats_df <- data.table::rbindlist(comb_stats, use.names = FALSE)
	colnames(stats_df) <- c("node","defining_muts","comp_sister1","comp_sister2","size1","size2","ratio_sizes","min_time_node","max_time1","max_time2","ratio_persist_time","logistic_growth", "logistic_growth_p", "homoplasies") #"diff_time_sister1","diff_time_sister2", "control_clade(s)" #"major_lineage"
	stats_df <- as.data.frame(stats_df)
	
	rm(comb_stats); gc()
	
	# Filter and write to CSV files
	# unfiltered output (just removing filtered based on overall min_descendants => -1 on all values)
	stats_df_unfilt <- stats_df[(stats_df$ratio_sizes > 0) & (stats_df$ratio_persist_time > 0),]
	rm(stats_df); gc()
	
	rm(rh, stimes, shs, nhs, mrst, tmrca, poedges, preedges, ndesc, Ndesc, tlsts, max_desc_time, min_desc_time, persistence_time, descendants, descendant_tips)
	gc()
	
	if(detailed_output) {
		tips_sisters1 <- pbmcapply::pbmclapply(seq_along(tgt_nodes), function(tp) { .get_comparison_sister_node(tgt_nodes[[tp]])[[5]] }, mc.cores = ncpu)
		tips_sisters2 <- pbmcapply::pbmclapply(seq_along(tgt_nodes), function(tp) { .get_comparison_sister_node(tgt_nodes[[tp]])[[6]] }, mc.cores = ncpu) 
		
		return(invisible(list(stats_df_unfilt, tgt_nodes, homoplasies, tips_sisters1, tips_sisters2, amd)))
	} else {
		return(invisible(list(stats_df_unfilt, tgt_nodes, homoplasies)))
	}
}


#' Extract mutations under multilevel selection for different cluster/quantile thresholds (0.25 to 25%) of statistics
#'
#' @description Detect nodes and their defining-homoplasies with the three statistics for sister clades within a small quantile
#'    threshold (0.25 to 25%). This means that number of descendants, persistence time or logistic growth rate of the numerator or target 
#'    is smaller than denominator or comparator, leading to the retrieval of a small (<1) value. This function also performs annotations by
#'    genomic region, different homoplasies at the same site, and within a 3-amino acid sliding window. It also provides detailed outputs if needed. 
#' 
#' @param stats_df_unfilt First list element returned by [mlsclust()]. A data.frame containing node/clade, defining_muts, 
#'    sister clades, ratio of sizes, ratio of persistence time, logistic growth rate, and if homoplasy.
#' @param tgt_nodes Second list element returned by [mlsclust()]. A vector with the node identifiers passing filtering criteria given by
#'    `min_descendants`, `max_descendants`, `min_cluster_age_yrs`, `max_date`, and `min_date` parameters.
#' @param homoplasies Third list element returned by [mlsclust()]. A data.frame (before matching with the quantile thresholds 
#'    of the statistics) containing defining_mut, node, major_lineage, Freq_homopl
#' @param output_dir Directory path where results will be saved
#' @param quantile_choice Number of subsets of equal size to partition statistics. Numeric value, default is 1/100, meaning 100 splits
#'    (one for each 1% of the statistics). E.g., 1, 2, 3, 4, 5, ..., 95, 96, 97, 98, 99, 100%
#' @param quantile_threshold_ratio_sizes Threshold of ratio sizes to be met for clade and respective defining-homoplasies to be 
#'    considered under multilevel selection
#' @param quantile_threshold_ratio_persist_time As above, but for the ratio of persistence time
#' @param quantile_threshold_logit_growth As above, but for the logistic growth rate
#' @param threshold_keep_lower Whether to keep values that are lower (default: TRUE) than the selected quantile_threshold or higher (FALSE)
#' @param compute_tree_annots Whether to plot segregating sites annotated in the time-scaled tree using \emph{ggtree} (default: FALSE)
#' @param plot_global_mut_freq Whether to use \emph{outbreakinfo} package to plot mutation prevalence worldwide. Homoplasies for all target nodes 
#'    passing filtering are plotted, not only the detected by clustering statistics (default: FALSE)
#' @param detailed_output Plot node-specific tables summarising (major) lineages, regions and sequences. Should only be set to TRUE if
#'    previous run of the [mlsclust()] also had this parameter set to TRUE
#' @param desc_sisters Tips descending for each sister of the clade (default: NULL). Should only be provided if `detailed_output==TRUE`. Needs to be 
#'    passed as a list of the returned elements 4 (tips_sisters1) and 5 (tips_sisters2) of [mlsclust()]. E.g. assuming `res` is the returned object of 
#'    [mlsclust()] when `detailed_output==TRUE`, then this parameter should be supplied as `desc_sisters=list(res[[4]], res[[5]])`
#' @param amd Filtered metadata matching tips in the tree (default: NULL). Should only be provided if `detailed_output==TRUE`. Passed as the 
#'    6th element returned by [mlsclust()], e.g. `amd=res[[6]]`
#'
#' @return Invisibly returns a contingency table (data.frame) containing the defining_mut, whether it is clustered (TFP) or not, and Freq of homoplasy.
#'    Also writes several tables to the supplied `output_dir` with: 
#'    * raw statistic outputs,
#'    * union and intersection of detected sites across stats, 
#'    * homoplasies across all target nodes, only detected, and not detected, 
#'    * homoplasies overlapping with dN/dS positively 
#'    selected sites,
#'    * homoplasies within a 3 amino acid sliding window, and
#'    * different homoplasies at the exact same sites
#'    
#'    The main returned file used for subsequent analyses is `clustered_all_df.csv`, which contains both TFP and non-TFP homoplasies and
#'    additional relevant columns for statistical analyses. 
#' @export
run_diff_thresholds <- function(stats_df_unfilt, tgt_nodes, homoplasies, output_dir=paste0("mlscluster-results-",Sys.Date()), quantile_choice=1/100, quantile_threshold_ratio_sizes=1, quantile_threshold_ratio_persist_time=1, quantile_threshold_logit_growth=1, threshold_keep_lower=TRUE, compute_tree_annots=FALSE, plot_global_mut_freq=FALSE, detailed_output=FALSE, desc_sisters, amd) {
	
	system(glue::glue("mkdir -p {output_dir}"))
	
	# Authenticate user in GISAID if wants global prevalence of identified mutations
	if(plot_global_mut_freq) {
		#if getting errors try 'sudo apt install libudunits2-dev libgdal-dev'
		outbreakinfo::authenticateUser()
	}
	
	stats_df_unfilt <- stats_df_unfilt[base::order(as.numeric(stats_df_unfilt$ratio_sizes), as.numeric(stats_df_unfilt$ratio_persist_time), as.numeric(stats_df_unfilt$logistic_growth), decreasing=TRUE),]
	utils::write.csv(stats_df_unfilt, file=glue::glue('{output_dir}/stats_unfiltered.csv'), quote=FALSE, row.names=FALSE)
	
	#quantile_options <- c(1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10)
	.extract_value_below_quantile_threshold <- function(df, var, quantile_choice_user, quantile_threshold) {
		if(quantile_choice_user >= 1/400) {
			quant <- stats::quantile(as.numeric(var), probs=seq(0,1,quantile_choice_user)) #na.rm=TRUE
			message("Quantiles")
			print(quant)
			value_chosen_quant <- unname(quant[names(quant) == quantile_threshold])
			value_chosen_quant <- round(value_chosen_quant, digits=5)
			message(glue::glue("Threshold value at quantile {quantile_threshold}"))
			message(value_chosen_quant)
			if(threshold_keep_lower)
				res_df <- df[(as.numeric(var) < value_chosen_quant),]
			else {
				res_df <- df[(as.numeric(var) > value_chosen_quant),]
			}
		}else {
			stop("Choices for quantile are: {1/2, ... , 1/400}")
		}
		return(res_df)
	}
	
	# ratio sizes threshold
	quantile_threshold_ratio_sizes <- paste0(quantile_threshold_ratio_sizes,"%")
	stats_df_rs_all <- .extract_value_below_quantile_threshold(stats_df_unfilt, stats_df_unfilt$ratio_sizes, quantile_choice, quantile_threshold_ratio_sizes)
	stats_df_rs <- stats_df_rs_all[, c(1:7,14)]
	utils::write.csv(stats_df_rs, file=glue::glue('{output_dir}/stats_ratio_size_threshold.csv'), quote=FALSE, row.names=FALSE)
	
	# ratio persistence time threshold
	quantile_threshold_ratio_persist_time <- paste0(quantile_threshold_ratio_persist_time,"%")
	stats_df_rpt_all <- .extract_value_below_quantile_threshold(stats_df_unfilt, stats_df_unfilt$ratio_persist_time, quantile_choice, quantile_threshold_ratio_persist_time)
	stats_df_rpt_all <- stats_df_rpt_all[base::order(as.numeric(stats_df_rpt_all$ratio_persist_time), decreasing=TRUE),]
	stats_df_rpt <- stats_df_rpt_all[, c(1:4,8:11,14)]
	utils::write.csv(stats_df_rpt, file=glue::glue('{output_dir}/stats_ratio_persist_time_threshold.csv'), quote=FALSE, row.names=FALSE)
	
	# logistic growth p-filtered
	quantile_threshold_logit_growth <- paste0(quantile_threshold_logit_growth,"%") #added
	stats_df_lg_p_all <- stats_df_unfilt[(stats_df_unfilt$logistic_growth_p > -1) & (stats_df_unfilt$logistic_growth_p <= 0.05),] #stats_df_lg_p <-
	stats_df_lg_p_all <- .extract_value_below_quantile_threshold(stats_df_lg_p_all, stats_df_lg_p_all$logistic_growth, quantile_choice, quantile_threshold_logit_growth) #added
	stats_df_lg_p_all <- stats_df_lg_p_all[base::order(as.numeric(stats_df_lg_p_all$logistic_growth), decreasing=TRUE),]
	stats_df_lg_p_all <- stats::na.omit(stats_df_lg_p_all)
	stats_df_lg_p <- stats_df_lg_p_all[, c(1,12:13,14)]
	utils::write.csv(stats_df_lg_p, file=glue::glue('{output_dir}/stats_logit_growth_p_threshold.csv'), quote=FALSE, row.names=FALSE)
	rm(stats_df_unfilt); gc()
	
	# Get intersection between applied filters
	stats_df_intersect <- base::merge(stats_df_rs, stats_df_rpt, by="node")
	stats_df_intersect <- base::merge(stats_df_intersect, stats_df_lg_p, by="node")
	stats_df_intersect <- stats_df_intersect[base::order(as.numeric(stats_df_intersect$ratio_sizes), as.numeric(stats_df_intersect$ratio_persist_time), as.numeric(stats_df_intersect$logistic_growth), decreasing=TRUE),]
	utils::write.csv(stats_df_intersect, file=glue::glue('{output_dir}/stats_intersection.csv'), quote=FALSE, row.names=FALSE)
	
	# Get union between applied filters to plot trees starting on node
	stats_df_rs_all <- as.data.frame(stats_df_rs_all); stats_df_rpt_all <- as.data.frame(stats_df_rpt_all); stats_df_lg_p_all <- as.data.frame(stats_df_lg_p_all)
	if(nrow(stats_df_rs_all) != 0)
		colnames(stats_df_rs_all) <- c("node","defining_muts","comp_sister1","comp_sister2","size1","size2","ratio_sizes","min_time_node","max_time1","max_time2","ratio_persist_time","logistic_growth","logistic_growth_p","homoplasies") #"major_lineage"
	if(nrow(stats_df_rpt_all) != 0)
		colnames(stats_df_rpt_all) <- c("node","defining_muts","comp_sister1","comp_sister2","size1","size2","ratio_sizes","min_time_node","max_time1","max_time2","ratio_persist_time","logistic_growth","logistic_growth_p","homoplasies") #"major_lineage"
	if(nrow(stats_df_lg_p_all) != 0)
		colnames(stats_df_lg_p_all) <- c("node","defining_muts","comp_sister1","comp_sister2","size1","size2","ratio_sizes","min_time_node","max_time1","max_time2","ratio_persist_time","logistic_growth","logistic_growth_p","homoplasies") #"major_lineage"
	stats_df_union <- dplyr::bind_rows(stats_df_rs_all, stats_df_rpt_all, stats_df_lg_p_all)
	stats_df_union <- stats_df_union[base::order(as.numeric(stats_df_union$ratio_sizes), as.numeric(stats_df_union$ratio_persist_time), as.numeric(stats_df_union$logistic_growth), decreasing=TRUE),]
	utils::write.csv(stats_df_union, file=glue::glue('{output_dir}/stats_union.csv'), quote=FALSE, row.names=FALSE)
	rm(stats_df_rs_all, stats_df_rpt_all, stats_df_lg_p_all); gc()
	
	# Annotate homoplasies in: (i) S regions of potential significance
	.annotate_s_homoplasies <- function(homopl_df) {
		regions_s <- utils::read.csv(system.file("extdata", "spike_regions.tsv", package="mlscluster"),sep="\t", header=T)
		# If has homoplasy in S
		homopl_df$s_mut <- ifelse( grepl(homopl_df$defining_mut, pattern ='^[S]:') , "yes", "no")
		# Get S mut coordinate only
		rgx_s <- regexpr('^[S]:[A-Z]{1}[0-9]{1,5}[A-Z]{1}', homopl_df$defining_mut)
		comm_coord <- rep(NA,length(homopl_df$defining_mut))
		comm_coord[rgx_s != -1] <- stringr::str_sub( regmatches(homopl_df$defining_mut, rgx_s) , 4, -2)
		comm_coord <- as.integer(comm_coord)
		comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("s_mut_coord")
		homopl_df <- cbind(homopl_df, comm_coord_df)
		data.table::setDT(homopl_df); data.table::setDT(regions_s)
		# Check if coordinate within boundaries of regions of interest (NTD, RBD, FCS)
		homopl_df[regions_s, on=c("s_mut_coord>=start", "s_mut_coord<=end"), s_mut_region_interest := region]
		
		homopl_df
	}
	
	# Annotate homoplasies in N(ucleocapsid) regions of potential significance
	.annotate_n_homoplasies <- function(homopl_df) {
		regions_n <- utils::read.csv(system.file("extdata", "n_regions.tsv", package="mlscluster"),sep="\t", header=T)
		# If has homoplasy in N
		homopl_df$n_mut <- ifelse( grepl(homopl_df$defining_mut, pattern ='^[N]:') , "yes", "no")
		# Get N mut coordinate only
		rgx_n <- regexpr('^[N]:[A-Z]{1}[0-9]{1,5}[A-Z]{1}', homopl_df$defining_mut)
		comm_coord <- rep(NA,length(homopl_df$defining_mut))
		comm_coord[rgx_n != -1] <- stringr::str_sub( regmatches(homopl_df$defining_mut, rgx_n) , 4, -2)
		comm_coord <- as.integer(comm_coord)
		comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("n_mut_coord")
		homopl_df <- cbind(homopl_df, comm_coord_df)
		data.table::setDT(homopl_df); data.table::setDT(regions_n)
		# Check if coordinate within boundaries of region of interest (Linker_domain)
		homopl_df[regions_n, on=c("n_mut_coord>=start", "n_mut_coord<=end"), n_mut_region_interest := region]
		
		homopl_df
	}
	
	# Annotate homoplasies in: (ii) different substitutions at the same site
	.annotate_diff_aa_mut_same_site <- function(homopl_df) {
		rgx_aadns <- regexpr('^[A-Za-z]:[A-Z]{1}[0-9]{1,5}[A-Z]{1}', homopl_df$defining_mut)
		comm_coord_aadns1 <- rep(NA,length(homopl_df$defining_mut))
		# All genomic regions + wt aa + position (without mutated aa)
		comm_coord_aadns1[rgx_aadns != -1] <- stringr::str_sub( regmatches(homopl_df$defining_mut, rgx_aadns) , 1, -2)
		comm_coord_aadns1_df <- as.data.frame(comm_coord_aadns1); colnames(comm_coord_aadns1_df) <- c("prot_coord_without_mut_site")
		homopl_df <- cbind(homopl_df, comm_coord_aadns1_df)
		homopl_df$prot_coord_without_mut_site <- stringr::str_sub( homopl_df$defining_mut, 1, -2)
		# Assign nodes with each prot_coord_without_mut_site
		homopl_aadns_nodes <- homopl_df %>% dplyr::group_by(prot_coord_without_mut_site) %>% dplyr::summarise(nodes_homopl = paste0(na.omit(node), collapse="|")) %>% dplyr::ungroup()
		homopl_aadns <- base::merge(homopl_df, homopl_aadns_nodes, by="prot_coord_without_mut_site")
		homopl_aadns$mut_site <- stringr::str_sub( homopl_aadns$defining_mut, -1, -1)
		# Filter repeated rows
		homopl_aadns <- homopl_aadns[!duplicated(homopl_aadns$defining_mut),]
		#homopl_aadns_freq <- homopl_aadns %>% dplyr::group_by(prot_coord_without_mut_site) %>% dplyr::summarise(Freq=dplyr::n(), diff_mutated_aa=paste0(na.omit(mut_site), collapse="|"))
		homopl_aadns_freq <- homopl_aadns %>% dplyr::group_by(prot_coord_without_mut_site) %>% dplyr::summarise(Freq_homopl=dplyr::n())
		homopl_aadns_final <- base::merge(homopl_aadns, homopl_aadns_freq, by="prot_coord_without_mut_site")
		homopl_aadns_final <- homopl_aadns_final %>% dplyr::select(nodes_homopl, prot_coord_without_mut_site, defining_mut, Freq_homopl.x, major_lineage) #diff_mutated_aa
		homopl_aadns_final <- homopl_aadns_final[!duplicated(homopl_aadns_final$prot_coord_without_mut_site),]
		homopl_aadns_final <- homopl_aadns_final[as.numeric(homopl_aadns_final$Freq_homopl.x) > 1,]
		homopl_aadns_final <- homopl_aadns_final %>% dplyr::select(nodes_homopl, defining_mut, Freq_homopl.x, major_lineage)
		
		homopl_aadns_final
	}
	
	# Annotate homoplasies in: (iii) moving window of neighbouring 3 residues
	.annotate_adjacent_muts_window_s3 <- function(homopl_df) {
		window_size <- 3
		
		prot_lengths <- utils::read.csv(system.file("extdata", "aa_length_prots.tsv", package="mlscluster"),sep="\t", header=T)
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
		prot_len_ids_df <- data.table::rbindlist(prot_len_ids)
		prot_windows_df <- data.table::rbindlist(prot_windows)
		prot_len_window_combined <- cbind(prot_len_ids_df, prot_windows_df)
		prot_len_window_combined <- data.frame(prot_len_window_combined, do.call(rbind,stringr::str_split(prot_len_window_combined$window,"-")))
		colnames(prot_len_window_combined) <- c("protein","window","window_pos1","window_pos2","window_pos3")
		prot_len_window_combined$window_pos1 <- as.integer(prot_len_window_combined$window_pos1)
		prot_len_window_combined$window_pos2 <- as.integer(prot_len_window_combined$window_pos2)
		prot_len_window_combined$window_pos3 <- as.integer(prot_len_window_combined$window_pos3)
		prot_len_window_combined <- prot_len_window_combined %>% dplyr::group_by(protein) %>% dplyr::mutate(window_id = rep(seq(dplyr::n()), each = window_size, length = dplyr::n())) %>% dplyr::ungroup()
		
		# Get protein name of homoplasies
		homopl_df$protein <- sub("\\:.*", "", homopl_df$defining_mut)
		homopl_df$protein <- toupper(homopl_df$prot)
		
		rgx_amw3 <- regexpr('^[A-Za-z]:[A-Z]{1}[0-9]{1,5}[A-Z]{1}', homopl_df$defining_mut)
		comm_coord <- rep(NA,length(homopl_df$defining_mut))
		# Mutation coordinate of homoplasies
		comm_coord[rgx_amw3 != -1] <- stringr::str_sub( regmatches(homopl_df$defining_mut, rgx_amw3) , 4, -2)
		comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("mut_coord")
		homopl_df <- cbind(homopl_df, comm_coord_df)
		homopl_df <- homopl_df[(homopl_df$protein != "SYNSNP"),]
		homopl_df$mut_coord <- as.integer(homopl_df$mut_coord)
		homopl_df <- stats::na.omit(homopl_df)
		
		homopl_df_merge <- homopl_df %>% dplyr::left_join(prot_len_window_combined, by="protein", relationship = "many-to-many")
		
		homopl_df_merge$window_pos1_match <- ifelse(homopl_df_merge$mut_coord == homopl_df_merge$window_pos1, TRUE, FALSE)
		homopl_df_merge$window_pos2_match <- ifelse(homopl_df_merge$mut_coord == homopl_df_merge$window_pos2, TRUE, FALSE)
		homopl_df_merge$window_pos3_match <- ifelse(homopl_df_merge$mut_coord == homopl_df_merge$window_pos3, TRUE, FALSE)
		
		homopl_df_matches <- homopl_df_merge[( (homopl_df_merge$window_pos1_match == TRUE) | (homopl_df_merge$window_pos2_match == TRUE) | (homopl_df_merge$window_pos3_match == TRUE)),]
		
		homopl_df_matches_nodes <- homopl_df_matches %>% dplyr::group_by(defining_mut) %>% dplyr::summarise(nodes_homopl = paste0(unique(node), collapse="|"), Freq_homopl=length(unique(node))) %>% dplyr::ungroup()
		homopl_df_matches_nodes <- homopl_df_matches_nodes[as.numeric(homopl_df_matches_nodes$Freq_homopl) > 1,]
		
		homopl_df_matches_neigh_mut <- homopl_df_matches %>% dplyr::inner_join(homopl_df_matches_nodes, by="defining_mut")
		homopl_df_matches_neigh_mut <- homopl_df_matches_neigh_mut %>% dplyr::select(nodes_homopl, protein, mut_coord, defining_mut, Freq_homopl.x, major_lineage, window, window_id)
		homopl_df_matches_neigh_mut <- homopl_df_matches_neigh_mut[base::order(homopl_df_matches_neigh_mut$protein, as.numeric(homopl_df_matches_neigh_mut$mut_coord)),]
		if(nrow(homopl_df_matches_neigh_mut) != 0) { #!is.null(homopl_df_matches_neigh_mut)
			homopl_df_matches_neigh_mut$window_id2 <- paste0(homopl_df_matches_neigh_mut$protein,":",homopl_df_matches_neigh_mut$window_id)
			homopl_df_matches_neigh_mut_final <- homopl_df_matches_neigh_mut[!duplicated(homopl_df_matches_neigh_mut[c("defining_mut","window_id2")]),] #4,8
			homopl_df_matches_neigh_mut_final <- homopl_df_matches_neigh_mut_final %>% dplyr::group_by(window_id2) %>% dplyr::filter(dplyr::n() != 1)
			homopl_df_matches_neigh_mut_final <- homopl_df_matches_neigh_mut_final[!duplicated(homopl_df_matches_neigh_mut_final$defining_mut),]
			homopl_df_matches_neigh_mut_final <- homopl_df_matches_neigh_mut_final %>% dplyr::select(nodes_homopl, defining_mut, Freq_homopl.x, major_lineage, window, window_id, window_id2)
			return(homopl_df_matches_neigh_mut_final)
		}else {
			homopl_df_matches_neigh_mut_final <- data.frame(nodes_homopl=c(NA),defining_mut=c(NA),Freq_homopl.x=c(NA),major_lineage=c(NA),window=c(NA),window_id=c(NA),window_id2=c(NA))
			return(homopl_df_matches_neigh_mut_final)
		}
	}
	
	# Flag know problematic sites listed here: https://github.com/W-L/ProblematicSites_SARS-CoV2
	.remove_known_problematic_sites <- function(homopl_df) {
		# Process problematic sites file
		probl_sites <- utils::read.csv(system.file("extdata", "problematic_sites_20221006.tsv", package="mlscluster"),sep="\t", header=T)
		# Split multiple potential mutated sites (separated by commas) for each site into unique row
		probl_sites <- probl_sites %>% dplyr::mutate(aa_alt = strsplit(as.character(aa_alt), ",")) %>% tidyr::unnest(aa_alt)
		probl_sites$aa_prot_site <- toupper( paste0(probl_sites$gene,":",probl_sites$aa_ref,probl_sites$aa_pos,probl_sites$aa_alt) )
		probl_sites <- probl_sites[probl_sites$filter == "mask",] #only removing sites flagged as 'mask' but not 'caution'
		
		# Remove sites listed as mask
		homopl_df <- dplyr::anti_join(homopl_df, probl_sites, by=c("defining_mut"="aa_prot_site"))
		homopl_df
	}
	
	# Return table with already described positively selected sites (https://observablehq.com/@spond/sars-cov-2-selection-countries) that should not appear here
	.sanity_check_positively_selected_sites <- function(homopl_df) {
		pos_sel_sites <- utils::read.csv(system.file("extdata", "positive_selected_sites_20221006.tsv", package="mlscluster"),sep="\t", header=T)
		pos_sel_sites$prot_site <- paste0(pos_sel_sites$protein,":",pos_sel_sites$site) 
		
		homopl_df$protein <- sub("\\:.*", "", homopl_df$defining_mut)
		homopl_df$protein <- toupper(homopl_df$prot)
		
		rgx_scpss <- regexpr(':[A-Z]{1}[0-9]{1,5}[A-Z]{1}', homopl_df$defining_mut)
		comm_coord <- rep(NA,length(homopl_df$defining_mut))
		# Mutation coordinate of homoplasies
		comm_coord[rgx_scpss != -1] <- stringr::str_sub( regmatches(homopl_df$defining_mut, rgx_scpss), 3, -2)
		comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("site")
		homopl_df <- cbind(homopl_df, comm_coord_df)
		homopl_df$prot_site <- paste0(homopl_df$protein,":",homopl_df$site)
		
		homopl_df_pss <- dplyr::inner_join(homopl_df, pos_sel_sites, by="prot_site")
		
		homopl_df_pss
	}
	
	# Homoplasies (all target nodes considered)
	homoplasies1_all_tgt_nodes <- data.frame(tgt_nodes); colnames(homoplasies1_all_tgt_nodes) <- c("node")
	#View(homoplasies1_all_tgt_nodes)
	#View(homoplasies[[1]])
	homoplasies1_all_tgt_nodes_df <- base::merge(homoplasies1_all_tgt_nodes, homoplasies, by="node", all.x=TRUE, all.y=FALSE)
	rm(homoplasies1_all_tgt_nodes); gc()
	homoplasies1_all_tgt_nodes_df <- stats::na.omit(homoplasies1_all_tgt_nodes_df)
	homoplasies1_all_tgt_nodes_df <- .remove_known_problematic_sites(homoplasies1_all_tgt_nodes_df)
	stats_df_homopl_nodes1 <- homoplasies1_all_tgt_nodes_df %>% dplyr::group_by(defining_mut) %>% dplyr::summarise(nodes_homopl = paste0(na.omit(node), collapse="|")) %>% dplyr::ungroup()
	stats_df_homopl1 <- base::merge(homoplasies1_all_tgt_nodes_df, stats_df_homopl_nodes1, by="defining_mut")
	rm(stats_df_homopl_nodes1); gc()
	stats_df_homopl1 <- stats_df_homopl1[!duplicated(stats_df_homopl1$defining_mut),]
	stats_df_homopl1 <- .annotate_s_homoplasies(stats_df_homopl1)
	stats_df_homopl1 <- .annotate_n_homoplasies(stats_df_homopl1)
	stats_df_homopl1 <- stats_df_homopl1 %>% dplyr::select(nodes_homopl, defining_mut, Freq_homopl, s_mut,	s_mut_region_interest, n_mut, n_mut_region_interest, major_lineage)
	stats_df_homopl1[is.na(stats_df_homopl1)] = ""
	stats_df_homopl1 <- stats_df_homopl1[base::order(stats_df_homopl1$s_mut, stats_df_homopl1$s_mut_region_interest, decreasing=TRUE),]
	stats_df_aadns1 <- .annotate_diff_aa_mut_same_site(homoplasies1_all_tgt_nodes_df)
	stats_df_aamw1 <- .annotate_adjacent_muts_window_s3(homoplasies1_all_tgt_nodes_df)
	rm(homoplasies1_all_tgt_nodes_df); gc()
	stats_df_scpss1 <- .sanity_check_positively_selected_sites(stats_df_homopl1)
	stats_df_scpss1 <- stats_df_scpss1[, c(1:8)]
	
	utils::write.csv(stats_df_homopl1, file=glue::glue('{output_dir}/homoplasy_DETAILS_all_tgt_nodes.csv'), quote=FALSE, row.names=FALSE)
	utils::write.csv(stats_df_aadns1, file=glue::glue('{output_dir}/homoplasy_SAME_SITE_muts_all_tgt_nodes.csv'), quote=FALSE, row.names=FALSE)
	utils::write.csv(stats_df_aamw1, file=glue::glue('{output_dir}/homoplasy_NEIGHBOUR_WINDOW_muts_all_tgt_nodes.csv'), quote=FALSE, row.names=FALSE)
	utils::write.csv(stats_df_scpss1, file=glue::glue('{output_dir}/homoplasy_KNOWN_POSITIVELY_SELECTED_all_tgt_nodes.csv'), quote=FALSE, row.names=FALSE)
	
	rm(stats_df_homopl1, stats_df_aadns1, stats_df_aamw1, stats_df_scpss1); gc()
	
	# Homoplasies (node detected by stat)
	stats_df_union$node <- as.integer(stats_df_union$node)
	homoplasies2_detect_df <- stats_df_union	%>% dplyr::left_join(homoplasies, by="node", relationship = "many-to-many")
	homoplasies2_detect_df <- stats::na.omit(homoplasies2_detect_df)
	homoplasies2_detect_df <- .remove_known_problematic_sites(homoplasies2_detect_df)
	stats_df_homopl_nodes2 <- homoplasies2_detect_df %>% dplyr::group_by(defining_mut) %>% dplyr::summarise(nodes_homopl = paste0(stats::na.omit(node), collapse="|")) %>% dplyr::ungroup()
	stats_df_homopl2 <- base::merge(homoplasies2_detect_df, stats_df_homopl_nodes2, by="defining_mut")
	rm(stats_df_homopl_nodes2); gc()
	stats_df_homopl2$Freq_homopl <- NULL
	stats_df_homopl2_freq <- data.table::setDT(stats_df_homopl2)[, .(Freq_homopl = .N), by = .(defining_mut)]
	stats_df_homopl2_freq_df <- base::merge(stats_df_homopl2, stats_df_homopl2_freq, by="defining_mut")
	gc(stats_df_homopl2, stats_df_homopl2_freq); gc()
	stats_df_homopl2_freq_df <- stats_df_homopl2_freq_df[as.numeric(stats_df_homopl2_freq_df$Freq_homopl) > 1,]
	stats_df_homopl2_freq_df <- .annotate_s_homoplasies(stats_df_homopl2_freq_df)
	stats_df_homopl2_freq_df <- stats_df_homopl2_freq_df[, c(2:3,1,16:17,18,20:21,4:15)]
	
	stats_df_homopl2_freq_df[is.na(stats_df_homopl2_freq_df)] = ""
	stats_df_homopl2_freq_df <- .annotate_n_homoplasies(stats_df_homopl2_freq_df)
	stats_df_homopl2_freq_df$is_clustered <- 1
	
	stats_df_aadns2 <- .annotate_diff_aa_mut_same_site(homoplasies2_detect_df)
	if(nrow(stats_df_aadns2) != 0) {
		stats_df_aadns2$is_clustered <- 1
	}
	stats_df_aamw2 <- .annotate_adjacent_muts_window_s3(homoplasies2_detect_df)
	if(nrow(stats_df_aamw2) != 0) {
		stats_df_aamw2$is_clustered <- 1
	}
	
	stats_df_scpss2 <- .sanity_check_positively_selected_sites(stats_df_homopl2_freq_df)
	stats_df_scpss2 <- stats_df_scpss2[, c(1:23)]
	stats_df_scpss2$is_clustered <- 1
	
	if(nrow(stats_df_homopl2_freq_df) != 0) {
		stats_df_homopl2_freq_df_cp <- stats_df_homopl2_freq_df %>% dplyr::select(nodes_homopl,defining_mut,Freq_homopl,is_clustered,s_mut_region_interest,n_mut_region_interest, major_lineage) #id_homopl_cat
		stats_df_homopl2_freq_df_cp <- stats_df_homopl2_freq_df_cp[!duplicated(stats_df_homopl2_freq_df_cp$defining_mut),]
	}else {
		stats_df_homopl2_freq_df_cp <- data.frame(nodes_homopl=c(NA),defining_mut=c(NA),Freq_homopl=c(NA),is_clustered=c(NA),s_mut_region_interest=c(NA),n_mut_region_interest=c(NA), major_lineage=c(NA))
	}
	
	if(nrow(stats_df_aadns2) != 0) {
		stats_df_aadns2_cp <- stats_df_aadns2 %>% dplyr::select(nodes_homopl,defining_mut,Freq_homopl.x,is_clustered, major_lineage) #id_homopl_cat
		stats_df_aadns2_cp <- stats_df_aadns2_cp[!duplicated(stats_df_aadns2_cp$defining_mut),]
		stats_df_aadns2_cp$s_mut_region_interest <- ""
		stats_df_aadns2_cp$n_mut_region_interest <- ""
		colnames(stats_df_aadns2_cp) <- c("nodes_homopl","defining_mut","Freq_homopl","is_clustered","major_lineage","s_mut_region_interest","n_mut_region_interest") #"id_homopl_cat"
	}
	
	stats_df_aamw2_cp <- data.frame()
	if(nrow(stats_df_aamw2) != 0) {
		stats_df_aamw2_cp <- stats_df_aamw2[, c(1:3,8,4)]
		stats_df_aamw2_cp <- stats_df_aamw2_cp[!duplicated(stats_df_aamw2_cp$defining_mut),]
		stats_df_aamw2_cp$s_mut_region_interest <- ""
		stats_df_aamw2_cp$n_mut_region_interest <- ""
		colnames(stats_df_aamw2_cp) <- c("nodes_homopl","defining_mut","Freq_homopl","is_clustered","major_lineage","s_mut_region_interest","n_mut_region_interest") #"id_homopl_cat"
	}else {
		stats_df_aamw2_cp <- data.frame(nodes_homopl=c(""),defining_mut=c(""),Freq_homopl.x=c(""),is_clustered=c(""),major_lineage=c(""),s_mut_region_interest=c(""),n_mut_region_interest=c(""))
	}
	
	utils::write.csv(stats_df_homopl2_freq_df, file=glue::glue('{output_dir}/homoplasy_DETAILS_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	utils::write.csv(stats_df_aadns2, file=glue::glue('{output_dir}/homoplasy_SAME_SITE_muts_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	utils::write.csv(stats_df_aamw2, file=glue::glue('{output_dir}/homoplasy_NEIGHBOUR_WINDOW_muts_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	utils::write.csv(stats_df_scpss2, file=glue::glue('{output_dir}/homoplasy_KNOWN_POSITIVELY_SELECTED_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	
	rm(stats_df_homopl2_freq_df, stats_df_aadns2, stats_df_aamw2); gc()
	
	# Homoplasies (node NOT detected by stat)
	homoplasies3_not_detect <- base::setdiff(tgt_nodes, stats_df_union$node)
	homoplasies3_not_detect <- data.frame(homoplasies3_not_detect); colnames(homoplasies3_not_detect) <- c("node")
	homoplasies3_not_detect_df <- homoplasies3_not_detect	%>% dplyr::left_join(homoplasies, by="node")
	rm(homoplasies3_not_detect); gc()
	homoplasies3_not_detect_df <- stats::na.omit(homoplasies3_not_detect_df)
	homoplasies3_not_detect_df <- .remove_known_problematic_sites(homoplasies3_not_detect_df)
	stats_df_homopl_nodes3 <- homoplasies3_not_detect_df %>% dplyr::group_by(defining_mut) %>% dplyr::summarise(nodes_homopl = paste0(stats::na.omit(node), collapse="|")) %>% dplyr::ungroup()
	stats_df_homopl3 <- base::merge(homoplasies3_not_detect_df, stats_df_homopl_nodes3, by="defining_mut")
	rm(stats_df_homopl_nodes3); gc()
	stats_df_homopl3$Freq_homopl <- NULL
	stats_df_homopl3_freq <- data.table::setDT(stats_df_homopl3)[, .(Freq_homopl = .N), by = .(defining_mut)]
	stats_df_homopl3_freq_df <- base::merge(stats_df_homopl3, stats_df_homopl3_freq, by="defining_mut")
	rm(stats_df_homopl3, stats_df_homopl3_freq); gc()
	stats_df_homopl3_freq_df <- stats_df_homopl3_freq_df[as.numeric(stats_df_homopl3_freq_df$Freq_homopl) > 1,]
	stats_df_homopl3_freq_df <- stats_df_homopl3_freq_df[!duplicated(stats_df_homopl3_freq_df$defining_mut),]
	stats_df_homopl3_freq_df <- .annotate_s_homoplasies(stats_df_homopl3_freq_df)
	stats_df_homopl3_freq_df <- .annotate_n_homoplasies(stats_df_homopl3_freq_df)
	stats_df_homopl3_freq_df <- stats_df_homopl3_freq_df %>% dplyr::select(nodes_homopl, defining_mut, Freq_homopl, s_mut, s_mut_region_interest, n_mut, n_mut_region_interest, major_lineage)
	stats_df_homopl3_freq_df[is.na(stats_df_homopl3_freq_df)] = ""
	stats_df_homopl3_freq_df <- stats_df_homopl3_freq_df[base::order(stats_df_homopl3_freq_df$s_mut, stats_df_homopl3_freq_df$s_mut_region_interest, decreasing=TRUE),]
	stats_df_homopl3_freq_df$is_clustered <- 0
	
	stats_df_aadns3 <- .annotate_diff_aa_mut_same_site(homoplasies3_not_detect_df)
	stats_df_aadns3$is_clustered <- 0
	
	stats_df_aamw3 <- .annotate_adjacent_muts_window_s3(homoplasies3_not_detect_df)
	stats_df_aamw3$is_clustered <- 0
	
	stats_df_scpss3 <- .sanity_check_positively_selected_sites(stats_df_homopl3_freq_df)
	stats_df_scpss3 <- stats_df_scpss3[, c(1:9)]
	stats_df_scpss3$is_clustered <- 0
	rm(homoplasies3_not_detect_df); gc()
	
	# Join friendly dfs for statistical analysis
	prot_lengths <- utils::read.csv(system.file("extdata", "aa_length_prots.tsv", package="mlscluster"),sep="\t", header=T)
	prot_lengths$protein <- toupper(prot_lengths$protein)
	
	if(nrow(stats_df_homopl3_freq_df) != 0) {
		stats_df_homopl3_freq_df_cp <- stats_df_homopl3_freq_df %>% dplyr::select(nodes_homopl,defining_mut,Freq_homopl,is_clustered,s_mut_region_interest,n_mut_region_interest, major_lineage) #id_homopl_cat
	}else {
		stats_df_homopl3_freq_df_cp <- data.frame(nodes_homopl=c(NA),defining_mut=c(NA),Freq_homopl=c(NA),is_clustered=c(NA),s_mut_region_interest=c(NA),n_mut_region_interest=c(NA), major_lineage=c(NA))
	}
	
	if(nrow(stats_df_aadns3) != 0) {
		stats_df_aadns3_cp <- stats_df_aadns3 %>% dplyr::select(nodes_homopl,defining_mut,Freq_homopl.x,is_clustered, major_lineage) #id_homopl_cat
		stats_df_aadns3_cp$s_mut_region_interest <- ""
		stats_df_aadns3_cp$n_mut_region_interest <- ""
		colnames(stats_df_aadns3_cp) <-  c("nodes_homopl","defining_mut","Freq_homopl","is_clustered","major_lineage","s_mut_region_interest","n_mut_region_interest") #"id_homopl_cat"
	}
	
	if(nrow(stats_df_aamw3) != 0) {
		stats_df_aamw3 <- as.data.frame(stats_df_aamw3)
		stats_df_aamw3_cp <- stats_df_aamw3[, c(1:3,8,4)] #6, 7:8
		stats_df_aamw3_cp$s_mut_region_interest <- ""
		stats_df_aamw3_cp$n_mut_region_interest <- ""
		colnames(stats_df_aamw3_cp) <- c("nodes_homopl","defining_mut","Freq_homopl","is_clustered","major_lineage","s_mut_region_interest","n_mut_region_interest") #"id_homopl_cat"
	}
	
	utils::write.csv(stats_df_homopl3_freq_df, file=glue::glue('{output_dir}/homoplasy_DETAILS_NOT_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	utils::write.csv(stats_df_aadns3, file=glue::glue('{output_dir}/homoplasy_SAME_SITE_muts_NOT_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	utils::write.csv(stats_df_aamw3, file=glue::glue('{output_dir}/homoplasy_NEIGHBOUR_WINDOW_muts_NOT_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	utils::write.csv(stats_df_scpss3, file=glue::glue('{output_dir}/homoplasy_KNOWN_POSITIVELY_SELECTED_NOT_detect_stats.csv'), quote=FALSE, row.names=FALSE)
	
	rm(stats_df_homopl3_freq_df, stats_df_aadns3, stats_df_aamw3); gc()
	
	clustered_df <- data.table::rbindlist(list(stats_df_homopl2_freq_df_cp, stats_df_homopl3_freq_df_cp, stats_df_aadns2_cp, stats_df_aadns3_cp, stats_df_aamw2_cp, stats_df_aamw3_cp), fill = TRUE)
	message("nrows clustered_df")
	message(nrow(clustered_df))
	clustered_df$protein <- sub("\\:.*", "", clustered_df$defining_mut)
	clustered_df$protein <- toupper(clustered_df$prot)
	clustered_df <- clustered_df %>% dplyr::inner_join(prot_lengths, by="protein")
	clustered_df$syn_non_syn <- ifelse(clustered_df$protein == "SYNSNP", "syn", "non-syn")
	clustered_df <- clustered_df %>% dplyr::select(defining_mut,major_lineage,Freq_homopl,nodes_homopl,is_clustered,syn_non_syn,s_mut_region_interest,n_mut_region_interest,protein,aa_length)
	clustered_df[clustered_df == ""] <- NA
	
	# Make known positively selection dfs compatible to create column stating whether this was independently found under selection
	stats_df_scpss2_cp <- stats_df_scpss2[!duplicated(stats_df_scpss2$defining_mut),]
	stats_df_scpss2_cp <- stats_df_scpss2_cp %>% dplyr::select(nodes_homopl,defining_mut,major_lineage,Freq_homopl,s_mut_region_interest,n_mut_region_interest,is_clustered) #id_homopl_cat,s_mut
	stats_df_scpss_join <- data.table::rbindlist(list(stats_df_scpss2_cp, stats_df_scpss3), fill=TRUE)
	
	clustered_df$indep_found_pos_selection <- ifelse(clustered_df$defining_mut %in% stats_df_scpss_join$defining_mut, "yes", "no")
	rm(stats_df_scpss2_cp, stats_df_scpss_join); gc()
	# remove duplicates coming from different homoplasy categories
	clustered_df <- clustered_df[!duplicated(clustered_df$defining_mut)]
	clustered_contingency_tab <- stats::xtabs(Freq_homopl ~ defining_mut+is_clustered, data=clustered_df) #id_homopl_cat
	clustered_contingency_tab <- as.data.frame(clustered_contingency_tab)
	clustered_contingency_tab <- clustered_contingency_tab[as.integer(clustered_contingency_tab$Freq) > 0,]
	utils::write.csv(clustered_df, file=glue::glue('{output_dir}/clustered_all_df.csv'), quote=FALSE, row.names=FALSE)
	utils::write.csv(clustered_contingency_tab, file=glue::glue('{output_dir}/clustered_contingency_df.csv'), quote=FALSE, row.names=FALSE)
	
	message(glue::glue("Nodes matching threshold {ifelse(threshold_keep_lower, '<', '>')} {quantile_threshold_ratio_sizes} quantile of ratio sizes: {nrow(stats_df_rs)}"))
	message(glue::glue("Nodes matching threshold {ifelse(threshold_keep_lower, '<', '>')} {quantile_threshold_ratio_persist_time} quantile of ratio persistence time: {nrow(stats_df_rpt)}"))
	message(glue::glue("Nodes matching threshold <= 0.05 of logistic regression p-value: {nrow(stats_df_lg_p)}"))
	message(glue::glue("Nodes matching all thresholds: {nrow(stats_df_intersect)}"))
	message(glue::glue("Nodes matching at least one threshold: {nrow(stats_df_union)}"))
	
	# Plot tree with defining mutations
	.cluster_tree <- function(tips) {
		
		tr <- ape::keep.tip(tre, tips)
		ggtr <- ggtree::ggtree(tr)
		mut_list <- strsplit( amd[match(tips, amd$sequence_name),]$mutations, split="\\|" )
		tab_mut <- table(do.call(c, mut_list))
		tab_mut <- tab_mut / Ntip(tr)
		shared_muts <- names(tab_mut)[tab_mut >= defining_mut_threshold]
		segregating_muts <- lapply(mut_list, function(x) base::setdiff(x, shared_muts) )
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
			aas <- aas[,base::order(sites)]
			aas[is.na(aas)] <- "X"
			sites <- sort(sites)
			aas2 <- phangorn::as.phyDat(aas)
			ap <- phangorn::ancestral.pars(tr, aas2, return="phyDat")
			ap1 <- as.character(ap)
			for(e in ape::postorder(tr)) {
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
		ggtr1 <- ggtr1 + ggplot2::geom_label( ggplot2::aes(x = branch, label = annot, size = 5)) #+ geom_tiplab(align=TRUE)
		ggtr2 <- ggtr1
		if(length(all_segregating) < 100) {
			suppressMessages( ggtr2 <- ggtree::gheatmap(ggtr1, as.data.frame(aas), width=0.66, offset=0.0005, colnames=FALSE, colnames_angle=-90, colnames_position="top", colnames_offset_y=-2) + ggplot2::theme(legend.position="none") )
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
	
	# summarise major lineage frequencies inside the node of interest
	.major_lineage_summary <- function(tips, max_rows=5) {
		if(is.null(tips)) {
			return("")
		}
		lineages <- amd$major_lineage[ match(tips, amd$sequence_name) ]
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
		
		for(i in 1:nrow(homopl_df)) {
			mut_freq <- outbreakinfo::getPrevalence(mutations = homopl_df[i,2], logInfo=FALSE)
			pl <- outbreakinfo::plotPrevalenceOverTime(mut_freq, title = glue::glue("Prevalence of {homopl_df[i,2]} worldwide"))
			ggplot2::ggsave(pl, file=glue::glue("{output_dir}/global_mut_freqs/{homopl_df[i,2]}_prevalence.jpg"), width=6, height=5)
		}
	}
	
	# Plot global prevalence of mutation for the homoplasies detected
	if(plot_global_mut_freq) {
		
		if(!dir.exists(glue::glue("{output_dir}/global_mut_freqs/")))
			suppressWarnings( dir.create(glue::glue("{output_dir}/global_mut_freqs/")) )
		
		.plot_global_mut_freqs_worldwide(stats_df_homopl1)
		.plot_global_mut_freqs_worldwide(stats_df_aamw1)
	}
	
	if(detailed_output) { # only plots node specific tables if change this flag to TRUE
		if(!is.null(desc_sisters)) {
			if(!dir.exists(glue::glue("{output_dir}/node_specific/")))
				suppressWarnings( dir.create(glue::glue("{output_dir}/node_specific/")) )
			
			message(glue::glue("Plotting trees + sequences + lineage + region summaries for nodes matching at least one of the thresholds..."))
			# Only plotting annotated trees, lineage and region summaries for nodes passing at least 1 of 3 criteria
			for(i in 1:nrow(stats_df_union)) { #length(tgt_nodes)
				tips <- list(desc_sisters[[1]][[i]], desc_sisters[[2]][[i]])
				tips <- unlist(tips)
				
				if(length(tips) > 0) {
					if(compute_tree_annots) { # tree annotations will only be computed if `write_node_specific_tables`==TRUE AND `compute_tree_annots`==TRUE
						ggtr <- .cluster_tree(tips)
						suppressMessages(ggplot2::ggsave(ggtr, file=glue::glue("{output_dir}/node_specific/{stats_df_union[i,1]}_clust_tree.pdf"),height=max(6, floor(length(tips) / 5 )), width = min(44, max(24, sqrt(length(tips)))), limitsize = FALSE  ))
					}
					tips_df <- amd[ amd$sequence_name %in% tips, ]
					utils::write.csv(tips_df, file=glue::glue('{output_dir}/node_specific/{stats_df_union[i,1]}_sequences.csv'), quote=FALSE, row.names=FALSE)
					lineage_summary <- .lineage_summary(tips=tips)
					utils::write.csv(lineage_summary, file=glue::glue('{output_dir}/node_specific/{stats_df_union[i,1]}_lineage_summary.csv'), quote=FALSE, row.names=FALSE)
					major_lineage_summary <- .major_lineage_summary(tips=tips)
					utils::write.csv(major_lineage_summary, file=glue::glue('{output_dir}/node_specific/{stats_df_union[i,1]}_MAJOR_lineage_summary.csv'), quote=FALSE, row.names=FALSE)
					region_summary <- .region_summary(tips=tips)
					utils::write.csv(region_summary, file=glue::glue('{output_dir}/node_specific/{stats_df_union[i,1]}_region_summary.csv'), quote=FALSE, row.names=FALSE)
				}
			}
		}
	}
	#message(paste0("Number of descendants vector length: ", length(ndesc)))
	
	invisible(clustered_contingency_tab)
}
