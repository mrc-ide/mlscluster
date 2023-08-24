## 1. Clustering algorithm run ##
sc2_md_curated <- readRDS("rds/sc2_md_curated.rds")
sc2_tre_curated <- readRDS("rds/sc2_tre_curated.rds")

source("R/utils.R")
source("R/mlsclust.R")
source("R/consts.R")
source("R/multiple_thresholds.R")
source("R/fdr.R")
source("R/selected_threshold.R")

NCPU <- 4
options(scipen=999)

###Testing only
#options(max.print=999999)

# Period 1: testing (~10 minutes: 6 for mlscluster and 4 for run_diff_thresholds)
start <- Sys.time()
res_p1 <- mlsclust(sc2_tre_curated, sc2_md_curated, min_descendants=10, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=as.Date("2019-12-30"), max_date=as.Date("2020-12-31"), branch_length_unit="days", defining_mut_threshold=0.75, root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995, ncpu=NCPU)
end <- Sys.time()
total_time_p1 <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time_p1,"mins"))

start <- Sys.time()
thr <- c(0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 10, 25)
quantl <- c(1/400, 1/200, 1/400, 1/100, 2/100, 3/100, 4/100, 5/100, 1/10, 1/4)
start <- Sys.time()
for(i in 1:length(thr)) {
	run_diff_thresholds(stats_df_unfilt=res_p1[[1]], tgt_nodes=res_p1[[2]], homoplasies=res_p1[[3]], output_dir=glue("results/TEST_01_sc2_root_to_dec2020/threshold_quantile_{thr[[i]]}/"), quantile_choice=quantl[i], quantile_threshold_ratio_sizes=thr[i], quantile_threshold_ratio_persist_time=thr[i], quantile_threshold_logit_growth=thr[i], threshold_keep_lower=TRUE, compute_tree_annots=FALSE, plot_global_mut_freq=FALSE, write_node_specific_tables=TRUE, desc_sisters=list(res_p1[[4]], res_p1[[5]]), amd=res_p1[[6]])
}
end <- Sys.time()
total_time_p1_t <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time_p1_t,"mins"))

# Period 2: june 2020 to mid-November 2021 -> testing ()
start <- Sys.time()
res_p2 <- mlsclust(sc2_tre_curated, sc2_md_curated, min_descendants=10, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=as.Date("2019-12-30"), max_date=as.Date("2021-11-15"), branch_length_unit="days", defining_mut_threshold=0.75, root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995, ncpu=NCPU)
end <- Sys.time()
total_time_p2 <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time_p2,"mins"))

start <- Sys.time()
for(i in 1:length(thr)) {
	run_diff_thresholds(stats_df_unfilt=res_p2[[1]], tgt_nodes=res_p2[[2]], homoplasies=res_p2[[3]], output_dir=glue("results/TEST_02_sc2_root_to_nov2021/threshold_quantile_{thr[[i]]}/"), quantile_choice=quantl[i], quantile_threshold_ratio_sizes=thr[i], quantile_threshold_ratio_persist_time=thr[i], quantile_threshold_logit_growth=thr[i], threshold_keep_lower=TRUE, compute_tree_annots=FALSE, plot_global_mut_freq=FALSE, write_node_specific_tables=TRUE, desc_sisters=list(res_p2[[4]], res_p2[[5]]), amd=res_p2[[6]])
}
end <- Sys.time()
total_time_p2_t <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time_p2_t,"mins"))

# Period 3: june 2020 to april 2022 -> testing ()
start <- Sys.time()
res_p3 <- mlsclust(sc2_tre_curated, sc2_md_curated, min_descendants=10, max_descendants=20e3, min_cluster_age_yrs=1/12, min_date=as.Date("2019-12-30"), max_date=as.Date("2022-04-30"), branch_length_unit="days", defining_mut_threshold=0.75, root_on_tip="Wuhan/WH04/2020", root_on_tip_sample_time=2019.995, ncpu=NCPU)
end <- Sys.time()
total_time_p3 <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time_p3,"mins"))

start <- Sys.time()
for(i in 1:length(thr)) {
	run_diff_thresholds(stats_df_unfilt=res_p3[[1]], tgt_nodes=res_p3[[2]], homoplasies=res_p3[[3]], output_dir=glue("results/TEST_03_sc2_whole_period/threshold_quantile_{thr[[i]]}/"), quantile_choice=quantl[i], quantile_threshold_ratio_sizes=thr[i], quantile_threshold_ratio_persist_time=thr[i], quantile_threshold_logit_growth=thr[i], threshold_keep_lower=TRUE, compute_tree_annots=FALSE, plot_global_mut_freq=FALSE, write_node_specific_tables=TRUE, desc_sisters=list(res_p3[[4]], res_p3[[5]]), amd=res_p3[[6]])
}
end <- Sys.time()
total_time_p2_t <- as.numeric (end - start, units = "mins")
message(paste("Total time elapsed:",total_time_p2_t,"mins"))
## End clustering algorithm run (test with 2020 data only) ##

## 2. Run for multiple thresholds ##
print("=========")
print("PERIOD 2 (PRE-OMICRON)")
print("=========")
# FOR PERIOD 2
# Expected paths containing results for each quantile threshold
table_names_periods <- c(2)
table_names_thresholds <- c(0.25,0.5,0.75,1,2,3,4,5,10,25)
table_combs <- do.call(paste, c(expand.grid(table_names_periods, table_names_thresholds), sep="_"))
out_folder <- "period2"
PERIOD_INTEREST <- 2
major_lineages <- c("EU1_B.1.177","Alpha_B.1.1.7","Delta_AY.4.*", "Delta_other","Other") #"Beta_B.1.351","Gamma_P.1"

system(glue("mkdir -p stat_results/plots_paper/"))
rmult_p2 <- stats_multiple_thresholds("results/02_sc2_root_to_nov2021", rm_freq_outliers=TRUE, "period2","FileS2.txt")
plot_tbc_stats("results/02_sc2_root_to_nov2021",out_folder)
stacked_nsites_genomic_region_mult_thresholds("stat_results/period2/genomewide_plot_non-syn/")
# Figure S10. spike for pre-Omicron period
fs10 <- ggarrange(rmult_p2[[2]][[2]][[10]] + theme(axis.text.x=element_blank(), axis.title.x=element_blank()), rmult_p2[[2]][[1]][[10]], align='h', nrow=2, ncol=1, labels=c('A', 'B'), legend="top", common.legend = FALSE) #common.legend = T
ggsave("stat_results/plots_paper/FigS10.png", units="px", width=2200, height=2000, dpi=300, bg="white")

# Figure S11. orf7a for pre-Omicron period
fs11 <- ggarrange(rmult_p2[[2]][[2]][[8]] + theme(axis.text.x=element_blank(), axis.title.x=element_blank()), rmult_p2[[2]][[1]][[8]], align='h', nrow=2, ncol=1, labels=c('A', 'B'), legend="top", common.legend = FALSE)
ggsave("stat_results/plots_paper/FigS11.png", plot=fs11, units="px", width=2200, height=2000, dpi=300, bg="white")

# Figure S12. orf3a for pre-Omicron period
fs12 <- ggarrange(rmult_p2[[2]][[2]][[6]] + theme(axis.text.x=element_blank(), axis.title.x=element_blank()), rmult_p2[[2]][[1]][[6]], align='h', nrow=2, ncol=1, labels=c('A', 'B'), legend="top", common.legend = FALSE)
ggsave("stat_results/plots_paper/FigS12.png", plot=fs12, units="px", width=2200, height=2000, dpi=300, bg="white")

# Figure S13. orf8 for pre-Omicron period
fs13 <- ggarrange(rmult_p2[[2]][[2]][[9]] + theme(axis.text.x=element_blank(), axis.title.x=element_blank()), rmult_p2[[2]][[1]][[9]], align='h', nrow=2, ncol=1, labels=c('A', 'B'), legend="top", common.legend = FALSE)
ggsave("stat_results/plots_paper/FigS13.png", plot=fs13, units="px", width=2200, height=2000, dpi=300, bg="white")

fs14 <- ggarrange(rmult_p2[[3]][[2]] + theme(axis.text.x=element_blank(), axis.title.x=element_blank()), rmult_p2[[3]][[1]], align='h', nrow=2, ncol=1, labels=c('A', 'B'), legend="top", common.legend = FALSE)
ggsave("stat_results/plots_paper/FigS14.png", plot=fs14, units="px", width=2200, height=2000, dpi=300, bg="white")

print("=========")
print("PERIOD 3 (INCLUDING OMICRON)")
print("=========")
# FOR PERIOD 3
table_names_periods <- c(3)
table_names_thresholds <- c(0.25,0.5,0.75,1,2,3,4,5,10,25)
table_combs <- do.call(paste, c(expand.grid(table_names_periods, table_names_thresholds), sep="_"))
out_folder <- "period3"
PERIOD_INTEREST <- 3
major_lineages <- c("EU1_B.1.177","Alpha_B.1.1.7","Delta_AY.4.*", "Delta_other","Omicron_BA.1.*","Omicron_BA.2.*","Other") #"Beta_B.1.351","Gamma_P.1","Omicron_BA.*",

rmult_p3 <- stats_multiple_thresholds("results/03_sc2_whole_period", rm_freq_outliers=TRUE, "period3","FileS3.txt")
s1 <- plot_tbc_stats("results/03_sc2_whole_period",out_folder)
ggsave("stat_results/plots_paper/FigS1.png", plot=s1, units="px", width=1500, height=2000, dpi=300, bg="white")

# Figs S3-S8
#major_lineages_supp_order <- c("Omicron_BA.1.*", "Delta_AY.4.*","Alpha_B.1.1.7","Delta_other","Other","Omicron_BA.2.*")
major_lineages_supp_indices <- c(5,3,2,4,7,6)
j <- 4
for(i in major_lineages_supp_indices) {
	print(i)
	ggsave(glue("stat_results/plots_paper/FigS{j}_{major_lineages[i]}.png"), plot=rmult_p3[[1]][[i]]+theme(legend.key.size = unit(0.5, 'cm')), units="px", width=2000, height=2200, dpi=300, bg="white")
	j <- j+1
}

sngrt_r3 <- stacked_nsites_genomic_region_mult_thresholds("stat_results/period3/genomewide_plot_non-syn/")

#saveRDS(sngrt_r3, "rds/sngrt_r3.rds")
## End run for multiple thresholds ##

## 3. Run FDR ##
# Extracting mutations from the entire dataset of >1.2 million sequences from England
cog_md_muts_p2 <- extract_muts_period(as.Date("2021-11-15"), "period2")
cog_md_muts_p3 <- extract_muts_period(as.Date("2022-04-30"), "period3")

# Computing codon positions of mutations
muts_match_codons_p2 <- compute_muts_match_codons(cog_md_muts_p2[[1]], cog_md_muts_p2[[2]])
muts_match_codons_p3 <- compute_muts_match_codons(cog_md_muts_p3[[1]], cog_md_muts_p3[[2]])

fdr_plots_p2 <- fdr_mlsclust("results/02_sc2_root_to_nov2021", muts_match_codons_p2[[1]], muts_match_codons_p2[[2]], muts_match_codons_p2[[3]], "period2")
fdr_plots_p3 <- fdr_mlsclust("results/03_sc2_whole_period", muts_match_codons_p3[[1]], muts_match_codons_p3[[2]], muts_match_codons_p3[[3]], "period3")

rem_labels <- function(ggplot_obj_list) {
	for(i in 1:length(ggplot_obj_list)) {
		if(i <= 20) {
			ggplot_obj_list[[i]] <- ggplot_obj_list[[i]] + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
		}
		if(!(i %in% seq(1,25,5))) {
			ggplot_obj_list[[i]] <- ggplot_obj_list[[i]] + theme(axis.text.y=element_blank(), axis.title.y=element_blank())
		}
	}
	ggplot_obj_list
}

fdr_plots_p2 <- rem_labels(fdr_plots_p2)
ggarrange(plotlist=fdr_plots_p2, ncol=5, nrow=5, common.legend = T)
ggsave(glue("stat_results/period2/fdr_codons/fdr_all_regions_p2.pdf"), width=10, height=10, dpi=600, bg="white")

fdr_plots_p3 <- rem_labels(fdr_plots_p3)
ggarrange(plotlist=fdr_plots_p3, ncol=5, nrow=5, common.legend = T)
ggsave(glue("stat_results/period3/fdr_codons/fdr_all_regions_p3.pdf"), width=10, height=10, dpi=600, bg="white")
## End run FDR ##

## 4. Run for the selected threshold (2% since FDR ~10% and epsilon ~2.5%) ##
print("=========")
print("PERIOD 2 (PRE-OMICRON)")
print("=========")
# FOR PERIOD 2
# Expected paths containing results for each quantile threshold
out_folder <- "period2_thr2"
PERIOD_INTEREST <- 2
THRESHOLD_INTEREST <- 2
major_lineages <- c("EU1_B.1.177","Alpha_B.1.1.7","Delta_AY.4.*", "Delta_other","Other")

r2 <- stats_selected_threshold("results/02_sc2_root_to_nov2021", rm_freq_outliers=TRUE, thr_index=5, out_folder) # NOTE: thr_index = 5 -> 2%
r2[[4]]
ggsave(glue("stat_results/plots_paper/FigS3.png"), units="px", width=2200, height=2500, dpi=300, bg="white")


print("=========")
print("PERIOD 3 (INCLUDING OMICRON)")
print("=========")
# FOR PERIOD 3
thr_pref <- "threshold_quantile"
out_folder <- "period3_thr2"
PERIOD_INTEREST <- 3
THRESHOLD_INTEREST <- 2
major_lineages <- c("EU1_B.1.177","Alpha_B.1.1.7","Delta_AY.4.*", "Delta_other","Omicron_BA.1.*","Omicron_BA.2.*","Other")

r3 <- stats_selected_threshold("results/03_sc2_whole_period", rm_freq_outliers=TRUE, thr_index=5, out_folder)

# Figure 3: stacked bar counts of unique homoplasies total & norm
r2[[2]][[1]] <- r2[[2]][[1]] + scale_x_continuous(limits=c(0, 110))
r3[[2]][[1]] <- r3[[2]][[1]] + scale_x_continuous(limits=c(0, 110))
rm_titles_all <- theme(axis.title.x=element_blank(), axis.title.y=element_blank())
rm_titles_y <- theme(axis.title.y=element_blank())
rm_axis_text_x <- theme(axis.text.x=element_blank())
my_legend <- get_legend(r3[[2]][[1]])
#ggarrange(legend_1, legend_2, legend_3, nrow=3)
f3 <- ggarrange(r2[[2]][[1]] + rm_titles_all + rm_axis_text_x + marg, r2[[2]][[2]] + rm_titles_all + rm_axis_text_x + marg, r3[[2]][[1]] + rm_titles_y + marg, r3[[2]][[2]] + rm_titles_y + marg, align='h', nrow=2, ncol=2, labels=c('A', 'B','C','D'), legend.grob=my_legend, legend="right") #common.legend = T
annotate_figure(f3, left = text_grob("Genomic region", rot = 90, vjust = 1, size=8))
ggsave("stat_results/plots_paper/Fig3.eps", device=cairo_ps, units="px", width=2100, height=1850, dpi=300, bg="white")

# Figure 4: Omicron results. A. comparison against positive selection (HyPhy) and B. identified sites by major lineage
r3[[4]]
ggsave("stat_results/plots_paper/Fig4.eps", device=cairo_ps, units="px", width=2100, height=1850, dpi=300, bg="white")

# Figure 5: spike-wide frequency of TFP-homoplasies & number of identified sites for each genomic position stacked by threshold
# sngrt_r3 variable from analysis_mult_thresholds.R
sngrt_r3 <- readRDS("rds/sngrt_r3.rds")
f5_configs <- theme(axis.title.x=element_text(size=8), axis.title.y=element_text(size=8), axis.text.y=element_text(size=6, color="black"), axis.text.x=element_text(size=7, color="black"))
ggarrange(r3[[3]] + f5_configs, sngrt_r3 + f5_configs, nrow=2, ncol=1, labels=c('A', 'B'), legend="right") #common.legend = T
ggsave("stat_results/plots_paper/Fig5.eps", units="px", device=cairo_ps, width=2000, height=1750, dpi=300, bg="white")
## End run for the selected threshold ##