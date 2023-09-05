# Threshold prefixes and paths
thr_pref <- "threshold_quantile"
thr <- c(0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 10, 25)
path_thresholds <- c(glue::glue("{thr_pref}_0.25"), glue::glue("{thr_pref}_0.5"), glue::glue("{thr_pref}_0.75"), glue::glue("{thr_pref}_1"), glue::glue("{thr_pref}_2"),
																					glue::glue("{thr_pref}_3"), glue::glue("{thr_pref}_4"), glue::glue("{thr_pref}_5"), glue::glue("{thr_pref}_10"), glue::glue("{thr_pref}_25"))

options(scipen=999)

# Palette for thresholds
pal_thresholds <- c("0.25"="#1984c5", "0.5"="#22a7f0", "0.75"="#63bff0", "1"="#a7d5ed", "2"="#cbfeff", "3"="#ffcccb", "4"="#e1a692", "5"="#de6e56", "10"="#e14b31", "25"="#c23728") # removed e2e2e2

# Palette for major_lineages
pal_lineages <- c("Alpha_B.1.1.7"="#fd7f6f", "Delta_AY.4.*"="#7eb0d5", "Delta_other"="#b2e061", "EU1_B.1.177"="#bd7ebe", "Omicron_BA.1.*"="#ffb55a", "Omicron_BA.2.*"="#ffcc00", "Other"="#fdcce5")

syn_p <- "SYNSNP:"
# Factors levels for syn mutations
lvls_syn_proteins <- c("NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NPS12","NSP13","NSP14","NSP15",
																							"NSP16","S","ORF3A","E","M","ORF6","ORF7A","ORF7B","ORF8","N","ORF10")

# Factors levels for non-syn mutations (with S and N annotations)
lvls_nonsyn_proteins_annots <- c("NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NSP12","NSP13","NSP14","NSP15","NSP16",
																																	"S","S:NTD","S:RBD","S:FCS","ORF3A","E","M","ORF6","ORF7A","ORF7B","ORF8","N","N:Linker_domain","ORF10")

# Factors levels for non-syn mutations (without S and N annotations)
lvls_nonsyn_proteins <- c("NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NSP12","NSP13","NSP14","NSP15","NSP16",
																										"S","ORF3A","E","M","ORF6","ORF7A","ORF7B","ORF8","N","ORF10")

# Factor levels for bubble plots: no NSPs
lvls_bubble <- c("ORF1AB","S","ORF3A","E","M","ORF6","ORF7A","ORF7B","ORF8","N","ORF10")


pal_id <- c("mlscluster"="#7c9885","Intersection"="#e1ad01", "HyPhy"="#033f63")

# Load config csv's
pos_sel_sites <- utils::read.csv(system.file("extdata", "positive_selected_sites_20221006.tsv", package="mlscluster"),sep="\t", header=T)
pos_sel_sites$prot_site <- paste0(pos_sel_sites$protein,":",pos_sel_sites$site)
annot_gen_range_positions <- utils::read.csv(system.file("extdata", "aa_ranges_stat_model.tsv", package="mlscluster"),sep="\t", header=T)
annot_gen_length_positions <- utils::read.csv(system.file("extdata", "aa_length_stat_model.tsv", package="mlscluster"),sep="\t", header=T)
syn_ranges_df <- utils::read.csv(system.file("extdata", "genomic_ranges_plot_syn.tsv", package="mlscluster"),sep="\t", header=T)

marg <- ggplot2::theme(plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, "cm"))