options(scipen=999)

#' Remove >99% quantile threshold of Freq_homopl since these sites are usually sequencing artifacts
#'
#' @keywords internal
#' @export
#' @noRd
remove_homopl_freq_outliers <- function(path_stats)  {
	clustered_dfs <- quant <- list()
	for(i in 1:length(path_thresholds)) {
		clustered_dfs[[i]] <- utils::read.csv(glue::glue("{path_stats}/{path_thresholds[i]}/clustered_all_df.csv"), header=T)
		clustered_dfs[[i]]$Freq_homopl <- as.integer(clustered_dfs[[i]]$Freq_homopl)
		quant[[i]] <- stats::quantile(clustered_dfs[[i]]$Freq_homopl, probs=seq(0,1,1/100))
		clustered_dfs[[i]] <- clustered_dfs[[i]][ clustered_dfs[[i]]$Freq_homopl <= unname(quant[[i]][names(quant[[i]]) == "99%"]), ]
		#print(max(clustered_dfs[[i]]$Freq_homopl))
	}
	return(clustered_dfs)
}

#' Include major_lineage column if it is missing and exclude recombinants
#'
#' @keywords internal
#' @export
#' @noRd
include_major_lineage_column <- function(md_df) {
	md_df$major_lineage = md_df$lineage
	md_df$major_lineage <- ifelse(grepl("^A[A-S]\\.|^AZ\\.", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^(AY\\.1$|AY\\.10|AY\\.11|AY\\.17|AY\\.19|AY\\.21|AY\\.120|AY\\.121|AY\\.122|AY\\.123|AY\\.124|AY\\.125|AY\\.126|AY\\.127|AY\\.128|AY\\.129|AY\\.13|AY\\.14|AY\\.16|AY\\.20|AY\\.22|AY\\.23|AY\\.24|AY\\.25|AY\\.26|AY\\.27|AY\\.28|AY\\.29|AY\\.3|AY\\.41|AY\\.42|AY\\.43|AY\\.44|AY\\.45|AY\\.46|AY\\.47|AY\\.5|AY\\.6|AY\\.7|AY\\.8|AY\\.9)", md_df$major_lineage), "Delta_other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^AY\\.4$", md_df$major_lineage), "Delta_AY.4.*", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^AY\\.4.", md_df$major_lineage), "Delta_AY.4.*", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^A\\.", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^A$", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("B\\.1\\.617", md_df$major_lineage), "Delta_other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("B\\.1\\.351", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B\\.1\\.1\\.7$", md_df$major_lineage), "Alpha_B.1.1.7", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B\\.1\\.177", md_df$major_lineage), "EU1_B.1.177", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B$", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B\\.1$", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B\\.1\\.1$", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^BA\\.1", md_df$major_lineage), "Omicron_BA.1.*", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^BA\\.2", md_df$major_lineage), "Omicron_BA.2.*", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^BA\\.3|^BA\\.4|^BA\\.5", md_df$major_lineage), "Omicron_BA.*", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^C|^W|^P\\.2", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^P\\.1", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^(AV\\.|AT\\.|BB\\.|D\\.|G\\.|L\\.|N\\.|M\\.|P\\.|Q\\.|R\\.|S\\.|U\\.|V\\.|Unassigned|Y\\.|Z\\.)", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B\\.1\\.1\\.[1-9]{1,3}", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B\\.1\\.[2-9]{1,3}", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^B\\.[2-9]{1,3}", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("B.18|B.1.104|B.1.105|B.1.110.2|B.1.110.3|B.1.111|B.1.115|B.1.117|B.1.12|B.1.128|B.1.13|B.1.131|B.1.142|B.1.146|B.1.147|B.1.149|B.1.153|B.1.157|B.1.160|B.1.160.|B.1.164|B.1.165|B.1.167|B.1.170|B.1.173|B.1.179|B.1.180|B.1.189|B.1.190|B.1.195|B.1.198|B.1.199", md_df$major_lineage), "Other", as.character(md_df$major_lineage))
	md_df$major_lineage <- ifelse(grepl("^X", md_df$major_lineage), "Recombinant", as.character(md_df$major_lineage))
	
	# UNCOMMENT FOR PLOTS
	# md_df_freq <- md_df %>% count(major_lineage, sort=T)
	# ggplot(data=md_df_freq, aes(x=n, y=major_lineage)) + geom_bar(stat="identity")
	# ggsave(glue::glue("rds/major_lineage_freqs.png"), width=8, height=5, dpi=600, bg="white")
	md_df_res = subset(md_df, major_lineage != "Recombinant")
	md_df_res <- md_df_res[!is.na(md_df_res$major_lineage),]
	# md_df_res_freq <- md_df_res %>% count(major_lineage, sort=T)
	# ggplot(data=md_df_res_freq, aes(x=n, y=major_lineage)) + geom_bar(stat="identity")
	# ggsave(glue::glue("rds/major_lineage_freqs_after_exclusions.png"), width=8, height=5, dpi=600, bg="white")
	
	return(md_df_res)
}
