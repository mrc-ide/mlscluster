library(glue)
library(stringr)
#load(file="rds/env.RData")
seq_muts_df_counts <- readRDS("rds/polymorphic_site_counts.rds")

c_path <- "results"
path_periods <- c(glue("{c_path}/01_sc2_root_to_dec2020"), glue("{c_path}/02_sc2_jan2021_to_may2021"),glue("{c_path}/03_sc2_jun2021_to_dec2021"), 
																		glue("{c_path}/04_sc2_jan2022_to_apr2022")) #,glue("{c_path}/05_sc2_whole_period")

thr_pref <- "threshold_quantile"
thr_path <- c(glue("{thr_pref}_0.25"), glue("{thr_pref}_0.5"), glue("{thr_pref}_0.75"), glue("{thr_pref}_1"), glue("{thr_pref}_2"),
														glue("{thr_pref}_3"), glue("{thr_pref}_4"), glue("{thr_pref}_5"), glue("{thr_pref}_10"), glue("{thr_pref}_25"), glue("{thr_pref}_50"))

sink(file = "rds/lm_output.txt")
clustered_dfs <- matrix(data=list(), nrow=length(path_periods), ncol=length(thr_path))
joined_mut_sites_clustered <- matrix(data=list(), nrow=length(path_periods), ncol=length(thr_path))
var_all <- mean_all <- poisson_models <- lr_models <- matrix(data=list(), nrow=length(path_periods), ncol=length(thr_path))
for(i in 1:length(path_periods)) { #length(path_periods)
	for(j in 1:length(thr_path)) {
		print("====================")
		print(glue("[{i},{j}]"))
		clustered_dfs[[i,j]] <- read.csv(glue("{path_periods[i]}/{thr_path[j]}/clustered_all_df.csv"), header=T)
		joined_mut_sites_clustered[[i,j]] <- clustered_dfs[[i,j]] %>% inner_join(seq_muts_df_counts, by=c("defining_mut"="mut_and_type"))
		# Extract genomic index position
		rgx_genomic_idx <- regexpr(':[A-Z]{1}[0-9]{1,5}[A-Z]{1}', joined_mut_sites_clustered[[i,j]]$defining_mut)
		comm_coord <- rep(NA,length(joined_mut_sites_clustered[[i,j]]$defining_mut))
		comm_coord[rgx_genomic_idx != -1] <- str_sub( regmatches(joined_mut_sites_clustered[[i,j]]$defining_mut, rgx_genomic_idx), 3, -2)
		comm_coord <- as.integer(comm_coord)
		comm_coord_df <- as.data.frame(comm_coord); colnames(comm_coord_df) <- c("genomic_index")
		joined_mut_sites_clustered[[i,j]] <- cbind(joined_mut_sites_clustered[[i,j]], comm_coord_df)
		#joined_mut_sites_clustered[[i,j]]$genomic_index <- as.integer(joined_mut_sites_clustered[[i,j]]$genomic_index)
		joined_mut_sites_clustered[[i,j]] <- joined_mut_sites_clustered[[i,j]] %>% select(defining_mut, Freq_homopl, is_clustered, syn_non_syn.x, protein.x, aa_length.x, indep_found_pos_selection, genomic_index) #s_mut_region_interest
		colnames(joined_mut_sites_clustered[[i,j]]) <- c("defining_mut","Freq_homopl","is_clustered","syn","protein","aa_length","indep_found_pos_selection","genomic_index") #s_mut_region_interest
		print(class(joined_mut_sites_clustered[[i,j]]$syn))
		#View(joined_mut_sites_clustered[[i,j]])
		# joined_mut_sites_clustered[[i,j]][joined_mut_sites_clustered[[i,j]] == "syn"] <- 1
		# joined_mut_sites_clustered[[i,j]][joined_mut_sites_clustered[[i,j]] == "non-syn"] <- 0
		joined_mut_sites_clustered[[i,j]]$syn <- ifelse(joined_mut_sites_clustered[[i,j]]$syn == as.character("non-syn"), 0, 1)
		joined_mut_sites_clustered[[i,j]]$indep_found_pos_selection <- ifelse(joined_mut_sites_clustered[[i,j]]$indep_found_pos_selection == as.character("yes"), 1, 0)
		
		# See if variance and mean approx equal: one of poisson regression assumptions
		# var_all[[i,j]] <- var(joined_mut_sites_clustered[[i,j]]$Freq, na.rm=T)
		# print(glue("Variance = {var_all[[i,j]]}"))
		# mean_all[[i,j]] <- mean(joined_mut_sites_clustered[[i,j]]$Freq, na.rm=T)
		# print(glue("Mean = {mean_all[[i,j]]}"))
		#hist(joined_mut_sites_clustered[[i,j]]$genomic_index)
		joined_mut_sites_clustered[[i,j]]$is_clustered <- as.factor(joined_mut_sites_clustered[[i,j]]$is_clustered)
		joined_mut_sites_clustered[[i,j]]$protein <- as.factor(joined_mut_sites_clustered[[i,j]]$protein)
		joined_mut_sites_clustered[[i,j]]$syn <- as.factor(joined_mut_sites_clustered[[i,j]]$syn)
		# TODO fix returning NAs: aa_length
		#joined_mut_sites_clustered[[i,j]]$s_mut_region_interest <- as.factor(joined_mut_sites_clustered[[i,j]]$s_mut_region_interest)
		joined_mut_sites_clustered[[i,j]]$indep_found_pos_selection <- as.factor(joined_mut_sites_clustered[[i,j]]$indep_found_pos_selection)
		poisson_models[[i,j]] <- glm(Freq_homopl ~ is_clustered + protein + aa_length + syn + indep_found_pos_selection, data=joined_mut_sites_clustered[[i,j]], na.action=na.omit, family=poisson(link = "log")) #defining_mut + genomic_index +s_mut_region_interest
		print("Summary poisson: ")
		print(summary(poisson_models[[i,j]]))
		lr_models[[i,j]] <- glm(is_clustered==1 ~ Freq_homopl + protein + aa_length + syn + indep_found_pos_selection, data=joined_mut_sites_clustered[[i,j]], na.action=na.omit, family=binomial(link="logit")) #defining_mut + genomic_index + s_mut_region_interest
		print("Summary logistic: ")
		print(summary(lr_models[[i,j]]))
		print("====================")
	}
}
sink(file = NULL)
#save.image("rds/env.RData")

# TODO below: diagnostics
# summary(poisson_models)
# exp(coef(poisson_models))
# predict(poisson_models, type="response") # predicted values
# residuals(poisson_models, type="deviance") # residuals
# plot(residuals(poisson_models, type="deviance")) # plot residuals

# exp(confint(poisson_models))
# To see which explanatory variables have an effect on response variable, look at the p-values
# If the Residual Deviance is greater than the degrees of freedom, then over-dispersion exists
# This means that the estimates are correct, but the standard errors (standard deviation) are wrong and unaccounted by the model.
