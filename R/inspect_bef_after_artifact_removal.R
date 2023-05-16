library(dplyr)

before <- read.csv("results/ok_bef_artifact_removal/02_sc2_root_to_nov2021/threshold_quantile_25/clustered_all_df.csv", header=T)
before <- before[before$is_clustered == 1,]
nrow(before[before$Freq_homopl >= 3,])

after <- read.csv("results/02_sc2_root_to_nov2021/threshold_quantile_25/clustered_all_df.csv", header=T)
after <- after[after$is_clustered == 1,]
nrow(after[after$Freq_homopl >= 3,])

aj <- before %>% anti_join(after, by="defining_mut")
#aj2 <- dplyr::setdiff(before, after)
