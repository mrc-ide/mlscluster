library(glue)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(viridis)

# Plot distribution of sequences included (fig. 2)
sc2_md_curated2 <- readRDS("rds/sc2_md_curated.rds")
sc2_md_curated2 <- include_major_lineage_column(sc2_md_curated2)
sc2_md_curated2$major_lineage <- ifelse(grepl("\\bOmicron_BA\\.\\*\\b", sc2_md_curated2$major_lineage), "Other", as.character(sc2_md_curated2$major_lineage))
#sc2_tre_curated2 <- readRDS("rds/sc2_tre_curated.rds")

# Palette for major_lineages
pal_lineages <- c("Alpha_B.1.1.7"="#fd7f6f", "Delta_AY.4.*"="#7eb0d5", "Delta_other"="#b2e061", "EU1_B.1.177"="#bd7ebe", "Omicron_BA.1.*"="#ffb55a", "Omicron_BA.2.*"="#ffee65", "Other"="#fdcce5") #removed "#beb9db", "#8bd3c7"

sc2_md_curated2$month <-as.Date(cut(sc2_md_curated2$sample_date, breaks = "1 month"))

df_lineages_per_month <- sc2_md_curated2 %>% group_by(month, major_lineage) %>% summarise(count_lin_month=n())
df_lineages_per_month <- df_lineages_per_month[!is.na(df_lineages_per_month$major_lineage),]

system("mkdir -p stat_results/plots_paper/england_maps/")

# Stacked bar chart of sequence distribution by sample date and lineage
pl1 <- df_lineages_per_month  %>% mutate(date=as.POSIXct(month)) %>% 
	ggplot(data=df_lineages_per_month, mapping=aes(fill=major_lineage, y=count_lin_month, x=month)) +  #mapping = aes(x = month, fill=major_lineage)
	geom_bar(position="stack", stat="identity", color="black")+ #width=25,
	scale_x_date(date_labels = "%b\n%Y", date_breaks = "2 months", limits = as.Date(c("2020-05-01","2022-05-01")), expand = c(0, 0))+
	scale_fill_manual(values=pal_lineages, name="Major lineage") + theme_minimal() + theme(legend.key.size = unit(0.3, 'cm'), legend.position="right") +
	xlab("Sample dates") + ylab("# sequences included")
ggsave("stat_results/plots_paper/england_maps/sequences_time_lineage.png", width=8, height=6.5, dpi=600, bg="white")

# Map of spatiotemporal distribution of samples in England across 4 waves and all period (1. June-Dec 2020 / Other+EU1, 2. Jan-May 2021 / Alpha, 3. June-Dec 2021 / Deltas, 4. Dec 2020-Apr 2021 / Omicrons, 5. All ~2 year period)
library(rgeos)
library(rgdal)
library(maptools)

# Remove UTLA duplicates and match with LTLA to allow plotting
sc2_md_curated2 <- sc2_md_curated2[!is.na(sc2_md_curated2$region),] #1275669
sc2_md_curated2 <- sc2_md_curated2[sc2_md_curated2$region != "",] # 20885 where region == ""
sort(unique(sc2_md_curated2$region)) #285 unique regions

sc2_md_curated2$region <- gsub("_", " ", sc2_md_curated2$region)
#sc2_md_curated2 <- sc2_md_curated2[sc2_md_curated2$region != "LONDON",] # 1 obs
#sc2_md_curated2 <- sc2_md_curated2[sc2_md_curated2$region != "HILLINGDON",] # 3 obs
# sc2_md_curated2 <- sc2_md_curated2[ sc2_md_curated2$region=="LONDON" ] <- "INNER LONDON"
# sc2_md_curated2 <- sc2_md_curated2[ sc2_md_curated2$region== "HILLINGDON" ] <- "OUTER LONDON"
sort(unique(sc2_md_curated2$region)) #265 unique regions

#lookup_ltla_utla_england <- read.csv("data/ltla_to_utla_lookup_england.csv", header=T) 
lookup_ltla_utla_adm2 <- read.csv("data/LAD_UTLA_adm2.csv", header=T) 
lookup_ltla_utla_adm2$adm2 <- gsub("_", " ", lookup_ltla_utla_adm2$adm2)
lookup_ltla_utla_adm2 <- lookup_ltla_utla_adm2[grepl("E", lookup_ltla_utla_adm2$UTLA_code),] # keep only England
sort(unique(lookup_ltla_utla_adm2$LAD_code)) #317 unique LTLAs

library(sp)
library(maptools)
library(raster)
#library(leaflet) #install.packages("leaflet")
library(broom)

shp <- getData('GADM', country='GBR', level = 2, path="data/") 
head(shp, 3)
shp_england <- shp[shp$NAME_1 == "England",]
sort(unique(shp_england$NAME_2))
shp_england$NAME_2 <- toupper(shp_england$NAME_2)
plot(shp_england)
shp_england_df <- tidy(shp_england, region="NAME_2") # convert to df with lat and long
length(unique(shp_england_df$id)) #118 adm2 regions from GADM

rm(shp, shp_england)
#print(length(unique(shp_england$id)) == length(unique(lookup_ltla_utla_england$LTLA22CD)))

sc2_md_curated3 <- sc2_md_curated2
#sc2_md_curated3$region[sc2_md_curated3$region == "GREATER MANCHESTER"] <- "MANCHESTER"
sc2_md_curated3$region[sc2_md_curated3$region == "COUNTY DURHAM"] <- "DURHAM"
sc2_md_curated3$region[sc2_md_curated3$region == "BERKSHIRE"] <- "WEST BERKSHIRE"
# extract matches between COG md region and adm2
sc2_md_curated3_test <- sc2_md_curated3 %>% inner_join(lookup_ltla_utla_adm2, by=c("region"="adm2"), multiple="first")
unique(sc2_md_curated3_test$region) #118 unique regions
#sc2_md_curated2_test$adm2 <- sc2_md_curated2_test$region
sc2_md_curated3_test <- sc2_md_curated3_test %>% dplyr::select(sequence_name, sample_date, lineage, major_lineage, region, mutations, month)
# extract NON-matches between COG md region and adm2
sc2_md_curated3_test_not <- sc2_md_curated3 %>% anti_join(lookup_ltla_utla_adm2, by=c("region"="adm2")) #multiple="first
sc2_md_curated3_test_not_c <- sc2_md_curated3_test_not %>% group_by(region) %>% summarise(n=n())
#sc2_md_curated2_test_not$region <- ifelse(sc2_md_curated2_test_not$region %in% lookup_ltla_utla_adm2$LAD_name, yes=, no=)
#%>% inner_join(lookup_ltla_utla_adm2, by=c("region"="LAD_name"), multiple="first")
# Extract other matches with LTLAs (commenting because this inflates DURHAM)
sc2_md_curated3_test_match_ltla <- sc2_md_curated3_test_not %>% inner_join(lookup_ltla_utla_adm2, by=c("region"="LAD_name"), multiple="first")
sc2_md_curated3_test_match_ltla$region <- sc2_md_curated3_test_match_ltla$adm2
unique(sc2_md_curated3_test_match_ltla$region)
sc2_md_curated3_test_match_ltla <- sc2_md_curated3_test_match_ltla %>% dplyr::select(sequence_name, sample_date, lineage, major_lineage, region, mutations, month)
sc2_md_curated3_test_match_ltla_c <- sc2_md_curated3_test_match_ltla %>% group_by(region) %>% summarise(n=n())
# Extract other matches with UTLAs
#sc2_md_curated2_test_match_utla <- sc2_md_curated2_test_not %>% inner_join(lookup_ltla_utla_adm2, by=c("region"="UTLA_name"), multiple="first")
#sc2_md_curated2_test_match_utla$region <- sc2_md_curated2_test_match_utla$adm2
#unique(sc2_md_curated2_test_match_utla$region)
#sc2_md_curated2_test_match_utla <- sc2_md_curated2_test_match_utla %>% dplyr::select(sequence_name, sample_date, lineage, major_lineage, region, mutations, month)
#sc2_md_curated3_all <- rbind(sc2_md_curated3_test, sc2_md_curated3_test_match_ltla) #sc2_md_curated2_test_match_utla
sc2_md_curated3_all <- sc2_md_curated3_test
sc2_md_curated3_all <- sc2_md_curated3_all[!duplicated(sc2_md_curated3_all$sequence_name), ]
#length(unique(sc2_md_curated2_all$region)) #119 (LIVERPOOL included here but not in GADM, ok since only 1639 obs)
#length(unique(shp_england$NAME_2))
#nrow(sc2_md_curated2) - nrow(sc2_md_curated2_all) #44513 without regions after matches + 20885 blank = 65398

london_ltla_to_adm2 <- lookup_ltla_utla_adm2$LAD_name[lookup_ltla_utla_adm2$adm2 == "GREATER LONDON"]

# Loading LTLA and UTLA data to plot proportion of cases instead of absolute genome counts
calculate_cases_england_regions <- function(file, period=5) {
	cases_region <- read.csv(file)
	cases_region$areaName <- toupper(cases_region$areaName)
	cases_region$areaName <- gsub("_", " ", cases_region$areaName)
	cases_region[ cases_region == "BRISTOL, CITY OF" ] <- "BRISTOL"
	cases_region[ cases_region == "COUNTY DURHAM" ] <- "DURHAM"
	cases_region[ cases_region == "KINGSTON UPON HULL, CITY OF" ] <- "KINGSTON UPON HULL"
	cases_region[ cases_region == "ST. HELENS" ] <- "SAINT HELENS"
	cases_region[ cases_region == "BEDFORD" ] <- "BEDFORDSHIRE"
	cases_region[ cases_region == "HEREFORDSHIRE, COUNTY OF" ] <- "HEREFORDSHIRE"
	cases_region[ cases_region == "HACKNEY AND CITY OF LONDON" ] <- "GREATER LONDON"
	cases_region$areaName <- ifelse(cases_region$areaName %in% london_ltla_to_adm2, yes="GREATER LONDON", no=as.character(cases_region$areaName))
	print("Unique areaNames")
	print(length(unique(cases_region$areaName))) #380 areaNames
	cases_region$month <- as.Date(cut(as.Date(cases_region$date), breaks = "1 month"))
	
	#cases_region <- cases_region %>% filter(between(month, as.Date('2020-06-01'), as.Date('2022-04-01'))) %>% group_by(areaName) %>% summarise(cases=sum(newCasesBySpecimenDate))
	
	if(period==1) { #June-Dec 2020 / Other+EU1
		cases_region <- cases_region %>% filter(between(month, as.Date('2020-06-01'), as.Date('2020-12-01'))) %>% group_by(areaName) %>% summarise(cases=sum(newCasesBySpecimenDate))
	}else if(period==2) { #Jan-May 2021 / Alpha
		cases_region <- cases_region %>% filter(between(month, as.Date('2021-01-01'), as.Date('2021-05-01'))) %>% group_by(areaName) %>% summarise(cases=sum(newCasesBySpecimenDate))
	}else if(period==3) { #June-Dec 2021  / Deltas
		cases_region <- cases_region %>% filter(between(month, as.Date('2021-06-01'), as.Date('2021-12-01'))) %>% group_by(areaName) %>% summarise(cases=sum(newCasesBySpecimenDate))
	}else if(period==4) { #Dec 2021-Apr 2022 / Omicrons 
		cases_region <- cases_region %>% filter(between(month, as.Date('2021-12-01'), as.Date('2022-04-01'))) %>% group_by(areaName) %>% summarise(cases=sum(newCasesBySpecimenDate))
	}else { # all period
		cases_region <- cases_region %>% filter(between(month, as.Date('2020-06-01'), as.Date('2022-04-01'))) %>% group_by(areaName) %>% summarise(cases=sum(newCasesBySpecimenDate))
	}
	
	# TODO solve dilema "POOLE" & "BOURNEMOUTH" & "CORNWALL" + "ISLES OF SCILLY" (duplicate rows)
	
	# Duplicate column
	cases_region <- cases_region %>% filter(areaName == "BOURNEMOUTH, CHRISTCHURCH AND POOLE") %>% bind_rows(cases_region, .) 
	# Change each to a different adm2 region (e.g. BOURNEMOUTH, CHRISTCHURCH AND POOLE -> 1. BOURNEMOUTH -> 2. POOLE)
	#cases_region %>% group_by(areaName) %>% filter(row_number()==1) %>% mutate(areaName=replace(areaName, "BOURNEMOUTH, CHRISTCHURCH AND POOLE", "BOURNEMOUTH"))
	cases_region <- cases_region %>% group_by(areaName) %>% mutate(areaName = ifelse(row_number() == 2, "POOLE", as.character(areaName)))
	cases_region[ cases_region == "BOURNEMOUTH, CHRISTCHURCH AND POOLE" ] <- "BOURNEMOUTH"
	cases_region <- cases_region %>% filter(areaName == "CORNWALL AND ISLES OF SCILLY") %>% bind_rows(cases_region, .)
	cases_region <- cases_region %>% group_by(areaName) %>% mutate(areaName = ifelse(row_number() == 2, "ISLES OF SCILLY", as.character(areaName)))
	cases_region[ cases_region == "CORNWALL AND ISLES OF SCILLY" ] <- "CORNWALL"
	
	print("head case counts")
	print(head(cases_region))
	
	return(cases_region)
}

cases_ltla <- cases_utla <- matches_cases_sequences_ltla <- matches_cases_sequences_utla <- matches_cases_sequences_all <- list()
for(i in 1:5) {
	cases_ltla[[i]] <- calculate_cases_england_regions("data/ltla_cases_2023-04-27.csv", period=i)
	cases_utla[[i]] <- calculate_cases_england_regions("data/utla_cases_2023-04-27.csv", period=i)
	
	matches_cases_sequences_ltla[[i]] <- cases_ltla[[i]] %>% inner_join(sc2_md_curated3_all, by=c("areaName"="region"), multiple="all")
	matches_cases_sequences_utla[[i]] <- cases_utla[[i]] %>% inner_join(sc2_md_curated3_all, by=c("areaName"="region"), multiple="all")
	
	matches_cases_sequences_all[[i]] <- rbind(matches_cases_sequences_ltla[[i]], matches_cases_sequences_utla[[i]])
	matches_cases_sequences_all[[i]] <- matches_cases_sequences_all[[i]][!duplicated(matches_cases_sequences_all[[i]]$sequence_name), ]
	
	print("diff between COG and LTLA/UTLA case counts regions")
	print(setdiff(unique(sc2_md_curated3_all$region), unique(matches_cases_sequences_all[[i]]$areaName)))
}

# Adjust breaks
library(scales)
trim_tails <- function(range = c(-Inf, Inf)) trans_new("trim_tails", 
	transform = function(x) {
		force(range)
		desired_breaks <- extended_breaks(n = 5)(x[x >= range[1] & x <= range[2]])
		break_increment <- diff(desired_breaks)[1]
		x[x < range[1]] <- range[1] - break_increment
		x[x > range[2]] <- range[2] + break_increment
		x
	},
	inverse = function(x) x,
	
	breaks = function(x) {
		force(range)
		extended_breaks(n = 7)(x)
	},
	format = function(x) {
		force(range)
		x[1] <- paste("<", range[1])
		x[length(x)] <- paste(">", range[2])
		x
})

plot_map_periods <- function(period=5, label="Jun 2020 to Apr 2022", matches_cases_sequences_df) {
	if(period==1) { #June-Dec 2020 / Other+EU1
		counts_seqs_all_map <- matches_cases_sequences_df %>% filter(between(month, as.Date('2020-06-01'), as.Date('2020-12-01'))) %>% group_by(areaName,cases) %>% summarise(n_seqs=n()) #%>% mutate(prop=n_seqs*100/cases)
	}else if(period==2) { #Jan-May 2021 / Alpha
		counts_seqs_all_map <- matches_cases_sequences_df %>% filter(between(month, as.Date('2021-01-01'), as.Date('2021-05-01'))) %>% group_by(areaName,cases) %>% summarise(n_seqs=n()) #%>% mutate(prop=n_seqs*100/cases)
	}else if(period==3) { #June-Dec 2021  / Deltas
		counts_seqs_all_map <- matches_cases_sequences_df %>% filter(between(month, as.Date('2021-06-01'), as.Date('2021-12-01'))) %>% group_by(areaName,cases) %>% summarise(n_seqs=n()) #%>% mutate(prop=n_seqs*100/cases)
	}else if(period==4) { #Dec 2021-Apr 2022 / Omicrons 
		counts_seqs_all_map <- matches_cases_sequences_df %>% filter(between(month, as.Date('2021-12-01'), as.Date('2022-04-01'))) %>% group_by(areaName,cases) %>% summarise(n_seqs=n()) #%>% mutate(prop=n_seqs*100/cases)
	}else { # all period
		counts_seqs_all_map <- matches_cases_sequences_df %>% group_by(areaName,cases) %>% summarise(n_seqs=n()) #%>% mutate(prop=n_seqs*100/cases)
	}
	#print(nrow(counts_seqs_all_map))
	
	counts_seqs_all_map$prop <- counts_seqs_all_map$n_seqs * 100 / counts_seqs_all_map$cases
	
	counts_seqs_all_map <- counts_seqs_all_map[order(counts_seqs_all_map$prop, decreasing = TRUE),]
	
	write.csv(counts_seqs_all_map, file=glue("stat_results/plots_paper/england_maps/period_{i}.csv"), quote=F, row.names=F)
	
	counts_seqs_shp <- counts_seqs_all_map %>% right_join(shp_england_df, by=c("areaName"="id"))
	counts_seqs_shp$n_seqs[ is.na(counts_seqs_shp$n_seqs) ] <- 0
	
	system("mkdir -p stat_results/plots_paper/england_maps/")
	#min_n <- 1 #min(counts_seqs_shp$prop)
	#max_n <- #max(counts_seqs_shp$prop)
	#n_breaks <- 7
	#break_every <- round((max(counts_seqs_shp$n_seqs)-min(counts_seqs_shp$n_seqs))/n_breaks, -3) # round to exact thousands
	print(glue("Period {period}"))
	print(mean(counts_seqs_shp$prop, na.rm=T))
	print(median(counts_seqs_shp$prop, na.rm=T))
	print(IQR(counts_seqs_shp$prop, na.rm=T))
 
	p <- ggplot(data = counts_seqs_shp, aes(x = long, y = lat, group = group, fill = prop)) + ggtitle(glue("{label} ({nrow(counts_seqs_all_map)} adm2 regions)")) +
		geom_polygon(color="grey70") + coord_equal() + theme_void() + theme(legend.position=c(0.1,0.45),plot.title = element_text(hjust = 0.5, size=9),legend.text=element_text(size=6), legend.title=element_text(size=7)) + #0.2,0.6
		scale_fill_continuous(type = "viridis", direction=-1, name="Proportion of cases sequenced (%)", trans = trim_tails(range = c(0,25)))#, breaks=seq(from=0, to=max_n, by=break_every), limits=c(0,max_n)) #na.value="white",
	ggsave(glue("stat_results/plots_paper/england_maps/period_{period}.png"), plot=p, width=5, height=6.5, dpi=600) #bg="white"
	return(p)
}

plot_maps_list <- list()
labels_england_plot <- c("Jun to Dec 2020", "Jan to May 2021", "Jun to Dec 2021", "Dec 2021 to Apr 2022", "Jun 2020 to Apr 2022")
for(i in 1:5) {
	plot_maps_list[[i]] <- plot_map_periods(period=i, labels_england_plot[i], matches_cases_sequences_all[[i]])
}

library(ggpubr)
plot_maps_list_1_4 <- plot_maps_list[c(1, 2, 3, 4)]
ggarrange(plotlist=plot_maps_list_1_4, ncol=2, nrow=2)
ggsave("stat_results/plots_paper/england_maps/maps_p1_p4.png", width=8, height=5.5, dpi=600, bg="white")

ggsave("stat_results/plots_paper/england_maps/map_complete_period.png", plot_maps_list[[5]], width=8, height=5.5, dpi=600, bg="white")

pl2 <- plot_maps_list[[5]]
ggarrange(pl1, pl2, ncol=1, nrow=2, labels=c("A","B"), align="h")
ggsave("stat_results/plots_paper/Fig2_seqDistr.pdf", width=6.5, height=8, dpi=600, bg="white")
