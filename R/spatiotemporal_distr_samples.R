library(glue)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(viridis)

# Plot distribution of sequences included (fig. 2)
sc2_md_curated2 <- readRDS("rds/sc2_md_curated.rds")
sc2_md_curated2 <- include_major_lineage_column(sc2_md_curated2)
sc2_md_curated2$major_lineage <- ifelse(grepl("\\bOmicron_BA\\.\\*\\b", sc2_md_curated2$major_lineage), "Other", as.character(sc2_md_curated2$major_lineage))
sc2_tre_curated2 <- readRDS("rds/sc2_tre_curated.rds")


sc2_md_curated2$month <-as.Date(cut(sc2_md_curated2$sample_date, breaks = "1 month"))

df_lineages_per_month <- sc2_md_curated2 %>% group_by(month, major_lineage) %>% summarise(count_lin_month=n())
df_lineages_per_month <- df_lineages_per_month[!is.na(df_lineages_per_month$major_lineage),]

system("mkdir -p stat_results/plots_paper/england_maps/")

# Stacked bar chart of sequence distribution by sample date and lineage
df_lineages_per_month  %>% mutate(date=as.POSIXct(month)) %>% 
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
sort(unique(lookup_ltla_utla_adm2$LAD_code))

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

#print(length(unique(shp_england$id)) == length(unique(lookup_ltla_utla_england$LTLA22CD)))

sc2_md_curated2_test <- sc2_md_curated2
# extract matches between COG md region and adm2
sc2_md_curated2_test <- sc2_md_curated2_test %>% inner_join(lookup_ltla_utla_adm2, by=c("region"="adm2"), multiple="first")
unique(sc2_md_curated2_test$region)
#sc2_md_curated2_test$adm2 <- sc2_md_curated2_test$region
sc2_md_curated2_test <- sc2_md_curated2_test %>% dplyr::select(sequence_name, sample_date, lineage, major_lineage, region, mutations, month)
# extract NON-matches between COG md region and adm2
sc2_md_curated2_test_not <- sc2_md_curated2 %>% anti_join(lookup_ltla_utla_adm2, by=c("region"="adm2")) #multiple="first
#sc2_md_curated2_test_not$region <- ifelse(sc2_md_curated2_test_not$region %in% lookup_ltla_utla_adm2$LAD_name, yes=, no=)
#%>% inner_join(lookup_ltla_utla_adm2, by=c("region"="LAD_name"), multiple="first")
# Extract other matches with LTLAs
sc2_md_curated2_test_match_ltla <- sc2_md_curated2_test_not %>% inner_join(lookup_ltla_utla_adm2, by=c("region"="LAD_name"), multiple="first")
sc2_md_curated2_test_match_ltla$region <- sc2_md_curated2_test_match_ltla$adm2
unique(sc2_md_curated2_test_match_ltla$region)
sc2_md_curated2_test_match_ltla <- sc2_md_curated2_test_match_ltla %>% dplyr::select(sequence_name, sample_date, lineage, major_lineage, region, mutations, month)
# Extract other matches with UTLAs
sc2_md_curated2_test_match_utla <- sc2_md_curated2_test_not %>% inner_join(lookup_ltla_utla_adm2, by=c("region"="UTLA_name"), multiple="first")
sc2_md_curated2_test_match_utla$region <- sc2_md_curated2_test_match_utla$adm2
unique(sc2_md_curated2_test_match_utla$region)
sc2_md_curated2_test_match_utla <- sc2_md_curated2_test_match_utla %>% dplyr::select(sequence_name, sample_date, lineage, major_lineage, region, mutations, month)
sc2_md_curated2_all <- rbind(sc2_md_curated2_test, sc2_md_curated2_test_match_ltla, sc2_md_curated2_test_match_utla)
sc2_md_curated2_all <- sc2_md_curated2_all[!duplicated(sc2_md_curated2_all$sequence_name), ]
#length(unique(sc2_md_curated2_all$region)) #119 (LIVERPOOL included here but not in GADM, ok since only 1639 obs)
#length(unique(shp_england$NAME_2))
#nrow(sc2_md_curated2) - nrow(sc2_md_curated2_all) #44513 without regions after matches + 20885 blank = 65398

plot_map_periods <- function(period=5, label="Jun 2020 to Apr 2022") {
	if(period==1) { #June-Dec 2020 / Other+EU1
		counts_seqs_all_map <- sc2_md_curated2_all %>% filter(between(month, as.Date('2020-06-01'), as.Date('2020-12-01'))) %>% group_by(region) %>% summarise(n_seqs=n())
	}else if(period==2) { #Jan-May 2021 / Alpha
		counts_seqs_all_map <- sc2_md_curated2_all %>% filter(between(month, as.Date('2021-01-01'), as.Date('2021-05-01'))) %>% group_by(region) %>% summarise(n_seqs=n())
	}else if(period==3) { #June-Dec 2021  / Deltas
		counts_seqs_all_map <- sc2_md_curated2_all %>% filter(between(month, as.Date('2021-06-01'), as.Date('2021-12-01'))) %>% group_by(region) %>% summarise(n_seqs=n())
	}else if(period==4) { #Dec 2021-Apr 2022 / Omicrons 
		counts_seqs_all_map <- sc2_md_curated2_all %>% filter(between(month, as.Date('2021-12-01'), as.Date('2022-04-01'))) %>% group_by(region) %>% summarise(n_seqs=n())
	}else { # all period
		counts_seqs_all_map <- sc2_md_curated2_all %>% group_by(region) %>% summarise(n_seqs=n())
	}
	print(nrow(counts_seqs_all_map))
	#View(counts_seqs_all_map)
	
	counts_seqs_shp <- counts_seqs_all_map %>% right_join(shp_england_df, by=c("region"="id"))
	#View(counts_seqs_shp)
	counts_seqs_shp$n_seqs[ is.na(counts_seqs_shp$n_seqs) ] <- 0
	
	system("mkdir -p stat_results/plots_paper/england_maps/")
	min_n <- min(counts_seqs_shp$n_seqs)
	max_n <- max(counts_seqs_shp$n_seqs)
	n_breaks <- 7
	break_every <- round((max(counts_seqs_shp$n_seqs)-min(counts_seqs_shp$n_seqs))/n_breaks, -3) # round to exact thousands
	p <- ggplot(data = counts_seqs_shp, aes(x = long, y = lat, group = group, fill = n_seqs)) + ggtitle(glue("{label} ({nrow(counts_seqs_all_map)} adm2 regions)")) +
		geom_polygon(color="grey70") + coord_equal() + theme_void() + theme(legend.position=c(0.2,0.6),plot.title = element_text(hjust = 0.5, size=9),legend.text=element_text(size=6), legend.title=element_text(size=7)) +
		scale_fill_continuous(type = "viridis", direction=-1, name="# sequences", breaks=seq(from=0, to=max_n, by=break_every), limits=c(0,max_n)) #na.value="white",
	ggsave(glue("stat_results/plots_paper/england_maps/period_{period}.png"), plot=p, width=5, height=6.5, dpi=600) #bg="white"
	return(p)
}

p1 <- plot_map_periods(period=1, "Jun to Dec 2020")
p2 <- plot_map_periods(period=2, "Jan to May 2021")
p3 <- plot_map_periods(period=3, "Jun to Dec 2021")
p4 <- plot_map_periods(period=4, "Dec 2021 to Apr 2022")
p5 <- plot_map_periods(period=5, "Jun 2020 to Apr 2022")

library(ggpubr)
ggarrange(p1, p2, p3, p4, ncol=2, nrow=2)
ggsave("stat_results/plots_paper/england_maps/maps_p1_p4.png", width=8, height=5.5, dpi=600, bg="white")