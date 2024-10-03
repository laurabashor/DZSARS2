### DZSARS2 NGS and RT-PCR results overview ###

### prep ###

# load packages

library(readxl)
library(openxlsx)
library(tidyverse)
library(lubridate)

# set working directory
setwd("")

# load list of best sequences chosen from the two technical replicates (includes replicate info)
best_names <- read.table("best_seqs_list.txt", col.names="names") %>%
  pull(names)

# load list of all samples with NGS data (no rep info)
NGS_names <- read.table("seqs_list.txt", col.names="names") %>%
  pull(names)

# load list of timeline dates for plotting
timeline_dates <- read.table("timeline_dates.txt", col.names="dates") %>%
  mutate(dates = as.Date(dates))

# load Pangolin table from nf-core/viralrecon MultiQC and clean up datasheet

df <- read_csv("summary_variants_metrics_mqc.csv") %>%
  rename_all(function(x) gsub(" ", "_", x)) %>% # get rid of spaces in column names
  rename_all(., .funs = tolower) %>% # make all lowercase too
  rename("name" = 1) %>% 
  separate_wider_delim(name, delim = "_", 
                       names = c("seq_id", "collection_date", "rep"),
                       too_few = "align_start",
                       cols_remove=F) %>%
  mutate(collection_date = ymd(collection_date),
         species = str_remove_all(seq_id, "[0-9]"),
         animal_number = str_remove_all(seq_id, "[a-z]"),
         ngs = ifelse((name %in% best_names), "seq_data", "no_seq_data"))


### RT-PCR data ###

# read in inconclusives

df_incon <- read_xlsx("DZ_incon_PCR_data.xlsx") %>%
  separate_wider_delim(Sample_ID, delim = "_", 
                       names = c("seq_id", "collection_date"),
                       too_few = "align_start",
                       cols_remove=F) %>%
  #separate(Sample_ID, c("seq_id", "collection_date")) %>%
  mutate(collection_date = ymd(collection_date),
         ngs = ifelse((Sample_ID %in% NGS_names), 
                      "seq_data", "no_seq_data"),
         dataset_ID = Sample_ID)

# read in positives and new animal ids and merge with inconclusives

# id datasheet
ids <- read_xlsx("DZ_new_ids.xlsx")

# combine everything
df_all_pcr <- read_xlsx("DZ_PCR_data_cleaned.xlsx") %>%
  mutate(ngs = ifelse((dataset_id %in% NGS_names), 
                      "seq_data", "no_seq_data")) %>%
  left_join(df_incon, by = c("collection_date", "seq_id", "ngs", c("dataset_id"="dataset_ID"))) %>%
  select(-c(Sample_ID, animal_id)) %>%
  left_join(ids, by = "seq_id") %>% select(animal_ID, everything()) %>%
  mutate(collection_date = ymd(collection_date),
         collection_date = as.character(collection_date),
         strain = paste0(animal_ID, "_", collection_date),
         strain = str_replace_all(strain, " ", ""),
         strain = str_replace_all(strain, "-", ""))

# write cleaned up datasheet

wb <- createWorkbook("all_DZ_metadata_cleaned.xlsx")
addWorksheet(wb, "all_metadata")
writeData(wb, "all_metadata", df_all_pcr, borders="all")
saveWorkbook(wb, "all_DZ_metadata_cleaned.xlsx", overwrite = TRUE)


### plotting timeline ###

timeline_df <- df_all_pcr %>%
  mutate(collection_date = as.Date(collection_date),
         animal_ID = factor(animal_ID, levels = c( "Tiger A","Tiger B", "Lion A" , "Lion B","Lion C" , 
                                                   "Lion D" , "Lion E" ,
                                                   "Lion F" , "Lion G","Lion H" ,"Lion I" ,"Lion J", "Lion K" ,
                                                   "Hyena A", "Hyena B","Hyena C", "Hyena D" )))

ggplot(timeline_df %>%
         drop_na(result) %>% 
         mutate(result = factor(result, 
                                levels= c("Positive", 
                                          "Inconclusive", 
                                          "Negative"))),
       aes(x=collection_date, 
           y=animal_ID)) + 
  geom_line(aes(group=animal_ID), alpha=0.5) +
  geom_point(data = timeline_df %>% filter(ngs == "seq_data"),
             aes(x=collection_date,
                 y=animal_ID, 
                 shape = ngs), fill=NA, size=5, stroke=1, color="#6666FF") +
  geom_point(aes(fill = result), #size=ngs
             shape=21, size=3) + 
  scale_fill_discrete(type = c("#FC6666", 
                               "#e8d04d", 
                               "lightgrey"),
                      na.value = NA) +
  # scale_size_manual(values = c(3, 4),
  #                       labels = c("No NGS data", "NGS data"))+
  scale_shape_manual(values = c("seq_data" = 21),
                     name = "",
                     labels = "Sequence data") +
  scale_x_date(
    date_labels = "%b %d",
    breaks = timeline_dates$dates) + #unique(timeline_df$collection_date))+
  labs(y="", x="", fill="qRT-PCR Result")+ 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60,
                                   hjust=1,size=12),
        axis.text.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14))

ggsave("fig1_timeline_with_NGS.pdf", width=9, height=6)
