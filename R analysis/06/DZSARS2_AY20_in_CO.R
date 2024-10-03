### Visualizing SARS2 AY20 lineage in CO

# load packages
library(tidyverse)
library(readxl)
library(openxlsx)
library(lubridate)

# set working directory 
setwd("")

#### look at case counts for AY20 and AY20+NA254V
# from 2021-06-01 to 2022-01-31 in CO
# metadata obtained from GISAID 2024-08-23

# load data for all AY20 seqs with NA254V mutation
# and pull IDs of AY20+NA254V seqs to filter them out of the overall AY20 data
NA254V_df <- read_tsv("01_input/ay20_NA254V.tsv")   %>%
  rename_with( ~ tolower(gsub(" ", "_", .x, fixed = TRUE))) %>%
  select(-submission_date) %>%
  mutate(lineage = "AY.20 + N A254V",
         collection_date = as_date(collection_date), #calculate epiyear and epiweek
         epiweek = epiweek(collection_date),
         epiyear = epiyear(collection_date))

NA254V_seqs <- NA254V_df %>% pull(accession_id)

# load all AY20 data from 2021-06-01 to 2022-01-31 in CO, filter out AY20+NA254V seqs then join with AY20+NA254V df
df <- read_tsv("01_input/all_ay20_Colorado_20210601_to_20220131_gisaid.tsv") %>%
  rename_with( ~ tolower(gsub(" ", "_", .x, fixed = TRUE))) %>%
  select(-submission_date) %>%
  filter(!accession_id %in% NA254V_seqs) %>% #filter out seqs with NA254V
  mutate(lineage = "AY.20",
         collection_date = mdy(collection_date), #calculate epiyear and epiweek
         epiweek = epiweek(collection_date),
         epiyear = epiyear(collection_date)) %>%
  full_join(NA254V_df)

# make dummy data for 0 cases in epiweek 2
epiweek2 <- data.frame(2, 2022, "AY.20", 0)
names(epiweek2) <- c("epiweek", "epiyear", "lineage", "n")

# group data by epiweek and add epiweek2
df_epiweek <- df %>%
  select(accession_id, collection_date, lineage, epiweek, epiyear) %>%
  group_by(epiweek, epiyear, lineage) %>%
  summarize(n=n()) %>%
  full_join(epiweek2) %>%
  arrange(epiyear)

# write epiweek data
write.xlsx(df_epiweek, "ay20_epiweek_data.xlsx")

# plot
ggplot(data = df_epiweek, aes(x=reorder(epiweek, epiyear), y=n, fill=lineage, group=lineage)) +
  geom_point(data = df_epiweek %>% filter(lineage == "AY.20"), 
             size=2.5, color = "darkgrey") +
  geom_line(data = df_epiweek %>% filter(lineage == "AY.20"), 
             linewidth=1, color = "darkgrey") +
  annotate("rect", xmin="40", xmax="43", ymin=-Inf, ymax=Inf, 
           alpha=0.2, fill="#6666FF") +
  annotate("text", x = "41", y = 43, label = "Oct 7th - Oct 30th", 
           color = "#6666FF", size=4) + #fontface="bold"
  geom_bar(data = df_epiweek %>% filter(lineage == "AY.20 + N A254V"),
           stat="identity", width=0.7) +
   scale_fill_manual(values =c("darkgrey", "#FEB853")) +
  #geom_vline(xintercept = c("40", "43"), linetype="dashed") +
  labs(x = "Epiweek (2021-2022)", 
       y = "Number of AY.20 sequences from humans in Colorado",
       fill = "SARS-CoV-2 lineage") +
  ylim(0,43)+
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

ggsave("AY20_NA254V.pdf", width=8.5, height=5)

