### analyze the variant output from the nfcore/viralrecon pipeline
## for the DZ SARS2 study 
# Laura Bashor

############ prep work ####################################

# load libraries
library(tidyverse)
library(openxlsx)
library(readxl)
library(data.table)

# set working directory
setwd()

# load data
# filter out sequences that were low coverage, would not be acceptable for downstream analyses

df <- read_csv("variants_long_table.csv") %>% 
  rename_all(., .funs = tolower) %>% # I don't like how iVar capitalizes everything
  rename("position"=pos, "dataset_ID"=sample, 
         "depth" = dp, "ref_depth" = ref_dp,
         "alt_depth" = alt_dp, "allele_frequency" = af) %>% # editing some column names
  filter(filter == "PASS") %>% # only variants that pass iVar fisher exact test
  select(-c(filter, caller)) %>% # don't need filter anymore, caller is always iVar
  filter(!dataset_ID %in% c("Lion10_20211115", "Lion10_20211115_R", "Lion11_20211104", "Lion11_20211104_R", 
                            "Lion1_20211102","Lion1_20211102_R", "Lion2_20211110", "Lion2_20211110_R", 
                            "Lion3_20211104","Lion3_20211104_R", # note Lion3_20211104 was bad, but R was good -- continuing w/o this sample for variant analysis
                            "Lion4_20211102", "Lion4_20211102_R","Lion5_20211102", "Lion5_20211102_R", 
                            "Lion5_20211104", "Lion5_20211104_R", "Lion6_20211029","Lion6_20211029_R", 
                            "Lion7_20211122", "Lion7_20211122_R")) %>% # note Lion7_20211122 was bad, but R was good -- continuing w/o this sample for variant analysis
  mutate(nucleotide_variant = paste0(ref, position, alt), #"nt "
         
# also, it calls the 5' and 3' UTR part of ORF1a and Spike respectively, weirdly
# change gene and effect for these regions
         gene = ifelse(position < 266, "5'UTR", ifelse(position > 29674, "3'UTR", 
                                                       gene)), 
         effect = ifelse(position < 266, "5'UTR", ifelse(position > 29674, "3'UTR", 
                                                       effect)))


# for some reason only orf1ab wasn't capitalized, fix that:
df$gene[df$gene == "orf1a"] <- "ORF1a"
df$gene[df$gene == "orf1b"] <- "ORF1b"


######### processing ##########

# make a new column labeling replicates (anything with _R is replicate 2, without it gets a 1)
df$rep = str_extract(df$dataset_ID, "_R$") #$ indicates end of string in regex
df$rep[is.na(df$rep)] <- 1
df$rep[df$rep == "_R"] <- 2
  
# now delete _R from sample names because that data is encoded in replicate column
df$dataset_ID <- gsub("_R$", "", df$dataset_ID)


# now average technical replicates for each sample
# what we want is to throw out any data where we didn't detect the variant in both replicates, whatever the reason for this may be

df2 <- df %>%
  mutate(rep = as.numeric(rep),
         variant = str_remove_all(hgvs_p_1letter, "p.")) %>% # get AA variant from hgvs.p annotation
  group_by(dataset_ID,  chrom, position, variant, ref, alt, nucleotide_variant, 
           gene, effect, hgvs_c, hgvs_p, hgvs_p_1letter, lineage) %>%
  summarise(rep_sum = sum(rep), # create a column summing the placeholder values for two technical replicates
            mean_freq = mean(allele_frequency),
            mean_depth = mean(depth),
            mean_ref_depth = mean(ref_depth), 
            mean_alt_depth = mean(alt_depth)) %>% 
  filter(rep_sum == 3) %>% # now use rep_sum to filter out any variants that weren't in both replicates (rep 1 and rep 2 should equal 3)
  ungroup() %>%
  select(-rep_sum) 

# sometimes we don't have an AA variant because it's an upstream or downstream gene variant
# replace the "." with the nucleotide variant 
df2$variant <- ifelse(df2$variant=="." , df2$nucleotide_variant, df2$variant)

##### Total variant observations: 2769


############ add additional gene/nsp annotations ######################

## load annotations and AY20 mutation list

# load in outbreak.info AY20 characteristic mutations
AY20 <- read_xlsx("AY20_variants_modified.xlsx") %>%
  select(-c(position, gene_or_product, notes)) %>%
  mutate(AY20 = T)

# add the AY20 and detailed gene product annotations to the big df

df3 <- df2 %>% 
  left_join(AY20, by = c("variant", "gene")) %>%
  replace_na(list(AY20 = FALSE)) %>%
  mutate(gene_or_product = case_when(
    between(position, 266, 805) ~ "leader_protein",
    between(position, 806, 2719) ~ "nsp2",
    between(position, 2720, 8554) ~ "nsp3",
    between(position, 8555, 10054) ~ "nsp4",
    between(position, 10055, 10972) ~ "3C-like proteinase",
    between(position, 10973, 11842) ~ "nsp6",
    between(position, 11843,	12091) ~ "nsp7",
    between(position, 12092,	12685) ~ "nsp8",
    between(position, 12686,	13024) ~ "nsp9",
    between(position, 13025,	13441) ~ "nsp10",
    between(position, 13442,	13468) ~ "RNA-dependent RNA polymerase",
    between(position, 13442,	13480) ~ "nsp11",
    between(position, 13468,	16236) ~ "RNA-dependent RNA polymerase",
    between(position, 16237,	18039) ~ "helicase",
    between(position, 18040,	19620) ~ "3'-to-5' exonuclease",
    between(position, 19621,	20658) ~ "endoRNAse",
    between(position, 20659,	21552) ~ "2'-O-ribose methyltransferase",
    TRUE ~ NA_character_
  ))

### get the metadata in order

# load in additional metadata
meta <- read_xlsx("DZ_data_cleaned.xlsx") %>%
  mutate(dataset_ID = dataset_id, animal_ID = animal_id) %>%
  select(animal_ID, sex, age, dataset_ID) %>%
  unique() 

# get info from dataset IDs and metadata
df4 <- df3 %>% 
  mutate(species = str_replace_all(dataset_ID, "[[_0-9]]", ""),
         collection_date = str_replace_all(dataset_ID, "(.*)_", ""),
         collection_date = ymd(collection_date), 
         seq_ID = str_replace_all(dataset_ID, "_(.*)", "")) %>%
  left_join(meta, by = "dataset_ID") %>% select(animal_ID, everything())
  

### calculate dpi: days post-first positive test

# pull out first test date for each animal
first <- df4 %>%
  group_by(seq_ID) %>%
  summarize(first_test = first(collection_date))

# then add a column with first test date
# and use it to add another column with days since first test (dpi)
# and make a nice variant name specifying the gene or UTR, or if it's upstream/downstream, just the nt change
df_final <- df4 %>%
  left_join(first, by = "seq_ID") %>%
  mutate(dpi = (collection_date - first_test),
         variant_name = ifelse(effect == "upstream_gene_variant", paste0("nt", variant),
                                ifelse(effect == "downstream_gene_variant", paste0("nt ", variant),
                                       paste0(gene, " ", variant)))) %>%
  select(-first_test)

rm(df, df2, AY20, df3, df4, first, meta) # clean up


# pull unique variants
unique_variants <- df_final %>% 
  select(position, variant_name, variant,ref, alt, gene, effect, hgvs_c, AY20) %>%
  unique()

# total unique variants: 63

########### write data ##############

# write unique variants for 
write_csv(unique_variants, "viralrecon_unique_variants_relative_to_MN908947.3.csv")

# write the variant table
wb <- createWorkbook("viralrecon_variant_table_processed.xlsx")
addWorksheet(wb, "variants")
writeData(wb, "variants", df_final, borders="all")
saveWorkbook(wb, "viralrecon_variant_table_processed.xlsx", overwrite = TRUE)

