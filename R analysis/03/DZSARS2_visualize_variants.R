### visualize SARS2 within-host variants ###
## relative to DZ tiger reference sequence ##

##### prep & load data #####

#load packages
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(readxl)
library(ragg)
library(png)

# set working directory
setwd("")

# load alternate variant names for within-host variants
DZref_names <- read_xlsx("DZ_reference_variant_names.xlsx") %>% pull(variant_name)


# load data frame of cleaned variant data, combined replicates
# and filter for only within-host (DZ tiger reference) variants
df <- read_xlsx("viralrecon_variant_table_processed.xlsx") %>%
  filter(variant_name %in% DZref_names) 

######### outbreak.info-style variant plot ##############

# prep colors

colfunc <- colorRampPalette(c("#e8d04d", "#1b7939"))
colfunc(9)

MUTATIONPALETTE <- c("#E8D04D" ,"#CEC54A", "#B4BA48" ,"#9BAF45", "#81A443", "#679940", "#4E8E3E","#34833B", "#1B7939" )
borderColour = "white" #"#555555"


# plot:

df_plot <- df %>%
  mutate(collection_date = as.Date(collection_date),
         nice_date = format(collection_date, "%b %d"),
         id_date = paste0(animal_ID, " ", nice_date),
         NS = ifelse(effect == "missense_variant", "nonsynonymous\nSNV", NA))

vertical_wh <- ggplot(df_plot, aes(x= fct_reorder(variant_name, position),
                    y = fct_reorder(nice_date, collection_date, .desc=T))) +
  geom_tile(aes(color="#dedede"), color=borderColour) +
  geom_tile(aes(fill = mean_freq),
            color = "lightgrey") +
  geom_point(aes(shape = NS), size=1)  +
  scale_fill_gradientn(colours = MUTATIONPALETTE, 
                       limits = c(0.25,1),
                       labels = scales::percent) +
  scale_x_discrete(position="top") +
  scale_color_manual(values=NA) + 
  scale_shape_manual(values=16,na.translate=F)+
  guides(color=guide_legend("No data", 
                            override.aes=list(color="white"))) +
  labs(fill="Allele\nfrequency", color = "No data", shape="") +
  theme_bw() +
  facet_grid(animal_ID ~ .,  
             scales = "free_y", 
             space = "free_y" , 
             switch="both") + 
  theme(axis.line = element_blank(),
        axis.title= element_blank(),
        axis.ticks=element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size=8),
        legend.title.align = 0,
        legend.spacing.x = unit(-0.01, "cm"),
        legend.margin=margin(0,0,0,0),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(angle = 60,
                                   vjust = 0,
                                   hjust=0,
                                   size=10),
        strip.text.y.left = element_text(size=9.5, angle = 0),
        strip.placement = "outside", 
        plot.margin=grid::unit(c(0.1,0.6,0.1,0.1), "in"))
  
vertical_wh

#### overview of variants and predicted effects ##################

# first, look at variant observation frequencies and effects 
# what is the distribution of allele frequencies? 

hist(df$mean_freq) # very skewed towards 100%

# do any columns have NA values? 
names(df)[sapply(df, anyNA)]
# just gene_or_product which I introduced

# how many variants observed were different effects?
df %>%
  count(effect) 

# effect                  n
# 1 3'UTR                  2
# 2 5'UTR                  2
# 3 missense_variant      87
# 4 synonymous_variant     3

# TOTAL OBSERVATIONS: 94 relative to tigers

# note that no AY20 variants remain after using tigers as DZ_reference
df %>%
  filter(AY20 == T) %>%
  count(effect) 


# count unique variants and effects 

# not observations of variants, but unique varaints
# need to make wide variant table with columns for individual animals
# remove unnecessary columns for this purpose

df_wide <- df %>%
  select(-c(chrom, hgvs_c, hgvs_p, hgvs_p_1letter, lineage, 
            mean_depth, mean_ref_depth, mean_alt_depth, 
            animal_ID, species, collection_date)) %>%
  pivot_wider(id_cols = c(position, variant, ref, alt, gene, effect,
                          AY20, gene_or_product), 
              names_from = dataset_ID, 
              values_from=c(mean_freq),
              names_sort = T)

### TOTAL UNIQUE VARIANTS: 14
### TOTAL OBSERVATIONS: 94


# variant effects:
df_wide %>%
  count(effect)

### EFFECTS: 

# effect                      n
# <chr>                   <int>
# 1 3'UTR                  1
# 2 5'UTR                  2
# 3 missense_variant       9
# 4 synonymous_variant     2


## plot substitution types ####

# paste together reference and variant bases
df$substitution <- paste(df$ref, df$alt, sep=">")
# factor
as.factor(df$substitution)

# classify as transitions or transversions
df <- df %>%
  mutate(sub_type = ifelse((substitution %in% c("C>T", "T>C", "A>G", "G>A")), 
                           "transition", ifelse((substitution %in% c("G>T", "C>G", "G>C","C>A", "A>T","T>A", "T>G")), 
                                                "transversion", "NA")))

# filter out insertions/deletions and count single nucleotide substitutions
substitution_count <- df %>% 
  filter(str_length(substitution) <= 3) %>%
  group_by(substitution, sub_type) %>%
  summarise(count = length(substitution))


# do it by species
species_substitution_count <- df %>% 
  filter(sub_type != "other", 
         str_length(substitution) <= 3) %>%
  group_by(dataset_ID, sub_type) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  complete(dataset_ID, sub_type) %>%
  mutate(species = str_replace_all(dataset_ID, "[[_0-9]]", "")) %>%
  replace(is.na(.), 0) 

# had to be a bit roundabout to include samples where 0 transitions were observed
# and make NAs 0s

# plot wh variant substitution types
whsub <- df %>% ggplot(aes(x = sub_type, fill=factor(species, levels = c("Lion", "Hyena")))) +
  geom_bar(position = "dodge") + 
  labs(x = "Substitution relative to tiger reference", 
       fill = "Species",
       y = "Number of observations") +
  scale_fill_manual(values = c("#FEB853", "#6666FF"),
                    labels = c("Lion", "Hyena")) + 
  scale_x_discrete(labels=c("conservative_inframe_deletion" = "Conservative\ninframe deletion",
                            "disruptive_inframe_deletion" = "Disruptive\ninframe deletion",
                            "missense_variant" = "nonsynonymous\nSNV",
                            "synonymous_variant" = "synonymous\nSNV",
                            "upstream_gene_variant" = "Upstream\ngene variant",
                            "downstream_gene_variant" = "Downstream\ngene variant",
                            "5'UTR" = "5'UTR\nnoncoding",
                            "3'UTR" = "3'UTR\nnoncoding")) +
  theme_bw()  +
  theme(axis.text=element_text(size=10),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        strip.text = element_text(size=12),
        legend.position = "none",
        plot.margin = margin(t=0.5, l=0, r=0.5, b=0.5, unit = "cm"))
whsub


### look at sub-consensus within-host variation #####

df_sub <- df %>% filter(mean_freq <= 0.5) 
# five variants found at <50% frequency

5/94 # 5.3%

df_under_99 <- df %>% filter(mean_freq < 0.99) 
# 14 variants found at <=99% frequency

14/94 # 14.9%


#### count variants ####

count_df <- df %>% 
  group_by(animal_ID, dataset_ID, collection_date, species, dpi) %>%
  summarize(number_of_variants=n()) %>%
  mutate(collection_date = as.Date(collection_date, format = "%Y%m%d"))


# get mean number of variants in each species

count_df %>%
  group_by(species) %>%
  summarize(mean(number_of_variants), 
            median(number_of_variants))

# species `mean(number_of_variants)` `median(number_of_variants)`
# 1 Hyena                         4.17                            5
# 2 Lion                          1.33                            1


##### look at distribution of # of wh variants ##########

ggplot(count_df, aes(x=number_of_variants, fill=species)) +
  geom_histogram(binwidth=1, position="dodge", 
                 color="white",  size=0.3) + 
  scale_fill_manual(values = c("#6666FF", "#FEB853")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6)) +
  labs(x="Within-host variants per sample (relative to tiger reference)", 
       y="Number of samples") +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        strip.text = element_text(size=12),
        legend.position="none")


###### plot wh variant counts: species and #variants ######

# plot number of variants by species
wh_count <- ggplot(count_df, aes(x=factor(species, levels = c("Lion", "Hyena")), y=number_of_variants))+
  geom_violin(outlier.shape=NA, lwd=0.3) +
  geom_jitter(aes(fill=species), position = position_jitter(width=0.2, height= 0.05), 
              alpha=0.7, size=3, shape=21, stroke=0.2) +
  labs(x="", y="Within-host variants per sample\n(relative to tiger reference)") +
  scale_fill_manual(values = c("#6666FF", "#FEB853")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12),
        axis.title.x=element_blank(),
        legend.position = "none",
        legend.background=element_rect(color="black"),
        legend.title = element_text(size=12),
        plot.margin = margin(t=0.5, l=0, r=0.5, b=0.5, unit = "cm"))

wh_count


# which animals had no change, just 1 wh variant the whole time?
count_df %>%
  group_by(animal_ID) %>%
  summarize(n = n(), 
            sum = sum(number_of_variants))

#HyenaD, LionB, LionC, LionD, LionJ


########## plot wh variant effects ##########

wh_effect <- df %>% ggplot(aes(x = effect, fill=factor(species, levels = c("Lion", "Hyena")))) +
  geom_bar(position = "dodge") + 
  labs(x = "Substitution relative to tiger reference", 
       fill = "Species",
       y = "Number of observations") +
  scale_fill_manual(values = c("#FEB853", "#6666FF"),
                    labels = c("Lion", "Hyena")) + 
  scale_x_discrete(labels=c("conservative_inframe_deletion" = "Conservative\ninframe deletion",
                            "disruptive_inframe_deletion" = "Disruptive\ninframe deletion",
                            "missense_variant" = "Nonsynonymous\nSNV",
                            "synonymous_variant" = "Synonymous\nSNV",
                            "upstream_gene_variant" = "Upstream\ngene variant",
                            "downstream_gene_variant" = "Downstream\ngene variant",
                            "5'UTR" = "5'UTR\nnoncoding",
                            "3'UTR" = "3'UTR\nnoncoding")) +
  theme_bw()  +
  theme(axis.text=element_text(size=10),
        axis.title.y=element_text(size=12),
        axis.title.x=element_blank(),
        strip.text = element_text(size=12),
        legend.position = "none",
        plot.margin = margin(t=0.5, l=0, r=0.5, b=0.5, unit = "cm"),
        axis.text.x = element_text(size = 10, angle = -50, vjust=0.1, hjust=0.5))
wh_effect



############# make final composite plot with schematic and other plots ###################

# load in genome schematic
genome_schematic <- readPNG("genome_schematic.png",
                            native = T)

wh_schematic <- ggplot() + 
  background_image(genome_schematic) + 
  # This ensures that the image leaves some space at the edges
  theme(plot.margin = margin(t=0.2, l=0.2, r=0.2, b=0.2, unit = "cm"))
wh_schematic 

#or leave space for genome schematic
space <- ggplot(df) + geom_blank() + theme_void()


ggarrange(ggarrange(vertical_wh, wh_schematic, ncol=1, heights = c(8,1), 
          widths = c(1, 0.9), labels = c("A","")), 
          ggarrange(wh_count, whsub, wh_effect, ncol=1, labels = c("B", "C", "D")), 
          ncol=2, widths=c(1.5,1))
ggsave("fig3_withinhost_composite.pdf", width=8.5, height=11)



