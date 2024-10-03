### DZ SARS2 virus population selection analysis ###

# Population-level and gene-product-level summary data output by SNPGenie
# Script for cleaning and combining SNPGenie data is selection_cleanup.R

# 63 virus population samples collected between Oct 7-Dec 13
# 16 individuals of 3 species: tiger (N=2), lion (N=11), hyena (N=3) 

#### Load packages and data ####

# load packages
library(tidyverse)
library(openxlsx)
library(readxl)
library(ggpubr)
library(scales)
library(lme4)
library(car)
library(performance)

# set working directory 
setwd("")

# load in cleaned, combined datasets
# filtering out samples that had [less than 80% Coverage > 10x] (this stat comes from the nf-core/viralrecon run)

# id datasheet
ids <- read_xlsx("DZ_new_ids.xlsx")

# population-level summary data
df_pop <- read_csv("snpgenie_population_summary_cleaned.csv") 

# gene product-level summary data
df_prod <- read_csv("snpgenie_gene_product_cleaned.csv")


#### tidying data for visualization and analysis ###

## population level prep

# get the mean of the two technical replicates for each animal
# calculate piN/piS ratio (name the column piN_piS) and add logical column
# also add animal ID and dates columns

pop <- df_pop %>%
  mutate(piN = as.numeric(piN), 
         piS = as.numeric(piS)) %>%
  select(dataset_ID, replicate, piN, piS, pi) %>%
  group_by(dataset_ID) %>%
  summarize(mean_piN = mean(piN), 
            mean_piS = mean(piS),
            mean_pi = mean(pi)) %>%
  separate(dataset_ID, remove = F, 
           c("seq_id", "collection_date"), sep = "_") %>%
  mutate(piN_piS = (mean_piN/mean_piS),
         piN_piS_greater1 = ifelse((piN_piS>1), TRUE, FALSE),
         piN_piS_less1 = ifelse((piN_piS<1), TRUE, FALSE),
         collection_date = as.Date(collection_date, format="%Y%m%d"),
         species = case_when(grepl("Tiger", dataset_ID) ~ "tiger",
                             grepl("Lion", dataset_ID) ~ "lion",
                             grepl("Hyena", dataset_ID) ~ "hyena")) %>%
  left_join(ids, by = "seq_id") %>% select(animal_ID, everything()) 

# remove negative controls from dataset
pop <- pop %>% 
  filter(!dataset_ID %in% c("NCon", "Negative_control", "NTC_DZ"))

# pull out first test date for each animal
first <- pop %>%
  group_by(animal_ID) %>%
  summarize(first_test = first(collection_date))

# then go back to pop, add a column with first test date
# and use it to add another column with days since first test (dpi)
pop <- pop %>%
  left_join(first, by = "animal_ID") %>%
  mutate(dpi = (collection_date - first_test)) %>%
  select(-first_test)

rm(first) # clean up


## gene product-level prep

# get the mean of the two technical replicates for each animal
# calculate piN/piS (name the column piN_piS) and add logical columns as well
# also need to clean up piN/piS ratio column, as when piN and piS are both zero, 
# ratio is NaN, and I want it to be zero

prod <- df_prod %>%
  mutate(piN = as.numeric(piN), 
         piS = as.numeric(piS)) %>%
  select(dataset_ID, product, replicate, piN, piS) %>%
  group_by(dataset_ID, product) %>%
  summarize(mean_piN = mean(piN), 
            mean_piS = mean(piS)) %>%
  separate(dataset_ID, remove = F, 
           c("seq_id", "collection_date"), sep = "_") %>%
  mutate(piN_piS = (mean_piN/mean_piS),
         log10_piN_piS = (log10(piN_piS)),
         piN_piS_greater1 = ifelse((piN_piS>1), TRUE, FALSE),
         piN_piS_less1 = ifelse((piN_piS<1), TRUE, FALSE),
         collection_date = as.Date(collection_date, format="%Y%m%d"),
         species = case_when(grepl("Tiger", dataset_ID) ~ "tiger",
                             grepl("Lion", dataset_ID) ~ "lion", 
                             grepl("Hyena", dataset_ID) ~ "hyena"),
         product = factor(product, levels = c("ORF1ab", "S", "ORF3a","E","M",
                                              "ORF6","ORF7a","ORF8","N","ORF10"))) %>%
  left_join(ids, by = "seq_id") %>% select(animal_ID, everything()) 

# modify NaN to 0 from 0/0 ratio
prod$piN_piS[is.na(prod$piN_piS)] <- 0

# confirm that NCon samples have no values before removing from dataset
prod <- prod %>% 
  filter(!dataset_ID %in% c("NCon", "Negative_control", "NTC_DZ"))

# pull out first test date for each animal
first <- prod %>%
  group_by(animal_ID) %>%
  summarize(first_test = first(collection_date)) 

# then go back to pop, add a column with first test date
# and use it to add another column with days since first test (dpi)
prod <- prod %>%
  left_join(first, by = "animal_ID") %>%
  mutate(dpi = (collection_date - first_test)) %>%  
  select(-first_test)

rm(first) # clean up



#### Stats for nucleotide diversity over time ####

# use a linear mixed model to show statistical difference over time and among species 
# accounting for repeated measures of the same individual over time

# for stats, outcome variable does not need to be normally distributed as a univariate variable
# BUT LME models assume that the residuals of the model are normally distributed

# mean nucleotide diversity = mean_pi

# first try, check species interaction *
#lmm_mod_sp <- lmer(mean_pi ~ dpi*species + (1|animal_ID), pop)
#Anova(lmm_mod_sp, type = "3")

# dpi is significant, species interaction is not so drop * to +
lmm_pi <- lmer(mean_pi ~ dpi + species + (1|animal_ID), pop)
Anova(lmm_pi, type = "3")
summary(lmm_pi)
lmm_pi

# notes:
# if there was an interaction, the slopes would not be the same by species
# no interaction, so the slopes are the same, and there are different intercepts
# whatever is on the right side of the | operator is a factor and referred to as a “grouping factor”
# the random effects part tells you how much variance you find
# among levels of the grouping factor, plus the residual variance
# should capture all influence of animal_ID on mean_dpi and species
# the variance for animalID = 3.139e-11 not that important: they don't explain a lot of variation. 
# we can take the variance for animalID and divide it by the total variance: 3.139e-11/(3.234e-09+3.139e-11)
# so the differences between individual animals explain ~1% of the variance that’s “left over” after the variance explained by our fixed effects.
# the fixed effect part is similar to a regular linear model

# check assumptions:
plot(lmm_pi) # fitted vs. residuals  
plot(lmm_pi, type=c("p","smooth"), col.line=1)
plot(lmm_pi,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1) # scale-location plot for the assumption of equal variance

qqnorm(resid(lmm_pi)) # Q-Q plot, looks good



#### Stats for mean synonymous nucleotide diversity ####

# try for species interaction first
# lmm_piS_sp <- lmer(mean_piS ~ dpi*species + (1|animal_ID), pop)
# Anova(lmm_piS_sp, type = "3") #nothing sig

# species interaction is not sig so drop * to +
lmm_piS <- lmer(mean_piS ~ dpi + species + (1|animal_ID), pop)
Anova(lmm_piS, type = "3")
lmm_piS
summary(lmm_piS)

# check assumptions
plot(lmm_piS) # fitted vs. residuals  
plot(lmm_piS, type=c("p","smooth"), col.line=1)
qqnorm(resid(lmm_piS)) # residuals of model, looks ok

###### make table for both models ######

sjPlot::tab_model(lmm_pi, lmm_piS,
                  show.ci = F,
                  digits = 6,
                  digits.re = 12, #11 digits
                  p.style = "numeric_stars",
                  show.re.var = TRUE,
                  pred.labels = c("(Intercept)", "Time (days after first positive test)", "Species [lion]", "Species [tiger]"),
                  dv.labels = c("Nucleotide diversity", "Synonymous nucleotide diversity"))
                 # file = "lmm_stats_table.html")


####### Plotting #######

# not normal distribution, left skewed
hist(pop$mean_pi) 
hist(pop$mean_piS)

# visualize piS (relative effective pop size) by individual with collection dates

# first, factor species and IDs for a nice order for plotting
pop <- pop %>% mutate(species = factor(species, levels = c("tiger", "lion", "hyena")),
                      animal_ID = factor(animal_ID, levels = c("Tiger A", "Tiger B",
                                       "Lion A" , "Lion B" ,
                                       "Lion C","Lion D","Lion E","Lion F",
                                       "Lion G", "Lion H", "Lion I","Lion J",
                                       "Lion K", "Hyena A" ,"Hyena B" ,"Hyena D" )))
  
Ne_plot <- ggplot(pop, aes(x=collection_date, 
                           y=mean_piS, fill=species)) +
  geom_line(aes(color=species), show.legend=F, alpha=0.8) +  
  geom_point(alpha=0.8, 
             pch=21, stroke= 0.2, 
             size=3, show.legend=F) +
  facet_wrap(~animal_ID) +
  labs(x= "Collection date", 
       y = expression("Synonymous nucleotide diversity"~(pi~S))) +
  scale_fill_manual(values = c(
    "#FC6666" ,"#FEB853","#6666FF")) +
  scale_color_manual(values = c(
    "#FC6666" ,"#FEB853","#6666FF")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60,
                                   hjust=1,size=12),
        axis.title = element_text(size=14),
        strip.text = element_text(size=12))

Ne_plot


#### Selection: piN/piS by species ####

# not normal distribution, left skewed
hist(pop$piN_piS)

# how many individuals had piN>piS overall? 
pop %>%
  group_by(species) %>%
  count(piN_piS_greater1 == T) 
# 31 samples piN>piS, 32 piS>piN, interesting
# all 6 hyenas had piN>piS overall


##### plot for nucleotide diversity, linear mixed effects model #####

# for linear mixed models, we plot marginal effects 
# show the slopes assoc. with species and dpi effects

# plot separate slopes for fixed effects
lmm_div_plot <- ggplot(pop, aes(x=dpi, y=mean_pi, fill=species)) +
  geom_point(size = 3, alpha = 0.8, pch=21, stroke=0.2) +
  geom_line(aes(y=predict(lmm_pi), group=animal_ID, color=species), 
            alpha=0.8, show.legend=F) + # lines show model predictions
  labs(x = "Days since first positive test", 
       y = expression("Within-host\nnucleotide diversity"~(pi)),
       fill = "Host species") +
  scale_fill_manual(values = c(
    "#FC6666" ,"#FEB853","#6666FF"),
    labels = c("Tiger", "Lion", "Hyena")) +
  scale_color_manual(values = c(
    "#FC6666" ,"#FEB853","#6666FF"),
    labels = c("Tiger", "Lion", "Hyena")) +
    theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        plot.margin= margin(0.1,0.5,0.1,0.5, "in")) 

lmm_div_plot


#### make composite plot Figure 4: expansion & diversification ####

ggarrange(lmm_div_plot, Ne_plot, ncol=1, labels="AUTO", 
          font.label = list(size = 16),heights = c(1.2,2))
ggsave("fig4_virus_populations.pdf", width=8.5, height=11)


####### Plot piN/piS genome level ############

# I really just want to know what piN/piS looks like at the final timepoint collected

# make a list with the last test date for each animal
last_list <- pop %>%
  group_by(animal_ID) %>%
  reframe(last_test = last(collection_date),
            last_test = format(last_test, "%Y%m%d"),
            dataset_ID = paste0(seq_id, 
                                sep = "_", last_test)) %>%
  pull(dataset_ID)

# then go back to pop df, pull out just the last test date for each individual
last <- pop %>% filter(dataset_ID %in% last_list) %>%
  mutate(type = "Full genome, last test")

hist(last$piN_piS) # still skewed


# all data
ggplot(pop, aes(x=species, y=piN_piS)) +
  geom_boxplot() +
  geom_point(alpha=0.5, show.legend=F) +
  geom_hline(yintercept=1, color="#e43537", 
             linewidth=0.8, linetype=2) +
  labs(x= "Species", y = expression(pi*N/pi*S~" ratio"))+
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))

# last test date
last_piNpiS_plot <- ggplot(last, aes(x=species, y=piN_piS)) +
  geom_jitter(aes(fill=species), position = position_jitter(0.1), 
              pch=21, stroke=0.2, 
              alpha=0.8, size=3, show.legend=F) +
   geom_hline(yintercept=1, color="black", 
             linewidth=0.8, linetype=2) +
  stat_summary(fun.y = "mean", fun.min = "mean", fun.max= "mean", 
               size= 0.3, geom = "crossbar", aes(color = species), show.legend=F)+
  labs(x= "", y = expression(pi*N/pi*S~" ratio"))+ #(whole genome)"))+
  scale_color_manual(values = c( "#FC6666" ,"#FEB853","#6666FF"
                                )) +
  scale_fill_manual(values = c( "#FC6666" ,"#FEB853","#6666FF"
  )) +
  scale_x_discrete(labels = c("Tiger", "Lion", "Hyena")) +
  scale_y_continuous(limits = c(0,3))+
  facet_wrap(~type) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank(),
        strip.text = element_text(size=14))
        
last_piNpiS_plot


#### Gene level plots ####

# first, factor species and IDs for a nice order for plotting
prod <- prod %>% mutate(species = factor(species, levels = c("tiger", "lion", "hyena")),
                      animal_ID = factor(animal_ID, levels = c("Tiger A", "Tiger B",
                                                               "Lion A" , "Lion B" ,
                                                               "Lion C","Lion D","Lion E","Lion F",
                                                               "Lion G", "Lion H", "Lion I","Lion J",
                                                               "Lion K", "Hyena A" ,"Hyena B" ,"Hyena D" )))
## all genes piN/piS ##

ggplot(prod %>% filter(!is.infinite(piN_piS )), 
       aes(x = product, y = piN_piS)) +
  geom_jitter(alpha=0.5, size=2, show.legend=F) + #aes(shape=species), 
  geom_hline(yintercept=1, color="#e43537", 
             size=0.8, linetype=2) +
  labs(x= "Gene", y = expression(pi*N/pi*S~" ratio")) +
  facet_wrap(~species, scales = "free_y") +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))

# look at a few specific genes:
# ORF1ab: more samples were piS>piN (42), but many were piN>piS (21), so not an overall pattern
prod %>%
  filter(product == "ORF1ab") %>%
  filter(!piN_piS_less1) 

# N: most samples were piN>piS (57),
prod %>%
  filter(product == "N") %>%
  filter(!piN_piS_greater1) 
# but 6 were piS>piN, two tigers, two lions, two hyenas

# S: most samples were piN>piS (45),
S <-prod %>%
  filter(product == "S") %>%
  filter(!piN_piS_greater1) 
# but 17 were piS>piN, two tigers, two lions, two hyenas
# including Tiger B 10/13, Tiger A 10/7, and Hyena B 10/28
# so kind of a mixed bag


##### composite plot for Figure 5 ####

p1 <- ggplot(prod %>% filter(product =="N"), aes(x=species, y = piN_piS)) +
  geom_jitter(aes(fill=species), alpha=0.8, shape=21, stroke=0.2,
              size=3, show.legend=F, height=0, width=0.2) + 
  stat_summary(fun.y = "mean", fun.min = "mean", fun.max= "mean", 
               size= 0.3, geom = "crossbar", aes(color = species), show.legend=F)+
  geom_hline(yintercept=1, color="black", size=0.8, linetype=2) + 
  scale_color_manual(values = c("#FC6666" ,"#FEB853","#6666FF")) +
  scale_fill_manual(values = c("#FC6666" ,"#FEB853","#6666FF")) +
  scale_x_discrete(labels = c("Tiger", "Lion", "Hyena"))+
  scale_y_continuous(labels = label_number(accuracy=1))+
  facet_wrap(~product)+ 
  labs(x= "", y = expression(pi*N/pi*S~" ratio")) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank(),
        strip.text = element_text(size=14))

p2 <- ggplot(prod %>% filter(product %in% c("ORF1ab", "S", "ORF3a")), 
             aes(x=species, y = piN_piS)) +
  geom_jitter(aes(fill=species), alpha=0.8, shape=21, stroke=0.2,
              size=3, show.legend=F, height=0, width=0.2) + 
  stat_summary(fun.y = "mean", fun.min = "mean", fun.max= "mean", 
               size= 0.3, geom = "crossbar", aes(color = species), show.legend=F)+
  geom_hline(yintercept=1, color="black", size=0.8, linetype=2) + 
  scale_color_manual(values = c("#FC6666" ,"#FEB853","#6666FF")) +
  scale_fill_manual(values = c("#FC6666" ,"#FEB853","#6666FF")) +
  scale_x_discrete(labels = c("Tiger", "Lion", "Hyena"))+
  scale_y_continuous(labels = label_number(accuracy=1), limits = c(0,15))+
  facet_wrap(~product, ncol=3)+ 
  labs(x= "", y = expression(pi*N/pi*S~" ratio")) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank(),
        strip.text = element_text(size=14))

p3 <- ggplot(prod %>% filter(!product %in% c("ORF1ab", "S", "ORF3a","N")), aes(x=species, y = piN_piS)) +
  geom_jitter(aes(fill=species), alpha=0.8, shape=21, stroke=0.2,
              size=3, show.legend=F, height=0, width=0.2) + 
  stat_summary(fun.y = "mean", fun.min = "mean", fun.max= "mean", 
               size= 0.3, geom = "crossbar", aes(color = species), show.legend=F)+
  geom_hline(yintercept=1, color="black", size=0.8, linetype=2) + 
  scale_color_manual(values = c("#FC6666" ,"#FEB853","#6666FF")) +
  scale_fill_manual(values = c("#FC6666" ,"#FEB853","#6666FF")) +
  scale_x_discrete(labels = c("Tiger", "Lion", "Hyena"))+
  scale_y_continuous(labels = label_number(accuracy=0.1))+
  facet_wrap(~product, ncol=3)+ 
  labs(x= "", y = expression(pi*N/pi*S~" ratio")) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_blank(),
        strip.text = element_text(size=14))


ggarrange(ggarrange(last_piNpiS_plot, p1, ncol=2, labels=c("A", "B")), p2, p3, ncol=1, 
          heights = c(1, 1, 2), labels=c("", "C", "D"), font.label = list(size = 16))
ggsave("fig5_selection.pdf", width=8.5, height=11)


#### Export excel sheets for supplemental tables ####
wb <- createWorkbook("supplemental_selection_tables.xlsx")
addWorksheet(wb, "population_level")
addWorksheet(wb, "gene_product_level")
writeData(wb, "population_level", pop, borders="all")
writeData(wb, "gene_product_level", prod,borders="all")
saveWorkbook(wb, "supplemental_selection_tables.xlsx", overwrite = TRUE)

