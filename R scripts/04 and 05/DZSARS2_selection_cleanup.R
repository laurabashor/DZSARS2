### cleanup of selection analysis dataset from SNPGenie ###

# load libraries
library(tidyverse)
library(openxlsx)

# load in selection analysis tables from SnpGenie analysis
setwd("")

dir<-""
subdir<-list.dirs(dir, recursive=F)

# get dir/file names
ff <- do.call(rbind, lapply(subdir, function(x) {
  ff <-
    list.files(x,
               "\\.txt$",
               include.dirs = FALSE,
               full.names = TRUE)
  data.frame(
    dir = basename(x),
    file = basename(ff),
    fullpath = ff,
    stringsAsFactors = F
  )
}))

# helper function from Mr Flick https://gist.github.com/MrFlick/11407464
read.stack <-
  function(files,
           ...,
           select = NULL,
           extra = NULL,
           reader = read.table) {
    dd <- data.frame()
    if (!is.null(extra)) {
      stopifnot(is.list(extra))
      stopifnot(all(sapply(extra, length) == length(files)))
    }
    for (i in 1:length(files)) {
      d <- reader(files[i], ...)
      if (!is.null(select)) {
        stopifnot(all(select %in% names(d)))
        d <- d[, select]
      }
      if (!is.null(extra)) {
        d <- do.call(cbind, c(list(d), lapply(extra, '[', i)))
      }
      if (nrow(dd) > 0) {
        dd <- rbind(dd, d)
      } else {
        dd <- d
      }
    }
    dd
  }

#### population summary ####
# read pop summary files in as a dataframe and clean up
ff_pop <- ff %>%
  filter(file == "population_summary.txt")

df_pop_working <- read.stack(ff_pop$fullpath, extra=list(file=ff_pop$file, dir=ff_pop$dir))

# helper function from https://stackoverflow.com/questions/32054368/use-first-row-data-as-column-names-in-r
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

df_pop_working <- header.true(df_pop_working)

names(df_pop_working)

df_pop_working$file <- sapply(strsplit(df_pop_working$file, "/\\s*"), tail, 1)
df_pop_working$file <- gsub("\\..*","", df_pop_working$file)

#helpful deletion of all the headers from each table which were introduced as rows
#function from https://stackoverflow.com/questions/39106128/delete-every-evenuneven-row-from-a-dataset
toDelete <- seq(1, nrow(df_pop_working), 2)
df_pop_working <- df_pop_working[ toDelete ,]

# final df_pop with replicates
df_pop <- df_pop_working %>%
  mutate(dataset_ID = file) %>%
  select(dataset_ID, everything()) %>%
  select(-c(population_summary.txt, 
            file,
            Hyena1_20211028_R.MN908947.3_SNPGenie_Results)) %>%
  filter(!dataset_ID %in% c("Lion1_20211102", "Lion1_20211102_R", # remove low-coverage samples that failed NGS QC
                           "Lion10_20211115", "Lion10_20211115_R",
                           "Lion11_20211104", "Lion11_20211104_R",
                           "Lion2_20211110", "Lion2_20211110_R",
                           "Lion3_20211104",#"Lion3_20211104_R", R was not good
                          "Lion4_20211102", "Lion4_20211102_R",
                           "Lion5_20211102", "Lion5_20211102_R",
                           "Lion5_20211104", "Lion5_20211104_R",
                           "Lion6_20211029","Lion6_20211029_R",
                           "Lion7_20211122")) #"Lion7_20211122_R")) #R was not good

# make a column indicating replicates
df_pop$replicate = str_extract(df_pop$dataset_ID, "_R")
df_pop$replicate[is.na(df_pop$replicate)] <- 1
df_pop$replicate[df_pop$replicate == "_R"] <- 2

# now create a column "ID_rep" that has the full ID with replicate
# and remove the replicate indicator from the dataset_ID column

df_pop <- df_pop %>% 
  mutate(ID_rep = dataset_ID)

df_pop$dataset_ID <- gsub("_R", "", df_pop$dataset_ID)


#### product results ####

# read product summary files in as a dataframe and clean up

ff_prod <- ff %>%
  filter(file == "product_results.txt")

df_prod_working <- read.stack(ff_prod$fullpath, extra=list(file=ff_prod$file, dir=ff_prod$dir))

df_prod_working <- header.true(df_prod_working)

names(df_prod_working)

df_prod_working$file <- sapply(strsplit(df_prod_working$file, "/\\s*"), tail, 1)
df_prod_working$file <- gsub("\\..*","", df_prod_working$file)

# final version of df_prod
df_prod <- df_prod_working %>%
  mutate(dataset_ID = file) %>%
  filter(dataset_ID != "file") %>%
  select(dataset_ID, everything()) %>%
  select(-c(product_results.txt, 
            Hyena1_20211028_R.MN908947.3_SNPGenie_Results)) %>%
  filter(!dataset_ID %in% c("Lion1_20211102", "Lion1_20211102_R", # remove low-coverage samples that failed NGS QC
                            "Lion10_20211115", "Lion10_20211115_R",
                            "Lion11_20211104", "Lion11_20211104_R",
                            "Lion2_20211110", "Lion2_20211110_R",
                            "Lion3_20211104",#"Lion3_20211104_R",
                            "Lion4_20211102", "Lion4_20211102_R",
                            "Lion5_20211102", "Lion5_20211102_R",
                            "Lion5_20211104", "Lion5_20211104_R",
                            "Lion6_20211029","Lion6_20211029_R",
                            "Lion7_20211122")) #"Lion7_20211122_R"))

# make a column indicating replicates
df_prod$replicate = str_extract(df_prod$dataset_ID, "_R")
df_prod$replicate[is.na(df_prod$replicate)] <- 1
df_prod$replicate[df_prod$replicate == "_R"] <- 2

# now create a column "ID_rep" that has the full ID with replicate
# and remove the replicate indicator from the dataset_ID column

df_prod <- df_prod %>% 
  mutate(ID_rep = dataset_ID) 

df_prod$dataset_ID <- gsub("_R", "", df_prod$dataset_ID)

# for some reason only orf1ab wasn't capitalized, fix that:
df_prod$product[df_prod$product == "orf1ab"] <- "ORF1ab"


#### save cleaned files ####

# save cleaned files
write_csv(df_pop, "snpgenie_population_summary_cleaned.csv")
write_csv(df_prod, "snpgenie_gene_product_cleaned.csv")

