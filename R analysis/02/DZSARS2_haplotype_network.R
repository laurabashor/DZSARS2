### DZSARS2 haplotype network ###

#load packages
library(pegas)
library(tidyverse)
library(RColorBrewer)

# set working directory
setwd("")

# read in data
data <- read.FASTA("DZ_latest_best_seqs_aligned.afa")

meta <- data.frame(names(data)) %>%
  rename("name" = 1) %>%
  mutate(
         species = str_replace_all(name, "[[_0-9]]", ""),
         date = str_replace_all(name, "(.*)_", ""),
         date = ymd(date),
         id = str_replace_all(name, "_(.*)", ""))
  
# how many haplotypes & who shares a haplotype
hap <- haplotype(data)
hap

# list which individuals have which haplotype
hapInfo <- stack(setNames(attr(hap,"index"),rownames(hap)))
head(hapInfo)

# index and merge with metadata
names(hapInfo) <- c("index","haplotype")

merged <- data.frame(cbind(hapInfo,meta[hapInfo$index,]))
head(merged)

# make a haplotype network
net <- haploNet(hap)
plot(net, size=attr(net,"freq"), labels=F) #*.1

# prep colors for species
colors <- c("#6666FF","#FEB853","#FC6666")

# plot colored by species
pie <- table(merged$haplotype,merged$species)
head(pie)

# set margins and character expansion (for legend text size)
#par(mar = c(bottom, left, top, right))
op <- par(mar = c(0, 0, 0, 0), cex = 3.5)

# plot and save to pdf
pdf("best_latest_hapnet.pdf", width=6, height=5)
plot(net, size=attr(net,"freq"), labels=F, #*.2
     pie=pie, bg = colors)
legend("bottomleft", colnames(pie), col=colors,
       pch=19, ncol=2, inset=0)
dev.off()

