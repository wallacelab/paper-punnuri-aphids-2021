#! /usr/bin/Rscript

library(dplyr)
library(ggplot2)
library(ggrepel)

# Plot a quick population structure plot showing resistance lines vs checks
setwd("~/Documents/Papers/1_Submitted/PUNNURI_Sorghum Aphid Resistance/GWAS_github/publication_figures/")
genos=readRDS("../farmcpu_pipeline/0_Inputs/281_Genotyp_Numeric.rds")
num_sites=10000 # Number of sites to subsample to before calculating distances
checks=c("PI257599", "PI533794", "PI656001")
best=c("PI276837", "PI597964", "PI656047", "PI659753")
set.seed(1) # For randomization

# Reformat genotypes
rownames(genos) = genos$taxa
genos$taxa=NULL
genos=as.matrix(genos)

# Make distance matrix of genotypes; randomly subsample to 10k loci to speed computation
targets = sample(1:ncol(genos), size=num_sites)
pcoa = genos[,targets] %>%
  dist() %>%
  cmdscale() %>%
  as.data.frame() %>%
  rename("PC1" = "V1", "PC2" = "V2")

# Assign whether each check or best hit
taxa=rownames(genos)
find_taxa = function(x){
  taxa[which(grepl(x, taxa))]
}
is_check = sapply(checks, find_taxa)
is_best = sapply(best, find_taxa)

# Update PCoA
pcoa = pcoa %>%
  mutate(type=ifelse(taxa %in% is_check, "check", ifelse(taxa %in% is_best, "best", "other")),
         name = ifelse(type == "other", "", taxa)) %>%
  mutate(name = sub(name, pattern="\\..+", repl="")) #%>%
  #mutate(type = factor(type, levels=c("best", "check", "other")))# Remove name bits after first period
names = pcoa %>%
  filter(name != "")


# Plot
myplot = ggplot(pcoa %>% arrange(desc(type))) + # Sort descending so best and checks on top
  aes(x=PC1, y=PC2, fill=type, label=name) +
  geom_point(size=3, shape=21, color="black", stroke=0.25) +
  geom_label_repel(data=names, size=1.5, 
                   fontface="bold", fill="white", label.padding=0.1)
ggsave(myplot, file="Figure - Population PCoA.png", width=4, height=3, dpi=300)

