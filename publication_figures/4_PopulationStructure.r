#! /usr/bin/Rscript

library(dplyr)
library(ggplot2)
library(ggrepel)
library(zoo)
# library(ggcorrplot)

# Plot a quick population structure plot showing resistance lines vs checks
setwd("~/Documents/Papers/1_Submitted/PUNNURI_Sorghum Aphid Resistance/GWAS_github/publication_figures/")
genos=readRDS("../farmcpu_pipeline/0_Inputs/281_Genotyp_Numeric.rds")
num_sites=10000 # Number of sites to subsample to before calculating distances
checks=c("PI257599", "PI533794", "PI656001")
best=c("PI276837", "PI597964", "PI656047", "PI659753")
set.seed(1) # For randomization
rmes1 = list(chr=6, start=324855, end=326854) # RMES1 locus

# Reformat genotypes
taxa=sub(genos$taxa, pattern="\\..+", repl="") # Reformat taxa names
rownames(genos) = taxa
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
is_check = rownames(pcoa) %in% checks
is_best = rownames(pcoa) %in%  best

# Update PCoA with checks vs controls
pcoa = pcoa %>%
  mutate(type=ifelse(is_check, "check", ifelse(is_best, "best", "other")),
         name = ifelse(type == "other", "", taxa)) %>%
  mutate(name = sub(name, pattern="\\..+", repl="")) #%>%
  #mutate(type = factor(type, levels=c("best", "check", "other")))# Remove name bits after first period
names = pcoa %>%
  filter(name != "")

# Plot PCoA
myplot = ggplot(pcoa %>% arrange(desc(type))) + # Sort descending so best and checks on top
  aes(x=PC1, y=PC2, fill=type, label=name) +
  geom_point(size=3, shape=21, color="black", stroke=0.25) +
  geom_label_repel(data=names, size=1.5, 
                   fontface="bold", fill="white", label.padding=0.1)
ggsave(myplot, file="Figure - Population PCoA.png", width=4, height=3, dpi=300)


# Calculate genetic similarity; as function in case want to do subsets elsewhere in genome
make_corplot = function(targets, mygenos){
  similar = matrix(NA, nrow=length(targets), ncol=length(targets), dimnames = list(targets, targets))
  subgenos = mygenos[targets,]
  for(i in 1:length(targets)){
    for(j in 1:i){
      similar[i,j] = similar[j,i] = sum(subgenos[i,] == subgenos[j,])/ncol(subgenos)
    }
  }
  ## Hierarchically cluster
  clust.order = hclust(as.dist(1-similar))$order
  similar=similar[clust.order, clust.order]
  
  cordata = reshape2::melt(similar) %>%
   mutate(Var1 = factor(Var1, levels=rownames(similar)), Var2=factor(Var2, levels=rownames(similar))) %>% # Put in right order
    filter(as.numeric(Var1) <= as.numeric(Var2))
  corplot = ggplot(cordata) +
    aes(x=Var1, y=Var2, fill=value, label=round(value, digits=2)) +
    geom_raster() +
    scale_y_discrete(limits=rev) +
    scale_fill_gradient(low="white", high="dodgerblue", na.value="red", limits=c(0.5,1)) +
    geom_text(fontface="bold", size=2.5) +
    theme(axis.title = element_blank(), legend.position="none", 
          axis.text=element_text(face="bold"), axis.text.x = element_text(angle=45, hjust=1),
          axis.ticks = element_blank(),
          axis.line = element_blank(), panel.grid=element_blank(),
          panel.background = element_blank())
}

targets = c(checks, best)
corplot = make_corplot(targets, genos) +
  ggtitle("Genetic Similarity - All")

# Plot genetic similarities
ggsave(corplot, file="Figure - Genetic Similarity.png", width=3, height=3, dpi=300)



# Plot genetic similarity over chromosome distance
span=50000
genos.chrom = as.integer(sub(colnames(genos), pattern="^S(.+)_.+", repl="\\1"))
genos.pos = as.integer(sub(colnames(genos), pattern="^S.+_(.+)", repl="\\1"))
focus = genos[targets,genos.chrom==rmes1$chr & genos.pos >= rmes1$start - span &
                genos.pos <= rmes1$end + span]

focus.pos = as.integer(sub(colnames(focus), pattern="^S.+_(.+)", repl="\\1"))
donor="PI257599"
# Similarity to PI257599, which donated RMES1
similar = matrix(NA, nrow=nrow(focus), ncol=ncol(focus), dimnames=dimnames(focus))
for(i in 1:nrow(similar)){
  similar[i,] = focus[i,] == focus[donor,]
}
# Create sliding window of similarity
winsize=4
similar.window = lapply(targets, function(t){
  data.frame(line=t, 
             start=rollapply(focus.pos, width=winsize, FUN=min),
             stop=rollapply(focus.pos, width=winsize, FUN=max),
             similarity=rollapply(similar[t,], width=winsize, FUN=mean, na.rm=TRUE))
}) %>%
  bind_rows() %>%
  mutate(pos=(start + stop)/2, line=factor(line, levels=c(checks, best))) 

# Plot
simplot = ggplot(similar.window) +
  aes(x=pos, y=similarity, color=line) +
  geom_line() +
  #geom_vline(xintercept=(rmes1$start + rmes1$end) / 2, size=1, alpha=0.5) + 
  annotate("rect", xmin=rmes1$start, xmax=rmes1$end, ymin=0, ymax=Inf, 
            color=NA, fill="black", alpha=0.4) + 
  annotate("text", x=(rmes1$start + rmes1$end)/2, y=0, hjust=0, label="RMES1",
           angle=90, size=2) +  
  facet_wrap(~line, ncol=1) +
  theme(legend.position = "none") +
  labs(x="Position (SBI-06)", y=paste("Similarity to",donor))
ggsave(simplot, file="Figure - Similarity around RMES1.png", width=4, height=8, dpi=300)
    



