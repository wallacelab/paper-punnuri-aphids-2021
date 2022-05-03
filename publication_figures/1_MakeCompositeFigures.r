#! /usr/bin/Rscript

library(gridExtra)
library(ggpubr)
library(grid)
library(png)
options(stringsAsFactors=F)
setwd('~/Documents/Papers/1_Submitted/PUNNURI_Sorghum Aphid Resistance/GWAS_github/publication_figures/')

# Arguments to pass along (easier to keep in a list like this)
args=list()
args$rmip_cutoff=0.05
args$perm_cutoff=0.01
args$winsize=5e6
args$step=5e5
args$offsets="sorghum_offsets.csv"
args$debug=FALSE

# Output figure size (in inches)
out.width=8
#out.height=3 # Depends on the specific plot
out.dpi=600
legend.file="wallace-legend.png"

# Load helper functions
source('0_MakeManhattanPlot.r')
offsets = load_offsets(args$offsets)

# Load key connecting human-fiendly phenotype names to ones used in pipeline
source("0_PhenotypeKey.r")


# Helper function to make plots for each phenotype
get_plots = function(traits, legend.height=0.04){
    file.info = pheno_name_key[traits]
    #glm_files = paste("Results_GLM/glm_results.",file.info,".csv.gz", sep='')
    glm_files = paste("IGNORE",file.info) # Ignoring GLM results now (reviewer comments)
    farm_files = paste("Results_FarmCPU/1e_",file.info,".p05.rmip05.csv", sep='')
    
    # Make plots
    myplots=list()
    for(i in 1:length(traits)){
        myplots[[i]] = plot_manhattan(traits[i], glm_files[i], farm_files[i], args)
    }
    
    # Arrange into grid
    chrom_label=text_grob("Chromosome", size = 12, face = "bold") 
    #left_label =text_grob(expression(bold(paste("-", log[10], " p-value (GLM)"))), size = 12, rot=90) 
    left_label =text_grob("RMIP (FarmCPU)", size = 12, face='bold', rot=90) 
    right_label =text_grob(expression(""), size = 12, rot=90) 
    legend=rasterGrob(readPNG(legend.file))
    grid.arrange(arrangeGrob(grobs=myplots, ncol=1), legend,
                 nrow=2, heights=c(1,legend.height), bottom=chrom_label, left=left_label, right=right_label)
}

# Aphid damage
png("Figure - Aphid Damage Manhattan plots.png", width=out.width, height=5.5, units="in", res=out.dpi)
    get_plots(c("AC1 (log)", "D1 (log)", "D2 (log)")) #,"Greenhouse Damage"))
dev.off()

# Drones
png("Figure - Drone Manhattan plots (no covariates).png", width=out.width, height=5.5, units="in", res=out.dpi)
    get_plots(c("NDRE", "NDVI", "SAVI"))
dev.off()
png("Figure - Drone Manhattan plots (flowering covariate).png", width=out.width, height=5.5, units="in", res=out.dpi)
    get_plots(c("NDRE (flowering covariate)", "NDVI (flowering covariate)", "SAVI (flowering covariate)"))
dev.off()

# Flowering
png("Figure - Flowering plot.png", width=out.width, height=2.5, units="in", res=out.dpi)
    get_plots(c("Flowering Time"), legend.height=0.1)
dev.off()




