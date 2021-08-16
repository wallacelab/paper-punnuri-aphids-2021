#! /usr/bin/Rscript

library(gridExtra)
library(ggpubr)
library(grid)
library(png)
options(stringsAsFactors=F)
setwd('~/Documents/Papers/PUNNURI_Sorghum Aphid Resistance/GWAS_github/publication_figures/')

# Arguments to pass along (easier to keep in a list like this)
args=list()
args$rmip_cutoff=0.05
args$perm_cutoff=0.01
args$winsize=5e6
args$step=5e5
args$offsets="sorghum_offsets.csv"
args$debug=FALSE  # Set to true to speed plotting during debugging of changes

# Output figure size (in inches)
out.width=8
out.dpi=600
legend.file="wallace-legend.png"

# Load helper functions
source('0_MakeManhattanPlot.r')
offsets = load_offsets(args$offsets)

# Phenotypes
pheno_name_key=c("Flowering Time" = "FLoweringtime",
                 "Greenhouse Damage" = "GHD",
                 "Plant Height (log)" = "Log_2019.2020_PH",
                 "D1 Aphid Count (log)" = "Log_D1_2019.2020AC",
                 "D1 Aphid Damage (log)" = "Log_D1_2019.2020",
                 "D2 Aphid Count (log)" = "Log_D2_2019.2020AC",
                 "D2 Aphid Damage (log)" = "Log_D2_2019.2020",
                 "NDRE" = "NDREAug25",
                 "NDVI" = "NDVIAug25",
                 "SAVI" = "SAVI_Aug25",
                 "NDRE (flowering covariate)" = "NDREAug25.with_flowering",
                 "NDVI (flowering covariate)" = "NDVIAug25.with_flowering",
                 "SAVI (flowering covariate)" = "SAVI_Aug25.with_flowering")



# Helper function to make plots for each phenotype
get_plots = function(traits){
    file.info = pheno_name_key[traits]
    glm_files = paste("Results_GLM/glm_results.",file.info,".csv.gz", sep='')
    farm_files = paste("Results_FarmCPU/1e_",file.info,".p05.rmip05.csv", sep='')
    
    # Make plots
    myplots=list()
    for(i in 1:length(traits)){
        myplots[[i]] = plot_manhattan(traits[i], glm_files[i], farm_files[i], args)
    }
    
    # Arrange into grid
    chrom_label=text_grob("Chromosome", size = 12, face = "bold") 
    left_label =text_grob(expression(bold(paste("-", log[10], " p-value (GLM)"))), size = 12, rot=90) 
    right_label =text_grob("RMIP (FarmCPU)", size = 12, face='bold', rot=-90) 
    legend=rasterGrob(readPNG(legend.file))
    grid.arrange(arrangeGrob(grobs=myplots, ncol=1), legend,
                 nrow=2, heights=c(10,1), bottom=chrom_label, left=left_label, right=right_label)
}


# Aphid damage
png("Figure - Aphid Damage Manhattan plots.png", width=out.width, height=8, units="in", res=out.dpi)
    get_plots(c("D1 Aphid Count (log)", "D1 Aphid Damage (log)", "D2 Aphid Damage (log)","Greenhouse Damage"))
dev.off()

# Drones
png("Figure - Drone Manhattan plots (no covariates).png", width=out.width, height=5.5, units="in", res=out.dpi)
    get_plots(c("NDRE", "NDVI", "SAVI"))
dev.off()
png("Figure - Drone Manhattan plots (flowering covariate).png", width=out.width, height=5.5, units="in", res=out.dpi)
    get_plots(c("NDRE (flowering covariate)", "NDVI (flowering covariate)", "SAVI (flowering covariate)"))
dev.off()

# Flowering + Height
png("Figure - Flowering and Height Manhattan plots.png", width=out.width, height=4, units="in", res=out.dpi)
    get_plots(c("Flowering Time", "Plant Height (log)"))
dev.off()

