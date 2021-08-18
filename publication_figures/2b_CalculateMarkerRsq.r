#! /usr/bin/Rscript

# Calculate the r-squared of markers for each trait (basically, 
#  how much phenotypic variance they explain, both individually and as a group)

# TODO: Confirm heritability correct with Addissu. Seems kind of high

library(rrBLUP)
library(GAPIT3)
options(stringsAsFactors=F)
setwd('~/Documents/Papers/PUNNURI_Sorghum Aphid Resistance/GWAS_github/publication_figures/')

# Parameters
glm_max_p=0.01 # Significant permutation p-value for GLM hits (max value)
farm_min_rmip=0.05 # Significant RMIP value for FarmCPU resampling analysis (min value)
num_pcs=3 # Number of principal coordinates to include in the model
source("0_PhenotypeKey.r") # Key connecting human-fiendly phenotype names to ones used in pipeline

# Input files
hits = read.csv("Table - Combined hits.full.csv")
genos = readRDS("../farmcpu_pipeline/0_Inputs/281_Genotyp_Numeric.rds")
phenos = read.csv("../farmcpu_pipeline/1_traits.phenos.csv")
covariates = read.csv("../farmcpu_pipeline/1_traits.covariates.csv")

# Add additional trait columns for drone data with flowering covariates
phenos$NDREAug25.with_flowering = phenos$NDREAug25
phenos$NDVIAug25.with_flowering = phenos$NDVIAug25
phenos$SAVI_Aug25.with_flowering = phenos$SAVI_Aug25

# Calculate kinship & principal components
snps = genos[,-1]
rownames(snps) = genos$taxa
kinship = GAPIT3::GAPIT.kinship.VanRaden(snps)
pca = prcomp(snps)$x[,1:num_pcs]
rownames(pca) = rownames(snps)

# Combine into a single data frame for rrBLUP, keeping the trait names
traits = names(phenos)[-1]
all.data = cbind(phenos, covariates, pca)

# Helper function to return the names of covariates depending on a trait name
get_covariates = function(trait){
    covariates="FLoweringtime" # Flowering time by default
    # Greenhouse damage has nSPAD reading as covariate
    if(trait == "GHD"){
        covariates="nSpad"
    }
    # nSPAD, flowering, and height have no covariates, as do basic drone data
    if(trait %in% c("nSpad", "FLoweringtime", "Log_2019.2020_PH", 
                    "NDVIAug25", "NDREAug25", "SAVI_Aug25")){
        covariates=NULL
    }
    return(covariates)
}

# Get residuals for each trait
models = lapply(traits, function(mytrait){
    # Get covariates
    mycov = get_covariates(mytrait)
    mycov = c(mycov, colnames(pca)) # Add PCs as covariates
    
    # Run model
    mymodel = kin.blup(all.data, geno="Taxa", pheno=mytrait, K=kinship)
    return(mymodel)
})
names(models) = traits

# Get narrow-sense heritability estimates
h2 = lapply(models, function(m){
    m$Vg / (m$Vg+m$Ve)
})

# Extract residuals
residuals = lapply(models, function(m){
    myresid = m$resid
    names(myresid) = m$pred
})

# Helper function to fit models of SNPs
get.r2 = function(myresiduals, mysnps){
    # If no SNPs, return missing
    if(length(mysnps) == 0){
        return(NA)
    }
    
    # Extract target SNPs
    target.snps = snps[,mysnps, drop=FALSE]
    if(!identical(names(myresiduals), rownames(target.snps))){
        stop("Residual and SNP names do not match")
    }
    
    # Merge into a single data frame
    mydata = data.frame(residuals=myresiduals, target.snps)
    
    # Fit basic linear regression model
    independents = paste(mysnps, collapse=" + ")
    myformula = paste("residuals", independents, sep=" ~ ")
    mymodel = lm(myformula, data=mydata)
    
    # Return R-squared value
    return(round(summary(mymodel)$r.squared, digits=3))
}

# Get variance explained for individual markers
hits$traitname = pheno_name_key[hits$Trait] # Convert back to raw phenotype label
hits$marker_r2 = NA
for(i in 1:nrow(hits)){
    myresid = residuals[[hits$traitname[i]]]
    mysnp = hits$Marker[i]
    hits$marker_r2[i] = get.r2(myresid, mysnp)
}

# Get variance explained by groups of SNPs
snp_groups = unique(hits[,c("Trait", "traitname")])
snp_groups$combined_r2 = snp_groups$rmip_r2 = snp_groups$glm_r2 = NA
for(i in 1:nrow(snp_groups)){
    # Pull out residuals and indices of trait
    myresid = residuals[[snp_groups$traitname[i]]]
    is_trait = hits$Trait == snp_groups$Trait[i]
    
    # GLM SNPs
    glm_snps = hits$Marker[ is_trait & hits$glm_perm_p <= glm_max_p]
    snp_groups$glm_r2[i] = get.r2(myresid, na.omit(glm_snps))
    
    # FarmCPU SNPs
    farm_snps = hits$Marker[ is_trait & hits$farm_RMIP >= farm_min_rmip]
    snp_groups$rmip_r2[i] = get.r2(myresid, na.omit(farm_snps))
    
    # Both together
    combined_snps = hits$Marker[is_trait]
    snp_groups$combined_r2[i] = get.r2(myresid, na.omit(combined_snps))
}

# Output
hits$traitname = NULL
snp_groups$traitname = NULL
write.csv(hits, file="Table - Combined hits.full.r_squared.csv", row.names=F)
write.csv(snp_groups, file="Table - Combined hits.r_squared_by_group.csv", row.names=F)