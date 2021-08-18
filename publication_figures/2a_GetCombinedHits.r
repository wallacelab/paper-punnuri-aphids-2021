#! /usr/bin/Rscript

options(stringsAsFactors=F)

# Arguments
setwd('/home/jgwall/Documents/Papers/PUNNURI_Sorghum Aphid Resistance/GWAS_github/publication_figures/')
args=list()
args$farmcpu = list.files(path = "Results_FarmCPU", pattern="1e_", full.names=TRUE)
args$glm =  list.files(path = "Results_GLM", pattern="glm_results", full.names=TRUE)
args$rmip_cutoff = 0.05
args$perm_cutoff = 0.05
args$outprefix = "Table - Combined hits"
farm_trait_pattern=".+1e_(.+)\\.p..\\.rmip.+.csv" # File pattern for extracting FarmCPU trait name

# Load key connecting human-fiendly phenotype names to ones used in pipeline
source("0_PhenotypeKey.r")

# Load GLM results
cat("Loading GLM results from", length(args$glm), "input files\n")
glm = lapply(args$glm, function(infile){
    mydata = read.csv(infile)
    mydata$Trait = sub(mydata$Trait, pattern="&", repl=".", fixed=T)
    
    # Need to add "with_flowering" to drone traits with flowering covariate 
    if(grepl(infile, pattern="with_flowering")){
        mydata$Trait = paste(mydata$Trait, ".with_flowering", sep="")
    }
    
    return(mydata[,c("Trait", "Marker", "Chr", "Pos", "p", "perm_p")])
})
glm = do.call(rbind, glm)

# Load FarmCPU results
cat("Loading FarmCPU results from", length(args$farmcpu), "input files\n")
farm = lapply(args$farmcpu, function(infile){
    mydata = read.csv(infile)
    
    # Get trait
    trait=sub(infile, pattern=farm_trait_pattern, repl="\\1")
    
    mydata$Trait=trait
    
    # Rename columns
    names(mydata)[names(mydata)=="chr"] = "Chr"
    names(mydata)[names(mydata)=="pos"] = "Pos"
    names(mydata)[names(mydata)=="snp"] = "Marker"
    
    # Return data
    return(mydata)
})
farm = do.call(rbind, farm)

# Subset to high-quality ones
glm.hits = subset(glm, glm$perm_p <= args$perm_cutoff)
farm.hits = subset(farm, farm$rmip >= args$rmip_cutoff)

# Get unique trait-SNP combinations
combined = rbind(glm.hits[,c("Trait", "Marker", "Chr", "Pos")],
                 farm.hits[,c("Trait", "Marker", "Chr", "Pos")])
combined = unique(combined)
combined = combined[order(combined$Trait, combined$Chr, combined$Pos),]

# Add in GLM p-values 
sitekey = paste(combined$Trait, combined$Chr, combined$Pos, combined$Marker)
glmkey = paste(glm$Trait, glm$Chr, glm$Pos, glm$Marker)
glm.match = match(sitekey, glmkey)
combined$glm_p = glm[glm.match, "p"]
combined$glm_perm_p = glm[glm.match, "perm_p"]

# Add RMIP values
farmkey = paste(farm$Trait, farm$Chr, farm$Pos, farm$Marker)
farm.match = match(sitekey, farmkey)
combined$farm_RMIP = farm[farm.match, "rmip"]

# Correct phenotype names
pheno_name_key.inverse = setNames(names(pheno_name_key), pheno_name_key) # Flip which is name and which is value
combined$Trait = pheno_name_key.inverse[combined$Trait]

# Write out full set
write.csv(combined, file=paste(args$outprefix, ".full.csv", sep=""), row.names=F)

# Write out only unique positions
spots = split(combined, paste(combined$Marker, combined$Chr, combined$Pos))
spots = lapply(spots, function(s){
    data.frame(Marker=unique(s$Marker), Chr=unique(s$Chr), Pos=unique(s$Pos), Traits=paste(s$Trait, collapse=";"))
})
spots=do.call(rbind, spots)
spots = spots[order(spots$Chr, spots$Pos),]
write.csv(spots, file=paste(args$outprefix, ".unique_locations.csv", sep=""), row.names=F)