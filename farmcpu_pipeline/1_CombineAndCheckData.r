#! /usr/bin/Rscript

# This is a basic template to use for making R scripts
library(argparse)
options(stringsAsFactors=F)

# Arguments
parser=ArgumentParser()
parser$add_argument("-g", "--genofile", help="Input genotype file (numeric matrix, RDS format)")
parser$add_argument("-p", "--phenofiles", nargs="*", help="Input files of phenotyes")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
parser$add_argument("-c", "--covariates", nargs="*", help="Which traits are the covariates")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Sorghum/Punurri_aphid_resistance/FarmCPU-RMIP-analysis/0_Inputs_corrected/')  # Debug working directory
# args=parser$parse_args(c("-p","aphid damage rating AV 2019&2020.txt", "Dron data .txt", "Greenhouse data .txt",
#                          "-g",'ImputedGAPIT.Genotype.Numerical.rds',
#                          "-c", "Log_2019.2020.Flowering.time", "-o",'99_tmp'))


# Load cat
cat("Loading and comparing genetic and phenotype data\n")
phenos = lapply(args$phenofiles, read.delim, row.names="Taxa")
genos = readRDS(args$genofile)

# Join phenos
for(i in 1:length(phenos)){
    phenos[[i]] = phenos[[i]][order(rownames(phenos[[i]])),, drop=F]
    
    if(!identical(rownames(phenos[[i]]), rownames(phenos[[1]]))){
        cat("\tWarning: Taxa names in file",i,"do not match the first file!")
    }
}
phenos = do.call(cbind, phenos)

# Compare genos and phenos
display = function(a, b, text){
    percent = round(length(a) / length(b) * 100, digits=1)
    cat("\t",length(a),"out of",length(b),text, "taxa are shared (",percent,"% )\n" )    
}
shared = intersect(rownames(phenos), genos$taxa)
display(shared, rownames(phenos), "trait")
display(shared, genos$taxa, "genotype")


# Write out combined phenotypes
covariate = subset(phenos, select = names(phenos) %in% args$covariates)
covariate = data.frame(Taxa=rownames(covariate), covariate)
write.csv(covariate, file=paste(args$outprefix, ".covariates.csv", sep=""), row.names = F)

#outphenos = subset(phenos, select = !names(phenos) %in% args$covariates)
outphenos = phenos  # Keeping all phenotypes
outphenos = data.frame(Taxa=rownames(outphenos), outphenos)
write.csv(outphenos, file=paste(args$outprefix, ".phenos.csv", sep=""), row.names = F)