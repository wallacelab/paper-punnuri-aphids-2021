#! /usr/bin/Rscript

# Determine the sample names to be used for genotypes

# Libraries
library(argparse)

# Command-line arguments
options(stringsAsFactors=F)
parser=ArgumentParser()
parser$add_argument("-i", "--infile", help="Input file of sample names (1 column) from Hu et al 2019")
parser$add_argument("-m", "--metadata", help="Hu et al's file of metadata (CSV format)")
parser$add_argument("-p", "--pattern", help="Pattern to search for in the 'Status' column to identify germplasm names")
parser$add_argument("-e", "--extra-samples", nargs="*", default="", help="Additional sample names to capture that don't match the supplied pattern")
parser$add_argument("-o", "--outfile", help="Output file")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Sorghum/Punnuri_SAP_Genotypes/2_MergeGenotypes/')
# args=parser$parse_args(c("-i","/home/jgwall/Projects/0_RawData/SorghumGermplasm/Hu_et_al_10k_GBS_genotypes/sample_names.txt", "-p", "SAP", "-o", "99_tmp.txt", "-m", "../0_BackgroundLit/Hu_et_al_Table_S1_Genotypes_and_origins.csv"))


# Load data
cat("Identifying sample names to pull out of public dataset\n")
samples=scan(args$infile, what=character())
metadata = read.csv(args$metadata, check.names=F)[,c("Accession_ID", "Status", "Country", "Races")]


# Subset metadata file to just the samples specified by args$pattern OR the ones supplied as extras with args$extra_samples 
tokeep = grepl(metadata$Status, pattern=args$pattern) # by default just do pattern search
if(length(args$extra_samples)>1 || args$extra_samples != "") {    # If extra samples provided, search for them, too
    tokeep = tokeep | (metadata$Accession_ID %in% args$extra_samples)
}
targets = subset(metadata, tokeep)

# Process GBS sample names to find the SAP sample names; since are separted by periods, just take everything before the first period
samples_trimmed = sapply(strsplit(samples, split="\\."), function(x){x[1]})

# Find matches
isfound = targets$Accession_ID %in% samples_trimmed
notfound = ! isfound
cat("\tFound",sum(isfound),"out of",nrow(targets),"target samples; failed to find",sum(notfound),"targets\n")
if(sum(notfound)>=1){
    cat("\t\tSamples not found are",sort(targets$Accession_ID[notfound]),"\n")
}

# Output original GBS names of samples
cat("Writing sample names to", args$outfile,"\n")
output = samples[samples_trimmed %in% targets$Accession_ID]
write(output, file=args$outfile)
