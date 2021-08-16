#! /usr/bin/Rscript

# Compile permutations and resamples to get RMIPs
library(argparse)
library(ggplot2)
options(stringsAsFactors=F)

# Arguments
parser=ArgumentParser()
parser$add_argument("-p", "--perms", nargs="*", help="List of input files of GAPIT results (permutations)")
parser$add_argument("-r", "--resamples", nargs="*", help="List of input files of GAPIT results (resample iterations)")
parser$add_argument("--perm-p-cutoff", type='double', default=0.01, help="P-value cutoff (determined from permutations) to be considered a hit in a resample run")
parser$add_argument("--rmip-cutoff", type='double', default=0.05, help="RMIP cutoff to be consider a genuine hit")
parser$add_argument("-o", "--outprefix", help="Output file prefix")
parser$add_argument("-t", "--trait", default="unnamed trait", help="Trait name")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/Sorghum/Punurri_aphid_resistance/FarmCPU-RMIP-analysis/')
# args=parser$parse_args(c("-p","1b_Permutations/Log_D1_2019.2020/perm1/GAPIT.FarmCPU.Log_D1_2019.2020.GWAS.Results.csv.gz",
#                        "1b_Permutations/Log_D1_2019.2020/perm1/GAPIT.FarmCPU.Log_D1_2019.2020.GWAS.Results.csv.gz",
#                        "-r", "1d_Resampling/Log_D1_2019.2020/resample1/GAPIT.FarmCPU.Log_D1_2019.2020.GWAS.Results.csv.gz",
#                        "1d_Resampling/Log_D1_2019.2020/resample1/GAPIT.FarmCPU.Log_D1_2019.2020.GWAS.Results.csv.gz",
#                        "-o", "99_tmp"))


cat("Compiling permutations and resamplings for",args$outprefix,"\n")

# Load permutations
cat("\tLoading permutations from",length(args$perms),"input files\n")
perms = lapply(args$perms, read.csv, nrow=1) # GAPIT best hit is always the top row
perms = sort(do.call(rbind, perms)$P.value)

# Load resamples; process at this stage so as to not save a bunch of unnecessary data in memory
cat("\tLoading resamples from",length(args$resamples),"input files\n")
resamples = lapply(args$resamples, function(infile){
    mydata = read.csv(infile)
    
    # # Get empirical p-value
    # mydata$empirical_p = sapply(mydata$P.value, function(p){
    #     sum(p>=perms)/length(perms)
    # })
    
    # Below is equivalent to above; need to determine which is faster
    mydata$empirical_p = findInterval(mydata$P.value, vec=perms) / length(perms)
    
    
    # Filter based on input arguments
    mydata=subset(mydata, mydata$empirical_p <= args$perm_p_cutoff)
    return(mydata[,c('SNP', 'Chromosome', 'Position', 'empirical_p')])
})

# Get RMIPs
cat("\tCalculating RMIPs\n")
resamples.all = do.call(rbind, resamples)
rmip = paste(resamples.all$SNP, resamples.all$Chromosome, resamples.all$Position, sep="###")
rmip = table(rmip) / length(args$resamples)
cat("\t\tTotal of", sum(rmip >= args$rmip_cutoff), "traits pass RMIP cutoff\n")

# Separate out for plotting
if(length(rmip)==0){
    toplot = data.frame(snp=character(), chr=numeric(), pos=numeric(), rmip=numeric())
}

toplot = strsplit(names(rmip), split="###")
toplot = data.frame(do.call(rbind, toplot))
names(toplot) = c("snp", "chr", "pos")
toplot$chr = as.numeric(toplot$chr)
toplot$pos = as.numeric(toplot$pos)
toplot$rmip = as.numeric(rmip)

    
# Graph results
toplot$significant = ifelse(toplot$rmip >= args$rmip_cutoff, yes="yes", no="no")
myplot = ggplot(toplot, mapping=aes(x=pos, y=rmip, color=significant)) +
    geom_point() + 
    facet_wrap(~chr, nrow=2) +
    ggtitle(args$trait) +
    scale_color_manual(values=c("yes"="red", "no"="black"))

# Output 
cat("\tWriting out results\n")
toplot = toplot[order(toplot$rmip, toplot$chr, toplot$pos),]
write.csv(toplot, file=paste(args$outprefix, ".csv", sep=""), row.names=F)
ggsave(myplot, file=paste(args$outprefix, ".png", sep=""), width=10, height=6, dpi=300)
