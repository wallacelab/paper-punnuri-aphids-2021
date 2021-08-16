#! /bin/bash

# Run FarmCPU in GAPIT with resampling and permutation analysis

# Parameters
num_pcs=3
num_perms=100
num_resamples=100
resample_fraction=0.9


# Directories
srcdir=0_Inputs
permdir=1b_Permutations
resampledir=1d_Resampling
resultdir=1e_Results
if [ ! -e $permdir ] ; then mkdir $permdir; fi
if [ ! -e $resampledir ] ; then mkdir $resampledir; fi
if [ ! -e $resultdir ] ; then mkdir $resultdir; fi

# Input files
map=$srcdir/Genotype.map.txt
genos=$srcdir/281_Genotyp_Numeric.rds
damage="$srcdir/281_Log_D1&D2_2019&2020.txt"
greenhouse="$srcdir/281_GHD.txt"
drone="$srcdir/281_drone data.txt"
aphids="$srcdir/281_aphid count.txt"
flowering="$srcdir/281_Flowering time.txt"
spad="$srcdir/GHD_nSpad.txt"
height="$srcdir/Plant height .txt"

# Covariates
cov_flower="FLoweringtime"
cov_spad="nSpad"
cov_height="Log_2019.2020_PH"


###########
# Setup
###########

# Combine and check
Rscript 1_CombineAndCheckData.r --genofile $genos --phenofiles "$flowering" "$spad" "$height" "$damage" "$aphids" "$greenhouse" "$drone" --covariates $cov_flower $cov_spad $cov_height --outprefix 1_traits 

# Get list of traits
head -n1 1_traits.phenos.csv | cut -f2- -d, | tr ',' '\n' | sed -e 's|"||g' > 1a_traitlist.txt

###########
# Permutations for p-value cutoffs
###########

phenofile=1_traits.phenos.csv
covfile=1_traits.covariates.csv

# Loop over traits
while read trait; do
        
    # Set up Covariates.
    case $trait in
        "GHD") covariates="$cov_spad" ;;  # SPAD for greenhouse damage
        "NDVIAug25") covariates="$cov_flower $cov_height" ;;  # Flowering + plant height for the drone data
        "NDREAug25") covariates="$cov_flower $cov_height" ;;
        "SAVI_Aug25") covariates="$cov_flower $cov_height" ;;
        "nSpad") covariates="NONE" ;;  # None for the covariates (SPAD, flowering, height) so can see where they hit on their own
        "FLoweringtime") covariates="NONE" ;;  
        "Log_2019.2020_PH") covariates="NONE" ;;  
        *) covariates="$cov_flower" ;;  # Just flowering for everything else (=aphid counts and plant height)
    esac
    echo "Running trait $trait with covariates $covariates"
    
 
    # Output directories
    mypermdir=$permdir/$trait
    myresampledir=$resampledir/$trait
    if [ ! -e $mypermdir ]; then mkdir $mypermdir; fi
    if [ ! -e $myresampledir ]; then mkdir $myresampledir; fi
    
    # Run permutations & resamplings
    Rscript 1b_ResampleFarmCpu.r --genofile $genos --phenofile $phenofile --map $map --trait $trait --num-pcs $num_pcs --covfile $covfile --covariates $covariates  \
        --permdir $mypermdir --resampledir $myresampledir --num-perms $num_perms --num-resamples $num_resamples --resample-fraction $resample_fraction

    # Save actual permutation results with p-values; remove all other intermediate files
    pigz $mypermdir/*/GAPIT*.GWAS.Results.csv
    rm $mypermdir/*/*.pdf $mypermdir/*/*.csv
    
    # Save actual resample results with p-values; remove all other intermediate files
    pigz $myresampledir/*/GAPIT*.GWAS.Results.csv
    rm $myresampledir/*/*.pdf $myresampledir/*/*.csv
    
    
    # Collect permuted p-values
    myperms=$permdir/$trait/*/*.GWAS.Results.csv.gz
    myresamples=$resampledir/$trait/*/*.GWAS.Results.csv.gz   
    Rscript 1e_GetRmips.r --perms $myperms --resamples $myresamples --perm-p-cutoff 0.05 --rmip-cutoff 0.05 -o $resultdir/1e_$trait.p05.rmip05 -t $trait
    
    #break
done < 1a_traitlist.txt


# TODO: Get linkage among all the RMIP hits and see if any should be combined into single hits?
