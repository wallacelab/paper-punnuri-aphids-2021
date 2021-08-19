#! /bin/bash

# Merge my SNPs with the existing ones from Hu et al 2019

# Add BedOps package to file path
BEDOPS=/home/jgwall/Software/FileTools/BedOps
export PATH="$PATH:$BEDOPS" 
TASSEL5="perl /home/jgwall/Software/TASSEL/tassel-5-standalone/run_pipeline.pl -Xms10g -Xmx40g"

# Raw data from elsewhere
genome=$HOME/Projects/0_RawData/SorghumGenome/3.1_phytozome/Sbicolor_313_v3.0.fa
src=$HOME/Projects/0_RawData/Sequencing/2019_11_20_Sorghum_GBS_for_FVSU/
public_snps=$HOME/Projects/0_RawData/SorghumGermplasm/Hu_et_al_10k_GBS_genotypes/SNPs_imp.recode.vcf.bgzip.gz
public_sample_names=$HOME/Projects/0_RawData/SorghumGermplasm/Hu_et_al_10k_GBS_genotypes/sample_names.txt
extra_samples="PI257599 PI535794 PI535783"  # Samples that are outside the SAP but want to incldue b/c they're also ones in this set; useful for sanity checks/positive controls
metadata=0_BackgroundLit/Hu_et_al_Table_S1_Genotypes_and_origins.csv

# Specific to this analysis
chrom_rename_key=0_chrom_rename_key.txt
new_snps=1_CallGenotypes/1d_snp_calls.renamed.vcf.gz
numcores=7 # Number of cores/threads to use

workdir=2_MergeGenotypes
if [ ! -e $workdir ]; then mkdir $workdir; fi


###################
# CONDA ENVIRONMENT
###################

# Conda environment names. (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_bcftools=bcftools-1.9
conda_beagle=beagle-4.1-21Jan17 # Not quite the same version of BEAGLE used by Hu et al (Theirs was 27 Jul 2016, but that's not in Bioconda)

# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED


################
# Extract genotypes from the Sorghum Association Panel
################

conda activate $conda_bcftools

# # Get the sample names needed for subsetting
# Rscript 2a_IdentifySapSampleNames.r -i $public_sample_names -m $metadata --pattern "SAP" -o $workdir/2_target_samples.txt --extra-samples $extra_samples
# zcat $new_snps | grep "^#CHROM" | tr '\t' '\n' | tail -n +10 > $workdir/2_new_samples.txt  

# # Subset master VCF file
# bcftools view --samples-file $workdir/2_target_samples.txt --output-type z --output-file $workdir/2a_sap.vcf.gz $public_snps
# bcftools index $workdir/2a_sap.vcf.gz

# # Sort new SNPs so TASSEL can deal with them without throwing an error
# $TASSEL5 -SortGenotypeFilePlugin -inputFile $new_snps -outputFile $workdir/2a_new_snps.fix_order.vcf.gz




#####################
# Imputation Pipeline with BEAGLE
#####################

conda activate $conda_bcftools

imputedir=$workdir/2b_Imputation_BEAGLE
if [ ! -e $imputedir ] ; then mkdir $imputedir; fi

# # Split genotypes by chromsome to reduce computational loads while imputing
# for chr in $(seq 1 10); do
#     
#     # Split reference and new panels by chromosomes
#     bcftools view -o $imputedir/2b_sap_snps.$chr.vcf.gz --output-type z $workdir/2a_sap.vcf.gz $chr       
#     $TASSEL5 -vcf $workdir/2a_new_snps.fix_order.vcf.gz -separate $chr -export $imputedir/2b_new_snps.$chr.hmp.txt.gz  
# 
#     # Merge
#     $TASSEL5 -fork1 -vcf $imputedir/2b_sap_snps.$chr.vcf.gz -fork2 -h $imputedir/2b_new_snps.$chr.hmp.txt.gz \
#         -combine3 -input1 -input2 -mergeGenotypeTables \
#         -export $imputedir/2c_merged_snps.chr$chr.vcf -exportType VCF  \
#         -runfork1 -runfork2
#     
# 
#     # Select only biallelic sites (=ones with exactly 2 alleles; invariant sites and multiallelic ones filtered out)
#     bgzip $imputedir/2c_merged_snps.chr$chr.vcf # Have to use blocked GZIP to index
#     tabix $imputedir/2c_merged_snps.chr$chr.vcf.gz  # Make index for subsetting
#     bcftools view --min-alleles 2 --max-alleles 2 --output-type z -o $imputedir/2d_merged_snps_standard.chr$chr.vcf.gz $imputedir/2c_merged_snps.chr$chr.vcf.gz
#     
# #     break
# done


# # Run imputation 
# conda activate $conda_beagle
# for chr in $(seq 1 10); do
# 
#     infile=$imputedir/2d_merged_snps_standard.chr$chr.vcf.gz
#     outprefix=$imputedir/2e_imputed_snps.merged.$chr
# 
#     beagle -Xms10g -Xmx40g gt=$infile out=$outprefix nthreads=$numcores window=2000 overlap=500
# 
#     # Quick check of how many / if any sites got changed 
#     python3 2e_CompareVcfFilesSanityCheck.py -i $outprefix.vcf.gz --reffile $infile -o $imputedir/2e_imputed_snps.merged.$chr.site_and_allele_report.txt
#     python3 2e_CompareVcfFilesGenotypeCalls.py -i $outprefix.vcf.gz --reffile $infile -o $imputedir/2e_imputed_snps.merged.$chr.genotype_report.txt 
# 
# #     break
# done


# Merge back together
conda activate $conda_bcftools
files=""    # Build files up in chromosome order (requirement for bcftools concat)
for chr in $(seq 1 10); do
    files="$files $imputedir/2e_imputed_snps.merged.$chr.vcf.gz"
    bcftools index $imputedir/2e_imputed_snps.merged.$chr.vcf.gz
done
bcftools concat  --output-type z -o $workdir/2f_sap.combined_imputed.vcf.gz $files 


