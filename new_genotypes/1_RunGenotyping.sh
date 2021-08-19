#! /bin/bash

# Get SNPs for each of the samples

# Add BedOps package to file path
BEDOPS=/home/jgwall/Software/FileTools/BedOps
export PATH="$PATH:$BEDOPS" 

genome=$HOME/Projects/0_RawData/SorghumGenome/3.1_phytozome/Sbicolor_313_v3.0.fa
src=$HOME/Projects/0_RawData/Sequencing/2019_11_20_Sorghum_GBS_for_FVSU/
public_snps=$HOME/Projects/0_RawData/SorghumGermplasm/Hu_et_al_10k_GBS_genotypes/SNPs_imp.recode.vcf.bgzip.gz
sample_name_key=0_sample_rename_key.txt
numcores=7 # Number of cores/threads to use

workdir=1_CallGenotypes
if [ ! -e $workdir ]; then mkdir $workdir; fi


###################
# CREATE CONDA ENVIRONMENT
###################

# Conda environment names. (These are handled with Miniconda and are meant to make things reproducible by making identical environments)
conda_bcftools=bcftools-1.9

# # Create BCFtools Conda environment (for reproducibility; only has to be done once, so leave commented out most of the time.)
# # note: this requires the bioconda channel to be included; see https://bioconda.github.io/user/install.html#set-up-channels for instructions
# conda create -n $conda_bcftools bcftools

# Load functions required to be able to activate conda environments within scripts.
. $(conda info --root)/etc/profile.d/conda.sh   # Loads the functions required to activate conda; KEEP THIS COMMAND UNCOMMENTED
conda activate $conda_bcftools

###################
# Run Genotyping
###################

# Index genome for alignment (bwa version 0.7.12-r1039)
bwa index $genome 

# Align each sample to the genome
for sample in $src/*.fastq.gz; do
    name=`basename $sample`
    name=${name/_R1*/}
    echo "Processing alignments for $name"
    
#     Align to genome
    bwa aln -t $numcores -f $workdir/1_$name.sai $genome $sample

    # Convert to SAM
    bwa samse $genome $workdir/1_$name.sai $sample > $workdir/1_$name.sam
    
#     Convert to BAM & sort
    samtools view -h -S -b $workdir/1_$name.sam > $workdir/1_$name.bam
    samtools sort -@ $numcores $workdir/1_$name.bam $workdir/1_$name.sorted
    samtools index $workdir/1_$name.sorted.bam

    rm $workdir/1_$name.sam $workdir/1_$name.sai $workdir/1_$name.bam 
#     break
done

# Index genome for samtools
samtools faidx $genome  

# Get list of locations to call at; have to tweak the format
zcat $public_snps | vcf2bed | cut -f1-3 > $workdir/1a_target_locations.bed
cat $workdir/1a_target_locations.bed | sed -r -e "s/^/Chr0/" | sed -e "s/Chr010/Chr10/" > $workdir/1a_target_locations.corrected.bed

# Call SNPs on all samples (skip indels with --skip-variants b/c there are no indels in the Hu et al data and it causes issues trying to merge & impute them)
max_depth=1000
min_base_quality=20 
bcftools mpileup -d $max_depth -f $genome --regions-file $workdir/1a_target_locations.corrected.bed -Q $min_base_quality --output-type z --output $workdir/1b_raw_probs.vcf.gz  $workdir/*.sorted.bam  
bcftools call --multiallelic-caller --threads $numcores --skip-variants indels --output $workdir/1c_snp_calls.vcf.gz --output-type z $workdir/1b_raw_probs.vcf.gz 

# Replace sample, site, and chromosome names throughout the VCF file.
python3 1d_RenameVcfSamplesAndSites.py -i $workdir/1c_snp_calls.vcf.gz -o $workdir/1d_snp_calls.renamed.vcf.gz --sample-keyfile $sample_name_key # --debug


# # Make keys to correct sample and site names for when merged with existing key

# DEPRECATED zcat $workdir/1c_snp_calls.vcf.gz | grep "^Chr" | sed -r "s|^Chr(..)\t([0-9]+)\t.+|Chr\1\t\2\tS\1_\2|" | sed "s|\tS0|\tS|"  > $workdir/1d.site-rename-key.txt  # Capture chrom + position to make a site name compatible with existing data


