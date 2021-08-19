# FarmCPU & plotting code for Punnuri _et al_ 2021
This respository contains the code for running a resmapling analysis of FarmCPU on sorghum aphid resistance for Punnuri et al 2021 (in preparation). It also contains the code to plot the results in publication-quality figures.

## farmcpu_pipeline
This folder contains the FarmCPU resampling pipeline used to generate resample model inclusion probabilities (RMIPs) for the analysis. It is designed to be run in bash using the master script "0_RunFarmCpu.sh" and requires R and its packages argparse, GAPIT3, and ggplot2, along with the pigz (parallel gzip) command-line program.

## new_genotypes
Most of the genotype data for this analysis came from [Hu et al. 2019](https://dx.doi.org/10.3835/plantgenome2018.06.0044). This folder contains the scripts used to call genotypes on 10 additional lines that were included but not part of the previous data. Running it requires downloading the [Sorghum v3.1.1 genome](https://phytozome-next.jgi.doe.gov/pz/portal.html#!info?alias=Org_Sbicolor) from Phytozome and the raw FASTQ data on the 10 lines. (TODO: add link to NCBI.)

## publication_figures
This folder contains a series of R scripts to generate figures and tables for the final publication. It is designed to be run as a series of R scripts.
* Scripts starting wtih "0_" are support scripts and are not meant to be run directly but called from other scripts
* Scripts starting with 1_ 2_, etc are the actual scripts used to create the figures and tables. The working directory in each will need to be udpated for your setup, but otherwise they should be ready to run.
