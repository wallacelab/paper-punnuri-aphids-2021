__author__ = 'jgwall'

import argparse
import gzip
import pandas as pd
import numpy as np

debug = False

chromID, posID, nameID, refID, altID = 0, 1, 2, 3, 4
firstSampleID = 9

def main():
    args = parse_args()
    print("Comparing genotype calls in",args.infile,"versus",args.reffile)

    ref = load_genotypes(args.reffile)
    tocheck = load_genotypes(args.infile)
    report = compare_genotypes(ref, tocheck, args.reffile, args.infile)

    print("Saving report to",args.outfile)
    open(args.outfile, "w").write(report)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="VCF file to perform sanity check on")
    parser.add_argument("-r", "--reffile", help="VCF file to use as the standard / to compare against")
    parser.add_argument("-o", "--outfile", help="File to write report to")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def load_genotypes(infile):
    print("\tLoading genotypes from", infile)

    # Go down past VCF header information
    IN = get_filehandle(infile, "rt")
    myline = IN.readline()
    while myline.startswith("##"):
        myline = IN.readline()      # myline should end with the header line starting with #CHROM

    # Build genotypes from VCF
    header= myline.lstrip("#").strip().split('\t')
    genos = pd.read_csv(IN, sep='\t', header=None, names=header, nrows = 1000 if debug else None)

    # Print progress
    print("\t\tVCF file has", len(genos),"sites and ", len(genos.columns) - firstSampleID, "samples")

    return genos


def compare_genotypes(ref, tocheck, reffile, tocheckfile):
    print("\tComparing sites and alleles")

    # String to store everything in (to be printed and written to a file)
    report = "Reference file: " + reffile + "\n"
    report += "\t" + str(len(ref)) + " sites and " + str(len(ref.columns) - firstSampleID) + " samples\n"
    report += "Query file: " + tocheckfile + "\n"
    report += "\t" + str(len(tocheck)) + " sites and " + str(len(tocheck.columns) - firstSampleID) + " samples\n"

    # Determine sites and samples in common
    ref_samples, tocheck_samples = set(ref.columns[firstSampleID:]), set(tocheck.columns[firstSampleID:])
    ref_sites, tocheck_sites = set(ref["ID"]) , set(tocheck["ID"])
    common_samples =  ref_samples & tocheck_samples
    common_sites = ref_sites & tocheck_sites
    unshared_samples = (ref_samples | tocheck_samples) - common_samples
    unshared_sites = (ref_sites | tocheck_sites) - common_sites
    report += "-- Commonalities --\n"
    report += str(len(common_samples)) + " samples in common\n"
    report += str(len(common_sites)) + " sites in common\n"

    # Get common sites and samples
    ref_calls = np.array(ref.loc[ [id in common_sites for id in ref['ID']] , common_samples])
    tocheck_calls = np.array(tocheck.loc[[id in common_sites for id in tocheck['ID']] , common_samples])
    report += "Shape of reference after culling is " + str(ref_calls.shape) + "\n"
    report += "Shape of query after culling is " + str(tocheck_calls.shape) + "\n"

    # Remove "|"  and "/" from genotypes (=phased and unphased) because otherwise everything is different
    for i in range(ref_calls.shape[0]):
        for j in range(ref_calls.shape[1]):
            ref_calls[i,j] =ref_calls[i,j].replace("|", "").replace("/","")
            tocheck_calls[i, j] = tocheck_calls[i, j].replace("|", "").replace("/", "")

    # Compare genotypes
    same = np.array(ref_calls) ==  np.array(tocheck_calls)
    same_count = np.sum(same)
    different_count = np.sum(~same)

    # Different genotypes
    old = ref_calls[~same].flatten()
    new = tocheck_calls[~same].flatten()
    different_genos = pd.DataFrame({"old":old, "new":new})

    report += "-- Genotypes --\n"
    report += "Genotypes the same after culling is " + str(same_count) + "\n"
    report += "Genotypes different after culling is " + str(different_count) + "\n"

    # Print out the report string
    print("### Genotype check report ###")
    print(report)

    # Add the list of unshared sites and samples
    report += "-- Unshared samples --\n"
    report += "\n".join(unshared_samples) + "\n"
    report += "-- Unshared sites --\n"
    report += "\n".join(unshared_sites) + "\n"
    report += "-- Different genotypes (unique combinations of changes) --\n"
    #report += str(different_genos.drop_duplicates())
    report += str(different_genos.groupby(['old','new']).size().reset_index(name='Count'))

    #print(different_genos)

    return(report)


def get_filehandle(file, mode):
    if file.endswith(".gz"):
        return gzip.open(file, mode)
    else:
        return open(file, mode)


if __name__ == '__main__': main()
