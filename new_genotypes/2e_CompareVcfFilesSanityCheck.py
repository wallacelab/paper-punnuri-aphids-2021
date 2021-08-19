__author__ = 'jgwall'

import argparse
import gzip

debug = False

chromID, posID, nameID, refID, altID = 0, 1, 2, 3, 4

def main():
    args = parse_args()
    print("Comparing sites and alleles calls in",args.infile,"versus",args.reffile)

    ref = load_sites(args.reffile)
    tocheck = load_sites(args.infile)
    report = compare_sites(ref, tocheck, args.reffile, args.infile)

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


def load_sites(infile):
    print("\tLoading sites and alleles from", infile)
    sitekey = dict()    # To store information
    for line in get_filehandle(infile, "rt"):
        if line.startswith("#"):
            continue
        data = line.strip().split('\t')

        # Extract needed info
        mychrom, mypos, myname = data[chromID], data[posID], data[nameID]
        myref, myalt = data[refID], data[altID]

        # Make a key
        mykey = "|".join([mychrom, mypos, myname])
        if mykey in sitekey:
            print("\t\tWARNING: Site key",mykey,"already exists and will be overwritten!")

        # Save data
        sitekey[mykey] = dict(ref=myref, alt=set(myalt.split(',')))
        if debug: print("Stored:", sitekey[mykey])

    print("\t\tLoaded",len(sitekey),"sites")
    return sitekey


def compare_sites(ref, tocheck, reffile, tocheckfile):
    print("\tComparing sites and alleles")

    # String to store everything in (to be printed and written to a file)
    report = "Reference file: " + reffile + "\n"
    report += "Query file: " + tocheckfile + "\n"

    # Compare sites
    check_keys , ref_keys = set(tocheck.keys()), set(ref.keys())
    notfound_1 = check_keys - ref_keys
    notfound_2 = ref_keys - check_keys
    shared_sites = ref_keys & check_keys
    report += "-- Sites --\n"
    report += str(len(notfound_1)) +  " sites are found in the query but not in the reference\n"
    report += str(len(notfound_2)) + " sites are found in the reference but not in the query\n"
    report += str(len(shared_sites)) + " sites are shared\n"

    # Compare alleles
    ref_different, alt_different, compatible, incompatible = 0, 0, 0, 0
    incompatible_sites = list()
    for mysite in shared_sites:
        ref1 = tocheck[mysite]['ref']
        ref2 = ref[mysite]['ref']

        alt1 = tocheck[mysite]['alt']
        alt2 = ref[mysite]['alt']

        # Check if reference alleles match
        if ref1 != ref2 : ref_different +=1
        # Check if alternate alleles match
        if len(alt1 - alt2) != 0 : alt_different +=1

        # Check if sites are compatible (=same collection of ref and alt alleles, even if in different order)
        check_alleles = alt1 | set(ref1)
        ref_alleles = alt2 | set(ref2)

        # Remove "no alternate" allele coding
        if "." in check_alleles: check_alleles.remove(".")
        if "." in ref_alleles:   ref_alleles.remove(".")
        #print(check_alleles, ref_alleles, alt1, alt2)

        # Calculate if are incompatible
        if len(ref_alleles - check_alleles)>0 & len(check_alleles - ref_alleles)>0:  # Only incompatible if both directions have leftover alleles
            incompatible +=1
            incompatible_sites.append(mysite + "\t" + ref1 + "|" + ",".join(alt1) + "\t" + ref2 + "|" + ",".join(alt2))
        else:
            compatible +=1

    # Add to report
    report += "-- Alleles --\n"
    report += str(ref_different) + " sites have different reference alleles\n"
    report += str(alt_different) + " sites have different alternate alleles\n"
    report += str(compatible) + " sites are compatible (=same collection of total alleles)\n"
    report += str(incompatible) + " sites are incompatible\n"


    # Print out the report string
    print("### Sanity check report ###")
    print(report)

    # Add the list of incompatible sites
    report += "-- Incompatible sites --\n"
    report += "\n".join(incompatible_sites)

    return(report)


def get_filehandle(file, mode):
    if file.endswith(".gz"):
        return gzip.open(file, mode)
    else:
        return open(file, mode)


if __name__ == '__main__': main()