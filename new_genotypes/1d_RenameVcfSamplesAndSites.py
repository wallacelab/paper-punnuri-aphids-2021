__author__ = 'jgwall'

import argparse
import gzip

debug = False

chromID, posID, nameID = 0, 1, 2    # Columns in VCF file with needed information


def main():
    args = parse_args()
    print("Renaming samples and sites in",args.infile)

    # Load keyfile
    samplekey = load_sample_keys(args.sample_keyfile)

    # Go through and update VCF file
    OUT = get_filehandle(args.outfile, "wt")
    changed=0
    for line in get_filehandle(args.infile, "rt"):

        # Skip header lines unchanged
        if line.startswith("##"):
            OUT.write(line)
            continue

        # Change sample names in data header
        elif line.startswith("#CHROM"):
            new_header = rename_samples(line, samplekey)
            OUT.write(new_header)

        # Generate site name based on chromosome and position
        else:
            data=line.strip().split('\t')
            chrom, pos = data[chromID], data[posID]
            
            # Rename chromosomes and sites
            newchrom = chrom.lstrip('Chr0')
            name = "S" + newchrom + "_" + pos
            data[nameID] = name
            data[chromID] = str(newchrom)
            
            # Write out new data
            newdata  ="\t".join(data) + "\n"  
            OUT.write(newdata)
            changed+=1

        if debug and changed % 100000 == 0:
            print("\tGenerated site names for", changed, "total sites so far")

    print("\tGenerated site names for",changed,"total sites")
    OUT.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-k", "--sample-keyfile", help="Two-column, tab-separated file of old and new sample names")

    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


def get_filehandle(file, mode):
    if file.endswith(".gz"):
        return gzip.open(file, mode)
    else:
        return open(file, mode)


def load_sample_keys(infile):
    print("\tLoading sample rename key from",infile)
    key = dict()

    # Go through file and get names
    for line in get_filehandle(infile, "rt"):
        old, new = line.strip().split('\t')
        if old in key:
            print("\tWARNING! Old file name",old,"already in key and will be overwritten!")
        key[old]=new

    # Return completed key
    print("\t\tLoaded",len(key),"names to change")
    return key


def rename_samples(line, samplekey):
    print("\tRenaming samples")
    data=line.strip().split('\t')

    # Check each field for a matching one
    changed=0
    for i in range(len(data)):
        if data[i] in samplekey:
            oldname = data[i]
            if debug:
                print("\t\tChanging",oldname,"to",samplekey[oldname])
            data[i] = samplekey[oldname]
            changed+=1
    print("\t\tUpdated a total of",changed,"sample names")

    # Regenerate header row and return
    return "\t".join(data) + '\n'

if __name__ == '__main__': main()
