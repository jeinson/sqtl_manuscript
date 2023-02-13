#!/gpfs/commons/home/jeinson/python3_env/bin/python

import argparse
import pysam
import pandas
import sys 

def flush_print(text):
    print(text)
    sys.stdout.flush()

def fun_dict_from_vcf_info(info):
    info = info.split(";")
    dict_info = {}
    for item in info:
        fields = item.split("=")
        if len(fields) == 2:
            tag = fields[0]
            value = fields[1]
            dict_info[tag] = value
    return dict_info

def main():

    # Arguments passed
    parser = argparse.ArgumentParser()

    # required
    parser.add_argument("--sQTLs", help="A file containing sQTLs variants which you want to find ancestral alleles")
    parser.add_argument("--anc_allele", help="Phased genotype VCF with all variants")

    global args
    args = parser.parse_args()

    # Open a tabix object for the ancestral allele file
    tabix_anc_allele = pysam.Tabixfile(args.anc_allele, "r")

    # Read in the eqtls as a pandas object
    sqtls = open(args.sQTLs, 'r')

    header = next(iter(sqtls)).strip().split("\t")

    # Add the new columns onto header, and flush print it
    header.extend(['ref_allele', 'alt_allele', 'alt_allele_freq', 'anc_allele', 'anc_allele_freq', 'hi_allele', 'li_allele', 'hi_allele_freq'])
    flush_print("\t".join(header))

    # Get ID index and slope index
    id_ix = header.index("ID")
    slope_ix = header.index("slope")

    chrom_old = "none"
    for line in sqtls:
        # Parse out the sQTL variant
        var = line.strip().split()[id_ix]

        chrom, pos, ref, alt = var.split("_")[0:4]

        # Get the gnomad MAF
        if chrom != chrom_old:
            gnomad_fp = "/gpfs/commons/datasets/gnomAD/3.0/gnomad.genomes.r3.0.sites." + chrom + ".vcf.bgz"
            # Open a tabix object for the gnomad file
            tabix_gnomad = pysam.Tabixfile(gnomad_fp, "r")
            chrom_old = chrom

        for record in tabix_gnomad.fetch(chrom, int(pos)-1, int(pos)):
            record = record.split("\t")
            if record[4] == alt:
                gnomad_af = float(fun_dict_from_vcf_info(record[7])['AF'])
            # else:
            #     import code; code.interact(local=locals()); sys.exit()
            

        # Parse out ancestral allele
        anc_allele = tabix_anc_allele.fetch(chrom[3:], int(pos)-1, int(pos))
        try:
            anc_allele = next(iter(anc_allele)).split("\t")[6]
        except:
            anc_allele = "NA"

        # Get the ancestral allele frequency
        if alt == anc_allele:
            anc_allele_freq = gnomad_af
        elif ref == anc_allele:
            anc_allele_freq = 1-gnomad_af
        else:
            anc_allele_freq = "NA"

        # Determine if the ref is the higher or lower included allele
        if float(line.strip().split()[slope_ix]) > 0:
            li_allele, hi_allele = ref, alt
            hi_allele_freq = gnomad_af
        else:
            hi_allele, li_allele = ref, alt  
            hi_allele_freq = 1-gnomad_af
        
        # Print the output
        flush_print(
            line.strip() + "\t" + "\t".join(map(str, [ref, alt, gnomad_af, anc_allele, anc_allele_freq, hi_allele, li_allele, hi_allele_freq]))
        )
        #import code; code.interact(local=locals()); sys.exit()

if __name__ == "__main__":
    main()