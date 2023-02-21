import argparse
import pysam
import numpy
import re
import sys

def main():

	parser = argparse.ArgumentParser()

	# required arguments
	parser.add_argument("psi_file", help="The raw output of the IPSA pipeline, where each exon is listed as its coordinate in the first column")
	parser.add_argument("annotation_gfx", help="The tabix indexed annotation file used in the IPSA run, but with more details")
	parser.add_argument("--contig", help="The chromosome number to use, for running this in parallel")

	#if  len(sys.argv) == 1:
		#parser.print_help(sys.stderr)
		#sys.exit(1)

	global args
	args, leftovers = parser.parse_known_args()


	# Get the tabix indexed annoatation file
	gfx_path = args.annotation_gfx
	tabix_gfx = pysam.Tabixfile(gfx_path, "r")

	psi_file_path = args.psi_file
	psi_file = open(psi_file_path)

	# Get contig
	contig = args.contig

	# Loop through the results file
	for line in psi_file:
		loc = line.rstrip().split("\t")[0].split("_")

		# Check if we're at the header line
		if not(bool(re.match("^chr", loc[0]))):
			continue
		
		if args.contig is not None:
			if not(bool(re.match(contig + "$", loc[0]))):
				continue

		matches = tabix_gfx.fetch(loc[0], int(loc[1]), int(loc[2]))
		for match in matches:
			match = match.split("\t")
			if (loc[1] == match[3]) & (loc[2] == match[4]):
				print(match[0] + ":" + match[3] + "-" + match[4] + "\t" + match[8])
				break # This fixes a weird issue where the are a few cases that double match





if __name__ == "__main__":
	main()
