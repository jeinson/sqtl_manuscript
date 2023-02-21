import argparse
import numpy as np
import re

def main():

    parser = argparse.ArgumentParser()

    # required arguments
    parser.add_argument("--input", help="The combined PSI file from")
    parser.add_argument("--output", help="The name of the output file to save")

    global args
    args = parser.parse_args()

    input_file = open(args.input, "r")
    output_file = open(args.output, "w")

    output_file.write(input_file.readline())
    next(input_file)

    for line in input_file:

        line_list = line.strip().split("\t")

        # Do some subsetting
        psi = np.array(line_list[1:])

        # if line_list[0] == "chr10_110503433_110503490_+":
        #     import code; code.interact(local=locals()); sys.exit()

        if sum(psi == "NA") > len(psi) / 2:
            continue

        if len(np.unique(psi[psi != "NA"])) < 5:
            continue
        
        output_file.write(line)

    input_file.close()
    output_file.close()

if __name__ == "__main__":
	main()
