# This script will read in a set of extracted barcodes, and a file of
# all possible barcodes, and attempt to best determine to which
# sequence an extracted barcode matches

# Usage:
#
# python3 map_barcodes.py mapping_file extracted_barcode_file output_file
#
# e.g.
#
# python3 map_barcodes.py BC_Genotype_Mapping.txt 0_pool_0/table_file.txt 0_pool_0/remapped_barcodes.txt
#
# the mapping file is expected to contain all the known BC1_BC2
# combinations in the second column. The output file is the same as
# the input file of extracted barcodes, except for those barcodes that
# are within a levenshtein distance of 5 of a known double barcode,
# the barcode is replaced with the correct version. Those at a greater
# distance are not written to the new file.

import sys
from Levenshtein import distance as levenshtein_distance

known_bc_file = sys.argv[1]
extracted_bc_file = sys.argv[2]
output_file = sys.argv[3]

# read in the known barcodes and store in a list

known_barcodes = []

with open (known_bc_file, "r+") as known:

    header = known.readline()
    
    for line in known:

        data = line.rstrip().split('\t')

        barcode = data[1]

        bc1, bc2 = barcode.split('_');
        
        known_barcodes.append([bc1, bc2])


# now read through the extracted barcodes and match them up should
# really be a separate function - will clean up code later to extract
# out common utilities

output_file_h = open(output_file, 'w+')

with open (extracted_bc_file, "r+") as extracted:

    distances = []
    passed = 0
    
    for line in extracted:

        data = line.rstrip().split('\t')

        bc1, bc2 = data[3], data[4]

        min_dist = 999
        best_bc1 = ''
        best_bc2 = ''

        # now let's find the best matching set of barcodes
        
        for known in known_barcodes:

            known_bc1 = known[0]
            known_bc2 = known[1]

            edit_dist = levenshtein_distance(bc1, known_bc1) + levenshtein_distance(bc2, known_bc2)

            if edit_dist < min_dist:

                min_dist = edit_dist
                best_bc1 = known_bc1
                best_bc2 = known_bc2

                if edit_dist == 0:

                    break

        if min_dist <= 5:
            
            # now print out the data, but with the barcodes replaced
            # with the ones to which they map

            passed += 1

            data[3] = best_bc1
            data[4] = best_bc2
            
            output_file_h.write("\t".join(data))
            output_file_h.write("\n")
            
        distances.append(min_dist)

output_file_h.close()
        
def frequency_histogram(distances):
    '''
    This method will make a histogram of the frequency distribution of
    Levenshtein distances. It's not currently used, but I didn't want
    to remove the code
    '''
    
    # now make a frequency graph

    import matplotlib.pyplot as plt

    bin_list = range(min(distances), max(distances))
        
    plt.hist(x=distances, bins='auto', color='#0504aa')
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Levenshtein Distance')
    plt.ylabel('Frequency')
    plt.yscale('log', nonposy='clip')
    plt.title('Frequency Distribution of Levenshtein Distances')
    
    plt.show()
