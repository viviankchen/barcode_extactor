# This is a modification of Lucas' original script, and assume that
# there is only a single entry in the provided file; it does no
# demultiplexing, as it assumes that demultiplexing is already done

# use as follows:
#
# python3 combine_UMIs.py meta_file.txt fastq_file BC_file
#
# e.g. python3 combine_UMIs.py 0_multiplex.txt N716S518_S301_L004_R1_001.assembled.fastq new_extraction.txt 

import sys
import os

def main():

    # should add some error checking to make sure we get input files
    
    meta_file = sys.argv[1]
    fastq_file = sys.argv[2]
    barcode_file = sys.argv[3]

    with open(meta_file, 'r') as meta_file_h:

        # create directory and file for each multiplex pairs based
        # on the multiplex files that contains: name of sample \t forward
        # multiplex \t Reverse multiplex with (reversed complemented)
        # sequence compare to the primer

        (condition, replicate, generation, multiplex_f, multiplex_r)  = meta_file_h.readline().rstrip('\n').split(',')
    
        multiplex = multiplex_f + '_' + multiplex_r

        dir_name = condition + '_' + replicate + '_' + generation

        os.system('mkdir %s' % dir_name)
    
        out_file_h = open('%s/table_file.txt' % dir_name, 'w')


    combine_UMIs(out_file_h, fastq_file, barcode_file)
    
    out_file_h.close()

def combine_UMIs(out_file_h, fastq_file, barcode_file):
    '''

    This function will take a fastq file and a file of extracted
    barcodes, and reunite those barcode sequences with their UMIs, by
    virtue of the sequence IDs that are common to both files

    '''

    line_num = -3
    count_unmatch_id = 0
    
    with open(fastq_file, 'r+') as fastq, open(barcode_file, 'r+') as cigar_break:

        # read in a line from the barcode file, then read through the
        # fastq file until we hit the same sequence id. Because the
        # barcode file is in the same order (though not all sequences
        # got an entry in the barcode file) we only need to read
        # forwards in the fastq file until we find the correct
        # sequence. We never need to go back to the beginning of
        # either file
        
        for line in cigar_break:

            # want to rename these variables to make consistent with
            # barcode extraction code
                
            (seq_id, barcode_f, q_bar_f, barcode_r, q_bar_r) = line.strip('\n').split('\t')

            while True:

                # read in the four lines of the fastq record
                    
                seq_id_fastq = fastq.readline().strip('\n').split(' ')[0]
                seq_id_fastq = seq_id_fastq[1:] # Removes the '@' from the sequence ID

                seq = fastq.readline().strip('\n')
                sign = fastq.readline().strip('\n')
                quality = fastq.readline().strip('\n')

                line_num += 4 # line number of the first line of this fastq record
                
                if seq_id == seq_id_fastq:

                    umi_f, umi_quality_f, umi_r, umi_quality_r = extract_umis(seq, quality)
                    
                    break

                else:
                    count_unmatch_id += 1 # not yet reached the correct fastq record

            numb_line = str(line_num)

            w = str.join('\t',(numb_line,seq_id,umi_f,barcode_f,barcode_r,umi_r,umi_quality_f,q_bar_f,q_bar_r,umi_quality_r + '\n'))

            out_file_h.write(w)


def extract_umis(seq, quality):

    # coordinates that we will need for extracting UMIs, which start
    # and end the read.

    # TODO:
    # Note, in reads that are not the expected length, it's possible
    # the UMIs don't start at the correct coordinate, or that some of
    # the UMI is missing. Currently, we are ignoring this

    umi_length = 8
    
    end_read = len(seq)
    umi_r_start = end_read - umi_length

    # now extract the UMIs and multiplex sequences,
    # and their qualities
                    
    umi_f = seq[:umi_length]
    umi_quality_f = quality[:umi_length]

    umi_r = seq[umi_r_start:end_read]
    umi_quality_r = quality[umi_r_start:end_read]

    return umi_f, umi_quality_f, umi_r, umi_quality_r
    
main()
