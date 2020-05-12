import re


def extract_barcodes(sam_file, out_file,
                     first_bc_ref_start=46, second_bc_ref_start=161,
                     first_bc_length=26, second_bc_length=26):

    '''
    Extracts two barcodes from the read on each line of
    a SAM file, based on the CIGAR string.

    Parameters
    ----------
    sam_file : str
        The name of a sam file
    out_file : str
        The desired name of the produced output file
    first_bc_ref_start : integer, optional
        Where in the reference sequence the first (leftmost) barcode starts (default is 46)
    second_bc_ref_start : integer, optional
        Where in the reference sequence the second (tightmost) barcode starts (default is 161)
    first_bc_length : integer, optional
        Length of the second barcode (default is 26)
    second_bc_length : integer, optional
        Length of the first barcode (default is 26)

    Note on barcodes
    ----------------

    Note, the first barcode (left most in the read) has historically
    been referred to as BC2, and has been used as the high complexity
    barcode, while the second barcode has been referred to BC1, and
    has been used as the low complexity barcode. While the name change
    is somewhat confusing, I did this to make the code generic for the
    extraction of any two barcodes from a read

    Returns
    -------
    Nothing

    Output
    ------

    Tab-delimited, with the following columns:

    seq_id, first_barcode_sequence, first_barcode_quality, second_barcode_sequence, second_barcode_quality

    where the qualities are the extracted quality scores for the
    residues in the barcodes.

    Also, it will print a report of the barcode extraction, indicating
    how many barcodes were extracted with the expected or unexpected
    lengths, and how many reads did not align

    Usage Examples
    --------------

    extract_barcodes("align.sam", "cigar_break.txt")
    extract_barcodes("align.sam", "cigar_break.txt", 46, 161, 26, 26)

    '''

    # variables we'll use to keep track of some statistic about barcode extraction
    
    count_no_alignment = 0
    count_diff_size = 0
    count_corr_size = 0
    
    out_file_h = open(out_file, 'w')

    with open(sam_file) as handle_sam:

        for line in handle_sam:

            read_list = line.rstrip().split('\t')

            if len(read_list) >= 11:

                seq_id, start, cigar, sequence, quality = read_list[0], read_list[3], read_list[5], read_list[9], read_list[10]

            pos = 0  # current position in the read

            # start of the read alignment
            start_read = int(start) - 1

            # where leftmost barcode starts and ends in the alignment
            first_bc_read_start = first_bc_ref_start - start_read
            first_bc_read_end = first_bc_read_start + first_bc_length

            # where rightmost barcode starts and ends in the alignment
            second_bc_read_start = second_bc_ref_start - start_read
            second_bc_read_end = second_bc_read_start + second_bc_length

            if cigar == '*':  # skip line if read wasn't aligned
                count_no_alignment += 1
                continue

            # determine the locations of the barcodes in the read, based on the CIGAR operations

            cigar_ops = re.compile("([M,I,D])").split(cigar)

            for i in range(0, len(cigar_ops) - 1, 2):

                operation = cigar_ops[i+1]
                op_size = int(cigar_ops[i])

                if operation == 'M':

                    # if it's matches or mismatches, simply advance
                    # our position along the read

                    pos += op_size

                else:

                    # use a lambda to encode that insertions
                    # require addition, while deletions require
                    # subtraction from the barcode start and ends

                    if operation == 'I':

                        def combine(a, b): return (a+b)

                    elif operation == 'D':

                        def combine(a, b): return (a-b)

                    # adjust barcode start and ends, depending on
                    # our current position in the read

                    if pos < first_bc_read_start:

                        first_bc_read_start = combine(first_bc_read_start, op_size)

                    if pos < first_bc_read_end:

                        first_bc_read_end = combine(first_bc_read_end, op_size)

                    if pos < second_bc_read_start:

                        second_bc_read_start = combine(second_bc_read_start, op_size)

                    if pos < second_bc_read_end:

                        second_bc_read_end = combine(second_bc_read_end, op_size)

                    # update the position in the read, based on
                    # the number of insertions/deletion

                    pos = combine(pos, op_size)

            # extract barcode sequences and their qualities

            first_bc_sequence = sequence[first_bc_read_start:first_bc_read_end]
            first_bc_quality = quality[first_bc_read_start:first_bc_read_end]

            second_bc_sequence = sequence[second_bc_read_start:second_bc_read_end]
            second_bc_quality = quality[second_bc_read_start:second_bc_read_end]

            if len(first_bc_sequence) != first_bc_length or len(second_bc_sequence) != second_bc_length:
                count_diff_size += 1
            else:
                count_corr_size += 1
            
            out_file_h.write("\t".join([seq_id, first_bc_sequence, first_bc_quality, second_bc_sequence, second_bc_quality]) + '\n')

    out_file_h.close()

    print('number of unaligned read with * = ' + str(count_no_alignment))
    print('number reads with at least one BC with the wrong size = ' + str(count_diff_size))
    print('number reads with both BC with the correct size = ' + str(count_corr_size))
