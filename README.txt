README

Vivian Chen's Notes on...

Gavin's Barcode Extraction

Order of operations:
1) barcode_extractor.py + barcode_runner.py
2) combine_UMIs.py
3) map_barcodes.py

Step 1: barcode_extractor.py

Parameters:
1) sam_file = Monica's sam_files....?
2) out_file = name...usually barcode or new_extraction.txt?
3) first_bc_ref_start : integer, optional. Where in the reference sequence the first (leftmost) barcode starts (default is 46)
4) second_bc_ref_start : integer, optional. Where in the reference sequence the second (rightmost) barcode starts (default is 161)
5) first_bc_length : integer, optional. Length of the second barcode (default is 26)
6) second_bc_length : integer, optional. Length of the first barcode (default is 26)

BC2 is the first barcode (leftmost)
BC1 is the second barcode (right most)

Output:
Tab-delimited, with the following columns:

seq_id, first_barcode_sequence, first_barcode_quality, second_barcode_sequence, second_barcode_quality

where the qualities are the extracted quality scores for the residues in the barcodes.

Also, it will print a report of the barcode extraction, indicating how many barcodes were extracted with the expected or unexpected lengths, and how many reads did not align





Step 2: combine_UMIs.py

Input files: 

1) meta_file = multiplex.txt file
2) fastq_file = fastq that matches the demultiplex of meta_file
3) barcode_file = new_extraction.txt name of output file

Questions:
Step 1) is there a reason why you separated the code rather than running it through main int he original .py file?

Step 1) so you used Monica's sam_file? The one with the weird ref? You don't do PEAR yourself.

Step 1) style question again. handle_sam.


Step 2) Gavin you create this multiplex variable, but what is happening...? I don't see if used again.

Step 2) style question, what does meta_file_h mean compared to meta_file what does the _h mean? Same with out_file_h.

