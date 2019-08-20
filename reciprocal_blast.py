#!/usr/bin/env python3.6

import sys, os, subprocess, re, glob
sys.path.append("/newhome/ab17362/Python_modules")
from my_module import *

#takes a fasta as input and loops through the blast databases in a directory, checking if the best hit is for each gene in the gene family
#corresponds to the representative of that species that is present in the fasta file

if len(sys.argv) != 4:
    print("Usage:  python3 reciprocal_blast.py fasta_file db_directory output directory\n\ndb_directory can be full or relative path. We're easy")
    exit()
else:
    seq_file = sys.argv[1]
    db_directory = sys.argv[2]
    output_directory = sys.argv[3]

seqs = get_file_data(seq_file)
seq_name = seq_file.split("/")[-1]
files = glob.glob(db_directory + "/*pin")
dbs = sorted(set([i.split("/")[-1].split(".", 1)[0] for i in files]))
fail = 0
for db in dbs:
    subprocess.run("blastp -query " + seq_file + " -db " + db_directory + "/"  + db + " -out " + output_directory + "/" + seq_name + "_vs_" + db + " -outfmt 6 -evalue 1e-09 -num_alignments 1 -max_hsps 1", shell = True)
