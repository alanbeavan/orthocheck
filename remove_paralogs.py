#!/usr/bin/env python3.6

import sys, re, glob
sys.path.append("/newhome/ab17362/Python_modules")
from my_module import *
from blasty_module import *



def check_args():
    """check the program has been called properly and assign variables"""
    if len(sys.argv) == 3:
        seqs_file = sys.argv[1]
        blast_dir = sys.argv[2]
        return seqs_file, blast_dir
    else:
        print("Usage: python remove_paralogs.py all_sequences_file blast_results_directoy\n\nfull or relative paths are fine\n")
        exit()


#initialise variables
seqs_file, blast_dir = check_args()
seqs = read_fasta(seqs_file)
files = glob.glob(blast_dir + "/*")
families = []
species_list = []


#Update the lists of families and species
for file in files:
    family = file.split("/")[-1].split(".")[0]
    if family not in families:
        families.append(family)
    sp = file.split("/")[1].split("_")[2]
    if sp not in species_list:
        species_list.append(sp)


#for each family
for family in families:
    #initialse results class
    current_results = Results(family, dict(), [], dict(), species_list)
    
    #loop through blast hit files
    r = re.compile(re.escape(family + ".fa"))
    for file in list(filter(r.search, files)):
        lines = get_file_data(file)
        species = file.split("_")[-1]
        sp = re.compile(re.escape(species + "_"))
        expected = ""
        query_hit_dict = {}
        for line in lines:
            fields = line.split("\t")
            query = fields[0]
            query_hit_dict[query] = fields[1]
            
            if re.search(sp, fields[0]):
                expected = query
        current_results.graph[expected] = query_hit_dict
    
    #get largest set of orthologs possible (this currently also writes them to files - I will refactor this later
    current_results.largest_ortholog_set(seqs, "putative_orthologs")




