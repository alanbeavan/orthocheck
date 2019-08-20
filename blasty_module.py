#!/usr/bin/env python3.6

import sys, re
sys.path.append("/newhome/ab17362/Python_modules")
from my_module import *

def seqs_are_different(hit, query, seqs):
    """for the identifier of the hit, the query and the database, check if the sequences the headers refer to are different"""
    r = re.compile(re.escape(hit))
    q = re.compile(re.escape(query))
    if len(list(filter(r.search, seqs.keys()))) == 1 and len(list(filter(q.search, seqs.keys()))) == 1:
        if seqs[list(filter(r.search, seqs.keys()))[0]] == seqs[list(filter(q.search, seqs.keys()))[0]]:
            return 0
        else:
            return 1
    elif len(list(filter(r.search, seqs.keys()))) >= 2 or len(list(filter(q.search, seqs.keys()))) >= 2:
            print("something is wrong - there are muliple headers in the reference fasta that match this hit or query")
            print("the query was " + query + "\nthe hit was " + hit)
    else:
        print("something is wrong - there are no headers in the reference fasta that match this hit or query")
        print("the query was " + query + "\nthe hit was " + hit)


class Results():
    """the blast results for a gene family describing the sequences hitting and some measure of how successful it is"""
    def __init__(self, name, hit_dict, hits, graph, species):
        """initialise results"""
        self.name = name
        self.hit_dict = hit_dict
        self.hits = hits
        self.graph = graph
        self.species = species

    def write_fasta(self, dir, seqs):
        """write a fasta file of all hits"""
        f = open(dir + "/" + self.name + "_all_blast_hits.fa", "w")
        for header in self.hits:
            r = re.compile(re.escape(header))
            f.write(">" + list(filter(r.search, seqs.keys()))[0] + "\n" + seqs[list(filter(r.search, seqs.keys()))[0]] + "\n")
        f.close()
            
    def is_fully_reciprocal(self, seqs):
        """confirm if the family is fully reciprocal"""
        for key, values in self.hit_dict.items():
            k = re.compile(re.escape(key))
            for value in values:
                v = re.compile(re.escape(value))
                if re.search(k, value) or re.search(v, key):
                    next
                elif seqs_are_different(value, key, seqs):
                    return 0
        return 1

    def add_to_list(self, filename):
        """add the family name to a file"""
        f = open(filename, "a")
        f.write(self.name + "\n")

    def largest_ortholog_set(self, seqs, dir):
        """return the largest possible set of reciprocal all v all genes from the family"""
        species = self.species.copy()
        coverred = []
        seq_names = {}
        for key in self.graph:
            #update the species coverred the the name of the orginal sequence
            coverred.append(key.split("_")[0].title())
            seq_names[key.split("_")[0].title()] = key

            #assess if the expected gene is matched by >50% of the hits
            #list of hits
            hits = []
            for value in self.graph[key].values():
                hits.append(value)
            
            #Check if the hits are reciprocal - adding one to the count if they are
            key_re = re.compile(re.escape(key))            
            count = 0
            for hit in hits:
                hit_re = re.compile(re.escape(hit))
                if re.search(hit_re, key) or re.search(key_re, hit):
                    count += 1
                elif not seqs_are_different(key, hit, seqs):
                    count += 1

            #Assess if the gene is matched reciprically by >50% of species
            target = len(self.species)/2
            print("species = " + key.split("_")[0].title())
            print("target = " + str(target))
            print("count = " + str(count))
            #ie. in confident families
            if count > target:  
                for hit in hits:
                    #if the sequence matches the expected - move on: recipricocity fulfilled
                    if re.search(hit_re, key) or re.search(key_re, hit):
                        next
                    #Otherwise the species that the query came from does not reciprically best hit with this species, so needs to be removed
                    else:
                        if hit.split("_")[0].title() in species:
                            print(hit + " removed because the hit didn't match the expected in a confident family")
                            species.remove(hit.split("_")[0].title())
                            if hit.split("_")[0].title() in seq_names:
                                del seq_names[hit.split("_")[0].title()]
            #if the expected was matched in fewer than 50% of hits, the species
            else:
                if key.split("_")[0].title() in species:
                    print(key + " removed because most genes didn't match it")
                    species.remove(key.split("_")[0].title())
                    if hit.split("_")[0].title() in seq_names:
                        del seq_names[hit.split("_")[0].title()] 

        print(coverred)
        #We now remove species from the results that were never in the original orthgroup
        if len(coverred) < len(self.species):
            for id in self.species:
                if id not in coverred:
                    if id in species:
                        print(id + " removed because it was not coverred in the results (ie nothing hit it and it wasn't in the original family")
                        species.remove(id)

        
        print(self.name + "\n" + " ".join(species) + "\n\n")
        #if the number of species satisfying reciprical all v all best hits is more than half the species, write the results
        if len(species) > len(self.species)/2:
            f = open(dir + "/" + self.name + "_putative_orthologs.fa", "w")
            for seq_name in seq_names:
                #look for the seq in the proteome and write it to a file
                to_search = seq_names[seq_name]
                r = re.compile(re.escape(to_search))
                if len(list(filter(r.search, seqs.keys()))) == 1: 
                    f.write(">" + seq_name + "\n" + seqs[list(filter(r.search, seqs.keys()))[0]] + "\n")
                else:
                    print("something is wrong - bad key: " + seq_names[seq_name])
            f.close()



    #I still need some sort of algorithm to decide which gene families to keep...
    #maybe
    #make max likelihood tree
    #decide if the genes are in-paralogs - ie. is there a bipartition that splits them from everything else with a reasonable level of support?
    #This wouldn't be for this program actually 
