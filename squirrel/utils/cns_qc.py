#!/usr/bin/env python3
import sys
import os
from squirrel.utils.log_colours import green,cyan
import select
from Bio import SeqIO
import collections
import csv
from squirrel.utils.config import *


def add_refs_to_input(input_fasta,assembly_refs,config):
    new_input_fasta = os.path.join(config[KEY_TEMPDIR], "qcref_combined.fasta")
    with open(new_input_fasta,"w") as fw:
        for record in SeqIO.parse(input_fasta,"fasta"):
            fw.write(f">{record.id}\n{record.seq}\n")
        for record in assembly_refs:
            fw.write(f">{record.id}\n{record.seq}\n")
    return new_input_fasta


def find_assembly_refs(cwd,assembly_refs,config):

    if not assembly_refs:
        sys.stderr.write(cyan(f'Error: no assembly references supplied.\nMust supply a reference file with one or more assembly references to do sequence QC.\n'))
        sys.exit(-1)

    refs = []
    ref_ids = []
    for record in SeqIO.parse(os.path.join(cwd,assembly_refs),"fasta"):
        refs.append(record)
        ref_ids.append(record.id)

    config[KEY_ASSEMBLY_REFERENCES] = ref_ids

    return refs


def get_seq_at_node(state_file,nodename):
    
    #returns a dict keys off 1-based positions in the genome
    #and the value is a list of tuples (id, base) at a given node
    # allows you to look up for a given site what the base is for a
    #given internal node or tip
    
    seq = ""
    with open(f"{state_file}","r") as f:
        for l in f:

            if not l.startswith("#"):
                try:
                    node,site,state,probA,probC,probG,probT = l.rstrip("\n").split("\t")
                except:
                    print(l)
                    break
                if node == nodename:
                    seq+=state

    return seq

def check_calls_to_reference(references):

    ref_vars = collections.defaultdict(set)
    ## get_reference_alleles
    for record in SeqIO.parse(alignment,"fasta"):
        for site in node_states:
            index = int(site)-1
            base = record.seq[index]
            if base in ["T","C","A","G"]:
                node_states[site].append((record.id,base))
            else:
                node_states[site].append((record.id,""))

def run_seq_qc(state_diff_file, aa_file, config):

    aligned_fasta = os.path.join(config[KEY_OUTDIR],config[KEY_OUTFILE])

    assembly_refs = config[KEY_ASSEMBLY_REFERENCES]
    print(green(f"{len(assembly_refs)} assembly references read in."))
    
    assembly_ref_shared = {}
    for i in assembly_refs:
        print(f"- {i}")
        assembly_ref_shared[i] = {}
    
    with open("shared_with_ref.csv","w") as fw: 
        with open(state_diff_file,"r") as f:
            reader = csv.DictReader(f)
            writer = csv.DictWriter(fw,header=reader.fieldnames,lineterminator="\n")
            writer.writeheader()
            for row in reader:
                root_var = reader["Node1"]
                ref_vars = {}
                for i in assembly_refs:
                    ref_vars[i] = row[i]

                for j in reader.fieldnames:
                    if j not in assembly_refs and j!= "site":
                        
                        if row[j] != root_var:
                            shared_with_ref = False
                            snp = row[j]
                            
                            for i in ref_var:
                                if snp == ref_vars[i]:
                                    shared_with_ref = True
                        
                        # if shared_with_ref:





