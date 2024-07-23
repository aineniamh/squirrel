#!/usr/bin/env python3
import sys
import os
from squirrel.utils.log_colours import green,cyan
import select
from Bio import SeqIO
import collections

from squirrel.utils.config import *

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
                


def run_seq_qc(seq_qc,cwd,assembly_refs,config):

    if seq_qc:

        aligned_fasta = os.path.join(config[KEY_OUTDIR],config[KEY_OUTFILE])

        refs = []
        for record in SeqIO.parse(os.path.join(cwd,assembly_refs),"fasta"):
            refs.append(record)

        print(green(f"{len(refs)} assembly references read in."))
        for i in refs:
            print(f"- {i.id}")

def add_refs_to_input(input_fasta,assembly_refs,config):
    new_input_fasta = os.path.join(KEY_TEMPDIR, "qcref_combined.fasta")
    with open(new_input_fasta,"w") as fw:
        for record in SeqIO.parse(input_fasta,"fasta"):
            fw.write(f">{record.id}\n{record.seq}\n")
        for record in SeqIO.parse(assembly_refs,"fasta"):
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

    config[KEY_ASSEMBLY_REFS] = ref_ids

    return refs

