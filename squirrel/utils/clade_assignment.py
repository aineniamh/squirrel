#!/usr/bin/env python3
import sys
import os
from squirrel.utils.log_colours import green,cyan
from Bio import SeqIO
import csv

from squirrel.utils.config import *

def get_clade_dict(panel):
    clade_dict = {"*":"unassigned"}
    for record in SeqIO.parse(panel, "fasta"):
        tokens = record.description.split(" ")
        clade = [i.split("=")[1] for i in tokens if i.startswith("clade")][0]            
        clade_dict[record.id] = clade.lower()
    return clade_dict

def parse_paf(paf_file,panel, clade_out, config):
    assigned_clade = {}
    clade_dict = get_clade_dict(panel)
    with open(paf_file, "r") as f:
        for line in f:
            tokens = line.rstrip("\n").split("\t")
            seq_name = tokens[0]
            ref_hit = tokens[5]
            clade = clade_dict[ref_hit]
            # print("####",seq_name, ref_hit, clade)

            assigned_clade[seq_name] = clade
    return assigned_clade

def annotate_fasta_with_clade(assigned_clade, input_fasta, output_fasta):
    with open(output_fasta,"w") as fw:
        for record in SeqIO.parse(input_fasta,"fasta"):
            clade = assigned_clade[record.id]
            fw.write(f">{record.description} clade={clade}\n{record.seq}\n")

def condense_clade_info(assigned_clade):
    clades = set()
    for i in assigned_clade.values():
        if i.startswith("cladeii"):
            clades.add("cladeii")
        elif i.startswith("cladei"):
            clades.add("cladei")
        else:
            clades.add("unassigned")
    return {KEY_ASSIGNED_CLADES: list(clades)}
