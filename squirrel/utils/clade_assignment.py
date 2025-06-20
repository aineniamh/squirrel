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
        clade_dict[record.id] = clade.lower().rstrip("ab")
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


def annotate_fasta_with_clade(clade_info, assigned_clade, input_fasta,output_report, tempdir):
    clades = clade_info[KEY_ASSIGNED_CLADES]
    file_handle_dict = {}

    with open(output_report, "w") as fw:
        writer = csv.DictWriter(fw, lineterminator="\n",fieldnames=["sequence_name","clade"])
        writer.writeheader()
        for clade in clades:
            file_handle_dict[clade] = open(os.path.join(tempdir,f"{clade}.fasta"),"w")

        for record in SeqIO.parse(input_fasta,"fasta"):
            clade = assigned_clade[record.id]
            file_handle_dict[clade].write(f">{record.description}\n{record.seq}\n")

            row = {"sequence_name":record.id,"clade":clade}
            writer.writerow(row)

        for clade in clades:
            file_handle_dict[clade].close()
        
