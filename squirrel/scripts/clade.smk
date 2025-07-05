import os
from squirrel.utils.config import *
from Bio import SeqIO
from Bio.Seq import Seq
import csv
from squirrel.utils.log_colours import green,cyan
import collections
import yaml

import squirrel.utils.clade_assignment as ca
import squirrel.pangolearn.pangolearn as pangolearn

rule all:
        input:
            os.path.join(config[KEY_TEMPDIR],"clades.yaml"),
            os.path.join(config[KEY_OUTDIR],"assignment_report.csv")

rule align_to_c2_reference:
    input:
        fasta = config[KEY_INPUT_FASTA],
        reference = config[KEY_ASSIGNMENT_REFERENCE]
    params:
        trim_start = 0,
        trim_end = config[KEY_TRIM_END],
        sam = os.path.join(config[KEY_TEMPDIR],"assignment_map.sam")
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"assignment_msa.aln.fasta")
    log:
        os.path.join(config[KEY_TEMPDIR], "logs/assignment_minimap2_sam.log")
    run:
        # the first line of this streams through the fasta and replaces '-' in sequences with empty strings
        # this could be replaced by a python script later
        #  {{ gsub(" ","_",$0); }} {{ gsub(",","_",$0); }}
        shell_command =  """ | awk '{{ if ($0 !~ /^>/) {{ gsub("-", "",$0); }} print $0; }}'   | \
            awk '{{ {{ gsub(" ", "_",$0); }} {{ gsub(",", "_",$0); }} print $0; }}'  | \
            minimap2 -a -x asm20 -rmq=no --junc-bonus=0 --for-only --sam-hit-only --secondary=no --score-N=0  -t  {workflow.cores} {input.reference:q} - -o {params.sam:q} &> {log:q}
            gofasta sam toMultiAlign \
                -s {params.sam:q} \
                -t {workflow.cores} \
                --reference {input.reference:q} \
                --trimstart {params.trim_start} \
                --trimend {params.trim_end} \
                --trim \
                --pad -o '{output.fasta}' &> {log:q}
            """

        shell_command = "cat \"{input.fasta}\"" + shell_command
        shell(shell_command)

rule pangolearn:
    input:
        fasta = rules.align_to_c2_reference.output.fasta,
        model = config[KEY_PLEARN_MODEL],
        header = config[KEY_PLEARN_HEADER],
        reference = config[KEY_ASSIGNMENT_REFERENCE]
    params:
        tempdir = config[KEY_TEMPDIR]
    output:
        csv = os.path.join(config[KEY_TEMPDIR],"assignment_report.csv"),
        clade_config = os.path.join(config[KEY_TEMPDIR],"clades.yaml")
    run:
        pangolearn.assign_lineage(input.header,input.model,input.reference,input.fasta,output.csv)

        clade_info,assigned_clade = ca.condense_clade_info(output.csv)

        ca.annotate_fasta_with_clade(clade_info, assigned_clade, input.fasta, params.tempdir)

        with open(output.clade_config, 'w') as f:
            yaml.dump(clade_info, f)
