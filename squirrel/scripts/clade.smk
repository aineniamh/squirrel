import os
from squirrel.utils.config import *
from Bio import SeqIO
from Bio.Seq import Seq
import csv
from squirrel.utils.log_colours import green,cyan
import collections
import yaml

import squirrel.utils.clade_assignment as ca

rule all:
        input:
            os.path.join(config[KEY_TEMPDIR],"clades.yaml"),
            os.path.join(config[KEY_OUTDIR],"assignment_report.csv")

rule map_seqs:
    input:
        ref = config[KEY_REFERENCE_PANEL],
        fasta = config[KEY_INPUT_FASTA]
    threads: workflow.cores
    log: os.path.join(config[KEY_TEMPDIR],"minimap2_clades.log")
    output:
        os.path.join(config[KEY_TEMPDIR],"mapped.paf")
    shell:
        """
        minimap2 -t {threads} -x asm20 --secondary=no --paf-no-hit \
        {input.ref:q} \
        {input.fasta:q} -o {output:q} &> {log:q}
        """

rule determine_clades:
    input:
        paf = rules.map_seqs.output[0],
        panel = config[KEY_REFERENCE_PANEL],
        fasta = config[KEY_INPUT_FASTA]
    params:
        tempdir = config[KEY_TEMPDIR]
    output:
        clade_config = os.path.join(config[KEY_TEMPDIR],"clades.yaml"),
        csv = os.path.join(config[KEY_OUTDIR],"assignment_report.csv")
    run:
        assigned_clade = ca.parse_paf(input.paf,input.panel, output.clade_config, config)

        clade_info = ca.condense_clade_info(assigned_clade)

        ca.annotate_fasta_with_clade(clade_info, assigned_clade, input.fasta, output.csv, params.tempdir)

        with open(output.clade_config, 'w') as f:
            yaml.dump(clade_info, f)
