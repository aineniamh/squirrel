import os
from squirrel.utils.config import *
from Bio import SeqIO
from Bio.Seq import Seq
import csv
from squirrel.utils.log_colours import green,cyan
import collections

from squirrel.utils.clade_assignment import parse_paf

rule all:
        input:
            os.path.join(config[KEY_TEMPDIR],"mapped.paf")

rule map_seqs:
    input:
        ref = config[KEY_REFERENCE_PANEL],
        fasta = config[KEY_INPUT_FASTA]
    threads: workflow.cores
    log: os.path.join(config[KEY_TEMPDIR],"minimap2_clades.log")
    output:
        paf = os.path.join(config[KEY_TEMPDIR],"mapped.paf")
    shell:
        """
        minimap2 -t {threads} -x asm20 --secondary=no --paf-no-hit \
        {input.ref:q} \
        {input.fastq:q} -o {output:q} &> {log:q}
        """

rule determine_clades:
    input:
        paf = rules.map_seqs.output.paf
    output:
        clade_config = os.path.join(config[KEY_TEMPDIR],"clades.yaml"),
        fasta = os.path.join(config[KEY_TEMPDIR],"input.clade_annotated.fasta")
    run:
        parse_paf()
    
