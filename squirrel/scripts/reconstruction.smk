import os
from squirrel.utils.config import *
from Bio import SeqIO
from Bio.Seq import Seq
import csv
from squirrel.utils.log_colours import green,cyan
import squirrel.utils.reconstruction_functions as recon

rule all:
    input:
        config[KEY_PHYLOGENY],
        config[KEY_PHYLOGENY_SVG]


rule iqtree:
    input:
        aln = config[KEY_OUTFILE]
    params:
        outgroup = config[KEY_OUTGROUP_STRING]
    output:
        tree = f"{config[KEY_OUTFILE]}.treefile"
    shell:
        """
        iqtree  -s {input.aln:q} \
                -m HKY \
                -czb \
                -blmin  0.0000000001 \
                -redo \
                -asr \
                -o '{params.outgroup}' 
        """

rule prune_outgroup:
    input:
        tree = rules.iqtree.output.tree
    params:
        outgroup = config[KEY_OUTGROUP_SENTENCE]
    output:
        tree = config[KEY_PHYLOGENY]
    shell:
        """
        jclusterfunk prune  -i {input.tree:q} \
                            -o {output.tree:q} \
                            -t '{params.outgroup}' 
        """


rule reconstruction_analysis:
    input:
        tree = rules.prune_outgroup.output.tree,
        alignment = config[KEY_OUTFILE]
    params:
        outdir = config[KEY_OUTDIR]
    output:
        tree = config[KEY_PHYLOGENY_SVG]
    run:
        directory = params.outdir
        width= 25

        height = recon.get_fig_height(input.alignment)

        recon.run_full_analysis(directory, input.alignment, input.tree,config,width,height)


