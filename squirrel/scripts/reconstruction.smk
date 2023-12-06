import os
from squirrel.utils.config import *
from Bio import SeqIO
from Bio.Seq import Seq
import csv
from squirrel.utils.log_colours import green,cyan


rule all:
    input:
        os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY])


rule iqtree:
    input:
        aln = config[KEY_OUTFILE]
    params:
        outgroup = config[KEY_OUTGROUP_STRING]
    output:
        tree = os.path.join(f"{config[KEY_OUTFILE]}.treefile")
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
        tree = os.path.join(config[KEY_PHYLOGENY])
    shell:
        """
        jclusterfunk prune  -i {input.tree:q} \
                            -o {output.tree:q} \
                            -t '{params.outgroup}' \
                            -f newick 
        """

