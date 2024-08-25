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
        outgroup = config[KEY_OUTGROUP_STRING],
        threads = config[KEY_PHYLO_THREADS]
    output:
        temp_aln = os.path.join(config[KEY_TEMPDIR],f"iqtree.fasta"),
        tree = os.path.join(config[KEY_TEMPDIR],f"iqtree.fasta.treefile")
    shell:
        """
        cp {input.aln:q} {output.temp_aln:q} && 
        iqtree  -s {output.temp_aln:q} \
                -m HKY \
                -czb \
                -nt {params.threads} \
                -blmin  0.0000000001 \
                -redo \
                -o '{params.outgroup}' 
        """

rule prune_outgroup:
    input:
        tree = rules.iqtree.output.tree
    params:
        treefile = f"iqtree.fasta.treefile",
        temp_outtree = "iqtree.pruned.tree",
        outdir = config[KEY_OUTDIR],
        tempdir = config[KEY_TEMPDIR],
        outtree = config[KEY_PHYLOGENY],
        outgroup = config[KEY_OUTGROUP_SENTENCE]
    output:
        tree = os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY])
    shell:
        """
        cd  {params.tempdir:q} &&
        jclusterfunk prune  -i "{params.treefile}" \
                            -o "{params.temp_outtree}" \
                            -t '{params.outgroup}' &&
        cp {params.temp_outtree:q} {output.tree:q}
        """
