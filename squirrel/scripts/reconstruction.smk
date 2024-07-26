import os
from squirrel.utils.config import *
from Bio import SeqIO
from Bio.Seq import Seq
import csv
from squirrel.utils.log_colours import green,cyan
import squirrel.utils.reconstruction_functions as recon

rule all:
    input:
        os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY]),
        os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY_SVG])


rule iqtree:
    input:
        aln = config[KEY_OUTFILE]
    params:
        outgroup = config[KEY_OUTGROUP_STRING],
        threads = config[KEY_PHYLO_THREADS]
    output:
        tree = f"{config[KEY_OUTFILE]}.treefile"
    shell:
        """
        iqtree  -s {input.aln:q} \
                -m HKY \
                --ufboot 1000 \
                -czb \
                -nt {params.threads} \
                -blmin  0.0000000001 \
                -redo \
                -asr \
                -o '{params.outgroup}' 
        """

rule prune_outgroup:
    input:
        tree = rules.iqtree.output.tree
    params:
        treefile = f"{config[KEY_OUTFILENAME]}.treefile",
        outdir = config[KEY_OUTDIR],
        outtree = config[KEY_PHYLOGENY],
        outgroup = config[KEY_OUTGROUP_SENTENCE]
    output:
        tree = os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY])
    shell:
        """
        cd  {params.outdir:q} &&
        jclusterfunk prune  -i "{params.treefile}" \
                            -o "{params.outtree}" \
                            -t '{params.outgroup}'
        """


rule reconstruction_analysis:
    input:
        tree = rules.prune_outgroup.output.tree,
        alignment = config[KEY_OUTFILE]
    params:
        outdir = config[KEY_OUTDIR]
    output:
        tree = os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY_SVG])
    run:
        directory = params.outdir
        width= 25

        height = recon.get_fig_height(input.alignment)

        recon.run_full_analysis(directory, input.alignment, input.tree,config,width,height)
