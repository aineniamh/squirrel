import os
from squirrel.utils.config import *
from Bio import SeqIO
from Bio.Seq import Seq
import csv
from squirrel.utils.log_colours import green,cyan
import squirrel.utils.reconstruction_functions as recon

if config[KEY_INTERACTIVE_SCRIPT]:
    rule all:
        input:
            os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY]),
            os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY_SVG]),
            os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY_INTERACTIVE]),
            os.path.join(config[KEY_OUTDIR],f"{config[KEY_PHYLOGENY]}.annotated.tree")
else:
    rule all:
        input:
            os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY]),
            os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY_SVG]),
            os.path.join(config[KEY_OUTDIR],f"{config[KEY_PHYLOGENY]}.annotated.tree")

rule iqtree:
    input:
        aln = config[KEY_OUTFILE]
    params:
        outgroup = config[KEY_OUTGROUP_STRING],
        threads = config[KEY_PHYLO_THREADS]
    output:
        temp_aln = os.path.join(config[KEY_TEMPDIR],f"iqtree.fasta"),
        tree = os.path.join(config[KEY_TEMPDIR],f"iqtree.fasta.treefile"),
        state_file = os.path.join(config[KEY_OUTDIR],f"{config[KEY_PHYLOGENY]}.state")
    shell:
        """
        cp {input.aln:q} {output.temp_aln:q} && 
        iqtree  -s {output.temp_aln:q} \
                -m HKY \
                -czb \
                -nt {params.threads} \
                -blmin  0.0000000001 \
                -redo \
                -asr \
                -o '{params.outgroup}' &&
        cp '{output.temp_aln}.state' {output.state_file:q}
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

rule reconstruction_analysis:
    input:
        tree = rules.prune_outgroup.output.tree,
        state_file = rules.iqtree.output.state_file,
        alignment = config[KEY_OUTFILE]
    params:
        outdir = config[KEY_OUTDIR]
    output:
        tree = os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY_SVG]),
        state_differences = os.path.join(config[KEY_OUTDIR],f"{config[KEY_PHYLOGENY]}.state_differences.csv"),
        amino_acids_out = os.path.join(config[KEY_OUTDIR],f"{config[KEY_PHYLOGENY]}.amino_acid.reconstruction.csv"),
        branch_snps_out = os.path.join(config[KEY_OUTDIR],f"{config[KEY_PHYLOGENY]}.branch_snps.reconstruction.csv"),
    run:
        directory = params.outdir
        point_style = config[KEY_POINT_STYLE]
        point_justify = config[KEY_POINT_JUSTIFY]
        # width= config[KEY_FIG_HEIGHT]
        # height = recon.get_fig_height(input.alignment)

        recon.run_full_analysis(directory, input.alignment, input.tree,input.state_file,config,point_style,point_justify,config[KEY_FIG_WIDTH],config[KEY_FIG_HEIGHT])

rule a3_annotation:
    input:
        branch_snps = rules.reconstruction_analysis.output.branch_snps_out,
        treefile = rules.prune_outgroup.output.tree
    output:
        outtree = os.path.join(config[KEY_OUTDIR],f"{config[KEY_PHYLOGENY]}.annotated.tree")
    run:
        recon.annotate_tree(output.outtree,input.branch_snps,input.treefile)

rule interactive_tree:
    input:
        tree = rules.prune_outgroup.output.tree,
        mutations = rules.reconstruction_analysis.output.amino_acids_out
    params:
        script_path = config[KEY_INTERACTIVE_SCRIPT]
    output:
        figure = os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY_INTERACTIVE])
    shell:
        """   
        if [ {params.script_path} ]; then
            Rscript {params.script_path} {input.tree} {input.mutations} {output.figure}
        fi
        """
