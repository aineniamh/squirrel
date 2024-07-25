#!/usr/bin/env python3
from squirrel.utils.log_colours import green,cyan

import squirrel.utils.custom_logger as custom_logger
from squirrel.utils.config import *
from squirrel.utils.initialising import *
import squirrel.utils.io_parsing as io
import squirrel.utils.cns_qc as qc
import squirrel.utils.reconstruction_functions as recon

import squirrel.utils.misc as misc
from squirrel import __version__
from . import _program

import os
import sys
import argparse
import snakemake

thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()


def main(sysargs = sys.argv[1:]):
    parser = argparse.ArgumentParser(prog = _program,
    description='squirrel: Some QUIck Rearranging to Resolve Evolutionary Links',
    usage='''squirrel <input> [options]''')

    io_group = parser.add_argument_group('Input-Output options')
    io_group.add_argument('input', nargs="*", help='Input fasta file of sequences to analyse.')
    io_group.add_argument('-o','--outdir', action="store",help="Output directory. Default: current working directory")
    io_group.add_argument('--outfile', action="store",help="Optional output file name. Default: <input>.aln.fasta")
    io_group.add_argument('--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")
    io_group.add_argument("--no-temp",action="store_true",help="Output all intermediate files, for dev purposes.")

    a_group = parser.add_argument_group("Pipeline options")
    a_group.add_argument("--seq-qc",action="store_true",help="Flag potentially problematic SNPs and sequences. Note that this will also run phylo mode, so you will need to specify both outgroup sequences and provide an assembly reference file. Default: don't run QC")
    a_group.add_argument("--assembly-refs",action="store",help="References to check for `calls to reference` against.")
    a_group.add_argument("--no-mask",action="store_true",help="Skip masking of repetitive regions. Default: masks repeat regions")
    a_group.add_argument("--no-itr-mask",action="store_true",help="Skip masking of end ITR. Default: masks ITR")
    a_group.add_argument("--extract-cds",action="store_true",help="Extract coding sequences based on coordinates in the reference")
    a_group.add_argument("--concatenate",action="store_true",help="Concatenate coding sequences for each genome, separated by `NNN`. Default: write out as separate records")
    a_group.add_argument("--clade",action="store",help="Specify whether the alignment is primarily for `cladei` or `cladeii` (will determine reference used for alignment). Default: `cladeii`", default="cladeii")
    a_group.add_argument("-p","--run-phylo",action="store_true",help="Run phylogenetic reconstruction pipeline")
    a_group.add_argument("--outgroups",action="store",help="Specify which MPXV outgroup(s) in the alignment to use in the phylogeny. These will get pruned out from the final tree.")

    m_group = parser.add_argument_group('Misc options')
    m_group.add_argument("-v","--version", action='version', version=f"squirrel {__version__}")
    m_group.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    m_group.add_argument("-t","--threads",action="store",default=1,type=int, help="Number of threads")


    if len(sysargs)<1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    # Initialise config dict
    config = setup_config_dict(cwd)

    get_datafiles(config,args.clade)

    config[KEY_OUTDIR] = io.set_up_outdir(args.outdir,cwd,config[KEY_OUTDIR])
    config[KEY_OUTFILE],config[KEY_CDS_OUTFILE],config[KEY_OUTFILENAME] = io.set_up_outfile(args.outfile,args.input, config[KEY_OUTFILE],config[KEY_OUTDIR])
    io.set_up_tempdir(args.tempdir,args.no_temp,cwd,config[KEY_OUTDIR], config)

    io.pipeline_options(args.no_mask, args.no_itr_mask, args.extract_cds, args.concatenate, config)

    config[KEY_INPUT_FASTA] = io.find_query_file(cwd, config[KEY_TEMPDIR], args.input)
    
    if args.seq_qc:
        assembly_refs = qc.find_assembly_refs(cwd,args.assembly_refs,config)
        args.run_phylo = True
        
        # config[KEY_INPUT_FASTA] = qc.add_refs_to_input(config[KEY_INPUT_FASTA],assembly_refs,config)

    config[KEY_FIG_HEIGHT] = recon.get_fig_height(config[KEY_INPUT_FASTA])
    io.phylo_options(args.run_phylo,args.outgroups,config[KEY_INPUT_FASTA],config)

    snakefile = get_snakefile(thisdir,"msa")

    status = misc.run_snakemake(config,snakefile,args.verbose,config)


    if status:

        if config[KEY_RUN_PHYLO]:
            phylo_snakefile = get_snakefile(thisdir,"reconstruction")
            phylo_stem = ".".join(config[KEY_OUTFILENAME].split(".")[:-1])
            phylo_stem=phylo_stem.split("/")[-1]

            config[KEY_PHYLOGENY] = f"{phylo_stem}.tree"
            config[KEY_PHYLOGENY_SVG] = f"{phylo_stem}.tree.svg"
            config[KEY_OUTGROUP_STRING] = ",".join(config[KEY_OUTGROUPS])
            config[KEY_OUTGROUP_SENTENCE] = " ".join(config[KEY_OUTGROUPS])

            status = misc.run_snakemake(config,phylo_snakefile,args.verbose,config)

            if status:
                print(green("Ancestral reconstruction & phylogenetics complete."))

                if args.seq_qc:

                    qc.check_for_snp_anomalies(config,config[KEY_FIG_HEIGHT])
        else:
            print(green("Alignment complete."))