#!/usr/bin/env python3
from squirrel.utils.log_colours import green,cyan

import squirrel.utils.custom_logger as custom_logger
from squirrel.utils.config import *
from squirrel.utils.initialising import *
import squirrel.utils.io_parsing as io
import squirrel.utils.cns_qc as qc
import squirrel.utils.clade_assignment as ca
import squirrel.utils.reconstruction_functions as recon
from squirrel.utils.make_report import *

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
    io_group.add_argument("--epi2me-outdir",action="store",help="Output directory for report paths in epi2me.")


    a_group = parser.add_argument_group("Alignment options")
    a_group.add_argument("-qc","--seq-qc",action="store_true",help="Flag potentially problematic SNPs and sequences. Default: don't run QC")
    a_group.add_argument("--assembly-refs",action="store",help="References to check for `calls to reference` against.")
    a_group.add_argument("--no-mask",action="store_true",help="Skip masking of repetitive regions. Default: masks repeat regions")
    a_group.add_argument("--no-itr-mask",action="store_true",help="Skip masking of end ITR. Default: masks ITR")
    a_group.add_argument("--additional-mask",action="store",help="Masking additional sites provided as a csv. Needs columns `Maximum` and `Minimum` in 1-base.")
    a_group.add_argument("--sequence-mask",action="store",help="Mask sites in specific sequences in the alignment as a csv, rather than the whole alignment column. Needs `sequence` and `site` (1-based) column.")
    a_group.add_argument("-ex","--exclude",action="store",help="Supply a csv file listing sequences that should be excluded from the analysis.")
    a_group.add_argument("--extract-cds",action="store_true",help="Extract coding sequences based on coordinates in the reference")
    a_group.add_argument("--concatenate",action="store_true",help="Concatenate coding sequences for each genome, separated by `NNN`. Default: write out as separate records")
    a_group.add_argument("--clade",action="store",help="Specify whether the alignment is primarily for `cladei` or `cladeii` (can also specify a or b, e.g. `cladeia`, `cladeiib`). `split` will separate input sequences into separate cladei and cladeii alignments. This will determine reference used for alignment, mask file and background set used if `--include-background` flag used in conjunction with the `--run-phylo` option. Default: `cladeii`")

    p_group = parser.add_argument_group("Phylo options")
    p_group.add_argument("-p","--run-phylo",action="store_true",help="Run phylogenetics pipeline")
    p_group.add_argument("-a","--run-apobec3-phylo",action="store_true",help="Run phylogenetics & APOBEC3-mutation reconstruction pipeline")
    p_group.add_argument("--outgroups",action="store",help="Specify which MPXV outgroup(s) in the alignment to use in the phylogeny. These will get pruned out from the final tree.")
    p_group.add_argument("-bg","--include-background",action="store_true",help="Include a default background set of sequences for the phylogenetics pipeline. The set will be determined by the `--clade` specified.")
    p_group.add_argument("-bf","--background-file",action="store",help="Include this additional FASTA file as background to the phylogenetics.")
    p_group.add_argument("-bm","--binary-partition-mask",action="store_true",help="Calculate and write binary partition mask")
    p_group.add_argument("--bm-separate-dimers",action="store_true",help="Write partition mask with 0 for non-apo, 1 for GA and 2 for TC target sites")

    pf_group = parser.add_argument_group("Tree figure options")
    pf_group.add_argument("-tfig","--tree-figure-only",action="store_true",help="Re-render tree figure custom height and width arguments. Requires: tree file, branch reconstruction file, height, width.")
    pf_group.add_argument('-tf',"--tree-file",action="store",help="Tree for re-rendering the figure.")
    pf_group.add_argument('-brf',"--branch-reconstruction-file",action="store",help="Reconstruction file for re-rendering the figure.")
    pf_group.add_argument("--fig-height",action="store",help="Overwrite tree figure default height.",type=int)
    pf_group.add_argument("--fig-width",action="store",help="Overwrite tree figure default width.",type=int)
    pf_group.add_argument("--point-style",action="store",help="Shape of points for apobec3 reconstruction figure. Options: circle, square. Default: circle")
    pf_group.add_argument("--point-justify",action="store",help="Justification of points for apobec3 reconstruction figure. Options: left, right. Default: left")
    pf_group.add_argument("--interactive-tree",action="store_true",help="Create a separate interactive file built from R functions")

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
    
    io.set_up_clade(args.clade,config)
    get_global_datafiles(config)
    config["version"] = __version__
    
    io.set_up_threads(args.threads,config)
    config[KEY_OUTDIR] = io.set_up_outdir(args.outdir,cwd,config[KEY_OUTDIR])
    io.set_up_epi2me_outdir(args.epi2me_outdir, config)
    io.parse_tf_options(args.tree_figure_only,args.tree_file,args.branch_reconstruction_file,args.fig_width,args.fig_height,args.point_style,args.point_justify,cwd,config)
    if args.tree_figure_only:
        new_tree = f"{config[KEY_TREE]}.rerender"
        outfile = os.path.join(config[KEY_OUTDIR],"")
        recon.make_reconstruction_tree_figure_w_labels(new_tree,config[KEY_BRANCH_RECONSTRUCTION],config[KEY_TREE],config[KEY_POINT_STYLE],config[KEY_POINT_JUSTIFY],config[KEY_FIG_WIDTH],config[KEY_FIG_HEIGHT])
        print(green("Success! New tree figure written."))
        sys.exit(0)
    
    config[KEY_OUTFILE_STEM],config[KEY_OUTDIR] = io.set_up_outfile_stem(args.outfile,cwd,args.input, config[KEY_OUTDIR])
    io.set_up_tempdir(args.tempdir,args.no_temp,cwd,config[KEY_OUTDIR], config)#outfile_arg, cwd, query_arg, outfile, outdir

    io.pipeline_options(args.no_mask, args.no_itr_mask, args.additional_mask,args.sequence_mask, args.extract_cds, args.concatenate,config[KEY_CLADE],cwd, config)

    config[KEY_INPUT_FASTA] = io.find_query_file(cwd, config[KEY_TEMPDIR], args.input)

    if args.background_file:
        config[KEY_INPUT_FASTA] = io.find_background_file(cwd,config[KEY_INPUT_FASTA],args.background_file,config)

    if args.exclude:
        config[KEY_INPUT_FASTA] = io.find_exclude_file(cwd,config[KEY_INPUT_FASTA],args.exclude,config)

    if config[KEY_CLADE] == "split":

        clade_snakefile = get_script(thisdir,"clade.smk")
        status = misc.run_snakemake(config,clade_snakefile,args.verbose,config)

        
        clade_config = misc.load_yaml(os.path.join(config[KEY_TEMPDIR],"clades.yaml"))
        config[KEY_ASSIGNED_CLADES] = clade_config[KEY_ASSIGNED_CLADES]
        print(green("Clades assigned:"))
        
        for clade in config[KEY_ASSIGNED_CLADES]:
            print(f"- {clade}")

        if "unassigned" in config[KEY_ASSIGNED_CLADES]:
            if len(config[KEY_ASSIGNED_CLADES])==1:
                sys.stderr.write(cyan(
                f'Error: no sequences had a clade successfully assigned.\nThis may be due to incomplete/ low coverage genome sequences.\n'))
                sys.exit(-1)
            else:
                print(cyan("Warning: some sequences were not successfully assigned a clade. \nThis may be due to incomplete/ low coverage genome sequences.\nSee assignment_report.csv for more information."))
                config[KEY_ASSIGNED_CLADES] = [i for i in config[KEY_ASSIGNED_CLADES] if not i=="unassigned"]
    else:
        config[KEY_ASSIGNED_CLADES] = [config[KEY_CLADE]]
    
    
    for clade in config[KEY_ASSIGNED_CLADES]:
        
        if len(config[KEY_ASSIGNED_CLADES]) >1:
            config[KEY_APPEND_CLADE_STR] = f".{clade}"

            print(f"Running analysis for {clade} sequences.")


        config[KEY_CLADE] = clade
        
        config[KEY_OUTFILENAME] = f"{config[KEY_OUTFILE_STEM]}{config[KEY_APPEND_CLADE_STR]}.aln.fasta"
        config[KEY_CDS_OUTFILENAME] = f"{config[KEY_OUTFILE_STEM]}{config[KEY_APPEND_CLADE_STR]}.cds.aln.fasta"


        config[KEY_OUTFILE] = os.path.join(config[KEY_OUTDIR],config[KEY_OUTFILENAME])
        config[KEY_CDS_OUTFILE] = os.path.join(config[KEY_OUTDIR],config[KEY_CDS_OUTFILENAME])
        #filter the input file

        get_datafiles(config)

        if args.clade == "split":
            config[KEY_INPUT_FASTA] = os.path.join(config[KEY_TEMPDIR],f"{clade}.fasta")
            
        if args.seq_qc:
            print(green("QC mode activated. Squirrel will flag:"))
            print("- Clumps of unique SNPs\n- SNPs adjacent to Ns\n- Sequences with high N content")
            config[KEY_SEQ_QC] = True
            config[KEY_EXCLUDE_FILE] = os.path.join(config[KEY_OUTDIR],f"{VALUE_EXCLUDE_FILE_STEM}{config[KEY_APPEND_CLADE_STR]}.csv")
            qc.check_flag_N_content(config[KEY_INPUT_FASTA],config[KEY_EXCLUDE_FILE],config)
        
        assembly_refs = []
        if args.seq_qc and args.run_apobec3_phylo:
            print("- Reversions to reference\n- Convergent mutations")
            assembly_refs = qc.find_assembly_refs(cwd,args.assembly_refs,config)
            # args.run_phylo = True

        # config[KEY_FIG_HEIGHT] = recon.get_fig_height(config[KEY_INPUT_FASTA])

        config[KEY_INPUT_FASTA] = io.phylo_options(args.run_phylo,args.run_apobec3_phylo,args.outgroups,args.include_background,args.binary_partition_mask,config[KEY_INPUT_FASTA],clade,config)
        snakefile = get_script(thisdir,"msa.smk")

        status = misc.run_snakemake(config,snakefile,args.verbose,config)

        if status:

            if config[KEY_RUN_PHYLO]:
                phylo_snakefile = get_script(thisdir,"phylo.smk")
                config[KEY_PHYLOGENY] = f"{config[KEY_OUTFILE_STEM]}{config[KEY_APPEND_CLADE_STR]}.tree"
                
                config[KEY_OUTGROUP_STRING] = ",".join(config[KEY_OUTGROUPS])
                config[KEY_OUTGROUP_SENTENCE] = " ".join(config[KEY_OUTGROUPS])

                config[KEY_PHYLOGENY_INTERACTIVE] = f"{config[KEY_OUTFILE_STEM]}.tree.interactive.html"
                
                if args.interactive_tree:
                    config[KEY_INTERACTIVE_SCRIPT] = get_script(thisdir, "interactive_tree.R")
                else: 
                    config[KEY_INTERACTIVE_SCRIPT] = ''

                if config[KEY_RUN_APOBEC3_PHYLO]:
                    config[KEY_PHYLOGENY_SVG] = f"{config[KEY_OUTFILE_STEM]}{config[KEY_APPEND_CLADE_STR]}.tree.svg"
                    phylo_snakefile = get_script(thisdir,"reconstruction.smk")

                status = misc.run_snakemake(config,phylo_snakefile,args.verbose,config)

                if status:
                    if config[KEY_RUN_APOBEC3_PHYLO]:
                        if args.binary_partition_mask:
                            outfile = os.path.join(config[KEY_OUTDIR],f"{config[KEY_OUTFILE_STEM]}{config[KEY_APPEND_CLADE_STR]}.binary_partition_mask.csv")
                            branch_reconstruction = os.path.join(config[KEY_OUTDIR],f"{config[KEY_OUTFILE_STEM]}{config[KEY_APPEND_CLADE_STR]}.tree.branch_snps.reconstruction.csv")
                            recon.find_binary_partition_mask(branch_reconstruction,args.bm_separate_dimers,config[KEY_REFERENCE_FASTA],outfile)
                            print(green(f"Binary partition mask string written to: "),outfile)
                        print(green("Ancestral reconstruction & phylogenetics complete."))
                    else:
                        print(green("Phylogenetics complete."))

            if args.seq_qc:
                mask_file = os.path.join(config[KEY_OUTDIR],f"{config[KEY_OUTFILE_STEM]}{config[KEY_APPEND_CLADE_STR]}.suggested_mask.csv")

                qc.check_for_snp_anomalies(assembly_refs,mask_file,config,config[KEY_FIG_HEIGHT])
                print(green("Flagged mutations writted to:"), f"{mask_file}")
            else:
                print(green("Alignment complete."))
                mask_file = ""
            
            # get the inputs for making the overall report

            
        else:

            sys.exit(-1)
            
    report =os.path.join(config[KEY_OUTDIR],f"{config[KEY_OUTFILE_STEM]}.report.html")
    make_output_report(report,config)

