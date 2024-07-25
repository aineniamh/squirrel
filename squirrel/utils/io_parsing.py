#!/usr/bin/env python3
import sys
import os
from squirrel.utils.log_colours import green,cyan
import select
from Bio import SeqIO

import tempfile
import shutil

from squirrel.utils.config import *


def set_up_outdir(outdir_arg,cwd,outdir):
    if outdir_arg:
        outdir = os.path.join(cwd, outdir_arg)
        if not os.path.exists(outdir):
            try:
                os.mkdir(outdir)
            except:
                sys.stderr.write(cyan(f'Error: cannot create directory:') + f"{outdir}")
                sys.exit(-1)
    return outdir


def set_up_outfile(outfile_arg,query_arg, outfile, outdir):
    if outfile_arg:
        outfile = os.path.join(outdir, outfile_arg)
        cds_outstr = ".".join(query_arg[0].split(".")[:-1]) + ".aln.cds.fasta"
        cds_outfile = os.path.join(outdir, cds_outstr)
        outfilename=outfile_arg
    elif query_arg:
        out_str = ".".join(query_arg[0].split(".")[:-1]) + ".aln.fasta"
        outfile = os.path.join(outdir, out_str)
        cds_outstr = ".".join(query_arg[0].split(".")[:-1]) + ".aln.cds.fasta"
        cds_outfile = os.path.join(outdir, cds_outstr)
        outfilename=out_str
    else:
        out_str = "sequences.aln.fasta"
        outfile = os.path.join(outdir, out_str)
        cds_outstr = "sequences.aln.cds.fasta"
        cds_outfile = os.path.join(outdir, cds_outstr)
        outfilename=out_str

    return outfile,cds_outfile,outfilename


def set_up_tempdir(tempdir_arg,no_temp_arg,cwd,outdir,config):

    if no_temp_arg:
        tempdir = outdir
        config[KEY_TEMPDIR] = tempdir
        print(green(f"\n--no-temp: ") + f"all intermediate files will be written to {outdir}\n")
    elif tempdir_arg:
        to_be_dir = os.path.join(cwd, tempdir_arg)
        try:
            if not os.path.exists(to_be_dir):
                os.mkdir(to_be_dir)
        except:
            sys.stderr.write(cyan(f'Error: cannot create temp directory {to_be_dir}.\n'))
            sys.exit(-1)
        tempdir = tempfile.mkdtemp(dir=to_be_dir)
        config[KEY_TEMPDIR] = tempdir
    else:
        tempdir = tempfile.mkdtemp()
        config[KEY_TEMPDIR] = tempdir
        try:
            if not os.path.exists(tempdir):
                os.mkdir(tempdir)
        except:
            sys.stderr.write(cyan(f'Error: cannot create temp directory {tempdir}.\n'))
            sys.exit(-1)
        
        try:
            with open(os.path.join(tempdir, "test.txt"),"w") as fw:
                fw.write("Test")
        except:
            sys.stderr.write(cyan(f'Error: cannot write to temp directory {tempdir}.\n'))
            sys.exit(-1)

def cleanup(no_temp,tempdir):
    if not no_temp:
        shutil.rmtree(tempdir)

def find_query_file(cwd, tempdir, query_arg):
    if len(query_arg) > 1:
        print(cyan(f"Error: Too many query (input) fasta files supplied: {query_arg}\nPlease supply one only."))
        sys.exit(-1)

    # find the query fasta
    try:
        if not os.path.exists(os.path.join(cwd, query_arg[0])):
            if select.select([sys.stdin,],[],[],0.0)[0]:
                query = os.path.join(tempdir, "stdin_query.fasta")
                with open(query,"w") as fw:
                    for l in sys.stdin:
                        l= l.rstrip("\n")
                        fw.write(l + '\n')
                
                print(green("Query:\t") + "reading from stdin.")
            elif not select.select([sys.stdin,],[],[],0.0)[0]:
                tried_path = os.path.join(cwd, query_arg[0])
                if tried_path.endswith("-"):
                    sys.stderr.write(cyan(
                        f'Error: cannot find query (input) fasta file using stdin.\n'))
                    sys.exit(-1)
                else:
                    sys.stderr.write(cyan(f'Error: cannot find query (input) fasta file at:') + f'{tried_path}\n')
                    sys.exit(-1)
        else:
            query = os.path.join(cwd, query_arg[0])
            print(green(f"Query file:\t") + f"{query}")
    except IndexError:
        sys.stderr.write(cyan(
            f'Error: input query fasta could not be detected from a filepath or through stdin.\n'))
        sys.exit(-1)

    return query


def find_additional_mask_file(cwd,additional_mask,config):

    path_to_try = os.path.join(cwd,additional_mask)
    try:
        with open(path_to_try,"r") as f:
            reader = csv.DictReader(f)
            header = reader.fieldnames
            for i in ["Maximum","Minimum"]:
                if i not in header:
                    sys.stderr.write(cyan(f'Error: additional mask file must contain columns `Maximum` and `Minimum`.\n'))
                    sys.exit(-1)
            for row in reader:
                try:
                    mx = int(row["Maximum"])
                    mn = int(row["Minimum"])
                except:
                    sys.stderr.write(cyan(f'Error: `Maximum` and `Minimum` columns must contain numeric values.\n'))
                    sys.exit(-1)
    except:
        sys.stderr.write(cyan(f'Error: cannot find additional mask file at: ') + f'{path_to_try}\n' + cyan('Please check file path and try again.\n'))
        sys.exit(-1)

    return path_to_try



def pipeline_options(no_mask, no_itr_mask, additional_mask,extract_cds,concatenate,cwd, config):
    config[KEY_NO_MASK] = no_mask
    if no_itr_mask:
        config[KEY_TRIM_END] = 197209
    
    if additional_mask:
        config[KEY_ADDITIONAL_MASK] = find_additional_mask_file(cwd,additional_mask,config)

    config[KEY_EXTRACT_CDS] = extract_cds
    config[KEY_CONCATENATE] = concatenate

def phylo_options(run_phylo,outgroups,alignment,config):
    config[KEY_RUN_PHYLO] = run_phylo

    if config[KEY_RUN_PHYLO]:

        if not outgroups:
            sys.stderr.write(cyan(
                        f'Error: must supply outgroup(s) for phylogenetics module.\n'))
            sys.exit(-1)
        if not type(outgroups) == list:
            outgroups = outgroups.split(",")
        config[KEY_OUTGROUPS] = outgroups
        
        seqs = SeqIO.index(alignment,"fasta")
        not_in = set()
        for outgroup in outgroups:
            if outgroup not in seqs:
                not_in.add(outgroup)

        if not_in:
            sys.stderr.write(cyan(
                        f'Error: outgroup(s) not found in input sequence file.\n'))
            for seq in not_in:
                sys.stderr.write(cyan(f"- {seq}\n"))
            sys.exit(-1)
