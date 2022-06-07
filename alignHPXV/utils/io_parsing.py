#!/usr/bin/env python3
import sys
import os
from alignHPXV.utils.log_colours import green,cyan
import select
from Bio import SeqIO

import tempfile
import shutil

from alignHPXV.utils.config import *


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
    else:
        out_str = ".".join(query_arg[0].split(".")[:-1]) + ".aln.fasta"
        outfile = os.path.join(outdir, out_str)
    return outfile


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