#!/usr/bin/env python3
import os
import sys
import itertools
import pkg_resources
from Bio import SeqIO

from squirrel.utils.log_colours import green,cyan
from squirrel.utils.config import *
from squirrel import __version__


def setup_config_dict(cwd):
    default_dict = {            

            KEY_INPUT_FASTA:None,
            KEY_OUTFILENAME:None,

            KEY_OUTDIR:cwd,
            KEY_OUTFILE:None,
            KEY_TEMPDIR:None,
            KEY_NO_TEMP:False,

            KEY_TRIM_END:190788,
            KEY_EXTRACT_CDS:False,
            KEY_CONCATENATE:False,

            KEY_VERBOSE: False,
            KEY_THREADS: 1
            }
    return default_dict

def get_snakefile(thisdir,filename):
    snakefile = ""
    # in this case now, the snakefile used should be the name of the analysis mode (i.e. pangolearn, usher or preprocessing)
    snakefile = os.path.join(thisdir, 'scripts',f'{filename}.smk')
    if not os.path.exists(snakefile):
        sys.stderr.write(cyan(f'Error: cannot find Snakefile at {snakefile}. Check installation\n'))
        sys.exit(-1)
    return snakefile


def package_data_check(filename,directory,key,config):
    try:
        package_datafile = os.path.join(directory,filename)
        data = pkg_resources.resource_filename('squirrel', package_datafile)
        config[key] = data
    except:
        sys.stderr.write(colour.cyan(f'Error: Missing package data.')+f'\n\t- {filename}\n')
        sys.exit(-1)

def get_datafiles(config,clade):
    fasta_filename = ""
    if clade.lower() == "cladei":
        fasta_filename = "NC_003310.fasta"
        mask_file = "to_mask.cladei.csv"
    elif clade.lower() == "cladeii":
        fasta_filename = "NC_063383.fasta"
        mask_file = "to_mask.cladeii.csv"
    else:
        sys.stderr.write(colour.cyan(f'Error: invalid clade specified. Please specify one of `cladei` or `cladeii`\n'))
        sys.exit(-1)


    resources = [
            {"key":KEY_REFERENCE_FASTA,
            "directory":"data",
            "filename":fasta_filename},
            {"key":KEY_TO_MASK,
            "directory":"data",
            "filename":mask_file},
            {"key":KEY_GENE_BOUNDARIES,
            "directory":"data",
            "filename":"gene_boundaries.csv"},
            {"key":KEY_GRANTHAM_SCORES,
            "directory":"data",
            "filename":"grantham_score.txt"}
            ]

    for resource in resources:
        package_data_check(resource["filename"],resource["directory"],resource["key"],config)