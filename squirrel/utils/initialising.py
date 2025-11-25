#!/usr/bin/env python3
import os
import sys
import itertools
import importlib.resources
from Bio import SeqIO

from squirrel.utils.log_colours import green,cyan
from squirrel.utils.config import *
from squirrel import __version__


def setup_config_dict(cwd):
    default_dict = {            

            KEY_INPUT_FASTA:None,
            KEY_OUTFILENAME:None,

            KEY_CLADE:"cladeii",
            KEY_OUTGROUPS:[],
            KEY_ASSIGNED_CLADES:[],
            KEY_APPEND_CLADE_STR:"",
            KEY_SPLIT_CLADE:False,

            KEY_OUTDIR:cwd,
            KEY_EPI2ME_OUTDIR:None,
            KEY_OUTFILE:None,
            KEY_OUTFILE_STEM:None,
            KEY_TEMPDIR:None,
            KEY_NO_TEMP:False,

            KEY_ASSEMBLY_REFERENCES:[],
            
            KEY_TREE:None,
            KEY_BRANCH_RECONSTRUCTION:None,
            KEY_FIG_HEIGHT:None,
            KEY_FIG_WIDTH:None,
            KEY_POINT_STYLE:"circle",
            KEY_POINT_JUSTIFY:"left",

            KEY_TRIM_END:VALUE_TRIM_END,
            KEY_EXTRACT_CDS:False,
            KEY_CONCATENATE:False,
            KEY_ADDITIONAL_MASK:None,
            KEY_SEQUENCE_MASK:None,
            KEY_SEQ_QC:False,
            KEY_RUN_PHYLO:False,
            KEY_RUN_APOBEC3_PHYLO:False,

            KEY_VERBOSE: False,
            KEY_THREADS: 1,
            KEY_PHYLO_THREADS: "AUTO",
            KEY_INCLUDE_BACKGROUND:False

            }
    return default_dict

# Necessary for scripts other than the python utilities (which should integrate more seamlessly through import)
def get_script(thisdir, filename):
    script_path = ""

    script_path = os.path.join(thisdir, 'scripts', filename)
    if not os.path.exists(script_path):
        sys.stderr.write(cyan(f'Error: cannot find file at {script_path}. Check installation\n'))
        sys.exit(-1)
    return script_path

def package_data_check(filename,directory,key,config):
    try:
        package_datafile = os.path.join(directory,filename)
        resource_path = importlib.resources.files('squirrel').joinpath(package_datafile)
        data = os.path.join(resource_path)
        config[key] = data
    except:
        sys.stderr.write(cyan(f'Error: Missing package data.')+f'\n\t- {filename}\n')
        sys.exit(-1)

def get_global_datafiles(config):

    resources = [
            {"key":KEY_REPORT_TEMPLATE,
            "directory":"data",
            "filename":"report.mako"},
            {"key":KEY_BACKGROUND_FASTA,
            "directory":"data",
            "filename":"background.fasta"},
            {"key":KEY_REFERENCE_PANEL,
            "directory":"data",
            "filename":"background.panel.fasta"},
            {"key":KEY_GRANTHAM_SCORES,
            "directory":"data",
            "filename":"grantham_score.txt"},
            {"key":KEY_ASSEMBLY_REFERENCES,
            "directory":"data",
            "filename":"ref_seq.fasta"},
            {"key":KEY_ASSIGNMENT_REFERENCE,
            "directory":"data",
            "filename":"NC_063383.fasta"},
            {"key":KEY_PLEARN_MODEL,
            "directory":"data",
            "filename":"randomForest_v2.joblib"},
            {"key":KEY_PLEARN_HEADER,
            "directory":"data",
            "filename":"randomForestHeaders_v2.joblib"}
            ]

    for resource in resources:
        package_data_check(resource["filename"],resource["directory"],resource["key"],config)

def get_datafiles(config):
    clade = config[KEY_CLADE].lower()
    config[KEY_CLADE] = clade
    fasta_filename = ""
    gene_boundaries_file = ""
    mask_file = ""
    if clade in ["cladei","cladeia","cladeib"]:
        fasta_filename = "NC_003310.fasta"
        mask_file = "to_mask.cladei.csv"
        gene_boundaries_file = "gene_boundaries.cladei.csv"
    elif clade in ["cladeii","cladeiia","cladeiib"]:
        fasta_filename = "NC_063383.fasta"
        mask_file = "to_mask.cladeii.csv"
        gene_boundaries_file = "gene_boundaries.cladeii.csv"
    elif clade=="variola":
        config[KEY_TRIM_END] = 184722 #184546
        fasta_filename = "NC_001611.fasta"
        mask_file = "to_mask.NC_001611.csv"
        gene_boundaries_file = "gene_boundaries.NC_001611.csv"
    elif clade=="split":
        pass
    else:
        sys.stderr.write(cyan(f'Error: invalid clade specified. Please specify one of `cladei`, `cladeia`,`cladeib`, `cladeii`, `cladeiia`,`cladeiib`\n'))
        sys.exit(-1)

    resources = [
            {"key":KEY_REPORT_TEMPLATE,
            "directory":"data",
            "filename":"report.mako"},
            {"key":KEY_BACKGROUND_FASTA,
            "directory":"data",
            "filename":"background.fasta"},
            {"key":KEY_REFERENCE_PANEL,
            "directory":"data",
            "filename":"background.panel.fasta"},
            {"key":KEY_REFERENCE_FASTA,
            "directory":"data",
            "filename":fasta_filename},
            {"key":KEY_TO_MASK,
            "directory":"data",
            "filename":mask_file},
            {"key":KEY_GENE_BOUNDARIES,
            "directory":"data",
            "filename":gene_boundaries_file},
            {"key":KEY_GRANTHAM_SCORES,
            "directory":"data",
            "filename":"grantham_score.txt"},
            {"key":KEY_ASSEMBLY_REFERENCES,
            "directory":"data",
            "filename":"ref_seq.fasta"}
            ]

    for resource in resources:
        package_data_check(resource["filename"],resource["directory"],resource["key"],config)