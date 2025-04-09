#!/usr/bin/env python3
import sys
import os
from squirrel.utils.log_colours import green,cyan
import select
from Bio import SeqIO
import csv

import tempfile
import shutil

from squirrel.utils.config import *

def set_up_threads(threads,config):
    if threads:
        try:
            tint = int(threads)
            config[KEY_THREADS] = threads
            config[KEY_PHYLO_THREADS] = threads
        except:
            sys.stderr.write(cyan(f'Error: threads specified must be an integer'))
            sys.exit(-1)


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

def set_up_outfile(outfile_arg,cwd,query_arg, outfile, outdir):
    outfile_stem = ""
    outfile_name = ""

    if outfile_arg:
        outfile_arg_dir = ""
        if "/" in outfile_arg:
            # figure out if a path was provided
            p = outfile_arg.split("/")
            outfile_arg_dir = "/".join(p[:-1])
            outfile_arg = p[-1]

        if "." in outfile_arg:
            # strip off extension if provided to get stem, keep for the name variable
            outfile_stem = ".".join(outfile_arg.split(".""")[:-1])
            outfile_name = outfile_arg
        else:
            # add a stem to the aln outfile
            outfile_stem = outfile_arg
            outfile_name = f"{outfile_stem}.aln.fasta"
            
        if outfile_arg_dir:
            # if the outfile provided was a path, join that path to the outdir
            
            outdir = os.path.join(outdir,outfile_arg_dir)
            try:
                if not os.path.exists(outdir):
                    os.mkdir(outdir)
            except:
                sys.stderr.write(cyan(f'Error: cannot create output directory for outfile {outdir}.\n'))
                sys.exit(-1)
    else:
        try:
            if not os.path.exists(os.path.join(cwd, query_arg[0])):
                outfile_stem = "sequences"
                outfile_name = "sequences.aln.fasta"
            else:
                # get the file name
                query_file = query_arg[0].split("/")[-1]
                
                # get the file stem & name
                outfile_stem = ".".join(query_file.split(".")[:-1])
                outfile_name = f'{outfile_stem}.aln.fasta'
        except IndexError:
            sys.stderr.write(cyan(
                f'Error: input query fasta could not be detected from a filepath or through stdin.\n'))
            sys.exit(-1)


    outfile = os.path.join(outdir, outfile_name)

    cds_outstr = f"{outfile_stem}.aln.cds.fasta"
    cds_outfile = os.path.join(outdir, cds_outstr)

    return outfile,cds_outfile,outfile_name,outfile_stem,outdir

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


def find_exclude_file(cwd,input_fasta,exclude_file,config):

    path_to_try = os.path.join(cwd,exclude_file)
    if not os.path.exists(path_to_try):
        sys.stderr.write(cyan(f'Error: cannot find exclude file at: ') + f'{path_to_try}\n' + cyan('Please check file path and try again.\n'))
        sys.exit(-1)

    to_exclude = set()
    with open(path_to_try,"r") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames

        if "name" not in header:
            sys.stderr.write(cyan(f'Error: exclude must contain column `name`.\n'))
            sys.exit(-1)
        for row in reader:
            if row["name"]:
                to_exclude.add(row["name"])
            else:
                sys.stderr.write(cyan(f'Error: empty lines or missing sequence names in the exclude file\n'))
                sys.exit(-1)
    print(green(f"Note: {len(to_exclude)} sequences to exclude"))

    new_input_fasta = os.path.join(config[KEY_TEMPDIR], "input.excluded.fasta")
    ex = 0
    with open(new_input_fasta,"w") as fw:
        i = 0
        ex +=1
        for record in SeqIO.parse(input_fasta,"fasta"):
            if record.description in to_exclude or record.id in to_exclude:
                ex +=1
            else:
                fw.write(f">{record.description}\n{record.seq}\n")
                i+=1
    
    print(green("Input FASTA file filtered by exclude file."))
    print("Number of excluded sequences:",ex)
    print("Number of sequences remaining in alignment:",i)

    return new_input_fasta
    
def find_additional_mask_file(cwd,additional_mask,config):

    path_to_try = os.path.join(cwd,additional_mask)
    if not os.path.exists(path_to_try):
        sys.stderr.write(cyan(f'Error: cannot find additional mask file at: ') + f'{path_to_try}\n' + cyan('Please check file path and try again.\n'))
        sys.exit(-1)

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

    return path_to_try

def find_sequence_mask_file(cwd,sequence_mask,config):

    path_to_try = os.path.join(cwd,sequence_mask)
    if not os.path.exists(path_to_try):
        sys.stderr.write(cyan(f'Error: cannot find sequence mask file at: ') + f'{path_to_try}\n' + cyan('Please check file path and try again.\n'))
        sys.exit(-1)

    with open(path_to_try,"r") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        for i in ["sequence","site"]:
            if i not in header:
                sys.stderr.write(cyan(f'Error: sequence mask file must contain columns `sequence` and `site`.\n'))
                sys.exit(-1)
        for row in reader:
            try:
                site = int(row["site"])
            except:
                sys.stderr.write(cyan(f'Error: `site` column must contain numeric value.\n'))
                sys.exit(-1)

    return path_to_try


def pipeline_options(no_mask, no_itr_mask, additional_mask,sequence_mask,extract_cds,concatenate,clade,cwd, config):
    config[KEY_NO_MASK] = no_mask
    if no_itr_mask:
        if clade.startswith("cladeii"):
            config[KEY_TRIM_END] = 197209
        elif clade.startswith("cladei"):
            config[KEY_TRIM_END] = 196858
        elif clade == "variola":
            config[KEY_TRIM_END] = 185578

    
    if additional_mask:
        config[KEY_ADDITIONAL_MASK] = find_additional_mask_file(cwd,additional_mask,config)

    if sequence_mask:
        config[KEY_SEQUENCE_MASK] = find_sequence_mask_file(cwd,sequence_mask,config)

    config[KEY_EXTRACT_CDS] = extract_cds
    config[KEY_CONCATENATE] = concatenate


def find_background_file(cwd,input_fasta,background_file,config):
    seqs = set()
    path_to_try = os.path.join(cwd,background_file)

    new_input_fasta = os.path.join(config[KEY_TEMPDIR], "input.custom_background.combined.fasta")
    with open(new_input_fasta,"w") as fw:
        i = 0
        for record in SeqIO.parse(input_fasta,"fasta"):
            seqs.add(record.description)
            fw.write(f">{record.description}\n{record.seq}\n")
            i+=1
        c = 0
        try:
            for record in SeqIO.parse(path_to_try,"fasta"):
                if record.description in seqs:
                    print(cyan("Ignoring duplicate seq in background:"),record.description)
                    seqs.add(record.description)
                else:
                    c +=1
                    fw.write(f">{record.description}\n{record.seq}\n")
        except:
            sys.stderr.write(cyan(f'Error: cannot find/parse background fasta file at: ') + f'{path_to_try}\n' + cyan('Please check file path and format.\n'))
            sys.exit(-1)
    
    print(green("Custom background combined with input FASTA file."))
    print("Number of input sequences:",i)
    print("Number of background sequences:",c)

    return new_input_fasta

def set_up_clade(clade,config):
    if clade:
        config[KEY_CLADE] = clade
        if config[KEY_CLADE] not in VALUE_VALID_CLADES:
            sys.stderr.write(cyan(
                f'Error: clade must be one of {VALUE_VALID_CLADES}.\n'))
            sys.exit(-1) 

def add_background_to_input(input_fasta,background,clade,config):
    in_name = input_fasta.rstrip("fasta").split("/")[-1]
    new_input_fasta = os.path.join(config[KEY_TEMPDIR], f"{in_name}.background_included.fasta")

    added = set()
    with open(new_input_fasta,"w") as fw:
        for record in SeqIO.parse(background,"fasta"):
            # include the outgroup seq
            if record.id in config[KEY_OUTGROUPS]:
                print("writing outgroup",record.id)
                fw.write(f">{record.id}\n{record.seq}\n")
                added.add(record.id)
            else:
                # include the relevant clade seqs
                c = ""
                for field in record.description.split(" "):
                    if field.startswith("clade"):
                        c = field.split("=")[1].lower()
                if c not in VALUE_VALID_CLADES:
                    sys.stderr.write(cyan(
                        f'Error: clade must be one of {VALUE_VALID_CLADES}.\n'))
                    sys.exit(-1)

                if clade == "cladei":
                    if c in ["cladei","cladeia","cladeib"]:
                        fw.write(f">{record.id}\n{record.seq}\n")
                        added.add(record.id)
                elif clade == "cladeii":
                    if c in ["cladeii","cladeiia","cladeiib"]:
                        fw.write(f">{record.id}\n{record.seq}\n")
                        added.add(record.id)
                else:
                    if c == clade:
                        fw.write(f">{record.id}\n{record.seq}\n")
                        added.add(record.id)

        for record in SeqIO.parse(input_fasta,"fasta"):
            for_iqtree = record.description.replace(" ","_")
            if for_iqtree in added:
                sys.stderr.write(cyan(f'Error: duplicate sequence name `{for_iqtree}` in background and supplied file.\nPlease modify sequence name and try again.\n'))
                sys.exit(-1)
            fw.write(f">{record.description}\n{record.seq}\n")
            

    return new_input_fasta

def parse_tf_options(tree_figure_only,tree_file,branch_reconstruction_file,width,height,point_style_arg,justify_arg,cwd,config):
    if point_style_arg:
        config[KEY_POINT_STYLE] = point_style_arg
    
    if  config[KEY_POINT_STYLE] not in ["circle","square"]:
        sys.stderr.write(cyan(f'Error: not a valid point style, please specify one of `circle` or `square`.\n'))
        sys.exit(-1)

    if justify_arg:
        config[KEY_POINT_JUSTIFY] = justify_arg
    
    if  config[KEY_POINT_JUSTIFY] not in ["left","right"]:
        sys.stderr.write(cyan(f'Error: not a valid point justification, please specify one of `left` or `right`.\n'))
        sys.exit(-1)

    if width:
        config[KEY_FIG_WIDTH] = width
    if height:
        config[KEY_FIG_HEIGHT] = height
    
    if tree_figure_only:
        if not tree_file or not branch_reconstruction_file:
            sys.stderr.write(cyan(f'Error: must supply a tree and branch reconstruction file alongside `-tfig/--tree-figure-only`.\n'))
            sys.exit(-1)

        tree = os.path.join(cwd, tree_file)
        if not os.path.exists(tree):
            sys.stderr.write(cyan(f'Error: cannot find treefile:') + f"{tree}")
            sys.exit(-1)
        config[KEY_TREE] = tree

        branch_reconstruction = os.path.join(cwd, branch_reconstruction_file)
        if not os.path.exists(branch_reconstruction):
            sys.stderr.write(cyan(f'Error: cannot find branch reconstruction file:') + f"{branch_reconstruction}")
            sys.exit(-1)
        config[KEY_BRANCH_RECONSTRUCTION] = branch_reconstruction


def phylo_options(run_phylo,run_apobec3_phylo,outgroups,include_background,binary_partition_mask,input_fasta,config):
    config[KEY_RUN_PHYLO] = run_phylo

    if run_apobec3_phylo:
        config[KEY_RUN_APOBEC3_PHYLO] = run_apobec3_phylo
        config[KEY_RUN_PHYLO] = True

    if binary_partition_mask and not run_apobec3_phylo:
        sys.stderr.write(cyan(
                        f'Error: binary partition mask file can only be written if APOBEC3 reconstruction mode is on.\n'))
        sys.exit(-1)


    if config[KEY_RUN_PHYLO]:

        if include_background:
            if outgroups:
                print(cyan('Note: outgroup(s) are selected automatically with using `--include-background` flag.\n'))
            config[KEY_INCLUDE_BACKGROUND] = "True"
            clade = config[KEY_CLADE]
            outgroups = OUTGROUP_DICT[clade]
            print(green("Outgroup selected:"),outgroups[0])

        if not outgroups:
            sys.stderr.write(cyan(
                        f'Error: must supply outgroup(s) for phylogenetics module.\n'))
            sys.exit(-1)
        if not type(outgroups) == list:
            outgroups = outgroups.split(",")
        config[KEY_OUTGROUPS] = outgroups
        

        if include_background:
            new_input_fasta = add_background_to_input(input_fasta,config[KEY_BACKGROUND_FASTA],config[KEY_CLADE],config)
            return new_input_fasta

        seqs = SeqIO.index(input_fasta,"fasta")
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

    
    return input_fasta
