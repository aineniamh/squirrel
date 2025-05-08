import os
from squirrel.utils.config import *
from squirrel.utils.log_colours import green,cyan
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import date
import datetime as dt

import csv

import baltic as bt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import collections
import pandas as pd
import matplotlib.patches as patches

import math
plt.switch_backend('Agg') 


mpl.rcParams.update({'font.size': 18})

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}
mpl.rcParams.update(new_rc_params)
plt.rcParams['font.family'] = 'Helvetica'

    
def get_node_states_all_sites(state_file,alignment):
    
    #returns a dict keys off 1-based positions in the genome
    #and the value is a list of tuples (id, base) at a given node
    # allows you to look up for a given site what the base is for a
    #given internal node or tip
    
    node_states = collections.defaultdict(list)
    c = 0
    
    ## first the reconstructed nodes
    with open(state_file,"r") as f:
        for l in f:

            if not l.startswith("#"):
                c+=1
                try:
                    node,site,state,probA,probC,probG,probT = l.rstrip("\n").split("\t")
                except:
                    print(l)
                    break
                if node != "Node":
                    if state not in ["N","-"]:
                        node_states[site].append((node,state))
                    else:
                        node_states[site].append((node,""))
    ## now the tips
    for record in SeqIO.parse(alignment,"fasta"):
        for site in node_states:
            index = int(site)-1
            base = record.seq[index]
            if base in ["T","C","A","G"]:
                node_states[site].append((record.id,base))
            else:
                node_states[site].append((record.id,""))
                
    return node_states

def get_header_str(dict_values):
    header_str = ""
    for i in sorted(dict_values, key = lambda i : i[0]):
        header_str += f"{i[0]},"
    header_str = header_str.rstrip(",")
    return header_str
    
    
def find_what_sites_vary_unambiguously(node_states,state_differences):
    header_str = get_header_str(node_states["1"])
    
    with open(state_differences,"w") as fw:
        fw.write(f"site,{header_str}\n")

        for site in node_states:
            info = node_states[site]
            
            # get the set of unique bases at a given site
            count = set([i[1] for i in info if i[1]])
            
            #if there's more than one
            if len(count)>1:
                
                #needs to be kep consistent with header str
                info = sorted(info, key = lambda i : i[0])
                base_str = ""
                for i in info:
                    base_str += f"{i[1]},"
                    
                base_str = base_str.rstrip(",")
                fw.write(f"{site},{base_str}\n")
    
def load_unambiguous_varying_sites(infile):
    node_states_diff = collections.defaultdict(dict)
    with open(infile,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            site = row["site"]
            for col in row:
                if col != "site":
                    node_states_diff[row["site"]][col] = row[col]
    return node_states_diff

def map_site_changes_to_branches(treefile, outfile,node_states,node_states_diff): 
    my_tree=bt.loadNewick(treefile,absoluteTime=False)
    last_node = ""
    current_node = ""

    with open(outfile,"w") as fw:
        fw.write("parent,child,site,snp,dimer\n")

        for k in my_tree.Objects:
            if k.branchType == 'leaf':
                current_node = k
                current_node.traits["label"]=k.name
            else:
                current_node = k

            if last_node:
                node_name = current_node.traits["label"]
                parent_name = current_node.parent.traits["label"]
                snps = []
                for site in node_states_diff:
                    node_base = node_states_diff[site][node_name]
                    parent_base = node_states_diff[site][parent_name]

                    if node_base != parent_base:
                        if node_base in ["A","C","G","T"] and parent_base in ["A","C","G","T"]:
                            snp = f"{parent_base}->{node_base}"
                            snps.append(snp)
                            if snp == "G->A":
                                dimer_site = f"{int(site)+1}"
                                dimer_base = ""

                                for i in node_states[dimer_site]:
                                    if i[0] == parent_name:
                                        dimer_base = i[1]
                                dimer = f"{parent_base}{dimer_base}"
                            elif snp == "C->T":
                                dimer_site = f"{int(site)-1}"
                                dimer_base = ""

                                for i in node_states[dimer_site]:
                                    if i[0] == parent_name:
                                        dimer_base = i[1]
                                dimer = f"{dimer_base}{parent_base}"
                            else:
                                dimer = ""
                            fw.write(f"{parent_name},{node_name},{site},{snp},{dimer}\n")

            last_node = current_node

def read_in_branch_snps(branch_snps):
    branch_snps_dict = collections.defaultdict(list)
    with open(branch_snps,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            branch = f"{row['parent']}_{row['child']}"
            branch_snps_dict[branch].append((row['site'],row['snp'],row['dimer'])) 
    return branch_snps_dict

def get_branch_snps_sites(branch_snps):
    all_snps = collections.Counter()
    branch_snps_dict = collections.defaultdict(list)
    with open(branch_snps,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            all_snps[int(row["site"])]+=1

            branch_snps_dict[int(row['site'])].append([row['parent'],row['child'],row['snp'],row['dimer']])
    
    homoplasies = {}
    for k in all_snps:
        if all_snps[k]>1:
            homoplasies[k] = all_snps[k]
            
    # print(len(homoplasies))
    # print(homoplasies)
    return branch_snps_dict,homoplasies
    

def get_acc_to_metadata_map(metadata):
    acc_dict = {}
    with open(metadata,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                acc_dict[row["accession"]] = row
            except:
                try:
                    acc_dict[row["\ufeffaccession"]] = row
                except:
                    pass
                
    return acc_dict

def get_fig_height(alignment):
    seqs = SeqIO.index(alignment,"fasta")

    height = 0.5*len(seqs)
    if height >15:
        return height
    else:
        return 15


def make_reconstruction_tree_figure_w_labels(outfile,branch_snps,treefile,point_style,justification,w=None,h=None):
    branch_snps_dict = read_in_branch_snps(branch_snps)
    
    my_tree=bt.loadNewick(treefile,absoluteTime=False)
    plt.switch_backend('Agg') 


    x_attr=lambda k: k.height ## x coordinate of branches will be absoluteTime attribute
    su_func=lambda k: 50-30*k.height/my_tree.treeHeight ## size of tips
    s_func=lambda k: 50-20*k.height/my_tree.treeHeight ## size of tips
    c_func=lambda k: "dimgrey"
    
    c_func=lambda k: 'dimgrey' 
    
    r2t = 200000*my_tree.treeHeight #rough number of snps root to tip
    increment = my_tree.treeHeight/(r2t*1.5) # divide the tree height by about twice the num of r2t snps
    # print(w,h)
    if w==None or w=="None":
        if r2t < 200:
            width = int(math.sqrt(r2t)*3)
        else:
            width = 25
    else:
        width = int(w)

    if h == None or h =="None":
        if my_tree.ySpan < 300:
            height = int(math.sqrt(my_tree.ySpan)*2)
        else:
            height = 40
    else:
        height = int(h)


    if not np.isfinite(height):
        height = 40
    if not np.isfinite(width):
        width = 25
    

    fig,ax = plt.subplots(figsize=(width,height),facecolor='w')

    my_tree.plotTree(ax,x_attr=x_attr) ## plot branches
    my_tree.plotPoints(ax,size=s_func,colour=c_func,x_attr=x_attr) ## plot circles at tips
    mpl.rcParams['font.family'] = 'sans-serif'
    
    target_func=lambda k: k.is_leaf() ## which branches will be annotated
    text_func=lambda k: k.name ## what text is plotted
    text_x_attr=lambda k: k.height+(increment*4) ## where x coordinate for text is

    my_tree.addText(ax,x_attr=text_x_attr,target=target_func,text=text_func) #

    
    for k in my_tree.Objects:
        current_node = k
        if k.branchType == 'leaf':
            current_node.traits["label"]=k.name

        node_name = current_node.traits["label"]
        try:
            parent_name = current_node.parent.traits["label"]
        except:
            continue
        branch_name= f"{parent_name}_{node_name}"

        if branch_name in branch_snps_dict:
            snps = []
#                 print(branch_name, len(branch_snps_dict[branch_name]))
            right_settings = {
                "apobec":(1,"#995E62"),
                "non_apobec":(2,"#D9B660")
            }
            left_settings = {
                "apobec":(2,"#995E62"),
                "non_apobec":(1,"#D9B660")
            }
            if justification == "right":
                setting_dict = right_settings
                snp_placement = current_node.height - increment
            else:
                setting_dict = left_settings
                snp_placement = current_node.parent.height + increment/2

            for s in branch_snps_dict[branch_name]:
                site,snp,dimer = s
                if snp == "G->A":
                    if dimer in ["GA"]:
                        snps.append(setting_dict["apobec"])
                    else:
                        snps.append(setting_dict["non_apobec"])
                elif snp == "C->T":
                    if dimer in ["TC"]:
                        snps.append(setting_dict["apobec"])
                    else:
                        snps.append(setting_dict["non_apobec"])
                else:
                    snps.append(setting_dict["non_apobec"])

            for snp in sorted(snps, key = lambda x : x[0]):
                
                if point_style == "circle":
                    plt.scatter([snp_placement],[k.y+0.5],color=snp[1],s=50)
                else:
                    rect = patches.Rectangle((snp_placement,k.y-0.5),increment/2,1,alpha=1, fill=True, edgecolor='none',facecolor=snp[1])
                    ax.add_patch(rect)

                if justification == "right":
                    snp_placement -= increment
                else:
                    snp_placement += increment

    [ax.spines[loc].set_visible(False) for loc in ['top','right','left']]
    ax.tick_params(axis='y',size=0)
    ax.tick_params(axis='x',size=0)

    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.margins(0.005,0.005,tight=True)
    plt.savefig(f"{outfile}.svg",bbox_inches='tight')
    plt.savefig(f"{outfile}.png",bbox_inches='tight'
                   );
    
def generate_reconstruction_files(alignment, state_out, state_differences):
    
    node_states = get_node_states_all_sites(state_out,alignment)
        
    find_what_sites_vary_unambiguously(node_states,state_differences)

    return node_states
    
def load_info(directory, alignment, treefile, state_out, state_differences, branch_snps_out, treefigureout,point_style,point_justify, node_states="",width=None,height=None):
    
    if not node_states:
        node_states = get_node_states_all_sites(state_out, alignment)

    node_states_diff = load_unambiguous_varying_sites(state_differences)

    map_site_changes_to_branches(treefile,
                                 branch_snps_out,
                                 node_states,
                                 node_states_diff)

    make_reconstruction_tree_figure_w_labels(treefigureout,
                                    branch_snps_out,
                                    treefile,
                                    point_style,
                                    point_justify,
                                    height,
                                    width)
    
def get_gene_boundaries(gene_boundaries_file):
    genes = {}
    gene_id = 0
    with open(gene_boundaries_file,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_id +=1
            name = f"{row['Name'].replace(' ','_')}_{gene_id}"
            start = int(row["Minimum"])
            end = int(row["Maximum"])+1
            length = int(row["Length"])
            direction = row["Direction"]
            genes[(start,end)]=(name,length,direction)
    return genes

def get_grantham_scores(grantham_scores_file):
    grantham_scores = {}

    with open(grantham_scores_file,"r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            for col in row:
                if col!="FIRST":
                    mutation = f"{row['FIRST']}{col}"

                    if row[col] != "0":
                        grantham_scores[mutation] = int(row[col])                    
                        grantham_scores[mutation[::-1]] = int(row[col])
    return grantham_scores
 
    
def categorise_amino_acid_mutation(aa1,aa2,grantham_scores):
    
    mutation_category = ""
    if aa1 == aa2:
        mutation_category = "synonymous"
    else:
        if aa2 == '*':
            mutation_category = "nonsense"
        else:
            mutation_category = "nonsynonymous"

    if f"{aa1}{aa2}" in grantham_scores:
        score = grantham_scores[f"{aa1}{aa2}"]
        if score < 51:
            prediction = "conservative"
        elif score <101:
            prediction = "moderately conservative"
        elif score <151:
            prediction = "moderately radical"
        else:
            prediction = "radical"
    else:
        score = "NA"
        prediction = "NA"
    
    return mutation_category,score,prediction

  
def get_aa_position(index,based=0):
    if based != 0:
        index -= based 

    remainder = index%3
    position_dict = {1:2,2:3,0:1}
    
    return position_dict[remainder]

def reverse_aa_position(start,end,site):
    index_dict = {}
    index = 0
    rev_position_dict = {1:2,2:3,0:1}
    for i in reversed(range(start,end)):
        remainder = index%3
        index_dict[i] = rev_position_dict[remainder]
        index +=1 

    position = index_dict[site]
    
    return position


def get_codon_indexes(aa_position,index):
    codon = []
    if aa_position == 1:
        codon = [index,index+1,index+2]
    elif aa_position == 2:
        codon = [index-1,index,index+1]
    elif aa_position == 3:
        codon =  [index-2,index-1,index]
    else:
        print("incorrect aa position")
        
    return codon

def get_codon_indexes_rev_strand(position,index):
    codon = []
    if position == 3:
        codon = [index,index+1,index+2]
    elif position == 2:
        codon = [index-1,index,index+1]
    elif position == 1:
        codon =  [index-2,index-1,index]
    return codon

def reconstruct_amino_acid_mutations(grantham_scores_file,gene_boundaries_file,branch_snps,node_states,outfile):
    branch_snps_dict,homoplasies = get_branch_snps_sites(branch_snps)
    genes = get_gene_boundaries(gene_boundaries_file)
    grantham_scores = get_grantham_scores(grantham_scores_file)
    
    fw = open(outfile,"w")
    fw.write("site,gene,direction,snp,dimer,apobec,aa_position,parent,parent_codon,parent_aa,")
    fw.write("child,child_codon,child_aa,mutation_category,score,prediction,homoplasy,occurrence\n")
    
    for site in branch_snps_dict:
        
        homoplasy = "False"
        occurrence = "1"
        if site in homoplasies:
            homoplasy = "True"
            occurrence = f"{homoplasies[site]}"
        site_found = False
        for gene in genes:
            start,end=gene

            if site in range(start,end):
                
                site_found = True
                name,length,direction = genes[gene]
                for site_snp in branch_snps_dict[site]:
                    parent,child,snp,dimer = site_snp

                    if direction == "forward":
                        aa_position = get_aa_position(site,gene[0])
                        codon_indexes = get_codon_indexes(aa_position,site)
                    else:
                        aa_position = reverse_aa_position(start,end,site)
                        codon_indexes = get_codon_indexes_rev_strand(aa_position,site)

                    parent_codon = []
                    child_codon = []

                    for base in codon_indexes:
                        reconstruction = node_states[f"{base}"]
                        for node in reconstruction:
                            if node[0] == parent:
                                parent_codon.append(node[1])

                            elif node[0] == child:
                                child_codon.append(node[1])

                    parent_codon = "".join(parent_codon)
                    child_codon = "".join(child_codon)

                    parent_codon = Seq(parent_codon)
                    child_codon = Seq(child_codon)

                    if direction == "reverse":
                        parent_codon= parent_codon.reverse_complement()
                        child_codon= child_codon.reverse_complement()

                    parent_aa = parent_codon.translate()
                    child_aa = child_codon.translate()

                    mutation_category,score,prediction = categorise_amino_acid_mutation(parent_aa,child_aa,grantham_scores)
                    apobec = "False"
                    if snp in ["C->T","G->A"] and dimer in ["GA","TC"]:
                        apobec = "True"

                    fw.write(f"{site},{name},{direction},{snp},{dimer},{apobec},{aa_position},{parent},{parent_codon},{parent_aa},{child},{child_codon},{child_aa},{mutation_category},{score},{prediction},{homoplasy},{occurrence}\n")                                    
        
        if not site_found:
            for snp_site in branch_snps_dict[site]:
                parent,child,snp,dimer = snp_site
                apobec = "False"
                if snp in ["C->T","G->A"] and dimer in ["GA","TC"]:
                    apobec = "True"
                fw.write(f"{site},NA,NA,{snp},{dimer},{apobec},NA,{parent},NA,NA,{child},NA,NA,intergenic,NA,NA,{homoplasy},{occurrence}\n")                                    
            
    fw.close()
            
def get_reconstruction_amino_acids(alignment,grantham_scores_file,gene_boundaries_file,branch_snps_out,state_out,amino_acids_out,node_states=""):
    if not node_states:
        node_states = get_node_states_all_sites(state_out,alignment)

    reconstruct_amino_acid_mutations(grantham_scores_file,gene_boundaries_file,branch_snps_out,
                                    node_states, amino_acids_out)
    
    
def year_fraction(date):
    start = dt.date(date.year, 1, 1).toordinal()
    year_length = dt.date(date.year+1, 1, 1).toordinal() - start
    return date.year + float(date.toordinal() - start) / year_length

def get_root_to_tip_counts(aa_reconstruction,state_diffs,root_to_tip_counts):

    site_info = {}
    with open(aa_reconstruction,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            site_info[row["site"]] = row

    # for i in site_info:
    #     print(i, site_info[i])

    seq_snps = collections.defaultdict(list)
    apobec_snps = collections.defaultdict(list)
    non_apo = collections.defaultdict(list)
    with open(state_diffs,"r") as f:
        reader = csv.DictReader(f)

        for row in reader:
            root_variant = row["Node1"]

            for seq in reader.fieldnames:
                if seq != "site" and not seq.startswith("Node"):
                    if row[seq] != root_variant:
                        if "" not in [row[seq],root_variant]:
                            try:
                                snp_info = site_info[row["site"]]
                                if snp_info["apobec"] == "True":
                                    apobec_snps[seq].append(row["site"])
                                else:
                                    non_apo[seq].append(row["site"])
                                seq_snps[seq].append(row["site"])
                            except:
                                pass
    s = []
    date_apo = {}
    with open(root_to_tip_counts,"w") as fw:
        fw.write("name,all_snps,apobec_snps,non_apobec_snps,date,decimal_year,precision\n")
        for i in seq_snps:
            datestring = i.split("|")[-1]
            
            precision = ""
            
            if len(datestring.split("-")) == 3:
                ddate = date.fromisoformat(datestring)
                odate = year_fraction(ddate)
                
                precision = "day"
                print(datestring,odate,precision)
            elif len(datestring.split("-")) == 2:
                if datestring.endswith("-2"):
                    datestring+="-15"
                else:
                    datestring+="-16"
                ddate = date.fromisoformat(datestring)
                odate = year_fraction(ddate)
                
                precision = "month"
                print(datestring,odate,precision)
            else:
                
                odate = float(datestring)+0.5
                precision="year"
                print(datestring,odate,precision)
            if i in apobec_snps:
                if i in non_apo:
                    fw.write(f"{i},{len(seq_snps[i])},{len(apobec_snps[i])},{len(non_apo[i])},{datestring},{odate},{precision}\n")
                else:
                    fw.write(f"{i},{len(seq_snps[i])},{len(apobec_snps[i])},0,{datestring},{odate},{precision}\n")
                date_apo[i] = [apobec_snps[i],datestring]
            else:
                if i in non_apo:
                    fw.write(f"{i},{len(seq_snps[i])},0,{len(non_apo[i])},{datestring},{odate},{precision}\n")
                else:
                    fw.write(f"{i},{len(seq_snps[i])},0,0,{datestring},{odate},{precision}\n")


def get_root_to_tip_counts_date_in(aa_reconstruction,state_diffs,root_to_tip_counts,date_dict):

    site_info = {}
    with open(aa_reconstruction,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            site_info[row["site"]] = row

    # for i in site_info:
    #     print(i, site_info[i])

    seq_snps = collections.defaultdict(list)
    apobec_snps = collections.defaultdict(list)
    non_apo = collections.defaultdict(list)
    with open(state_diffs,"r") as f:
        reader = csv.DictReader(f)

        for row in reader:
            root_variant = row["Node1"]

            for seq in reader.fieldnames:
                if seq != "site" and not seq.startswith("Node"):
                    if row[seq] != root_variant:
                        if "" not in [row[seq],root_variant]:
                            try:
                                snp_info = site_info[row["site"]]
                                if snp_info["apobec"] == "True":
                                    apobec_snps[seq].append(row["site"])
                                else:
                                    non_apo[seq].append(row["site"])
                                seq_snps[seq].append(row["site"])
                            except:
                                pass
    s = []
    date_apo = {}
    fw2 = open(APO_out,"w")
    fw2.write("name,date\n")
    with open(root_to_tip_counts,"w") as fw:
        fw.write("name,all_snps,apobec_snps,non_apobec_snps,date,decimal_year,precision\n")
        for i in seq_snps:
            datestring = date_dict[i]
            
            precision = ""
            
            if len(datestring.split("-")) == 3:
                ddate = date.fromisoformat(datestring)
                odate = year_fraction(ddate)
                
                precision = "day"
                print(datestring,odate,precision)
            elif len(datestring.split("-")) == 2:
                if datestring.endswith("-2"):
                    datestring+="-15"
                else:
                    datestring+="-16"
                ddate = date.fromisoformat(datestring)
                odate = year_fraction(ddate)
                
                precision = "month"
                print(datestring,odate,precision)
            else:
                
                odate = float(datestring)+0.5
                precision="year"
                print(datestring,odate,precision)
            if i in apobec_snps:
                if i in non_apo:
                    fw.write(f"{i},{len(seq_snps[i])},{len(apobec_snps[i])},{len(non_apo[i])},{datestring},{odate},{precision}\n")
                else:
                    fw.write(f"{i},{len(seq_snps[i])},{len(apobec_snps[i])},0,{datestring},{odate},{precision}\n")
                date_apo[i] = [apobec_snps[i],datestring]
            else:
                if i in non_apo:
                    fw.write(f"{i},{len(seq_snps[i])},0,{len(non_apo[i])},{datestring},{odate},{precision}\n")
                else:
                    fw.write(f"{i},{len(seq_snps[i])},0,0,{datestring},{odate},{precision}\n")


def run_full_analysis(directory, alignment, treefile,state_file,config,point_style,point_justify,width,height):

    state_out = state_file
    state_differences = f"{treefile}.state_differences.csv"
    branch_snps_out = f"{treefile}.branch_snps.reconstruction.csv"
    amino_acids_out= f"{treefile}.amino_acid.reconstruction.csv"

    node_states= generate_reconstruction_files(alignment,
                                  state_out,
                                  state_differences)

    tree_fig = f"{treefile}"
    load_info(directory,alignment,treefile,state_out,state_differences,branch_snps_out,tree_fig,point_style,point_justify,node_states,width,height)


    grantham_scores_file = config[KEY_GRANTHAM_SCORES]
    gene_boundaries_file = config[KEY_GENE_BOUNDARIES]
    get_reconstruction_amino_acids(alignment,grantham_scores_file,gene_boundaries_file,branch_snps_out,state_out,amino_acids_out,node_states)

def find_binary_partition_mask(branch_reconstruction,sep_status,reference,outfile):

    apobec_variable_ga = set()
    apobec_variable_tc = set()
    non_apobec_variable =  set()

    if sep_status:
        tc_mask = "2"
    else:
        tc_mask = "1"

    with open(branch_reconstruction,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["dimer"] == "GA":
                apobec_variable_ga.add(int(row['site']))
            elif row["dimer"] == "TC":
                apobec_variable_tc.add(int(row['site']))
            else:
                non_apobec_variable.add(int(row['site']))

    ref = str(SeqIO.read(reference, 'fasta').seq)
    pos = np.arange(len(ref))
    
    apo_ga_keep_0index = set()
    apo_tc_keep_0index = set()
    for i in pos:
        if ref[i:i+2]=="GA":
            if i+1 not in non_apobec_variable:
                apo_ga_keep_0index.add(i)
        elif ref[i:i+2]=="TC":
            if i+2 not in non_apobec_variable:
                apo_tc_keep_0index.add(i+1)

    for i in apobec_variable_ga:
        apo_ga_keep_0index.add(i-1)

    for i in apobec_variable_tc:
        apo_tc_keep_0index.add(i-1)

    apo_masked = 0
    non_apo_masked = 0
    tc_masked = 0
    ga_masked = 0
    mask_string = ""
    
    for i in range(len(ref)):
        if i not in apo_ga_keep_0index and i not in apo_tc_keep_0index:
            mask_string+="0"
            non_apo_masked +=1
        elif i in apo_ga_keep_0index:
            mask_string+="1"
            apo_masked +=1
            ga_masked +=1
        elif i in apo_tc_keep_0index:
            mask_string+=tc_mask
            apo_masked +=1
            tc_masked +=1
    print("TC sites",tc_masked)
    print("GA sites",ga_masked)
    print("All APOBEC3 sites",apo_masked)
    print("Non APOBEC3 sites",non_apo_masked)

    with open(outfile,"w") as fw:
        fw.write(mask_string + "\n")
    