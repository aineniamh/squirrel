#!/usr/bin/env python3
import os
import collections
from squirrel.utils.log_colours import green,cyan
from Bio import SeqIO
import csv

from squirrel.utils.config import *

import re

def read_fasta(path):
    """
    Read sequences from a FASTA file.
    
    Args:
        path (str): Path to the FASTA file.
    
    Returns:
        dict: Dictionary mapping sequence IDs to sequence strings.
    """
    sequences = {}
    for record in SeqIO.parse(path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences
    
def cigar_from_alignment(aligned_seq, gap='-'):
    """
    Generate a CIGAR string from an aligned sequence.
    
    Converts an aligned sequence into CIGAR format where 'M' represents
    a match (or mismatch) and 'I' represents an insertion (gap).
    
    Args:
        aligned_seq (str): The aligned sequence string.
        gap (str): The gap character used in the alignment. Defaults to '-'.
    
    Returns:
        str: CIGAR string representation of the alignment.
    """
    cigar = []
    curr_op = None
    curr_len = 0

    for base in aligned_seq:
        op = 'I' if base == gap else 'M'

        if op == curr_op:
            curr_len += 1
        else:
            if curr_op is not None:
                cigar.append(f"{curr_len}{curr_op}")
            curr_op = op
            curr_len = 1

    cigar.append(f"{curr_len}{curr_op}")
    return ''.join(cigar)

def build_coordinate_maps(aligned_seq, gap='-'):
    """
    Build bidirectional coordinate mapping between genome and alignment positions.
    
    Creates two mapping dictionaries that convert between genomic coordinates
    (ungapped) and alignment coordinates (gapped).
    
    Args:
        aligned_seq (str): The aligned sequence string.
        gap (str): The gap character used in the alignment. Defaults to '-'.
    
    Returns:
        tuple: A tuple of two dictionaries:
            - genome_to_aln (dict): Maps genome positions to alignment positions.
            - aln_to_genome (dict): Maps alignment positions to genome positions,
              with None for gap positions.
    """
    genome_pos = 0
    aln_pos = 0

    genome_to_aln = {}
    aln_to_genome = {}

    for base in aligned_seq:
        aln_pos += 1

        if base == gap:
            aln_to_genome[aln_pos] = None
        else:
            genome_pos += 1
            genome_to_aln[genome_pos] = aln_pos
            aln_to_genome[aln_pos] = genome_pos

    return genome_to_aln, aln_to_genome

def get_coordinates(path):
    """
    Process a multiple sequence alignment file and extract coordinate mappings.
    
    Reads a FASTA alignment file and generates CIGAR strings and bidirectional
    coordinate mappings for each sequence in the alignment.
    
    Args:
        path (str): Path to the FASTA alignment file.
    
    Returns:
        dict: Dictionary mapping sequence IDs to coordinate information dictionaries,
            each containing:
            - aligned_seq (str): The aligned sequence.
            - cigar (str): CIGAR string representation.
            - genome_to_aln (dict): Genome to alignment position mapping.
            - aln_to_genome (dict): Alignment to genome position mapping.
    
    Raises:
        ValueError: If sequences in the alignment have different lengths.
    """
    msa = read_fasta(path)
    results = {}

    aln_lengths = {len(seq) for seq in msa.values()}
    if len(aln_lengths) != 1:
        raise ValueError("Alignment sequences are not the same length")

    for genome_id, aligned_seq in msa.items():
        cigar = cigar_from_alignment(aligned_seq)
        g2a, a2g = build_coordinate_maps(aligned_seq)

        results[genome_id] = {
            "aligned_seq": aligned_seq,
            "cigar": cigar,
            "genome_to_aln": g2a,
            "aln_to_genome": a2g,
        }

    return results


def maps_from_cigar(cigar):
    genome_to_aln = {}
    aln_to_genome = {}

    g_pos = 0
    a_pos = 0

    for length, op in re.findall(r'(\d+)([MI])', cigar):
        length = int(length)

        if op == 'M':
            for _ in range(length):
                g_pos += 1
                a_pos += 1
                genome_to_aln[g_pos] = a_pos
                aln_to_genome[a_pos] = g_pos

        elif op == 'I':
            for _ in range(length):
                a_pos += 1
                aln_to_genome[a_pos] = None

    return genome_to_aln, aln_to_genome

def cigar_dict():
    cigars = {'NC_063383': '86I621M6I134M3I935M2I2907M1I1653M2I277M1I439M28I17M1I450M3I1166M14I664M1I300M2I868M65I2863M15I3109M4I2498M1953I443M2I3123M6I703M24I2380M3I4213M15I17980M5I32M3I158M253I3482M1I10624M1I27820M9I17811M1I24767M9I4108M27I15M6I6851M15I2628M25I28M1I25M1I839M1I495M3I853M1I2288M59I17M43I46M7I4869M8I306M1I197M2I2806M1I1210M3I3204M30I2868M7I2088M1I43M1I1379M4I1329M14I2225M1I536M1I774M1I3383M18I193M2I114M2I1071M1I9M1I2081M2I5765M1I88M1I775M1I406M5I229M1I56M1I545M1I879M2I895M2I1656M1I2903M2I942M3I134M6I621M86I', 
            'PX667572': '679M10I1644M3I2454M83I1655M2I615M6I120M4I258M1I646M1I276M2I92M4I48M1I22M1I432M1I78M1I290M1I240M1I3819M1I90M1I2193M2I319M1I744M6I1441M3I8234M1I28M1I4249M19I3622M3I981M1I214M1I5663M1I2856M1I4576M3I81M1I52M1I929M1I126M2I721M3I2334M2I74M2I20599M1I61655M1I266M15I71M10I319M453I5511M3I678M9I2638M1I720M1I2433M1I83M3I131M2I770M14I96M1I1375M1I2195M1I591M5I22M1I1222M3I2083M2I1064M1I799M1I163M2264I4560M8I1886M2I1606M2I54M1I4208M1I1158M1I1215M32I192M2I241M1I779M12I3817M1I2783M9I2438M3I3653M1I562M1I1977M1I460M2I196M1I2040M83I2451M3I1641M10I685M'}
    return cigars

def get_coordinate_map():
    """
    Extract coordinate mappings from predefined CIGAR strings.
    
    Retrieves a dictionary of CIGAR strings for reference sequences and
    converts each CIGAR string into bidirectional coordinate mappings
    between genomic and alignment positions.
    
    Returns:
        dict: Dictionary mapping sequence IDs to coordinate information dictionaries,
            each containing:
            - cigar (str): CIGAR string representation.
            - genome_to_aln (dict): Maps genome positions to alignment positions.
            - aln_to_genome (dict): Maps alignment positions to genome positions,
              with None for gap positions.
    """
    cigar_dict = cigar_dict()

    results = {}

    for genome_id, cigar in cigar_dict.items():
        g2a, a2g = maps_from_cigar(cigar)

        results[genome_id] = {
            "cigar": cigar,
            "genome_to_aln": g2a,
            "aln_to_genome": a2g,
        }

    return results

