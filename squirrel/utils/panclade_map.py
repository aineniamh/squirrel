#!/usr/bin/env python3
import os
import collections
from squirrel.utils.log_colours import green,cyan
from Bio import SeqIO
import csv

from squirrel.utils.config import *

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

