#!/usr/bin/env python3

KEY_INPUT_FASTA="fasta"
KEY_REFERENCE_FASTA = "reference_fasta"
KEY_BACKGROUND_FASTA = "background_fasta"
KEY_TO_MASK = "to_mask"
KEY_GENE_BOUNDARIES = "gene_boundaries"
KEY_OUTDIR = "outdir"
KEY_OUTFILE = "alignment_file"
KEY_OUTFILE_STEM = "outfile_stem"
KEY_OUTFILENAME = "outfilename"
KEY_TEMPDIR = "tempdir"
KEY_CDS_OUTFILE = "cds_outfile"
KEY_CONCATENATE = "concatenate"
KEY_NO_TEMP="no_temp"
KEY_VERBOSE="verbose"
KEY_THREADS = "threads"
KEY_PHYLO_THREADS = "phylo_threads"
KEY_NO_MASK="no_mask"
KEY_ADDITIONAL_MASK="additional_mask"
KEY_TRIM_END="trim_end"
KEY_EXTRACT_CDS="extract_cds"
KEY_SEQ_QC = "seq_qc"
KEY_ASSEMBLY_REFERENCES = "assembly_references"

KEY_CLADE = "clade"
KEY_RUN_PHYLO="run_phylo"
KEY_RUN_APOBEC3_PHYLO = "run_apobec3_phylo"
KEY_OUTGROUPS="outgroups"
KEY_PHYLOGENY="phylogeny"
KEY_PHYLOGENY_SVG="phylogeny_svg"
KEY_INCLUDE_BACKGROUND = "include_background"
KEY_BACKGROUND_FILE = "background_file"
KEY_OUTGROUP_STRING="outgroup_string"
KEY_OUTGROUP_SENTENCE="outgroup_sentence"
KEY_GRANTHAM_SCORES="grantham_scores"
KEY_REPORT_TEMPLATE = "report_template"

KEY_FIG_HEIGHT = "fig_height"

VALUE_VALID_CLADES = ["cladei","cladeia","cladeib","cladeii","cladeiia","cladeiib"]
OUTGROUP_DICT = {
    "cladei":["KJ642615|human|Nigeria||1978"],
    "cladeia":["KJ642615|human|Nigeria||1978"],
    "cladeib":["KJ642613|human|DRC|Equateur|1970-09-01"],
    "cladeii":["KJ642613|human|DRC|Equateur|1970-09-01"],
    "cladeiia":["KJ642613|human|DRC|Equateur|1970-09-01"],
    "cladeiib":["KJ642615|human|Nigeria||1978"]
}