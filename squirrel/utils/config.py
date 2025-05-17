#!/usr/bin/env python3

KEY_INPUT_FASTA="fasta"
KEY_REFERENCE_FASTA = "reference_fasta"
KEY_REFERENCE_PANEL = "reference_panel"
KEY_BACKGROUND_FASTA = "background_fasta"
KEY_TO_MASK = "to_mask"
KEY_GENE_BOUNDARIES = "gene_boundaries"
KEY_OUTDIR = "outdir"
KEY_OUTFILE = "alignment_file"
KEY_EXCLUDE_FILE = "exclude_file"
KEY_OUTFILE_STEM = "outfile_stem"
KEY_OUTFILENAME = "outfilename"
KEY_TEMPDIR = "tempdir"
KEY_CDS_OUTFILE = "cds_outfile"
KEY_CDS_OUTFILENAME = "cds_outfilename"
KEY_CONCATENATE = "concatenate"
KEY_NO_TEMP="no_temp"
KEY_VERBOSE="verbose"
KEY_THREADS = "threads"
KEY_PHYLO_THREADS = "phylo_threads"
KEY_NO_MASK="no_mask"
KEY_ADDITIONAL_MASK="additional_mask"
KEY_SEQUENCE_MASK="sequence_mask"
KEY_TRIM_END="trim_end"
KEY_EXTRACT_CDS="extract_cds"
KEY_SEQ_QC = "seq_qc"
KEY_ASSEMBLY_REFERENCES = "assembly_references"

KEY_CLADE = "clade"
KEY_ASSIGNED_CLADES = "assigned_clades"
KEY_SPLIT_CLADE="split_clade"
KEY_APPEND_CLADE_STR = "append_clade_str"
KEY_RUN_PHYLO="run_phylo"
KEY_RUN_APOBEC3_PHYLO = "run_apobec3_phylo"
KEY_OUTGROUPS="outgroups"
KEY_PHYLOGENY="phylogeny"
KEY_PHYLOGENY_SVG="phylogeny_svg"
KEY_PHYLOGENY_INTERACTIVE="phylogeny_interactive"
KEY_INCLUDE_BACKGROUND = "include_background"
KEY_BACKGROUND_FILE = "background_file"
KEY_OUTGROUP_STRING="outgroup_string"
KEY_OUTGROUP_SENTENCE="outgroup_sentence"
KEY_GRANTHAM_SCORES="grantham_scores"
KEY_REPORT_TEMPLATE = "report_template"
KEY_INTERACTIVE_SCRIPT = "interactive_script"

KEY_TREE = "tree"
KEY_BRANCH_RECONSTRUCTION = "branch_reconstruction"
KEY_ASR_TREE = "asr_tree"
KEY_ASR_STATE = "asr_state"
KEY_ASR_ALIGNMENT = "asr_alignment"

KEY_FIG_HEIGHT = "fig_height"
KEY_FIG_WIDTH = "fig_width"
KEY_POINT_STYLE = "point_style"
KEY_POINT_JUSTIFY = "point_justify"

VALUE_TRIM_END = 190788
VALUE_VALID_CLADES = ["cladei","cladeia","cladeib","cladeii","cladeiia","cladeiib","variola","split"]
OUTGROUP_DICT = {
    "variola":["KJ642615|human|Nigeria||1978"],
    "cladei":["KJ642615|human|Nigeria||1978"],
    "cladeia":["KJ642615|human|Nigeria||1978"],
    "cladeib":["KJ642613|human|DRC|Equateur|1970-09-01"],
    "cladeii":["KJ642613|human|DRC|Equateur|1970-09-01"],
    "cladeiia":["KJ642613|human|DRC|Equateur|1970-09-01"],
    "cladeiib":["KJ642615|human|Nigeria||1978"],
    "split":["KJ642615|human|Nigeria||1978"]
}

VALUE_OUTFILE_STEM = "sequences"
VALUE_EXCLUDE_FILE_STEM = "suggested_to_exclude"