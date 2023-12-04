import os
from squirrel.utils.config import *
from Bio import SeqIO
from Bio.Seq import Seq
import csv
from squirrel.utils.log_colours import green,cyan


rule all:
    input:
        os.path.join(config[KEY_OUTDIR],config[KEY_PHYLOGENY])



