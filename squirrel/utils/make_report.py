import os
from Bio import SeqIO
import gzip
import collections
import yaml
import csv
from mako.lookup import TemplateLookup
import datetime as dt
from datetime import date

from mako.template import Template
from mako.runtime import Context
from mako.exceptions import RichTraceback
from io import StringIO

from squirrel import __version__
from squirrel.utils.log_colours import green,cyan
from squirrel.utils.config import *


def get_tree_svg(tree_image_file):
    svg = ""
    with open(tree_image_file,"r") as f:
        for l in f:
            l = l.rstrip("\n")
            svg+=f"{l}\n"
    return svg

config = {}

def make_output_report(report_to_generate,mask_file,config):
    #need to call this multiple times if there are multiple reports wanted
    
    data_for_report = {}
    if config[KEY_RUN_PHYLO]:
        tree_image_file = os.path.join(config[KEY_PHYLOGENY_SVG])
        data_for_report["phylo_svg_string"] = get_tree_svg(tree_image_file)
    else:
        data_for_report["phylo_svg_string"] = ""

    if config[KEY_SEQ_QC]:
        rows = []
        with open(mask_file,"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                rows.append(row)
        data_for_report["mask_csv"]=rows
    else:
        data_for_report["mask_csv"]=[]


    template_dir = os.path.abspath(os.path.dirname(config[KEY_REPORT_TEMPLATE]))
    mylookup = TemplateLookup(directories=[template_dir]) #absolute or relative works

    mytemplate = Template(filename=config[KEY_REPORT_TEMPLATE], lookup=mylookup)
    buf = StringIO()

    ctx = Context(buf, 
                    date = date.today(),
                    version = __version__,
                    data_for_report = data_for_report,
                    config=config)

    try:
        mytemplate.render_context(ctx)
    except:
        traceback = RichTraceback()
        for (filename, lineno, function, line) in traceback.traceback:
            print("File %s, line %s, in %s" % (filename, lineno, function))
            print(line, "\n")
        print("%s: %s" % (str(traceback.error.__class__.__name__), traceback.error))

    with open(report_to_generate, 'w') as fw:
        print(green("Generating: ") + f"{report_to_generate}")
        fw.write(buf.getvalue())


# make_output_report(output,config,data_for_report)