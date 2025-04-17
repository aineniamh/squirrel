#!/usr/bin/env python3
import os
import csv
import sys
import datetime as dt

import snakemake
from squirrel.utils.config import *
from squirrel.utils.log_colours import green,cyan,red
import yaml

import squirrel.utils.custom_logger as custom_logger

def run_snakemake(snake_config,snakefile,verbose,config):
    if verbose:
        print(red("\n**** CONFIG ****"))
        for k in sorted(config):
            print(green(f" - {k}: ") + f"{config[k]}")
        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                    workdir=config[KEY_TEMPDIR], config=snake_config, cores=config[KEY_THREADS],lock=False
                                    )
    else:
        logger = custom_logger.Logger()
        status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True, force_incomplete=True,
                                    workdir=config[KEY_TEMPDIR], config=snake_config, cores=config[KEY_THREADS],lock=False,
                                    quiet=True,log_handler=logger.log_handler
                                    )
    return status


def load_yaml(yamlfile):
    input_config = ""
    with open(yamlfile, "r") as f:
        try:
            input_config = yaml.load(f, Loader=yaml.FullLoader)
        except:
            sys.stderr.write(cyan(f'Error: failed to read config file. Ensure your file in correct yaml format.\n'))
            sys.exit(-1)
    return input_config