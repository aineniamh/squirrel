#!/usr/bin/env python3
import os
import csv
import sys
import datetime as dt

import snakemake
from squirrel.utils.config import *
from squirrel.utils.log_colours import green,cyan

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