#!/usr/bin/env python3
import os
import csv
import sys
import datetime as dt
from pathlib import Path

from squirrel.utils.config import *
from squirrel.utils.log_colours import green,cyan,red
import yaml


from pathlib import Path

from snakemake.api import (
    OutputSettings,
    ResourceSettings,
    WorkflowSettings,
    SnakemakeApi,
    StorageSettings,
    ConfigSettings,
    DAGSettings,
    ExecutionSettings
)

import squirrel.utils.custom_logger as custom_logger

# def run_snakemake(snake_config,snakefile,verbose,config):
#     if verbose:
#         print(red("\n**** CONFIG ****"))
#         for k in sorted(config):
#             print(green(f" - {k}: ") + f"{config[k]}")

#         status = api.SnakemakeApi(snakefile, printshellcmds=True, 
#         forceall=True, f
#         orce_incomplete=True,
#                                     workdir=config[KEY_TEMPDIR], 
#                                     config=snake_config, cores=config[KEY_THREADS],lock=False
#                                     )
#     else:
#         logger = custom_logger.Logger()
#         status = api.SnakemakeApi(snakefile, printshellcmds=False, forceall=True, force_incomplete=True,
#                                     workdir=config[KEY_TEMPDIR], config=snake_config, cores=config[KEY_THREADS],lock=False,
#                                     quiet=True,log_handler=logger.log_handler
#                                     )
#     return status


def run_snakemake(snake_config,my_snakefile,verbose,config):
    with SnakemakeApi(
        OutputSettings(
            printshellcmds=True,
            quiet=False,
            verbose=True,
            log_handler_settings=custom_logger,
        )
    ) as snakemake_api:
        try:
        
            workflow_api = snakemake_api.workflow(
                resource_settings=ResourceSettings(
                    cores=config[KEY_THREADS]
                    ),
                config_settings=ConfigSettings(
                    config=config
                ),
                snakefile=Path(my_snakefile),
                workdir=Path(config[KEY_TEMPDIR]),
                )
            dag_api = workflow_api.dag(
                dag_settings=DAGSettings(
                    forceall=True,
                    force_incomplete=True
                ),
                    )
            dag_api.execute_workflow(
                execution_settings=ExecutionSettings(
                    lock=False,
                    keep_incomplete=True
                ),
            )


        except Exception as e:
            snakemake_api.print_exception(e)
            return False
    
    return True





def load_yaml(yamlfile):
    input_config = ""
    with open(yamlfile, "r") as f:
        try:
            input_config = yaml.load(f, Loader=yaml.FullLoader)
        except:
            sys.stderr.write(cyan(f'Error: failed to read config file. Ensure your file in correct yaml format.\n'))
            sys.exit(-1)
    return input_config