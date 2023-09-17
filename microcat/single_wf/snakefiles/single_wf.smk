#!/usr/bin/env snakemake

import sys
from snakemake.utils import min_version
import os
import numpy as np
import pandas

import microcat
MICROCAT_DIR = microcat.__path__[0]

wildcard_constraints:
    # Patient = "[a-zA-Z0-9_]+", # Any alphanumeric characters and underscore
    # tissue = "S[0-9]+",  # S followed by any number
    lane = "L[0-9]{3}",  # L followed by exactly 3 numbers
    plate = "P[0-9]{3}",  # L followed by exactly 3 numbers
    library = "[0-9]{3}"  # Exactly 3 numbers


min_version("7.0")
shell.executable("bash")

class ansitxt:
    RED = '\033[31m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def warning(msg):
    print(f"\n{ansitxt.BOLD}{ansitxt.RED}{msg}{ansitxt.ENDC}\n",file=sys.stderr)


PLATFORM = None

if config["params"]["host"]["starsolo"]["do"]:
    if  config["params"]["host"]["starsolo"]["soloType"]=="CB_UMI_Simple":
        PLATFORM = "lane"
    elif config["params"]["host"]["starsolo"]["soloType"]=="SmartSeq":
        PLATFORM = "plate"
    elif config["params"]["host"]["starsolo"]["soloType"]=="CB_UMI_Complex":
        PLATFORM = "lane"
    else:
        raise ValueError("Platform must be either 'CB_UMI' or 'smartseq'")
elif config["params"]["host"]["cellranger"]["do"]:
    PLATFORM = "lane"
else:
    raise ValueError("Platform must be either 'lane' or 'plate'")

if config["params"]["begin"] == "host":
    try:
        SAMPLES = microcat.parse_samples(config["params"]["samples"],platform = PLATFORM)
        SAMPLES_ID_LIST = SAMPLES.index.get_level_values("sample_id").unique()
    except FileNotFoundError:
        warning(f"ERROR: the samples file does not exist. Please see the README file for details. Quitting now.")
        sys.exit(1)
    include: "../rules/host.smk"
    # include: "../rules/ERCC.smk"
    include: "../rules/classfier.smk"

elif config["params"]["begin"] == "classifier":
    try:
        SAMPLES = microcat.parse_bam_samples(config["params"]["samples"],platform = PLATFORM)
        SAMPLES_ID_LIST = SAMPLES.index.get_level_values("sample_id").unique()
    except FileNotFoundError:
        warning(f"ERROR: the samples file does not exist. Please see the README file for details. Quitting now.")
        sys.exit(1)
    
    include: "../rules/bam_host.smk"
    # include: "../rules/ERCC.smk"
    include: "../rules/classfier.smk"

rule all:
    input:
        rules.host_all.input,
        rules.classifier_all.input