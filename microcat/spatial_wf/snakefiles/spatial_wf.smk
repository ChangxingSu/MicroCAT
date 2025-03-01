#!/usr/bin/env snakemake

import sys
from snakemake.utils import min_version
import os
import numpy as np
import pandas

import microcat
MICROCAT_DIR = microcat.__path__[0]

wildcard_constraints:
    lane = "L[0-9]{3}",  # L followed by exactly 3 numbers
    section = "SEC[0-9]{3}", # SEC followed by exactly 3 numbers
    library = "[0-9]{3}"  # Exactly 3 numbers

min_version("7.0")
shell.executable("bash")

class ansitxt:
    RED = '\033[31m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def warning(msg):
    print(f"\n{ansitxt.BOLD}{ansitxt.RED}{msg}{ansitxt.ENDC}\n",file=sys.stderr)

if config["params"]["begin"] == "host":
    try:
        SAMPLES = microcat.parse_spatial_samples(config["params"]["samples"])
        SAMPLES_ID_LIST = SAMPLES.index.get_level_values("sample_id").unique()
    except FileNotFoundError as e:
        if "File not found" in str(e):
            warning(f"ERROR: {e}. Please see the README file for details.")
            sys.exit(1)

    include: "../rules/visium.smk"
    include: "../rules/krakenuniq.smk"
rule all:
    input:
        rules.spaceranger_all.input,
        rules.krakuniq_all.input