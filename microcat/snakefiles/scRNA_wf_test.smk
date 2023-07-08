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
    if "tenX" in config["params"]["host"]["starsolo"]["assay"]:
        PLATFORM = "tenX"
    elif "Smartseq" in config["params"]["host"]["starsolo"]["assay"]:
        PLATFORM = "smartseq"
    elif "Seq-well" in config["params"]["host"]["starsolo"]["assay"]:
        PLATFORM = "tenX"
    else:
        raise ValueError("Platform must be either 'tenX' or 'smartseq'")
elif config["params"]["host"]["cellranger"]["do"]:
    PLATFORM = "tenX"
else:
    raise ValueError("Platform must be either 'tenX' or 'smartseq'")




try:
    SAMPLES = microcat.parse_samples(config["params"]["samples"],platform = PLATFORM)
    SAMPLES_ID_LIST = SAMPLES.index.get_level_values("sample_id").unique()
except FileNotFoundError:
    warning(f"ERROR: the samples file does not exist. Please see the README file for details. Quitting now.")
    sys.exit(1)

def aggregate_input(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the scatter step
    '''
    checkpoint_output = checkpoints.split_starsolo_BAM_by_RG.get(**wildcards).output[0]
    return expand('unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam',
        sample=SAMPLES_ID_LIST)

rule generate_pe_manifest_file:
    input:
        config["params"]["samples"],
    output:
        PE_MANIFEST_FILE = os.path.join("data", "manifest.tsv"),
    script:
        "../scripts/generate_PE_manifest_file.py"            



rule starsolo_smartseq_count:
    # Input files
    input:
        # Path to the input manifest file
        manifest = "/home/microcat-sucx/project/microcat_smartseq/data/manifest.tsv",
        # fastqs_dir=lambda wildcards: microcat.get_fastqs_dir(SAMPLES,wildcards),
    output:
        # Path to the output features.tsv file
        features_file = os.path.join(
            config["output"]["host"],
            "starsolo_count/features.tsv"),
        # Path to the output matrix.mtx file
        matrix_file = os.path.join(
            config["output"]["host"],
            "starsolo_count/matrix.mtx"),
        # Path to the output barcodes.tsv file
        barcodes_file = os.path.join(
            config["output"]["host"],
            "starsolo_count/barcodes.tsv"),
        mapped_bam_file = os.path.join(
            config["output"]["host"],
            "starsolo_count/Aligned_out.bam")
    params:
        # Path to the output directory
        starsolo_out = os.path.join(
            config["output"]["host"],
            "starsolo_count/"),
        # Path to the STAR index directory
        reference = config["params"]["host"]["starsolo"]["reference"],
        SAMattrRGline = microcat.get_SAMattrRGline_from_manifest(os.path.join("data", "manifest.tsv")),
        # Additional parameters for STAR
        variousParams = config["params"]["host"]["starsolo"]["variousParams"],
        # Number of threads for STAR
        threads = config["params"]["host"]["starsolo"]["threads"]
    log:
        os.path.join(config["logs"]["host"],
                    "starsolo/starsolo_count_smartseq2.log")
    benchmark:
        os.path.join(config["benchmarks"]["host"],
                    "starsolo/starsolo_count_smartseq2.benchmark")
    conda:
        config["envs"]["star"]
    shell:
        '''
        mkdir -p {params.starsolo_out}; 
        cd {params.starsolo_out} ;
        STAR \
        --soloType SmartSeq \
        --genomeDir {params.reference} \
        --readFilesManifest {input.manifest} \
        --runThreadN {params.threads} \
        --soloUMIdedup Exact \
        --soloStrand Unstranded \
        --outSAMtype BAM Unsorted\
        --outSAMattrRGline {params.SAMattrRGline} \
        --readFilesCommand zcat \
        --outSAMunmapped Within \
        --quantMode GeneCounts \
        {params.variousParams}  \
        2>&1 | tee ../../../{log} ;
        pwd ;\
        cd ../../../;\
        ln -sr "{params.starsolo_out}/Solo.out/Gene/filtered/features.tsv" "{output.features_file}" ;\
        ln -sr "{params.starsolo_out}/Solo.out/Gene/filtered/matrix.mtx" "{output.matrix_file}" ; \
        ln -sr "{params.starsolo_out}/Solo.out/Gene/filtered/barcodes.tsv" "{output.barcodes_file}" ;\
        mv "{params.starsolo_out}/Aligned.out.bam" "{output.mapped_bam_file}";\
        '''


rule starsolo_smartseq_extracted:
    input:
        mapped_bam_file = os.path.join(
            config["output"]["host"],
            "starsolo_count/Aligned_out.bam"),
    output:
        unmapped_bam_unsorted_file = os.path.join(
            config["output"]["host"],
            "starsolo_count/Aligned_out_unmapped.bam"),
    params:
        threads=16
    conda:
        config["envs"]["star"]
    shell:
        '''
        samtools view --threads  {params.threads}  -b -f 4   {input.mapped_bam_file}  >  {output.unmapped_bam_unsorted_file}
        '''

rule starsolo_smartseq_unmapped_sorted_bam:
    input:
        unmapped_bam_unsorted_file = os.path.join(
            config["output"]["host"],
            "starsolo_count/Aligned_out_unmapped.bam")
    output:
        unmapped_sorted_bam = os.path.join(
            config["output"]["host"],
            "starsolo_count/Aligned_out_unmapped_RGsorted.bam"),
    params:
        threads=40,
        tag="RG"
    log:
        os.path.join(config["logs"]["host"],
                    "starsolo/unmapped_sorted_bam.log")
    benchmark:
        os.path.join(config["benchmarks"]["host"],
                    "starsolo/unmapped_sorted_bam.benchmark")
    conda:
        config["envs"]["star"]
    shell:
        '''
        samtools sort -@ {params.threads} -t {params.tag} -o {output.unmapped_sorted_bam}  {input.unmapped_bam_unsorted_file};
        '''
rule split_starsolo_BAM_by_RG:
    input:
        unmapped_sorted_bam = os.path.join(
                config["output"]["host"],
                "starsolo_count/Aligned_out_unmapped_RGsorted.bam"),
    output:
        unmapped_bam_sorted_file =os.path.join(
            config["output"]["host"],
            "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
    # params:
    #     SampleID="{sample}",
    # shell:
    #     '''
    #     python /data/project/host-microbiome/microcat/microcat/scripts/split_Starsolo_BAM_by_RG.py \
    #     --bam_path {input.unmapped_sorted_bam} \
    #     --tag {params.SampleID} \
    #     --output_bam {output.unmapped_bam_sorted_file} 
    #     '''
    script:
        "../src/split_PathSeq_BAM_by_CB_UB.py"


rule all:
    input:
        os.path.join(
            config["output"]["host"],
            "starsolo_count/Aligned_out.bam")
        # expand(os.path.join(
        #     config["output"]["host"],
        #     "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),sample=SAMPLES_ID_LIST)
# rule all:
#     input:
#         os.path.join(
#             config["output"]["host"],
#             "starsolo_count/Aligned_out_unmapped.bam")


# rule all:
#     input:
#         # expand(os.path.join(config["output"]["host"],"unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),sample=SAMPLES_ID_LIST)
#         os.path.join(
#             config["output"]["host"],
#             "starsolo_count/Aligned_out_unmapped_RGsorted.bam")
        # os.path.join(
        #     config["output"]["host"],
        #     "starsolo_count/Aligned_out_unmapped.bam")
        # os.path.join(
        #         config["output"]["host"],
        #         "starsolo_count/features.tsv"),
        # # Path to the output matrix.mtx file
        # os.path.join(
        #         config["output"]["host"],
        #         "starsolo_count/matrix.mtx"),
        # # Path to the output barcodes.tsv file
        # os.path.join(
        #         config["output"]["host"],
        #         "starsolo_count/barcodes.tsv"),
        # os.path.join(
        #         config["output"]["host"],
        #         "starsolo_count/Aligned_out.bam")
        # os.path.join("data", "manifest.tsv"),
#         rules.classifier_all.input
