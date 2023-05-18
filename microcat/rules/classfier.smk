import glob
import os

# def get_barcodes(wildcards):
#     checkpoint_output = checkpoints.cellranger_unmapped_demultiplex.get(sample=wildcards.sample).output[0]
#     # Get all barcodes by parsing file names
#     barcodes = [os.path.basename(x).split("_")[1] for x in glob(os.path.join(checkpoint_output, "CB_*.bam"))]

#     # Construct fastq file paths
#     fastq_files = []
#     for barcode in barcodes:
#         fastq_files.extend([
#             os.path.join(
#                 config["output"]["host"],
#                 "cellranger_count",
#                 wildcards.sample,
#                 "unmapped_bam_CB_demultiplex",
#                 f"CB_{barcode}_R1.fastq"),
#             os.path.join(
#                 config["output"]["host"],
#                 "cellranger_count",
#                 wildcards.sample,
#                 "unmapped_bam_CB_demultiplex",
#                 f"CB_{barcode}_R2.fastq")
#         ])

#     return fastq_files

def get_CB_bam_files(wildcards):
    bam_dir = os.path.join(
        config["output"]["host"],
        "cellranger_count",
        wildcards.sample,
        "unmapped_bam_CB_demultiplex"
    )
    return glob.glob(os.path.join(bam_dir, "CB_*.bam"))


def aggregate_CB_bam_output(wildcards):
    demultiplex_output = checkpoints.cellranger_unmapped_demultiplex.get(**wildcards).output.unmapped_bam_CB_demultiplex_dir
    Barcode_list, = glob_wildcards(os.path.join(demultiplex_output, "CB_{barcode}.bam"))
    return expand(os.path.join(
        config["output"]["host"],
        "cellranger_count/{sample}//unmapped_bam_CB_demultiplex/CB_{barcode}.bam"),
        sample=wildcards.sample,
        barcode=Barcode_list)



# rule paired_bam_to_fastq:
#     input:
#         # expand("{sample_dir}/unmapped_bam_CB_demultiplex/CB_{barcode}.bam", 
#         #        sample_dir=os.path.join(config["output"]["host"], "cellranger_count/{sample}"), 
#         #        barcode=get_barcodes(wildcards.sample))
#         get_CB_bam_files
#     output:
#         # unmapped_fastq_1 = expand("{sample_dir}/unmapped_bam_CB_demultiplex/CB_{barcode}_R1.fastq", 
#         #                           sample_dir=os.path.join(config["output"]["host"], "cellranger_count/{sample}"), 
#         #                           barcode=get_barcodes(wildcards.sample)),
#         # unmapped_fastq_2 = expand("{sample_dir}/unmapped_bam_CB_demultiplex/CB_{barcode}_R2.fastq", 
#         #                           sample_dir=os.path.join(config["output"]["host"], "cellranger_count/{sample}"), 
#         #                           barcode=get_barcodes(wildcards.sample))
#         unmapped_fastq_1 = lambda wildcards: expand("{sample_dir}/unmapped_bam_CB_demultiplex/CB_{barcode}_R1.fastq", 
#                                                     sample_dir=os.path.join(config["output"]["host"], "cellranger_count/{sample}"), 
#                                                     barcode=get_barcodes(wildcards)),
#         unmapped_fastq_2 = lambda wildcards: expand("{sample_dir}/unmapped_bam_CB_demultiplex/CB_{barcode}_R2.fastq", 
#                                                     sample_dir=os.path.join(config["output"]["host"], "cellranger_count/{sample}"), 
#                                                     barcode=get_barcodes(wildcards))
#     threads:
#         16
#     priority: 11
#     shell:
#         '''
#         samtools fastq --threads {threads}   {input}  -1 {output.unmapped_fastq_1} -2  {output.unmapped_fastq_2}
#         '''


if config["params"]["classifier"]["kraken2uniq"]["do"]:
    rule kraken2uniq_cell_classified:
        input:
            r1 = expand(os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R1.fastq"), barcode=get_barcodes(wildcards.sample)),
            r2 = expand(os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R2.fastq"), barcode=get_barcodes(wildcards.sample))
        output:
            krak2_classified_output_r1 = expand(os.path.join(
                config["output"]["classifier"],
                "classified_output/{sample}/cell_level/{sample}_{barcode}_kraken2_classified_1.fq"), barcode=get_barcodes(wildcards.sample)),
            krak2_classified_output_r2 = expand(os.path.join(
                config["output"]["classifier"],
                "classified_output/{sample}/cell_level/{sample}_{barcode}_kraken2_classified_2.fq"), barcode=get_barcodes(wildcards.sample)),
            krak2_unclassified_output_r1 = expand(os.path.join(
                config["output"]["classifier"],
                "unclassified_output/{sample}/cell_level/{sample}_{barcode}_kraken2_unclassified_1.fq"), barcode=get_barcodes(wildcards.sample)),
            krak2_unclassified_output_r2 = expand(os.path.join(
                config["output"]["classifier"],
                "unclassified_output/{sample}/cell_level/{sample}_{barcode}_kraken2_unclassified_2.fq"), barcode=get_barcodes(wildcards.sample)),
            krak2_output = expand(os.path.join(
                config["output"]["classifier"],
                "kraken2_output/{sample}/cell_level/{sample}_{barcode}_kraken2_output.txt"), barcode=get_barcodes(wildcards.sample)),
            krak2_report = expand(os.path.join(
                config["output"]["classifier"],
                "kraken2_report/custom/{sample}/cell_level/{sample}_{barcode}_kraken2_report.txt"), barcode=get_barcodes(wildcards.sample)),
            krak2_std_report=expand(os.path.join(
                config["output"]["classifier"],
                "kraken2_report/standard/{sample}/cell_level/{sample}_{barcode}_kraken2_std_report.txt"), barcode=get_barcodes(wildcards.sample)),
            krak2_mpa_report=expand(os.path.join(
                config["output"]["classifier"],
                "kraken2_report/mpa/{sample}/cell_level/{sample}_{barcode}_kraken2_mpa_report.txt"), barcode=get_barcodes(wildcards.sample))
        params:
            database = config["params"]["classifier"]["kraken2uniq"]["kraken2_database"],
            #since kraken2uniq acquire specific input fomrat "#fq",so we put it it params
            krak2_classified_output=expand(os.path.join(
                config["output"]["classifier"],
                "classified_output/{sample}/cell_level/{sample}_{barcode}_kraken2_classified#.fq"), barcode=get_barcodes(wildcards.sample)),
            krak2_unclassified_output=expand(os.path.join(
                config["output"]["classifier"],
                "unclassified_output/{sample}/cell_level/{sample}_{barcode}_kraken2_unclassified#.fq"), barcode=get_barcodes(wildcards.sample)),
            threads=config["params"]["classifier"]["kraken2uniq"]["threads"]
        benchmark:
            expand(os.path.join(config["benchmarks"]["classifier"],
                                "kraken2uniq/{sample}/cell_level/{sample}_{barcode}_kraken2uniq_classifier_benchmark.log"), barcode=get_barcodes(wildcards.sample))
        log:
            expand(os.path.join(config["logs"]["classifier"],
                                "kraken2uniq/{sample}/cell_level/{sample}_{barcode}_kraken2uniq_classifier.log"), barcode=get_barcodes(wildcards.sample))
        conda:
            config["envs"]["kraken2"]
        shell:
            '''
            kraken2 --db {params.database} \
            --threads {params.threads} \
            --classified-out {params.krak2_classified_output}\
            --unclassified-out {params.krak2_unclassified_output}\
            --output {output.krak2_output} \
            --report {output.krak2_report} \
            --report-minimizer-data \
            {input.r1} {input.r2} \
            --use-names \
            --paired \
            --memory-mapping \
            2> >(tee {log})
            cut -f 1-3,6-8 {output.krak2_report} > {output.krak2_std_report}
            python /data/scRNA/test-10x/scripts/kraken2mpa.py -r {output.krak2_std_report} -o {output.krak2_mpa_report}
            '''
    # rule kraken2uniq_classified:
    #     input:
    #         r1 = os.path.join(
    #             config["output"]["host"],
    #             "starsolo_count/unmapped_fastq/{sample}_1.fastq"),
    #         r2 = os.path.join(
    #             config["output"]["host"],
    #             "starsolo_count/unmapped_fastq/{sample}_2.fastq")
    #     output:
    #         krak2_classified_output_r1 = os.path.join(
    #             config["output"]["classifier"],
    #             "classified_output/{sample}/{sample}_kraken2_classified_1.fq"),
    #         krak2_classified_output_r2 = os.path.join(
    #             config["output"]["classifier"],
    #             "classified_output/{sample}/{sample}_kraken2_classified_2.fq"),
    #         krak2_unclassified_output_r1 = os.path.join(
    #             config["output"]["classifier"],
    #             "unclassified_output/{sample}/{sample}_kraken2_unclassified_1.fq"),
    #         krak2_unclassified_output_r2 = os.path.join(
    #             config["output"]["classifier"],
    #             "unclassified_output/{sample}/{sample}_kraken2_unclassified_2.fq"),
    #         krak2_output = os.path.join(
    #             config["output"]["classifier"],
    #             "kraken2_output/{sample}/{sample}_kraken2_output.txt"),
    #         krak2_report = os.path.join(
    #             config["output"]["classifier"],
    #             "kraken2_report/custom/{sample}/{sample}_kraken2_report.txt"),
    #         krak2_std_report=os.path.join(
    #             config["output"]["classifier"],
    #             "kraken2_report/standard/{sample}/{sample}_kraken2_std_report.txt"),
    #         krak2_mpa_report=os.path.join(
    #             config["output"]["classifier"],
    #             "kraken2_report/mpa/{sample}/{sample}_kraken2_mpa_report.txt")
    #     params:
    #         database = config["params"]["classifier"]["kraken2uniq"]["kraken2_database"],
    #         #since kraken2uniq acquire specific input fomrat "#fq",so we put it it params
    #         krak2_classified_output=os.path.join(
    #             config["output"]["classifier"],
    #             "classified_output/{sample}/{sample}_kraken2_classified#.fq"),
    #         krak2_unclassified_output=os.path.join(
    #             config["output"]["classifier"],
    #             "unclassified_output/{sample}/{sample}_kraken2_unclassified#.fq")
    #     threads: 
    #         config["params"]["classifier"]["kraken2uniq"]["threads"]
    #     benchmark:
    #         os.path.join(config["benchmarks"]["classifier"],
    #                     "kraken2uniq/{sample}_kraken2uniq_classifier_benchmark.log")
    #     log:
    #         os.path.join(config["logs"]["classifier"],
    #                     "kraken2uniq/{sample}_kraken2uniq_classifier.log")
    #     message:
    #         "Classifier: Performing Taxonomic Classifcation of Sample {sample} with kraken2uniq."
    #     conda:
    #         config["envs"]["kraken2"]
    #     shell:
    #         '''
    #         kraken2 --db {params.database} \
    #         --threads 8 \
    #         --classified-out {params.krak2_classified_output}\
    #         --unclassified-out {params.krak2_unclassified_output}\
    #         --output {output.krak2_output} \
    #         --report {output.krak2_report} \
    #         --report-minimizer-data \
    #         {input.r1} {input.r2} \
    #         --use-names \
    #         --paired \
    #         --memory-mapping \
    #         2> >(tee {log})
    #         cut -f 1-3,6-8 {output.krak2_report} > {output.krak2_std_report}
    #         python /data/scRNA/test-10x/scripts/kraken2mpa.py -r {output.krak2_std_report} -o {output.krak2_mpa_report}
    #         '''
else:
    rule kraken2uniq_classified:
        input:    

if config["params"]["classifier"]["krakenuniq"]["do"]:
    # rule krakenuniq_classifier:
    #     input:
    #         r1 = os.path.join(
    #             config["output"]["host"],
    #             "starsolo_count/unmapped_fastq/{sample}_1.fastq"),
    #         r2 = os.path.join(
    #             config["output"]["host"],
    #             "starsolo_count/unmapped_fastq/{sample}_2.fastq")
    #     output:
    #         krakenuniq_output = os.path.join(
    #             config["output"]["classifier"],
    #             "krakenuniq_output/{sample}/{sample}_krakenuniq_output.txt"),
    #         krakenuniq_report = os.path.join(
    #             config["output"]["classifier"],
    #             "krakenuniq_report/custom/{sample}/{sample}_krakenuniq_report.txt"),
    #         krakenuniq_classified_output = os.path.join(
    #             config["output"]["classifier"],
    #             "krakenuniq_classified_output/{sample}/{sample}_krakenuniq_classified.fq"),
    #         krakenuniq_unclassified_output = os.path.join(
    #             config["output"]["classifier"],
    #             "krakenuniq_classified_output/{sample}/{sample}_krakenuniq_classified.fq")
    #     params:
    #         database = config["params"]["classifier"]["krakenuniq"]["krakenuniq_database"],
    #         threads=config["params"]["classifier"]["krakenuniq"]["threads"],
    #         estimate_precision=config["params"]["classifier"]["krakenuniq"]["estimate_precision"]
    #     conda:
    #         config["envs"]["krakenuniq"]
    #     message:
    #         "Classifier: Performing Taxonomic Classifcation of Sample {sample} with krakenuniq."
    #     shell:
    #         '''
    #         krakenuniq --db {params.database} \
    #         --threads {params.threads} \
    #         --hll-precision {params.estimate_precision} \
    #         --classified-out {params.krakenuniq_classified_output}\
    #         --unclassified-out {params.krakenuniq_unclassified_output}\
    #         --output {output.krakenuniq_output} \
    #         --report-file {output.krakenuniq_report} \
    #         {input.r1} {input.r2} \
    #         --paired \
    #         --preload \
    #         --check-names \
    #         2>&1 | tee {log})
    #         '''
    rule krakenuniq_cell_level_classifier:
        input:
            r1 = expand(os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R1.fastq"), barcode=get_barcodes(wildcards.sample)),
            r2 = expand(os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R2.fastq"), barcode=get_barcodes(wildcards.sample))
        output:
            krakenuniq_output = expand(os.path.join(
                config["output"]["classifier"],
                "krakenuniq_output/{sample}/cell_level/{sample}_{barcode}_krakenuniq_output.txt"), barcode=get_barcodes(wildcards.sample)),
            krakenuniq_report = expand(os.path.join(
                config["output"]["classifier"],
                "krakenuniq_report/custom/{sample}/cell_level/{sample}_{barcode}_krakenuniq_report.txt"), barcode=get_barcodes(wildcards.sample)),
            krakenuniq_classified_output = expand(os.path.join(
                config["output"]["classifier"],
                "krakenuniq_classified_output/{sample}/cell_level/{sample}_{barcode}_krakenuniq_classified.fq"), barcode=get_barcodes(wildcards.sample)),
            krakenuniq_unclassified_output = expand(os.path.join(
                config["output"]["classifier"],
                "krakenuniq_classified_output/{sample}/cell_level/{sample}_{barcode}_krakenuniq_unclassified.fq"), barcode=get_barcodes(wildcards.sample))
        params:
            database = config["params"]["classifier"]["krakenuniq"]["krakenuniq_database"],
            threads=config["params"]["classifier"]["krakenuniq"]["threads"],
            estimate_precision=config["params"]["classifier"]["krakenuniq"]["estimate_precision"]
        benchmark:
            expand(os.path.join(config["benchmarks"]["classifier"],
                        "krakenuniq/{sample}/cell_level/{sample}_{barcode}_krakenuniq_classifier_benchmark.log"), barcode=get_barcodes(wildcards.sample))
        log:
            expand(os.path.join(config["logs"]["classifier"],
                        "krakenuniq/{sample}/cell_level/{sample}_{barcode}_krakenuniq_classifier.log"), barcode=get_barcodes(wildcards.sample))
        conda:
            config["envs"]["krakenuniq"]
        shell:
            '''
            krakenuniq --db {params.database} \
            --threads {params.threads} \
            --hll-precision {params.estimate_precision} \
            --classified-out {params.krakenuniq_classified_output}\
            --unclassified-out {params.krakenuniq_unclassified_output}\
            --output {output.krakenuniq_output} \
            --report-file {output.krakenuniq_report} \
            {input.r1} {input.r2} \
            --paired \
            --preload \
            --check-names \
            2>&1 | tee {log})
            '''

else:
    rule krakenuniq_classifier:
        input:    

# rule classifier_all:
#     input:
#         rules.kraken2uniq_classifier.input,
#         rules.krakenuniq_classifier.input,
rule classifier_all:
    input:
        expand([
            os.path.join(
            config["output"]["host"],
            "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R1.fastq"),
            os.path.join(
            config["output"]["host"],
            "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R2.fastq")],
                sample=glob_wildcards(os.path.join(config["output"]["host"], "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}.bam")).sample,
                barcode=glob_wildcards(os.path.join(config["output"]["host"], "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}.bam")).barcode
                )