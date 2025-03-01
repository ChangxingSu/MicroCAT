
KUNIQ_OUTPUT = os.path.join(
            config["output"]["classifier"],
            "rmhost_krakenuniq_output/{sample}/{sample}_krakenuniq_output.txt")
KUNIQ_CLASSIFIED_OUTPUT = os.path.join(
            config["output"]["classifier"],
            "rmhost_krakenuniq_classified_output/{sample}/{sample}_krakenuniq_classified.fq"),
KUNIQ_REPORT = os.path.join(
            config["output"]["classifier"],
            "rmhost_krakenuniq_report/{sample}/{sample}_krakenuniq_report.txt")
KUNIQ_TAXDB = os.path.join(
            config["params"]["classifier"]["krakenuniq"]["krakenuniq_database"],
            "taxDB")
KUNIQ_GENUS_COUNT = os.path.join(
            config["output"]["classifier"],
            "rmhost_krakenuniq_report/{sample}/counts/{sample}_genus.tsv")
KUNIQ_SPECIES_COUNT = os.path.join(
            config["output"]["classifier"],
            "rmhost_krakenuniq_report/{sample}/counts/{sample}_species.tsv")
KUNIQ_GENUS_TAXA_COUNT = os.path.join(
            config["output"]["classifier"],
            "rmhost_krakenuniq_report/{sample}/counts/{sample}_taxa_genus.tsv")
KUNIQ_SPECIES_TAXA_COUNT = os.path.join(
            config["output"]["classifier"],
            "rmhost_krakenuniq_report/{sample}/counts/{sample}_taxa_species.tsv")      

rule paired_bam_to_fastq:
    input:
        unmapped_bam_sorted_file =os.path.join(
        config["output"]["host"],
        "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
    output:
        unmapped_fastq = os.path.join(
            config["output"]["host"],
            "unmapped_host/{sample}/{sample}_unmappped2human_bam.fastq"),
        unmapped_r1_fastq = os.path.join(
            config["output"]["host"],
            "unmapped_host/{sample}/{sample}_unmappped2human_bam_r1.fastq"),
        unmapped_r2_fastq = os.path.join(
            config["output"]["host"],
            "unmapped_host/{sample}/{sample}_unmappped2human_bam_r2.fastq")
    # log:
    #     os.path.join(config["logs"]["host"],
    #                 "bam2fastq/{sample}_bam_convert_fastq.log")
    params:
        bam2fastq_script = config["scripts"]["bam2fastq"],
    threads:
        config["resources"]["paired_bam_to_fastq"]["threads"]
    resources:
        mem_mb=config["resources"]["paired_bam_to_fastq"]["mem_mb"]
    priority: 11
    conda:
        "MicroCAT"
    shell:
        '''
        bash {params.bam2fastq_script} {input.unmapped_bam_sorted_file} {output.unmapped_r1_fastq} {output.unmapped_r2_fastq} {output.unmapped_fastq} {threads}
        '''

rule krakenuniq_simple_classifier:
    input:
        unmapped_fastq = os.path.join(
            config["output"]["host"],
            "unmapped_host/{sample}/{sample}_unmappped2human_bam.fastq"),
        unmapped_r1_fastq = os.path.join(
            config["output"]["host"],
            "unmapped_host/{sample}/{sample}_unmappped2human_bam_r1.fastq"),
        unmapped_r2_fastq = os.path.join(
            config["output"]["host"],
            "unmapped_host/{sample}/{sample}_unmappped2human_bam_r2.fastq")
    output:
        krakenuniq_output = KUNIQ_OUTPUT,
        krakenuniq_report = KUNIQ_REPORT,
    params:
        database = config["params"]["classifier"]["krakenuniq"]["krakenuniq_database"],
        estimate_precision=config["params"]["classifier"]["krakenuniq"]["estimate_precision"],
        variousParams = config["params"]["classifier"]["krakenuniq"]["variousParams"],
    benchmark:
        os.path.join(config["benchmarks"]["classifier"],
                    "rmhost_krakenuniq/{sample}_kraken2uniq_classifier_benchmark.tsv")
    log:
        os.path.join(config["logs"]["classifier"],
                    "rmhost_krakenuniq/{sample}_krakenuniq_classifier.log")
    threads:
        40
    resources:
        mem_mb= 140000,
    conda:
        "krakuniq"
    # message:
    #     "classifier: Performing Taxonomic Classifcation of Sample {sample} with krakenuniq."
    shell:
        '''
        if [ -s "{input.unmapped_fastq}" ]; then
            krakenuniq --db {params.database} \
            --threads {threads} \
            --hll-precision {params.estimate_precision} \
            --output {output.krakenuniq_output} \
            --report-file {output.krakenuniq_report} \
            {input.unmapped_fastq}  \
            --preload-size 140G \
            {params.variousParams} \
            2>&1 | tee {log}
        else
            krakenuniq --db {params.database} \
            --threads {threads} \
            --hll-precision {params.estimate_precision} \
            --output {output.krakenuniq_output} \
            --report-file {output.krakenuniq_report} \
            {input.unmapped_r1_fastq} {input.unmapped_r2_fastq}\
            --preload-size 140G \
            --paired \
            {params.variousParams} \
            2>&1 | tee {log}
        fi
        '''


rule krakuniq2sc:
    input:
        krakenuniq_output = KUNIQ_OUTPUT,
        krakenuiq_report = KUNIQ_REPORT,
        taxdb = KUNIQ_TAXDB,
        cb_bam = os.path.join(
        config["output"]["host"],
        "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
    output:
        kuniq_genus_tsv = KUNIQ_GENUS_COUNT,
        kuniq_species_tsv = KUNIQ_SPECIES_COUNT
    log:
        os.path.join(config["logs"]["classifier"],
                    "rmhost_krakenuniq/{sample}_krakenuniq2sc.log")
    threads:
        40
    resources:
        mem_mb = 120000
    shell:
        '''
        python /scratch/2025-01-05/med-sucx/microcat_manuscript/benchmarks/scripts/kuniq2sc_arg.py \
        --kraken_output {input.krakenuniq_output} \
        --kraken_report {input.krakenuiq_report} \
        --taxdb {input.taxdb} \
        --cb_bam {input.cb_bam} \
        --genus_tsv {output.kuniq_genus_tsv} \
        --species_tsv {output.kuniq_species_tsv} \
        --log_file {log}
        '''

rule krakuniq_all:
    input:
        expand(KUNIQ_SPECIES_COUNT, sample=SAMPLES_ID_LIST),

