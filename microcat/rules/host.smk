import pandas as pd
import glob
import os
## beta test
sys.path.append('/data/med-sucx/host-microbiome/MicroCAT/microcat/')
import sample


def gather_fastq_files(wildcards):
    from snakemake.io import Wildcards
    fastq_files = []
    for sample in SAMPLES:
        wildcards_sample = Wildcards(fromdict={"sample": sample['sample']})
        barcode_dict = checkpoints.cellranger_unmapped_demultiplex.get(**wildcards_sample).output
        for barcode in barcode_dict.keys():
            fastq_files.extend([
                os.path.join(
                    config["output"]["host"],
                    f"cellranger_count/{sample['sample']}/unmapped_bam_CB_demultiplex/CB_{barcode}_R1.fastq"),
                os.path.join(
                    config["output"]["host"],
                    f"cellranger_count/{sample['sample']}/unmapped_bam_CB_demultiplex/CB_{barcode}_R2.fastq")
            ])
    return fastq_files






# rule starsolo_count_microbiome:
#     input:
#         # unmapped bam file from 10x
#         unmapped_bam_file = os.path.join(
#             config["output"]["host"],
#             "cellranger_run/{sample}/{sample}_unmappped2human_bam.bam")
#     output:
#         # Path to the output features.tsv file
#         features_file = os.path.join(
#             config["output"]["host"],
#             "starsolo_count_Mycobacterium_canettii/{sample}/{sample}_features.tsv"),
#         matrix_file = os.path.join(
#             config["output"]["host"],
#             "starsolo_count_Mycobacterium_canettii/{sample}/{sample}_matrix.mtx"),
#         barcodes_file = os.path.join(
#             config["output"]["host"],
#             "starsolo_count_Mycobacterium_canettii/{sample}/{sample}_barcodes.tsv"),
#     params:
#         starsolo_out = os.path.join(
#             config["output"]["host"],
#             "starsolo_count_Mycobacterium_canettii/"),
#         reference = config["params"]["host"]["starsolo"]["reference"],
#         soloCBwhitelist=config["params"]["host"]["starsolo"]["soloCBwhitelist"],
#         soloUMIlen=config["params"]["host"]["starsolo"]["soloUMIlen"],
#         soloType = config["params"]["host"]["starsolo"]["soloType"],
#         variousParams = config["params"]["host"]["starsolo"]["variousParams"],
#         threads = config["params"]["host"]["starsolo"]["threads"]
#     log:
#         os.path.join(config["logs"]["host"],
#                     "starsolo_count_Mycobacterium_canettii/{sample}_starsolo_count_Mycobacterium_canettii.log")
#     benchmark:
#         os.path.join(config["benchmarks"]["host"],
#                     "starsolo_count_Mycobacterium_canettii/{sample}_starsolo_count_Mycobacterium_canettii.benchmark")
#     shell:
#         '''
#         mkdir -p {params.starsolo_out}; 
#         cd {params.starsolo_out} ;
#         STAR \
#         --soloType {params.soloType} \
#         --soloCBwhitelist {params.soloCBwhitelist} \
#         --soloUMIstart 17 \
#         --soloUMIlen {params.soloUMIlen} \
#         --soloCBlen 16 \
#         --soloCBstart 1 \
#         --genomeDir {params.reference} \
#         --readFilesIn ../../../{input.unmapped_bam_file} \
#         --readFilesType SAM SE \
#         --soloInputSAMattrBarcodeSeq CR UR \
#         --soloInputSAMattrBarcodeQual CY UY \
#         --runThreadN {params.threads} \
#         --clipAdapterType CellRanger4 \
#         --outFilterScoreMin 30 \
#         --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
#         --soloUMIfiltering MultiGeneUMI_CR \
#         --outSAMtype BAM SortedByCoordinate\
#         --outSAMattributes CR UR CY UY CB UB \
#         --readFilesCommand samtools view -F 0x100\
#         --soloUMIdedup 1MM_CR \
#         --outFileNamePrefix ./{wildcards.sample}/\
#         {params.variousParams}  \
#         2>&1 | tee ../../../{log} ;
#         pwd ;\
#         cd ../../../;\
#         ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/features.tsv" "{output.features_file}" ;
#         ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/matrix.mtx" "{output.matrix_file}" ; 
#         ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/barcodes.tsv" "{output.barcodes_file}" ;\
#         '''


if config["params"]["host"]["starsolo"]["do"]:

    if "tenX" in config["params"]["host"]["starsolo"]["assay"]:
        # This will be executed if the string "tenX" is in the assay parameter

        if config["params"]["host"]["starsolo"]["assay"]=="tenX_v3":
            rule starsolo_10x_count:
                input:
                    # Directory containing input fastq files
                    fastqs_dir=lambda wildcards: sample.get_fastqs_dir(SAMPLES,wildcards),
                output:
                    # Path to the output features.tsv file
                    features_file = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/{sample}_features.tsv"),
                    matrix_file = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/{sample}_matrix.mtx"),
                    barcodes_file = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/{sample}_barcodes.tsv"),
                    # Path to the output unmapped bam
                    mapped_bam_file = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam")
                params:
                    barcode_reads = lambda wildcards: sample.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                    cdna_reads = lambda wildcards: sample.get_starsolo_sample_id(SAMPLES, wildcards, "fq2"),
                    starsolo_out = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/"),
                    reference = config["params"]["host"]["starsolo"]["reference"],
                    variousParams = config["params"]["host"]["starsolo"]["variousParams"],
                    threads = config["params"]["host"]["starsolo"]["threads"]
                log:
                    os.path.join(config["logs"]["host"],
                                "starsolo/{sample}_starsolo_count.log")
                benchmark:
                    os.path.join(config["benchmarks"]["host"],
                                "starsolo/{sample}_starsolo_count.benchmark")
                message: "Executing starsolo with {params.threads} threads on the following files {wildcards.sample}.Library with 10x 3' V3"
                shell:
                    '''
                    mkdir -p {params.starsolo_out}; 
                    cd {params.starsolo_out} ;
                    STAR \
                    --soloType CB_UMI_Simple \
                    --soloCBwhitelist ../data/3M-february-2018.txt \
                    --soloCBstart 1 \
                    --soloCBlen 16 \
                    --soloUMIstart 17 \
                    --soloUMIlen 12 \
                    --soloBarcodeReadLength 150 \
                    --genomeDir {params.reference} \
                    --readFilesIn {params.cdna_reads} {params.barcode_reads} \
                    --runThreadN {params.threads} \
                    --clipAdapterType CellRanger4 \
                    --outFilterScoreMin 30 \
                    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
                    --soloUMIfiltering MultiGeneUMI_CR \
                    --outSAMtype BAM Unsorted\
                    --outSAMattrRGline ID:{wildcards.sample} PL:illumina SM:{wildcards.sample} LB:tenX_v3 \
                    --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN \
                    --readFilesCommand zcat \
                    --soloUMIdedup 1MM_CR \
                    --outSAMunmapped Within \
                    --outFileNamePrefix ./{wildcards.sample}/\
                    {params.variousParams}  \
                    2>&1 | tee ../../../{log} ;
                    pwd ;\
                    cd ../../../;\
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/features.tsv" "{output.features_file}" ;
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/matrix.mtx" "{output.matrix_file}" ; 
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/barcodes.tsv" "{output.barcodes_file}" ;\
                    mv "{params.starsolo_out}/{wildcards.sample}/Aligned.out.bam" "{output.mapped_bam_file}";\
                    '''        
        if config["params"]["host"]["starsolo"]["assay"]=="tenX_v1":
            rule starsolo_10x_count:
                input:
                    # Directory containing input fastq files
                    fastqs_dir=lambda wildcards: sample.get_fastqs_dir(SAMPLES,wildcards),
                output:
                    # Path to the output features.tsv file
                    features_file = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/{sample}_features.tsv"),
                    matrix_file = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/{sample}_matrix.mtx"),
                    barcodes_file = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/{sample}_barcodes.tsv"),
                    # Path to the output unmapped bam
                    mapped_bam_file = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam")
                params:
                    barcode_reads = lambda wildcards: sample.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                    cdna_reads = lambda wildcards: sample.get_starsolo_sample_id(SAMPLES, wildcards, "fq2"),
                    starsolo_out = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/"),
                    reference = config["params"]["host"]["starsolo"]["reference"],
                    variousParams = config["params"]["host"]["starsolo"]["variousParams"],
                    threads = config["params"]["host"]["starsolo"]["threads"]
                log:
                    os.path.join(config["logs"]["host"],
                                "starsolo/{sample}_starsolo_count.log")
                benchmark:
                    os.path.join(config["benchmarks"]["host"],
                                "starsolo/{sample}_starsolo_count.benchmark")
                message: "Executing starsolo with {params.threads} threads on the following files {wildcards.sample}.Library with 10x 3' V1"
                shell:
                    '''
                    mkdir -p {params.starsolo_out}; 
                    cd {params.starsolo_out} ;
                    STAR \
                    --soloType CB_UMI_Simple \
                    --soloCBwhitelist ../data/737K-april-2014_rc.txt \
                    --soloCBstart 1 \
                    --soloCBlen 16 \
                    --soloUMIstart 17 \
                    --soloUMIlen 10 \
                    --soloBarcodeReadLength 150 \
                    --genomeDir {params.reference} \
                    --readFilesIn {input.cdna_reads} {input.barcode_reads} \
                    --runThreadN {params.threads} \
                    --clipAdapterType CellRanger4 \
                    --outFilterScoreMin 30 \
                    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
                    --soloUMIfiltering MultiGeneUMI_CR \
                    --outSAMtype BAM Unsorted\
                    --outSAMattrRGline ID:{wildcards.sample} PL:illumina SM:{wildcards.sample} LB:tenX_v3 \
                    --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN \
                    --readFilesCommand zcat \
                    --soloUMIdedup 1MM_CR \
                    --outSAMunmapped Within \
                    --outFileNamePrefix ./{wildcards.sample}/\
                    {params.variousParams}  \
                    2>&1 | tee ../../../{log} ;
                    pwd ;\
                    cd ../../../;\
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/features.tsv" "{output.features_file}" ;
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/matrix.mtx" "{output.matrix_file}" ; 
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/barcodes.tsv" "{output.barcodes_file}" ;\
                    mv "{params.starsolo_out}/{wildcards.sample}/Aligned.sortedByCoord.out.bam" "{output.mapped_bam_file}";
                    '''   

        if config["params"]["host"]["starsolo"]["assay"]=="tenX_v2":
            rule starsolo_10x_count:
                input:
                    # Directory containing input fastq files
                    fastqs_dir=lambda wildcards: sample.get_fastqs_dir(SAMPLES,wildcards),
                output:
                    # Path to the output features.tsv file
                    features_file = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/{sample}_features.tsv"),
                    matrix_file = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/{sample}_matrix.mtx"),
                    barcodes_file = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/{sample}_barcodes.tsv"),
                    # Path to the output unmapped bam
                    mapped_bam_file = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam")
                params:
                    starsolo_out = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/"),
                    reference = config["params"]["host"]["starsolo"]["reference"],
                    variousParams = config["params"]["host"]["starsolo"]["variousParams"],
                    threads = config["params"]["host"]["starsolo"]["threads"],
                    barcode_reads = lambda wildcards: sample.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                    cdna_reads = lambda wildcards: sample.get_starsolo_sample_id(SAMPLES, wildcards, "fq2")
                log:
                    os.path.join(config["logs"]["host"],
                                "starsolo/{sample}_starsolo_count.log")
                benchmark:
                    os.path.join(config["benchmarks"]["host"],
                                "starsolo/{sample}_starsolo_count.benchmark")
                resources:
                    threads=20,      # This rule needs 30 threads
                    mem_mb=100000  # This rule needs 100 GB of memory
                message: "Executing starsolo with {params.threads} threads on the following files {wildcards.sample}.Library with 10x 3' V2."
                shell:
                    '''
                    mkdir -p {params.starsolo_out}; 
                    cd {params.starsolo_out} ;
                    STAR \
                    --soloType CB_UMI_Simple \
                    --soloCBwhitelist /data/project/host-microbiome/microcat/microcat/data/737K-august-2016.txt \
                    --soloCBstart 1 \
                    --soloCBlen 16 \
                    --soloUMIstart 17 \
                    --soloUMIlen 10 \
                    --genomeDir {params.reference} \
                    --readFilesIn {params.cdna_reads}  {params.barcode_reads} \
                    --runThreadN {params.threads} \
                    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
                    --soloUMIfiltering MultiGeneUMI_CR \
                    --outSAMtype BAM SortedByCoordinate\
                    --outFilterScoreMin 30  \
                    --outSAMattrRGline ID:{wildcards.sample} PL:illumina SM:{wildcards.sample} LB:tenX_v2 \
                    --outSAMattributes CR UR CY UY CB UB\
                    --readFilesCommand zcat \
                    --soloUMIdedup 1MM_CR \
                    --outSAMunmapped Within \
                    --outFileNamePrefix ./{wildcards.sample}/\
                    {params.variousParams}  \
                    2>&1 | tee ../../../{log} ;
                    pwd ;\
                    cd ../../../;\
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/features.tsv" "{output.features_file}" ;
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/matrix.mtx" "{output.matrix_file}" ; 
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/barcodes.tsv" "{output.barcodes_file}" ;\
                    mv "{params.starsolo_out}/{wildcards.sample}/Aligned.sortedByCoord.out.bam" "{output.mapped_bam_file}";\
                    '''  
        rule starsolo_10x_unmapped_extracted_sorted:
            input:
                mapped_bam_file = os.path.join(
                    config["output"]["host"],
                    "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam")
            output:
                unmapped_bam_sorted_file = os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
            params:
                threads=16,
                unmapped_bam_unsorted_file = os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByCoord_unmapped_out.bam")
            ## because bam is sorted by Coord,it's necessary to sort it by read name
            shell:
                '''
                samtools view --threads  {params.threads}  -b -f 4   {input.mapped_bam_file}  >  {params.unmapped_bam_unsorted_file};\
                samtools sort -n  --threads  {params.threads} {params.unmapped_bam_unsorted_file} -o {output.unmapped_bam_sorted_file}
                '''
        # rule starsolo_10x_unmapped_sorted_bam:
        #     input:
        #         unmapped_bam_unsorted_file = os.path.join(
        #             config["output"]["host"],
        #             "unmapped_host/{sample}/Aligned_sortedByCoord_unmapped_out.bam")
        #     output:
        #         unmapped_sorted_bam = os.path.join(
        #             config["output"]["host"],
        #             "unmapped_host/{sample}/Aligned.out.unmapped.CBsorted.bam"),
        #     params:
        #         threads=40,
        #         tag="CB"
        #     log:
        #         os.path.join(config["logs"]["host"],
        #                     "starsolo/{sample}/unmapped_sorted_bam.log")
        #     benchmark:
        #         os.path.join(config["benchmarks"]["host"],
        #                     "starsolo/{sample}/unmapped_sorted_bam.benchmark")
        #     shell:
        #         '''
        #         samtools sort -@ {params.threads} -t {params.tag} -o {output.unmapped_sorted_bam}  {input.unmapped_bam_unsorted_file};
        #         '''
        # rule starsolo_10X_demultiplex_bam_by_cell_barcode:
        #     input:
        #         unmapped_sorted_bam = os.path.join(
        #             config["output"]["host"],
        #             "unmapped_host/{sample}/Aligned.out.unmapped.RGsorted.bam")
        #     output:
        #         unmapped_bam_demultiplex_dir = directory(os.path.join(
        #             config["output"]["host"],
        #             "unmapped_host/{sample}/unmapped_bam_CB_demultiplex/"))
        #     params:
        #         threads = 40, # Number of threads
        #         tag="CB"
        #     log:
        #         os.path.join(
        #             config["logs"]["host"],
        #             "starsolo_count/{sample}/demultiplex_bam_by_read_group.log")
        #     benchmark:
        #         os.path.join(
        #             config["benchmarks"]["host"], 
        #             "starsolo_count/{sample}/demultiplex_bam_by_read_group.benchmark")
        #     shell:
        #         """
        #         python /data/project/host-microbiome/microcat/microcat/scripts/spilt_bam_by_tag.py --tag {params.tag} --bam_path {input.unmapped_sorted_bam} --output_dir {output.unmapped_bam_demultiplex_dir}  --log_file {log}
        #         """
        rule starsolo_10X_all:
            input:
                # expand(os.path.join(
                #     config["output"]["classifier"],
                #     "rmhost_classified_qc/kmer_UMI/{sample}/{sample}_kraken2_sckmer.txt"),sample=SAMPLES_ID_LIST),
                # expand(os.path.join(config["output"]["classifier"],
                #     "microbiome_matrix_build/{sample}/data.txt"),sample=SAMPLES_ID_LIST)
                expand(os.path.join(config["output"]["host"],"unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),sample=SAMPLES_ID_LIST)
                # expand(os.path.join(
                #     config["output"]["classifier"],
                #     "rmhost_classified_qc/kmer_UMI/{sample}/{sample}_kraken2_sckmer_correlation_test.txt"),sample=SAMPLES_ID_LIST) 
        
    else:
        rule starsolo_10X_all:
            input: 


    if config["params"]["host"]["starsolo"]["assay"]=="SmartSeq2" or config["params"]["host"]["starsolo"]["assay"]=="SmartSeq":

        # rule run_fastp:
        #     group:
        #         "run_fastp"
        #     conda:
        #         join(ENV_DIR, "fastp.yml")
        #     input:
        #         FASTQ1_FILE,
        #         FASTQ2_FILE
        #     output:
        #         TRIMMED_FASTQ1_FILE,
        #         TRIMMED_FASTQ2_FILE,
        #         TRIMMED_UNPAIRED_FILE,
        #         FAILED_READS_FILE,
        #         FASTP_JSON_REPORT,
        #         FASTP_HTML_REPORT
        #     threads:
        #         6
        #     benchmark:
        #         "benchmarks/{patient}-{sample}-{cell}.run_fastp.benchmark.txt"
        #     shell:
        #         "fastp -w {threads} "
        #         "--unqualified_percent_limit 40 " # filter reads where 40% of bases have phred quality < 15
        #         "--cut_tail " # use defaults --cut_window_size 4 --cut_mean_quality 20
        #         "--low_complexity_filter " # filter reads with less than 30% complexity (30% of the bases are different from the preceeding base)
        #         "--trim_poly_x " # trim poly X's - useful for trimming polyA tails
        #         "-i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} "
        #         "--unpaired1 {output[2]} --unpaired2 {output[2]} --failed_out {output[3]} "
        #         "-j {output[4]} -h {output[5]} "
        #         + config["params"]["fastp"]

        # rule generate_pe_manifest_file:

        # input:
        #     config["units"]
        # output:
        #     PE_MANIFEST_FILE
        
        # script:
        #     "../src/generate_pe_manifest_file.py"
        
        rule starsolo_smartseq_count:
            # Input files
            input:
                # Path to the input manifest file
                manifest = config["params"]["host"]["starsolo"]["manifest"],
                fastqs_dir=lambda wildcards: sample.get_fastqs_dir(SAMPLES,wildcards),
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
                # Path to the output unmapped fastq file for read1
                # ummapped_fastq_1 = os.path.join(
                #     config["output"]["host"],
                #     "starsolo_count/Unmapped.out.mate1"),
                # # Path to the output unmapped fastq file for read2
                # ummapped_fastq_2 = os.path.join(
                #     config["output"]["host"],
                #     "starsolo_count/Unmapped.out.mate2"),
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
                # Type of sequencing library
                soloType = config["params"]["host"]["starsolo"]["soloType"],
                SAMattrRGline = sample.get_SAMattrRGline_from_manifest(config["params"]["host"]["starsolo"]["manifest"]),
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
        rule starsolo_smartseq_unmapped_extracted_sorted:
            input:
                mapped_bam_file = os.path.join(
                    config["output"]["host"],
                    "starsolo_count/Aligned_out.bam")
            output:
                unmapped_bam_unsorted_file = os.path.join(
                    config["output"]["host"],
                    "starsolo_count/Aligned_out_unmapped.bam")
            params:
                threads=16
            shell:
                '''
                samtools view --threads  {params.threads}  -b -f 4   {input.mapped_bam_file}  >  {output.unmapped_bam_unsorted_file};\
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
            shell:
                '''
                samtools sort -@ {params.threads} -t {params.tag} -o {output.unmapped_sorted_bam}  {input.unmapped_bam_unsorted_file};
                '''

        rule starsolo_smartseq_demultiplex_bam_by_read_group:
            input:
                unmapped_sorted_bam = os.path.join(
                    config["output"]["host"],
                    "starsolo_count/Aligned_out_unmapped_RGsorted.bam")
            output:
                unmapped_bam_demultiplex_dir = directory(os.path.join(
                    config["output"]["host"],
                    "unmapped_host/"))
            params:
                threads = 40, # Number of threads
                tag="RG"
            log:
                os.path.join(
                    config["logs"]["host"],
                    "starsolo_count/demultiplex_bam_by_read_group.log")
            benchmark:
                os.path.join(
                    config["benchmarks"]["host"], 
                    "starsolo_count/demultiplex_bam_by_read_group.benchmark")
            shell:
                """
                python /data/project/host-microbiome/microcat/microcat/scripts/spilt_bam_by_tag.py --tag {params.tag} --bam_path {input.unmapped_sorted_bam} --output_dir {output.unmapped_bam_demultiplex_dir}  --log_file {log}
                """
        rule starsolo_smartseq_all:
            input:
                directory(os.path.join(
                    config["output"]["host"],
                    "unmapped_host/"))
                # os.path.join(
                #     config["output"]["host"],
                #     "starsolo_count/Aligned_out_unmapped_RGsorted.bam")
        
    else:
        rule starsolo_smartseq_all:
            input: 
    #ALL Input
    rule starsolo_all:
        input: 
            rules.starsolo_smartseq_all.input,
            rules.starsolo_10X_all.input

else:
    rule starsolo_all:
        input:

if config["params"]["host"]["cellranger"]["do"]:
# expected input format for FASTQ file
# cellranger call to process the raw samples
    rule cellranger_count:
        input:
            # fastqs_dir = config["params"]["data_dir"],
            # r1 = lambda wildcards: get_sample_id(SAMPLES, wildcards, "fq1"),
            # r2 = lambda wildcards: get_sample_id(SAMPLES, wildcards, "fq2")
            fastqs_dir=lambda wildcards: sample.get_fastqs_dir(SAMPLES,wildcards),
        output:
            features_file = os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/{sample}_features.tsv"),
            matrix_file = os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/{sample}_matrix.mtx"),
            barcodes_file = os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/{sample}_barcodes.tsv"),
            mapped_bam_file = os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/{sample}_mappped2human_bam.bam"),
            mapped_bam_index_file = os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/{sample}_mappped2human_bam.bam.bai")
        priority: 10
        params:
            cr_out = os.path.join(
                config["output"]["host"],
                "cellranger_count/"),
            reference = config["params"]["host"]["cellranger"]["reference"],
            # local_cores = config["params"]["host"]["cellranger"]["local_cores"],
            metrics_summary = os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/{sample}.metrics_summary.csv"),
            web_summary = os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/{sample}.web_summary.html"),
            SampleID="{sample}",
            variousParams = config["params"]["host"]["cellranger"]["variousParams"],
        # resources:
        #     mem_mb=config["tools"]["cellranger_count"]["mem_mb"],
        #     runtime=config["tools"]["cellranger_count"]["runtime"],
        threads: 
            config["params"]["host"]["cellranger"]["threads"]
        resources:
            mem_mb=102400,
            disk_mb=10000
        log:
            os.path.join(config["logs"]["host"],
                        "cellranger/{sample}_cellranger_count.log")
        benchmark:
            os.path.join(config["benchmarks"]["host"],
                        "cellranger/{sample}_cellranger_count.benchmark")
        # NOTE: cellranger count function cannot specify the output directory, the output is the path you call it from.
        # Therefore, a subshell is used here.
        shell:
            '''
            cd {params.cr_out}  
            cellranger count \
            --id={params.SampleID} \
            --sample={params.SampleID}  \
            --transcriptome={params.reference} \
            --localcores={threads} \
            --fastqs={input.fastqs_dir} \
            --nosecondary \
            {params.variousParams} \
            2>&1 | tee ../../../{log} ;  
            cd ../../../;
            gunzip {params.cr_out}{params.SampleID}/outs/filtered_feature_bc_matrix/features.tsv.gz ; 
            gunzip {params.cr_out}{params.SampleID}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ; 
            gunzip {params.cr_out}{params.SampleID}/outs/filtered_feature_bc_matrix/matrix.mtx.gz ; 
            ln -sr "{params.cr_out}{params.SampleID}/outs/filtered_feature_bc_matrix/features.tsv" "{output.features_file}"; 
            ln -sr "{params.cr_out}{params.SampleID}/outs/filtered_feature_bc_matrix/matrix.mtx" "{output.matrix_file}"; 
            ln -sr "{params.cr_out}{params.SampleID}/outs/filtered_feature_bc_matrix/barcodes.tsv" "{output.barcodes_file}" ; 
            ln -sr "{params.cr_out}{params.SampleID}/outs/web_summary.html" "{params.web_summary}" ; 
            ln -sr "{params.cr_out}{params.SampleID}/outs/metrics_summary.csv" "{params.metrics_summary}";
            ln -sr "{params.cr_out}{params.SampleID}/outs/possorted_genome_bam.bam" "{output.mapped_bam_file}";
            ln -sr "{params.cr_out}{params.SampleID}/outs/possorted_genome_bam.bam.bai" "{output.mapped_bam_index_file}";
            '''
    # rule cellranger_unmapped_extracted:
    #     input:
    #         mapped_bam_file = os.path.join(
    #         config["output"]["host"],
    #         "cellranger_count/{sample}/{sample}_mappped2human_bam.bam")
    #     output:
    #         unmapped_bam_unsorted_file = os.path.join(
    #         config["output"]["host"],
    #         "unmapped_host/{sample}/{sample}_unmappped2human_unsorted_bam.bam")
    #     shell:
    #         '''
    #         samtools view --threads  {threads}  -b -f 4   {input.mapped_bam_file}  >  {output.unmapped_bam_unsorted_file};
    #         '''
    rule cellranger_unmapped_extracted:
        input:
            mapped_bam_file = os.path.join(
            config["output"]["host"],
            "cellranger_count/{sample}/{sample}_mappped2human_bam.bam")
        output:
            unmapped_bam_unsorted_file = os.path.join(
            config["output"]["host"],
            "cellranger_count/{sample}/{sample}_unmappped2human_sorted_bam.bam")
        threads: 
            16
        resources:
            mem_mb=16384
        shell:
            '''
            samtools view --threads  {threads}  -b -f 4   {input.mapped_bam_file}  >  {output.unmapped_bam_unsorted_file};
            '''
    rule cellranger_unmapped_sorted:
        input:
            unmapped_bam_unsorted_file = os.path.join(
            config["output"]["host"],
            "cellranger_count/{sample}/{sample}_unmappped2human_sorted_bam.bam")
        output:
            unmapped_bam_CBsorted_file = os.path.join(
            config["output"]["host"],
            "cellranger_count/{sample}/{sample}_unmappped2human_CB_sorted_bam.bam")
        params:
            threads=16,
            unmapped_bam_unsorted_file = os.path.join(
                config["output"]["host"],
                "unmapped_host/{sample}/Aligned_sortedByCoord_unmapped_out.bam")
        ## because bam is sorted by Coord,it's necessary to sort it by read name
        shell:
            '''
            samtools view --threads  {params.threads}  -b -f 4   {input.mapped_bam_file}  >  {params.unmapped_bam_unsorted_file};\
            samtools sort -n  --threads  {params.threads} {params.unmapped_bam_unsorted_file} -o {output.unmapped_bam_sorted_file}
            '''
    # rule cellranger_unmapped_sorted:
    #     input:
    #         unmapped_bam_unsorted_file = os.path.join(
    #         config["output"]["host"],
    #         "unmapped_host/{sample}/{sample}_unmappped2human_unsorted_bam.bam")
    #     output:
    #         unmapped_bam_CBsorted_file = os.path.join(
    #         config["output"]["host"],
    #         "unmapped_host/{sample}/{sample}_unmappped2human_CB_sorted_bam.bam")
    #     params:
    #         threads=40,
    #         tag="CB"
    #     shell:
    #         '''
    #         samtools sort -@ {params.threads} -t {params.tag} -o {output.unmapped_bam_CBsorted_file}  {input.unmapped_bam_unsorted_file};
    #         '''
    # rule cellranger_unmapped_sorted:
    #     input:
    #         unmapped_bam_unsorted_file = os.path.join(
    #         config["output"]["host"],
    #         "unmapped_host/{sample}/{sample}_unmappped2human_unsorted_bam.bam")
    #     output:
    #         unmapped_bam_CBsorted_file = os.path.join(
    #         config["output"]["host"],
    #         "unmapped_host/{sample}/{sample}_unmappped2human_CB_sorted_bam.bam")
    #     params:
    #         threads=40,
    #         tag="CB"
    #     shell:
    #         '''
    #         samtools sort -@ {params.threads} -t {params.tag} -o {output.unmapped_bam_CBsorted_file}  {input.unmapped_bam_unsorted_file};
    #         '''
    #since output barcode.bam is unknown, here we use checkpoint
    # checkpoint cellranger_unmapped_demultiplex:
    #     input:
    #         unmapped_bam_CBsorted_file = os.path.join(
    #         config["output"]["host"],
    #         "cellranger_count/{sample}/{sample}_unmappped2human_CB_sorted_bam.bam"),
    #     output:
    #         unmapped_bam_CB_demultiplex_dir = directory(os.path.join(
    #             config["output"]["host"],
    #             "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/"))
    #     params:
    #         threads = 40, # Number of threads
    #         tag="CB"
    #     log:
    #         os.path.join(
    #             config["logs"]["host"],
    #             "cellranger_count/{sample}/cellranger_unmapped_demultiplex_by_CB.log")
    #     benchmark:
    #         os.path.join(
    #             config["benchmarks"]["host"], 
    #             "cellranger_count/{sample}/cellranger_unmapped_demultiplex_by_CB.benchmark")
    #     shell:
    #         """
    #         python /data/project/host-microbiome/microcat/microcat/scripts/spilt_bam_by_tag.py --tag {params.tag} \
    #         --bam_path {input.unmapped_bam_CBsorted_file} \
    #         --output_dir {output.unmapped_bam_CB_demultiplex_dir}  \
    #         --log_file {log};
    #         """
    # rule paired_bam_to_fastq:
    #     input:
    #         bam_file=aggregate_CB_bam_output("{sample}"),
    #     output:
    #         unmapped_fastq_1 = os.path.join(
    #         config["output"]["host"],
    #         "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R1.fastq"),
    #         unmapped_fastq_2 = os.path.join(
    #         config["output"]["host"],
    #         "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R2.fastq")
    #     log:
    #         os.path.join(
    #             config["logs"]["host"],
    #             "cellranger_count/{sample}/cell_level/CB_{barcode}_paired_bam_to_fastq.log")
    #     benchmark:
    #         os.path.join(
    #             config["benchmarks"]["host"], 
    #             "cellranger_count/{sample}/cell_level/CB_{barcode}_paired_bam_to_fastq.benchmark")
    #     threads:
    #         16
    #     priority: 11
    #     shell:
    #         '''
    #         samtools fastq --threads {threads} {input.bam_file} -1 {output.unmapped_fastq_1} -2 {output.unmapped_fastq_2}
    #         '''

    rule cellranger_all:
        input:
            expand(os.path.join(
            config["output"]["host"],
            "cellranger_count/{sample}/{sample}_unmappped2human_sorted_bam.bam"),sample=SAMPLES_ID_LIST)

else:
    rule cellranger_all:
        input:

rule host_all:
    input:
        rules.starsolo_all.input,
        rules.cellranger_all.input,

    



