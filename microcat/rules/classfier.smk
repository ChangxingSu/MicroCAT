# import glob
# import os
# ## beta test
# sys.path.append('/data/project/host-microbiome/microcat/microcat/')
# import sample
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

# def get_CB_bam_files(wildcards):
#     bam_dir = os.path.join(
#         config["output"]["host"],
#         "cellranger_count",
#         wildcards.sample,
#         "unmapped_bam_CB_demultiplex"
#     )
#     return glob.glob(os.path.join(bam_dir, "CB_*.bam"))


# def aggregate_CB_bam_output(wildcards):
#     demultiplex_output = checkpoints.cellranger_unmapped_demultiplex.get(**wildcards).output.unmapped_bam_CB_demultiplex_dir
#     Barcode_list, = glob_wildcards(os.path.join(demultiplex_output, "CB_{barcode}.bam"))
#     return expand(os.path.join(
#         config["output"]["host"],
#         "cellranger_count/{sample}//unmapped_bam_CB_demultiplex/CB_{barcode}.bam"),
#         sample=wildcards.sample,
#         barcode=Barcode_list)



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

rule paired_bam_to_fastq:
    input:
        unmapped_bam_sorted_file =os.path.join(
        config["output"]["host"],
        "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
    output:
        unmapped_fastq = os.path.join(
            config["output"]["host"],
            "unmapped_host/{sample}/{sample}_unmappped2human_bam.fastq")
    # threads:
    #     16
    resources:
        mem_mb=100000
    priority: 11
    conda:
        config["envs"]["star"]
    shell:
        '''
        samtools fastq --threads {threads}  -n {input.unmapped_bam_sorted_file}  > {output.unmapped_fastq}
        '''


if config["params"]["classifier"]["kraken2uniq"]["do"]:
    rule kraken2uniq_classified:
        input:
            unmapped_fastq = os.path.join(
                config["output"]["host"],
                "unmapped_host/{sample}/{sample}_unmappped2human_bam.fastq")
        output:
            krak2_classified_output_fq = os.path.join(
                config["output"]["classifier"],
                "rmhost_classified_output/{sample}/{sample}_kraken2_classified.fq"),
            krak2_unclassified_output_fq = os.path.join(
                config["output"]["classifier"],
                "rmhost_unclassified_output/{sample}/{sample}_kraken2_unclassified.fq"),
            krak2_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_output/{sample}/{sample}_kraken2_output.txt"),
            krak2_report = os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_report/custom/{sample}/{sample}_kraken2_report.txt"),
            krak2_std_report=os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_report/standard/{sample}/{sample}_kraken2_std_report.txt"),
            krak2_mpa_report=os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_report/mpa/{sample}/{sample}_kraken2_mpa_report.txt")
        params:
            database = config["params"]["classifier"]["kraken2uniq"]["kraken2_database"],
            kraken2mpa_script = config["scripts"]["kraken2mpa"]
            #since kraken2 acquire specific input fomrat "#fq",so we put it it params
            # krak2_classified_output=os.path.join(
            #     config["output"]["classifier"],
            #     "classified_output/{sample}/{sample}_kraken2_classified#.fq"),
            # krak2_unclassified_output=os.path.join(
            #     config["output"]["classifier"],
            #     "unclassified_output/{sample}/{sample}_kraken2_unclassified#.fq")

        resources:
            mem_mb=config["params"]["classifier"]["kraken2uniq"]["mem_mb"]
        priority: 12
        threads: 
            config["params"]["classifier"]["kraken2uniq"]["threads"]
        benchmark:
            os.path.join(config["benchmarks"]["classifier"],
                        "rmhost_kraken2uniq/{sample}_kraken2uniq_classifier_benchmark.log")
        log:
            os.path.join(config["logs"]["classifier"],
                        "rmhost_kraken2uniq/{sample}_kraken2uniq_classifier.log")
        conda:
            config["envs"]["kraken2"]
        shell:
            '''
            kraken2 --db {params.database} \
            --threads 20 \
            --classified-out {output.krak2_classified_output_fq}\
            --unclassified-out {output.krak2_unclassified_output_fq}\
            --output {output.krak2_output} \
            --report {output.krak2_report} \
            --report-minimizer-data \
            {input.unmapped_fastq} \
            --use-names \
            --memory-mapping \
            2>&1 | tee {log};\
            cut -f 1-3,6-8 {output.krak2_report} > {output.krak2_std_report};\
            python {params.kraken2mpa_script} -r {output.krak2_std_report} -o {output.krak2_mpa_report};
            '''
    rule extract_kraken2_reads:
        input:
            krak2_classified_output_fq = os.path.join(
                config["output"]["classifier"],
                "rmhost_classified_output/{sample}/{sample}_kraken2_classified.fq"),
            krak2_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_output/{sample}/{sample}_kraken2_output.txt"),
            krak2_report = os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_report/custom/{sample}/{sample}_kraken2_report.txt"),
        output:
            krak2_extracted_output_fq = os.path.join(
                config["output"]["classifier"],
                "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_extracted_classified.fq"),
        log:
            os.path.join(config["logs"]["classifier"],
                        "rmhost_kraken2uniq_extracted/{sample}_kraken2uniq_classifier_extracted.log")
        benchmark:
            os.path.join(config["benchmarks"]["classifier"],
                            "rmhost_kraken2uniq/{sample}_kraken2uniq_classifier_extracted_benchmark.log")
        params:
            extract_kraken_reads_script = config["scripts"]["extract_kraken_reads"]
        conda:
            config["envs"]["kmer_qc"]
        priority: 
            13
        shell:
            '''
            python {params.extract_kraken_reads_script} \
            -k {input.krak2_output} \
            -s1 {input.krak2_classified_output_fq} \
            --report {input.krak2_report}\
            -o {output.krak2_extracted_output_fq} \
            --taxid 9606 \
            --exclude \
            --include-parents \
            2>&1 | tee {log};
            '''

    rule extract_kraken2_report:
        input:
            krak2_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_output/{sample}/{sample}_kraken2_output.txt"),
            krak2_report = os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_report/custom/{sample}/{sample}_kraken2_report.txt"),
            krak2_mpa_report=os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_report/mpa/{sample}/{sample}_kraken2_mpa_report.txt")
        output:
            krak2_extracted_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_extracted_classified_output.txt"),
        log:
            os.path.join(config["logs"]["classifier"],
                        "rmhost_kraken2uniq_extracted/{sample}_kraken2uniq_classifier_report_extracted.log")
        params:
            extract_microbiome_output_script = config["scripts"]["extract_microbiome_output"]
        threads:
            8
        priority: 
            14
        conda:
            config["envs"]["kmer_qc"]
        shell:
            '''
            pwd
            Rscript {params.extract_microbiome_output_script} \
            --output_file {input.krak2_output} \
            --kraken_report {input.krak2_report} \
            --mpa_report {input.krak2_mpa_report} \
            --extract_file {output.krak2_extracted_output}\
            --cores {threads} \
            --ntaxid 6000 \
            2>&1 | tee {log};
            '''

    rule k_mer_test:
        input:
            krak2_extracted_output_fq = os.path.join(
                config["output"]["classifier"],
                "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_extracted_classified.fq"),
            krak2_extracted_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_extracted_classified_output.txt"),
            krak2_report = os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_report/custom/{sample}/{sample}_kraken2_report.txt"),
            krak2_mpa_report=os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_report/mpa/{sample}/{sample}_kraken2_mpa_report.txt")
        output:
            krak2_sckmer_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_classified_qc/kmer_UMI/{sample}/{sample}_kraken2_sckmer.txt"),
            krak2_sckmer_test_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_classified_qc/kmer_UMI/{sample}/{sample}_kraken2_sckmer_correlation_test.txt"),
        log:
            os.path.join(config["logs"]["classifier"],
                        "classified_qc/kmer_UMI/{sample}_kraken2uniq_sckmer.log")
        priority: 
            15
        params:
            cb_len = config["params"]["classifier"]["sckmer"]["cb_len"],
            umi_len = config["params"]["classifier"]["sckmer"]["umi_len"],
            min_frac = config["params"]["classifier"]["sckmer"]["min_frac"],
            kmer_len = config["params"]["classifier"]["sckmer"]["kmer_len"],
            nsample = config["params"]["classifier"]["sckmer"]["nsample"],
            sckmer_unpaired_script= config["scripts"]["sckmer_unpaired"]
        conda:
            config["envs"]["kmer_qc"]
        shell:
            '''
            Rscript {params.sckmer_unpaired_script} \
            --fa1 {input.krak2_extracted_output_fq} \
            --microbiome_output_file {input.krak2_extracted_output} \
            --cb_len {params.cb_len} \
            --umi_len {params.umi_len} \
            --kraken_report {input.krak2_report} \
            --mpa_report {input.krak2_mpa_report} \
            --min_frac {params.min_frac} \
            --kmer_len {params.kmer_len} \
            --nsample {params.nsample} \
            --kmer_file {output.krak2_sckmer_output} \
            --kmer_test_file {output.krak2_sckmer_test_output} \
            2>&1 | tee {log};
            '''
    rule sample_denosing:
        input:
            krak2_sckmer_test_output = os.path.join(
                                config["output"]["classifier"],
                                "rmhost_classified_qc/kmer_UMI/{sample}/{sample}_kraken2_sckmer_correlation_test.txt"),
        output:
            krak2_sample_denosing = os.path.join(
                                config["output"]["classifier"],
                                "rmhost_classified_qc/kmer_UMI/{sample}/{sample}_sample_denosing.txt"),
        priority: 
            16
        params:
            krak2_report_dir = os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_report/custom/"),
            SampleID="{sample}",
            sample_denosing_script= config["scripts"]["sample_denosing"]
        conda:
            config["envs"]["kmer_qc"]
        shell:
            '''
            Rscript  {params.sample_denosing_script}\
            --path {params.krak2_report_dir} \
            --kmer_data {input.krak2_sckmer_test_output} \
            --out_path {output.krak2_sample_denosing} \
            --sample_names {params.SampleID} \
            --min_reads 2 
            '''

    rule krak2_output_denosing:
        input:
            krak2_extracted_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_extracted_classified_output.txt"),
            krak2_sample_denosing = os.path.join(
                config["output"]["classifier"],
                "rmhost_classified_qc/kmer_UMI/{sample}/{sample}_sample_denosing.txt"),
        output:
            krak2_output_denosing = os.path.join(
                config["output"]["classifier"],
                "rmhost_classified_qc/kmer_UMI/{sample}/{sample}_kraken2_output_denosing.txt"),
        conda:
            config["envs"]["kmer_qc"]
        priority: 
            17
        params:
            krak2_output_denosing_script = config["scripts"]["krak2_output_denosing"]
        shell:
            '''
            Rscript {params.krak2_output_denosing_script} \
            --output_file {input.krak2_extracted_output} \
            --taxa {input.krak2_sample_denosing} \
            --out_krak2_denosing {output.krak2_output_denosing}
            '''
    rule krak2_matrix_build:
        input:
            krak2_output_denosing = os.path.join(
                config["output"]["classifier"],
                "rmhost_classified_qc/kmer_UMI/{sample}/{sample}_kraken2_output_denosing.txt"),
            unmapped_bam_sorted_file =os.path.join(
                config["output"]["host"],
                "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
        output:
            krak2_output_denosing_label = os.path.join(config["output"]["classifier"],"microbiome_matrix_build/{sample}/data.txt")
        params:
            database = config["params"]["classifier"]["kraken2uniq"]["kraken2_database"],
            threads = 20,
            matrix_outdir = os.path.join(config["output"]["classifier"],
                "microbiome_matrix_build/{sample}/"),
            kraken2sc_script = config["scripts"]["kraken2sc"]
        priority: 
            18
        log:
            os.path.join(config["logs"]["classifier"],
                        "microbiome_matrix_build/{sample}_matrix.log")
        conda:
            config["envs"]["kmer_python"]
        shell:
            '''
            python {params.kraken2sc_script} \
            --bam {input.unmapped_bam_sorted_file} \
            --kraken_output {input.krak2_output_denosing}  \
            --dbfile {params.database} \
            --log_file {log} \
            --processors {params.threads} \
            --outdir {params.matrix_outdir}
            '''
            
    rule kraken2uniq_classified_all:
        input:
            expand(os.path.join(config["output"]["classifier"],"microbiome_matrix_build/{sample}/data.txt"),sample=SAMPLES_ID_LIST)
else:
    rule kraken2uniq_classified_all:
        input:    

if config["params"]["classifier"]["krakenuniq"]["do"]:
    rule krakenuniq_classifier:
        input:
            unmapped_fastq = os.path.join(
                config["output"]["host"],
                "unmapped_host/{sample}/{sample}_unmappped2human_bam.fastq")
        output:
            krakenuniq_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_krakenuniq_output/{sample}/{sample}_krakenuniq_output.txt"),
            krakenuniq_report = os.path.join(
                config["output"]["classifier"],
                "rmhost_krakenuniq_report/custom/{sample}/{sample}_krakenuniq_report.txt"),
            krakenuniq_classified_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_krakenuniq_classified_output/{sample}/{sample}_krakenuniq_classified.fq"),
            krakenuniq_unclassified_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_krakenuniq_unclassified_output/{sample}/{sample}_krakenuniq_unclassified.fq")
        params:
            database = config["params"]["classifier"]["krakenuniq"]["krakenuniq_database"],
            threads=config["params"]["classifier"]["krakenuniq"]["threads"],
            estimate_precision=config["params"]["classifier"]["krakenuniq"]["estimate_precision"]
        benchmark:
            os.path.join(config["benchmarks"]["classifier"],
                        "rmhost_krakenuniq/{sample}_kraken2uniq_classifier_benchmark.tsv")
        log:
            os.path.join(config["logs"]["classifier"],
                        "rmhost_krakenuniq/{sample}_krakenuniq_classifier.log")
        resources:
            mem_mb=config["params"]["classifier"]["krakenuniq"]["mem_mb"]
        conda:
            config["envs"]["krakenuniq"]
        message:
            "Classifier: Performing Taxonomic Classifcation of Sample {sample} with krakenuniq."
        shell:
            '''
            krakenuniq --db {params.database} \
            --threads {params.threads} \
            --hll-precision {params.estimate_precision} \
            --classified-out {output.krakenuniq_classified_output}\
            --unclassified-out {output.krakenuniq_unclassified_output}\
            --output {output.krakenuniq_output} \
            --report-file {output.krakenuniq_report} \
            {input.unmapped_fastq}  \
            --preload \
            2>&1 | tee {log}
            '''
    # rule krakenuniq_cell_level_classifier:
    #     input:
    #         r1 = expand(os.path.join(
    #             config["output"]["host"],
    #             "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R1.fastq"), barcode=get_barcodes(wildcards.sample)),
    #         r2 = expand(os.path.join(
    #             config["output"]["host"],
    #             "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R2.fastq"), barcode=get_barcodes(wildcards.sample))
    #     output:
    #         krakenuniq_output = expand(os.path.join(
    #             config["output"]["classifier"],
    #             "krakenuniq_output/{sample}/cell_level/{sample}_{barcode}_krakenuniq_output.txt"), barcode=get_barcodes(wildcards.sample)),
    #         krakenuniq_report = expand(os.path.join(
    #             config["output"]["classifier"],
    #             "krakenuniq_report/custom/{sample}/cell_level/{sample}_{barcode}_krakenuniq_report.txt"), barcode=get_barcodes(wildcards.sample)),
    #         krakenuniq_classified_output = expand(os.path.join(
    #             config["output"]["classifier"],
    #             "krakenuniq_classified_output/{sample}/cell_level/{sample}_{barcode}_krakenuniq_classified.fq"), barcode=get_barcodes(wildcards.sample)),
    #         krakenuniq_unclassified_output = expand(os.path.join(
    #             config["output"]["classifier"],
    #             "krakenuniq_classified_output/{sample}/cell_level/{sample}_{barcode}_krakenuniq_unclassified.fq"), barcode=get_barcodes(wildcards.sample))
    #     params:
    #         database = config["params"]["classifier"]["krakenuniq"]["krakenuniq_database"],
    #         threads=config["params"]["classifier"]["krakenuniq"]["threads"],
    #         estimate_precision=config["params"]["classifier"]["krakenuniq"]["estimate_precision"]
    #     benchmark:
    #         expand(os.path.join(config["benchmarks"]["classifier"],
    #                     "krakenuniq/{sample}/cell_level/{sample}_{barcode}_krakenuniq_classifier_benchmark.log"), barcode=get_barcodes(wildcards.sample))
    #     log:
    #         expand(os.path.join(config["logs"]["classifier"],
    #                     "krakenuniq/{sample}/cell_level/{sample}_{barcode}_krakenuniq_classifier.log"), barcode=get_barcodes(wildcards.sample))
    #     conda:
    #         config["envs"]["krakenuniq"]
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
    rule krakenuniq_classified_all:
        input:   
            expand(os.path.join(
                config["output"]["classifier"],
                "rmhost_krakenuniq_output/{sample}/{sample}_krakenuniq_output.txt"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                config["output"]["classifier"],
                "rmhost_krakenuniq_report/custom/{sample}/{sample}_krakenuniq_report.txt"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                config["output"]["classifier"],
                "rmhost_krakenuniq_classified_output/{sample}/{sample}_krakenuniq_classified.fq"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                config["output"]["classifier"],
                "rmhost_krakenuniq_unclassified_output/{sample}/{sample}_krakenuniq_unclassified.fq"),sample=SAMPLES_ID_LIST)

else:
    rule krakenuniq_classified_all:
        input:    

if config["params"]["classifier"]["pathseq"]["do"]:
    rule pathseq_classified:
        input:
            unmapped_bam_sorted_file =os.path.join(
                config["output"]["host"],
                "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
        output:
            pathseq_classified_bam_file = os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_classified.bam"),
            pathseq_output = os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_classified.txt"),
            filter_metrics = os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_filter_metrics.txt"),
            score_metrics = os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_score_metrics.txt"),
        params:
            host_bwa_image = config["params"]["classifier"]["pathseq"]["host_bwa_image"],
            microbe_bwa_image = config["params"]["classifier"]["pathseq"]["microbe_bwa_image"],
            microbe_dict_file = config["params"]["classifier"]["pathseq"]["microbe_dict"],
            host_hss_file = config["params"]["classifier"]["pathseq"]["host_bfi"],
            taxonomy_db = config["params"]["classifier"]["pathseq"]["taxonomy_db"],
            pathseq_output_dir = os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/")
        resources:
            mem_mb=config["params"]["classifier"]["pathseq"]["mem_mb"]
        priority: 12
        threads: 
            config["params"]["classifier"]["pathseq"]["threads"]
        conda:
            config["envs"]["pathseq"]
        benchmark:
            os.path.join(config["benchmarks"]["classifier"],
                        "rmhost_pathseq/{sample}_pathseq_classifier_benchmark.log")
        log:
            os.path.join(config["logs"]["classifier"],
                        "rmhost_pathseq/{sample}_pathseq_classifier.log")
        shell:
            '''
            mkdir -p {params.pathseq_output_dir};\
            gatk PathSeqPipelineSpark \
            --filter-duplicates false \
            --min-score-identity .7 \
            --input {input.unmapped_bam_sorted_file} \
            --filter-bwa-image {params.host_bwa_image} \
            --kmer-file {params.host_hss_file} \
            --microbe-bwa-image {params.microbe_bwa_image} \
            --microbe-dict {params.microbe_dict_file} \
            --taxonomy-file {params.taxonomy_db} \
            --output {output.pathseq_classified_bam_file}\
            --scores-output {output.pathseq_output}\
            --filter-metrics {output.filter_metrics}\
            --score-metrics {output.score_metrics}\
            --java-options "-Xmx200g" \
            2>&1 | tee {log}\
            '''
    rule pathseq_extract_paired_bam:
        input:
            pathseq_classified_bam_file = os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_classified.bam"),
        output:
            pathseq_classified_paired_bam_file = temp(os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_paired_classified.bam"))
        params:
            threads=8
        conda:
            config["envs"]["star"]
        shell:
            '''
            samtools view --threads {params.threads} -h -b -f 1 -o {output.pathseq_classified_paired_bam_file} {input.pathseq_classified_bam_file}
            '''
    rule pathseq_sort_extract_paired_bam:
        input:
            pathseq_classified_paired_bam_file = temp(os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_paired_classified.bam")),
        output:
            pathseq_classified_paired_sorted_bam_file = temp(os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_sorted_paired_classified.bam"))
        params:
            threads = 8
        conda:
            config["envs"]["star"]
        shell:
            '''
            samtools sort --threads {params.threads} -n -o {output.pathseq_classified_paired_sorted_bam_file} {input.pathseq_classified_paired_bam_file} 
            '''
    rule pathseq_extract_unpaired_bam:
        input:
            pathseq_classified_bam_file = os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_classified.bam"),
        output:
            pathseq_classified_unpaired_bam_file = temp(os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_unpaired_classified.bam"))
        params:
            threads=8
        conda:
            config["envs"]["star"]
        shell:
            '''
            samtools view --threads {params.threads} -h -b -F 1 -o {output.pathseq_classified_unpaired_bam_file} {input.pathseq_classified_bam_file}
            '''

    rule pathseq_score_cell_BAM:
        input:
            pathseq_classified_paired_sorted_bam_file = os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_sorted_paired_classified.bam"),
            pathseq_classified_unpaired_bam_file = os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_unpaired_classified.bam")
        output:
            pathseq_classified_score_output = os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_output.txt")
        params:
            taxonomy_db = config["params"]["classifier"]["pathseq"]["taxonomy_db"],
            pathseqscore_other_params = config["params"]["classifier"]["pathseqscore"] 
        resources:
            mem_mb=16000
        conda:
            config["envs"]["pathseq"]
        benchmark:
            os.path.join(config["benchmarks"]["classifier"],
                        "rmhost_pathseq_score/{sample}_pathseq_classifier_score_benchmark.log")
        log:
            os.path.join(config["logs"]["classifier"],
                        "rmhost_pathseq_score/{sample}_pathseq_classifier_score.log")
        shell:
            '''
            gatk PathSeqScoreSpark \
            --min-score-identity .7 \
            --unpaired-input {input.pathseq_classified_unpaired_bam_file} \
            --paired-input {input.pathseq_classified_paired_sorted_bam_file}\
            --taxonomy-file {params.taxonomy_db} \
            --scores-output {output.pathseq_classified_score_output} \
            --java-options "-Xmx15g -Xms15G" \
            --conf spark.port.maxRetries=64 \
            {params.pathseqscore_other_params}\
            2>&1 | tee {log}; \
            '''
    rule pathseq_INVADESEQ:
        input:
            unmapped_bam_sorted_file =os.path.join(
                config["output"]["host"],
                "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
            features_file = os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/{sample}_features.tsv"),
            pathseq_classified_bam_file = os.path.join(
                            config["output"]["classifier"],
                            "pathseq_classified_output/{sample}/{sample}_pathseq_classified.bam"),
            pathseq_output = os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_classified.txt")
        output:
            filtered_matrix_readname = os.path.join(
                config["output"]["classifier"],
                "pathseq_final_output/{sample}/{sample}_filtered_matrix_readname.txt"),
            unmap_cbub_bam = os.path.join(
                config["output"]["classifier"],
                "pathseq_final_output/{sample}/{sample}_pathseq_unmap_cbub.bam"),
            unmap_cbub_fasta = os.path.join(
                config["output"]["classifier"],
                "pathseq_final_output/{sample}/{sample}_pathseq_unmap_cbub.fasta"),
            filtered_matrix_list= os.path.join(
                config["output"]["classifier"],
                "pathseq_final_output/{sample}/{sample}_pathseq_filtered_matrix_list.txt"),
            matrix_readnamepath = os.path.join(
                    config["output"]["classifier"],
                    "pathseq_final_output/{sample}/{sample}_filtered_matrix.readnamepath"),
            genus_cell = os.path.join(
                    config["output"]["classifier"],
                    "pathseq_final_output/{sample}/{sample}_genus_cell.txt"),
            filtered_matrix_genus_csv = os.path.join(
                    config["output"]["classifier"],
                    "pathseq_final_output/{sample}/{sample}_filtered_matrix_genus.csv"),
            filtered_matrix_validate = os.path.join(
                    config["output"]["classifier"],
                    "pathseq_final_output/{sample}/{sample}_filtered_matrix.validate.csv")
        conda:
            config["envs"]["kmer_python"]
        params:
            SampleID="{sample}",
            INVADEseq_script = config["scripts"]["INVADEseq"]
        shell:
            '''
            python {params.INVADEseq_script} \
            {input.unmapped_bam_sorted_file} \
            {params.SampleID} \
            {input.features_file} \
            {input.pathseq_classified_bam_file}\
            {input.pathseq_output} \
            {output.filtered_matrix_readname} \
            {output.unmap_cbub_bam} \
            {output.unmap_cbub_fasta} \
            {output.filtered_matrix_list} \
            {output.matrix_readnamepath} \
            {output.genus_cell} \
            {output.filtered_matrix_genus_csv} \
            {output.filtered_matrix_validate}
            '''
    rule pathseq_classified_all:
        input:   
            # expand(os.path.join(
            #     config["output"]["classifier"],
            #     "pathseq_classified_output/{sample}/{sample}_pathseq_classified.bam"),sample=SAMPLES_ID_LIST),
            # expand(os.path.join(
            #     config["output"]["classifier"],
            #     "pathseq_classified_output/{sample}/{sample}_pathseq_classified.txt"),sample=SAMPLES_ID_LIST),
            # expand(os.path.join(
            #     config["output"]["classifier"],
            #     "pathseq_classified_output/{sample}/{sample}_pathseq_filter_metrics.txt"),sample=SAMPLES_ID_LIST),
            # expand(os.path.join(
            #     config["output"]["classifier"],
            #     "pathseq_classified_output/{sample}/{sample}_pathseq_score_metrics.txt"),sample=SAMPLES_ID_LIST)
            # expand(os.path.join(
            #     config["output"]["classifier"],
            #     "pathseq_classified_output/{sample}/{sample}_pathseq_output.txt"),sample=SAMPLES_ID_LIST)
            expand(os.path.join(
                config["output"]["classifier"],
                "pathseq_final_output/{sample}/{sample}_filtered_matrix_readname.txt"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                config["output"]["classifier"],
                "pathseq_final_output/{sample}/{sample}_pathseq_unmap_cbub.bam"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                config["output"]["classifier"],
                "pathseq_final_output/{sample}/{sample}_pathseq_unmap_cbub.fasta"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                config["output"]["classifier"],
                "pathseq_final_output/{sample}/{sample}_pathseq_filtered_matrix_list.txt"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                    config["output"]["classifier"],
                    "pathseq_final_output/{sample}/{sample}_filtered_matrix.readnamepath"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                    config["output"]["classifier"],
                    "pathseq_final_output/{sample}/{sample}_genus_cell.txt"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                    config["output"]["classifier"],
                    "pathseq_final_output/{sample}/{sample}_filtered_matrix_genus.csv"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                    config["output"]["classifier"],
                    "pathseq_final_output/{sample}/{sample}_filtered_matrix.validate.csv"),sample=SAMPLES_ID_LIST)
else:
    rule pathseq_classified_all:
        input:    

if config["params"]["classifier"]["metaphlan"]["do"]:
    rule metaphlan_classified:
        input:  
            unmapped_fastq = os.path.join(
                            config["output"]["host"],
                            "unmapped_host/{sample}/{sample}_unmappped2human_bam.fastq")
        output:
            mpa_bowtie2_out=os.path.join(
                config["output"]["classifier"],
                "metaphlan4_classified_output/{sample}/{sample}_metphlan4_classifier_bowtie2.bz2"),
            mpa_profile_out=os.path.join(
                config["output"]["classifier"],
                "metaphlan4_classified_output/{sample}/{sample}_metphlan4_classifier_profile.txt"),
        params:
            threads = config["params"]["classifier"]["metaphlan"]["threads"],
            sequence_type = config["params"]["classifier"]["metaphlan"]["sequence_type"],
            bowtie2db = config["params"]["classifier"]["metaphlan"]["bowtie2db"],
            db_index = config["params"]["classifier"]["metaphlan"]["db_index"],
            analysis_type = config["params"]["classifier"]["metaphlan"]["analysis_type"],
            metaphlan_other_params = config["params"]["classifier"]["metaphlan"]["metaphlan_other_params"] 
        benchmark:
            os.path.join(config["benchmarks"]["classifier"],
                        "rmhost_metaphlan_classifier/{sample}/{sample}_metaphalan4_classifier_benchmark.log")
        log:
            os.path.join(config["logs"]["classifier"],
                        "rmhost_metaphlan_classifier/{sample}/{sample}_metaphalan4_classifier.log")
        conda:
            config["envs"]["metaphlan"]
        resources:
            mem_mb=40000
        shell:
            '''
            metaphlan {input.unmapped_fastq} \
            -t {params.analysis_type} \
            --bowtie2out {output.mpa_bowtie2_out} \
            -o {output.mpa_profile_out} \
            --unclassified_estimation \
            --nproc {params.threads} \
            --input_type {params.sequence_type} \
            --bowtie2db {params.bowtie2db}  \
            --index {params.db_index} \
            {params.metaphlan_other_params}\
            2>&1 | tee {log}; \
            '''

    # rule mergeprofiles:
    #     input: 
    #         expand(os.path.join(
    #             config["output"]["classifier"],
    #             "metaphlan4_classified_output/{sample}/{sample}_metphlan4_classifier_profile.txt"), sample=SAMPLES_ID_LIST)
    #     output: 
    #         merged_abundance_table = os.path.join(
    #             config["output"]["classifier"],
    #             "metaphlan4_classified_output/merged_abundance_table.txt"),
    #         merged_species_abundance_table = os.path.join(
    #             config["output"]["classifier"],
    #             "metaphlan4_classified_output/merged_abundance_table_species.txt")
    #     params: 
    #         profiles=config["output_dir"]+"/metaphlan/*_profile.txt"
    #     conda: "utils/envs/metaphlan4.yaml"
    #     shell: """
    #         python utils/merge_metaphlan_tables.py {params.profiles} > {output.o1}
    #         grep -E "(s__)|(^ID)|(clade_name)|(UNKNOWN)|(UNCLASSIFIED)" {output.o1} | grep -v "t__"  > {output.o2}
    #         """
    rule metaphlan_classified_all:
        input:
            expand(os.path.join(
                config["output"]["classifier"],
                "metaphlan4_classified_output/{sample}/{sample}_metphlan4_classifier_bowtie2.bz2"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                config["output"]["classifier"],
                "metaphlan4_classified_output/{sample}/{sample}_metphlan4_classifier_profile.txt"),sample=SAMPLES_ID_LIST)
else:
    rule metaphlan_classified_all:
        input:    

rule classifier_all:
    input:
        rules.kraken2uniq_classified_all.input,
        rules.krakenuniq_classified_all.input,
        rules.pathseq_classified_all.input,
        rules.metaphlan_classified_all.input