import glob
def get_whitelist(wildcards):
    if config["params"]["begin"] == "host":
        if config["params"]["host"]["spaceranger"]["do"]:
            whitelist_file = os.path.join(
                config["output"]["host"],
                "spaceranger_count",wildcards.sample,f"{wildcards.sample}_barcodes.tsv")
            
            return whitelist_file

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
        config["envs"]["star"]
    shell:
        '''
        bash {params.bam2fastq_script} {input.unmapped_bam_sorted_file} {output.unmapped_r1_fastq} {output.unmapped_r2_fastq} {output.unmapped_fastq} {threads}
        '''

rule kraken2uniq_classified:
    input:
        unmapped_fastq = temp(os.path.join(
            config["output"]["host"],
            "unmapped_host/{sample}/{sample}_unmappped2human_bam.fastq")),
        unmapped_r1_fastq = temp(os.path.join(
            config["output"]["host"],
            "unmapped_host/{sample}/{sample}_unmappped2human_bam_r1.fastq")),
        unmapped_r2_fastq = temp(os.path.join(
            config["output"]["host"],
            "unmapped_host/{sample}/{sample}_unmappped2human_bam_r2.fastq"))
    output:
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
        minimum_hit = config["params"]["classifier"]["kraken2uniq"]["minimum_hit"],
        confidence = config["params"]["classifier"]["kraken2uniq"]["confidence"],
        kraken2mpa_script = config["scripts"]["kraken2mpa"],
        variousParams = config["params"]["classifier"]["kraken2uniq"]["variousParams"],
        #since kraken2 acquire specific input fomrat "#fq",so we put it it params
        # krak2_classified_output_fq_pair=os.path.join(
        #     config["output"]["classifier"],
        #     "classified_output/{sample}/{sample}_kraken2_classified#.fq"),
        # krak2_unclassified_output_fq_pair=os.path.join(
        #     config["output"]["classifier"],
        #     "unclassified_output/{sample}/{sample}_kraken2_unclassified#.fq"),
        # krak2_classified_output_fq = os.path.join(
        #     config["output"]["classifier"],
        #     "rmhost_classified_output/{sample}/{sample}_kraken2_classified.fq"),
        # krak2_unclassified_output_fq = os.path.join(
        #     config["output"]["classifier"],
        #     "rmhost_unclassified_output/{sample}/{sample}_kraken2_unclassified.fq"),
    resources:
        mem_mb=config["resources"]["kraken2uniq"]["mem_mb"]
    priority: 12
    threads: 
        config["resources"]["kraken2uniq"]["threads"]
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
        if [ -s "{input.unmapped_fastq}" ]; then
            kraken2 --db {params.database} \
            --threads {threads} \
            --output {output.krak2_output} \
            --report {output.krak2_report} \
            --report-minimizer-data \
            --minimum-hit-groups {params.minimum_hit} \
            --confidence {params.confidence} \
            {input.unmapped_fastq} \
            --use-names \
            {params.variousParams} \
            2>&1 | tee {log};\
        else
            kraken2 --db {params.database} \
            --threads {threads} \
            --output {output.krak2_output} \
            --report {output.krak2_report} \
            --minimum-hit-groups {params.minimum_hit}\
            --confidence {params.confidence} \
            --report-minimizer-data \
            {input.unmapped_r1_fastq} {input.unmapped_r2_fastq}\
            --use-names \
            --paired \
            {params.variousParams} \
            2>&1 | tee {log};\
        fi

        cut -f 1-3,6-8 {output.krak2_report} > {output.krak2_std_report};\
        python {params.kraken2mpa_script} -r {output.krak2_std_report} -o {output.krak2_mpa_report};
        '''
if config["params"]["host"]["starsolo"]["soloType"]=="SmartSeq":
    rule extract_kraken2_classified_bam:
        input:
            krak2_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_output/{sample}/{sample}_kraken2_output.txt"),
            krak2_report = os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_report/custom/{sample}/{sample}_kraken2_report.txt"),
            krak2_mpa_report=os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_report/mpa/{sample}/{sample}_kraken2_mpa_report.txt"),
            unmapped_bam_sorted_file =os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
        output:
            krak2_extracted_bam = os.path.join(
                config["output"]["classifier"],
                "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_extracted_classified.bam"),
            krak2_extracted_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_extracted_classified_output.txt"),                
        log:
            os.path.join(config["logs"]["classifier"],
                        "rmhost_kraken2uniq_extracted/{sample}_kraken2uniq_classifier_bam_extracted.log")
        params:
            extract_kraken_bam_script = config["scripts"]["extract_kraken_bam"],
            ktaxonomy_file = os.path.join(
                config["params"]["classifier"]["kraken2uniq"]["kraken2_database"],
                "ktaxonomy.tsv"),
            barcode_tag = ("CB") if PLATFORM == "lane" else "RG"
        resources:
            # mem_mb=config["resources"]["extract_kraken2_classified_bam"]["mem_mb"]
            mem_mb=(
                lambda wildcards: os.stat(os.path.join(
                config["output"]["host"],
                f"unmapped_host/{wildcards.sample}/Aligned_sortedByName_unmapped_out.bam")).st_size  /1024 /1024  * 1.5
            ),
        threads: 
            config["resources"]["extract_kraken2_classified_bam"]["threads"]
        priority: 
            14
        conda:
            config["envs"]["kmer_python"]
        shell:
            '''
            python {params.extract_kraken_bam_script} \
            --krak_output_file {input.krak2_output} \
            --kraken_report {input.krak2_report} \
            --ktaxonomy {params.ktaxonomy_file} \
            --extracted_bam_file {output.krak2_extracted_bam}\
            --input_bam_file {input.unmapped_bam_sorted_file} \
            --barcode_tag {params.barcode_tag} \
            --extracted_output_file {output.krak2_extracted_output} \
            --log_file {log}
            '''
else:
    rule extract_kraken2_classified_bam:
        input:
            krak2_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_output/{sample}/{sample}_kraken2_output.txt"),
            krak2_report = os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_report/custom/{sample}/{sample}_kraken2_report.txt"),
            krak2_mpa_report=os.path.join(
                config["output"]["classifier"],
                "rmhost_kraken2_report/mpa/{sample}/{sample}_kraken2_mpa_report.txt"),
            unmapped_bam_sorted_file =os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
            whitelist = get_whitelist
        output:
            krak2_extracted_bam = os.path.join(
                config["output"]["classifier"],
                "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_extracted_classified.bam"),
            krak2_extracted_output = os.path.join(
                config["output"]["classifier"],
                "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_extracted_classified_output.txt"),                
        log:
            os.path.join(config["logs"]["classifier"],
                        "rmhost_kraken2uniq_extracted/{sample}_kraken2uniq_classifier_bam_extracted.log")
        params:
            extract_kraken_bam_script = config["scripts"]["extract_kraken_bam"],
            ktaxonomy_file = os.path.join(
                config["params"]["classifier"]["kraken2uniq"]["kraken2_database"],
                "ktaxonomy.tsv"),
            barcode_tag = ("CB") if PLATFORM == "lane" else "RG"
        resources:
            # mem_mb=config["resources"]["extract_kraken2_classified_bam"]["mem_mb"]
            mem_mb=(
                lambda wildcards: os.stat(os.path.join(
                config["output"]["host"],
                f"unmapped_host/{wildcards.sample}/Aligned_sortedByName_unmapped_out.bam")).st_size  /1024 /1024  * 1.5
            ),
        threads: 
            config["resources"]["extract_kraken2_classified_bam"]["threads"]
        priority: 
            14
        conda:
            config["envs"]["kmer_python"]
        shell:
            '''
            python {params.extract_kraken_bam_script} \
            --krak_output_file {input.krak2_output} \
            --kraken_report {input.krak2_report} \
            --ktaxonomy {params.ktaxonomy_file} \
            --extracted_bam_file {output.krak2_extracted_bam}\
            --input_bam_file {input.unmapped_bam_sorted_file} \
            --barcode_tag {params.barcode_tag} \
            --whitelist {input.whitelist} \
            --extracted_output_file {output.krak2_extracted_output} \
            --log_file {log}
            '''
rule krak_sample_denosing:
    input:
        krak2_extracted_output = os.path.join(
            config["output"]["classifier"],
            "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_extracted_classified_output.txt"),        
        krak2_report = os.path.join(
            config["output"]["classifier"],
            "rmhost_kraken2_report/custom/{sample}/{sample}_kraken2_report.txt"),
    output:
        krak_sample_denosing_result = os.path.join(
            config["output"]["classifier"],
            "rmhost_kraken2_qc/{sample}/{sample}_krak_sample_denosing.txt"),
        krak_sample_raw_result = os.path.join(
            config["output"]["classifier"],
            "rmhost_kraken2_qc/{sample}/{sample}_krak_sample_raw.txt"),
    resources:
        mem_mb=(
            lambda wildcards, input: os.stat(input.krak2_extracted_output).st_size  /1024 /1024  * 2
        ),
    threads: 
        config["resources"]["krak_sample_denosing"]["threads"]
    log:
        os.path.join(config["logs"]["classifier"],
                    "classified_qc/{sample}/{sample}_krak_sample_denosing.log")
    priority: 
        15
    params:
        krak_sample_denosing_script= config["scripts"]["krak_sample_denosing"],
        min_read_fraction = 0.15,
        inspect_file = os.path.join(config["params"]["classifier"]["kraken2uniq"]["kraken2_database"],"inspect.txt"),
        ktaxonomy_file = os.path.join(config["params"]["classifier"]["kraken2uniq"]         ["kraken2_database"],"ktaxonomy.tsv"),
        barcode_tag = ("CB") if PLATFORM == "lane" else "RG"
    conda:
        config["envs"]["kmer_python"]
    benchmark:
        os.path.join(config["benchmarks"]["classifier"],
                    "rmhost_kraken2_qc/{sample}_sample_denosing_benchmark.tsv")
    shell:
        '''
        python {params.krak_sample_denosing_script} \
        --krak_report {input.krak2_report} \
        --krak_output {input.krak2_extracted_output} \
        --ktaxonomy {params.ktaxonomy_file}\
        --inspect {params.inspect_file} \
        --min_read_fraction {params.min_read_fraction} \
        --qc_output_file {output.krak_sample_denosing_result} \
        --raw_qc_output_file {output.krak_sample_raw_result} \
        --log_file {log};
        '''
rule krak_study_denosing:
    input:
        krak_sample_denosing_result_list = expand(os.path.join(
            config["output"]["classifier"],
            "rmhost_kraken2_qc/{sample}/{sample}_krak_sample_denosing.txt"),sample=SAMPLES_ID_LIST),
        krak_sample_raw_result_list = expand(os.path.join(
            config["output"]["classifier"],
            "rmhost_kraken2_qc/{sample}/{sample}_krak_sample_raw.txt"),sample=SAMPLES_ID_LIST)
    output:
        candidate_species =  os.path.join(
            config["output"]["classifier"],
            "rmhost_kraken2_qc/study/krak_candidate_species.txt"),
    priority: 
        15
    log:
        os.path.join(config["logs"]["classifier"],
                    "classified_qc/study/krak_candidate_species.log")
    params:
        krak_sample_denosing_output_dir = os.path.join(
            config["output"]["classifier"],
            "rmhost_kraken2_qc/"),
        min_reads = config["params"]["classifier"]["krak_study_denosing"]["min_reads"],
        min_uniq = config["params"]["classifier"]["krak_study_denosing"]["min_uniq"],
        krak_study_denosing_script= config["scripts"]["krak_study_denosing"]
    conda:
        config["envs"]["kmer_python"]
    shell:
        '''
        python  {params.krak_study_denosing_script}\
        --file_list {input.krak_sample_denosing_result_list} \
        --out_path {output.candidate_species} \
        --raw_file_list {input.krak_sample_raw_result_list} \
        --min_reads {params.min_reads} \
        --min_uniq {params.min_uniq} \
        --log_file {log}
        '''
rule extract_kraken2_denosied_fastq:
    input:
        krak2_output = os.path.join(
            config["output"]["classifier"],
            "rmhost_kraken2_output/{sample}/{sample}_kraken2_output.txt"),
        unmapped_bam_sorted_file =os.path.join(
            config["output"]["host"],
            "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
        candidate_species =  os.path.join(
                        config["output"]["classifier"],
                        "rmhost_kraken2_qc/study/krak_candidate_species.txt")
    output:
        krak_screened_fastq = temp(os.path.join(
            config["output"]["classifier"],
            "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_screen.fastq")),
        krak_screened_r1_fastq = temp(os.path.join(
            config["output"]["classifier"],
            "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_screen_r1.fastq")),
        krak_screened_r2_fastq = temp(os.path.join(
            config["output"]["classifier"],
            "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_screen_r2.fastq"))
    log:
        os.path.join(config["logs"]["classifier"],
                    "rmhost_kraken2uniq_extracted/{sample}_kraken2uniq_classifier_fastq_extracted.log")
    params:
        extract_kraken_fastq_script = config["scripts"]["extract_kraken_fastq"],
        SampleID="{sample}",
        ktaxonomy_file = os.path.join(
            config["params"]["classifier"]["kraken2uniq"]["kraken2_database"],
            "ktaxonomy.tsv"),
    benchmark:
        os.path.join(config["benchmarks"]["classifier"],
                    "rmhost_kraken2_qc/{sample}_sample_extraced_denosing_benchmark.tsv")
    resources:
        mem_mb=config["resources"]["extract_kraken2_classified_bam"]["mem_mb"]
    threads: 
        config["resources"]["extract_kraken2_classified_bam"]["threads"]
    priority: 
        14
    conda:
        config["envs"]["kmer_python"]
    shell:
        '''
        python {params.extract_kraken_fastq_script} \
        --candidate {input.candidate_species} \
        --krak_output_file {input.krak2_output} \
        --sample_name {params.SampleID} \
        --ktaxonomy {params.ktaxonomy_file} \
        --fastq_r1 {output.krak_screened_r1_fastq}\
        --fastq_r2 {output.krak_screened_r2_fastq}\
        --fastq {output.krak_screened_fastq}\
        --input_bam_file {input.unmapped_bam_sorted_file} \
        --processes {threads} \
        --log_file {log}
        '''
rule kraken2uniq_classified_all:
    input:
        expand(os.path.join(
            config["output"]["classifier"],
            "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_screen.fastq"),sample=SAMPLES_ID_LIST),
        expand(os.path.join(
            config["output"]["classifier"],
            "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_screen_r1.fastq"),sample=SAMPLES_ID_LIST),
        expand(os.path.join(
            config["output"]["classifier"],
            "rmhost_extracted_classified_output/{sample}/{sample}_kraken2_screen_r2.fastq"),sample=SAMPLES_ID_LIST)

rule classifier_all:
    input:
        rules.kraken2uniq_classified_all.input,