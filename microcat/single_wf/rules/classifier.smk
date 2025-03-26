import glob
def get_whitelist(wildcards):
    if config["params"]["begin"] == "host":
        if config["params"]["host"]["starsolo"]["do"]:
            if config["params"]["host"]["starsolo"]["soloType"]=="CB_UMI_Complex" or config["params"]["host"]["starsolo"]["soloType"]=="CB_UMI_Simple":
                # whitelist_file = os.path.join(
                #                 config["output"]["host"],
                #                 "starsolo_count/",wildcards.sample,"Solo.out/Gene/filtered/barcodes.tsv.gz")
                whitelist_file = os.path.join(
                    config["output"]["host"],
                    "starsolo_count/",wildcards.sample,f"{wildcards.sample}_barcodes.tsv")
                # if os.path.exists(whitelist_file):
                return whitelist_file
                # else:
                #     whitelist_file = os.path.join(
                #                 config["output"]["host"],
                #                 "starsolo_count/",wildcards.sample,"Solo.out/Gene/filtered/barcodes.tsv")
                #     if os.path.exists(whitelist_file):
                #         return whitelist_file
                #     else:
                #         raise(ValueError("No whitelist file found for sample %s" % wildcards.sample))

            elif config["params"]["host"]["starsolo"]["soloType"]=="SmartSeq":
                # whitelist_file = os.path.join(
                #                 config["output"]["host"],
                #                 "starsolo_count_SE/",wildcards.sample,"Solo.out/Gene/raw/barcodes.tsv.gz")
                whitelist_file = glob.glob(os.path.join(
                                    config["output"]["host"],
                                    "starsolo_count_SE/", wildcards.sample, "Solo.out/Gene/raw/barcodes.tsv*"))
                if whitelist_file:
                    return whitelist_file[0]
                else:
                    # whitelist_file = os.path.join(
                    #                 config["output"]["host"],
                    #                 "starsolo_count_PE/",wildcards.sample,"Solo.out/Gene/raw/barcodes.tsv.gz")
                    whitelist_file = glob.glob(os.path.join(
                                        config["output"]["host"],
                                        "starsolo_count_PE/", wildcards.sample, "Solo.out/Gene/raw/barcodes.tsv*"))
                    if whitelist_file:
                        return whitelist_file[0]
                
                raise(ValueError("No whitelist file found for sample %s" % wildcards.sample))
        elif config["params"]["host"]["cellranger"]["do"]:
            # whitelist_file = os.path.join(
            #                 config["output"]["host"],
            #                 "cellranger_count/",wildcards.sample,"outs/raw_feature_bc_matrix/barcodes.tsv.gz")
            whitelist_file = os.path.join(
                config["output"]["host"],
                "cellranger_count",wildcards.sample,f"{wildcards.sample}_barcodes.tsv")
            
            return whitelist_file
            # if os.path.exists(whitelist_file):
            #     return whitelist_file
            # else:
            #     whitelist_file = os.path.join(
            #                 config["output"]["host"],
            #                 "cellranger_count/",wildcards.sample,"outs/raw_feature_bc_matrix/barcodes.tsv")
            #     if os.path.exists(whitelist_file):
            #         return whitelist_file
            #     else:
            #         raise(ValueError("No whitelist file found for sample %s" % wildcards.sample))
    if config["params"]["begin"] == "classifier":
        mtx_path = microcat.get_starsolo_sample_id(SAMPLES, wildcards, "mtx")
        whitelist_file = glob.glob(os.path.join(mtx_path, "barcodes.tsv*"))
        if whitelist_file is not None:
            return whitelist_file[0]
        else:
            raise(ValueError("No whitelist file found for sample %s" % wildcards.sample))
        # if os.path.exists(whitelist_file):
        #     return whitelist_file
    raise(ValueError("No whitelist file found for sample %s" % wildcards.sample))

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
# rule bwa_remove_host:
#     input:
#         unmapped_fastq = os.path.join(
#             config["output"]["host"],
#             "unmapped_host/{sample}/{sample}_unmappped2human_bam.fastq"),
#         unmapped_r1_fastq = os.path.join(
#             config["output"]["host"],
#             "unmapped_host/{sample}/{sample}_unmappped2human_bam_r1.fastq"),
#         unmapped_r2_fastq = os.path.join(
#             config["output"]["host"],
#             "unmapped_host/{sample}/{sample}_unmappped2human_bam_r2.fastq")
#     output:
#         rmhost_fastq = os.path.join(
#             config["output"]["host"],
#             "unmapped_host/{sample}/{sample}_bwa_rmhost.fastq"),
#         rmhost_r1_fastq = os.path.join(
#             config["output"]["host"],
#             "unmapped_host/{sample}/{sample}_bwa_rmhost_r1.fastq"),
#         rmhost_r2_fastq = os.path.join(
#             config["output"]["host"],
#             "unmapped_host/{sample}/{sample}_bwa_rmhost_r2.fastq"),
#     params:
#         host_index = config["params"]["align"]["bwa2"]["host_db"],
#         rmhost_bam = os.path.join(
#             config["output"]["host"],
#             "unmapped_host/{sample}/{sample}_rmhost.bam")
#     log:
#         os.path.join(config["logs"]["classifier"],
#                     "bwa_rmhost/{sample}_rmhost.log")
#     threads:
#         config["resources"]["rmhost"]["threads"]
#     resources:
#         mem_mb=config["resources"]["rmhost"]["mem_mb"]
#     priority: 11
#     conda:
#         config["envs"]["star"]
#     shell:
#         '''
#         if [ -s "{input.unmapped_fastq}" ]; then
#             bwa-mem2 mem -t {threads} {params.host_index} {input.unmapped_fastq} | samtools view --threads {threads} -b -f 4 > {params.rmhost_bam} ;\
#             samtools fastq --threads {threads} -n {params.rmhost_bam} > {output.rmhost_fastq};  \
#             touch {output.rmhost_r1_fastq} ;\
#             touch {output.rmhost_r2_fastq} ;
#         else
#             bwa-mem2 mem -t {threads} -p {params.host_index} {input.unmapped_r1_fastq} {input.unmapped_r2_fastq} | samtools view --threads {threads} -b -f12 > {params.rmhost_bam} ;\
#             samtools fastq --threads {threads}  -n {params.rmhost_bam} -1 {output.rmhost_r1_fastq} -2 {output.rmhost_r2_fastq}; \
#             touch {output.rmhost_fastq};
#         fi
#         '''
if config["params"]["classifier"]["kraken2uniq"]["do"]:
    rule kraken2uniq_classified:
        input:
            # rmhost_fastq = os.path.join(
            #     config["output"]["host"],
            #     "unmapped_host/{sample}/{sample}_bwa_rmhost.fastq"),
            # rmhost_r1_fastq = os.path.join(
            #     config["output"]["host"],
            #     "unmapped_host/{sample}/{sample}_bwa_rmhost_r1.fastq"),
            # rmhost_r2_fastq = os.path.join(
            #     config["output"]["host"],
            #     "unmapped_host/{sample}/{sample}_bwa_rmhost_r2.fastq"),
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
            estimate_precision=config["params"]["classifier"]["krakenuniq"]["estimate_precision"],
            variousParams = config["params"]["classifier"]["krakenuniq"]["variousParams"],
        benchmark:
            os.path.join(config["benchmarks"]["classifier"],
                        "rmhost_krakenuniq/{sample}_kraken2uniq_classifier_benchmark.tsv")
        log:
            os.path.join(config["logs"]["classifier"],
                        "rmhost_krakenuniq/{sample}_krakenuniq_classifier.log")
        threads:
            config["resources"]["krakenuniq"]["threads"]
        resources:
            mem_mb=config["resources"]["krakenuniq"]["mem_mb"]
        conda:
            config["envs"]["krakenuniq"]
        # message:
        #     "classifier: Performing Taxonomic Classifcation of Sample {sample} with krakenuniq."
        shell:
            '''
            krakenuniq --db {params.database} \
            --threads {threads} \
            --hll-precision {params.estimate_precision} \
            --classified-out {output.krakenuniq_classified_output}\
            --unclassified-out {output.krakenuniq_unclassified_output}\
            --output {output.krakenuniq_output} \
            --report-file {output.krakenuniq_report} \
            {input.unmapped_fastq}  \
            --preload \
            {params.variousParams} \
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
                "pathseq_classified_output/{sample}/"),
            variousParams = config["params"]["classifier"]["pathseq"]["variousParams"],
        resources:
            mem_mb=config["resources"]["pathseq"]["mem_mb"]
        priority: 12
        threads: 
            config["resources"]["pathseq"]["threads"]
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
            {params.variousParams} \
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
        threads:
            8
        conda:
            config["envs"]["star"]
        shell:
            '''
            samtools view --threads {threads} -h -b -f 1 -o {output.pathseq_classified_paired_bam_file} {input.pathseq_classified_bam_file}
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
        threads:
            8
        conda:
            config["envs"]["star"]
        shell:
            '''
            samtools sort --threads {threads} -n -o {output.pathseq_classified_paired_sorted_bam_file} {input.pathseq_classified_paired_bam_file} 
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
        threads:
            8
        resources:
            mem_mb=config["resources"]["samtools_extract"]["mem_mb"]
        conda:
            config["envs"]["star"]
        shell:
            '''
            samtools view --threads {threads} -h -b -F 1 -o {output.pathseq_classified_unpaired_bam_file} {input.pathseq_classified_bam_file}
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
    # rule pathseq_INVADESEQ:
    #     input:
    #         unmapped_bam_sorted_file =os.path.join(
    #             config["output"]["host"],
    #             "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
    #         features_file = os.path.join(
    #             config["output"]["host"],
    #             "cellranger_count/{sample}/{sample}_features.tsv"),
    #         pathseq_classified_bam_file = os.path.join(
    #                         config["output"]["classifier"],
    #                         "pathseq_classified_output/{sample}/{sample}_pathseq_classified.bam"),
    #         pathseq_output = os.path.join(
    #             config["output"]["classifier"],
    #             "pathseq_classified_output/{sample}/{sample}_pathseq_classified.txt")
    #     output:
    #         filtered_matrix_readname = os.path.join(
    #             config["output"]["classifier"],
    #             "pathseq_final_output/{sample}/{sample}_filtered_matrix_readname.txt"),
    #         unmap_cbub_bam = os.path.join(
    #             config["output"]["classifier"],
    #             "pathseq_final_output/{sample}/{sample}_pathseq_unmap_cbub.bam"),
    #         unmap_cbub_fasta = os.path.join(
    #             config["output"]["classifier"],
    #             "pathseq_final_output/{sample}/{sample}_pathseq_unmap_cbub.fasta"),
    #         filtered_matrix_list= os.path.join(
    #             config["output"]["classifier"],
    #             "pathseq_final_output/{sample}/{sample}_pathseq_filtered_matrix_list.txt"),
    #         matrix_readnamepath = os.path.join(
    #                 config["output"]["classifier"],
    #                 "pathseq_final_output/{sample}/{sample}_filtered_matrix.readnamepath"),
    #         genus_cell = os.path.join(
    #                 config["output"]["classifier"],
    #                 "pathseq_final_output/{sample}/{sample}_genus_cell.txt"),
    #         filtered_matrix_genus_csv = os.path.join(
    #                 config["output"]["classifier"],
    #                 "pathseq_final_output/{sample}/{sample}_filtered_matrix_genus.csv"),
    #         filtered_matrix_validate = os.path.join(
    #                 config["output"]["classifier"],
    #                 "pathseq_final_output/{sample}/{sample}_filtered_matrix.validate.csv")
    #     conda:
    #         config["envs"]["kmer_python"]
    #     params:
    #         SampleID="{sample}",
    #         INVADEseq_script = config["scripts"]["INVADEseq"]
    #     shell:
    #         '''
    #         python {params.INVADEseq_script} \
    #         {input.unmapped_bam_sorted_file} \
    #         {params.SampleID} \
    #         {input.features_file} \
    #         {input.pathseq_classified_bam_file}\
    #         {input.pathseq_output} \
    #         {output.filtered_matrix_readname} \
    #         {output.unmap_cbub_bam} \
    #         {output.unmap_cbub_fasta} \
    #         {output.filtered_matrix_list} \
    #         {output.matrix_readnamepath} \
    #         {output.genus_cell} \
    #         {output.filtered_matrix_genus_csv} \
    #         {output.filtered_matrix_validate}
    #         '''
    rule pathseq_classified_all:
        input:   
            expand(os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_classified.bam"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_classified.txt"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_filter_metrics.txt"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_score_metrics.txt"),sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                config["output"]["classifier"],
                "pathseq_classified_output/{sample}/{sample}_pathseq_output.txt"),sample=SAMPLES_ID_LIST)
            # expand(os.path.join(
            #     config["output"]["classifier"],
            #     "pathseq_final_output/{sample}/{sample}_filtered_matrix_readname.txt"),sample=SAMPLES_ID_LIST),
            # expand(os.path.join(
            #     config["output"]["classifier"],
            #     "pathseq_final_output/{sample}/{sample}_pathseq_unmap_cbub.bam"),sample=SAMPLES_ID_LIST),
            # expand(os.path.join(
            #     config["output"]["classifier"],
            #     "pathseq_final_output/{sample}/{sample}_pathseq_unmap_cbub.fasta"),sample=SAMPLES_ID_LIST),
            # expand(os.path.join(
            #     config["output"]["classifier"],
            #     "pathseq_final_output/{sample}/{sample}_pathseq_filtered_matrix_list.txt"),sample=SAMPLES_ID_LIST),
            # expand(os.path.join(
            #         config["output"]["classifier"],
            #         "pathseq_final_output/{sample}/{sample}_filtered_matrix.readnamepath"),sample=SAMPLES_ID_LIST),
            # expand(os.path.join(
            #         config["output"]["classifier"],
            #         "pathseq_final_output/{sample}/{sample}_genus_cell.txt"),sample=SAMPLES_ID_LIST),
            # expand(os.path.join(
            #         config["output"]["classifier"],
            #         "pathseq_final_output/{sample}/{sample}_filtered_matrix_genus.csv"),sample=SAMPLES_ID_LIST),
            # expand(os.path.join(
            #         config["output"]["classifier"],
            #         "pathseq_final_output/{sample}/{sample}_filtered_matrix.validate.csv"),sample=SAMPLES_ID_LIST)
else:
    rule pathseq_classified_all:
        input:    

if config["params"]["classifier"]["metaphlan4"]["do"]:
    rule metaphlan_classified:
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
            mpa_bowtie2_out=os.path.join(
                config["output"]["classifier"],
                "metaphlan4_classified_output/{sample}/{sample}_metphlan4_classifier_bowtie2.bz2"),
            mpa_profile_out=os.path.join(
                config["output"]["classifier"],
                "metaphlan4_classified_output/{sample}/{sample}_metphlan4_classifier_profile.txt"),
        params:
            sequence_type = config["params"]["classifier"]["metaphlan4"]["sequence_type"],
            bowtie2db = config["params"]["classifier"]["metaphlan4"]["bowtie2db"],
            db_index = config["params"]["classifier"]["metaphlan4"]["db_index"],
            analysis_type = config["params"]["classifier"]["metaphlan4"]["analysis_type"],
            variousParams = config["params"]["classifier"]["metaphlan4"]["variousParams"] 
        threads:
            config["resources"]["metaphlan4"]["threads"]
        benchmark:
            os.path.join(config["benchmarks"]["classifier"],
                        "rmhost_metaphlan_classifier/{sample}/{sample}_metaphalan4_classifier_benchmark.log")
        log:
            os.path.join(config["logs"]["classifier"],
                        "rmhost_metaphlan_classifier/{sample}/{sample}_metaphalan4_classifier.log")
        conda:
            config["envs"]["metaphlan"]
        resources:
            mem_mb = config["resources"]["metaphlan4"]["mem_mb"]
        shell:
            '''
            if [ -s "{input.unmapped_fastq}" ]; then
                metaphlan {input.unmapped_fastq} \
                -t {params.analysis_type} \
                --bowtie2out {output.mpa_bowtie2_out} \
                -o {output.mpa_profile_out} \
                --unclassified_estimation \
                --nproc {threads} \
                --input_type {params.sequence_type} \
                --bowtie2db {params.bowtie2db}  \
                --index {params.db_index} \
                {params.variousParams}\
                2>&1 | tee {log}; \
            else
                metaphlan {input.unmapped_r1_fastq} {input.unmapped_r2_fastq} \
                -t {params.analysis_type} \
                --bowtie2out {output.mpa_bowtie2_out} \
                -o {output.mpa_profile_out} \
                --unclassified_estimation \
                --nproc {threads} \
                --input_type {params.sequence_type} \
                --bowtie2db {params.bowtie2db}  \
                --index {params.db_index} \
                {params.variousParams}\
                2>&1 | tee {log}; \
            fi
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
        # rules.kraken2uniq_classified_all.input,
        rules.krakenuniq_classified_all.input,
        rules.pathseq_classified_all.input,
        rules.metaphlan_classified_all.input