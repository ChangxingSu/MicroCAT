localrules: generate_se_manifest_file, generate_pe_manifest_file

if config["params"]["host"]["starsolo"]["do"]:
    if config["params"]["host"]["starsolo"]["soloType"]=="CB_UMI_Simple":
        # Auto detect 10x Genomics
        if config["params"]["host"]["starsolo"]["name"]=="tenX_AUTO":
            if config["params"]["host"]["starsolo"]["variousParams"]:
                rule starsolo_CB_UMI_Simple_count:
                    input:
                        # Directory containing input fastq files
                        fastqs_dir=lambda wildcards: microcat.get_fastqs_dir(SAMPLES,wildcards),
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
                            "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam"),
                        starsolo_count_report = os.path.join(
                            config["output"]["host"],
                            "starsolo_count/{sample}/Log.final.out")
                    params:
                        barcode_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                        cdna_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq2"),
                        starsolo_out = os.path.join(
                            config["output"]["host"],
                            "starsolo_count/"),
                        starsolo_10X_auto = config["scripts"]["starsolo_10X_auto"],
                        reference = config["params"]["host"]["starsolo"]["reference"],
                        barcode_data_dir = config["datas"]["barcode_list_dirs"]["tenX"],
                        variousParams = config["params"]["host"]["starsolo"]["variousParams"],
                        outSAMattributes = config["params"]["host"]["starsolo"]["outSAMattributes"],
                        soloUMIdedup = config["params"]["host"]["starsolo"]["algorithm"]["soloUMIdedup"],
                        soloCBmatchWLtype = config["params"]["host"]["starsolo"]["algorithm"]["soloCBmatchWLtype"],
                        soloUMIfiltering = config["params"]["host"]["starsolo"]["algorithm"]["soloUMIfiltering"],
                        soloCellFilter = config["params"]["host"]["starsolo"]["algorithm"]["soloCellFilter"],
                        clipAdapterType = config["params"]["host"]["starsolo"]["algorithm"]["clipAdapterType"],
                        outFilterScoreMin = config["params"]["host"]["starsolo"]["algorithm"]["outFilterScoreMin"],
                        soloMultiMappers = config["params"]["host"]["starsolo"]["algorithm"]["soloMultiMappers"],
                    log:
                        os.path.join(config["logs"]["host"],
                                    "starsolo/{sample}_starsolo_count.log")
                    benchmark:
                        os.path.join(config["benchmarks"]["host"],
                                    "starsolo/{sample}_starsolo_count.benchmark")
                    conda:
                        config["envs"]["star"]
                    resources:
                        mem_mb=config["resources"]["starsolo"]["mem_mb"]
                    threads:
                        config["resources"]["starsolo"]["threads"]
                    shell:
                        '''
                            mkdir -p {params.starsolo_out}; 
                            cd {params.starsolo_out} ;
                            
                            bash {params.starsolo_10X_auto} \
                            --barcode_reads {params.barcode_reads} \
                            --cdna_reads {params.cdna_reads} \
                            --barcode_data_dir {params.barcode_data_dir} \
                            --sample {wildcards.sample} \
                            --threads {threads} \
                            --reference {params.reference} \
                            --soloUMIdedup {params.soloUMIdedup} \
                            --soloCBmatchWLtype {params.soloCBmatchWLtype} \
                            --soloUMIfiltering {params.soloUMIfiltering} \
                            --soloCellFilter {params.soloCellFilter} \
                            --outFilterScoreMin {params.outFilterScoreMin} \
                            --soloMultiMappers {params.soloMultiMappers} \
                            --clipAdapterType {params.clipAdapterType}\
                            --variousParams '{params.variousParams}'\
                            2>&1 | tee ../../../{log} ;
                        '''
            else:
                rule starsolo_CB_UMI_Simple_count:
                    input:
                        # Directory containing input fastq files
                        fastqs_dir=lambda wildcards: microcat.get_fastqs_dir(SAMPLES,wildcards),
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
                            "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam"),
                        starsolo_count_report = os.path.join(
                            config["output"]["host"],
                            "starsolo_count/{sample}/Log.final.out")
                    params:
                        barcode_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                        cdna_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq2"),
                        starsolo_out = os.path.join(
                            config["output"]["host"],
                            "starsolo_count/"),
                        starsolo_10X_auto = config["scripts"]["starsolo_10X_auto"],
                        reference = config["params"]["host"]["starsolo"]["reference"],
                        barcode_data_dir = config["datas"]["barcode_list_dirs"]["tenX"],
                        outSAMattributes = config["params"]["host"]["starsolo"]["outSAMattributes"],
                        soloUMIdedup = config["params"]["host"]["starsolo"]["algorithm"]["soloUMIdedup"],
                        soloCBmatchWLtype = config["params"]["host"]["starsolo"]["algorithm"]["soloCBmatchWLtype"],
                        soloUMIfiltering = config["params"]["host"]["starsolo"]["algorithm"]["soloUMIfiltering"],
                        soloCellFilter = config["params"]["host"]["starsolo"]["algorithm"]["soloCellFilter"],
                        clipAdapterType = config["params"]["host"]["starsolo"]["algorithm"]["clipAdapterType"],
                        outFilterScoreMin = config["params"]["host"]["starsolo"]["algorithm"]["outFilterScoreMin"],
                        soloMultiMappers = config["params"]["host"]["starsolo"]["algorithm"]["soloMultiMappers"],
                    log:
                        os.path.join(config["logs"]["host"],
                                    "{sample}/{sample}_starsolo_count.log")
                    benchmark:
                        os.path.join(config["benchmarks"]["host"],
                                    "{sample}/{sample}_starsolo_count.benchmark")
                    conda:
                        config["envs"]["star"]
                    resources:
                        mem_mb=config["resources"]["starsolo"]["mem_mb"]
                    threads:
                        config["resources"]["starsolo"]["threads"]
                    shell:
                        '''
                            mkdir -p {params.starsolo_out}; 
                            cd {params.starsolo_out} ;
                            
                            bash {params.starsolo_10X_auto} \
                            --barcode_reads {params.barcode_reads} \
                            --cdna_reads {params.cdna_reads} \
                            --barcode_data_dir {params.barcode_data_dir} \
                            --sample {wildcards.sample} \
                            --threads {threads} \
                            --reference {params.reference} \
                            --soloUMIdedup {params.soloUMIdedup} \
                            --soloCBmatchWLtype {params.soloCBmatchWLtype} \
                            --soloUMIfiltering {params.soloUMIfiltering} \
                            --soloCellFilter {params.soloCellFilter} \
                            --outFilterScoreMin {params.outFilterScoreMin} \
                            --soloMultiMappers {params.soloMultiMappers} \
                            --clipAdapterType {params.clipAdapterType}\
                            2>&1 | tee ../../../{log} ;
                        '''
        else:
            rule starsolo_CB_UMI_Simple_count:
                input:
                    # Directory containing input fastq files
                    fastqs_dir=lambda wildcards: microcat.get_fastqs_dir(SAMPLES,wildcards),
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
                        "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam"),
                    starsolo_count_report = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/Log.final.out")
                params:
                    barcode_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                    cdna_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq2"),
                    starsolo_out = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/"),
                    reference = config["params"]["host"]["starsolo"]["reference"],
                    chemistry = config["params"]["host"]["starsolo"]["name"],
                    description = config["params"]["host"]["starsolo"]["description"],
                    soloCBlen = config["params"]["host"]["starsolo"]["barcode"]["soloCBlen"],
                    soloCBstart = config["params"]["host"]["starsolo"]["barcode"]["soloCBstart"],
                    soloBarcodeMate = config["params"]["host"]["starsolo"]["barcode"]["soloBarcodeMate"],
                    soloUMIstart = config["params"]["host"]["starsolo"]["barcode"]["soloUMIstart"],
                    soloUMIlen = config["params"]["host"]["starsolo"]["barcode"]["soloUMIlen"],
                    clip5pNbases = config["params"]["host"]["starsolo"]["barcode"]["clip5pNbases"],
                    clip3pNbases = config["params"]["host"]["starsolo"]["barcode"]["clip3pNbases"],
                    soloUMIdedup = config["params"]["host"]["starsolo"]["algorithm"]["soloUMIdedup"],
                    soloCBmatchWLtype = config["params"]["host"]["starsolo"]["algorithm"]["soloCBmatchWLtype"],
                    soloUMIfiltering = config["params"]["host"]["starsolo"]["algorithm"]["soloUMIfiltering"],
                    soloCellFilter = config["params"]["host"]["starsolo"]["algorithm"]["soloCellFilter"],
                    clipAdapterType = config["params"]["host"]["starsolo"]["algorithm"]["clipAdapterType"],
                    outFilterScoreMin = config["params"]["host"]["starsolo"]["algorithm"]["outFilterScoreMin"],
                    soloMultiMappers = config["params"]["host"]["starsolo"]["algorithm"]["soloMultiMappers"],
                    barcode_list = "None" if config["params"]["host"]["starsolo"]["barcode"]["soloCBwhitelist"] == "None"
                                            else os.path.join(config["datas"]["barcode_list_dirs"]["tenX"],
                                                config["params"]["host"]["starsolo"]["barcode"]["soloCBwhitelist"]),
                    outSAMattributes = config["params"]["host"]["starsolo"]["outSAMattributes"],
                    outSAMtype = config["params"]["host"]["starsolo"]["outSAMtype"],
                    variousParams = config["params"]["host"]["starsolo"]["variousParams"],
                resources:
                    mem_mb=config["resources"]["starsolo"]["mem_mb"]
                threads:
                    config["resources"]["starsolo"]["threads"]
                log:
                    os.path.join(config["logs"]["host"],
                                "{sample}/{sample}_starsolo_count.log")
                benchmark:
                    os.path.join(config["benchmarks"]["host"],
                                "{sample}/{sample}_starsolo_count.benchmark")
                conda:
                    config["envs"]["star"]
                message: "Executing starsolo with {threads} threads on the following files {wildcards.sample}.Library with {params.description}"
                shell:
                    '''
                    if echo {params.cdna_reads} | grep -q "\.gz" ; then
                        file_command='--readFilesCommand zcat'
                    else
                        file_command=''
                    fi

                    mkdir -p {params.starsolo_out}; 
                    cd {params.starsolo_out} ;
                    STAR \
                    $file_command  \
                    --soloType CB_UMI_Simple \
                    --soloCBwhitelist {params.barcode_list} --soloCBstart {params.soloCBstart} --soloCBlen {params.soloCBlen} \
                    --soloUMIstart {params.soloUMIstart} --soloUMIlen {params.soloUMIlen} \
                    --genomeDir {params.reference} \
                    --readFilesIn {params.cdna_reads} {params.barcode_reads} \
                    --runThreadN {threads} \
                    --clipAdapterType {params.clipAdapterType} --outFilterScoreMin {params.outFilterScoreMin} --soloCBmatchWLtype {params.soloCBmatchWLtype} \
                    --soloUMIfiltering {params.soloUMIfiltering} --soloUMIdedup {params.soloUMIdedup} \
                    --outSAMtype {params.outSAMtype}\
                    --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} LB:{params.chemistry} \
                    --outSAMattributes {params.outSAMattributes} \
                    --outSAMunmapped Within \
                    --outFileNamePrefix ./{wildcards.sample}/\
                    {params.variousParams} \
                    2>&1 | tee ../../../{log} ;
                    pwd ;\
                    cd ../../../;\
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/features.tsv" "{output.features_file}" ;
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/matrix.mtx" "{output.matrix_file}" ; 
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/barcodes.tsv" "{output.barcodes_file}" ;\
                    mv "{params.starsolo_out}/{wildcards.sample}/Aligned.sortedByCoord.out.bam" "{output.mapped_bam_file}";
                    '''
        rule starsolo_CB_UMI_Simple_unmapped_extracted_sorted:
            input:
                mapped_bam_file = os.path.join(
                    config["output"]["host"],
                    "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam")
            output:
                unmapped_bam_sorted_file = os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
                unmapped_bam_sorted_index = os.path.join(
                        config["output"]["host"],
                        "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bai")
            params:
                unmapped_bam_unsorted_file = os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByCoord_unmapped_out.bam")
            ## because bam is sorted by Coord,it's necessary to sort it by read name
            conda:
                config["envs"]["star"]
            threads:
                config["resources"]["samtools_extract"]["threads"]
            resources:
                mem_mb=config["resources"]["samtools_extract"]["mem_mb"],
            log:
                os.path.join(config["logs"]["host"],
                            "{sample}/unmapped_extracted_sorted_bam.log")
            benchmark:
                os.path.join(config["benchmarks"]["host"],
                            "{sample}/unmapped_extracted_sorted_bam.benchmark")
            shell:
                '''
                samtools view --threads  {threads}  -b -f 4  {input.mapped_bam_file}  >  {params.unmapped_bam_unsorted_file};\
                samtools sort -n  --threads  {threads} {params.unmapped_bam_unsorted_file} -o {output.unmapped_bam_sorted_file};\
                samtools index -@  {threads} {output.unmapped_bam_sorted_file} -o {output.unmapped_bam_sorted_index};\
                rm -rf {params.unmapped_bam_unsorted_file};
                '''
        rule starsolo_CB_UMI_Simple_all:
            input:
                expand(os.path.join(config["output"]["host"],"unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),sample=SAMPLES_ID_LIST),
                expand(os.path.join(
                        config["output"]["host"],
                        "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bai"), sample=SAMPLES_ID_LIST)
        
    else:
        rule starsolo_CB_UMI_Simple_all:
            input: 

    if config["params"]["host"]["starsolo"]["soloType"]=="SmartSeq":
        rule generate_pe_manifest_file:
            input:
                config["params"]["samples"],
            output:
                PE_MANIFEST_FILE = os.path.join("data", "{sample}-pe-manifest.tsv")
            script:
                "../scripts/generate_PE_manifest_file.py"      

        rule generate_se_manifest_file:
            input:
                config["params"]["samples"],
            output:
                PE_MANIFEST_FILE = os.path.join("data", "{sample}-se-manifest.tsv")
            script:
                "../scripts/generate_SE_manifest_file.py"

        rule starsolo_smartseq_PE_count:
            # Input files
            input:
                # Path to the input manifest file
                manifest = os.path.join("data", "{sample}-pe-manifest.tsv"),
            output:
                # # Path to the output features.tsv file
                # features_file = os.path.join(
                #     config["output"]["host"],
                #     "starsolo_count/{sample}/features.tsv"),
                # # Path to the output matrix.mtx file
                # matrix_file = os.path.join(
                #     config["output"]["host"],
                #     "starsolo_count/{sample}/matrix.mtx"),
                # # Path to the output barcodes.tsv file
                # barcodes_file = os.path.join(
                #     config["output"]["host"],
                #     "starsolo_count/{sample}/barcodes.tsv"),
                mapped_pe_bam_file = os.path.join(
                    config["output"]["host"],
                    "starsolo_count_PE/{sample}/Aligned_out.bam"),
                starsolo_count_report = os.path.join(
                    config["output"]["host"],
                    "starsolo_count_PE/{sample}/Log.final.out")
            params:
                cdna_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                # Path to the output directory
                starsolo_out = os.path.join(
                    config["output"]["host"],
                    "starsolo_count_PE/"),
                # Path to the STAR index directory
                reference = config["params"]["host"]["starsolo"]["reference"],
                # Type of sequencing library
                soloType = config["params"]["host"]["starsolo"]["soloType"],
                soloUMIdedup = config["params"]["host"]["starsolo"]["algorithm"]["soloUMIdedup"],
                soloCellFilter = config["params"]["host"]["starsolo"]["algorithm"]["soloCellFilter"],
                soloStrand = config["params"]["host"]["starsolo"]["soloStrand"],
                # SAMattrRGline = microcat.get_SAMattrRGline_from_manifest(config["params"]["host"]["starsolo"]["manifest"]),
                SAMattrRGline = lambda wildcards: microcat.get_SAMattrRGline_by_sample(SAMPLES, wildcards),
                # Additional parameters for STAR
                variousParams = config["params"]["host"]["starsolo"]["variousParams"],
            log:
                os.path.join(config["logs"]["host"],
                            "{sample}/starsolo_count_smartseq2_PE.log")
            benchmark:
                os.path.join(config["benchmarks"]["host"],
                            "{sample}/starsolo_count_smartseq2_PE.benchmark")
            conda:
                config["envs"]["star"]
            resources:
                mem_mb=config["resources"]["starsolo"]["mem_mb"]
            threads:
                # Number of threads for STAR
                config["resources"]["starsolo"]["threads"]
            shell:
                '''
                if [ -s "{input.manifest}" ]; then
                    if echo {params.cdna_reads} | grep -q "\.gz" ; then
                        file_command='--readFilesCommand zcat'
                    else
                        file_command=''
                    fi

                    mkdir -p {params.starsolo_out}; 
                    cd {params.starsolo_out} ;
                    STAR \
                    --soloType SmartSeq \
                    --genomeDir {params.reference} \
                    --readFilesManifest ../../../{input.manifest} \
                    --runThreadN {threads} \
                    --soloUMIdedup {params.soloUMIdedup} \
                    --soloStrand {params.soloStrand} \
                    --soloCellFilter {params.soloCellFilter}\
                    --outSAMtype BAM Unsorted\
                    $file_command \
                    --outSAMunmapped Within \
                    --outSAMattrRGline {params.SAMattrRGline}\
                    --outFileNamePrefix ./{wildcards.sample}/\
                    {params.variousParams} \
                    2>&1 | tee ../../../{log} ;
                    pwd ;\
                    cd ../../../;\
                    mv "{params.starsolo_out}/{wildcards.sample}/Aligned.out.bam" "{output.mapped_pe_bam_file}";\
                else
                    touch "{output.mapped_pe_bam_file}"
                    touch "{output.starsolo_count_report}"
                fi
                '''
        rule starsolo_smartseq_SE_count:
            # Input files
            input:
                # Path to the input manifest file
                manifest = os.path.join("data", "{sample}-se-manifest.tsv"),
            output:
                # # Path to the output features.tsv file
                # features_file = os.path.join(
                #     config["output"]["host"],
                #     "starsolo_count/{sample}/features.tsv"),
                # # Path to the output matrix.mtx file
                # matrix_file = os.path.join(
                #     config["output"]["host"],
                #     "starsolo_count/{sample}/matrix.mtx"),
                # # Path to the output barcodes.tsv file
                # barcodes_file = os.path.join(
                #     config["output"]["host"],
                #     "starsolo_count/{sample}/barcodes.tsv"),
                mapped_se_bam_file = os.path.join(
                    config["output"]["host"],
                    "starsolo_count_SE/{sample}/Aligned_out.bam"),
                starsolo_count_report = os.path.join(
                    config["output"]["host"],
                    "starsolo_count_SE/{sample}/Log.final.out")
            params:
                cdna_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                # Path to the output directory
                starsolo_out = os.path.join(
                    config["output"]["host"],
                    "starsolo_count_SE/"),
                # Path to the STAR index directory
                reference = config["params"]["host"]["starsolo"]["reference"],
                # Type of sequencing library
                soloType = config["params"]["host"]["starsolo"]["soloType"],
                soloUMIdedup = config["params"]["host"]["starsolo"]["algorithm"]["soloUMIdedup"],
                soloCellFilter = config["params"]["host"]["starsolo"]["algorithm"]["soloCellFilter"],
                soloStrand = config["params"]["host"]["starsolo"]["soloStrand"],
                # SAMattrRGline = microcat.get_SAMattrRGline_from_manifest(config["params"]["host"]["starsolo"]["manifest"]),
                SAMattrRGline = lambda wildcards: microcat.get_SAMattrRGline_by_sample(SAMPLES, wildcards),
                # Additional parameters for STAR
                variousParams = config["params"]["host"]["starsolo"]["variousParams"],
            log:
                os.path.join(config["logs"]["host"],
                            "starsolo/{sample}/starsolo_count_smartseq2_SE.log")
            benchmark:
                os.path.join(config["benchmarks"]["host"],
                            "starsolo/{sample}/starsolo_count_smartseq2_SE.benchmark")
            conda:
                config["envs"]["star"]
            resources:
                mem_mb=config["resources"]["starsolo"]["mem_mb"]
            threads:
                # Number of threads for STAR
                config["resources"]["starsolo"]["threads"]
            shell:
                '''
                if [ -s "{input.manifest}" ]; then
                    if echo {params.cdna_reads} | grep -q "\.gz" ; then
                        file_command='--readFilesCommand zcat'
                    else
                        file_command=''
                    fi

                    mkdir -p {params.starsolo_out}
                    cd {params.starsolo_out}
                    STAR \
                    --soloType SmartSeq \
                    --genomeDir {params.reference} \
                    --readFilesManifest ../../../{input.manifest} \
                    --runThreadN {threads} \
                    --soloUMIdedup {params.soloUMIdedup} \
                    --soloStrand {params.soloStrand} \
                    --soloCellFilter {params.soloCellFilter} \
                    --outSAMtype BAM Unsorted \
                    $file_command \
                    --outSAMunmapped Within \
                    --outSAMattrRGline {params.SAMattrRGline} \
                    --outFileNamePrefix ./{wildcards.sample}/ \
                    {params.variousParams} \
                    2>&1 | tee ../../../{log}
                    pwd
                    cd ../../../
                    mv "{params.starsolo_out}/{wildcards.sample}/Aligned.out.bam" "{output.mapped_se_bam_file}"
                else
                    touch "{output.mapped_se_bam_file}"
                    touch "{output.starsolo_count_report}"
                fi
                '''
        rule starsolo_smartseq_combined_bam:
            input:
                mapped_se_bam_file = os.path.join(
                    config["output"]["host"],
                    "starsolo_count_SE/{sample}/Aligned_out.bam"),
                mapped_pe_bam_file = os.path.join(
                    config["output"]["host"],
                    "starsolo_count_PE/{sample}/Aligned_out.bam")
            output:
                unmapped_sorted_bam_file = os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
                unmapped_bam_sorted_index = os.path.join(
                        config["output"]["host"],
                        "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bai"),
            params:
                unmapped_bam_unsorted_file = temp(os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_out_unmapped.bam")),
                unmapped_se_bam_unsorted_file = temp(os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_out_unmapped_SE.bam")),
                unmapped_pe_bam_unsorted_file = temp(os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_out_unmapped_PE.bam"))
            log:
                os.path.join(config["logs"]["host"],
                            "{sample}/unmapped_sorted_bam.log")
            benchmark:
                os.path.join(config["benchmarks"]["host"],
                            "{sample}/unmapped_sorted_bam.benchmark")
            threads:
                config["resources"]["samtools_extract"]["threads"]
            conda:
                config["envs"]["star"]
            resources:
                mem_mb=config["resources"]["samtools_extract"]["mem_mb"],
            shell:
                '''
                if [ -s "{input.mapped_se_bam_file}" ] && [ -s "{input.mapped_pe_bam_file}" ]; then
                    samtools view --threads {threads} -b -f 4 {input.mapped_se_bam_file} > {params.unmapped_se_bam_unsorted_file}
                    samtools view --threads {threads} -b -f 4 {input.mapped_pe_bam_file} > {params.unmapped_pe_bam_unsorted_file}
                    samtools merge -@ {threads} {params.unmapped_bam_unsorted_file} {params.unmapped_pe_bam_unsorted_file} {params.unmapped_se_bam_unsorted_file}
                    rm -rf {params.unmapped_se_bam_unsorted_file}
                    rm -rf {params.unmapped_se_bam_unsorted_file}
                elif [ -s "{input.mapped_pe_bam_file}" ]; then
                    samtools view --threads {threads} -b -f 4 {input.mapped_pe_bam_file} > {params.unmapped_bam_unsorted_file}
                elif [ -s "{input.mapped_se_bam_file}" ]; then
                    samtools view --threads {threads} -b -f 4 {input.mapped_se_bam_file} > {params.unmapped_bam_unsorted_file}
                else
                    echo "Both mapped BAM files are empty! Exiting..."
                    exit 1
                fi

                samtools sort -n  --threads  {threads} {params.unmapped_bam_unsorted_file} -o {output.unmapped_sorted_bam_file};\
                samtools index -@  {threads} {output.unmapped_sorted_bam_file} -o {output.unmapped_bam_sorted_index};\
                
                rm -rf {params.unmapped_bam_unsorted_file}
                '''
        # rule starsolo_smartseq_extracted_sorted:
        #     input:
        #         mapped_bam_file = os.path.join(
        #             config["output"]["host"],
        #             "starsolo_count/{sample}/Aligned_out.bam")
        #     output:
        #         unmapped_sorted_bam_file = os.path.join(
        #             config["output"]["host"],
        #             "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
        #         unmapped_bam_sorted_index = os.path.join(
        #                 config["output"]["host"],
        #                 "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bai")
        #     params:
        #         unmapped_bam_unsorted_file = temp(os.path.join(
        #             config["output"]["host"],
        #             "unmapped_host/{sample}/Aligned_out_unmapped.bam"))
        #     log:
        #         os.path.join(config["logs"]["host"],
        #                     "starsolo/{sample}/unmapped_sorted_bam.log")
        #     benchmark:
        #         os.path.join(config["benchmarks"]["host"],
        #                     "starsolo/{sample}/unmapped_sorted_bam.benchmark")
        #     threads:
        #         config["resources"]["samtools_extract"]["threads"]
        #     conda:
        #         config["envs"]["star"]
        #     resources:
        #         mem_mb=config["resources"]["samtools_extract"]["mem_mb"],
        #     shell:
        #         '''
        #         samtools view --threads  {threads}  -b -f 4   {input.mapped_bam_file}  >  {params.unmapped_bam_unsorted_file};\
        #         samtools sort -n  --threads  {threads} {params.unmapped_bam_unsorted_file} -o {output.unmapped_sorted_bam_file};\
        #         samtools index -@  {threads} {output.unmapped_sorted_bam_file} -o {output.unmapped_bam_sorted_index};\
        #         rm -rf {params.unmapped_bam_unsorted_file}
        #         '''

        # rule starsolo_smartseq_unmapped_sorted_bam:
        #     input:
        #         unmapped_bam_unsorted_file = os.path.join(
        #             config["output"]["host"],
        #             "unmapped_host/{sample}/Aligned_out_unmapped.bam")
        #     output:
        #         unmapped_sorted_bam_file = os.path.join(
        #             config["output"]["host"],
        #             "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
        #         unmapped_bam_sorted_index = os.path.join(
        #                 config["output"]["host"],
        #                 "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bai")
        #     params:
        #         tag="RG"
        #     log:
        #         os.path.join(config["logs"]["host"],
        #                     "starsolo/{sample}/unmapped_sorted_bam.log")
        #     benchmark:
        #         os.path.join(config["benchmarks"]["host"],
        #                     "starsolo/{sample}/unmapped_sorted_bam.benchmark")
        #     threads:
        #         config["resources"]["samtools_extract"]["threads"]
        #     resources:
        #         mem_mb=config["resources"]["samtools_extract"]["mem_mb"],
        #     conda:
        #         config["envs"]["star"]
        #     shell:
        #         '''
        #         samtools sort -n  --threads  {threads} {input.unmapped_bam_unsorted_file} -o {output.unmapped_sorted_bam_file};\
        #         samtools index -@  {threads} {output.unmapped_sorted_bam_file} -o {output.unmapped_bam_sorted_index};\
        #         '''
        # checkpoint starsolo_smartseq_demultiplex_bam_by_read_group:
        #     input:
        #         unmapped_sorted_bam = os.path.join(
        #             config["output"]["host"],
        #             "starsolo_count/Aligned_out_unmapped_RGsorted.bam")
        #     output:
        #         unmapped_bam_demultiplex_dir = directory(os.path.join(
        #             config["output"]["host"],
        #             "unmapped_host/"))
        #     params:
        #         threads = 40, # Number of threads
        #         tag="RG"
        #     conda:
        #         config["envs"]["star"]
        #     log:
        #         os.path.join(
        #             config["logs"]["host"],
        #             "starsolo_count/demultiplex_bam_by_read_group.log")
        #     benchmark:
        #         os.path.join(
        #             config["benchmarks"]["host"], 
        #             "starsolo_count/demultiplex_bam_by_read_group.benchmark")
        #     shell:
        #         """
        #         python /data/project/host-microbiome/microcat/microcat/scripts/spilt_bam_by_tag.py --tag {params.tag} --bam_path {input.unmapped_sorted_bam} --output_dir {output.unmapped_bam_demultiplex_dir}  --log_file {log}
        #         """
        # split the PathSeq BAM into one BAM per cell barcode
        # rule split_starsolo_BAM_by_RG:
        #     input:
        #         unmapped_sorted_bam = os.path.join(
        #                 config["output"]["host"],
        #                 "starsolo_count/Aligned_out_unmapped_RGsorted.bam")
        #     output:
        #         unmapped_bam_sorted_file =os.path.join(
        #             config["output"]["host"],
        #             "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
        #     params:
        #         SampleID="{sample}",
        #     shell:
        #         '''
        #         python /data/project/host-microbiome/microcat/microcat/scripts/split_Starsolo_BAM_by_RG.py \
        #         --bam_path {input.unmapped_sorted_bam} \
        #         --tag {params.SampleID} \
        #         --output_bam {output.unmapped_bam_sorted_file} 
        #         '''

        rule starsolo_SmartSeq_all:
            input:
                expand(os.path.join(
                        config["output"]["host"],
                        "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"), sample=SAMPLES_ID_LIST),
                expand(os.path.join(
                        config["output"]["host"],
                        "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bai"), sample=SAMPLES_ID_LIST)
        
    else:
        rule starsolo_SmartSeq_all:
            input: 


    if config["params"]["host"]["starsolo"]["soloType"]=="CB_UMI_Complex":
            rule starsolo_CB_UMI_Complex_count:
                input:
                    # Directory containing input fastq files
                    fastqs_dir=lambda wildcards: microcat.get_fastqs_dir(SAMPLES,wildcards),
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
                        "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam"),
                    starsolo_count_report = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/Log.final.out")
                params:
                    barcode_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                    cdna_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq2"),
                    starsolo_out = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/"),
                    reference = config["params"]["host"]["starsolo"]["reference"],
                    chemistry = config["params"]["host"]["starsolo"]["name"],
                    soloCBlen = config["params"]["host"]["starsolo"]["barcode"]["soloCBlen"],
                    soloCBposition = config["params"]["host"]["starsolo"]["barcode"]["soloCBposition"],
                    soloBarcodeMate = config["params"]["host"]["starsolo"]["barcode"]["soloBarcodeMate"],
                    soloUMIposition = config["params"]["host"]["starsolo"]["barcode"]["soloUMIposition"],
                    clip5pNbases = config["params"]["host"]["starsolo"]["barcode"]["clip5pNbases"],
                    clip3pNbases = config["params"]["host"]["starsolo"]["barcode"]["clip3pNbases"],
                    soloUMIdedup = config["params"]["host"]["starsolo"]["algorithm"]["soloUMIdedup"],
                    soloCBmatchWLtype = config["params"]["host"]["starsolo"]["algorithm"]["soloCBmatchWLtype"],
                    soloUMIfiltering = config["params"]["host"]["starsolo"]["algorithm"]["soloUMIfiltering"],
                    soloAdapterSequence = config["params"]["host"]["starsolo"]["soloAdapterSequence"],
                    soloCellFilter = config["params"]["host"]["starsolo"]["algorithm"]["soloCellFilter"],
                    clipAdapterType = config["params"]["host"]["starsolo"]["algorithm"]["clipAdapterType"],
                    outFilterScoreMin = config["params"]["host"]["starsolo"]["algorithm"]["outFilterScoreMin"],
                    soloMultiMappers = config["params"]["host"]["starsolo"]["algorithm"]["soloMultiMappers"],
                    soloAdapterMismatchesNmax = config["params"]["host"]["starsolo"]["algorithm"]["soloAdapterMismatchesNmax"],
                    soloStrand = config["params"]["host"]["starsolo"]["soloStrand"],
                    barcode_list_1 =  os.path.join(config["datas"]["barcode_list_dirs"]["tenX"],
                                                config["params"]["host"]["starsolo"]["barcode"]["soloCBwhitelist"]["barcode_1"]),
                    barcode_list_2 =  os.path.join(config["datas"]["barcode_list_dirs"]["tenX"],
                                                config["params"]["host"]["starsolo"]["barcode"]["soloCBwhitelist"]["barcode_2"]),
                    outSAMattributes = config["params"]["host"]["starsolo"]["outSAMattributes"],
                    outSAMtype = config["params"]["host"]["starsolo"]["outSAMtype"],
                    variousParams = config["params"]["host"]["starsolo"]["variousParams"],
                    description = config["params"]["host"]["starsolo"]["description"],
                resources:
                    mem_mb=config["resources"]["starsolo"]["mem_mb"]
                threads:
                    config["resources"]["starsolo"]["threads"]
                log:
                    os.path.join(config["logs"]["host"],
                                "{sample}/{sample}_starsolo_count.log")
                benchmark:
                    os.path.join(config["benchmarks"]["host"],
                                "{sample}/{sample}_starsolo_count.benchmark")
                conda:
                    config["envs"]["star"]
                message: "Executing starsolo with {threads} threads on the following files {wildcards.sample}.Library with {params.description}"
                shell:
                    '''
                    if echo {params.cdna_reads} | grep -q "\.gz" ; then
                        file_command='--readFilesCommand zcat'
                    else
                        file_command=''
                    fi

                    mkdir -p {params.starsolo_out}; 
                    cd {params.starsolo_out} ;
                    STAR \
                    $file_command  \
                    --soloType CB_UMI_Complex \
                    --soloCBwhitelist {params.barcode_list_1} {params.barcode_list_2}  \
                    --soloCBposition {params.soloCBposition}  \
                    --soloUMIposition {params.soloUMIposition} \
                    --genomeDir {params.reference} \
                    --readFilesIn {params.cdna_reads} {params.barcode_reads} \
                    --runThreadN {threads} \
                    --soloAdapterMismatchesNmax {params.soloAdapterMismatchesNmax} \
                    --soloAdapterSequence {params.soloAdapterSequence} \
                    --soloCBmatchWLtype {params.soloCBmatchWLtype} \
                    --soloCellFilter {params.soloCellFilter} \
                    --soloStrand {params.soloStrand} \
                    --outSAMtype {params.outSAMtype}\
                    --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} LB:{params.chemistry} \
                    --outSAMattributes {params.outSAMattributes} \
                    --outSAMunmapped Within \
                    --outFileNamePrefix ./{wildcards.sample}/\
                    {params.variousParams}  \
                    2>&1 | tee ../../../{log} ;
                    pwd ;\
                    cd ../../../;\
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/raw/features.tsv" "{output.features_file}" ;\
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/raw/matrix.mtx" "{output.matrix_file}" ; \
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/raw/barcodes.tsv" "{output.barcodes_file}" ;\
                    mv "{params.starsolo_out}/{wildcards.sample}/Aligned.sortedByCoord.out.bam" "{output.mapped_bam_file}";
                    '''
            rule starsolo_CB_UMI_Complex_unmapped_extracted_sorted:
                input:
                    mapped_bam_file = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam")
                output:
                    unmapped_bam_sorted_file = os.path.join(
                        config["output"]["host"],
                        "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
                    unmapped_bam_sorted_index = os.path.join(
                            config["output"]["host"],
                            "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bai")
                params:
                    unmapped_bam_unsorted_file = os.path.join(
                        config["output"]["host"],
                        "unmapped_host/{sample}/Aligned_sortedByCoord_unmapped_out.bam")
                ## because bam is sorted by Coord,it's necessary to sort it by read name
                conda:
                    config["envs"]["star"]
                threads:
                    config["resources"]["samtools_extract"]["threads"]
                resources:
                    mem_mb=config["resources"]["samtools_extract"]["mem_mb"],
                log:
                    os.path.join(config["logs"]["host"],
                                "{sample}/unmapped_extracted_sorted_bam.log")
                benchmark:
                    os.path.join(config["benchmarks"]["host"],
                                "{sample}/unmapped_extracted_sorted_bam.benchmark")
                shell:
                    '''
                    samtools view --threads  {threads}  -b -f 4  {input.mapped_bam_file}  >  {params.unmapped_bam_unsorted_file};\
                    samtools sort -n  --threads  {threads} {params.unmapped_bam_unsorted_file} -o {output.unmapped_bam_sorted_file};\
                    samtools index -@  {threads} {output.unmapped_bam_sorted_file} -o {output.unmapped_bam_sorted_index};\
                    rm -rf {params.unmapped_bam_unsorted_file};
                    '''
            rule starsolo_CB_UMI_Complex_all:
                input:
                    expand(os.path.join(config["output"]["host"],"unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),sample=SAMPLES_ID_LIST),
                    expand(os.path.join(
                            config["output"]["host"],
                            "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bai"), sample=SAMPLES_ID_LIST)
    else:
        rule starsolo_CB_UMI_Complex_all:
            input: 
    rule starsolo_all:
        input: 
            rules.starsolo_SmartSeq_all.input,
            rules.starsolo_CB_UMI_Simple_all.input,
            rules.starsolo_CB_UMI_Complex_all.input

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
            fastqs_dir=lambda wildcards: microcat.get_fastqs_dir(SAMPLES,wildcards),
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
                "cellranger_count/{sample}/{sample}_mappped2human_bam.bam.bai"),
            metrics_summary = os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/{sample}.metrics_summary.csv"),
            web_summary = os.path.join(
                config["output"]["host"],
                "cellranger_count/{sample}/{sample}.web_summary.html"),
        priority: 10
        params:
            cr_out = os.path.join(
                config["output"]["host"],
                "cellranger_count/"),
            reference = config["params"]["host"]["cellranger"]["reference"],
            # local_cores = config["params"]["host"]["cellranger"]["local_cores"],
            SampleID="{sample}",
            variousParams = config["params"]["host"]["cellranger"]["variousParams"],
        # resources:
        #     mem_mb=config["tools"]["cellranger_count"]["mem_mb"],
        #     runtime=config["tools"]["cellranger_count"]["runtime"],
        threads: 
            config["resources"]["cellranger"]["threads"]
        resources:
            mem_mb=config["resources"]["cellranger"]["mem_mb"]
        conda:
            config["envs"]["star"]
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
            ln -sr "{params.cr_out}{params.SampleID}/outs/web_summary.html" "{output.web_summary}" ; 
            ln -sr "{params.cr_out}{params.SampleID}/outs/metrics_summary.csv" "{output.metrics_summary}";
            ln -sr "{params.cr_out}{params.SampleID}/outs/possorted_genome_bam.bam" "{output.mapped_bam_file}";
            ln -sr "{params.cr_out}{params.SampleID}/outs/possorted_genome_bam.bam.bai" "{output.mapped_bam_index_file}";
            '''
    rule cellranger_unmapped_extracted_sorted:
        input:
            # unmapped_bam_unsorted_file = os.path.join(
            # config["output"]["host"],
            # "cellranger_count/{sample}/{sample}_unmappped2human_sorted_bam.bam")
            mapped_bam_file = os.path.join(
            config["output"]["host"],
            "cellranger_count/{sample}/{sample}_mappped2human_bam.bam")
        output:
            unmapped_bam_sorted_file = os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
            unmapped_bam_sorted_index = os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bai")
        params:
            unmapped_bam_unsorted_file = os.path.join(
                config["output"]["host"],
                "unmapped_host/{sample}/Aligned_sortedByCoord_unmapped_out.bam")
        ## because bam is sorted by Coord,it's necessary to sort it by read name
        resources:
            mem_mb=config["resources"]["samtools_extract"]["mem_mb"]
        threads:
            config["resources"]["samtools_extract"]["threads"]
        conda:
            config["envs"]["star"]
        shell:
            '''
            samtools view --threads  {threads}  -b -f 4   {input.mapped_bam_file}  >  {params.unmapped_bam_unsorted_file};\
            samtools sort -n  --threads  {threads} {params.unmapped_bam_unsorted_file} -o {output.unmapped_bam_sorted_file};\
            samtools index -@  {threads} {output.unmapped_bam_sorted_file} -o {output.unmapped_bam_sorted_index};\
            rm -rf {params.unmapped_bam_unsorted_file} 
            '''

    rule cellranger_all:
        input:
            expand(os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"), sample=SAMPLES_ID_LIST),
            expand(os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bai"), sample=SAMPLES_ID_LIST)

else:
    rule cellranger_all:
        input:

rule host_all:
    input:
        rules.starsolo_all.input,
        rules.cellranger_all.input,