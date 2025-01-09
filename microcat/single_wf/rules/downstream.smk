if config["params"]["begin"] == "host":
    if config["params"]["host"]["starsolo"]["do"]:
        if config["params"]["host"]["cellbender"]["do"]:
            if config["params"]["host"]["cellbender"]["gpu"]:
                rule cellbender_filter:
                    input:
                        unmapped_bam_sorted_file =os.path.join(
                            config["output"]["host"],
                            "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
                    output:
                        filtered_hdf5 = os.path.join(
                                config["output"]["profile"],
                                "{sample}/cellbender/filterd_feature_bc_matrix.h5") 
                    params:
                        fpr_cutoff = config["params"]["host"]["cellbender"]["fpr"],
                        raw_mtx_dir = os.path.join(
                                config["output"]["host"],
                                "starsolo_count/{sample}/Solo.out/Gene/raw"),
                        variousParams = config["params"]["host"]["cellbender"]["variousParams"],
                    log:
                        os.path.join(config["logs"]["profile"],
                                    "cellbender/{sample}_create_hdf5.log")
                    threads: 
                        config["resources"]["cellbender"]["threads"]
                    benchmark:
                        os.path.join(config["benchmarks"]["profile"],
                                    "cellbender/{sample}_cellbender.benchmark")
                    shell:
                        '''
                        cellbender remove-background \
                        --cuda \
                        --input {params.raw_mtx_dir} \
                        --output {output.filtered_hdf5} \
                        --fpr {params.fpr_cutoff}
                        '''
            else:
                rule cellbender_filter:
                    input:
                        unmapped_bam_sorted_file =os.path.join(
                            config["output"]["host"],
                            "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
                    output:
                        filtered_hdf5 = os.path.join(
                                config["output"]["profile"],
                                "{sample}/cellbender/filterd_feature_bc_matrix.h5") 
                    params:
                        fpr_cutoff = config["params"]["host"]["cellbender"]["fpr"],
                        raw_mtx_dir = os.path.join(
                                config["output"]["host"],
                                "starsolo_count/{sample}/Solo.out/Gene/raw"),
                        variousParams = config["params"]["host"]["cellbender"]["variousParams"],
                    log:
                        os.path.join(config["logs"]["profile"],
                                    "cellbender/{sample}_create_hdf5.log")
                    threads: 
                        config["resources"]["cellbender"]["threads"]
                    benchmark:
                        os.path.join(config["benchmarks"]["profile"],
                                    "cellbender/{sample}_cellbender.benchmark")
                    shell:
                        '''
                        cellbender remove-background \
                        --cpu-threads {threads} \
                        --input {params.raw_mtx_dir} \
                        --output {output.filtered_hdf5} \
                        --fpr {params.fpr_cutoff}
                        '''

            rule leiden_pre_cluster:
                input:
                    filtered_hdf5 = os.path.join(
                            config["output"]["profile"],
                            "{sample}/cellbender/filterd_feature_bc_matrix.h5") 
                output:
                    ledian_cluster = os.path.join(
                            config["output"]["profile"],
                            "{sample}/cellbender/leiden_cluster.tsv") 
                log:
                    os.path.join(config["logs"]["profile"],
                                "{sample}/{sample}_leidan_cellbender_clsuter.log")
                params:
                    ledian_cluster_noncellbender_script = config["scripts"]["leiden_pre_cluster"],
                resources:
                    mem_mb=config["resources"]["samtools_extract"]["mem_mb"],    
                shell:
                    '''
                    python {params.ledian_cluster_noncellbender_script} \
                    --input_hdf5 {input.filtered_hdf5} \
                    --output_cluster {output.ledian_cluster}\
                    '''    
        elif config["params"]["host"]["starsolo"]["soloType"]=="CB_UMI_Complex" or config["params"]["host"]["starsolo"]["soloType"]=="CB_UMI_Simple":
            rule leiden_pre_cluster:
                input:
                    unmapped_bam_sorted_file =os.path.join(
                        config["output"]["host"],
                        "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
                output:
                    ledian_cluster = os.path.join(
                            config["output"]["profile"],
                            "{sample}/cellbender/leiden_cluster.tsv") 
                log:
                    os.path.join(config["logs"]["profile"],
                                "{sample}/{sample}_leidan_cellbender_clsuter.log")
                params:
                    filter_mtx_dir = os.path.join(
                                        config["output"]["host"],
                                        "starsolo_count/{sample}/Solo.out/Gene/filtered"),
                    ledian_cluster_noncellbender_script = config["scripts"]["ledian_cluster_noncellbender"],
                resources:
                    mem_mb=config["resources"]["samtools_extract"]["mem_mb"],
                shell:
                    '''
                    python {params.ledian_cluster_noncellbender_script} \
                    -i {params.filter_mtx_dir} \
                    --output_cluster {output.ledian_cluster}\
                    &> {log}

                    '''
        else:
            rule leiden_pre_cluster:
                input:
                    unmapped_bam_sorted_file =os.path.join(
                        config["output"]["host"],
                        "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
                output:
                    ledian_cluster = os.path.join(
                            config["output"]["profile"],
                            "{sample}/cellbender/leiden_cluster.tsv")
                params:
                    SE_filter_mtx_dir = os.path.join(
                                        config["output"]["host"],
                                        "starsolo_count_SE/{sample}/Solo.out/Gene/raw"),
                    PE_filter_mtx_dir = os.path.join(
                                        config["output"]["host"],
                                        "starsolo_count_PE/{sample}/Solo.out/Gene/raw"),             
                    ledian_cluster_noncellbender_script = config["scripts"]["ledian_cluster_noncellbender"]
                log:
                    os.path.join(config["logs"]["profile"],
                                "{sample}/{sample}_leidan_cellbender_clsuter.log")
                resources:
                    mem_mb=config["resources"]["samtools_extract"]["mem_mb"],
                shell:
                    '''
                    if [ -d "{params.SE_filter_mtx_dir}" ] && [ "$(ls -A {params.SE_filter_mtx_dir})" ]; then
                        python {params.ledian_cluster_noncellbender_script} \
                        -i {params.SE_filter_mtx_dir} \
                        --output_cluster {output.ledian_cluster}\
                        &> {log}
                    elif [ -d "{params.PE_filter_mtx_dir}" ] && [ "$(ls -A {params.PE_filter_mtx_dir})" ]; then
                        python {params.ledian_cluster_noncellbender_script} \
                        -i {params.PE_filter_mtx_dir} \
                        --output_cluster {output.ledian_cluster}\
                        &> {log}
                    else
                        echo "No mtx files found"
                        exit 1
                    fi
                    '''

        rule starsolo_downstream_all:
            input:
                expand(os.path.join(
                        config["output"]["profile"],
                        "{sample}/cellbender/leiden_cluster.tsv"), sample=SAMPLES_ID_LIST),
            
    else:
        rule starsolo_downstream_all:
            input: 

    if config["params"]["host"]["cellranger"]["do"]:
        if config["params"]["host"]["cellbender"]["do"]:
            if config["params"]["host"]["cellbender"]["gpu"]:
                rule cellbender_filter:
                    input:
                        unmapped_bam_sorted_file =os.path.join(
                                config["output"]["host"],
                                "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
                    output:
                        filtered_hdf5 = os.path.join(
                                config["output"]["profile"],
                                "{sample}/cellbender/filterd_feature_bc_matrix.h5") 
                    params:
                        fpr_cutoff = config["params"]["host"]["cellbender"]["fpr"],
                        variousParams = config["params"]["host"]["cellbender"]["variousParams"],
                        raw_hdf5 = os.path.join(
                            config["output"]["host"],
                            "cellranger_count/{sample}/outs/raw_feature_bc_matrix.h5")
                    log:
                        os.path.join(config["logs"]["profile"],
                                    "cellbender/{sample}_create_hdf5.log")
                    threads: 
                        config["resources"]["cellbender"]["threads"]
                    benchmark:
                        os.path.join(config["benchmarks"]["profile"],
                                    "cellbender/{sample}_cellbender.benchmark")
                    shell:
                        '''
                        cellbender remove-background \
                        --cuda \
                        --input {params.raw_hdf5} \
                        --output {output.filtered_hdf5} \
                        --fpr {params.fpr_cutoff}
                        '''
            else:
                rule cellbender_cellranger_filter:
                    input:
                        unmapped_bam_sorted_file =os.path.join(
                                config["output"]["host"],
                                "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
                    output:
                        filtered_hdf5 = os.path.join(
                                config["output"]["profile"],
                                "{sample}/cellbender/filterd_feature_bc_matrix.h5") 
                    params:
                        fpr_cutoff = config["params"]["host"]["cellbender"]["fpr"],
                        variousParams = config["params"]["host"]["cellbender"]["variousParams"],
                        raw_hdf5 = os.path.join(
                            config["output"]["host"],
                            "cellranger_count/{sample}/outs/raw_feature_bc_matrix.h5")
                    threads: 
                        config["resources"]["cellbender"]["threads"]
                    log:
                        os.path.join(config["logs"]["profile"],
                                    "cellbender/{sample}_create_hdf5.log")
                    benchmark:
                        os.path.join(config["benchmarks"]["profile"],
                                    "cellbender/{sample}_cellbender.benchmark")
                    shell:
                        '''
                        cellbender remove-background \
                        --cpu-threads {threads} \
                        --input {params.raw_hdf5} \
                        --output {output.filtered_hdf5} \
                        --fpr {params.fpr_cutoff}
                        '''
            rule leiden_pre_cluster:
                input:
                    filtered_hdf5 = os.path.join(
                            config["output"]["profile"],
                            "{sample}/cellbender/filterd_feature_bc_matrix.h5") 
                output:
                    ledian_cluster = os.path.join(
                            config["output"]["profile"],
                            "{sample}/cellbender/leiden_cluster.tsv") 
                log:
                    os.path.join(config["logs"]["profile"],
                                "{sample}/{sample}_leidan_cellbender_clsuter.log")
                params:
                    leiden_pre_cluster_script = config["scripts"]["leiden_pre_cluster"], 
                resources:
                    mem_mb=config["resources"]["samtools_extract"]["mem_mb"],
                shell:
                    '''
                    python {params.leiden_pre_cluster_script} \
                    --input_hdf5 {input.filtered_hdf5} \
                    --output_cluster {output.ledian_cluster}\
                    '''
        else:
            rule leiden_pre_cluster:
                input:
                    # genes_file = os.path.join(
                    #                     config["output"]["host"],
                    #                     "starsolo_count/{sample}/{sample}_features.tsv"),
                    # matrix_file = os.path.join(
                    #                     config["output"]["host"],
                    #                     "starsolo_count/{sample}/{sample}_matrix.mtx"),
                    # barcodes_file = os.path.join(
                    #                     config["output"]["host"],
                    #                     "starsolo_count/{sample}/{sample}_barcodes.tsv"),
                    unmapped_bam_sorted_file =os.path.join(
                        config["output"]["host"],
                        "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
                output:
                    ledian_cluster = os.path.join(
                            config["output"]["profile"],
                            "{sample}/cellbender/leiden_cluster.tsv") 
                log:
                    os.path.join(config["logs"]["profile"],
                                "{sample}/{sample}_leidan_cellbender_clsuter.log")
                params:
                    filter_mtx_dir = os.path.join(
                                        config["output"]["host"],
                                        "cellranger_count/{sample}/outs/filtered_feature_bc_matrix"),
                    ledian_cluster_noncellbender_script = config["scripts"]["ledian_cluster_noncellbender"],
                resources:
                    mem_mb=config["resources"]["samtools_extract"]["mem_mb"],
                shell:
                    '''
                    python {params.ledian_cluster_noncellbender_script} \
                    -i {params.filter_mtx_dir} \
                    --output_cluster {output.ledian_cluster}\
                    &> {log}
                    '''
        rule cellranger_downstream_all:
            input:
                expand(os.path.join(
                        config["output"]["profile"],
                        "{sample}/cellbender/leiden_cluster.tsv"), sample=SAMPLES_ID_LIST),

    else:
        rule cellranger_downstream_all:
            input: 

    rule downstream_all:
        input:
            rules.starsolo_downstream_all.input,
            rules.cellranger_downstream_all.input

if config["params"]["begin"] == "classifier":        
    if config["params"]["host"]["cellbender"]["do"]:
        if config["params"]["host"]["cellbender"]["gpu"]:
            rule cellbender_filter:
                input:
                    unmapped_bam_sorted_file =os.path.join(
                            config["output"]["host"],
                            "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
                    mtx_file = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "mtx"),
                output:
                    filtered_hdf5 = os.path.join(
                            config["output"]["profile"],
                            "{sample}/cellbender/filterd_feature_bc_matrix.h5") 
                params:
                    fpr_cutoff = config["params"]["host"]["cellbender"]["fpr"],
                    variousParams = config["params"]["host"]["cellbender"]["variousParams"],
                    raw_hdf5 = os.path.join(
                        config["output"]["host"],
                        "cellranger_count/{sample}/outs/raw_feature_bc_matrix.h5")
                log:
                    os.path.join(config["logs"]["profile"],
                                "cellbender/{sample}_create_hdf5.log")
                threads: 
                    config["resources"]["cellbender"]["threads"]
                benchmark:
                    os.path.join(config["benchmarks"]["profile"],
                                "cellbender/{sample}_cellbender.benchmark")
                shell:
                    '''
                    cellbender remove-background \
                    --cuda \
                    --input {input.mtx_file} \
                    --output {output.filtered_hdf5} \
                    --fpr {params.fpr_cutoff}
                    '''
        else:
            rule cellbender_cellranger_filter:
                input:
                    unmapped_bam_sorted_file =os.path.join(
                            config["output"]["host"],
                            "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
                    mtx_file = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "mtx"),
                output:
                    filtered_hdf5 = os.path.join(
                            config["output"]["profile"],
                            "{sample}/cellbender/filterd_feature_bc_matrix.h5") 
                params:
                    fpr_cutoff = config["params"]["host"]["cellbender"]["fpr"],
                    variousParams = config["params"]["host"]["cellbender"]["variousParams"],
                    raw_hdf5 = os.path.join(
                        config["output"]["host"],
                        "cellranger_count/{sample}/outs/raw_feature_bc_matrix.h5")
                threads: 
                    config["resources"]["cellbender"]["threads"]
                log:
                    os.path.join(config["logs"]["profile"],
                                "cellbender/{sample}_create_hdf5.log")
                benchmark:
                    os.path.join(config["benchmarks"]["profile"],
                                "cellbender/{sample}_cellbender.benchmark")
                shell:
                    '''
                    cellbender remove-background \
                    --cpu-threads {threads} \
                    --input {input.mtx_file} \
                    --output {output.filtered_hdf5} \
                    --fpr {params.fpr_cutoff}
                    '''
        rule leiden_pre_cluster:
            input:
                filtered_hdf5 = os.path.join(
                        config["output"]["profile"],
                        "{sample}/cellbender/filterd_feature_bc_matrix.h5") 
            output:
                ledian_cluster = os.path.join(
                        config["output"]["profile"],
                        "{sample}/cellbender/leiden_cluster.tsv") 
            log:
                os.path.join(config["logs"]["profile"],
                            "{sample}/{sample}_leidan_cellbender_clsuter.log")
            params:
                leiden_pre_cluster_script = config["scripts"]["leiden_pre_cluster"], 
            resources:
                mem_mb=config["resources"]["samtools_extract"]["mem_mb"],
            shell:
                '''
                python {params.leiden_pre_cluster_script} \
                --input_hdf5 {input.filtered_hdf5} \
                --output_cluster {output.ledian_cluster}\
                '''
    else:
        rule leiden_pre_cluster:
            input:
                unmapped_bam_sorted_file =os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
                mtx_file = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "mtx"),
            output:
                ledian_cluster = os.path.join(
                        config["output"]["profile"],
                        "{sample}/cellbender/leiden_cluster.tsv") 
            log:
                os.path.join(config["logs"]["profile"],
                            "{sample}/{sample}_leidan_cellbender_clsuter.log")
            params:
                ledian_cluster_noncellbender_script = config["scripts"]["ledian_cluster_noncellbender"],
            resources:
                mem_mb=config["resources"]["samtools_extract"]["mem_mb"],
            shell:
                '''
                python {params.ledian_cluster_noncellbender_script} \
                -i {input.mtx_file} \
                --output_cluster {output.ledian_cluster}\
                &> {log}
                '''
    rule downstream_all:
        input:
            expand(os.path.join(
                                config["output"]["profile"],
                                "{sample}/cellbender/leiden_cluster.tsv"), sample=SAMPLES_ID_LIST)