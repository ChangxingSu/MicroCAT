if config["params"]["host"]["starsolo"]["do"]:
    if config["params"]["host"]["cellbender"]["do"]:
        if config["params"]["host"]["cellbender"]["gpu"] == "off" :
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
                    --cuda \
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
            shell:
                '''
                python {params.ledian_cluster_noncellbender_script} \
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
                                    "starsolo_count/{sample}/Solo.out/Gene/filtered"),
                ledian_cluster_noncellbender_script = config["scripts"]["ledian_cluster_noncellbender"],
            # conda:
            #     "../envs/create_hdf5.yaml"
            shell:
                '''
                python {params.ledian_cluster_noncellbender_script} \
                -i {params.filter_mtx_dir} \
                --output_cluster {output.ledian_cluster}\
                &> {log}
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
        if config["params"]["host"]["cellbender"]["gpu"] == "off" :
            rule cellbender_cellranger_filter:
                input:
                    unmapped_bam_sorted_file =os.path.join(
                            config["output"]["host"],
                            "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
                    raw_hdf5 = os.path.join(
                        config["output"]["host"],
                        "cellranger_count/{sample}/outs/raw_feature_bc_matrix.h5")
                output:
                    filtered_hdf5 = os.path.join(
                            config["output"]["profile"],
                            "{sample}/cellbender/filterd_feature_bc_matrix.h5") 
                params:
                    fpr_cutoff = config["params"]["host"]["cellbender"]["fpr"],
                    variousParams = config["params"]["host"]["cellbender"]["variousParams"],
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
                    --input {input.raw_hdf5} \
                    --output {output.filtered_hdf5} \
                    --fpr {params.fpr_cutoff}
                    '''
        else:
            rule cellbender_filter:
                input:
                    unmapped_bam_sorted_file =os.path.join(
                            config["output"]["host"],
                            "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
                    raw_hdf5 = os.path.join(
                        config["output"]["host"],
                        "cellranger_count/{sample}/outs/raw_feature_bc_matrix.h5")
                output:
                    filtered_hdf5 = os.path.join(
                            config["output"]["profile"],
                            "{sample}/cellbender/filterd_feature_bc_matrix.h5") 
                params:
                    fpr_cutoff = config["params"]["host"]["cellbender"]["fpr"],
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
            # conda:
            #     "../envs/create_hdf5.yaml"
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

# cellbender remove-background \
# --cuda \
# --input /data/comics-sucx/raw_feature_bc_matrix.h5 \
# --output /data/comics-sucx/microcat_debug/microcat_singlecell/cellbender_feature_bc_matrix.h5 \
# --fpr 0.01

# cellbender remove-background \
# --cuda \
# --input /data/comics-sucx/Fn_Ec_16S_S4.h5 \
# --output /data/comics-sucx/microcat_debug/microcat_singlecell/fnec/cellbender_feature_bc_matrix.h5 \
# --fpr 0.01


# /data/comics-sucx/microcat_debug/microcat_singlecell/fnec/fnec_groups.tsv


# python /data/project/host-microbiome/MicroCAT/microcat/single_wf/scripts/krak_sample_denosing.py                 --krak_report results/03.classifier/rmhost_kraken2_report/custom/4_HT29Cells_BMixB_Fn_Ec_16S_S4/4_HT29Cells_BMixB_Fn_Ec_16S_S4_kraken2_report.txt                 --krak_output results/03.classifier/rmhost_kraken2_output/4_HT29Cells_BMixB_Fn_Ec_16S_S4/4_HT29Cells_BMixB_Fn_Ec_16S_S4_kraken2_output.txt                 --krak_mpa_report results/03.classifier/rmhost_kraken2_report/mpa/4_HT29Cells_BMixB_Fn_Ec_16S_S4/4_HT29Cells_BMixB_Fn_Ec_16S_S4_kraken2_mpa_report.txt                 --bam results/03.classifier/rmhost_extracted_classified_output/4_HT29Cells_BMixB_Fn_Ec_16S_S4/4_HT29Cells_BMixB_Fn_Ec_16S_S4_kraken2_extracted_classified.bam                 --ktaxonomy /data/database/kraken2uniq_database/k2_pluspf_16gb_20231009/ktaxonomy.tsv                --inspect /data/database/kraken2uniq_database/k2_pluspf_16gb_20231009/inspect.txt                 --min_frac 0.5                 --min_entropy 1.2                 --min_dust 0.08                --qc_output_file results/03.classifier/rmhost_classified_qc/4_HT29Cells_BMixB_Fn_Ec_16S_S4/4_HT29Cells_BMixB_Fn_Ec_16S_S4_krak_sample_denosing.txt                 --raw_qc_output_file results/03.classifier/rmhost_classified_qc/4_HT29Cells_BMixB_Fn_Ec_16S_S4/4_HT29Cells_BMixB_Fn_Ec_16S_S4_krak_sample_raw.txt                 --barcode_tag CB       --cluster   /data/microcat_debug/microcat_singlecell/fnec/fnec_groups.tsv  --log_file logs/03.classifier/classified_qc/4_HT29Cells_BMixB_Fn_Ec_16S_S4/4_HT29Cells_BMixB_Fn_Ec_16S_S4_krak_sample_denosing.log;


# python /data/project/host-microbiome/MicroCAT/microcat/single_wf/scripts/krak_sample_denosing.py                 --krak_report results/03.classifier/rmhost_kraken2_report/custom/1_HT29_Cells_no_bacteria_16S_S1/1_HT29_Cells_no_bacteria_16S_S1_kraken2_report.txt                 --krak_output results/03.classifier/rmhost_kraken2_output/1_HT29_Cells_no_bacteria_16S_S1/1_HT29_Cells_no_bacteria_16S_S1_kraken2_output.txt                 --krak_mpa_report results/03.classifier/rmhost_kraken2_report/mpa/1_HT29_Cells_no_bacteria_16S_S1/1_HT29_Cells_no_bacteria_16S_S1_kraken2_mpa_report.txt                 --bam results/03.classifier/rmhost_extracted_classified_output/1_HT29_Cells_no_bacteria_16S_S1/1_HT29_Cells_no_bacteria_16S_S1_kraken2_extracted_classified.bam                 --ktaxonomy /data/database/kraken2uniq_database/k2_pluspf_16gb_20231009/ktaxonomy.tsv                --inspect /data/database/kraken2uniq_database/k2_pluspf_16gb_20231009/inspect.txt                 --min_frac 0.5                 --min_entropy 1.2                 --min_dust 0.08                --qc_output_file results/03.classifier/rmhost_classified_qc/1_HT29_Cells_no_bacteria_16S_S1/1_HT29_Cells_no_bacteria_16S_S1_krak_sample_denosing.txt                 --raw_qc_output_file results/03.classifier/rmhost_classified_qc/1_HT29_Cells_no_bacteria_16S_S1/1_HT29_Cells_no_bacteria_16S_S1_krak_sample_raw.txt    --cluster  /data/microcat_debug/microcat_singlecell/16s_louvain_groups.tsv            --barcode_tag CB                 --log_file logs/03.classifier/classified_qc/1_HT29_Cells_no_bacteria_16S_S1/1_HT29_Cells_no_bacteria_16S_S1_krak_sample_denosing.log;


# python  /data/project/host-microbiome/MicroCAT/microcat/single_wf/scripts/krak_study_denosing_copy.py            --file_list results/03.classifier/rmhost_classified_qc/4_HT29Cells_BMixB_Fn_Ec_16S_S4/4_HT29Cells_BMixB_Fn_Ec_16S_S4_krak_sample_denosing.txt results/03.classifier/rmhost_classified_qc/2_HT29_COCA36F3Fn_GEX_S2/2_HT29_COCA36F3Fn_GEX_S2_krak_sample_denosing.txt results/03.classifier/rmhost_classified_qc/1_HT29_Cells_no_bacteria_GEX_S1/1_HT29_Cells_no_bacteria_GEX_S1_krak_sample_denosing.txt results/03.classifier/rmhost_classified_qc/1_HT29_Cells_no_bacteria_16S_S1/1_HT29_Cells_no_bacteria_16S_S1_krak_sample_denosing.txt results/03.classifier/rmhost_classified_qc/2_HT29_COCA36F3Fn_16S_S2/2_HT29_COCA36F3Fn_16S_S2_krak_sample_denosing.txt \
# --raw_file_list /data/scRNA_analysis/benchmark/Galeano2022_HT29/results/03.classifier/rmhost_classified_qc/1_HT29_Cells_no_bacteria_16S_S1/1_HT29_Cells_no_bacteria_16S_S1_krak_sample_raw.txt /data/scRNA_analysis/benchmark/Galeano2022_HT29/results/03.classifier/rmhost_classified_qc/1_HT29_Cells_no_bacteria_GEX_S1/1_HT29_Cells_no_bacteria_GEX_S1_krak_sample_raw.txt /data/scRNA_analysis/benchmark/Galeano2022_HT29/results/03.classifier/rmhost_classified_qc/2_HT29_COCA36F3Fn_16S_S2/2_HT29_COCA36F3Fn_16S_S2_krak_sample_raw.txt /data/scRNA_analysis/benchmark/Galeano2022_HT29/results/03.classifier/rmhost_classified_qc/2_HT29_COCA36F3Fn_GEX_S2/2_HT29_COCA36F3Fn_GEX_S2_krak_sample_raw.txt \
# --out_path results/03.classifier/rmhost_classified_qc/2_HT29_COCA36F3Fn_16S_S2/2_HT29_COCA36F3Fn_16S_S2_krak_study_denosing.txt /data/scRNA_analysis/benchmark/Galeano2022_HT29/results/03.classifier/rmhost_classified_qc/4_HT29Cells_BMixB_Fn_Ec_16S_S4/4_HT29Cells_BMixB_Fn_Ec_16S_S4_krak_sample_raw.txt             --sample_name 2_HT29_COCA36F3Fn_16S_S2             --min_reads 2             --min_uniq 2             --log_file logs/03.classifier/classified_qc/2_HT29_COCA36F3Fn_16S_S2/2_HT29_COCA36F3Fn_16S_S2_krak_study_denosing.log