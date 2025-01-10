import os
import subprocess

# Rule to prepare and format FASTQ files for spaceranger processing
rule format_fastq_for_spaceranger:
    input:
        config["params"]["samples"],
    output:
        os.path.join(config["output"]["raw"], "reads/{sample}/{sample}_summary.json")
    script:
        "../scripts/preprocess_raw.py"

# Main rule to run spaceranger count for spatial transcriptomics analysis
rule spaceranger_count:
        input:
            sample_summary = os.path.join(config["output"]["raw"], "reads/{sample}/{sample}_summary.json"),
        output:
            # Basic matrix files for gene expression data
            features_file = os.path.join(
                config["output"]["host"],
                "spaceranger_count/{sample}/{sample}_features.tsv"),
            matrix_file = os.path.join(
                config["output"]["host"], 
                "spaceranger_count/{sample}/{sample}_matrix.mtx"),
            barcodes_file = os.path.join(
                config["output"]["host"],
                "spaceranger_count/{sample}/{sample}_barcodes.tsv"),
            
            # Report and statistics files
            metrics_summary = os.path.join(
                config["output"]["host"],
                "spaceranger_count/{sample}/{sample}.metrics_summary.csv"),
            web_summary = os.path.join(
                config["output"]["host"],
                "spaceranger_count/{sample}/{sample}.web_summary.html"),
                
            # Alignment files
            mapped_bam = os.path.join(
                config["output"]["host"],
                "spaceranger_count/{sample}/{sample}_possorted_genome_bam.bam"),
            mapped_bam_index = os.path.join(
                config["output"]["host"],
                "spaceranger_count/{sample}/{sample}_possorted_genome_bam.bam.bai"),
                
            # Spatial-specific files for Visium data
            spatial_scale = os.path.join(
                config["output"]["host"],
                "spaceranger_count/{sample}/scalefactors_json.json"),
            tissue_positions = os.path.join(
                config["output"]["host"],
                "spaceranger_count/{sample}/tissue_positions.csv"),
            spatial_enrichment = os.path.join(
                config["output"]["host"],
                "spaceranger_count/{sample}/spatial_enrichment.csv"),
        params:
            sr_out = os.path.join(
                config["output"]["host"],
                "spaceranger_count/"),
            reference = config["params"]["host"]["spaceranger"]["reference"],
            fastqs_dir = os.path.abspath(os.path.join(config["output"]["raw"], "reads/{sample}/")),
            SampleID="{sample}",
            variousParams = config["params"]["host"]["spaceranger"]["variousParams"],
        threads:
            config["resources"]["spaceranger"]["threads"]
        resources:
            mem_mb=config["resources"]["spaceranger"]["mem_mb"]
        log:
            os.path.join(config["logs"]["host"],
                        "spaceranger/{sample}_spaceranger_count.log")
        benchmark:
            os.path.join(config["benchmarks"]["host"],
                        "spaceranger/{sample}_spaceranger_count.benchmark")
        shell:
            """
            # Run spaceranger count with specified parameters
            spaceranger count \
            --id={params.SampleID} \
            --transcriptome={params.reference} \
            --fastqs={params.fastqs_dir} \
            --sample={params.SampleID} \
            --localcores={threads} \
            --localmem=$(({resources.mem_mb}/1024)) \
            --unknown-slide \
            --reorient-images \
            {params.variousParams} \
            2>&1 | tee ../../../{log} ;  
            
            cd ../../../;
            
            # Copy and decompress filtered feature matrix files
            cp {params.sr_out}{params.SampleID}/outs/filtered_feature_bc_matrix/features.tsv.gz "{params.sr_out}{params.SampleID}/outs/features.tsv.gz";
            cp {params.sr_out}{params.SampleID}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz "{params.sr_out}{params.SampleID}/outs/barcodes.tsv.gz"; 
            cp {params.sr_out}{params.SampleID}/outs/filtered_feature_bc_matrix/matrix.mtx.gz "{params.sr_out}{params.SampleID}/outs/matrix.mtx.gz"; 
            
            # Decompress the copied files
            gunzip "{params.sr_out}{params.SampleID}/outs/matrix.mtx.gz"; 
            gunzip "{params.sr_out}{params.SampleID}/outs/features.tsv.gz";
            gunzip "{params.sr_out}{params.SampleID}/outs/barcodes.tsv.gz"; 
            
            # Move decompressed feature matrix files to output locations
            mv "{params.sr_out}{params.SampleID}/outs/features.tsv" "{output.features_file}"; 
            mv "{params.sr_out}{params.SampleID}/outs/matrix.mtx" "{output.matrix_file}"; 
            
            # Create symbolic links for output files
            ln -sr "{params.sr_out}{params.SampleID}/outs/web_summary.html" "{output.web_summary}";
            ln -sr "{params.sr_out}{params.SampleID}/outs/metrics_summary.csv" "{output.metrics_summary}";
            ln -sr "{params.sr_out}{params.SampleID}/outs/possorted_genome_bam.bam" "{output.mapped_bam}";
            ln -sr "{params.sr_out}{params.SampleID}/outs/possorted_genome_bam.bam.bai" "{output.mapped_bam_index}";
            ln -sr "{params.sr_out}{params.SampleID}/outs/spatial/scalefactors_json.json" "{output.spatial_scale}";
            ln -sr "{params.sr_out}{params.SampleID}/outs/spatial/tissue_positions.csv" "{output.tissue_positions}";
            ln -sr "{params.sr_out}{params.SampleID}/outs/spatial/spatial_enrichment.csv" "{output.spatial_enrichment}";
            
            # Move the decompressed barcodes file
            mv "{params.sr_out}{params.SampleID}/outs/barcodes.tsv" "{output.barcodes_file}";
            """

# Rule to extract and sort unmapped reads from spaceranger BAM output
rule spaceranger_unmapped_extracted_sorted:
    input:
        mapped_bam = os.path.join(config["output"]["host"], "spaceranger_count/{sample}/{sample}_possorted_genome_bam.bam")
    output:
        unmapped_bam_sorted_file = os.path.join(config["output"]["host"], "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),
        unmapped_bam_sorted_index = os.path.join(config["output"]["host"], "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam.bai"),
    params:
        unmapped_bam_unsorted_file = os.path.join(
            config["output"]["host"],
            "unmapped_host/{sample}/Aligned_sortedByCoord_unmapped_out.bam")
    resources:
        mem_mb=config["resources"]["samtools_extract"]["mem_mb"]
    threads:
        config["resources"]["samtools_extract"]["threads"]
    run:    
        # Extract unmapped reads from the BAM file
        shell(
            'samtools view --threads {threads} -b -f 4 {input.mapped_bam_file} > {params.unmapped_bam_unsorted_file}'
        )

        # Check if the output contains any reads by examining the first 10 lines
        result = subprocess.run(f'samtools view {params.unmapped_bam_unsorted_file} | head -n 10', shell=True, capture_output=True)
        head_output_str = result.stdout.decode('utf-8')
        line_count = len(head_output_str.strip().split('\n'))

        # Raise error if no unmapped reads were found
        if line_count == 0:
            raise ValueError(f"Error: The unmapped BAM unsorted file for sample {wildcards.sample} is empty. Please check your data.")
            
        # Sort the unmapped reads by name and create index
        shell(
            '''
            samtools sort -n --threads {threads} {params.unmapped_bam_unsorted_file} -o {output.unmapped_bam_sorted_file};\
            samtools index -@ {threads} {output.unmapped_bam_sorted_file} -o {output.unmapped_bam_sorted_index};\
            rm -rf {params.unmapped_bam_unsorted_file};
            '''
        )

# Rule to aggregate all outputs
rule spaceranger_all:
        input:
            expand(os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"), sample=SAMPLES_ID_LIST), # Unmapped reads
            expand(os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bai"), sample=SAMPLES_ID_LIST) # Unmapped reads index
