
rule raw_prepare_bam:
    input:
        # Directory containing input fastq files
        bam_dir=lambda wildcards: microcat.get_fastqs_dir(SAMPLES,wildcards),
    params:
        bam_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "bam"),
        sample_id = "{sample}",
    output:
        mapped_bam_file = os.path.join(
            config["output"]["host"],
            "/{sample}/Aligned_sortedByCoord_out.bam")
    shell:
        '''
        ln -sr "{params.bam_reads}" "{output.mapped_bam_file}";
        '''    

rule bam_unmapped_extracted_sorted:
    input:
        mapped_bam_file = os.path.join(
            config["output"]["host"],
            "{sample}/Aligned_sortedByCoord_out.bam")
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
                    "starsolo/{sample}/unmapped_extracted_sorted_bam.log")
    benchmark:
        os.path.join(config["benchmarks"]["host"],
                    "starsolo/{sample}/unmapped_extracted_sorted_bam.benchmark")
    shell:
        '''
        samtools view --threads  {threads}  -b -f 4  {input.mapped_bam_file}  >  {params.unmapped_bam_unsorted_file};\
        samtools sort -n  --threads  {threads} {params.unmapped_bam_unsorted_file} -o {output.unmapped_bam_sorted_file};\
        samtools index -@  {threads} {output.unmapped_bam_sorted_file} -o {output.unmapped_bam_sorted_index};\
        rm -rf {params.unmapped_bam_unsorted_file};
        '''

rule host_all:
    input:
        expand(os.path.join(config["output"]["host"],"unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),sample=SAMPLES_ID_LIST),
        expand(os.path.join(
                config["output"]["host"],
                "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bai"), sample=SAMPLES_ID_LIST)