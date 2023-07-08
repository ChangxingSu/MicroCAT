# import pandas as pd
# import glob
# import os
# ## beta test

# import sample


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

# def aggregate_RG_bam_output(wildcards):
#     checkpoint_output = checkpoints.starsolo_smartseq_demultiplex_bam_by_read_group.get(**wildcards).output.unmapped_bam_demultiplex_dir

#     sample_folders = glob_wildcards(os.path.join(config["output"]["host"], "unmapped_host/{sample}"))
#     samples = sample_folders.sample

#     return expand(
#         os.path.join(
#             config["output"]["host"],
#             "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"
#         ),
#         sample=samples
#     )


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

        if config["params"]["host"]["starsolo"]["assay"]=="tenX_auto":
            rule starsolo_10x_count:
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
                        "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam")
                resources:
                    mem_mb=100000  # This rule needs 100 GB of memory
                params:
                    barcode_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                    cdna_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq2"),
                    starsolo_out = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/"),
                    reference = config["params"]["host"]["starsolo"]["reference"],
                    variousParams = config["params"]["host"]["starsolo"]["variousParams"],
                    threads = config["params"]["host"]["starsolo"]["threads"],
                    barcode_data_dir = config["datas"]["barcode_list_dirs"]["tenX"]
                log:
                    os.path.join(config["logs"]["host"],
                                "starsolo/{sample}_starsolo_count.log")
                benchmark:
                    os.path.join(config["benchmarks"]["host"],
                                "starsolo/{sample}_starsolo_count.benchmark")
                conda:
                    config["envs"]["star"]
                message: "Executing starsolo with {params.threads} threads on the following files {wildcards.sample}.Library with 10x 3' V3"
                shell:
                    '''
                    mkdir -p {params.starsolo_out}; 
                    cd {params.starsolo_out} ;

                    if echo {params.cdna_reads} | grep -q "\.gz" ; then
                        file_command='--readFilesCommand zcat'
                    else
                        file_command=''
                    fi

                    ## also define one file from R1/R2; we choose the largest one, because sometimes there are tiny files from trial runs
                    R1F=`echo {params.barcode_reads} | tr ',' ' ' | xargs ls -s | tail -n1 | awk '{{print $2}}'`
                    R2F=`echo {params.cdna_reads} | tr ',' ' ' | xargs ls -s | tail -n1 | awk '{{print $2}}'`

                    ## let's see if the files are archived or not. Gzip is tprinthe most common, but bgzip archives should work too since they are gzip-compatible.
                    GZIP=""
                    BC=""
                    CHEMISTRY=""
                    NBC1=""
                    NBC2=""
                    NBC3=""
                    NBCA=""
                    R1LEN=""
                    R2LEN=""
                    R1DIS=""
                    WL="{params.barcode_data_dir}"
                    
                    echo $R2F 
                    
                    ## randomly subsample 200k reads - let's hope there are at least this many (there should be):
                    seqtk sample -s100 $R1F 200000 > {wildcards.sample}.test.R1.fastq 
                    seqtk sample -s100 $R2F 200000 > {wildcards.sample}.test.R2.fastq 
                    wait
                    
                    NBC1=`cat {wildcards.sample}.test.R1.fastq | awk 'NR%4==2' | grep -F -f $WL/737K-april-2014_rc.txt | wc -l`
                    NBC2=`cat {wildcards.sample}.test.R1.fastq | awk 'NR%4==2' | grep -F -f $WL/737K-august-2016.txt | wc -l`
                    NBC3=`cat {wildcards.sample}.test.R1.fastq | awk 'NR%4==2' | grep -F -f $WL/3M-february-2018.txt | wc -l`
                    NBCA=`cat {wildcards.sample}.test.R1.fastq | awk 'NR%4==2' | grep -F -f $WL/737K-arc-v1.txt | wc -l`
                    R1LEN=`cat {wildcards.sample}.test.R1.fastq | awk 'NR%4==2' | awk '{{sum+=length($0)}} END {{printf "%d\\n",sum/NR+0.5}}'`
                    R2LEN=`cat {wildcards.sample}.test.R2.fastq | awk 'NR%4==2' | awk '{{sum+=length($0)}} END {{printf "%d\\n",sum/NR+0.5}}'`
                    R1DIS=`cat {wildcards.sample}.test.R1.fastq | awk 'NR%4==2' | awk '{{print length($0)}}' | sort | uniq -c | wc -l`

                    ## elucidate the right barcode whitelist to use. Grepping out N saves us some trouble. Note the special list for multiome experiments (737K-arc-v1.txt):
                    ## 80k (out of 200,000) is an empirical number - I've seen <50% barcodes matched to the whitelist, but a number that's < 40% suggests something is very wrong
                    max_count=0
                    max_count_variable=""

                    if (( $NBC3 > 80000 )) 
                    then
                        if (( $NBC3 > max_count ))
                        then
                            max_count=$NBC3
                            max_count_variable="NBC3"
                            BC=$WL/3M-february-2018.txt
                            CHEMISTRY="V3"
                        fi
                    fi

                    if (( $NBC2 > 80000 ))
                    then
                        if (( $NBC2 > max_count ))
                        then
                            max_count=$NBC2
                            max_count_variable="NBC2"
                            BC=$WL/737K-august-2016.txt
                            CHEMISTRY="V2"
                        fi
                    fi

                    if (( $NBCA > 80000 ))
                    then
                        if (( $NBCA > max_count ))
                        then
                            max_count=$NBCA
                            max_count_variable="NBCA"
                            BC=$WL/737K-arc-v1.txt
                            CHEMISTRY="V1"
                        fi
                    fi

                    if (( $NBC1 > 80000 )) 
                    then
                        if (( $NBC1 > max_count ))
                        then
                            max_count=$NBC1
                            max_count_variable="NBC1"
                            CHEMISTRY="multiome"
                            BC=$WL/737K-april-2014_rc.txt
                        fi
                    fi

                    if [[ -n "$max_count_variable" ]]
                    then
                        echo "The variable with the highest count is $max_count_variable"
                        # Do something with the selected variables BC and CHEMISTRY
                    else
                        >&2 echo "ERROR: No whitelist has matched a random selection of 200,000 barcodes! Match counts: $$NBC1 (v1), $$NBC2 (v2), $$NBC3 (v3), $$NBCA (multiome)."
                        exit 1
                    fi

                    ## check read lengths, fail if something funky is going on: 
                    PAIRED=False
                    UMILEN=""
                    CBLEN=""
                    if (( $R1DIS > 1 && $R1LEN <= 30 ))
                    then 
                        >&2 echo "ERROR: Read 1 (barcode) has varying length; possibly someone thought it's a good idea to quality-trim it. Please check the fastq files."
                        exit 1
                    elif (( $R1LEN < 24 )) 
                    then
                        >&2 echo "ERROR: Read 1 (barcode) is less than 24 bp in length. Please check the fastq files."
                        exit 1
                    elif (( $R2LEN < 40 )) 
                    then
                        >&2 echo "ERROR: Read 2 (biological read) is less than 40 bp in length. Please check the fastq files."
                        exit 1
                    fi


                    ## assign the necessary variables for barcode/UMI length/paired-end processing. 
                    ## scripts was changed to not rely on read length for the UMIs because of the epic Hassan case
                    # (v2 16bp barcodes + 10bp UMIs were sequenced to 28bp, effectively removing the effects of the UMIs)
                    if (( $R1LEN > 50 )) 
                    then
                        PAIRED=True
                    fi

                    if [[ $CHEMISTRY == "V3" || $CHEMISTRY == "multiome" ]] 
                    then 
                        CBLEN=16
                        UMILEN=12
                    elif [[ $CHEMISTRY == "V2" ]] 
                    then
                        CBLEN=16
                        UMILEN=10
                    elif [[ $CHEMISTRY == "V1" ]] 
                    then
                        CBLEN=14
                        UMILEN=10
                    fi 

                    ## finally, see if you have 5' or 3' experiment. I don't know and easier way than to run a test alignment:  
                    STRAND=Forward

                    STAR --runThreadN {params.threads} --genomeDir {params.reference} --readFilesIn {wildcards.sample}.test.R2.fastq {wildcards.sample}.test.R1.fastq --runDirPerm All_RWX --outSAMtype None \
                        --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) \
                        --soloUMIlen $UMILEN --soloStrand Forward \
                        --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
                        --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
                        --soloFeatures Gene GeneFull --soloOutFileNames {wildcards.sample}_test_strand/ features.tsv barcodes.tsv matrix.mtx &> /dev/null 

                    ## the following is needed in case of bad samples: when a low fraction of reads come from mRNA, experiment will look falsely reverse-stranded
                    UNIQFRQ=`grep "Reads Mapped to Genome: Unique," {wildcards.sample}_test_strand/GeneFull/Summary.csv | awk -F "," '{{print $2}}'`
                    GENEPCT=`grep "Reads Mapped to GeneFull: Unique GeneFull" {wildcards.sample}_test_strand/GeneFull/Summary.csv | awk -F "," -v v=$UNIQFRQ '{{printf "%d\\n",$2*100/v}}'`

                    ## this percentage is very empirical, but was found to work in 99% of cases. 
                    ## any 10x 3' run with GENEPCT < 35%, and any 5' run with GENEPCT > 35% are 
                    ## *extremely* strange and need to be carefully evaluated
                    if (( $GENEPCT < 35 )) 
                    then
                        STRAND=Reverse
                    fi

                    ## finally, if paired-end experiment turned out to be 3' (yes, they do exist!), process it as single-end: 
                    if [[ $STRAND == "Forward" && $PAIRED == "True" ]]
                    then
                        PAIRED=False
                    fi

                    echo "Done setting up the STARsolo run; here are final processing options:" >> ../../../{log}
                    echo "=============================================================================" >> ../../../{log}
                    echo "Paired-end mode: $PAIRED" >> ../../../{log}
                    echo "Strand (Forward = 3', Reverse = 5'): $STRAND, %reads same strand as gene: $GENEPCT" >> ../../../{log}
                    echo "CB whitelist: $BC, matches out of 200,000: $NBC3 (v3), $NBC2 (v2), $NBC1 (v1), $NBCA (multiome) " >> ../../../{log}
                    echo "CB length: $CBLEN" >> ../../../{log}
                    echo "UMI length: $UMILEN" >> ../../../{log}
                    echo "-----------------------------------------------------------------------------" >> ../../../{log}
                    echo "Read 1 files: {params.barcode_reads}" >> ../../../{log}
                    echo "-----------------------------------------------------------------------------" >> ../../../{log}
                    echo "Read 2 files: {params.cdna_reads}" >> ../../../{log}
                    echo "-----------------------------------------------------------------------------" >> ../../../{log}


                    if [[ $PAIRED == "True" ]]
                    then
                        ## note the R1/R2 order of input fastq reads and --soloStrand Forward for 5' paired-end experiment
                        STAR --runThreadN {params.threads} --genomeDir {params.reference} --readFilesIn {params.barcode_reads} {params.cdna_reads} --runDirPerm All_RWX $file_command \
                        --soloBarcodeMate 1 --clip5pNbases 39 0 \
                        --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloCBstart 1 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand Forward \
                        --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
                        --soloCellFilter EmptyDrops_CR --outFilterScoreMin 30 \
                        --soloMultiMappers EM \
                        --outSAMtype BAM SortedByCoordinate \
                        --outSAMattrRGline ID:{wildcards.sample} PL:illumina SM:{wildcards.sample} LB:$CHEMISTRY \
                        --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN \
                        --outFileNamePrefix ./{wildcards.sample}/\
                        {params.variousParams}  \
                        2>&1 | tee ../../../{log} ;
                    else 
                        STAR --runThreadN {params.threads} --genomeDir {params.reference} --readFilesIn {params.cdna_reads} {params.barcode_reads} --runDirPerm All_RWX $file_command \
                        --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand $STRAND \
                        --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
                        --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
                        --soloMultiMappers EM \
                        --outSAMtype BAM SortedByCoordinate \
                        --outSAMattrRGline ID:{wildcards.sample} PL:illumina SM:{wildcards.sample} LB:$CHEMISTRY \
                        --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN \
                        --outFileNamePrefix ./{wildcards.sample}/\
                        {params.variousParams}  \
                        2>&1 | tee ../../../{log} ;
                    fi

                    rm -rf {wildcards.sample}.test.R1.fastq
                    rm -rf {wildcards.sample}.test.R2.fastq
                    rm -rf {wildcards.sample}_test_strand
                    
                    cd ../../../;\
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/features.tsv" "{output.features_file}" ;
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/matrix.mtx" "{output.matrix_file}" ; 
                    ln -sr "{params.starsolo_out}/{wildcards.sample}/Solo.out/Gene/filtered/barcodes.tsv" "{output.barcodes_file}" ;\
                    mv "{params.starsolo_out}/{wildcards.sample}/Aligned.sortedByCoord.out.bam" "{output.mapped_bam_file}";\


                    '''

        if config["params"]["host"]["starsolo"]["assay"]=="tenX_v3":
            rule starsolo_10x_count:
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
                        "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam")
                resources:
                    mem_mb=100000  # This rule needs 100 GB of memory
                params:
                    barcode_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                    cdna_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq2"),
                    starsolo_out = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/"),
                    reference = config["params"]["host"]["starsolo"]["reference"],
                    variousParams = config["params"]["host"]["starsolo"]["variousParams"],
                    threads = config["params"]["host"]["starsolo"]["threads"],
                    barcode_list =  os.path.join(config["datas"]["barcode_list_dirs"]["tenX"],"3M-february-2018.txt")
                log:
                    os.path.join(config["logs"]["host"],
                                "starsolo/{sample}_starsolo_count.log")
                benchmark:
                    os.path.join(config["benchmarks"]["host"],
                                "starsolo/{sample}_starsolo_count.benchmark")
                conda:
                    config["envs"]["star"]
                message: "Executing starsolo with {params.threads} threads on the following files {wildcards.sample}.Library with 10x 3' V3"
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
                    --soloType CB_UMI_Simple \
                    --soloCBwhitelist {params.barcode_list} \
                    --soloCBstart 1 \
                    --soloCBlen 16 \
                    --soloUMIstart 17 \
                    --soloUMIlen 12 \
                    --genomeDir {params.reference} \
                    --readFilesIn {params.cdna_reads} {params.barcode_reads} \
                    --runThreadN {params.threads} \
                    --clipAdapterType CellRanger4 \
                    --outFilterScoreMin 30 \
                    $file_command  \
                    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
                    --soloUMIfiltering MultiGeneUMI_CR \
                    --outSAMtype BAM SortedByCoordinate\
                    --outSAMattrRGline ID:{wildcards.sample} PL:illumina SM:{wildcards.sample} LB:tenX_v3 \
                    --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN \
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
        if config["params"]["host"]["starsolo"]["assay"]=="tenX_v1":
            rule starsolo_10x_count:
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
                        "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam")
                params:
                    barcode_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                    cdna_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq2"),
                    starsolo_out = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/"),
                    reference = config["params"]["host"]["starsolo"]["reference"],
                    variousParams = config["params"]["host"]["starsolo"]["variousParams"],
                    threads = config["params"]["host"]["starsolo"]["threads"],
                    barcode_list =  os.path.join(config["datas"]["barcode_list_dirs"]["tenX"],"737K-april-2014_rc.txt")
                log:
                    os.path.join(config["logs"]["host"],
                                "starsolo/{sample}_starsolo_count.log")
                benchmark:
                    os.path.join(config["benchmarks"]["host"],
                                "starsolo/{sample}_starsolo_count.benchmark")
                conda:
                    config["envs"]["star"]
                message: "Executing starsolo with {params.threads} threads on the following files {wildcards.sample}.Library with 10x 3' V1"
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
                    --soloType CB_UMI_Simple \
                    --soloCBwhitelist {params.barcode_list} \
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
                    $file_command  \
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
                        "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam")
                params:
                    starsolo_out = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/"),
                    reference = config["params"]["host"]["starsolo"]["reference"],
                    variousParams = config["params"]["host"]["starsolo"]["variousParams"],
                    threads = config["params"]["host"]["starsolo"]["threads"],
                    barcode_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                    cdna_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq2"),
                    barcode_list =  os.path.join(config["datas"]["barcode_list_dirs"]["tenX"],"737K-august-2016.txt")
                log:
                    os.path.join(config["logs"]["host"],
                                "starsolo/{sample}_starsolo_count.log")
                benchmark:
                    os.path.join(config["benchmarks"]["host"],
                                "starsolo/{sample}_starsolo_count.benchmark")
                conda:
                    config["envs"]["star"]
                resources:
                    threads=20,      # This rule needs 30 threads
                    mem_mb=100000  # This rule needs 100 GB of memory
                message: "Executing starsolo with {params.threads} threads on the following files {wildcards.sample}.Library with 10x 3' V2."
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
                    --soloType CB_UMI_Simple \
                    --soloCBwhitelist {params.barcode_list} \
                    --soloCBstart 1 \
                    --soloCBlen 16 \
                    --soloUMIstart 17 \
                    --soloUMIlen 10 \
                    --genomeDir {params.reference} \
                    --readFilesIn {params.cdna_reads}  {params.barcode_reads} \
                    --runThreadN {params.threads} \
                    --clipAdapterType CellRanger4 \
                    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
                    --soloUMIfiltering MultiGeneUMI_CR \
                    --outSAMtype BAM SortedByCoordinate\
                    --outFilterScoreMin 30  \
                    --outSAMattrRGline ID:{wildcards.sample} PL:illumina SM:{wildcards.sample} LB:tenX_v2 \
                    --outSAMattributes CR UR CY UY CB UB\
                    $file_command \
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
                    mv "{params.starsolo_out}/{wildcards.sample}/Aligned.sortedByCoord.out.bam" "{output.mapped_bam_file}"\
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
            conda:
                config["envs"]["star"]
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

    if config["params"]["host"]["starsolo"]["assay"]=="Cel-Seq2":
            rule starsolo_Cel_count:
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
                        "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam")
                params:
                    starsolo_out = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/"),
                    reference = config["params"]["host"]["starsolo"]["reference"],
                    variousParams = config["params"]["host"]["starsolo"]["variousParams"],
                    threads = config["params"]["host"]["starsolo"]["threads"],
                    barcode_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                    cdna_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq2"),
                    barcode_list =  os.path.join(config["datas"]["barcode_list_dirs"]["tenX"],"737K-august-2016.txt")
                log:
                    os.path.join(config["logs"]["host"],
                                "starsolo/{sample}_starsolo_count.log")
                benchmark:
                    os.path.join(config["benchmarks"]["host"],
                                "starsolo/{sample}_starsolo_count.benchmark")
                conda:
                    config["envs"]["star"]
                resources:
                    threads=20,      # This rule needs 30 threads
                    mem_mb=100000  # This rule needs 100 GB of memory
                message: "Executing starsolo with {params.threads} threads on the following files {wildcards.sample}.Library with 10x 3' V2."
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
                    --soloType CB_UMI_Simple \
                    --soloCBwhitelist {params.barcode_list} \
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
                    --soloUMIdedup 1MM_CR \
                    $file_command \
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

        rule generate_pe_manifest_file:
            input:
                config["params"]["samples"],
                fastqs_dir=lambda wildcards: microcat.get_fastqs_dir(SAMPLES,wildcards),
            output:
                PE_MANIFEST_FILE = "./manifest.tsv"
            script:
                "../scripts/generate_PE_manifest_file.py"
        
        rule starsolo_smartseq_count:
            # Input files
            input:
                # Path to the input manifest file
                manifest = "./manifest.tsv",
                fastqs_dir=lambda wildcards: microcat.get_fastqs_dir(SAMPLES,wildcards),
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
                SAMattrRGline = microcat.get_SAMattrRGline_from_manifest(config["params"]["host"]["starsolo"]["manifest"]),
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
            conda:
                config["envs"]["star"]
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
            conda:
                config["envs"]["star"]
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
            conda:
                config["envs"]["star"]
            shell:
                '''
                samtools sort -@ {params.threads} -t {params.tag} -o {output.unmapped_sorted_bam}  {input.unmapped_bam_unsorted_file};
                '''

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
        rule split_starsolo_BAM_by_RG:
            input:
                unmapped_sorted_bam = os.path.join(
                        config["output"]["host"],
                        "starsolo_count/Aligned_out_unmapped_RGsorted.bam")
            output:
                unmapped_bam_sorted_file =os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
            params:
                SampleID="{sample}",
            shell:
                '''
                python /data/project/host-microbiome/microcat/microcat/scripts/split_Starsolo_BAM_by_RG.py \
                --bam_path {input.unmapped_sorted_bam} \
                --tag {params.SampleID} \
                --output_bam {output.unmapped_bam_sorted_file} 
                '''

        rule starsolo_smartseq_all:
            input:
                # directory(os.path.join(
                #     config["output"]["host"],
                #     "unmapped_host/"))
                expand(os.path.join(config["output"]["host"],"unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),sample=SAMPLES_ID_LIST)
        
    else:
        rule starsolo_smartseq_all:
            input: 
    #ALL Input

    if config["params"]["host"]["starsolo"]["assay"]=="Seq-well":
        rule starsolo_seqwell_count:
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
                    "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam")
            params:
                starsolo_out = os.path.join(
                    config["output"]["host"],
                    "starsolo_count/"),
                reference = config["params"]["host"]["starsolo"]["reference"],
                variousParams = config["params"]["host"]["starsolo"]["variousParams"],
                threads = config["params"]["host"]["starsolo"]["threads"],
                barcode_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                cdna_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq2"),
            log:
                os.path.join(config["logs"]["host"],
                            "starsolo/{sample}_starsolo_count.log")
            benchmark:
                os.path.join(config["benchmarks"]["host"],
                            "starsolo/{sample}_starsolo_count.benchmark")
            conda:
                config["envs"]["star"]
            resources:
                threads=20,      # This rule needs 30 threads
                mem_mb=100000  # This rule needs 100 GB of memory
            message: "Executing starsolo with {params.threads} threads on the following files {wildcards.sample}.Library with 10x 3' V2."
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
                --soloType CB_UMI_Simple \
                --soloCBwhitelist None \
                --soloCBstart 1 \
                --soloCBlen 8 \
                --soloUMIstart 8 \
                --soloUMIlen 12 \
                --genomeDir {params.reference} \
                --readFilesIn {params.cdna_reads}  {params.barcode_reads} \
                --runThreadN {params.threads} \
                --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
                --soloUMIfiltering MultiGeneUMI_CR \
                --outSAMtype BAM SortedByCoordinate\
                --outFilterScoreMin 30  \
                --outSAMattrRGline ID:{wildcards.sample} PL:illumina SM:{wildcards.sample} LB:Seq-Well \
                --outSAMattributes CR UR CY UY CB UB\
                $file_command \
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
        rule starsolo_seqwell_unmapped_extracted_sorted:
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
            conda:
                config["envs"]["star"]
            shell:
                '''
                samtools view --threads  {params.threads}  -b -f 4   {input.mapped_bam_file}  >  {params.unmapped_bam_unsorted_file};\
                samtools sort -n  --threads  {params.threads} {params.unmapped_bam_unsorted_file} -o {output.unmapped_bam_sorted_file}
                '''
        rule starsolo_seqwell_all:
            input:
                expand(os.path.join(config["output"]["host"],"unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),sample=SAMPLES_ID_LIST)
    else:
        rule starsolo_seqwell_all:
            input: 


    if config["params"]["host"]["starsolo"]["assay"]=="Seq-well":
        rule starsolo_seqwell_count:
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
                    "starsolo_count/{sample}/Aligned_sortedByCoord_out.bam")
            params:
                starsolo_out = os.path.join(
                    config["output"]["host"],
                    "starsolo_count/"),
                reference = config["params"]["host"]["starsolo"]["reference"],
                variousParams = config["params"]["host"]["starsolo"]["variousParams"],
                threads = config["params"]["host"]["starsolo"]["threads"],
                barcode_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq1"),
                cdna_reads = lambda wildcards: microcat.get_starsolo_sample_id(SAMPLES, wildcards, "fq2"),
            log:
                os.path.join(config["logs"]["host"],
                            "starsolo/{sample}_starsolo_count.log")
            benchmark:
                os.path.join(config["benchmarks"]["host"],
                            "starsolo/{sample}_starsolo_count.benchmark")
            conda:
                config["envs"]["star"]
            resources:
                threads=20,      # This rule needs 30 threads
                mem_mb=100000  # This rule needs 100 GB of memory
            message: "Executing starsolo with {params.threads} threads on the following files {wildcards.sample}.Library with 10x 3' V2."
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
                --soloType CB_UMI_Simple \
                --soloCBwhitelist None \
                --soloCBstart 1 \
                --soloCBlen 8 \
                --soloUMIstart 8 \
                --soloUMIlen 12 \
                --genomeDir {params.reference} \
                --readFilesIn {params.cdna_reads}  {params.barcode_reads} \
                --runThreadN {params.threads} \
                --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
                --soloUMIfiltering MultiGeneUMI_CR \
                --outSAMtype BAM SortedByCoordinate\
                --outFilterScoreMin 30  \
                --outSAMattrRGline ID:{wildcards.sample} PL:illumina SM:{wildcards.sample} LB:Seq-Well \
                --outSAMattributes CR UR CY UY CB UB\
                $file_command \
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


    rule starsolo_all:
        input: 
            rules.starsolo_smartseq_all.input,
            rules.starsolo_10X_all.input,
            rules.starsolo_seqwell_all.input,

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
            ln -sr "{params.cr_out}{params.SampleID}/outs/web_summary.html" "{params.web_summary}" ; 
            ln -sr "{params.cr_out}{params.SampleID}/outs/metrics_summary.csv" "{params.metrics_summary}";
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
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam")
        params:
            threads=16,
            unmapped_bam_unsorted_file = os.path.join(
                config["output"]["host"],
                "unmapped_host/{sample}/Aligned_sortedByCoord_unmapped_out.bam")
        ## because bam is sorted by Coord,it's necessary to sort it by read name
        conda:
            config["envs"]["star"]
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
            expand( os.path.join(
                    config["output"]["host"],
                    "unmapped_host/{sample}/Aligned_sortedByName_unmapped_out.bam"),sample=SAMPLES_ID_LIST)

else:
    rule cellranger_all:
        input:

rule host_all:
    input:
        rules.starsolo_all.input,
        rules.cellranger_all.input,



