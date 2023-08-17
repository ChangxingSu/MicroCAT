#!/bin/bash 
# default
barcode_reads=""
cdna_reads=""
barcode_data_dir=""
sample=""
threads=""
reference=""
soloUMIdedup=""
soloCBmatchWLtype=""
soloUMIfiltering=""
soloCellFilter=""
outFilterScoreMin=""
soloMultiMappers=""
clipAdapterType=""
variousParams=""
file_command=""

# 解析命令行选项
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --barcode_reads)
      barcode_reads=$2
      shift 2
      ;;
    --cdna_reads)
      cdna_reads=$2
      shift 2
      ;;
    --barcode_data_dir)
      barcode_data_dir=$2
      shift 2
      ;;
    --sample)
      sample=$2
      shift 2
      ;;
    --threads)
      threads=$2
      shift 2
      ;;
    --reference)
      reference=$2
      shift 2
      ;;
    --soloUMIdedup)
      soloUMIdedup=$2
      shift 2
      ;;
    --soloCBmatchWLtype)
      soloCBmatchWLtype=$2
      shift 2
      ;;
    --soloUMIfiltering)
      soloUMIfiltering=$2
      shift 2
      ;;
    --soloCellFilter)
      soloCellFilter=$2
      shift 2
      ;;
    --outFilterScoreMin)
      outFilterScoreMin=$2
      shift 2
      ;;
    --soloMultiMappers)
      soloMultiMappers=$2
      shift 2
      ;;
    --clipAdapterType)
      clipAdapterType=$2
      shift 2
      ;;
    --variousParams)
      variousParams=$2
      shift 2
      ;;
    *)
      echo "Invalid args: $1" >&2
      exit 1
      ;;
  esac
done

if [ -z "$variousParams" ]; then
  variousParams=""
fi

green=$(tput setaf 2)
blue=$(tput setaf 4)
reset=$(tput sgr0)
red=$(tput setaf 1)

echo "Running preflight checks (please wait)..."


if echo $cdna_reads | grep -q "\.gz" ; then
    file_command='--readFilesCommand zcat'
else
    file_command=''
fi

## also define one file from R1/R2; we choose the largest one, because sometimes there are tiny files from trial runs
R1F=`echo $barcode_reads | tr ',' ' ' | xargs ls -s | tail -n1 | awk '{{print $2}}'`
R2F=`echo $cdna_reads | tr ',' ' ' | xargs ls -s | tail -n1 | awk '{{print $2}}'`




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
WL=$barcode_data_dir


current_date=$(date +%Y-%m-%d)
current_time=$(date +%H:%M:%S)
state="[runtime] (run)"
message="Subsample 200k reads for test barcode"
alignment=100

printf "\e[32m%s %s\e[0m %s \e[34m%-${alignment}s\e[0m\n" "$current_date" "$current_time" "$state" "$message"

# echo -e "\e[32m$current_date $current_time\e[0m [runtime] (run:local)       Subsample 200k reads"

## randomly subsample 200k reads - let's hope there are at least this many (there should be):
# seqtk sample -s100 $R1F 200000 > {$sample}.test.R1.fastq 
# seqtk sample -s100 $R2F 200000 > {$sample}.test.R2.fastq 
seqtk sample -s100 $R1F 200000 > {$sample}.test.R1.fastq $
seqtk sample -s100 $R2F 200000 > {$sample}.test.R2.fastq $
wait

current_date=$(date +%Y-%m-%d)
current_time=$(date +%H:%M:%S)
state="[runtime] (ready)"
message="Successfully subsample 200k reads for test barcode"

printf "%s %s %s %-${alignment}s\n" "${green}${current_date}" "${current_time}" "${state}" "${blue}${message}${reset}"

# echo "$current_date $current_time [runtime] (ready)       Successfully subsample 200k reads"

NBC1=`cat {$sample}.test.R1.fastq | awk 'NR%4==2' | grep -F -f $WL/737K-april-2014_rc.txt | wc -l`
NBC2=`cat {$sample}.test.R1.fastq | awk 'NR%4==2' | grep -F -f $WL/737K-august-2016.txt | wc -l`
NBC3=`cat {$sample}.test.R1.fastq | awk 'NR%4==2' | grep -F -f $WL/3M-february-2018.txt | wc -l`
NBCA=`cat {$sample}.test.R1.fastq | awk 'NR%4==2' | grep -F -f $WL/737K-arc-v1.txt | wc -l`
R1LEN=`cat {$sample}.test.R1.fastq | awk 'NR%4==2' | awk '{sum+=length($0)} END {printf "%d\\n",sum/NR+0.5}'`
R2LEN=`cat {$sample}.test.R2.fastq | awk 'NR%4==2' | awk '{sum+=length($0)} END {printf "%d\\n",sum/NR+0.5}'`
R1DIS=`cat {$sample}.test.R1.fastq | awk 'NR%4==2' | awk '{print length($0)}' | sort | uniq -c | wc -l`

## elucidate the right barcode whitelist to use. Grepping out N saves us some trouble. Note the special list for multiome experiments (737K-arc-v1.txt):
## 80k (out of 200,000) is an empirical number - I've seen <50% barcodes matched to the whitelist, but a number that's < 40% suggests something is very wrong
max_count=0
max_count_variable=""

if (( $NBC3 > 80000 )) 
then
    if (( $NBC3 > max_count ))
    then
        max_count=$NBC3
        max_count_variable="10X v3"
        BC=$WL/3M-february-2018.txt
        CHEMISTRY="tenx_v3"
    fi
fi

if (( $NBC2 > 80000 ))
then
    if (( $NBC2 > max_count ))
    then
        max_count=$NBC2
        max_count_variable="10X v2"
        BC=$WL/737K-august-2016.txt
        CHEMISTRY="tenx_v2"
    fi
fi

if (( $NBCA > 80000 ))
then
    if (( $NBCA > max_count ))
    then
        max_count=$NBCA
        max_count_variable="10X Multiome"
        BC=$WL/737K-arc-v1.txt
        CHEMISTRY="Multiome"
    fi
fi

if (( $NBC1 > 80000 )) 
then
    if (( $NBC1 > max_count ))
    then
        max_count=$NBC1
        max_count_variable="10X v1"
        CHEMISTRY="tenx_v1"
        BC=$WL/737K-april-2014_rc.txt
    fi
fi

current_date=$(date +%Y-%m-%d)
current_time=$(date +%H:%M:%S)
if [[ -n "$max_count_variable" ]]
then
    # echo "The variable with the highest count is $max_count_variable"
    state="[runtime] (ready)"
    message="Successfully detected fastq barcode chemistry: $max_count_variable"
    printf "%s %s %s %-${alignment}s\n" "${green}${current_date}" "${current_time}" "${state}" "${blue}${message}${reset}"

else
    >&2 echo "${red}ERROR${reset}: No whitelist has matched a random selection of 200,000 barcodes! Match counts: $$NBC1 (v1), $$NBC2 (v2), $$NBC3 (v3), $$NBCA (multiome)."
    exit 1
fi

## check read lengths, fail if something funky is going on: 
PAIRED=False
UMILEN=""
CBLEN=""
if (( $R1DIS > 1 && $R1LEN <= 30 ))
then 
    >&2 echo "${red}ERROR${reset}: Read 1 (barcode) has varying length; possibly someone thought it's a good idea to quality-trim it. Please check the fastq files."
    exit 1
elif (( $R1LEN < 24 )) 
then
    >&2 echo "${red}ERROR${reset}: Read 1 (barcode) is less than 24 bp in length. Please check the fastq files."
    exit 1
elif (( $R2LEN < 40 )) 
then
    >&2 echo "${red}ERROR${reset}: Read 2 (biological read) is less than 40 bp in length. Please check the fastq files."
    exit 1
fi


## assign the necessary variables for barcode/UMI length/paired-end processing. 
## scripts was changed to not rely on read length for the UMIs because of the epic Hassan case
# (v2 16bp barcodes + 10bp UMIs were sequenced to 28bp, effectively removing the effects of the UMIs)
if (( $R1LEN > 50 )) 
then
    PAIRED=True
fi

if [[ $CHEMISTRY == "tenx_v3" || $CHEMISTRY == "Multiome" ]] 
then 
    CBLEN=16
    UMILEN=12
elif [[ $CHEMISTRY == "tenx_v2" ]] 
then
    CBLEN=16
    UMILEN=10
elif [[ $CHEMISTRY == "tenx_v1" ]] 
then
    CBLEN=14
    UMILEN=10
fi 

## finally, see if you have 5' or 3' experiment. I don't know and easier way than to run a test alignment:  
STRAND=Forward

report_dir=""
report_dir="${sample}_test_strand"

current_date=$(date +%Y-%m-%d)
current_time=$(date +%H:%M:%S)
state="[runtime] (run)"
message="Detected sequencing strand"

green=$(tput setaf 2)
blue=$(tput setaf 4)
reset=$(tput sgr0)

printf "%s %s %s %-${alignment}s\n" "${green}${current_date}" "${current_time}" "${state}" "${blue}${message}${reset}"

STAR --runThreadN $threads --genomeDir $reference --readFilesIn {$sample}.test.R2.fastq {$sample}.test.R1.fastq --runDirPerm All_RWX --outSAMtype None \
    --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) \
    --soloUMIlen $UMILEN --soloStrand Forward \
    --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
    --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
    --soloFeatures Gene GeneFull --soloOutFileNames $report_dir/ features.tsv barcodes.tsv matrix.mtx

## the following is needed in case of bad samples: when a low fraction of reads come from mRNA, experiment will look falsely reverse-stranded
report="${sample}_test_strand/GeneFull/Summary.csv"

UNIQFRQ=`grep "Reads Mapped to Genome: Unique," $report | awk -F "," '{print $2}'`
GENEPCT=`grep "Reads Mapped to GeneFull: Unique GeneFull" $report | awk -F "," -v v=$UNIQFRQ '{printf "%d\n",$2*100/v}'`

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

current_date=$(date +%Y-%m-%d)
current_time=$(date +%H:%M:%S)
state="[runtime] (ready)"
message="Successfully detected sequencing strand"

printf "%s %s %s %-${alignment}s\n" "${green}${current_date}" "${current_time}" "${state}" "${blue}${message}${reset}"

# calculate length
max_length=$(echo "$STRAND" | awk '{print length($0)}')
if [ ${#GENEPCT} -gt $max_length ]; then
  max_length=${#GENEPCT}
fi

# print as table
printf "%-25s %-${max_length}s\n" "Paired-end mode:" "$PAIRED"
printf "%-25s %-${max_length}s\n" "Strand (Forward = 3', Reverse = 5'):" "$STRAND"
printf "%-25s %-${max_length}s\n" "%reads same strand as gene:" "$GENEPCT"
printf "%-25s %s\n" "CB whitelist:" "$BC"
printf "%-25s %s\n" "matches out of 200,000:" "$NBC3 (v3), $NBC2 (v2), $NBC1 (v1), $NBCA (multiome)"
printf "%-25s %s\n" "CB length:" "$CBLEN"
printf "%-25s %s\n" "UMI length:" "$UMILEN"
echo "-----------------------------------------------------------------------------"
printf "%-25s %s\n" "Read 1 files:" "$barcode_reads"
echo "-----------------------------------------------------------------------------"
printf "%-25s %s\n" "Read 2 files:" "$cdna_reads"
echo "-----------------------------------------------------------------------------"



current_date=$(date +%Y-%m-%d)
current_time=$(date +%H:%M:%S)
state="[runtime] (run)"
message="Start to run starsolo"

printf "%s %s %s %-${alignment}s\n" "${green}${current_date}" "${current_time}" "${state}" "${blue}${message}${reset}"


if [[ $PAIRED == "True" ]]
then
    ## note the R1/R2 order of input fastq reads and --soloStrand Forward for 5' paired-end experiment
    STAR --runThreadN $threads --genomeDir $reference --readFilesIn $barcode_reads $cdna_reads --runDirPerm All_RWX $file_command \
    --soloBarcodeMate 1 --clip5pNbases 39 0 \
    --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloCBstart 1 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand Forward \
    --soloUMIdedup $soloUMIdedup --soloCBmatchWLtype $soloCBmatchWLtype --soloUMIfiltering $soloUMIfiltering \
    --soloCellFilter $soloCellFilter --outFilterScoreMin $outFilterScoreMin \
    --soloMultiMappers $soloMultiMappers \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattrRGline ID:$sample PL:illumina SM:$sample LB:$CHEMISTRY \
    --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN \
    --outFileNamePrefix ./$sample/\
    $variousParams  
else 
    STAR --runThreadN $threads --genomeDir $reference --readFilesIn $cdna_reads $barcode_reads --runDirPerm All_RWX $file_command \
    --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand $STRAND \
    --soloUMIdedup $soloUMIdedup --soloCBmatchWLtype $soloCBmatchWLtype --soloUMIfiltering $soloUMIfiltering \
    --soloCellFilter $soloCellFilter --outFilterScoreMin $outFilterScoreMin \
    --soloMultiMappers $soloMultiMappers \
    --clipAdapterType CellRanger4 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattrRGline ID:$sample PL:illumina SM:$sample LB:$CHEMISTRY \
    --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN \
    --outFileNamePrefix ./$sample/ \
    $variousParams

fi



current_date=$(date +%Y-%m-%d)
current_time=$(date +%H:%M:%S)
state="[runtime] (ready)"
message="Finished work"

printf "%s %s %s %-${alignment}s\n" "${green}${current_date}" "${current_time}" "${state}" "${blue}${message}${reset}"


rm -rf {$sample}.test.R1.fastq
rm -rf {$sample}.test.R2.fastq
rm -rf $report_dir

sleep 15

cd $sample/

current_date=$(date +%Y-%m-%d)
current_time=$(date +%H:%M:%S)
state="[runtime] (run)"
message="Proprecess output"

printf "%s %s %s %-${alignment}s\n" "${green}${current_date}" "${current_time}" "${state}" "${blue}${message}${reset}"


features_file="${sample}_features.tsv"
matrix_file="${sample}_matrix.mtx"
barcodes_file="${sample}_barcodes.tsv"
mapped_bam_file="Aligned_sortedByCoord_out.bam"
ln -sr "Solo.out/Gene/raw/features.tsv" "$features_file" ;
ln -sr "Solo.out/Gene/raw/matrix.mtx" "$matrix_file" ; 
ln -sr "Solo.out/Gene/raw/barcodes.tsv" "$barcodes_file" ;\
mv "Aligned.sortedByCoord.out.bam" "$mapped_bam_file";\


current_date=$(date +%Y-%m-%d)
current_time=$(date +%H:%M:%S)
state="[runtime] (ready)"
message="Sucessfully proprecessed output"

printf "%s %s %s %-${alignment}s\n" "${green}${current_date}" "${current_time}" "${state}" "${blue}${message}${reset}"
