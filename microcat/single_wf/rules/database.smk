

rule download_bacteria_fna_refseq:
    params:
        
    output:
        temp(BACTERIA_FNA_FILE)
    shell:
        "wget -O - {params.url}/bacteria.{wildcards.fn}.1.genomic.fna.gz | "
        "gunz
        

rule zip:

    shell:
    '''
    python /data/project/host-microbiome/kraken_metaphlan_comparison/database_building_scripts/zip.py --directory /data/database/RefseqV218_download/bacteria --output /data/database/RefseqV218_download/bacteria.fna --num_workers 40 --log /data/database/RefseqV218_download/add_bacteria.log
    '''


    python /data/project/host-microbiome/kraken_metaphlan_comparison/database_building_scripts/zip.py --directory /data/database/RefseqV218_download/viral --output /data/database/RefseqV218_download/viral.fna --num_workers 40 --log /data/database/RefseqV218_download/add_viral.log

    python /data/project/host-microbiome/kraken_metaphlan_comparison/database_building_scripts/zip.py --directory /data/database/RefseqV218_download/archaea --output /data/database/RefseqV218_download/archaea.fna --num_workers 40 --log /data/database/RefseqV218_download/add_archaea.log

    python /data/project/host-microbiome/kraken_metaphlan_comparison/database_building_scripts/zip.py --directory /data/database/RefseqV218_download/fungi --output /data/database/RefseqV218_download/fungi.fna --num_workers 40 --log /data/database/RefseqV218_download/add_fungi.log



python /data/project/host-microbiome/kraken_metaphlan_comparison/database_building_scripts/zip.py --directory /data/database/RefseqV218_download/fungi --output /data/database/RefseqV218_download/protozoa.fna --num_workers 40 --log /data/database/RefseqV218_download/add_protozoa.log

python /data/project/host-microbiome/kraken_metaphlan_comparison/database_building_scripts/get_ncbi_other_domains.py --domain protozoa --complete True --folder  /data/database/RefseqV218_download/protozoa --download_genomes True --log_file /data/database/RefseqV218_download/run_protozoa_download.log --processors 40


    samtools faidx microbiome_RefseqV218_Compelete-Chromon.fasta

gunzip -c ./*.fna.gz >>../microbiome_RefseqV218_Compelete-Chromon.fna

    awk '!/^>/ { n=length($0); if (n < l) { l=n } } END { print l }' microbiome_RefseqV218_Compelete-Chromon.fna >test.log


    seqkit seq -w 80 microbiome_RefseqV218_Compelete-Chromon.fna > microbiome_RefseqV218_Compelete-Chromon.fna

for i in $(ls *.fna);
do
  cat ${i}
  # 在每一个fasta输出之后再输出一个空行
  echo
done > .microbiome_RefseqV218_Compelete-Chromon.fna

cat .microbiome_RefseqV218_Compelete-Chromon.fna | grep "NZ_CP109831.1"

rule BuildBwaMemIndexImage:
    input:
        fasta_file = HOST_FASTA,
        fai_file = rules.IndexFasta.output.fai_file,
        gatk4_jar_override = GATK4_JAR_OVERRIDE
    output:
        img_file = "{fasta_file}.img"
    shell:
        """
        set -e
        export GATK_LOCAL_JAR={input.gatk4_jar_override}
        gatk BwaMemIndexImageCreator -I {input.fasta_file}
        """

        gatk BwaMemIndexImageCreator -I /data/database/RefseqV218_download/microbiome_RefseqV218_Compelete-Chromon.fasta --output /data/database/pathseq_database/refv218/microbiome_RefseqV218_Compelete-Chromon.fasta.img

ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/RefSeq-release218.catalog.gz

gatk PathSeqBuildReferenceTaxonomy \
    -R /data/database/RefseqV218_download/microbiome_RefseqV218_Compelete-Chromon.fasta \
    --refseq-catalog /data/database/pathseq_database/refv218/RefSeq-release218.catalog.gz \
    --min-non-virus-contig-length 2000\
    --tax-dump /data/database/pathseq_database/refv218/taxdump.tar.gz \
    -O /data/database/pathseq_database/refv218/microbiome_RefseqV218_Compelete-Chromon.db


ln -s
gatk CreateSequenceDictionary R=/data/database/RefseqV218_download/microbiome_RefseqV218_Compelete-Chromon.fasta O=/data/database/pathseq_database/refv218/microbiome_RefseqV218_Compelete-Chromon.dict