

rule download_bacteria_fna_refseq:
    params:
        
    output:
        temp(BACTERIA_FNA_FILE)
    shell:
        "wget -O - {params.url}/bacteria.{wildcards.fn}.1.genomic.fna.gz | "
        "gunzip -c > {output}"