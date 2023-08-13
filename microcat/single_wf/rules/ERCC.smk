ERCC92_URL = "https://tools.thermofisher.com/content/sfs/manuals/ERCC92.zip"
ERCC92_FA = os.path.join("data", "ERCC92.fa")
ERCC92_GTF = os.path.join("data", "ERCC92.gtf")

localrules: download_ERCC92_spike_ins, combine_human_ERCC92_GTF, combine_human_ERCC92_fasta

rule download_ERCC92_spike_ins:
    params:
        url=ERCC92_URL,
        odir="data"
    output:
        zip=temp("ERCC92.zip"),
        gtf=ERCC92_GTF,
        fa=ERCC92_FA
    shell:
        "wget {params.url} && unzip {output.zip} -d {params.odir}"

        
rule combine_human_ERCC92_fasta:
    input:
        human_genome=config["human_ref"]["genome"],
        ercc92_fa=ERCC92_FA
    output:
        combined_genome=config["ref"]["genome"]
    shell:
        "cat {input.human_genome} {input.ercc92_fa} > {output.combined_genome}"

rule combine_human_ERCC92_GTF:
    input:
        human_annotation=config["human_ref"]["annotation"],
        ercc92_gtf=ERCC92_GTF
    output:
        combined_annotation=config["ref"]["annotation"]
    shell:
        "cat {input.human_annotation} {input.ercc92_gtf} > {output.combined_annotation}"

