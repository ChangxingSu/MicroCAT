rule saw_count:
    input:
        image = lambda wildcards: microcat.get_img(SAMPLES, wildcards),
        chip_mask = config["params"]["host"]["saw"]["chip_mask"],
        reference = config["params"]["host"]["saw"]["reference"],
    params:
        variousParams = config["params"]["host"]["saw"]["variousParams"],
        SampleID="{sample}",
        reference = config["params"]["host"]["saw"]["reference"],
        fastqs_dir = os.path.abspath(os.path.join(config["output"]["raw"], "reads/{sample}/")),
        output_dir = os.path.abspath(os.path.join(config["output"]["host"], "spaceranger_count/{sample}")),
        kit_version = config["params"]["host"]["saw"]["kit_version"],
        sequencing_type = config["params"]["host"]["saw"]["sequencing_type"],
    shell:
        """
        saw count \    
        --id={params.SampleID} \
        --kit-version={params.kit_version} \
        --sequencing-type={params.sequencing_type} \
        --chip-mask={input.chip_mask} \
        --organism={params.organism} \
        --tissue={params.tissue} \
        --fastqs={params.fastqs} \
        --reference={input.reference} \
        --image={input.image}
        """