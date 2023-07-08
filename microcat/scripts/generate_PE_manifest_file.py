import microcat


PLATFORM = None

if config["params"]["host"]["starsolo"]["do"]:
    if "tenX" in config["params"]["host"]["starsolo"]["assay"]:
        PLATFORM = "tenX"
    elif "Smartseq" in config["params"]["host"]["starsolo"]["assay"]:
        PLATFORM = "smartseq"
    elif "Seq-well" in config["params"]["host"]["starsolo"]["assay"]:
        PLATFORM = "tenX"
    else:
        raise ValueError("Platform must be either 'tenX' or 'smartseq'")
elif config["params"]["host"]["cellranger"]["do"]:
    PLATFORM = "tenX"
else:
    raise ValueError("Platform must be either 'tenX' or 'smartseq'")

    
try:
    SAMPLES = microcat.parse_samples(snakemake.input[0],platform = PLATFORM)
except FileNotFoundError:
    warning(f"ERROR: the samples file does not exist. Please see the README file for details. Quitting now.")
    sys.exit(1)


SAMPLES = SAMPLES.reset_index()
# Extract required columns from the parsed samples DataFrame
manifest_df = SAMPLES[['fq1', 'fq2', 'patient']]

# Rename the columns to match the desired format
# manifest_df.columns = ['Read1-file-name', 'Read2-file-name', 'Cell-id']

# Save the manifest DataFrame to a TSV file
manifest_df.to_csv(snakemake.output[0], sep='\t', index=False)