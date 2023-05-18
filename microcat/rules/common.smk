import pandas

class ansitxt:
    RED = '\033[31m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def warning(msg):
    print(f"\n{ansitxt.BOLD}{ansitxt.RED}{msg}{ansitxt.ENDC}\n",file=sys.stderr)


def parse_samples(samples_tsv):
    # Load the sample.tsv file into a pandas DataFrame
    df = pandas.read_csv(samples_tsv, sep='\t')
    
    # Split the 'id' column into 'sample' and 'lane' columns at the last underscore '_'
    # The resulting 'sample' and 'lane' columns are added to the DataFrame
    df[['sample', 'SX', 'lane']] = df['id'].str.rsplit(pat='_', n=2, expand=True)

    # Create a 'fastqs_dir' column that contains the directory of the fastq files
    # This is done by removing the filename from the 'fq1' column
    df['fastqs_dir'] = df['fq1'].apply(lambda x: '/'.join(x.split('/')[:-1]))
    
    # If the 'lane' column contains any NaN values (due to the absence of a lane in 'id'), replace them with 'no_lane'
    df['lane'] = df['lane'].fillna('no_lane')  
    
    # Return a list of dictionaries where each dictionary contains the 'sample', 'lane', and 'fastqs_dir' for a sample
    return df[['sample', 'lane', 'fastqs_dir']].to_dict(orient='records')

def get_fastqs_dir(wildcards):

    # Iterate over each sample in the global '_samples' list
    for s in _samples:
        # If the 'sample' field of the current sample dictionary matches the 'sample' wildcard
        if s['sample'] == sample_name:
            # Return the 'fastqs_dir' field of the matching sample
            return s['fastqs_dir']

    # If no matching sample is found in the '_samples' list, raise a ValueError with an appropriate message
    raise ValueError(f"No fastqs_dir found for sample: {sample_name}")

    
def get_starsolo_sample_id(samples, wildcards, read):
    sample_reads = [s[read] for s in samples if s['sample'] == wildcards.sample]
    if read == "fq1":
        return ','.join(sorted(sample_reads))
    elif read == "fq2":
        return ' '.join(sorted(sample_reads))