import pysam
import argparse
import logging
import os
from tqdm import tqdm

# Parse command line arguments
parser = argparse.ArgumentParser(description="Split BAM file by tag")
parser.add_argument('--tag', choices=['CB', 'RG'], required=True, help='Tag to use for splitting')
parser.add_argument('--bam_path', required=True, help='Path to input BAM file')
parser.add_argument('--output_dir', required=True, help='Path to output directory')
parser.add_argument('--log_file', default='logfile.log', help='Path to log file')
args = parser.parse_args()

# Setup logging
logging.basicConfig(filename=args.log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


# Log which tag is being processed
logging.info('Processing tag: %s', args.tag)


# Create output directory if it doesn't exist
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# Variable to hold tag index
tag_hold = 'unset'
# Read in unsplit file and loop reads by line
samfile = pysam.AlignmentFile(args.bam_path, "rb")

logging.info('Starting processing of file %s', args.bam_path)

# Wrap the iterator with tqdm to show progress
for read in tqdm(samfile, unit="read"):
    # Tag itr for current read
    try:
        tag_itr = read.get_tag(args.tag)
    except KeyError:
        # if not tag, exclude it 
        continue
        
    # If change in tag or first line; open new file  
    if tag_itr != tag_hold:
        # Close previous split file, only if not first read in file
        if tag_hold != 'unset':
            split_file.close()
        tag_hold = tag_itr

        output_subdir = os.path.join(args.output_dir, "{}".format(tag_itr))
        if not os.path.exists(output_subdir):
            os.makedirs(output_subdir,exist_ok=True)
        
        split_file = pysam.AlignmentFile(os.path.join(output_subdir, "Aligned_sortedByName_unmapped_out.bam".format(args.tag, tag_itr)), "wb", template=samfile)
        # split_file = pysam.AlignmentFile(os.path.join(args.output_dir, "{}_{}.bam".format(args.tag, tag_itr)), "wb", template=samfile)
        logging.info('Opened new file for writing: %s', split_file.filename.decode())

    # Write read with same tag to file
    split_file.write(read)

# Close the last split file and the input file
split_file.close()
samfile.close()

logging.info('Finished processing of file %s', args.bam_path)
