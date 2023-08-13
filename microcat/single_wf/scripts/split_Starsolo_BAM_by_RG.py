import pysam
from collections import defaultdict
import os
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description="Split BAM file by RG")
parser.add_argument('--tag', required=True, help='Tag to use for splitting')
parser.add_argument('--bam_path', required=True, help='Path to input BAM file')
parser.add_argument('--output_bam', required=True, help='Path to output directory')
# parser.add_argument('--log_file', default='logfile.log', help='Path to log file')
args = parser.parse_args()


# suppress incorrect error warning - https://github.com/pysam-developers/pysam/issues/939
save = pysam.set_verbosity(0)
# load and iterate through the PathSeq BAM file
starsolo_bam = pysam.AlignmentFile(args.bam_path, mode="rb")
# set verbosity back to original setting
pysam.set_verbosity(save)

output = []
UMI_dict = defaultdict(list)
# seg is an AlignedSegment object
for seg in starsolo_bam.fetch(until_eof=True):
    # not all records will have the CB tag and the UB tag - they should now
    if seg.has_tag("RG"):
        if (seg.get_tag(tag="RG") == args.tag):
            UMI_dict.append(seg)

output_subdir = os.path.join(args.output_dir, "{}".format(tag))
if not os.path.exists(output_subdir):
    os.makedirs(output_subdir,exist_ok=True)
        
# split_file = pysam.AlignmentFile(os.path.join(output_subdir, "Aligned_sortedByName_unmapped_out.bam".format(args.tag, tag_itr)), "wb", template=samfile)

read_group_bam = pysam.AlignmentFile(args.output_bam, mode="wb", template=starsolo_bam)
for UMI in UMI_dict:
    #print(UMI)
    # keep one read per UMI - the read with the highest mapping quality
    # UMI_reads = UMI_dict[UMI]
    # UMI_read = UMI_reads[0]
    # #print(UMI_read)
    # for read in UMI_reads:
    #     #print(read)
    #     if read.mapping_quality > UMI_read.mapping_quality:
    #         UMI_read = read
    read_group_bam.write(UMI)
read_group_bam.close()

starsolo_bam.close()


# suppress incorrect error warning - https://github.com/pysam-developers/pysam/issues/939
save = pysam.set_verbosity(0)
# load and iterate through the PathSeq BAM file
STAR_BAM = pysam.AlignmentFile(snakemake.input[0], mode="rb")
# set verbosity back to original setting
pysam.set_verbosity(save)

output = []
UMI_dict = defaultdict(list)
# seg is an AlignedSegment object
for seg in STAR_BAM.fetch(until_eof=True):
    # not all records will have the CB tag and the UB tag - they should now
    if seg.has_tag("CB") and seg.has_tag("UB"):
        if (seg.get_tag(tag="CB") == snakemake.wildcards["cell"]):
            UMI_dict[seg.get_tag(tag="UB")].append(seg)


CELL_BAM = pysam.AlignmentFile(snakemake.output[0], mode="wb", template=STAR_BAM)
for UMI in UMI_dict:
    #print(UMI)
    # keep one read per UMI - the read with the highest mapping quality
    UMI_reads = UMI_dict[UMI]
    UMI_read = UMI_reads[0]
    #print(UMI_read)
    for read in UMI_reads:
        #print(read)
        if read.mapping_quality > UMI_read.mapping_quality:
            UMI_read = read
    CELL_BAM.write(UMI_read)
CELL_BAM.close()

STAR_BAM.close()