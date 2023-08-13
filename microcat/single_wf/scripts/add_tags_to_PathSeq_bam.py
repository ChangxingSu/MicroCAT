import pysam
import argparse
import logging

# 解析命令行参数
parser = argparse.ArgumentParser(description="Script to process BAM files")
parser.add_argument("--classified_bam",dest='classified_bam', required=True,
        help='BAM file containing the tax ID')
parser.add_argument("--host_bam",dest='host_bam', required=True,
        help='BAM file containing the scRNA Tag')
parser.add_argument("--output",dest='output', required=True,
        help='Output BAM file containing the scRNA Tag and tax ID')
parser.add_argument("--log", help="Path to log file", default="script.log")
args = parser.parse_args()

# 配置日志记录
logging.basicConfig(filename=args.log, level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

# 记录脚本开始
logging.info("Start to add tags to PathSeq classified bam")

# Load and index the CellRanger BAM file
with pysam.AlignmentFile(args.host_bam, mode="rb") as cr_bam:
    cr_idx = pysam.IndexedReads(cr_bam)
    cr_idx.build()
    logging.info("CellRanger BAM file loaded and indexed")

    # Load and iterate through the PathSeq BAM file
    with pysam.AlignmentFile(args.classified_bam, mode="rb") as pathseq_bam:
        output = []
        for seg in pathseq_bam.fetch(until_eof=True):
            cr_list = list(cr_idx.find(seg.query_name))
            if cr_list and cr_list[0].has_tag("CB") and cr_list[0].has_tag("UB"):
                CB = cr_list[0].get_tag(tag="CB")
                UB = cr_list[0].get_tag(tag="UB")
                seg.set_tag("CB", CB, "Z")
                seg.set_tag("UB", UB, "Z")
            output.append(seg)
        logging.info("PathSeq BAM file processed")
        # Write all PathSeq alignments to the output BAM file
        with pysam.AlignmentFile(args.output, mode="wb", template=pathseq_bam) as all_pathseq_bam:
            for seg in output:
                all_pathseq_bam.write(seg)
            logging.info("Output BAM file created")


# 记录脚本结束
logging.info("Script finished")
