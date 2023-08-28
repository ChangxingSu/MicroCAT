import os
from pathlib import Path
import pysam

class BamWriter:
    def __init__(self, alignment, barcodes, prefix, output_dir):
        self.alignment = alignment
        self.prefix = prefix
        self.barcodes = set(barcodes)
        self.output_dir = output_dir
        self._out_files = {}

    def write_record_to_barcode(self, rec, barcode):
        if barcode not in self.barcodes:
            return
        if barcode not in self._out_files:
            self._open_file_for_barcode(barcode)
        self._out_files[barcode].write(rec)

    def _open_file_for_barcode(self, barcode):
        output_path = os.path.join(self.output_dir, f"{self.prefix}_{barcode}.bam")
        self._out_files[barcode] = pysam.AlignmentFile(
            output_path, "wb", template=self.alignment
        )

def main(input_bam, barcodes_file, output_prefix, contigs, barcode_tag, output_dir):
    alignment = pysam.AlignmentFile(input_bam)
    if Path(barcodes_file).is_file():
        with open(barcodes_file, "r") as fh:
            barcodes = [l.rstrip() for l in fh.readlines()]
    else:
        barcodes = [barcodes_file]
        print(f"Extracting single barcode: {barcodes}")
    writer = BamWriter(alignment=alignment, barcodes=barcodes, prefix=output_prefix, output_dir=output_dir)
    
    # if contigs == ".":
    #     print("Extracting reads from all contigs")
    #     recs = [alignment.fetch()]
    # else:
    #     if "-" in contigs:
    #         start, end = contigs.split("-")
    #         print(f"Extracting reads from contigs {start} to {end}")
    #         recs = (alignment.fetch(str(contig)) for contig in range(int(start), int(end) + 1))
    #     elif "," in contigs:
    #         contigs = contigs.split(",")
    #         print(f"Extracting reads from contigs {contigs}")
    #         recs = (alignment.fetch(str(contig)) for contig in contigs)
    #     else:
    #         print(f"Extracting reads for one contig: {contigs}")
    #         recs = (alignment.fetch(c) for c in [contigs])
    for sread in alignment:
        try:
            barcode = sread.get_tag(barcode_tag)
            print(barcode)
            writer.write_record_to_barcode(rec=sread, barcode=barcode)
        except KeyError:
            pass
            


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Split a 10x barcoded sequencing file into barcode-specific BAMs")
    parser.add_argument("--input_bam", help="Input BAM file")
    parser.add_argument("--barcodes_file", help="File containing barcodes, or a single barcode")
    parser.add_argument("--output_prefix", help="Output file prefix for barcode-specific BAM files")
    parser.add_argument("--output_dir", help="Output directory for saving barcode-specific BAM files")
    parser.add_argument("--contigs", help="Contigs to extract reads from: '.' for all contigs, 'chr1' for the contig 'chr1', '1-5' for chromosomes 1, 2, 3, 4, and 5, or '1,2,3' for contigs 1, 2, and 3")
    parser.add_argument("--barcode_tag", default="CB", help="Barcode tag to use for extracting barcodes")

    args = parser.parse_args()

    main(
        input_bam=args.input_bam,
        barcodes_file=args.barcodes_file,
        output_prefix=args.output_prefix,
        output_dir=args.output_dir,
        contigs=args.contigs,
        barcode_tag=args.barcode_tag,
    )
