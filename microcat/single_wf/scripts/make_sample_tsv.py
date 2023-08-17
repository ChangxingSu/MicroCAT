import re
import os
import csv

fastq_dir = '/data/scRNA_data/Saliba2016'
fastq_files = os.listdir(fastq_dir)

prefixes = set()
for f in fastq_files:
    m = re.match(r'(.+)_R(\d)_001.fastq.gz', f)
    if m:
        prefixes.add(m.group(1))

rows = []
rows.append(['id', 'fq1', 'fq2'])

for prefix in prefixes:
    r1 = f'{prefix}_R1_001.fastq.gz'
    r2 = f'{prefix}_R2_001.fastq.gz'

    fq1_abs_path = os.path.abspath(os.path.join(fastq_dir, r1))
    fq2_abs_path = os.path.abspath(os.path.join(fastq_dir, r2))

    if r1 in fastq_files and r2 in fastq_files:
        rows.append([prefix, fq1_abs_path, fq2_abs_path])

with open('/home/microcat-sucx/project/scRNA-analysis/benchmark/Saliba2016/sample.tsv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(rows)
