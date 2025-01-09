import pandas as pd
import warnings
import os
import json
import subprocess
import shutil


try:
    samples_df = pd.read_csv(snakemake.input[0], sep="\t")
except FileNotFoundError:
    warnings.warn(f"ERROR: the samples file does not exist. Please see the README file for details. Quitting now.")
    sys.exit(1)


# 创建一个空字典以存储处理信息
processing_info = []

if not set(['id', 'fq1', 'fq2']).issubset(samples_df.columns):
    raise ValueError("Columns 'id', 'fq1', 'fq2' must exist in the sample.tsv")

# Extract library, lane, and plate from id
samples_df[['patient_tissue_lane_plate', 'library']] = samples_df['id'].str.rsplit("_", n=1, expand=True)
samples_df['is_lane'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1].startswith("L"))
samples_df.loc[samples_df['is_lane'], 'lane'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1])
samples_df['patient_tissue'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: '_'.join(x.split('_')[:-1]))
samples_df = samples_df.drop(columns=['patient_tissue_lane_plate'])
# samples_df = samples_df.loc[(samples_df.plate == snakemake.wildcards["plate"]) & (samples_df["patient_tissue"] == snakemake.wildcards["sample"])]
samples_df = samples_df.loc[(samples_df["patient_tissue"] == snakemake.wildcards["sample"])]


samples_df = samples_df.reset_index()

# # Extract required columns from the parsed samples DataFrame
# manifest_df = samples_df[['fq1', 'fq2', 'cell']]

# # Filter out rows where fq2 is NaN
# manifest_df = manifest_df[manifest_df['fq2'].notna()]
output = str(snakemake.output[0])
outdir = os.path.dirname(output)
# temp file
# execute(f'''rm -rf {outdir}''')
# execute(f'''mkdir -p {outdir}''')
# 如果目录存在，先删除
if os.path.exists(outdir):
    shutil.rmtree(outdir)

# 创建目录
os.makedirs(outdir)
# 遍历样本数据，根据patient_tissue重新命名文件并创建软链接
for index, sample in samples_df.iterrows():
    fq1 = sample['fq1']
    fq2 = sample['fq2']
    PATIENT_TISSUE = sample['patient_tissue']
    LANE = sample['lane']
    LIBRARY = sample['library']


    # 检查文件后缀是否为fastq
    if fq1.endswith(".fastq") and fq2.endswith(".fastq"):
        TYPE = 'fastq'
        # 构建新的文件名
        new_fq1 = os.path.join(outdir, f'{PATIENT_TISSUE}_S1_{LANE}_R1_{LIBRARY}.fastq')
        new_fq2 = os.path.join(outdir, f'{PATIENT_TISSUE}_S1_{LANE}_R2_{LIBRARY}.fastq')


        # 创建软链接
        # execute(f'''ln -s {fq1} {new_fq1} ''')
        # execute(f'''ln -s {fq2} {new_fq2} ''')
        subprocess.call(['ln', '-s', fq1, new_fq1])
        subprocess.call(['ln', '-s', fq2, new_fq2])
        # 构建处理信息字典
        processing_info.append({
            'file_type': TYPE,
            'Lane': LANE,
            'Library': LIBRARY,
            'Original_fq1': fq1,
            'New_fq1': new_fq1,
            'Original_fq2': fq2,
            'New_fq2': new_fq2
        })
    elif fq1.endswith(".fastq.gz") and fq2.endswith(".fastq.gz"):
        TYPE = 'fastq.gz'
        # 构建新的文件名
        new_fq1 = os.path.join(outdir, f'{PATIENT_TISSUE}_S1_{LANE}_R1_{LIBRARY}.fastq.gz')
        new_fq2 = os.path.join(outdir, f'{PATIENT_TISSUE}_S1_{LANE}_R2_{LIBRARY}.fastq.gz')

        # 创建软链接
        subprocess.call(['ln', '-s', fq1, new_fq1])
        subprocess.call(['ln', '-s', fq2, new_fq2])

        # 构建处理信息字典
        processing_info.append({
            'file_type': TYPE,
            'Lane': LANE,
            'Library': LIBRARY,
            'Original_fq1': fq1,
            'New_fq1': new_fq1,
            'Original_fq2': fq2,
            'New_fq2': new_fq2
        })
    else:
        print(f"Skipping files {fq1} and {fq2} with incorrect extensions.")

# 将处理信息写入JSON文件
with open(snakemake.output[0], 'w') as json_file:
    json.dump(processing_info, json_file, indent=4)