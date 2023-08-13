import os
import microcat
import requests
import gzip

MICROCAT_DIR = microcat.__path__[0]


def download_file(url, destination):
    response = requests.get(url)
    with open(destination, 'wb') as file:
        file.write(response.content)


def extract_gzip_file(gzip_file, destination):
    with gzip.open(gzip_file, 'rb') as f_in:
        with open(destination, 'wb') as f_out:
            f_out.write(f_in.read())


# 创建 data 文件夹
data_folder = os.path.join(MICROCAT_DIR, 'data')
if not os.path.exists(data_folder):
    os.makedirs(data_folder)


# 设置 tenX_v1 文件路径和下载链接
tenX_v1_barcode_url = 'https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-april-2014_rc.txt'
tenX_v1_barcode_filename = '737K-april-2014_rc.txt'
tenX_v1_barcode_save_path = os.path.join(data_folder, tenX_v1_barcode_filename)

# 检查文件是否已经存在，如果不存在则下载
if not os.path.exists(tenX_v1_barcode_save_path):
    download_file(tenX_v1_barcode_url, tenX_v1_barcode_save_path)


# 设置 tenX_v2 文件路径和下载链接
tenX_v2_barcode_url = 'https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt'
tenX_v2_barcode_filename = '737K-august-2016.txt'

tenX_v2_barcode_save_path = os.path.join(data_folder, tenX_v2_barcode_filename)

# 检查文件是否已经存在，如果不存在则下载
if not os.path.exists(tenX_v2_barcode_save_path):
    download_file(tenX_v2_barcode_url, tenX_v2_barcode_save_path)


# 设置 tenX_v3 文件路径和下载链接
tenX_v3_barcode_url = 'https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz'
tenX_v3_barcode_filename = '3M-february-2018.txt'

tenX_v3_barcode_save_path = os.path.join(data_folder, tenX_v3_barcode_filename)

# 检查文件是否已经存在，如果不存在则下载并解压缩
if not os.path.exists(tenX_v3_barcode_save_path):
    tenX_v3_barcode_gz_save_path = tenX_v3_barcode_save_path + '.gz'
    download_file(tenX_v3_barcode_url, tenX_v3_barcode_gz_save_path)
    extract_gzip_file(tenX_v3_barcode_gz_save_path, tenX_v3_barcode_save_path)
    os.remove(tenX_v3_barcode_gz_save_path)
