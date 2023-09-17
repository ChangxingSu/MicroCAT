import os
import microcat
import subprocess
import gzip
from ruamel.yaml import YAML

MICROCAT_DIR = microcat.__path__[0]
CONFIGS_DIR = os.path.join(MICROCAT_DIR,'single_wf','config')
DEFAULT_DATA_DIR = os.path.join(MICROCAT_DIR,'single_wf', 'data')


def download_file(url, destination):
    command = f"wget -O {destination} {url}"
    subprocess.call(command, shell=True)


def extract_gzip_file(gzip_file, destination):
    with gzip.open(gzip_file, 'rb') as f_in:
        with open(destination, 'wb') as f_out:
            f_out.write(f_in.read())


def update_config_file(cell_white_list_path):
    config_file_path = os.path.join(CONFIGS_DIR, 'config.yaml')
    conf = microcat.parse_yaml(config_file_path)
    # Update the cell whitelist location in the config.yaml file here
    # Modify according to your specific requirements
    conf["datas"]["barcode_list_dirs"]["tenX"] = cell_white_list_path
    # Update the configuration file
    yaml = YAML()
    with open(config_file_path, 'w') as f_out:
        yaml.dump(conf, f_out)


# Create the data folder
data_folder = os.path.join(MICROCAT_DIR,'single_wf','data')
if not os.path.exists(data_folder):
    os.makedirs(data_folder)


print("Please choose an option:")
print("1. Download files to the default location")
print("2. Download files to a custom location and update config.yaml")
print("3. Input the cell whitelist location directly and update config.yaml")
print("4. Abort the script")

choice = input("Enter the number corresponding to your choice: ")

if choice == '1':
    # Download files to the default location
    data_save_path = DEFAULT_DATA_DIR

elif choice == '2':
    # Download files to a custom location
    custom_data_dir = input("Enter the custom download path: ")
    if not os.path.exists(custom_data_dir):
        os.makedirs(custom_data_dir)
    data_save_path = custom_data_dir

elif choice == '3':
    # Input the cell whitelist location directly
    cell_white_list_path = input("Enter the cell whitelist location: ")
    update_config_file(cell_white_list_path)
    exit()

elif choice == '4':
    # Abort the script
    exit()

else:
    print("Invalid choice!")
    exit()

tenX_v1_barcode_url = 'https://teichlab.github.io/scg_lib_structs/data/737K-april-2014_rc.txt.gz'
tenX_v1_barcode_filename = '737K-april-2014_rc.txt'
tenX_v1_barcode_save_path = os.path.join(data_save_path, tenX_v1_barcode_filename)


# Download the file if it doesn't exist
if not os.path.exists(tenX_v1_barcode_save_path):
    tenX_v1_barcode_gz_save_path = tenX_v1_barcode_save_path + '.gz'
    print(f"Downloading {tenX_v1_barcode_filename}.gz")
    download_file(tenX_v1_barcode_url, tenX_v1_barcode_gz_save_path)
    print(f"Extracting {tenX_v1_barcode_filename}.gz")
    extract_gzip_file(tenX_v1_barcode_gz_save_path, tenX_v1_barcode_save_path)
    os.remove(tenX_v1_barcode_gz_save_path)

tenX_v2_barcode_url = 'https://teichlab.github.io/scg_lib_structs/data/737K-august-2016.txt.gz'
tenX_v2_barcode_filename = '737K-august-2016.txt'
tenX_v2_barcode_save_path = os.path.join(data_save_path, tenX_v2_barcode_filename)

# Download the file if it doesn't exist
if not os.path.exists(tenX_v2_barcode_save_path):
    tenX_v2_barcode_gz_save_path = tenX_v2_barcode_save_path + '.gz'
    print(f"Downloading {tenX_v2_barcode_filename}.gz")
    download_file(tenX_v2_barcode_url, tenX_v2_barcode_gz_save_path)
    print(f"Extracting {tenX_v2_barcode_filename}.gz")
    extract_gzip_file(tenX_v2_barcode_gz_save_path, tenX_v2_barcode_save_path)
    os.remove(tenX_v2_barcode_gz_save_path)

tenX_v3_barcode_url = 'https://teichlab.github.io/scg_lib_structs/data/3M-february-2018.txt.gz'
tenX_v3_barcode_filename = '3M-february-2018.txt'
tenX_v3_barcode_save_path = os.path.join(data_save_path, tenX_v3_barcode_filename)

# Download and extract the file if it doesn't exist
if not os.path.exists(tenX_v3_barcode_save_path):
    tenX_v3_barcode_gz_save_path = tenX_v3_barcode_save_path + '.gz'
    print(f"Downloading {tenX_v3_barcode_filename}.gz")
    download_file(tenX_v3_barcode_url, tenX_v3_barcode_gz_save_path)
    print(f"Extracting {tenX_v3_barcode_filename}.gz")
    extract_gzip_file(tenX_v3_barcode_gz_save_path, tenX_v3_barcode_save_path)
    os.remove(tenX_v3_barcode_gz_save_path)


tenX_multiome_barcode_url = 'https://woldlab.caltech.edu/~diane/genome/737K-arc-v1.txt.gz'
tenX_multiome_barcode_filename = '737K-arc-v1.txt'
tenX_multiome_barcode_save_path = os.path.join(data_save_path, tenX_multiome_barcode_filename)


# Download the file if it doesn't exist
if not os.path.exists(tenX_multiome_barcode_save_path):
    tenX_multiome_barcode_gz_save_path = tenX_multiome_barcode_save_path + '.gz'
    print(f"Downloading {tenX_multiome_barcode_filename}.gz")
    download_file(tenX_multiome_barcode_url, tenX_multiome_barcode_gz_save_path)
    print(f"Extracting {tenX_multiome_barcode_filename}.gz")
    extract_gzip_file(tenX_multiome_barcode_gz_save_path, tenX_multiome_barcode_save_path)
    os.remove(tenX_multiome_barcode_gz_save_path)
    
# Update the config.yaml file
update_config_file(data_save_path)