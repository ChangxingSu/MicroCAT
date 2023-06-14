import os
import microcat
import requests
import gzip

MICROCAT_DIR = microcat.__path__[0]
CONFIGS_DIR = os.path.join(MICROCAT_DIR, 'config')
DEFAULT_DATA_DIR = os.path.join(MICROCAT_DIR, 'data')


def download_file(url, destination):
    response = requests.get(url)
    with open(destination, 'wb') as file:
        file.write(response.content)


def extract_gzip_file(gzip_file, destination):
    with gzip.open(gzip_file, 'rb') as f_in:
        with open(destination, 'wb') as f_out:
            f_out.write(f_in.read())


def update_config_file(cell_white_list_path):
    config_file_path = os.path.join(CONFIGS_DIR, 'config.yaml')
    conf = microcat.parse_yaml(config_file_path)
    # Update the cell whitelist location in the config.yaml file here
    # Modify according to your specific requirements
    conf["params"]["data_dir"] = cell_white_list


# Create the data folder
data_folder = os.path.join(MICROCAT_DIR, 'data')
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

# Set the paths and URLs for the files
tenX_v1_barcode_url = 'https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-april-2014_rc.txt'
tenX_v1_barcode_filename = '737K-april-2014_rc.txt'
tenX_v1_barcode_save_path = os.path.join(data_save_path, tenX_v1_barcode_filename)

# Download the file if it doesn't exist
if not os.path.exists(tenX_v1_barcode_save_path):
    download_file(tenX_v1_barcode_url, tenX_v1_barcode_save_path)


tenX_v2_barcode_url = 'https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt'
tenX_v2_barcode_filename = '737K-august-2016.txt'
tenX_v2_barcode_save_path = os.path.join(data_save_path, tenX_v2_barcode_filename)

# Download the file if it doesn't exist
if not os.path.exists(tenX_v2_barcode_save_path):
    download_file(tenX_v2_barcode_url, tenX_v2_barcode_save_path)


tenX_v3_barcode_url = 'https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz'
tenX_v3_barcode_filename = '3M-february-2018.txt'
tenX_v3_barcode_save_path = os.path.join(data_save_path, tenX_v3_barcode_filename)

# Download and extract the file if it doesn't exist
if not os.path.exists(tenX_v3_barcode_save_path):
    tenX_v3_barcode_gz_save_path = tenX_v3_barcode_save_path + '.gz'
    download_file(tenX_v3_barcode_url, tenX_v3_barcode_gz_save_path)
    extract_gzip_file(tenX_v3_barcode_gz_save_path, tenX_v3_barcode_save_path)
    os.remove(tenX_v3_barcode_gz_save_path)

# Update the config.yaml file
update_config_file(data_save_path)