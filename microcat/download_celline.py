import os
import microcat
import subprocess
import gzip
from ruamel.yaml import YAML

MICROCAT_DIR = microcat.__path__[0]
CONFIGS_DIR = os.path.join(MICROCAT_DIR,'single_wf','config')
DEFAULT_DATA_DIR = os.path.join(MICROCAT_DIR,'single_wf', 'data')



# Create the data folder
if not os.path.exists(DEFAULT_DATA_DIR):
    os.makedirs(DEFAULT_DATA_DIR)

cell_line_save_path = os.path.join(DEFAULT_DATA_DIR, "sahmi_cell_line.xlsx")
url = "https://github.com/sjdlabgroup/SAHMI/raw/main/Table%20S4.xlsx"
command = f"wget -O {cell_line_save_path} {url}"

subprocess.call(command, shell=True)
print("Updating config")
config_file_path = os.path.join(CONFIGS_DIR, 'config.yaml')
conf = microcat.parse_yaml(config_file_path)
# Update the cell whitelist location in the config.yaml file here
# Modify according to your specific requirements
config["params"]["classifier"]["krak_study_denosing"]["cell_line"] = cell_line_save_path
# Update the configuration file
yaml = YAML()
with open(config_file_path, 'w') as f_out:
    yaml.dump(conf, f_out)

print("Finishing updating config")