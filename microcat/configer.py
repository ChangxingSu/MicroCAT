#!/usr/bin/env python

import argparse
import os
import sys
import shutil  #shutil provides a higher-level interface for file operations
from ruamel.yaml import YAML

# parse yaml style
def parse_yaml(yaml_file):
    """
    Read and parse a YAML file.

    :param yaml_file: Path to the YAML file
    :return: Parsed YAML data
    """
    yaml = YAML()
    with open(yaml_file, "r") as f:
        return yaml.load(f)


def update_config(yaml_file_old, yaml_file_new, yaml_content, remove=True):
    """
    Update the configuration file.

    :param yaml_file_old: Path to the old YAML file
    :param yaml_file_new: Path to the new YAML file
    :param yaml_content: YAML content to be written
    :param remove: Whether to remove the old YAML file
    """
    yaml = YAML()
    yaml.default_flow_style = False
    if remove:
        os.remove(yaml_file_old)
    with open(yaml_file_new, "w") as f:
        yaml.dump(yaml_content, f)


class MicrocatConfig:
    """
    config project directory
    """
    sub_dirs = [
        "envs",
        "profiles",
        "results",
        "figures",
        "notebooks",
        "logs"
    ]

    def __init__(self, work_dir,config_type):
        if config_type == "single":
            config_subdir = "single_wf"
        elif config_type == "bulk":
            config_subdir = "bulk_wf"
        elif config_type == "spatial":
            config_subdir = "spatial_wf"
        elif config_type == "multi":
            config_subdir = "multi_wf"
        else:
            raise ValueError("Invalid config_type")
        self.work_dir = os.path.realpath(work_dir)
        self.microcat_dir = os.path.dirname(os.path.abspath(__file__))
        self.config_file = os.path.join(self.microcat_dir,config_subdir,"config", "config.yaml")
        self.envs_dir = os.path.join(self.microcat_dir,config_subdir, "envs")
        self.profiles_dir = os.path.join(self.microcat_dir,"profiles")
        self.new_config_file = os.path.join(self.work_dir, "config.yaml")

    def __str__(self):
        message = """
            ███╗   ███╗██╗ ██████╗██████╗  ██████╗  ██████╗ █████╗ ████████╗
            ████╗ ████║██║██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔══██╗╚══██╔══╝
            ██╔████╔██║██║██║     ██████╔╝██║   ██║██║     ███████║   ██║   
            ██║╚██╔╝██║██║██║     ██╔══██╗██║   ██║██║     ██╔══██║   ██║   
            ██║ ╚═╝ ██║██║╚██████╗██║  ██║╚██████╔╝╚██████╗██║  ██║   ██║   
            ╚═╝     ╚═╝╚═╝ ╚═════╝╚═╝  ╚═╝ ╚═════╝  ╚═════╝╚═╝  ╚═╝   ╚═╝   
    Microbiome Identification in Cell Resolution from Omics-Computational Analysis Toolbox\n
Thanks for using microcat.
A microbiome identification project has been created at %s
if you want to create fresh conda environments:
        microcat run-local --conda-create-envs-only
        microcat run-remote --conda-create-envs-only
if you have environments:
        microcat run-local --help
        microcat run-remote --help
""" % (
            self.work_dir
        )

        return message

    def create_dirs(self):
        """
        create project directory
        """
        # If the project path does not exist, create a project directory
        if not os.path.exists(self.work_dir):
            os.mkdir(self.work_dir)
        
        # Iterate through the list of subdirectories, creating a path for each subdirectory
        for sub_dir in MicrocatConfig.sub_dirs:
            os.makedirs(os.path.join(self.work_dir, sub_dir), exist_ok=True)
        
        # Iterate through the files in the envs directory and copy them into the envs subdirectory of the project directory
        overwrite_all = None
        for i in os.listdir(self.envs_dir):
            dest_file = os.path.join(self.work_dir, "envs", i)
            # If the target file already exists, output a prompt message
            if os.path.exists(dest_file):
                if overwrite_all is None:
                    print(f"Warning: The file '{dest_file}' already exists.")
                    proceed = input("Do you want to overwrite it? (y/n/all/quit): ").lower()
                    if proceed == 'n':
                        print("Skip updating this file.")
                        continue
                    elif proceed == 'all':
                        overwrite_all = True
                    elif proceed == 'quit':
                        print("Aborted.")
                        sys.exit(1)
                elif overwrite_all is False:
                    print("Skip updating this file.")
                    continue
            # Copy the source file to the target file
            shutil.copyfile(os.path.join(self.envs_dir, i), dest_file)
        
        if os.path.exists(self.profiles_dir):
            overwrite_all = None
            # Iterate through the subdirectories of the profiles directory and copy them into the profiles subdirectory of the project directory
            for i in os.listdir(self.profiles_dir):
                dest_dir = os.path.join(self.work_dir, "profiles", i)

                # If the target directory already exists, output a prompt message and terminate the program
                if os.path.exists(dest_dir):
                    if overwrite_all is None:
                        print(f"{dest_dir} already exists.")
                        proceed = input("Do you want to overwrite it? (y/n/all/quit): ").lower()
                        if proceed == 'n':
                            print("Skip updating this file.")
                            continue
                        elif proceed == 'all':
                            overwrite_all = True
                        elif proceed == 'quit':
                            print("Aborted.")
                            sys.exit(1)
                    elif overwrite_all is False:
                        print("Skip updating this dir.")
                        continue
                # Otherwise copy the source directory to the target directory
                else:
                    shutil.copytree(os.path.join(self.profiles_dir, i), dest_dir)
        else:
            # Handle the case when 'profiles' directory does not exist
            print("\033[93mWARNING: The 'profiles' directory does not exist in the original folder.")
            print("You can only run the program locally using the run-local option.")
            print("If you want to run it on a cluster, please refer to the documentation for the 'cluster run' section.\033[0m")
            print("Skipping the profiles directory.")

    def get_config(self):
        """
        Reads the default configuration file config.yaml and returns the parsed configuration
        """
        return parse_yaml(self.config_file)


# https://github.com/Ecogenomics/CheckM/blob/master/checkm/customHelpFormatter.py
class custom_help_formatter(argparse.HelpFormatter):
    """Provide a customized format for help output.
    http://stackoverflow.com/questions/9642692/argparse-help-without-duplicate-allcaps
    """

    def _split_lines(self, text, width):
        return text.splitlines()

    def _get_help_string(self, action):
        h = action.help
        if "%(default)" not in action.help:
            if (
                action.default != ""
                and action.default != []
                and action.default != None
                and action.default != False
            ):
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [
                        argparse.OPTIONAL, argparse.ZERO_OR_MORE]

                    if action.option_strings or action.nargs in defaulting_nargs:
                        if "\n" in h:
                            lines = h.splitlines()
                            lines[0] += " (default: %(default)s)"
                            h = "\n".join(lines)
                        else:
                            h += " (default: %(default)s)"
            return h

    def _fill_text(self, text, width, indent):
        return "".join([indent + line for line in text.splitlines(True)])

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            (metavar,) = self._metavar_formatter(action, default)(1)
            return metavar

        else:
            parts = []

            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append(option_string)

                return "%s %s" % (", ".join(parts), args_string)

            return ", ".join(parts)

    def _get_default_metavar_for_optional(self, action):
        return action.dest.upper()

    def _get_default_metavar_for_positional(self, action):
        return action.dest