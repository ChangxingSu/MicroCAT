#!/usr/bin/env python

import click
import os
import microcat
import os
import subprocess
import sys
import textwrap
from io import StringIO 
import pandas as pd
import psutil
import json
import re
import sys
import csv
import yaml
from click_help_colors import HelpColorsGroup, HelpColorsCommand
from microcat import MicrocatConfig
from microcat import update_config
from microcat import parse_yaml


DEFAULT_CONFIG = "./config.yaml"

WORKFLOWS_SCRNA = [
    "host_all",
    "kraken2uniq_classified_all",
    "krakenuniq_classified_all",
    "pathseq_classified_all",
    "metaphlan_classified_all",
    "classifier_all",
    "all"]

WORKFLOWS_BULK = [
    "host_all",
    "kraken2uniq_classified_all",
    "classifier_all",
    "all"]

WORKFLOWS_SPATIAL = [
    "host_all",
    "kraken2uniq_classified_all",
    "classifier_all",
    "all"]

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'],ignore_unknown_options=True,allow_extra_args=True)

def read_config(config_path):
    # Return the appropriate configuration dictionary based on the workflow value
    with open(config_path, "r") as file:
        config_data = yaml.safe_load(file)

    workflow = config_data.get("workflow")
    if workflow == "single":
        return WORKFLOWS_SCRNA
    elif workflow == "bulk":
        return WORKFLOWS_BULK
    elif workflow == "spatial":
        return WORKFLOWS_SPATIAL
    elif workflow == "multi":
        return WORKFLOWS_MULTI
    else:
        raise ValueError("Invalid workflow value in config file.")

def validate_task(ctx, param, value):
    config_path = ctx.params.get("config", DEFAULT_CONFIG)  # Use the default config if not provided
    config_data = read_config(config_path)
    available_tasks = list(config_data)  # Get the list of tasks from the dictionary

    if value not in available_tasks:
        raise click.BadParameter(f"'{value}' is not a valid task. Choose from: {', '.join(available_tasks)}")

    return value

def version_callback(ctx, param, value):
    """Callback for supplying version information"""
    if not value or ctx.resilient_parsing:
        return

    from snakemake import __version__ as smv
    from microcat import __version__ as mcv

    versions = {
        "Snakemake": smv,
        "MicroCAT": mcv
    }
    versions_string = "\n".join(f"{tt:<18}: {vv}" for tt, vv in versions.items())
    click.echo(versions_string)
    ctx.exit()

def abort_if_false(ctx, param, value):
    if not value:
        ctx.abort()


# ##############################################################################
# Microcat Main function Init

@click.group("microcat",cls=HelpColorsGroup,
            context_settings=CONTEXT_SETTINGS,
            help_headers_color='yellow',
            help_options_color='green')
@click.option("-v",'--version',
    is_flag=True,
    expose_value=False,
    is_eager=True,
    help="Show the version and exit.",
    callback=version_callback,
)
def microcat():
    """
            ███╗   ███╗██╗ ██████╗██████╗  ██████╗  ██████╗ █████╗ ████████╗
            ████╗ ████║██║██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔══██╗╚══██╔══╝
            ██╔████╔██║██║██║     ██████╔╝██║   ██║██║     ███████║   ██║   
            ██║╚██╔╝██║██║██║     ██╔══██╗██║   ██║██║     ██╔══██║   ██║   
            ██║ ╚═╝ ██║██║╚██████╗██║  ██║╚██████╔╝╚██████╗██║  ██║   ██║   
            ╚═╝     ╚═╝╚═╝ ╚═════╝╚═╝  ╚═╝ ╚═════╝  ╚═════╝╚═╝  ╚═╝   ╚═╝   
    Microbiome Identification upon Cell Resolution from Omics-Computational Analysis Toolbox
    """
    pass




@microcat.command("path")
def path():
    """
    Print out microcat install path
    """
    import microcat
    MICROCAT_DIR = microcat.__path__[0]
    click.echo(MICROCAT_DIR)

# ##############################################################################
# Utils module

@microcat.group("download")
def download():
    """
    Download necessary files for running microcat
    """
    pass

@download.command("cellline")
def download_celline():
    """
    Download cell line per million microbiome reads (rpmm) percentile data for single analysis pipeline
    """

    script_path = os.path.join(os.path.dirname(__file__),"download_celline.py")  
    command = ["python", script_path]

    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError as e:
        click.secho(f"Error: Failed to run download_celline.py ({e})",fg="red")



@download.command("whitelist")
def download_barcode():
    """
    Download barcode whitelist for single analysis pipeline
    """
    script_path = os.path.join(os.path.dirname(__file__),"prepare.py")  
    command = ["python", script_path]

    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError as e:
        click.secho(f"Error: Failed to run prepare.py ({e})",fg="red")

@download.command("profile")
@click.option(
    "--cluster",
    default=None,
    type=click.Choice(["slurm", "sge", "lsf"]),
    help="Cluster workflow manager engine, now support generic",
)
def download_profile(cluster):
    """
    Download profile config from Github

    $ microcat download profile --cluster lsf
    \n
    $ microcat download profile --cluster slurm
    \n
    $ microcat download profile --cluster sge
    """
    # create configuration directory that snakemake searches for profiles
    microcat_dir = os.path.dirname(os.path.abspath(__file__))
    profile_dir = os.path.join(microcat_dir, "profiles")

    # Create the profile directory if it doesn't exist
    os.makedirs(profile_dir, exist_ok=True)

    # Use the cluster variable to form the template string
    template = f"gh:Snakemake-Profiles/{cluster}"

    click.echo("Downloading profile config from github")
    try:
        # Execute cookiecutter command using subprocess
        subprocess.run(["cookiecutter", "--output-dir", profile_dir, template])
    except subprocess.CalledProcessError as e:
        click.secho(f"Error: Failed to download ({e})", fg="red")


# ##############################################################################
# Init module

@microcat.group("init", context_settings=CONTEXT_SETTINGS)
def init():
    """
    Init microcat style analysis project

    """
    pass

# ##############################################################################
## Single Init
@init.command("single")
@click.option("-b",'--begin',
              type=click.Choice(['host','classifier','denosing'], case_sensitive=False),
              default='host', 
              show_default=True,
              help='Pipeline starting point', 
              show_choices=True)
@click.option('--host',
              default='starsolo',
              show_default=True,
              help='Which hoster used',
              show_choices=True,
              type=click.Choice(['starsolo', 'cellranger'], case_sensitive=False),
              prompt='Select what host process software you use')
@click.option('--chemistry', type=click.Choice(['smartseq', 'smartseq2', 'tenx_3pv1', 'tenx_3pv2',
                                                'tenx_3pv3', 'seqwell', 'tenx_auto', 'dropseq',
                                                'tenx_multiome', 'tenx_5ppe', 'seqwell', 'indrop_v1','indrop_v2','celseq2']),
              default="tenx_auto",
              show_default=True,
              help='Sequencing chemistry option, required when host is starsolo',
              prompt='Select sequencing chemistry you use')
@click.option('--classifier',
              default='kraken2uniq',
              show_default=True,
              help='Which classifier used',
              show_choices=True,
              type=click.Choice(['kraken2uniq', 'krakenuniq', 'pathseq', 'metaphlan'], case_sensitive=False),
              prompt='Select classifier software you use')
@click.option("-s", "--samples",
              type=click.Path(exists=True),
              prompt='Give sample tsv path',
              help="Samples list, tsv format required.")
@click.option(
    "-d",
    "--workdir",
    metavar="WORKDIR",
    type=str,
    default="./",
    help="Project workdir")
def single_init(
    begin,
    workdir,
    host,
    chemistry,
    classifier,
    samples):
    """
    Single-cell RNA-seq microbiome mining pipeline init.\n

    Init 10x Genomics Example:
    \n
    $ microcat init single --host starsolo --chemistry tenx_auto --classifier kraken2uniq
    \n
    Init Smartseq2 Example:
    \n
    $ microcat init single --host starsolo --chemistry smartseq2 --classifier kraken2uniq
    """
    if workdir:
    # Create a MicrocatConfig object using the provided working directory
        project = MicrocatConfig(workdir,config_type="single")
    
        # Check if the working directory already exists
        if os.path.exists(workdir):
            click.secho(f"WARNING: The working directory '{workdir}' already exists.",fg="yellow")
            proceed = input("Do you want to proceed? (y/n): ").lower()
            if proceed != 'y':
                print("Aborted.")
                sys.exit(1)

        # Print the project structure and create the necessary subdirectories
        print(project.__str__())
        project.create_dirs()

        # Get the default configuration
        conf = project.get_config()

        # Update environment configuration file paths
        for env_name in conf["envs"]:
            conf["envs"][env_name] = os.path.join(os.path.realpath(workdir), f"envs/{env_name}.yaml")


        for script_path in conf["scripts"]:
            origin_path = conf["scripts"][script_path]
            conf["scripts"][script_path] = os.path.join(os.path.dirname(__file__),"single_wf",f"{origin_path}")

        # Get the single cell chemistry defination
        with open(os.path.join(os.path.join(os.path.dirname(__file__),"chemistry_defs.json"))) as file:
            CHEMISTRY_DEFS = json.load(file)

        # conf["params"]["simulate"]["do"] = False
        conf["params"]["begin"] = begin
        for hoster_ in ["starsolo","cellranger"]:
            if hoster_ == host:
                conf["params"]["host"][hoster_]["do"] = True
                if hoster_ == "starsolo":
                    if chemistry in CHEMISTRY_DEFS:
                        chem_params = CHEMISTRY_DEFS[chemistry]

                        # 更新host.starsolo.barcode下的参数
                        barcode_params = chem_params.get("barcode")

                        if barcode_params:
                            for key, value in barcode_params[0].items():
                                conf["params"]["host"]["starsolo"]["barcode"][key] = value

                        # 更新host.starsolo.algorithm下的参数
                        algorithm_params = chem_params.get("algorithm")
                        if algorithm_params:
                            for key, value in algorithm_params[0].items():
                                conf["params"]["host"]["starsolo"]["algorithm"][key] = value

                        # 更新host.starsolo下的其他参数
                        for key, value in chem_params.items():
                            if key not in ["barcode", "algorithm"]:
                                conf["params"]["host"]["starsolo"][key] = value
            else:
                conf["params"]["host"][hoster_]["do"] = False

        for classifier_ in ["kraken2uniq","krakenuniq","pathseq","metaphlan4"]:
            if classifier_ in classifier:
                conf["params"]["classifier"][classifier_]["do"] = True
            else:
                conf["params"]["classifier"][classifier_]["do"] = False

        # Add the user-supplied samples table to the configuration
        if samples:
            conf["params"]["samples"] = os.path.abspath(samples)
        else:
            click.secho("ERROR:Please supply samples table",fg="red")
            sys.exit(-1)

        # Update the configuration file
        update_config(
            project.config_file, project.new_config_file, conf, remove=False
        )

        click.secho("NOTE: Congfig.yaml reset to default values.", fg='green')

    else:
        # If the user didn't provide a working directory, print an error message and exit
        click.secho("ERROR:Please supply a workdir!",fg='red')
        sys.exit(-1)

@init.command("bulk")
@click.option('--begin',
              type=click.Choice(['simulate','trimming','host','classifier','denosing'], case_sensitive=False),
              default='trimming', 
              show_default=True,
              help='Pipeline starting point', 
              show_choices=True,
              prompt='Your name please')
@click.option('--host',
              default='starsolo',
              show_default=True,
              help='Which hoster used',
              show_choices=True,
              type=click.Choice(['starsolo', 'cellranger'], case_sensitive=False),
              prompt='Select what host you use')
@click.option('--chemistry', type=click.Choice(['smartseq', 'smartseq2', 'tenx_3pv1', 'tenx_3pv2',
                                                'tenx_3pv3', 'seqwell', 'tenx_auto', 'dropseq',
                                                'tenx_multiome', 'tenx_5ppe', 'seqwell','indrop_v1','celseq2']),
              default=None, help='Sequencing chemistry option, required when host is starsolo')
@click.option('--classifier',
              default='pathseq',
              show_default=True,
              help='Which classifier used',
              show_choices=True,
              type=click.Choice(['kraken2uniq', 'krakenuniq', 'pathseq', 'metaphlan'], case_sensitive=False))
@click.option("-s", "--samples",
              type=click.Path(exists=True),
              prompt='Select what host you use',
              help="Samples list, tsv format required.")
@click.option(
    "-d",
    "--workdir",
    metavar="WORKDIR",
    type=str,
    default="./",
    help="Project workdir")
def bulk_init(
    workdir,
    begin,
    host,
    chemistry,
    classifier,
    samples):
    """
    Bulk seq microbiome mining pipeline init.
    """
    if workdir:
    # Create a MicrocatConfig object using the provided working directory
        project = MicrocatConfig(workdir,config_type="bulk")
    
        # Check if the working directory already exists
        if os.path.exists(workdir):
            print(f"WARNING: The working directory '{workdir}' already exists.")
            proceed = input("Do you want to proceed? (y/n): ").lower()
            if proceed != 'y':
                print("Aborted.")
                sys.exit(1)

        # Print the project structure and create the necessary subdirectories
        print(project.__str__())
        project.create_dirs()

        # Get the default configuration
        conf = project.get_config()

        # Update environment configuration file paths
        for env_name in conf["envs"]:
            conf["envs"][env_name] = os.path.join(os.path.realpath(workdir), f"envs/{env_name}.yaml")


        for script_path in conf["scripts"]:
            origin_path = conf["scripts"][script_path]
            conf["scripts"][script_path] = os.path.join(os.path.dirname(__file__),f"{origin_path}")

        # Get the single cell chemistry defination
        with open(os.path.join(os.path.join(os.path.dirname(__file__),"chemistry_defs.json"))) as file:
            CHEMISTRY_DEFS = json.load(file)

        conf["params"]["simulate"]["do"] = False
        conf["params"]["begin"] = begin
        for hoster_ in ["starsolo","cellranger"]:
            if hoster_ == host:
                conf["params"]["host"][hoster_]["do"] = True
                if hoster_ == "starsolo":
                    if chemistry in CHEMISTRY_DEFS:
                        chem_params = CHEMISTRY_DEFS[chemistry]

                        # 更新host.starsolo.barcode下的参数
                        barcode_params = chem_params.get("barcode")

                        if barcode_params:
                            for key, value in barcode_params[0].items():
                                conf["params"]["host"]["starsolo"]["barcode"][key] = value

                        # 更新host.starsolo.algorithm下的参数
                        algorithm_params = chem_params.get("algorithm")
                        if algorithm_params:
                            for key, value in algorithm_params[0].items():
                                conf["params"]["host"]["starsolo"]["algorithm"][key] = value

                        # 更新host.starsolo下的其他参数
                        for key, value in chem_params.items():
                            if key not in ["barcode", "algorithm"]:
                                conf["params"]["host"]["starsolo"][key] = value
            else:
                conf["params"]["host"][hoster_]["do"] = False

        for classifier_ in ["kraken2uniq","krakenuniq","pathseq","metaphlan"]:
            if classifier_ in classifier:
                conf["params"]["classifier"][classifier_]["do"] = True
            else:
                conf["params"]["classifier"][classifier_]["do"] = False

        # Add the user-supplied samples table to the configuration
        if samples:
            conf["params"]["samples"] = os.path.abspath(samples)
        else:
            click.secho("ERROR:Please supply samples table!",fg='red')
            sys.exit(-1)

        # Update the configuration file
        update_config(
            project.config_file, project.new_config_file, conf, remove=False
        )

        click.secho("NOTE: Congfig.yaml reset to default values.", fg='green')

    else:
        # If the user didn't provide a working directory, print an error message and exit
        click.secho("ERROR:Please supply a workdir!",fg='red')
        sys.exit(-1)


@init.command("spatial")
@click.option('--begin',
              type=click.Choice(['simulate','trimming','host','classifier','denosing'], case_sensitive=False),
              default='trimming', 
              show_default=True,
              help='Pipeline starting point', 
              show_choices=True,
              prompt='Your name please')
@click.option('--host',
              default='starsolo',
              show_default=True,
              help='Which hoster used',
              show_choices=True,
              type=click.Choice(['starsolo', 'cellranger'], case_sensitive=False),
              prompt='Select what host you use')
@click.option('--chemistry', type=click.Choice(['smartseq', 'smartseq2', 'tenx_3pv1', 'tenx_3pv2',
                                                'tenx_3pv3', 'seqwell', 'tenx_auto', 'dropseq',
                                                'tenx_multiome', 'tenx_5ppe', 'seqwell', 'celseq2']),
              default=None, help='Sequencing chemistry option, required when host is starsolo')
@click.option('--classifier',
              default='pathseq',
              show_default=True,
              help='Which classifier used',
              show_choices=True,
              type=click.Choice(['kraken2uniq', 'krakenuniq', 'pathseq', 'metaphlan'], case_sensitive=False))
@click.option("-s", "--samples",
              type=click.Path(exists=True),
              prompt='Select what host you use',
              help="Samples list, tsv format required.")
def spatial_init(hash_type,host,chemistry,classifier,samples):
    """Single-cell RNA-seq microbiome mining pipeline init."""
    click.echo(hash_type)
    click.secho('Hello World!', fg='green')
    click.echo('Initialized the database')

@init.command("multi")
@click.option('--begin',
              type=click.Choice(['simulate','trimming','host','classifier','denosing'], case_sensitive=False),
              default='trimming', 
              show_default=True,
              help='Pipeline starting point', 
              show_choices=True,
              prompt='Your name please')
@click.option('--host',
              default='starsolo',
              show_default=True,
              help='Which hoster used',
              show_choices=True,
              type=click.Choice(['starsolo', 'cellranger'], case_sensitive=False),
              prompt='Select what host you use')
@click.option('--chemistry', type=click.Choice(['smartseq', 'smartseq2', 'tenx_3pv1', 'tenx_3pv2',
                                                'tenx_3pv3', 'seqwell', 'tenx_auto', 'dropseq',
                                                'tenx_multiome', 'tenx_5ppe', 'seqwell', 'celseq2']),
              default=None, help='Sequencing chemistry option, required when host is starsolo')
@click.option('--classifier',
              default='pathseq',
              show_default=True,
              help='Which classifier used',
              show_choices=True,
              type=click.Choice(['kraken2uniq', 'krakenuniq', 'pathseq', 'metaphlan'], case_sensitive=False))
@click.option("-s", "--samples",
              type=click.Path(exists=True),
              prompt='Select what host you use',
              help="Samples list, tsv format required.")
def multi_init(hash_type,host,chemistry,classifier,samples):
    """Single-cell RNA-seq microbiome mining pipeline init."""
    click.echo(hash_type)
    click.secho('Hello World!', fg='green')
    click.echo('Initialized the database')



# ##############################################################################
# Config module

@microcat.command("config", context_settings=CONTEXT_SETTINGS)
@click.option("--star_ref",
              type=click.Path(exists=True),
              help="STAR reference genome for aligning bulk sequencing data")
@click.option("--starsolo_ref",
              type=click.Path(exists=True),
              help="STARsolo reference genome for aligning single and spatial sequencing data")
@click.option("--cellranger_ref",
              type=click.Path(exists=True),
              help="Cellranger reference genome for aligning single sequencing data")
@click.option("--krak2_ref",
              type=click.Path(exists=True),
              help="Kraken2 reference genome for aligning microbiome data")
@click.option("--krakuniq_ref",
              type=click.Path(exists=True),
              help="Samples list, tsv format required.")
@click.option("--pathseq_ref",
              type=click.Path(exists=True),
              help="Samples list, tsv format required.")
@click.option("--metaphlan4_ref",
              type=click.Path(exists=True),
              help="Samples list, tsv format required.")
@click.option("--generic_profile",
              type=click.Choice(["lsf","slurm","pbs"]),
              multiple=False,
              help="Samples list, tsv format required.")
@click.option("-e","--edit",
              default=None,
              type=click.Choice(["single", "bulk", "spatial","multi"]),
              multiple=False,
              help="Edit config file")
def config(
    star_ref,
    starsolo_ref,
    krak2_ref,
    krakuniq_ref,
    pathseq_ref,
    metaphlan4_ref,
    cellranger_ref,
    generic_profile,
    edit):
    """
    Quickly adjust microcat's default configurations
    """

    microcat_dir = os.path.dirname(os.path.abspath(__file__))
    # bulk_config_file = os.path.join(microcat_dir,"bulk_wf","config", "config.yaml")
    # multi_config_file = os.path.join(microcat_dir,"multi_wf","config", "config.yaml")
    single_config_file = os.path.join(microcat_dir,"single_wf","config", "config.yaml")
    # spatial_config_file = os.path.join(microcat_dir,"spatial_wf","config", "config.yaml")

    edited = False  # Flag to check if the user has edited any config file

    if edit:
        if edit == "single":
            if os.path.exists(single_config_file):
                click.edit(filename=single_config_file)
                edited = True
            else:
                click.echo(f"Config file not found: {single_config_file}")
        elif edit == "bulk":
            if os.path.exists(bulk_config_file):
                click.edit(filename=bulk_config_file)
                edited = True
            else:
                click.echo(f"Config file not found: {bulk_config_file}")
        if edit == "spatial":
            if os.path.exists(spatial_config_file):
                click.edit(filename=spatial_config_file)
                edited = True
            else:
                click.echo(f"Config file not found: {spatial_config_file}")
        if edit == "multi":
            if os.path.exists(multi_config_file):
                click.edit(filename=multi_config_file)
                edited = True
            else:
                click.echo(f"Config file not found: {multi_config_file}")

    # If a specific configuration file is specified and has not been edited yet, update it
    if cellranger_ref is not None and not edited:
        # Single have cellranger refernence 
        single_config = parse_yaml(single_config_file)
        single_config["params"]["host"]["cellranger"]["reference"] = cellranger_ref
        update_config(single_config_file, single_config_file, single_config)

        click.echo(f"Cellranger Refernce location updated as '{cellranger_ref}' path.")
    if starsolo_ref is not None and not edited:
        # Spatial and single both have starsolo refernence 
        single_config = parse_yaml(single_config_file)
        single_config["params"]["host"]["starsolo"]["reference"] = starsolo_ref
        update_config(single_config_file, single_config_file, single_config)

        # spatial_config = parse_yaml(spatial_config_file)
        # spatial_config["host"]["starsolo"]["reference"] = starsolo_ref
        # update_config(spatial_config_file, spatial_config_file, spatial_config)

        click.echo(f"Starsolo Refernce location updated as '{starsolo_ref}' path.")
    if krak2_ref is not None and not edited:
        # Spatial and single both have starsolo refernence 
        single_config = parse_yaml(single_config_file)
        single_config["params"]["classifier"]["kraken2uniq"]["kraken2_database"] = krak2_ref
        update_config(single_config_file, single_config_file, single_config)

        # spatial_config = parse_yaml(spatial_config_file)
        # spatial_config["params"]["classifier"]["kraken2uniq"]["kraken2_database"] = krak2_ref
        # update_config(spatial_config_file, spatial_config_file, spatial_config)

        # bulk_config = parse_yaml(bulk_config_file)
        # bulk_config["params"]["classifier"]["kraken2uniq"]["kraken2_database"] = krak2_ref
        # update_config(bulk_config_file, bulk_config_file, bulk_config)

        click.echo(f"Krake2nuniq Refernce location updated as '{krak2_ref}' path.")

    if krakuniq_ref is not None and not edited:
        # Spatial and single both have starsolo refernence 
        single_config = parse_yaml(single_config_file)
        single_config["params"]["classifier"]["krakenuniq"]["krakenuniq_database"] = krakuniq_ref
        update_config(single_config_file, single_config_file, single_config)

        # spatial_config = parse_yaml(spatial_config_file)
        # spatial_config["params"]["classifier"]["krakenuniq"]["krakenuniq_database"] = krakuniq_ref
        # update_config(spatial_config_file, spatial_config_file, spatial_config)

        # bulk_config = parse_yaml(bulk_config_file)
        # bulk_config["params"]["classifier"]["krakenuniq"]["krakenuniq_database"] = krakuniq_ref
        # update_config(bulk_config_file, bulk_config_file, bulk_config)

        click.echo(f"Krakenuniq Refernce location updated as '{krakuniq_ref}' path.")
    if generic_profile is not None :
        profile_path = os.path.join(microcat_dir,"profiles", "generic")
        profile_conf_file = os.path.join(profile_path, "config.yaml")
        submmit_conf_file = os.path.join(profile_path, "key_mapping.yaml")

        profile_conf =  parse_yaml(profile_conf_file)
        submmit_conf = parse_yaml(submmit_conf_file)
        
        submmit_conf['system'] = generic_profile
        profile_conf['cluster-status'] = generic_profile+"_status.py"
        
        update_config(profile_conf_file, profile_conf_file, submmit_conf)
        update_config(submmit_conf_file, submmit_conf_file, profile_conf)
        click.echo(f"Cluster Configuration updated as '{generic_profile}' system.")


# ##############################################################################
# Debug module
@microcat.command("debug", context_settings=CONTEXT_SETTINGS)
@click.option(
    "--config",
    type=str,
    default=DEFAULT_CONFIG,
    help="Path of config.yaml",
    show_default=True
)
@click.option("--task",
              default="all",
              callback=validate_task)
@click.pass_context
def run_debug(
    ctx,
    config,
    task):
    """
    Execute the analysis workflow on debug mode.
    """
    # Retrieve the 'config' value from the 'context' object
    conf = parse_yaml(config)
    # Get the workflow type
    if conf["workflow"] == "single":
        snakefile = os.path.join(os.path.dirname(__file__),"single_wf/snakefiles/single_wf.smk")
    elif conf["workflow"] == "bulk":
        snakefile = os.path.join(os.path.dirname(__file__),"bulk_wf/snakefiles/bulk_wf.smk")
    elif conf["workflow"] == "spatial":
        snakefile = os.path.join(os.path.dirname(__file__),"spatial_wf/snakefiles/spatial_wf.smk")
    elif conf["workflow"] == "multi":
        snakefile = os.path.join(os.path.dirname(__file__),"multi_wf/snakefiles/multi_wf.smk")
    else:
        raise ValueError("Invalid config_type")
    
    unknown = click.get_current_context().args

    # Check if the sample list is provided, exit if not
    if not os.path.exists(conf["params"]["samples"]):
        click.secho("ERROR:Please specific samples list on init step or change config.yaml manualy.",fg="red")
        sys.exit(1)

    # Prepare the command list for running Snakemake
    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--configfile",
        config,
        "--until",
        task,
        "--dry-run"
    ] + unknown

    # Convert the command list to a string and print it
    cmd_str = " ".join(cmd).strip()
    click.echo("Running microcat %s:\n%s" % (conf["workflow"], cmd_str))

    # Execute the Snakemake command and capture the output
    env = os.environ.copy()
    proc = subprocess.Popen(
        cmd_str,
        shell=True,
        stdout=sys.stdout,
        stderr=sys.stderr,
        env=env,
    )
    proc.communicate()

    # Print the actual executed command
    click.echo(f'''\nReal running cmd:\n{cmd_str}''')

# ##############################################################################
# Run local module
@microcat.command("run-local", context_settings=CONTEXT_SETTINGS)
@click.option(
    "--cores",
    type=int,
    default=60,
    help="All job cores, available on 'run-local'",
)
@click.option(
    "--config",
    type=str,
    default=DEFAULT_CONFIG,
    help="Path of config.yaml",
)
@click.option(
    "--wait",
    type=int,
    default=60,
    help="Wait given seconds",
)
@click.option(
    "--use-conda",
    default=False,
    is_flag=True,
    help="Use conda environment",
)
@click.option(
    "--list",
    default=False,
    is_flag=True,
    help="List",
)
@click.option(
    "--conda-prefix",
    default="~/.conda/envs",
    help="Conda environment prefix",
)
@click.option(
    "--conda-create-envs-only",
    default=False,
    is_flag=True,
    help="Conda create environments only",
)
@click.option(
    "--use-singularity",
    default=False,
    is_flag=True,
    help="Use a singularity container",
)
@click.option(
    "--singularity-prefix",
    default="",
    help="Singularity images prefix",
)
@click.option(
    "--jobs",
    type=int,
    default=30,
    help="Job numbers",
)
@click.option("--task",
              default="all",
              callback=validate_task) 
@click.pass_context
def run_local(
    ctx,
    config,
    cores,
    task,
    use_conda,
    list,
    conda_prefix,
    conda_create_envs_only,
    use_singularity,
    singularity_prefix,
    jobs,
    wait):
    """
    Execute the analysis workflow on local computer mode
    """
    from snakemake import __version__ as smv
    from microcat import __version__ as mcv

    versions = {
        "Snakemake": smv,
        "MicroCAT": mcv
    }
    versions_string = "\n".join(f"{tt:<18}: {vv}" for tt, vv in versions.items())
    click.echo("---------------------------")
    click.echo("Version inforamtion for microcat:")
    click.echo(versions_string)
    click.echo("---------------------------")
    # Parse the YAML configuration file
    conf = parse_yaml(config)
    # Get the workflow type
    if conf["workflow"] == "single":
        snakefile = os.path.join(os.path.dirname(__file__),"single_wf/snakefiles/single_wf.smk")
    elif conf["workflow"] == "bulk":
        snakefile = os.path.join(os.path.dirname(__file__),"bulk_wf/snakefiles/bulk_wf.smk")
    elif conf["workflow"] == "spatial":
        snakefile = os.path.join(os.path.dirname(__file__),"spatial_wf/snakefiles/spatial_wf.smk")
    elif conf["workflow"] == "multi":
        snakefile = os.path.join(os.path.dirname(__file__),"multi_wf/snakefiles/multi_wf.smk")
    else:
        raise ValueError("Invalid config_type")
    
    unknown = click.get_current_context().args

    # Check if the sample list is provided, exit if not
    if not os.path.exists(conf["params"]["samples"]):
        click.secho("ERROR:Please specific samples list on init step or change config.yaml manualy.",fg="red")
        sys.exit(1)

    # Prepare the command list for running Snakemake
    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--configfile",
        config,
        "--until",
        task
    ] + unknown

    # Add specific flags to the command based on the input arguments
    if "--touch" in unknown:
        pass
    elif conda_create_envs_only:
        cmd += ["--use-conda",
                "--conda-create-envs-only",
                "--cores 4"]
        if conda_prefix is not None:
            cmd += ["--conda-prefix", conda_prefix]
    else:
        cmd += [
            "--rerun-incomplete",
            "--keep-going",
            "--printshellcmds",
            "--show-failed-logs",
            "--reason",
        ]

        # Add flags for using conda environments
        if use_conda:
            cmd += ["--use-conda"]
            if conda_prefix is not None:
                cmd += ["--conda-prefix", conda_prefix]
        
        # Add flags for listing tasks
        if list:
            cmd += ["--list"]

        # Get the system's available resources
        num_cores = psutil.cpu_count(logical=False)
        mem_usage = round((float(psutil.virtual_memory().total)/1024/1024),2)
        # Adjust the resource parameters based on available resources
        cores_limit = int(num_cores * 0.80)
        mem_limit = int(mem_usage * 0.80)

        # Compare the requested resources with the limits
        if cores > num_cores:
            click.secho("WARNING: Maximum allowable cores is exceeded. Automatically set to 75% of the current available cores ({0}).".format(cores_limit),fg="yellow")
            cores = cores_limit
        # Initialize the resources variable before conditional assignment
        resources = {}
        # Check if "--resources" is present in unknown arguments
        
        resources_index = None
        for index, arg in enumerate(unknown):
            if arg == "--resources":
                resources_index = index
                break

        if resources_index is not None:
            resources_args = unknown[resources_index + 1]
            # resources = {}

            # Split the resources_args into key-value pairs
            resources_items = resources_args.split()
            
            # Find the index of the next parameter after resources_args
            next_parameter_index = resources_index + 2 if resources_index + 2 < len(unknown) else None
            
            # If there is a next parameter and it doesn't start with "--", include it in resources_items
            if next_parameter_index is not None and not unknown[next_parameter_index].startswith("--"):
                resources_items.append(unknown[next_parameter_index])

            # Extract key-value pairs from resources_items
            for item in resources_items:
                key, value = item.split("=")
                resources[key] = value

                # Get user setting mem_mb
                if key == "mem_mb":
                    user_mem_mb = value

            if "mem_mb" not in resources:
                user_mem_mb = None
        else:
            user_mem_mb = None

        # Check if mem_mb_value exceeds mem_usage or is None, and adjust it if necessary
        if user_mem_mb is None or (int(user_mem_mb) > mem_limit or int(user_mem_mb) <= 0):
            click.secho("WARNING: Maximum allowable memory is not set or exceeded. Automatically set to 75% of the current memory ({0}MB).".format(mem_limit),fg="yellow")
            user_mem_mb = mem_limit
            resources["mem_mb"] = user_mem_mb

        # Create a string representation of the resources dictionary
        resources_str = " ".join(["{0}={1}".format(key, value) for key, value in resources.items()])
        cmd += ["--cores", str(cores),
                "--jobs", str(jobs),
                "--resources", str(resources_str)]
    # Convert the command list to a string and print it
    cmd_str = " ".join(cmd).strip()
    click.echo("Running microcat %s:\n%s" % (conf["workflow"], cmd_str))

    # Execute the Snakemake command and capture the output
    env = os.environ.copy()
    proc = subprocess.Popen(
        cmd_str,
        shell=True,
        stdout=sys.stdout,
        stderr=sys.stderr,
        env=env,
    )
    proc.communicate()

    # Print the actual executed command
    click.echo(f'''\nReal running cmd:\n{cmd_str}''')



@microcat.command("run-remote", context_settings=CONTEXT_SETTINGS)
@click.option(
    "--jobs",
    type=int,
    default=30,
    help="Cluster job numbers, available on '--run-remote'",
)
@click.option(
    "--cluster-engine",
    default="generic",
    type=click.Choice(["generic","slurm", "sge", "lsf","custom"]),
    help="Cluster workflow manager engine, now support generic",
)
@click.option(
    "--local-cores",
    type=int,
    default=8,
    help="Local cores, available on '--run-remote'",
)
@click.option(
    "--use-conda",
    default=False,
    is_flag=True,
    help="Use conda environment",
)
@click.option(
    "--config",
    type=str,
    default=DEFAULT_CONFIG,
    help="Path of config.yaml",
)
@click.option(
    "--conda-prefix",
    default="~/.conda/envs",
    help="Conda environment prefix",
)
@click.option(
    "--conda-create-envs-only",
    default=False,
    is_flag=True,
    help="Conda create environments only",
)
@click.option(
    "--use-singularity",
    default=False,
    is_flag=True,
    help="Use a singularity container",
)
@click.option(
    "--list",
    default=False,
    is_flag=True,
    help="List rule",
)
@click.option(
    "--singularity-prefix",
    default="",
    help="Singularity images prefix",
)
# @click.argument("task",callback=validate_task, metavar="[TASK]") 
@click.option("--task",
              default="all",
              callback=validate_task) 
@click.pass_context
def run_remote(
    ctx,
    config,
    task,
    use_conda,
    conda_prefix,
    conda_create_envs_only,
    use_singularity,
    singularity_prefix,
    jobs,
    list,
    local_cores,
    cluster_engine
    ):
    """
    Execute the analysis workflow on remote cluster mode
    """

    # Parse the YAML configuration file
    conf = parse_yaml(config)
    # Get the workflow type
    if conf["workflow"] == "single":
        snakefile = os.path.join(os.path.dirname(__file__),"single_wf/snakefiles/single_wf.smk")
    elif conf["workflow"] == "bulk":
        snakefile = os.path.join(os.path.dirname(__file__),"bulk_wf/snakefiles/bulk_wf.smk")
    elif conf["workflow"] == "spatial":
        snakefile = os.path.join(os.path.dirname(__file__),"spatial_wf/snakefiles/spatial_wf.smk")
    elif conf["workflow"] == "multi":
        snakefile = os.path.join(os.path.dirname(__file__),"multi_wf/snakefiles/multi_wf.smk")
    else:
        raise ValueError("Invalid config_type")
    
    unknown = click.get_current_context().args

    # Check if the sample list is provided, exit if not
    if not os.path.exists(conf["params"]["samples"]):
        click.secho("ERROR:Please specific samples list on init step or change config.yaml manualy.",fg="red")
        sys.exit(1)

    # Prepare the command list for running Snakemake
    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--configfile",
        config,
        "--until",
        task
    ] + unknown

    # Add specific flags to the command based on the input arguments
    if "--touch" in unknown:
        pass
    elif conda_create_envs_only:
        cmd += ["--use-conda",
                "--conda-create-envs-only",
                "--cores", 
                str(local_cores)]
        if conda_prefix is not None:
            cmd += ["--conda-prefix", conda_prefix]
    else:
        cmd += [
            "--rerun-incomplete",
            "--keep-going",
            "--printshellcmds",
            "--reason",
        ]

        # Add flags for using conda environments
        if use_conda:
            cmd += ["--use-conda"]
            if conda_prefix is not None:
                cmd += ["--conda-prefix", conda_prefix]
        
        # Add flags for listing tasks
        if list:
            cmd += ["--list"]


        profile_path = os.path.join("./profiles", cluster_engine)
        if cluster_engine == "generic":
            profile_conf = parse_yaml(os.path.join(profile_path, "config.yaml"))
            submmit_conf = parse_yaml(os.path.join(profile_path, "key_mapping.yaml"))
            cluster_setting = submmit_conf['system']
            cluster_status_setting = profile_conf['cluster-status']
            if cluster_status_setting != cluster_setting + "_status.py":
                click.secho(f"ERROR：Cluster system and Cluster system status script dont match！Please check profile.\n"
                            f"cluster: {cluster_setting}\n"
                            f"cluster-status: {cluster_status_setting}", fg="red")
            else:
                click.secho(f"Cluster system is {cluster_setting}",fg="green")
                queue_file = os.path.join(profile_path, "queue.tsv")   
                if os.path.exists(queue_file):
                    click.secho("Snakemake will automatically selects the optimal queue for you.",fg="green")
                    df = pd.read_csv(queue_file, sep="\t")
                    click.echo(df.to_string(index=False))
                else:
                    click.secho("Snakemake will selects the queue with config",fg="green")
        else:
            if not os.path.exists(profile_path):
                click.secho(f"ERROR: Cluster profile {profile_path} don't exist!", fg="red")
                sys.exit(1)

        cmd += ["--profile", profile_path,
                "--local-cores", str(local_cores),
                "--jobs", str(jobs)]
    
    # Convert the command list to a string and print it
    cmd_str = " ".join(cmd).strip()
    click.echo("Running microcat %s:\n%s" % (conf["workflow"], cmd_str))

    # Execute the Snakemake command and capture the output
    env = os.environ.copy()
    proc = subprocess.Popen(
        cmd_str,
        shell=True,
        stdout=sys.stdout,
        stderr=sys.stderr,
        env=env,
    )
    proc.communicate()

    # Print the actual executed command
    click.echo(f'''\nReal running cmd:\n{cmd_str}''')



if __name__ == '__main__':
    microcat()