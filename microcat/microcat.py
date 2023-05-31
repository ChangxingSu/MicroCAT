#!/usr/bin/env python
import argparse
import os
import subprocess
import sys
import textwrap
from io import StringIO
import configer
import pandas as pd

WORKFLOWS_MAG = [
    "host_all",
    "check_all",]

def update_config_tools(conf, host, classifier,assay=None):
    conf["params"]["simulate"]["do"] = False
    # conf["params"]["begin"] = begin

    # for trimmer_ in ["sickle", "fastp"]:
    #     if trimmer_ == trimmer:
    #         conf["params"]["trimming"][trimmer_]["do"] = True
    #     else:
    #         conf["params"]["trimming"][trimmer_]["do"] = False

    for hoster_ in ["starsolo","cellranger"]:
        if hoster_ == host:
            conf["params"]["host"][hoster_]["do"] = True
            if hoster_ == "starsolo":
                conf["params"]["host"][hoster_]["assay"] = assay
        else:
            conf["params"]["host"][hoster_]["do"] = False

    for classifier_ in ["kraken2uniq","krakenuniq","pathseq"]:
        if classifier_ in classifier:
            conf["params"]["classifier"][classifier_]["do"] = True
        else:
            conf["params"]["classifier"][classifier_]["do"] = False

    # if begin == "simulate":
    #     conf["params"]["simulate"]["do"] = True
    # elif begin == "rmhost":
    #     conf["params"]["trimming"][trimmer]["do"] = False
    # elif (begin == "assembly") or (begin == "binning"):
    #     conf["params"]["raw"]["save_reads"] = True
    #     conf["params"]["raw"]["fastqc"]["do"] = False
    #     conf["params"]["qcreport"]["do"] = False

    #     conf["params"]["trimming"][trimmer]["do"] = False
    #     conf["params"]["rmhost"][host]["do"] = False
    return conf


def init(args, unknown):

    # Check if the user provided a working directory
    if args.workdir:
        # Create a MicrocatConfig object using the provided working directory
        project = configer.MicrocatConfig(args.workdir)


        # Check if the working directory already exists
        if os.path.exists(args.workdir):
            print(f"Warning: The working directory '{args.workdir}' already exists.")
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
            conf["envs"][env_name] = os.path.join(os.path.realpath(args.workdir), f"envs/{env_name}.yaml")

        # Update the configuration with the selected tools
        conf = update_config_tools(
            conf,args.host,args.classifier,args.assay
        )

        # Add the user-supplied samples table to the configuration
        if args.samples:
            conf["params"]["samples"] = args.samples
        else:
            print("Please supply samples table")
            sys.exit(-1)

        # Update the configuration file
        configer.update_config(
            project.config_file, project.new_config_file, conf, remove=False
        )

    else:
        # If the user didn't provide a working directory, print an error message and exit
        print("Please supply a workdir!")
        sys.exit(-1)


def run_snakemake(args, unknown, snakefile, workflow):
    """
    Use subprocess.Popen to run the MicroCAT workflow.

    Args:
        args (argparse.Namespace): An object containing parsed command-line arguments.
    """
    # Parse the YAML configuration file
    conf = configer.parse_yaml(args.config)

    # Check if the sample list is provided, exit if not
    if not os.path.exists(conf["params"]["samples"]):
        print("Please specific samples list on init step or change config.yaml manualy")
        sys.exit(1)

    # Prepare the command list for running Snakemake
    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--configfile",
        args.config,
        "--cores",
        str(args.cores),
    ] + unknown

    # Add specific flags to the command based on the input arguments
    if "--touch" in unknown:
        pass
    elif args.conda_create_envs_only:
        cmd += ["--use-conda", "--conda-create-envs-only"]
        if args.conda_prefix is not None:
            cmd += ["--conda-prefix", args.conda_prefix]
    else:
        cmd += [
            "--rerun-incomplete",
            "--keep-going",
            "--printshellcmds",
            "--reason",
        ]

        # Add flags for using conda environments
        if args.use_conda:
            cmd += ["--use-conda"]
            if args.conda_prefix is not None:
                cmd += ["--conda-prefix", args.conda_prefix]
        
        # Add flags for listing tasks
        if args.list:
            cmd += ["--list"]

        # Add flags for running tasks locally
        elif args.run_local:
            cmd += ["--cores", str(args.cores),
                    "--local-cores", str(args.local_cores),
                    "--jobs", str(args.jobs)]
        elif args.run_remote:
            profile_path = os.path.join("./profiles", args.cluster_engine)
            cmd += ["--profile", profile_path,
                    "--local-cores", str(args.local_cores),
                    "--jobs", str(args.jobs)]
        
        # Add flags for running tasks remotely
        elif args.debug:
            cmd += ["--debug-dag"]
        # Add flags for a dry run
        else:
            cmd += ["--dry-run"]

        # Add --dry-run flag if it's specified and not already in the command list
        if args.dry_run and ("--dry-run" not in cmd):
            cmd += ["--dry-run"]
    
    # Convert the command list to a string and print it
    cmd_str = " ".join(cmd).strip()
    print("Running microcat %s:\n%s" % (workflow, cmd_str))

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
    print(f'''\nReal running cmd:\n{cmd_str}''')


def bulk_wf(args, unknown):
    print("bulk")
    snakefile = os.path.join(os.path.dirname(__file__), "snakefiles/bulk_wf.smk")
    run_snakemake(args, unknown, snakefile, "bulk_wf")

def scRNA_wf(args, unknown):
    snakefile = os.path.join(os.path.dirname(__file__), "snakefiles/scRNA_wf.smk")
    run_snakemake(args, unknown, snakefile, "scRNA_wf")

def spatial_wf(args, unknown):
    snakefile = os.path.join(os.path.dirname(__file__), "snakefiles/spatial_wf.smk")
    run_snakemake(args, unknown, snakefile, "spatial_wf")

def main():

    # Banner text and program information
    banner = '''
            ███╗   ███╗██╗ ██████╗██████╗  ██████╗  ██████╗ █████╗ ████████╗
            ████╗ ████║██║██╔════╝██╔══██╗██╔═══██╗██╔════╝██╔══██╗╚══██╔══╝
            ██╔████╔██║██║██║     ██████╔╝██║   ██║██║     ███████║   ██║   
            ██║╚██╔╝██║██║██║     ██╔══██╗██║   ██║██║     ██╔══██║   ██║   
            ██║ ╚═╝ ██║██║╚██████╗██║  ██║╚██████╔╝╚██████╗██║  ██║   ██║   
            ╚═╝     ╚═╝╚═╝ ╚═════╝╚═╝  ╚═╝ ╚═════╝  ╚═════╝╚═╝  ╚═╝   ╚═╝   
              *         *     *        *   *     _      *      *      *
                   *         *                  / )        *     *
              *            *     (\__/) *       ( (  *       *      *
                         ,-.,-.,)     (.,-.,-.,-.) ).,-.,-.
                |  |  |  | @|  ={      }= | @|  / / | @|o |
    j__j__j__j__j__j__j__j__j__j_)     `-------/ /__j__j__j__j__j__j__j__j_
    ____________________________(               /__________________________
                   |    |  | @| \              || o|O | @|
                   |    |o |  |,'\       ,   ,'"|  |  |  |  hjw
                   |vV\ vV\|/vV|`-'\  ,---\   | \Vv\hjwVv\//v
                                  _) )  `. \ /
                                 (__/      ) )
                                          (_/

    Microbiome Identification in Cell Resolution from Omics-Computational Analysis Toolbox
    
    Version: dev
    Manual: https://github.com/ChangxingSu/microcat/wiki
    Source code: https://github.com/ChangxingSu/microcat
    '''

    # Create the main parser object with a banner and program name
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(banner),
        prog="microcat",
        epilog="Contact with cheonghing.so@gmail.com"
    )

    # Add the "version" argument to the microcat
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        default=False,
        help="print software version and exit",
    )
    # Create a sub-parser object for the common arguments
    common_parser = argparse.ArgumentParser(add_help=False)
    ## Add the "workdir" argument to the common sub-parser
    common_parser.add_argument(
        "-d",
        "--workdir",
        metavar="WORKDIR",
        type=str,
        default="./",
        help="project workdir",
    )
    ## Add the "check_samples" argument to the common sub-parser
    common_parser.add_argument(
        "--check-samples",
        dest="check_samples",
        default=False,
        action="store_true",
        help="check samples, default: False",
    )
    ###############################################################
    # Create a sub-parser object for the "run" command
    run_parser = argparse.ArgumentParser(add_help=False)
    ## Add the "config" argument to the "run" sub-parser
    run_parser.add_argument(
        "--config",
        type=str,
        default="./config.yaml",
        help="config.yaml",
    )
    ## run local (no cluster)
    ### Add the "run-local" argument to the "run" sub-parser
    run_parser.add_argument(
        "--run-local",
        default=False,
        dest="run_local",
        action="store_true",
        help="run pipeline on local computer",
    )
    #### Add the "cores" argument to the "run" sub-parser
    run_parser.add_argument(
        "--cores",
        type=int,
        default=60,
        help="all job cores, available on '--run-local'"
    )

    ##############################################################
    ## run remote cluster
    ## More detail in https://snakemake.readthedocs.io/en/stable/executing/cluster.html
    ### Add the "run-remote" argument to the "run" sub-parser
    run_parser.add_argument(
        "--run-remote",
        default=False,
        dest="run_remote",
        action="store_true",
        help="run pipeline on remote cluster",
    )
    ### Add the "local-cores" argument to the "run" sub-parser
    run_parser.add_argument(
        "--local-cores",
        type=int,
        dest="local_cores",
        default=8,
        help="local job cores, available on '--run-remote'"
    )
    ### Add the "jobs" argument to the "run" sub-parser
    run_parser.add_argument(
        "--jobs",
        type=int,
        default=30,
        help="cluster job numbers, available on '--run-remote'"
    )
    ## Add the "cluster-engine" argument to the "run" sub-parser
    ### TODO: support sge(qsub) and slurm(sbatch)
    run_parser.add_argument(
        "--cluster-engine",
        default="bsub",
        choices=["slurm", "sge", "lsf"],
        help="cluster workflow manager engine, now only support lsf(bsub)"
    )
    ## Add the "list" argument to the "run" sub-parser
    run_parser.add_argument(
        "--list",
        default=False,
        action="store_true",
        help="list pipeline rules",
    )
    ## Add the "debug" argument to the "run" sub-parser
    run_parser.add_argument(
        "--debug",
        default=False,
        action="store_true",
        help="debug pipeline",
    )
    ## Add the "dry-run" argument to the "run" sub-parser
    run_parser.add_argument(
        "--dry-run",
        default=False,
        dest="dry_run",
        action="store_true",
        help="dry run pipeline",
    )
    ## Add the "wait" argument to the "run" sub-parser
    run_parser.add_argument("--wait", type=int, default=60, help="wait given seconds")

    ############################################################################
    ## Add the conda related argument to the "run" sub-parser
    ## When --use-conda is activated, Snakemake will automatically create 
    ## software environments for any used wrapper
    ## Moer detail in https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management
    ### Add the "use-conda" argument to the "run" sub-parser
    run_parser.add_argument(
        "--use-conda",
        default=False,
        dest="use_conda",
        action="store_true",
        help="use conda environment",
    )
    ### 
    run_parser.add_argument(
        "--conda-prefix",
        default="~/.conda/envs",
        dest="conda_prefix",
        help="conda environment prefix",
    )
    run_parser.add_argument(
        "--conda-create-envs-only",
        default=False,
        dest="conda_create_envs_only",
        action="store_true",
        help="conda create environments only",
    )
    ################################################################
    ## Running jobs in containers
    run_parser.add_argument(
        "--use-singularity",
        default=False,
        dest="use_singularity",
        action="store_true",
        help="use a singularity container",
    )
    run_parser.add_argument(
        "--singularity-prefix",
        default="",
        dest="singularity_prefix",
        help="singularity images prefix",
    )
    ###############################################################
    # Create a sub-parser object
    subparsers = parser.add_subparsers(title="available subcommands", metavar="")
    ##  
    parser_init = subparsers.add_parser(
        "init",
        formatter_class=configer.custom_help_formatter,
        parents=[common_parser],# add common parser
        prog="microcat init",
        help="init project",
    )
    parser_init.add_argument(
            "-s",
            "--samples",
            type=str,
            default=None,
            help="""desired input:
    samples list, tsv format required.
    """,
        )
    ## add init project contain work 
    parser_init.add_argument(
        "-b",
        "--begin",
        type=str,
        default="trimming",
        choices=["simulate", "trimming", "host", "classifier", "denosing"],
        help="pipeline starting point",
    )

    parser_init.add_argument(
        "--trimmer",
        type=str,
        default="fastp",
        required=False,
        choices=["sickle", "fastp", "trimmomatic"],
        help="which trimmer used",
    )
    parser_init.add_argument(
        "--host",
        type=str,
        default="starsolo",
        required=False,
        choices=["starsolo","cellranger"],
        help="which rmhoster used",
    )
    parser_init.add_argument(
        "--assay",
        type=str,
        default=None,
        choices=["Smartseq", "Smartseq2", "tenX_v3", "tenX_v2", "tenX_v1"],
        help="Sequencing assay option, required when host is starsolo",
    )
    
    parser_init.add_argument(
        "--classifier",
        nargs="+",
        required=False,
        default="kraken2uniq",
        choices=["kraken2uniq","krakenuniq","pathseq"],
        help="wchich classifier used",
    )

    parser_init.set_defaults(func=init)
    


    parser_bulk_wf = subparsers.add_parser(
            "bulk_wf",
            formatter_class=configer.custom_help_formatter,
            parents=[common_parser, run_parser],
            prog="smart bulk_wf",
            help="bulk rna seq microbiome mining pipeline",
        )
    # parser_bulk_wf.add_argument(
    #     "task",
    #     metavar="TASK",
    #     nargs="?",
    #     type=str,
    #     default="all",
    #     choices=[WORKFLOWS_MAG],
    #     help="pipeline end point. Allowed values are " + ", ".join(WORKFLOWS_MAG),
    # )
    parser_bulk_wf.set_defaults(func=bulk_wf)

    parser_scrna_wf = subparsers.add_parser(
            "scrna_wf",
            formatter_class=configer.custom_help_formatter,
            parents=[common_parser, run_parser],
            prog="smart scrna_wf",
            help="single cell rna seq microbiome mining pipeline",
        )

    parser_scrna_wf.add_argument(
            "--platform",
            type=str,
            default="10x",
            choices=["10x", "smart-seq2"],
            help="single cell sequencing platform,support 10x and smart-seq2",
        )
    # parser_scrna_wf.add_argument(
    #     "task",
    #     metavar="TASK",
    #     nargs="?",
    #     type=str,
    #     default="all",
    #     choices=WORKFLOWS_MAG,
    #     help="pipeline end point. Allowed values are " + ", ".join(WORKFLOWS_MAG),
    # )
    parser_scrna_wf.set_defaults(func=scRNA_wf)

    args, unknown = parser.parse_known_args()
        
    if hasattr(args, 'host') and args.host == "starsolo" and args.assay is None:
        parser_init.error("--assay option is required when host is starsolo")

    try:
        if args.version:
            print("microcat version %s" % __version__)
            sys.exit(0)
        args.func(args, unknown)
    except AttributeError as e:
        print(e)
        parser.print_help()

if __name__ == "__main__":
    main()