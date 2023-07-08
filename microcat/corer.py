#!/usr/bin/env python

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import textwrap
from io import StringIO
import pandas as pd

import microcat
import re
import sys
from typing import TYPE_CHECKING, Callable, ClassVar, Iterable, Iterator

# rich is only used to display help. It is imported inside the functions in order
# not to add delays to command line tools that use this formatter.
if TYPE_CHECKING:
    from argparse import Action, _ArgumentGroup
    from typing_extensions import Self

    from rich.console import Console, ConsoleOptions, RenderableType, RenderResult
    from rich.containers import Lines
    from rich.style import StyleType
    from rich.text import Span, Text

__all__ = [
    "RichHelpFormatter",
    "RawDescriptionRichHelpFormatter",
    "RawTextRichHelpFormatter",
    "ArgumentDefaultsRichHelpFormatter",
    "MetavarTypeRichHelpFormatter",
]


class RichHelpFormatter(argparse.HelpFormatter):
    """An argparse HelpFormatter class that renders using rich."""

    group_name_formatter: ClassVar[Callable[[str], str]] = str.title
    """A function that formats group names. Defaults to ``str.title``."""
    styles: ClassVar[dict[str, StyleType]] = {
        "argparse.args": "dark_cyan",
        "argparse.groups": "yellow",
        "argparse.help": "default",
        "argparse.metavar": "cyan",
        "argparse.syntax": "bold",
        "argparse.text": "default",
        "argparse.prog": "grey50",
    }
    """A dict of rich styles to control the formatter styles.

    The following styles are used:

    - ``argparse.args``: for positional-arguments and --options (e.g "--help")
    - ``argparse.groups``: for group names (e.g. "positional arguments")
    - ``argparse.help``: for argument's help text (e.g. "show this help message and exit")
    - ``argparse.metavar``: for meta variables (e.g. "FILE" in "--file FILE")
    - ``argparse.prog``: for %(prog)s in the usage (e.g. "foo" in "Usage: foo [options]")
    - ``argparse.syntax``: for highlights of back-tick quoted text (e.g. "``` `some text` ```"),
    - ``argparse.text``: for the descriptions and epilog (e.g. "A foo program")
    """
    highlights: ClassVar[list[str]] = [
        r"(?:^|\s)(?P<args>-{1,2}[\w]+[\w-]*)",  # highlight --words-with-dashes as args
        r"`(?P<syntax>[^`]*)`",  # highlight `text in backquotes` as syntax
    ]
    """A list of regex patterns to highlight in the help text.

    It is used in the description, epilog, groups descriptions, and arguments' help. By default,
    it highlights ``--words-with-dashes`` with the `argparse.args` style and
    ``` `text in backquotes` ``` with the `argparse.syntax` style.

    To disable highlighting, clear this list (``RichHelpFormatter.highlights.clear()``).
    """
    usage_markup: ClassVar[bool] = False
    """If True, render the usage string passed to ``ArgumentParser(usage=...)`` as markup.

    Defaults to ``False`` meaning the text of the usage will be printed verbatim.

    Note that the auto-generated usage string is always colored.
    """

    _root_section: _Section
    _current_section: _Section

    def __init__(
        self,
        prog: str,
        indent_increment: int = 2,
        max_help_position: int = 24,
        width: int | None = None,
    ) -> None:
        super().__init__(prog, indent_increment, max_help_position, width)
        self._console: Console | None = None

        # https://docs.python.org/3/library/stdtypes.html#printf-style-string-formatting
        self._printf_style_pattern = re.compile(
            r"""
            %                               # Percent character
            (?:\((?P<mapping>[^)]*)\))?     # Mapping key
            (?P<flag>[#0\-+ ])?             # Conversion Flags
            (?P<width>\*|\d+)?              # Minimum field width
            (?P<precision>\.(?:\*?|\d*))?   # Precision
            [hlL]?                          # Length modifier (ignored)
            (?P<format>[diouxXeEfFgGcrsa%]) # Conversion type
            """,
            re.VERBOSE,
        )

    @property
    def console(self) -> Console:  # deprecate?
        if self._console is None:
            from rich.console import Console
            from rich.theme import Theme

            self._console = Console(theme=Theme(self.styles))
        return self._console

    @console.setter
    def console(self, console: Console) -> None:  # is this needed?
        self._console = console

    class _Section(argparse.HelpFormatter._Section):  # type: ignore[misc]
        def __init__(
            self, formatter: RichHelpFormatter, parent: Self | None, heading: str | None = None
        ) -> None:
            if heading is not None:
                heading = f"{type(formatter).group_name_formatter(heading)}:"
            super().__init__(formatter, parent, heading)
            self.rich_items: list[RenderableType] = []
            self.rich_actions: list[tuple[Text, Text | None]] = []
            if parent is not None:
                parent.rich_items.append(self)
            if TYPE_CHECKING:  # already assigned in super().__init__ but not present in typeshed
                self.formatter = formatter
                self.heading = heading
                self.parent = parent

        def __rich_console__(self, console: Console, options: ConsoleOptions) -> RenderResult:
            from rich.text import Text

            # empty section
            if not self.rich_items and not self.rich_actions:
                return
            # root section
            if self is self.formatter._root_section:
                yield from self.rich_items
                return
            # group section
            help_pos = min(self.formatter._action_max_length + 2, self.formatter._max_help_position)
            help_width = max(self.formatter._width - help_pos, 11)
            if self.heading:
                yield Text(self.heading, style="argparse.groups")
            yield from self.rich_items  # (optional) group description
            indent = Text(" " * help_pos)
            for action_header, action_help in self.rich_actions:
                if not action_help:
                    yield action_header  # no help, yield the header and finish
                    continue
                action_help_lines = self.formatter._rich_split_lines(action_help, help_width)
                if len(action_header) > help_pos - 2:
                    yield action_header  # the header is too long, put it on its own line
                    action_header = indent
                action_header.set_length(help_pos)
                action_help_lines[0].rstrip()
                yield action_header + action_help_lines[0]
                for line in action_help_lines[1:]:
                    line.rstrip()
                    yield indent + line
            yield "\n"

    def add_text(self, text: str | None) -> None:
        if text is not argparse.SUPPRESS and text is not None:
            self._current_section.rich_items.append(self._rich_format_text(text))

    def add_usage(
        self,
        usage: str | None,
        actions: Iterable[Action],
        groups: Iterable[_ArgumentGroup],
        prefix: str | None = None,
    ) -> None:
        if usage is argparse.SUPPRESS:
            return
        from rich.text import Span, Text

        if prefix is None:
            prefix = self._format_usage(usage="", actions=(), groups=(), prefix=None).rstrip("\n")
        prefix_end = ": " if prefix.endswith(": ") else ""
        prefix = prefix[: len(prefix) - len(prefix_end)]
        prefix = type(self).group_name_formatter(prefix) + prefix_end

        usage_spans = [Span(0, len(prefix.rstrip()), "argparse.groups")]
        usage_text = self._format_usage(usage, actions, groups, prefix=prefix)
        if usage is None:  # get colour spans for generated usage
            prog = f"{self._prog}"
            if actions:
                prog_start = usage_text.index(prog, len(prefix))
                usage_spans.append(Span(prog_start, prog_start + len(prog), "argparse.prog"))
            actions_start = len(prefix) + len(prog) + 1
            try:
                spans = list(self._rich_usage_spans(usage_text, actions_start, actions=actions))
            except ValueError:
                spans = []
            usage_spans.extend(spans)
            rich_usage = Text(usage_text)
        elif self.usage_markup:  # treat user provided usage as markup
            usage_spans.extend(self._rich_prog_spans(prefix + Text.from_markup(usage).plain))
            rich_usage = Text.from_markup(usage_text)
            usage_spans.extend(rich_usage.spans)
            rich_usage.spans.clear()
        else:  # treat user provided usage as plain text
            usage_spans.extend(self._rich_prog_spans(prefix + usage))
            rich_usage = Text(usage_text)
        rich_usage.spans.extend(usage_spans)
        self._root_section.rich_items.append(rich_usage)

    def add_argument(self, action: Action) -> None:
        super().add_argument(action)
        if action.help is not argparse.SUPPRESS:
            self._current_section.rich_actions.extend(self._rich_format_action(action))

    def format_help(self) -> str:
        with self.console.capture() as capture:
            self.console.print(self._root_section, highlight=False, soft_wrap=True)
        help = capture.get()
        if help:
            help = self._long_break_matcher.sub("\n\n", help).rstrip() + "\n"
        return help

    # ===============
    # Utility methods
    # ===============
    def _rich_prog_spans(self, usage: str) -> Iterator[Span]:
        from rich.text import Span

        if "%(prog)" not in usage:
            return
        params = {"prog": self._prog}
        formatted_usage = ""
        last = 0
        for m in self._printf_style_pattern.finditer(usage):
            start, end = m.span()
            formatted_usage += usage[last:start]
            sub = usage[start:end] % params
            prog_start = len(formatted_usage)
            prog_end = prog_start + len(sub)
            formatted_usage += sub
            last = end
            yield Span(prog_start, prog_end, "argparse.prog")

    def _rich_usage_spans(self, text: str, start: int, actions: Iterable[Action]) -> Iterator[Span]:
        from rich.text import Span

        options: list[Action] = []
        positionals: list[Action] = []
        for action in actions:
            if action.help is not argparse.SUPPRESS:
                options.append(action) if action.option_strings else positionals.append(action)
        pos = start
        for action in options:  # start with the options
            if sys.version_info >= (3, 9):  # pragma: >=3.9 cover
                usage = action.format_usage()
                if isinstance(action, argparse.BooleanOptionalAction):
                    for option_string in action.option_strings:
                        start = text.index(option_string, pos)
                        end = start + len(option_string)
                        yield Span(start, end, "argparse.args")
                        pos = end + 1
                    continue
            else:  # pragma: <3.9 cover
                usage = action.option_strings[0]
            start = text.index(usage, pos)
            end = start + len(usage)
            yield Span(start, end, "argparse.args")
            if action.nargs != 0:
                metavar = self._format_args(action, self._get_default_metavar_for_optional(action))
                start = text.index(metavar, end)
                end = start + len(metavar)
                yield Span(start, end, "argparse.metavar")
            pos = end + 1
        for action in positionals:  # positionals come at the end
            usage = self._format_args(action, self._get_default_metavar_for_positional(action))
            start = text.index(usage, pos)
            end = start + len(usage)
            yield Span(start, end, "argparse.args")
            pos = end + 1

    def _rich_whitespace_sub(self, text: Text) -> Text:
        # do this `self._whitespace_matcher.sub(' ', text).strip()` but text is Text
        spans = [m.span() for m in self._whitespace_matcher.finditer(text.plain)]
        for start, end in reversed(spans):
            if end - start > 1:  # slow path
                space = text[start : start + 1]
                space.plain = " "
                text = text[:start] + space + text[end:]
            else:  # performance shortcut
                text.plain = text.plain[:start] + " " + text.plain[end:]
        # Text has no strip method yet
        lstrip_at = len(text.plain) - len(text.plain.lstrip())
        if lstrip_at:
            text = text[lstrip_at:]
        text.rstrip()
        return text

    # =====================================
    # Rich version of HelpFormatter methods
    # =====================================
    def _rich_expand_help(self, action: Action) -> Text:
        from rich.markup import escape
        from rich.text import Text

        params = dict(vars(action), prog=self._prog)
        for name in list(params):
            if params[name] is argparse.SUPPRESS:
                del params[name]
            elif hasattr(params[name], "__name__"):
                params[name] = params[name].__name__
        if params.get("choices") is not None:
            params["choices"] = ", ".join([str(c) for c in params["choices"]])
        help_string = self._get_help_string(action)
        assert help_string is not None
        help_string % params  # pyright: ignore[reportUnusedExpression] # raise ValueError if needed
        parts = []
        last = 0
        for m in self._printf_style_pattern.finditer(help_string):
            start, end = m.span()
            parts.append(help_string[last:start])
            parts.append(escape(help_string[start:end] % params))
            last = end
        parts.append(help_string[last:])
        rich_help = Text.from_markup("".join(parts), style="argparse.help")
        for highlight in self.highlights:
            rich_help.highlight_regex(highlight, style_prefix="argparse.")
        return rich_help

    def _rich_format_text(self, text: str) -> Text:
        from rich.markup import escape
        from rich.text import Text

        if "%(prog)" in text:
            text = text % {"prog": escape(self._prog)}
        rich_text = Text.from_markup(text, style="argparse.text")
        for highlight in self.highlights:
            rich_text.highlight_regex(highlight, style_prefix="argparse.")
        text_width = max(self._width - self._current_indent * 2, 11)
        indent = Text(" " * self._current_indent)
        return self._rich_fill_text(rich_text, text_width, indent)

    def _rich_format_action(self, action: Action) -> Iterator[tuple[Text, Text | None]]:
        header = self._rich_format_action_invocation(action)
        header.pad_left(self._current_indent)
        help = self._rich_expand_help(action) if action.help and action.help.strip() else None
        yield header, help
        for subaction in self._iter_indented_subactions(action):
            yield from self._rich_format_action(subaction)

    def _rich_format_action_invocation(self, action: Action) -> Text:
        from rich.text import Text

        if not action.option_strings:
            return Text().append(self._format_action_invocation(action), style="argparse.args")
        else:
            action_header = Text(", ").join(Text(o, "argparse.args") for o in action.option_strings)
            if action.nargs != 0:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                action_header.append_tokens(((" ", None), (args_string, "argparse.metavar")))
            return action_header

    def _rich_split_lines(self, text: Text, width: int) -> Lines:
        return self._rich_whitespace_sub(text).wrap(self.console, width)

    def _rich_fill_text(self, text: Text, width: int, indent: Text) -> Text:
        lines = self._rich_whitespace_sub(text).wrap(self.console, width)
        return type(text)("\n").join(indent + line for line in lines) + "\n\n"


class RawDescriptionRichHelpFormatter(RichHelpFormatter):
    """Rich help message formatter which retains any formatting in descriptions."""

    def _rich_fill_text(self, text: Text, width: int, indent: Text) -> Text:
        return type(text)("\n").join(indent + line for line in text.split()) + "\n\n"


class RawTextRichHelpFormatter(RawDescriptionRichHelpFormatter):
    """Rich help message formatter which retains formatting of all help text."""

    def _rich_split_lines(self, text: Text, width: int) -> Lines:
        return text.split()


class ArgumentDefaultsRichHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, RichHelpFormatter):
    """Rich help message formatter which adds default values to argument help."""


class MetavarTypeRichHelpFormatter(argparse.MetavarTypeHelpFormatter, RichHelpFormatter):
    """Rich help message formatter which uses the argument 'type' as the default
    metavar value (instead of the argument 'dest').
    """

def add_parser(*args, **kwds):
    kwds.setdefault("formatter_class", parser.formatter_class)
    return subparsers.add_parser(*args, **kwds)

WORKFLOWS_MAG = [
    "host_all",
    "check_all",]

WORKFLOWS_SCRNA = [
    "host_all",
    "kraken2uniq_classified_all",
    "krakenuniq_classified_all",
    "pathseq_classified_all",
    "metaphlan_classified_all",
    "classifier_all",]

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

    for classifier_ in ["kraken2uniq","krakenuniq","pathseq","metaphlan"]:
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
        project = microcat.MicrocatConfig(args.workdir)


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


        for script_path in conf["scripts"]:
            origin_path = conf["scripts"][script_path]
            conf["scripts"][script_path] = os.path.join(os.path.dirname(__file__),f"{origin_path}")

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
        microcat.update_config(
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
    conf = microcat.parse_yaml(args.config)

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
        "--until",
        args.task
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

def prepare(args, unknown):

    script_path = os.path.join(os.path.dirname(__file__),"prepare.py")  
    command = ["python", script_path]
    command.extend(unknown)  

    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError as e:
        print(f"Error: Failed to run prepare.py ({e})")


def main():

    # Banner text and program information
    banner = '''
    Microbiome Identification in Cell Resolution from Omics-Computational Analysis Toolbox
    '''

    # Create the main parser object with a banner and program name
    parser = argparse.ArgumentParser(
        formatter_class=RichHelpFormatter,
        description=textwrap.dedent(banner),
        prog="microcat",
        epilog="Contact with changxingsu42@gmail.com"
    )

    # Add the "version" argument to the microcat
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        default=False,
        help="Print software version and exit",
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
        help="Project workdir",
    )
    ## Add the "check_samples" argument to the common sub-parser
    common_parser.add_argument(
        "--check-samples",
        dest="check_samples",
        default=False,
        action="store_true",
        help="Check samples, default: False",
    )
    ###############################################################
    # Create a sub-parser object for the "run" command
    run_parser = argparse.ArgumentParser(add_help=False)
    ## Add the "config" argument to the "run" sub-parser
    run_parser.add_argument(
        "--config",
        type=str,
        default="./config.yaml",
        help="Path of config.yaml",
    )
    ## run local (no cluster)
    ### Add the "run-local" argument to the "run" sub-parser
    run_parser.add_argument(
        "--run-local",
        default=False,
        dest="run_local",
        action="store_true",
        help="Run pipeline on local computer",
    )
    #### Add the "cores" argument to the "run" sub-parser
    run_parser.add_argument(
        "--cores",
        type=int,
        default=60,
        help="All job cores, available on '--run-local'"
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
        help="Run pipeline on remote cluster",
    )
    ### Add the "local-cores" argument to the "run" sub-parser
    run_parser.add_argument(
        "--local-cores",
        type=int,
        dest="local_cores",
        default=8,
        help="Local job cores, available on '--run-remote'"
    )
    ### Add the "jobs" argument to the "run" sub-parser
    run_parser.add_argument(
        "--jobs",
        type=int,
        default=30,
        help="Cluster job numbers, available on '--run-remote'"
    )
    ## Add the "cluster-engine" argument to the "run" sub-parser
    ### TODO: support sge(qsub) and slurm(sbatch)
    run_parser.add_argument(
        "--cluster-engine",
        default="bsub",
        choices=["slurm", "sge", "lsf"],
        help="Cluster workflow manager engine, now only support lsf(bsub)"
    )
    ## Add the "list" argument to the "run" sub-parser
    run_parser.add_argument(
        "--list",
        default=False,
        action="store_true",
        help="List pipeline rules",
    )
    ## Add the "debug" argument to the "run" sub-parser
    run_parser.add_argument(
        "--debug",
        default=False,
        action="store_true",
        help="Debug pipeline",
    )
    ## Add the "dry-run" argument to the "run" sub-parser
    run_parser.add_argument(
        "--dry-run",
        default=False,
        dest="dry_run",
        action="store_true",
        help="Dry run pipeline",
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
        help="Use conda environment",
    )
    ### 
    run_parser.add_argument(
        "--conda-prefix",
        default="~/.conda/envs",
        dest="conda_prefix",
        help="Conda environment prefix",
    )
    run_parser.add_argument(
        "--conda-create-envs-only",
        default=False,
        dest="conda_create_envs_only",
        action="store_true",
        help="Conda create environments only",
    )
    ################################################################
    ## Running jobs in containers
    run_parser.add_argument(
        "--use-singularity",
        default=False,
        dest="use_singularity",
        action="store_true",
        help="Use a singularity container",
    )
    run_parser.add_argument(
        "--singularity-prefix",
        default="",
        dest="singularity_prefix",
        help="Singularity images prefix",
    )
    ###############################################################
    # Create a sub-parser object
    subparsers = parser.add_subparsers(title="available subcommands", metavar="")
    ##  
    parser_init = subparsers.add_parser(
        "init",
        formatter_class=RichHelpFormatter,
        parents=[common_parser],# add common parser
        prog="microcat init",
        help="Init project",
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
        help="Pipeline starting point",
    )

    parser_init.add_argument(
        "--trimmer",
        type=str,
        default="fastp",
        required=False,
        choices=["sickle", "fastp", "trimmomatic"],
        help="Which trimmer used",
    )
    parser_init.add_argument(
        "--host",
        type=str,
        default="starsolo",
        required=False,
        choices=["starsolo","cellranger"],
        help="Which rmhoster used",
    )
    parser_init.add_argument(
        "--assay",
        type=str,
        default=None,
        choices=["Smartseq", "Smartseq2", "tenX_v3", "tenX_v2", "tenX_v1","Seq-well","tenX_auto"],
        help="Sequencing assay option, required when host is starsolo",
    )
    
    parser_init.add_argument(
        "--classifier",
        nargs="+",
        required=False,
        default="kraken2uniq",
        choices=["kraken2uniq","krakenuniq","pathseq","metaphlan"],
        help="Which classifier used",
    )

    parser_init.set_defaults(func=init)
    


    parser_bulk_wf = subparsers.add_parser(
            "bulk_wf",
            formatter_class=RichHelpFormatter,
            parents=[common_parser, run_parser],
            prog="microcat bulk_wf",
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
            formatter_class=RichHelpFormatter,
            parents=[common_parser, run_parser],
            prog="microcat scrna_wf",
            help="single cell rna seq microbiome mining pipeline",
        )

    parser_prepare = subparsers.add_parser(
            "prepare",
            formatter_class=RichHelpFormatter,
            parents=[common_parser],
            prog="microcat prepare",
            help="Prepare for microcat",
        )
    parser_prepare.set_defaults(func=prepare)

    # parser_scrna_wf.add_argument(
    #         "--platform",
    #         type=str,
    #         default="10x",
    #         choices=["10x", "smart-seq2"],
    #         help="single cell sequencing platform,support 10x and smart-seq2",
    #     )
    parser_scrna_wf.add_argument(
        "task",
        metavar="TASK",
        nargs="?",
        type=str,
        default="all",
        choices=WORKFLOWS_SCRNA,
        help="pipeline end point. Allowed values are " + ", ".join(WORKFLOWS_MAG),
    )
    parser_scrna_wf.set_defaults(func=scRNA_wf)

    args, unknown = parser.parse_known_args()
        
    if hasattr(args, 'host') and args.host == "starsolo" and args.assay is None:
        parser_init.error("--assay option is required when host is starsolo")

    try:
        if args.version:
            print("microcat version %s" % microcat.__version__)
            sys.exit(0)
        args.func(args, unknown)
    except AttributeError as e:
        print(e)
        parser.print_help()

if __name__ == "__main__":
    main()