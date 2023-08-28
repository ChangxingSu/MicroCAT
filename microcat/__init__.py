#!/usr/bin/env python
from microcat.configer import MicrocatConfig
from microcat.configer import parse_yaml
from microcat.configer import update_config
# from microcat.configer import custom_help_formatter

from microcat.sample import HEADERS
from microcat.sample import parse_samples
from microcat.sample import parse_bam_samples


from microcat.sample import get_starsolo_sample_id
from microcat.sample import get_sample_id
from microcat.sample import get_fastqs_dir
from microcat.sample import get_samples_bax

from microcat.sample import get_samples_id_by_tissue
from microcat.sample import get_samples_id_by_patient
from microcat.sample import get_samples_id_by_lane
from microcat.sample import get_samples_id_by_library
from microcat.sample import get_tissue_by_patient
from microcat.sample import get_SAMattrRGline_by_sample
from microcat.sample import get_SAMattrRGline_from_manifest

from microcat.__about__ import __version__, __author__


name = "microcat"
