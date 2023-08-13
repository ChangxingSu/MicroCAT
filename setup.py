#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import sys
from shutil import rmtree

from setuptools import find_packages, setup, Command

exec(open("microcat/__about__.py").read())

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()


# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join('README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

packages = ["microcat"]

package_data = {
    "microcat": [
        "microcat/single_wf/config/*.yaml",
        "microcat/single_wf/envs/*.yaml",
        "microcat/single_wf/snakefiles/*.smk",
        "microcat/single_wf/rules/*.smk",
        "microcat/single_wf/scripts/*.py",
        "microcat/single_wf/scripts/*.R",
        "microcat/bulk_wf/config/*.yaml",
        "microcat/bulk_wf/envs/*.yaml",
        "microcat/bulk_wf/snakefiles/*.smk",
        "microcat/bulk_wf/rules/*.smk",
        "microcat/bulk_wf/scripts/*.py",
        "microcat/bulk_wf/scripts/*.R",
        "microcat/spatial_wf/config/*.yaml",
        "microcat/spatial_wf/envs/*.yaml",
        "microcat/spatial_wf/snakefiles/*.smk",
        "microcat/spatial_wf/rules/*.smk",
        "microcat/spatial_wf/scripts/*.py",
        "microcat/spatial_wf/scripts/*.R",
        "microcat/multi_wf/config/*.yaml",
        "microcat/multi_wf/envs/*.yaml",
        "microcat/multi_wf/snakefiles/*.smk",
        "microcat/multi_wf/rules/*.smk",
        "microcat/multi_wf/scripts/*.py",
        "microcat/multi_wf/scripts/*.R",
        "microcat/profiles/*",
        "microcat/*.py",
    ]
}

data_files = [(".", ["LICENSE", "README.md"])]

entry_points = {"console_scripts": ["microcat=microcat.cli:microcat"]}

requires = [
    req.strip()
    for req in open("requirements.txt", "r").readlines()
    if not req.startswith("#")
]

classifiers = [
    "Development Status :: 3 - Alpha",
    "Environment :: Console",
    "Intended Audience :: Developers",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Natural Language :: English",
    "Operating System :: OS Independent",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="microcat",
    version=__version__,
    author=__author__,
    author_email="changxingsu42@gmail.com",
    url="https://github.com/ChangxingSu/MicroCAT",
    description="a computational toolbox to identificated microbiome from Omics",
    long_description_content_type="text/markdown",
    long_description=long_description,
    entry_points=entry_points,
    packages=packages,
    package_data=package_data,
    data_files=data_files,
    include_package_data=True,
    install_requires=requires,
    license="GPLv3+",
    classifiers=classifiers,
)

