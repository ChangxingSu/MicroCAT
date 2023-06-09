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
        "microcat/config/*.yaml",
        "microcat/envs/*.yaml",
        "microcat/snakefiles/*.smk",
        "microcat/rules/*.smk",
        "microcat/scripts/*.py",
        "microcat/scripts/*.R",
        "microcat/*.py",
    ]
}

data_files = [(".", ["LICENSE", "README.md"])]

entry_points = {"console_scripts": ["microcat=microcat.corer:main"]}

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

