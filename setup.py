#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import sys
from shutil import rmtree
from setuptools import find_packages, setup, Command

# read version and author information
exec(open("microcat/__about__.py").read())

# publish command
if sys.argv[-1] == "publish":
    os.system("python setup.py sdist bdist_wheel")
    os.system("twine upload dist/*")
    sys.exit()

# read README
try:
    with io.open(os.path.join('README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# base packages
packages = find_packages(include=['microcat', 'microcat.*'])

# include files
package_data = {
    "microcat": [
        # Workflow files
        "*/config/*.yaml",
        "*/envs/*.yaml",
        "*/snakefiles/*.smk",
        "*/rules/*.smk",
        "*/scripts/*.*",
        
        # Other files
        "profiles/*",
        "*.py",
    ]
}

# base dependencies
requires = [
    req.strip()
    for req in open("requirements.txt", "r").readlines()
    if not req.startswith("#")
]

# test dependencies
test_requires = [
    "pytest>=6.0",
    "pytest-cov>=2.0",
    "pytest-mock>=3.6",
]

lint_requires = [
    "black>=21.0",
    "flake8>=3.9",
    "isort>=5.0",
    "mypy>=0.910",
]

docs_requires = [
    "jupyter-book>=0.15.1",
    "jupytext>=1.14.5",
    "sphinx-copybutton>=0.5.2",
    "sphinx-design>=0.4.1",
    "sphinx-togglebutton>=0.3.2",
    "sphinx-thebe>=0.2.1",
    "sphinx-external-toc>=0.3.1",
    "sphinx-book-theme>=1.0.1",
    "sphinx-comments>=0.0.3",
    "sphinx-autodoc2>=0.4.2",
    "myst-nb>=0.17.2",
]

dev_requires = test_requires + lint_requires + docs_requires + [
    "pre-commit>=2.15",
    "bump2version>=1.0",
    "twine>=3.4",
    "build>=0.7",
]

setup(
    name="microcat",
    version=__version__,
    author=__author__,
    author_email="changxingsu42@gmail.com",
    url="https://github.com/ChangxingSu/MicroCAT",
    description="A computational toolbox to identify microbiome from Omics data",
    long_description_content_type="text/markdown",
    long_description=long_description,
    entry_points={"console_scripts": ["microcat=microcat.cli:microcat"]},
    packages=packages,
    package_data=package_data,
    data_files=[(".", ["LICENSE", "README.md"])],
    include_package_data=True,
    install_requires=requires,
    extras_require={
        "test": test_requires,
        "lint": lint_requires,
        "docs": docs_requires,
        "dev": dev_requires,
    },
    python_requires=">=3.9",
    license="GPLv3+",
    keywords="bioinformatics microbiome single-cell spatial-transcriptomics",
    project_urls={
        "Source": "https://github.com/ChangxingSu/MicroCAT",
        "Issue Tracker": "https://github.com/ChangxingSu/MicroCAT/issues",
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
