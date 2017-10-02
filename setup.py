#!/usr/bin/env python
# coding=utf-8
from setuptools import setup, find_packages


def requirements_file_to_list(fn="requirements.txt"):
    """read a requirements file and create a list that can be used in setup.

    """
    with open(fn, 'r') as f:
        return [x.rstrip() for x in list(f) if x and not x.startswith('#')]


setup(
    name="uvplot",
    version="0.1.0",
    packages=find_packages(),
    install_requires=requirements_file_to_list(),
    author="Marco Tazzari",
    author_email="mtazzari@ast.cam.ac.uk",
    description="Utilities for plotting interferometric visibilities.",
    long_description=open('README.rst').read(),
    package_data={u'': [u'LICENSE']},
    license="GPLv3",
    url="tbd",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
    ]
)