#!/usr/bin/env python
# coding=utf-8
from setuptools import setup, find_packages


setup(
    name="uvplot",
    version="0.2.1",
    packages=find_packages(),
    author="Marco Tazzari",
    author_email="mtazzari@ast.cam.ac.uk",
    description="Utilities for plotting interferometric visibilities.",
    long_description=open('README.rst').read(),
    install_requires=["numpy>=1.9", "matplotlib"],
    license="LGPLv3",
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