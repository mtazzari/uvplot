#!/usr/bin/env python
# coding=utf-8
from setuptools import setup

# read version number
exec(open("uvplot/_version.py").read())

setup(
    name="uvplot",
    version=__version__,
    packages=['uvplot'],
    author="Marco Tazzari",
    author_email="mtazzari@gmail.com",
    description="Utilities for handling and plotting interferometric visibilities.",
    long_description=open('README.rst').read(),
    install_requires=["numpy>=1.9", "matplotlib"],
    license="LGPLv3",
    url="https://github.com/mtazzari/uvplot",
    download_url='https://github.com/mtazzari/uvplot/archive/{}.tar.gz'.format(__version__),
    keywords = ['science', 'astronomy', 'interferometry'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
    ]
)
