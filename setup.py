#!/usr/bin/env python
# coding=utf-8
from setuptools import setup

# read version number
exec(open("uvplot/_version.py", "r").read())

# use README as long_description (for PyPI)
try:
    long_description = open("README.rst", "r", encoding="utf-8").read()
except TypeError:
    # under Python 2.7 the `encoding` keyword doesn't exist.
    print(
        "DEPRECATION: Python 2.7 will reach the end of its life on January 1st, 2020. "
        "Please upgrade your Python as Python 2.7 won't be maintained after that date. "
        "A future version of uvplot will drop support for Python 2.7. "
    )
    long_description = open("README.rst", "r").read()


setup(
    name="uvplot",
    version=__version__,
    packages=["uvplot"],
    author="Marco Tazzari",
    author_email="mtazzari@gmail.com",
    description="Utilities for handling and plotting interferometric visibilities.",
    long_description=long_description,
    install_requires=["numpy>=1.9", "matplotlib"],
    license="LGPLv3",
    url="https://github.com/mtazzari/uvplot",
    download_url="https://github.com/mtazzari/uvplot/archive/{}.tar.gz".format(
        __version__
    ),
    keywords=["science", "astronomy", "interferometry"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
    ],
)
