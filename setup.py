from __future__ import absolute_import
from __future__ import print_function

import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup


def read(*names, **kwargs):
    with io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ) as fh:
        return fh.read()

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="swat-pytools",
    version="0.0.1",
    author="J. Sebastian Hernandez-Suarez",
    author_email="herna505@msu.edu",
    description="A Python wrapper for executing and calibrating SWAT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jshernandezs/swat-pytools",
    project_urls={
        "Bug Tracker": "https://github.com/jshernandezs/swat-pytools/issues",
    },
    install_requires=['autograd>=1.3', 'hvwfg>=1.0.2', 'matplotlib>=3.4',
    'numpy>=1.15', 'pandas>=1.1', 'pymoo==0.5.0', 'scipy>=1.1', 'dask>=2021.10',
    'dask-jobqueue>=0.7'],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: GNU General Public License",
        "Operating System :: Unix/macOS",
    ],
    package_dir={"": "src"},
    packages=find_packages("src"),
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    python_requires=">=3.7",
    entry_points={
        'console_scripts': [
            'api = api.api:main',
        ]
    },
)
