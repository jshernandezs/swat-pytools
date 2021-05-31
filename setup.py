import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
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
    'numpy>=1.15', 'pandas>=1.1', 'pymoo>=0.4.2', 'scipy>=1.1'],
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
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.7",
)
