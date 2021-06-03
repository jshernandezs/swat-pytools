# SWAT-pytools

A Python wrapper for executing and calibrating the [Soil and Water Assessment Tool (SWAT)](https://swat.tamu.edu/) in Unix/macOS systems. A module for computing the 171 hydrologic indices reported by Henriksen et al. (2006) and the Magnificent Seven indices proposed by Archfield et al. (2014) is also included here. For single- and multi-objective calibration, we use the [pymoo](https://pymoo.org/) framework.

## Installation

First, you need to have a Python 3 enviroment installed. We recommend using [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [anaconda](https://www.anaconda.com/). Check if conda is available in the command line:

 ```bash
 conda --version
 ```
 
 Create a new Python environment with Numpy, Pandas, and Matplotlib preinstalled and activate it:
 
 ```bash
 conda create -n swatpy -y python==3.7 numpy pandas matplotlib
 conda activate swatpy
 ```
 
 If you are familiar with conda, you might prefer to use an existing environment that you already created.
 
 Now, you can install the current SWAT-pytools version. `cd` to the directory where you normally save your repositories and execute these lines in a terminal:
 
 ```bash
 git clone https://github.com/jshernandezs/swat-pytools
 cd swat-pytools
 pip install .
 ```

## Getting started

Once you download this package, you will find the following local directory tree:

<pre>
swat-pytools
├── resources
│   ├── csv_files
│   ├── figFiles
│   ├── Models
│   └── Observed
├── src
│   └── swat_utilities
└── tests

</pre>

In the `resources` directory, you must place all the necessary files for executing and calibrating SWAT, such as the model executables, zip files containing input SWAT text files (`resources/Models` directory), csv files containing observed time series (`resources/Observed` directory), and other optional files. Regarding SWAT executables, revisions 622, 670, and 681 are available for Unix in this repository, only version 622 is available for macOS.

### Executing SWAT with Python

In this example, we are using the SWAT 622 version to run a model of the Honeyoey Creek - Pine Creek Watershed located in Michigan, US. The SWAT input text files, which normally are generated inside the `TxtInOut` folder when using ArcSWAT, are put together in a zip file that will be handled by the Python wrapper. In this case, we are using the `Honeyoy_Model.zip` file placed in the `resources/Observed` directory.

We assume that a new Python script is created in the `test` directory. 

First, we import the libraries that we are going to use:

```python
import os
from swat_utilities.swat_config import ModelSetup
```
Then, we define a variable with the path to the zip file containing the SWAT text files:

```python
model_file_path = os.path.abspath('../resources/Models/Honeyoy_Model.zip')
```
Now, we create the model object using the `ModelSetup` class:

```python
swat_model = ModelSetup(model_file_path)
```

We must indicate which SWAT version we are going to use, which is the property `swat_exec_name` of the model object created above:

```python
swat_model.swat_exec_name = 'SWAT_Rev622'
```
To see the execution progress when running the model, we set the property 'verbose_swat' as `True`:

```python
swat_model.verbose_swat = True
```
To run the model, we need to execute the ```prepare_swat()``` method followed by the ```run_swat()``` method as follows:

```python
swat_model.prepare_swat()
swat_model.run_swat()
```
By default, the results are stored in the `/tmp/output_swat/New_SWAT` directory. The user can modify the output directory and the model folder containing the input and output SWAT text files by using the ```output_dir``` and ```new_model_name``` properties of the `swat_model` object, respectively.
