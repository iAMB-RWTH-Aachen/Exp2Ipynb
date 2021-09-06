# Exp2Ipynb

Last Version: 2021/06/17<br>
Author: Ulf Liebal<br>
Contact: ulf.liebal@rwth-aachen.de<br>
License: see LICENSE file<br>



The Exp2Ipynb workflow facilitates analysis of promoter libraries using Jupyter Notebooks. To use it, clone the repository.

## Installation
To use the Exp2Ipynb workflow, a Python environment is required. For Windows, Apple and Linus systems, Anaconda is a suitable Python program. To run Jupyter notebooks, the Jupyter environment needs to be available.

Requirements:
- [Python](https://www.anaconda.com/products/individual-d)
- [Jupyter](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html)
- More detailed requirements can be found in the [requirements.txt](./requirements.txt) file

Install Exp2Ipynb by downloading or cloning the GitHub project. For example, in the git command line type: 

`git clone https://github.com/iAMB-RWTH-Aachen/Exp2Ipynb.git`

Once installed locally, navigate into the `Exp2Ipynb` and start the Jupyter environment (in Linux command line `jupyter lab` or `jupyter notebook`, in Windows-Anaconda click the Jupyter Lab icon). In the Jupyter environment, navigate into the `Analysis` folder and open the Notebooks files.

## Overview

The analysis is distributed in different notebooks in the sub-folder `Analysis`:

| *.ipynb | Description |
| ------ | ------ |
| 0-Worflow | Workflow set-up. Parameters are defined for the machine learning type, threshold for positional analysis, and output file names and types. |
| 1-Statistical-Analysis | Statistical analysis of sequence and reporter. The notebook reports the unique promoter number, sequence and position sampling diversity, and reporter cross-correlation. |
| 2-Regressor-Training | Machine learning training. The data is separated into training and test sets and and trained to the defined machine learning tool. |
| 3-Regressor-Performance | Performance evaluation of machine learning. The machine learning regressor is loaded and evaluated based on cross validation and feature importance. |
| 4-Promoter-Prediction | Activity prediction of defined sequences. The activity of single sequences can be assessed as well as to identify sequences with defined activity. |

## Data preparation and statistical analysis

The data input is a comma-separated-value file (`csv`) located in the `Analysis` directory with at least three columns: (i) an identifier column, (ii) a sequence column, and (iii) an expression value column, with header names. Optional columns are expression values of the sequences in other organisms or with a different reporter system and the standard deviation with the replicate numbers. The sequence column only accepts DNA abbreviations (A, C, G, T) with an identical length for each sequence. The output file names and figure file types can be defined. The column names for data import are defined in a separate configuration file `config.txt` and the notebook `0-Workflow` guides through its construction

The statistical analysis is conducted in the notebook `1-Statistical-Analysis`. In addition to a single expression value, the standard deviation and replicate number can be provided. Optionally, outliers in the original data set can be removed from further analysis. Machine learning performance is improved if replicates are available and the workflow enables the re-generation of replicates based on mean and standard deviation. The replicates are calculated from Python numpy random normal function and are valid for normal-distributed data while adding a reasonable prediction bias.

## Predictor training and performance

The workflow supports classification and regression. In regression, the expression values for sequences are quantitatively predicted but require high data quality (sufficient sample size and position entropy). Classification provides qualitative predictions (*e.g.*, low-medium-high), with more reliable predictions even for small data sets. The implemented machine learning routines are random forest (RF), gradient boosting trees (GB), and support vector machines (SVM). The input features are nucleotides on each position in one-hot encoding plus the overall sequence GC-content. The predicted target variable is the expression strength. The feature size depends on the sequence length, typically ranging between 40-100 nucleotides or 160-400 features. If a classification is chosen, the output is binned according to a parameter provided by the user (`Response_Value` in `config.txt`). If `Response_Value=1`, the original activity values are used; if the `Response_Value=0`, the data is centered with zero mean and unit variance and for larger values, bins are generated as equal-sized buckets (python pandas `qcut`), and the bin label is used as the target prediction. The data is split into training and test set, with a default ratio of 9:1 and a grid search on the training set identifies the optimal hyper-parameter.

## Synthetic sequence identification by optimization

The ability to predict new sequences with defined activities is required for effective strain engineering. The predictors resulting from the machine-learning training are used to search sequence and expression with a genetic algorithm based on the Python framework [DEAP](https://pypi.org/project/deap/). The genetic optimization algorithm is used as it can be easily applied to a wide variety of machine learning models. Sequence identification for a regressor is conducted with a search for the sequence with the closest expression. The classification requires that the predicted expression lies within the target expression class while the sequence distance to the reference sequences is minimized. The sequences already present in the library are excluded from the search. 
