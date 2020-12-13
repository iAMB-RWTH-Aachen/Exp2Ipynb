ExpressionExpertIpynb

The workflow facilitates analysis of promoter libraries using Jupyter Notebooks. Clone the repository or directly use it in the cloud with binder by clicking the following binder symbol.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/uliebal/exp2ipynb/master)

The analysis is distributed in different notebooks:

| *.ipynb | Description |
| ------ | ------ |
| 0-Worflow | Workflow set-up. Parameters are defined for the machine learning type, threshold for positional analysis, and output file names and types. |
| 1-Statistical-Analysis | Statistical analysis of sequence and reporter. The notebook reports the unique promoter number, sequence and position sampling diversity, and reporter cross-correlation. |
| 2-Regressor-Training | Machine learning training. The data is separated into training and test sets and and trained to the defined machine learning tool. |
| 3-Regressor-Performance | Performance evaluation of machine learning. The machine learning regressor is loaded and evaluated based on cross validation and feature importance. |
| 4-Exploration-Space | Sampling from the predictable sequence space. The machine learning regressor is used to predict sequences that are covered by the library sampling space. |
| 5-Promoter-Prediction | Activity prediction of defined sequences. The activity of single sequences can be assessed as well as to identify sequences with defined activity. |


Last Version: 2020/12/13<br>
Author: Ulf Liebal<br>
Contact: ulf.liebal@rwth-aachen.de<br>
License: see LICENSE file<br>
