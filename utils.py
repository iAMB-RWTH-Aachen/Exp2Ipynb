import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import os
import json

def save_model_to_json(filepath, filename, clf, scalers):
    """
    Write a json file to be used in Maingo / Melon
    return None
    """
    prediction_parameters = dict()
    prediction_parameters["rho"] = np.asscalar(clf.intercept_)
    prediction_parameters["support_vectors"] = clf.support_vectors_.tolist()
    prediction_parameters["dual_coefficients"] = clf.dual_coef_.ravel().tolist()
    prediction_parameters["kernel_parameters"] = [clf._gamma]
    prediction_parameters["kernel_function"] = clf.kernel

    prediction_parameters["scaling"] = dict()
    for var, scaler in scalers.items():
        if scaler is not None:
            if isinstance(scaler, MinMaxScaler):
                prediction_parameters["scaling"][var] = {
                    "scaler": "MinMax",
                    "min": scaler.data_min_.tolist(),
                    "max": scaler.data_max_.tolist()
                }    
            elif isinstance(scaler, StandardScaler):
                prediction_parameters["scaling"][var] = {
                    "scaler": "Standard",
                    "mean": scaler.mean_.tolist(),
                    "stddev": np.sqrt(scaler.var_).tolist()
                }
            else:
                raise TypeError("Can only write parameters for scalers of type 'MinMaxScaler' and 'StandardScaler' to JSON file.")
        else:
            prediction_parameters["scaling"][var] = {
                    "scaler": "None"
                }   

    if not os.path.exists(filepath):
        os.makedirs(filepath)
    
    with open(os.path.join(filepath,filename), 'w') as outfile:
        json.dump(prediction_parameters, outfile)


def plot_results_2D(clf, X, y, y_pred):
    lw = 2

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(15, 10), sharey=True)
    axes.plot(X, y_pred, color='m', lw=lw,
                    label='SVR prediction')
    axes.plot(X, np.sin(X).ravel(), color='b', lw=lw,
                    label='Original function')
    axes.scatter(X[clf.support_], y[clf.support_], facecolor="none",
                        edgecolor='m', s=50,
                        label='Support vectors')
    axes.scatter(X[np.setdiff1d(np.arange(len(X)), clf.support_)],
                        y[np.setdiff1d(np.arange(len(X)), clf.support_)],
                        facecolor="none", edgecolor="k", s=50,
                        label='other training data')
    axes.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
                    ncol=1, fancybox=True, shadow=True)

    fig.text(0.5, 0.04, 'data', ha='center', va='center')
    fig.text(0.06, 0.5, 'target', ha='center', va='center', rotation='vertical')
    fig.suptitle("Support Vector Regression", fontsize=14)
    plt.show()
