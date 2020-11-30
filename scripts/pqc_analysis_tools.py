"""This script contains a number of functions which are helpful for the PQC
analysis.
"""

import glob
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

from analysis_pqc import *


__all__ = [
    'find_most_recent_file',
    'find_all_files_from_path',
    'assign_label',
    'read_json_file',
    'units',
    'normalise_parameter',
    'plot_curve',
    'fit_curve'
]

def find_all_files_from_path(path, test, whitelist=None, blacklist=None, single=False):
    """
    returns a list of measurements for a given test
    optionally filters also whitelist an blacklist (all whitelists must match and no blacklist match)
    eg for forward right poly vdp:
      whitlelist=["PQCFlutesRight","polyslicon"] and
      blacklist=["reverse"]
    if single is set to true, we expect only one return value (not a list) and it throws an error if there is more
    """
    
    path_folder = path
    filedir =[]

    files = glob.glob(os.path.join(path_folder, "*.json"))
    files.sort(key=os.path.getmtime)
    
    for f in files:
        #the replace is necessary for van_der_pauw/van-der-pauw
        segments = [v.lower().replace("-", "_") for v in f.split('_')]
        if (test in segments) and \
           (blacklist is None or not any(e.lower() in segments for e in blacklist)) and \
           (whitelist is None or all(e .lower() in segments for e in whitelist)):
            filedir.append(f)
    if single:
        if len(filedir) > 1:
            #pass
            print("Warning: more than one measurement available, taking the most recent one!")
        if len(filedir) == 0:
            return None
        return filedir[-1]  # we sorted it first for time

    return np.sort(filedir)



def find_most_recent_file(path, test):
    """This function takes the pathfile and the name of test as inputs and
    returns the most recent file which corresponds to the selected test.
    """
    path_folder = path
    filedir = []

    # get all json files
    files = glob.glob(os.path.join(path_folder, "*.json"))
    # sort the folders accordind to datetime
    files.sort(key=os.path.getmtime)

    for f in files:
        if test in [v.lower() for v in f.split('_')]:
            filedir.append(f)

    # get the last, i.e the most recent file
    return filedir[-1]


def assign_label(path, test, vdp=False):
    """This function assigns the ID to the variable lbl which is used in the
    plots. if vdp is set to true, only the vdp info is extracted
    """
    file = path
    # un-comment this in case you want to test a specific file. You have to
    # assign it as a path through the terminal
    # file = path
    #print(path)
    lbl_list =[1, 2, 6, 8, 9]
    if vdp:
        lbl_list = [10, 11, 12]
    file = file.split(os.sep)[-1]
    lbl = '_'.join([file.split('_')[i] for i in lbl_list])
    return lbl


def read_json_file(path, test, parameter):
    """This function reads the json file and returns an array of the parameter
    you need.
    """
    
    file = path  #un-comment this out in case you want to test a specific file. You have to assign it as a path through the terminal
    try:
        with open(file) as f:
            a = json.load(f)
        return np.array([i for i in a['series'][parameter]])
    except:
        return []


def units(data, unit):
    """This function converts the unit scale to the correct one."""
    x = max(abs(data))
    numer = 0
    unit_scale = ['P', 'T', 'G', 'M', 'k', '', 'm', '$\mu$', 'n', 'p', 'f', 'a', 'z', 'y']
    i = -1
    lower_limit = 1e-24

    max_scale = 1e+18
    while max_scale> lower_limit:
        previous_max = max_scale
        max_scale = max_scale/1000
        i +=1
        if x >= max_scale and x< previous_max:
            numerator = max_scale
            string = '{}{}'.format(unit_scale[i], unit)

    return numerator, string


def normalise_parameter(parameter, unit):
    """This function scales the data and return the latter and the units in the
    correct form.
    """
    denominator, unit = units(parameter, unit)
    x = np.array([j / denominator for j in parameter])
    return x, unit


def plot_curve(ax, x, y, title, xlabel, ylabel, legend, annotate, x_loc, y_loc):
    """This function plots the x,y data and sets the desired labes into the axis
    etc. ax is an axes object and is necessary to create two or more subplots.
    In that case we want to overlay the curve and the fit.
    """

    #plt.figure()
    ax.plot(x, y, '-o', ms=3, label=legend)
    plt.annotate(annotate, (x_loc, y_loc), xycoords='figure fraction', color='black', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='upper left')
    plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    plt.tight_layout()
   # plt.show()


def fit_curve(ax, x, y1, y2=None, color='r'):
    """This function returns the object which corresponds to the fit function.
    You can plot 2 fit functions with respect to the same x, but if you want
    just one, then you can omit y2.

    This should be modified according to what the user wants.
    """
    #ax = plt.gca()
    ax.plot(x, y1, '--'+color)
    if y2:
      ax.plot(x, y2, '--'+color)

    return ax
