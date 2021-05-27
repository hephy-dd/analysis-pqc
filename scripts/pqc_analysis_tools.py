"""This script contains a number of functions which are helpful for the PQC
analysis.
"""

import glob
import json
import os
import dateutil.parser as timestamp_parser

import matplotlib.pyplot as plt
import numpy as np

__all__ = [
    'find_most_recent_file',
    'find_all_files_from_path',
    'get_timestamp',
    'assign_label',
    'read_json_file',
    'units',
    'normalise_parameter',
    'plot_curve',
    'fit_curve'
]


def find_all_files_from_path(path, test=None, *, whitelist=None, blacklist=None, pattern='*.json'):
    """
    returns a list of measurements for a given test
    optionally filters also with whitelist and blacklist (all whitelists must match and no blacklist match)
    eg for forward right poly vdp:
      whitlelist=["PQCFlutesRight","polyslicon"] and
      blacklist=["reverse"]
    """
    filedir = []

    filenames = glob.glob(os.path.join(path, pattern))
    filenames.sort()

    for filename in filenames:
        # the replace is necessary for van_der_pauw/van-der-pauw
        segments = [v.lower().replace("-", "_") for v in filename.split('_')]
        if (test is None or test in segments) and \
           (blacklist is None or not any(e.lower() in segments for e in blacklist)) and \
           (whitelist is None or all(e.lower() in segments for e in whitelist)):
            filedir.append(filename)

    return np.sort(filedir)


def find_most_recent_file(path, test=None, *, whitelist=None, blacklist=None):
    """This function takes the pathfile and the name of test as inputs and
    returns the most recent file which corresponds to the selected test. This is obtained via the name,
    not with the timestamp, so it gives a warning if the names don'z match for all candidats
    """
    all_files = find_all_files_from_path(path, test, whitelist=whitelist, blacklist=blacklist)
    if len(all_files) == 0:
        return None

    all_files.sort()  # the date should be the only difference, so we can sort

    # when we mix up measurements this gives a warning here (if the name is not equal (except for the timestamp))
    basename_first = '_'.join(all_files[-1].split('_')[0:-1])
    for f in all_files:
        basename = '_'.join(f.split('_')[0:-1])
        if basename != basename_first:
            print("Warning: heterogenous naming: {}".format(os.path.basename(basename)))
            print("   Info: first one was:       {}".format(os.path.basename(basename_first)))
            break  # we only want to show this once

    return all_files[-1]


def assign_label(path, test, vdp=False):
    """This function assigns the ID to the variable lbl which is used in the
    plots. if vdp is set to true, only the vdp info is extracted
    """
    file = path
    # un-comment this in case you want to test a specific file. You have to
    # assign it as a path through the terminal
    # file = path
    # print(path)
    lbl_list = [1, 2, 6, 8, 9]
    if vdp:
        lbl_list = [10, 11, 12, 13]
    basename = os.path.basename(file)
    try:
        lbl = '_'.join([basename.split('_')[i] for i in lbl_list])
    except IndexError:
        return path
    return lbl


def read_json_file(filename):
    """Return a PQC JSON formatted file as dictionary containing numpy arrays.

    >>> series = read_json_file('sample.json').get('series')
    >>> series.get('voltage')
    array([0.0, 0.1, 0.2, 0.3])
    """
    data = {"series": {}}
    try:
        with open(filename) as f:
            data = json.load(f)
        # convert to numpy arrays
        series = data.get('series', {})
        for k, v in series.items():
            series[k] = np.array(v)
    except Exception:
        raise RuntimeError(f"Failed to parse JSON file: {filename}")
    return data


def get_timestamp(filename):
    ret = None
    with open(filename) as json_file:
        data = json.load(json_file)
        ret = timestamp_parser.parse(data['meta']['start_timestamp'])
    return ret


def units(data, unit):
    """This function converts the unit scale to the correct one."""
    x = max(abs(data))
    unit_scale = ['P', 'T', 'G', 'M', 'k', '', 'm', '$\mu$', 'n', 'p', 'f', 'a', 'z', 'y']
    i = -1
    lower_limit = 1e-24

    max_scale = 1e+18
    while max_scale > lower_limit:
        previous_max = max_scale
        max_scale = max_scale / 1000
        i += 1
        if x >= max_scale and x < previous_max:
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


def plot_curve(ax, x, y, title, xlabel, ylabel, legend=None, annotate=None, x_loc=0, y_loc=0):
    """This function plots the x,y data and sets the desired labes into the axis
    etc. ax is an axes object and is necessary to create two or more subplots.
    In that case we want to overlay the curve and the fit.
    """

    # plt.figure()
    ax.plot(x, y, '-o', ms=3, label=legend)
    if annotate:
        plt.annotate(annotate, (x_loc, y_loc), xycoords='figure fraction', color='black', bbox=dict(facecolor='deepskyblue', alpha=0.75), horizontalalignment='right', verticalalignment='top')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if legend:
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
