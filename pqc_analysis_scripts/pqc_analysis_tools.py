import glob
import json
import matplotlib.pyplot as plt
from analysis_pqc import *
import numpy as np
import os
import sys


### This script contains a number of functions which are helpful for the PQC analysis.####


def find_most_recent_file(path, test):

    # This function takes the pathfile and the name of test as inputs and returns the
    # most recent file which corresponds to the selected test

  path_folder = path

  filedir = []

  files = glob.glob(path_folder + os.sep + "*.json")    # get all json files
  files.sort(key=os.path.getmtime)             #sort the folders accordind to datetime

  for f in files:
      if test in [v.lower() for v in f.split('_')]:
          filedir.append(f)

  fil = filedir[-1]               #get the last, i.e the most recent file

  return fil




def assign_label(path, test):

    # This function assigns the ID to the variable lbl which is used in the plots

    file = find_most_recent_file(path, test)
    # file = path   #un-comment this out in case you want to test a specific file. You have to assign it as a path through the terminal

    file.split('\\')[-1]
    lbl = '_'.join(file.split('_')[6:9])
    return lbl


def read_json_file(path, test, parameter):

    # This function reads the json file and returns an array of the parameter you need

    file = find_most_recent_file(path, test)
    #file = path  #un-comment this out in case you want to test a specific file. You have to assign it as a path through the terminal

    with open(file) as f:
       a = json.load(f)
       x = np.array([i for i in a['series'][parameter]])

    return x


def units(data, unit):

    # This function converts the unit scale to the correct one

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

    # This function scales the data and return the latter and the units in the correct form

    denominator, unit = units(parameter, unit)
    x = np.array([j / denominator for j in parameter])
    return x, unit



def plot_curve(ax, x, y, title, xlabel, ylabel, legend, annotate, x_loc, y_loc):

    # This function plots the x,y data and sets the desired labes into the axis etc
    # ax is an axes object and is necessary to create two or more subplots. In that case
    # we want to overlay the curve and the fit

    ax.plot(x, y, '-o', ms=3, label=legend)
    plt.annotate(annotate, (x_loc, y_loc), xycoords='figure fraction', color='black', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='upper left')
    plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    plt.tight_layout()
    plt.show()


def fit_curve(ax, x, y1, y2):

    # This function returns the object which corresponds to the fit function.
    # You can plot 2 fit functions with respect to the same x, but if you want just one,
    # then you can type y2=0.
    # This should be modified according to what the user wants.

    ax = plt.gca()
    ax.plot(x, y1, '--r')
    if y2 is not 0:
      ax.plot(x, y2, '--r')

    return ax
