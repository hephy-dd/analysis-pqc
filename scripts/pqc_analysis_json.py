#!/usr/bin/env python3

"""This script contains a number of functions which can be used to analyse the
PQC measurements.

By default it finds the most recent file of the path and the test (e.g iv) that
you have given as an input in the terminal and analyses the data.

If you select 'all' as an input test, then it finds the most recent file of each
individual possible test and analyses it.

You may want to test a particular file (which might not be the most recent one).
In that case, you can just modify the corresponding functions in the
pqc_analysis_tools.py script and give the file (and the path) as an input in the
terminal.
"""

import argparse
import math
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

from analysis_pqc import *
from pqc_analysis_tools import *

print_results = 1


'''
functions = {
      'iv': analyse_iv_data,
      'cv': analyse_cv_data,
      'fet': analyse_fet_data,
      'gcd': analyse_gcd_data,
      'mos': analyse_mos_data,
      'linewidth': analyse_linewidth_data,
      'van_der_pauw': analyse_van_der_pauw_data, 
      'breakdown': analyse_breakdown_data}
'''

def analyse_iv_data(path):
    test = 'iv'

    v = abs(read_json_file(path, test, 'voltage'))
    i_tot = read_json_file(path, test, 'current_elm')
    i = abs(read_json_file(path, test, 'current_hvsrc'))
    temp = read_json_file(path, test, 'temperature_chuck')
    humidity = read_json_file(path, test, 'humidity_box')

    x_loc = 0.3
    y_loc = 0.65

    v_norm, v_unit = normalise_parameter(v, 'V')
    i_norm, i_unit = normalise_parameter(i, 'A')

    v_max, i_max, i_800, i_600, status = analyse_iv(v_norm, i_norm)

    lbl = assign_label(path, test)

    annotate = 'I$_{max}$' + ': {} {} @ {}{}'.format(round(i_max, 2), i_unit, max(v_norm),v_unit) + '\n\nT$_{avg}$' + ': {:0.2f} $^\circ C$'.format(np.mean(temp)) \
          + '\n\n H$_{avg}$:' + '{:0.2f}'.format(np.mean(humidity)) + r'$\%$'

    fig,ax = plt.subplots(1,1)
    plot_curve(ax, v_norm, i_norm, 'IV Curve', 'Reverse Bias Voltage [{}]'.format(v_unit), 'Current [{}]'.format(i_unit), lbl, annotate, x_loc, y_loc)

    if print_results:
        print('%s:  IV:\ti_600: %.3f uA\ti_800: %.3f uA' % (lbl, i_600, i_800))

    return i_600, i_800


def analyse_cv_data(path):
    test = 'cv'

    v = abs(read_json_file(path, test, 'voltage_hvsrc'))
    i = read_json_file(path, test, 'current_hvsrc')
    c = read_json_file(path, test, 'capacitance')
    c2 = read_json_file(path, test, 'capacitance2')
    r = read_json_file(path, test, 'resistance')
    temp = read_json_file(path, test, 'temperature_chuck')
    humidity = read_json_file(path, test, 'humidity_box')

    x_loc = 0.3
    y_loc = 0.65

    v_norm, v_unit = normalise_parameter(v, 'V')
    c_norm, c_unit = normalise_parameter(c, 'F')  # then we get this is pF, we don't want that
    
    inv_c2 = 1/c**2
    
    lbl = assign_label(path, test)
    if "Flute_1" in path:
    	area = 1.56e-6  # m^2, quarter
    elif "Flute_3" in path:
    	area = 6.25e-6  # m^2, half (but without rounded edges)
    else:
    	area = 1
    	print("WARNING: clould not determine flute number - area dependent values will be wrong!")

    v_dep1, v_dep2, rho, conc, a_rise, b_rise, v_rise, a_const, b_const, v_const, spl_dev, status = analyse_cv(v, c, area=area, cut_param= 0.008)

    
    annotate = 'V$_{{fd}}}}$: {} V\n\nT$_{{avg}}$: {} \u00B0C\nH$_{{avg}}$: {}'.format(v_dep2, round(np.mean(temp),2), round(np.mean(humidity),2)) + r'$\%$'

    #fig1, ax1 = plt.subplots(1, 1)
    #plot_curve(ax1, v_norm, c_norm, 'CV curve', 'Voltage[{}]'.format(v_unit), 'Capacitance [{}]'.format(c_unit), lbl, annotate, x_loc, y_loc)

    fig2, ax2 = plt.subplots(1,1)
    #ax2b = ax2.twinx()
    fit_curve(ax2, v_rise, a_rise * v_rise+ b_rise, color='ro')
    fit_curve(ax2, v_const, a_const * v_const+ b_const, color='kx')
    #fit_curve(ax2b, v_norm, spl_dev, color='mx')
    plot_curve(ax2, v, 1./c**2, 'Full Depletion Voltage Estimation', 'Voltage[{}]'.format(v_unit), '1/C$^{2}$ [F$^{-2}$]', lbl, '', 0, 0 )

    if print_results:
    	#print(f"{lbl}: CV: v_fd: {}")
        print('%s: \tCV: v_fd: %.2e V\trho: %.2e kOhm\tconc: %.2e cm^-3' % (lbl, v_dep2, rho*1e-3, conc*1e-6))
   
    return v_dep2 


def analyse_mos_data(path):
    test = 'mos'

    v = read_json_file(path, test, 'voltage_hvsrc')
    i = read_json_file(path, test, 'current_hvsrc')
    c = read_json_file(path, test, 'capacitance')
    c2 = read_json_file(path, test, 'capacitance2')
    r = read_json_file(path, test, 'resistance')
    temp = read_json_file(path, test, 'temperature_chuck')
    humidity = read_json_file(path, test, 'humidity_box')

    v_norm, v_unit = normalise_parameter(v, 'V')
    c_norm, c_unit = normalise_parameter(c, 'F')

    v_fb1, v_fb2, c_acc, c_inv, t_ox, n_ox, a_acc, b_acc, v_acc, a_dep, b_dep, v_dep, a_inv, b_inv, v_inv,  spl_dev, status = analyse_mos(v_norm, c_norm)
    lbl = assign_label(path, test)


    fit_acc = [a_acc*x + b_acc for x in v_norm]
    fit_dep = [a_dep*i+b_dep for i in v_norm]
    annotate = 'V$_{{fb}}$: {} V (via intersection)\nV$_{{fb}}$: {} V (via inflection)\n\nt$_{{ox}}$: {} m\nn$_{{ox}}$: {} cm$^{{-2}}$'.format(round(v_fb2,2), round(v_fb1,2), round(t_ox, 2), round(n_ox, 2))

    x_loc = 0.15
    y_loc = 0.145

    fig, ax = plt.subplots(1,1)
    plt.ylim(0, 100)
    fit_curve(ax, v_norm, fit_acc, fit_dep)
    plt.axvline(x=v_fb2, color='black', linestyle='dashed')
    plot_curve(ax, v_norm,c_norm, 'CV Curve', 'Voltage [V]', 'Capacitance [{}]'.format(c_unit), lbl, annotate, x_loc, y_loc)

    if print_results:
        print('%s: \tMOS: v_fb2: %.2e V\tt_ox: %.2e um\tn_ox: %.2e cm^-3' % (lbl, v_fb2, t_ox, n_ox))
 

    return v_fb1, v_fb2, t_ox, n_ox

     
def analyse_gcd_data(path):
    test = 'gcd'

    v = read_json_file(path, test, 'voltage')
    i_em = read_json_file(path, test, 'current_elm')
    i_src = read_json_file(path, test, 'current_vsrc')
    i_hvsrc = read_json_file(path, test, 'current_hvsrc')


    i_em_norm, i_em_unit = normalise_parameter(i_em, 'A')

    lbl = assign_label(path, test)

    i_surf, i_bulk, i_acc, i_dep, i_inv, v_acc, v_dep, v_inv, spl_dev, status = analyse_gcd(v,i_em_norm)

    fig, ax = plt.subplots(1,1)
    plot_curve(ax, v, i_em_norm, 'I-V Curve GCD', 'Voltage [V]', 'Current [{}]'.format(i_em_unit), lbl, '', 0, 0)

    if print_results:
        print('%s: \tGCD: i_surf: %.2e A\t i_bulk: %.2e A' % (lbl, i_surf, i_bulk))
  
    return i_surf, i_bulk



def analyse_fet_data(path):
    test = 'fet'

    v = read_json_file(path, test, 'voltage')
    i_em = read_json_file(path, test, 'current_elm')
    i_src = read_json_file(path, test, 'current_vsrc')
    i_hvsrc = read_json_file(path, test, 'current_hvsrc')


    i_em_norm, i_em_unit = normalise_parameter(i_em, 'A')

    v_th, a, b, spl_dev, status = analyse_fet(v, i_em_norm)

    fit  = [a*i +b for i in v]

    lbl = assign_label(path, test)

    fig,ax1 = plt.subplots()
    lns1 = ax1.plot(v,i_em_norm, ls='', marker='s', ms=3, label='transfer characteristics')
    ax1.set_xlabel('V$_{GS}$ [V]')
    ax1.set_ylabel('I$_{D}$ [A]')
    ax1.set_ylabel(r'I$_\mathrm{D}$ [A]')
    ax2 = ax1.twinx()
    lns2 = ax2.plot(v, spl_dev, ls=' ', marker='s', ms=3, color='tab:orange', label="transconductance")
    ax2.tick_params(axis='y', labelcolor='tab:orange')
    ax2.set_ylabel(r'g$_\mathrm{m}$ [S]', color='tab:orange')
    lns3 = ax1.plot(v, fit, '--r', label="tangent")
    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    plt.legend(lns, labs, loc='upper left')
    plt.show()

    if print_results:
       print('%s: \tnFet: v_th: %.2e V' % (lbl, v_th))
 
    return v_th


def analyse_van_der_pauw_data(path):
    test = 'van-der-pauw'

    v = read_json_file(path, test, 'voltage_vsrc')
    i = abs(read_json_file(path, test, 'current'))

    lbl = assign_label(path, test)
    r_sheet, a, b, x_fit, spl_dev, status = analyse_van_der_pauw(i, v)
    fit = [a*x +b for x in x_fit]

    fig, ax = plt.subplots(1,1)
    fit_curve(ax, x_fit, fit, 0)
    plot_curve(ax, i, v, 'IV Curve', 'Current', 'Voltage', lbl, '', 0, 0)
    
    if print_results:
       print('%s: \tvan der Pauw: r_sheet: %.2e Ohm/sq' % (lbl, r_sheet))
 

    return r_sheet


def analyse_linewidth_data(path):
    test = 'linewidth'

    v = read_json_file(path, test, 'voltage_vsrc')
    i = read_json_file(path, test, 'current')

    i_norm, i_unit = normalise_parameter(i, 'A')

    lbl = assign_label(path, test)
    t_line, a, b, x_fit, spl_dev, status = analyse_linewidth(i_norm, v, r_sheet=-1, cut_param=0.01, debug=0)

    fit = [a*x +b for x in x_fit]

    fig, ax = plt.subplots(1, 1)
    fit_curve(ax, x_fit, fit, 0)
    plot_curve(ax, i_norm, v, 'IV Curve', 'Current [{}]'.format(i_unit), 'Voltage [V]', lbl, '', 0, 0)

    if print_results:
        print('%s: \tLinewidth: %.2e um' % (lbl, t_line))
 
    return t_line
  

def analyse_cbkr_data(path, r_sheet=-1):
    test = 'cbkr'

    v = read_json_file(path, test, 'voltage_vsrc')
    i = read_json_file(path, test, 'current')

    i_norm, i_unit = normalise_parameter(i, 'A')

    lbl = assign_label(path, test)

    r_contact, a, b, x_fit, spl_dev, status = analyse_cbkr(i_norm, v, r_sheet, cut_param=0.01, debug=0)
    fit = [a*x +b for x in x_fit]

    fig, ax = plt.subplots(1, 1)
    fit_curve(ax, x_fit, fit, 0)
    plot_curve(ax, i_norm, v, 'IV Curve', 'Current [{}]'.format(i_unit), 'Voltage [V]', lbl, '', 0, 0)
   
    if print_results:
       print('%s: \tcbkr: r_contact: %.2e Ohm' % (lbl, r_contact))
 

    return r_contact

   
def analyse_contact_data(path):
    test= 'cbkr'

    v = read_json_file(path, test, 'voltage_vsrc')
    i = read_json_file(path, test, 'current')

    i_norm, i_unit = normalise_parameter(i, 'A')

    lbl = assign_label(path, test)
    r_contact, a, b, x_fit, spl_dev, status = analyse_contact(i, v, cut_param=0.01, debug=0)

    fit = [a*x+b for x in x_fit]
    
    if print_results:
       print('%s: \tcontact: r_contact: %.2e Ohm' % (lbl, r_contact))
 
    return r_contact



def analyse_meander_data(path):
    test = 'meander'

    v = read_json_file(path, test, 'voltage_vsrc')
    i = read_json_file(path, test, 'current')

    i_norm, i_unit = normalise_parameter(i, 'A')

    lbl = assign_label(path, test)

    rho_sq, status = analyse_meander(i_norm, v, debug=0)
    
    if print_results:
       print('%s: \tMeander: rho_sq: %.2e' % (lbl, rho_sq))
   

    return rho_sq



def analyse_breakdown_data(path):
    test = 'breakdown'

    v = read_json_file(path, test, 'voltage')
    i = read_json_file(path, test, 'current_hvsrc')
    i_elm = read_json_file(path, test, 'current_elm')
    temp = read_json_file(path, test, 'temperature_chuck')
    humidity = read_json_file(path, test, 'humidity_box')

    i_elm_norm, i_elm_unit = normalise_parameter(i_elm, 'A')

    lbl = assign_label(path, test)
    x_loc = 0.3
    y_loc = 0.5

    v_bd, status = analyse_breakdown(v, i_elm_norm, debug=0)

    fig, ax = plt.subplots(1,1)
    annotate = 'V$_{{bd}}$: {} V \n\nT$_{{avg}}$ : {} \u00B0C \nH$_{{avg}}$: {} $\%$ '.format(v_bd, round(np.mean(temp),2), round(np.mean(humidity),2))
    plot_curve(ax, v, i_elm_norm, 'IV Curve', 'Voltage [V]', 'Current [{}]'.format(i_elm_unit), lbl, annotate, x_loc, y_loc)


    if print_results:
       print('%s: \tBreakdown: v_bd: %.2e V' % (lbl, v_bd))
   
    return v_bd


functions = {
        'iv': analyse_iv_data,
        'cv': analyse_cv_data,
        'fet': analyse_fet_data,
        'gcd': analyse_gcd_data,
        'mos': analyse_mos_data,
        'linewidth': analyse_linewidth_data,
        'van_der_pauw': analyse_van_der_pauw_data,
        'breakdown': analyse_breakdown_data
    }


def analyse_file(path, test):
    '''
    functions = {
        'iv': analyse_iv_data,
        'cv': analyse_cv_data,
        'fet': analyse_fet_data,
        'gcd': analyse_gcd_data,
        'mos': analyse_mos_data,
        'linewidth': analyse_linewidth_data,
        'van_der_pauw': analyse_van_der_pauw_data,
        'breakdown': analyse_breakdown_data
    }
    '''
    if test == 'all':
        for f in functions.values():
            print(f)
            f(path)
    elif test in functions:
        functions.get(test)(path)
    else:
        raise ValueError(f"no such test: '{test}'")



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    parser.add_argument('test')
    return parser.parse_args()


def main():
    args = parse_args()
    
    if args.test == 'all':
        tests = functions.keys()
    else:
        tests = [args.test]

    for test in tests:
        print(test)
        filedir = find_all_files_from_path(args.path, test)
        for f in filedir:
         #  print(f)
           analyse_file(f, test)
        plt.show()

if __name__ =="__main__":
    main()
