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
from scipy.stats import linregress
from analysis_pqc import *
from pqc_analysis_tools import *

print_results = 1


def analyse_iv_data(path, plot=True):
    test = 'iv'

    v = abs(read_json_file(path, test, 'voltage'))
    i_tot =abs( read_json_file(path, test, 'current_elm'))
    i = abs(read_json_file(path, test, 'current_hvsrc'))
    temp = read_json_file(path, test, 'temperature_chuck')
    humidity = read_json_file(path, test, 'humidity_box')

    x_loc = 0.3
    y_loc = 0.65

    v_norm, v_unit = normalise_parameter(v, 'V')
    i_norm, i_unit = normalise_parameter(i_tot, 'A')

    v_max, i_max, i_800, i_600, status = analyse_iv(v, i_tot)

    lbl = assign_label(path, test)

    annotate = 'I$_{max}$' + ': {:.3g} A @ {} V'.format(i_max, v_max) + '\n\nT$_{avg}$' + ': {:0.2f} $^\circ C$'.format(np.mean(temp)) \
          + '\n\n H$_{avg}$:' + '{:0.2f}'.format(np.mean(humidity)) + r'$\%$'

    
    if plot:
      fig,ax = plt.subplots(1,1)
      plot_curve(ax, v, i_tot, 'IV Curve', 'Reverse Bias Voltage [V]', 'Current [A]', lbl, annotate, x_loc, y_loc)
      plt.yscale('log')
    if print_results:
        print('%s:  IV:\ti_600: %.3f uA\ti_800: %.3f uA' % (lbl, i_600, i_800))

    return i_600, i_800



def analyse_cv_data(path, plot=True):
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
    c_norm, c_unit = normalise_parameter(c, 'F')

    inv_c2 = 1/c**2

   
    v_dep1, v_dep2, rho, conc, a_rise, b_rise, v_rise, a_const, b_const, v_const, spl_dev, status = analyse_cv(v, c)


    lbl = assign_label(path, test)
    annotate = 'V$_{{fd}}}}$: {} V\n\nT$_{{avg}}$: {} \u00B0C\nH$_{{avg}}$: {}'.format(v_dep2, round(np.mean(temp),2), round(np.mean(humidity),2)) + r'$\%$'


    if plot:
      fig2, ax2 = plt.subplots(1,1)
      fit_curve(ax2, v_rise, a_rise * v_rise+ b_rise, 0)
      fit_curve(ax2, v_const, a_const*v_const + b_const, 0)
      plot_curve(ax2, v, inv_c2, 'Full Depletion Voltage Estimation', 'Voltage[{}]'.format(v_unit), '1/C$^{2}$ [F$^{-2}$]', lbl, '', 0, 0 )
   
    if print_results:
    	
        print('%s: \tCV: v_fd: %.2e V\trho: %.2e Ohm\tconc: %.2e cm^-3' % (lbl, v_dep2, rho*1e-3, conc*1e-6))
   
    return v_dep2 



def analyse_mos_data(path, plot=True):
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

    v_fb1, v_fb2, c_acc, c_inv, t_ox, n_ox, Q_ox, a_acc, b_acc, v_acc, a_dep, b_dep, v_dep, a_inv, b_inv, v_inv,  spl_dev, status = analyse_mos(v, c)
    lbl = assign_label(path, test)


    fit_acc = [a_acc*x + b_acc for x in v]
    fit_dep = [a_dep*i+b_dep for i in v]
    annotate = 'V$_{{fb}}$: {} V (via intersection)\nV$_{{fb}}$: {} V (via inflection)\n\nCacc : {:.3g} F \nt$_{{ox}}$: {:.3g} m\nn$_{{ox}}$: {:.3g} cm$^{{-2}}$'.format(round(v_fb2,2), round(v_fb1,2), np.mean(c_acc), t_ox, n_ox)

    x_loc = 0.15
    y_loc = 0.145


    if plot:
      fig, ax = plt.subplots(1,1)
      plt.ylim(0, 1.5*np.max(c))
      fit_curve(ax, v, fit_acc, fit_dep)
      plt.axvline(x=v_fb2, color='black', linestyle='dashed')
      plot_curve(ax, v, c, 'CV Curve', 'Voltage [V]', 'Capacitance [F]', lbl, annotate, x_loc, y_loc)

    if print_results:

        print('%s: MOS: v_fb2: %.2e V \t t_ox: %.2e um \t n_ox: %.2e cm^-2 \t Mean Cacc: %.2e F' % (lbl, v_fb2, t_ox, n_ox, np.mean(c_acc)))      
  
  
    return v_fb1, v_fb2, t_ox, n_ox, Q_ox, np.mean(c_acc)

     

def analyse_gcd_data(path, plot=True):
  try:
    test = 'gcd'

    v = read_json_file(path, test, 'voltage')
    i_em = read_json_file(path, test, 'current_elm')
    i_src = read_json_file(path, test, 'current_vsrc')
    i_hvsrc = read_json_file(path, test, 'current_hvsrc')


    i_em_norm, i_em_unit = normalise_parameter(i_em, 'A')

    lbl = assign_label(path, test)

    i_surf, i_bulk, i_acc, i_dep, i_inv, v_acc, v_dep, v_inv, v_trans, a_acc, b_acc, a_dep, b_dep, a_inv, b_inv, v_fb2, v_fb3, spl_dev, status, idep_mean, iinv_mean, s0 = analyse_gcd(v,i_em_norm)

   
    if plot:
      v_draw = np.array([i for i in v if -7.5<i<0])
      v_draw2 = np.array([i for i in v if -2.5<i<v_inv[5]])
      fit_trans = np.array([a_acc*x + b_acc for x in v_trans])
      fit_dep = np.array([a_dep*x + b_dep for x in v_draw])
      fit_inv = np.array([a_inv*x + b_inv for x in v_draw2])
      annotate = 'I$_{{surf}}$ : {} pA \nI$_{{bulk}}$: {} pA'.format(round(i_surf,2), round(i_bulk,2))
      fig, ax = plt.subplots(1,1)
      plot_curve(ax, v, i_em_norm, 'I-V Curve GCD', 'Voltage [V]', 'Current [pA]', lbl, annotate, 0.15, 0.5)
      #plt.vlines(v_fb2, np.min(i_em_norm), np.max(i_em_norm), color='k', linestyle='dashed')
      fit_curve(ax, v_trans, fit_trans, 0)
      fit_curve(ax, v_draw, fit_dep, 0)
      fit_curve(ax, v_draw2, fit_inv, 0)
      plt.vlines(v_fb2, np.min(i_em_norm), np.max(i_em_norm), color='k', linestyle='dashed')
      plt.vlines(v_fb3, np.min(i_em_norm), np.max(i_em_norm), color='k', linestyle='dashed')

    if print_results:

        print('%s: GCD: i_surf: %.2e A\t i_bulk: %.2e A\t vfb_acc: %.2e \t vfb_inv: %.2e' % (lbl, i_surf, i_bulk, np.mean(v_acc), np.mean(v_inv)))

         
    return i_surf, i_bulk, v_fb2, v_fb3, s0
  
  except:
      return 0
   




def analyse_fet_data(path, plot=True):
    test = 'fet'

    v = read_json_file(path, test, 'voltage')
    i_em = read_json_file(path, test, 'current_elm')
    i_src = read_json_file(path, test, 'current_vsrc')
    i_hvsrc = read_json_file(path, test, 'current_hvsrc')


    i_em_norm, i_em_unit = normalise_parameter(i_em, 'A')

    v_th, a, b, spl_dev, status = analyse_fet(v, i_em)#_em_norm)

    fit  = [a*i +b for i in v]

    lbl = assign_label(path, test)


    if plot:
      fig,ax1 = plt.subplots()
      lns1 = ax1.plot(v,i_em, ls='', marker='s', ms=3, label='transfer characteristics')
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
      plt.annotate('V_th: {} V'.format(round(v_th,2)), (0.4, 0.15), xycoords='figure fraction', color='black', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5')) 
      plt.legend(lns, labs, loc='upper left')
      plt.show()

    if print_results:
       print('%s: \tnFet: v_th: %.2e V' % (lbl, v_th))
 
    return v_th



def analyse_van_der_pauw_data(path, plot=True):
   test = 'van-der-pauw'
    
   if find_json_parameter(path, test, 'voltage_vsrc'):

        v = read_json_file(path, test, 'voltage_vsrc')
   elif find_json_parameter(path, test, 'voltage'):
        v = read_json_file(path, test, 'voltage')

    
   i = read_json_file(path, test, 'current')


   lbl = assign_label(path, test)
   i_norm, i_unit = normalise_parameter(i, 'A')
   v_norm, v_unit = normalise_parameter(v, 'V')
   
   if len(i)!=0:
     r_sheet, a, b, x_fit, spl_dev, status = analyse_van_der_pauw(i,v)
     fit = [a*x +b for x in x_fit]

    
     if plot:
       fig, ax = plt.subplots(1,1)
       #annotate = 'Rsheet: {} Ohm/sq'.format(round(r_sheet,2)) 
       fit_curve(ax, x_fit, fit, 0)
       plot_curve(ax, i, v, 'IV Curve', 'Current [{}]'.format(i_unit), 'Voltage [{}]'.format(v_unit), lbl, '', 0, 0)  #annotate, 0.25, 0.7)
       plt.show()
   
     if print_results:
        print('%s: van der Pauw: r_sheet: %.2e Ohm/sq' % (lbl, r_sheet)) 

        return r_sheet, lbl
  
     else:
        print('No parameter found')



def analyse_linewidth_data(path, r_sheet, plot=True):
    test = 'linewidth'

    v = read_json_file(path, test, 'voltage_vsrc')
    i = read_json_file(path, test, 'current')

    i_norm, i_unit = normalise_parameter(i, 'A')
    v_norm, v_unit = normalise_parameter(v, 'V')

    lbl = assign_label(path, test)

    print(lbl)

    t_line, a, b, x_fit, spl_dev, status = analyse_linewidth(i, v, r_sheet=r_sheet, cut_param=0.01, debug=0)

       
    fit = [a*x +b for x in x_fit]


    if plot:
      fig, ax = plt.subplots(1, 1)
      fit_curve(ax, x_fit, fit, 0)
      plot_curve(ax, i_norm, v_norm, 'IV Curve', 'Current [{}]'.format(i_unit), 'Voltage [{}]'.format(v_unit), lbl, '', 0, 0)

    if print_results:
        print('%s: \tLinewidth: %.2e um' % (lbl, t_line))
 
    return t_line
  


def analyse_cross_data(path):

    test = 'bulk_cross'

    v = read_json_file(path, test, 'voltage_vsrc')
    i = read_json_file(path, test, 'current')

    i_norm, i_unit = normalise_parameter(i, 'A')
    v_norm, v_unit = normalise_parameter(v, 'V')

    lbl = assign_label(path, test)

    r_sheet, a, b, x_fit, spl_dev, status = analyse_cross(i,v)

    if print_results:
        print('%s: Bulk cross: r_sheet: %.2e Ohm/sq' % (lbl, r_sheet))


    return r_sheet, lbl



def analyse_cbkr_data(path, r_sheet, plot=True):
    test = 'cbkr'

    v = read_json_file(path, test, 'voltage_vsrc')
    i = read_json_file(path, test, 'current')

    i_norm, i_unit = normalise_parameter(i, 'A')

    lbl = assign_label(path, test)

   
    r_contact, a, b, x_fit, spl_dev, status = analyse_cbkr(i, v, r_sheet, cut_param=0.01, debug=0)
    fit = [a*x +b for x in x_fit]

    if plot:
      fig, ax = plt.subplots(1, 1)
      fit_curve(ax, x_fit, fit, 0)
      plot_curve(ax, i, v, 'IV Curve', 'Current [{}]'.format(i_unit), 'Voltage [V]', lbl, '', 0, 0)
   
    if print_results:
 
       print('%s: \tcbkr: r_contact: %.2e Ohm' % (lbl, r_contact))

 

    return r_contact

   
def analyse_contact_data(path):
    test= 'contact'

    v = read_json_file(path, test, 'voltage_vsrc')
    i = read_json_file(path, test, 'current')

   # i_norm, i_unit = normalise_parameter(i, 'A')

    lbl = assign_label(path, test)
    r_contact, a, b, x_fit, spl_dev, status = analyse_contact(i, v, cut_param=0.01, debug=0)

    fit = [a*x+b for x in x_fit]
    
    if print_results:
       print('%s: \tcontact: r_contact: %.2e Ohm' % (lbl, r_contact))
 
    return r_contact



def analyse_meander_data(path):
    test = 'meander'

    if find_json_parameter(path, test, 'voltage'):
        v = read_json_file(path, test, 'voltage')
    else:   
        v = read_json_file(path, test, 'voltage_vsrc')
    
    if find_json_parameter(path, test, 'current_elm'):
         i = read_json_file(path, test, 'current_elm')
    else:
        i = read_json_file(path, test, 'current')

    i_norm, i_unit = normalise_parameter(i, 'A')

    lbl = assign_label(path, test)

    print(lbl)
    rho_sq, status = analyse_meander(i, v, path, debug=0)
    
    if print_results:
       print('%s: \tMeander: rho_sq: %.2e' % (lbl, rho_sq))
   

    return rho_sq



def analyse_breakdown_data(path, plot=True):
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

    if plot:
      fig, ax = plt.subplots(1,1)
      annotate = 'V$_{{bd}}$: {} V \n\nT$_{{avg}}$ : {} \u00B0C \nH$_{{avg}}$: {} $\%$ '.format(v_bd, round(np.mean(temp),2), round(np.mean(humidity),2))
      plot_curve(ax, v, i_elm_norm, 'IV Curve', 'Voltage [V]', 'Current [{}]'.format(i_elm_unit), lbl, annotate, x_loc, y_loc)


    if print_results:
       print('%s: \tBreakdown: v_bd: %.2e V' % (lbl, v_bd))
   
    return v_bd


def analyse_capacitor_data(path, plot=True):

    test = 'capacitor'

    v = read_json_file(path, test, 'voltage_hvsrc')
    c = read_json_file(path, test, 'capacitance')
    
    c_norm, c_unit = normalise_parameter(c, 'F')

    c_mean, c_median, d,  status = analyse_capacitor(v, c)

    lbl = assign_label(path, test)


    if plot:
      fig, ax = plt.subplots(1,1)
      annotate = 'C$_{{median}}$: {:.3g} F \n \n Oxide thickness: {:.3g} m'.format(c_median, d)
      plot_curve(ax, v, c_norm, 'CV Curve', 'Voltage [V]', 'Capacitance [{}]'.format(c_unit), lbl, annotate, 0.3, 0.6)

    if print_results:
        print('%s: Capacitor: Mean_capacitance: %.2e F \t Median: %.2e F \t Oxide thicknes: %.2e' % (lbl, c_mean, c_median, d))

    return c_mean, c_median





functions = {
        'iv': analyse_iv_data,
        'cv': analyse_cv_data,
        'fet': analyse_fet_data,
        'gcd': analyse_gcd_data,
        'mos': analyse_mos_data,
        'linewidth': analyse_linewidth_data,
        'van-der-pauw': analyse_van_der_pauw_data,
        'bulk': analyse_cross_data,
        'cbkr': analyse_cbkr_data,
        'contact': analyse_contact_data,
        'meander': analyse_meander_data,
        'breakdown': analyse_breakdown_data,
        'capacitor': analyse_capacitor_data
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
        'bulk_cross: analyse_cross_data,
        'breakdown': analyse_breakdown_data
    }
    '''
    if test == 'all':
        for f in functions.values():
            print(f)
            f(path, True)
    elif test in functions:
       #v, i, lbl = 
       try:
           functions.get(test)(path)
       except (ValueError, TypeError, IndexError):
           print("Error in datafile")
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

    j=-1
    colors = ['c', 'b','g','r','m', 'y', 'black', 'limegreen', 'slategrey', 'darkblue', 'gold', 'brown', 'aqua', 'indigo', 'thistle', 'pink']
    for test in tests:
        print(test)
        filedir = find_all_files_from_path(args.path, test)
        for f in filedir:
          j +=1  
            
          analyse_file(f, test)
         
          plt.show()

if __name__ =="__main__":
    main()
