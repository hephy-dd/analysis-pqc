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
import glob

import matplotlib.pyplot as plt
import numpy as np

from analysis_pqc import *
from pqc_analysis_tools import *

print_results = 1


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

    v_max, i_max, i_800, i_600, status = analyse_iv(v, i)

    lbl = assign_label(path, test)

    annotate = 'I$_{max}$' + ': {} {} @ {}{}'.format(round(i_max, 2), i_unit, max(v_norm),v_unit) + '\n\nT$_{avg}$' + ': {:0.2f} $^\circ C$'.format(np.mean(temp)) \
          + '\n\n H$_{avg}$:' + '{:0.2f}'.format(np.mean(humidity)) + r'$\%$'

    fig,ax = plt.subplots(1,1)
    plot_curve(ax, v_norm, i_norm, 'IV Curve', 'Reverse Bias Voltage [{}]'.format(v_unit), 'Current [{}]'.format(i_unit), lbl, annotate, x_loc, y_loc)

    if print_results:
        print('%s:  IV:\ti_600: %.3f uA\ti_800: %.3f uA' % (lbl, i_600*1e6, i_800*1e6))

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
        print('%s: \tCV: v_fd: %.2e V\trho: %.2e Ohm\tconc: %.2e cm^-3' % (lbl, v_dep2, rho, conc*1e-6))
   
    return v_dep2 


def analyse_mos_data(path, plotResults=True, printResults=print_results):
    test = 'mos'

    v = read_json_file(path, test, 'voltage_hvsrc')
    i = read_json_file(path, test, 'current_hvsrc')
    c = read_json_file(path, test, 'capacitance')
    c2 = read_json_file(path, test, 'capacitance2')
    r = read_json_file(path, test, 'resistance')
    temp = read_json_file(path, test, 'temperature_chuck')
    humidity = read_json_file(path, test, 'humidity_box')
    
    if(len(v) == 0):
        return np.nan, np.nan, np.nan, np.nan, np.nan

    #v_norm, v_unit = normalise_parameter(v, 'V')
    #c_norm, c_unit = normalise_parameter(c, 'F')

    v_fb1, v_fb2, c_acc, c_inv, t_ox, n_ox, a_acc, b_acc, v_acc, a_dep, b_dep, v_dep, a_inv, b_inv, v_inv,  spl_dev, status = analyse_mos(v, c)
    lbl = assign_label(path, test)
    c_acc_m = np.mean(c_acc)

    fit_acc = [a_acc*x + b_acc for x in v]
    fit_dep = [a_dep*i+b_dep for i in v]
    annotate = 'V$_{{fb}}$: {} V (via intersection)\nV$_{{fb}}$: {} V (via inflection)\n\nt$_{{ox}}$: {} m\nn$_{{ox}}$: {} cm$^{{-2}}$'.format(round(v_fb2,2), round(v_fb1,2), round(t_ox, 2), round(n_ox, 2))

    x_loc = 0.15
    y_loc = 0.145
    
    if plotResults:
        fig, ax = plt.subplots(1,1)
        #plt.ylim(0, 100)
        fit_curve(ax, v, fit_acc, fit_dep)
        plt.axvline(x=v_fb2, color='black', linestyle='dashed')
        plot_curve(ax, v, c, 'CV Curve', 'Voltage [V]', 'Capacitance [F]', lbl, annotate, x_loc, y_loc)

    if printResults:
        print('%s: \tMOS: v_fb2: %.2e V\tc_acc: %.2e F\tt_ox: %.2e nm\tn_ox: %.2e cm^-2' % (lbl, v_fb2, c_acc_m, t_ox*1e3, n_ox))
         
    return v_fb1, v_fb2, t_ox, n_ox, c_acc_m

     
def analyse_gcd_data(path, plotResults=True, printResults=print_results):
    test = 'gcd'

    v = read_json_file(path, test, 'voltage')
    i_em = read_json_file(path, test, 'current_elm')
    i_src = read_json_file(path, test, 'current_vsrc')
    i_hvsrc = read_json_file(path, test, 'current_hvsrc')

    if(len(v) == 0):
        return np.nan, np.nan

    lbl = assign_label(path, test)

    i_surf, i_bulk, i_acc, i_dep, i_inv, v_acc, v_dep, v_inv, spl_dev, status = analyse_gcd(v,i_em)
    
    if plotResults:
        fig, ax = plt.subplots(1,1)
        plot_curve(ax, v, i_em, 'I-V Curve GCD', 'Voltage [V]', 'Current [{}]'.format("A"), lbl, '', 0, 0)

    if printResults:
        print('%s: \tGCD: i_surf: %.2e A\t i_bulk: %.2e A' % (lbl, i_surf, i_bulk))
  
    return i_surf, i_bulk



def analyse_fet_data(path, plotResults=True, printResults=print_results):
    test = 'fet'

    v = read_json_file(path, test, 'voltage')
    i_em = read_json_file(path, test, 'current_elm')
    i_src = read_json_file(path, test, 'current_vsrc')
    i_hvsrc = read_json_file(path, test, 'current_hvsrc')
    
    if(len(v) == 0):
        return np.nan

    v_th, a, b, spl_dev, status = analyse_fet(v, i_em)

    fit  = [a*i +b for i in v]

    lbl = assign_label(path, test)
    
    if plotResults:
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
        plt.legend(lns, labs, loc='upper left')
        plt.show()

    if printResults:
       print('%s: \tnFet: v_th: %.2e V' % (lbl, v_th))
 
    return v_th


def analyse_van_der_pauw_data(path, printResults=print_results, plotResults=True):
    test = 'van-der-pauw'

    v = read_json_file(path, test, 'voltage_vsrc')
    i = read_json_file(path, test, 'current')
    if(len(v) == 0):
        return np.nan, 0
    
    lbl = assign_label(path, test)
    lbl_vdp = assign_label(path, test, vdp=True)
    r_sheet, a, b, x_fit, spl_dev, status, r_value = analyse_van_der_pauw(i, v)
    if(abs(r_value) < 0.9):  # r_value is the correlation coefficient
        a = b = np.nan  # the quality of the fit is too bad, we don't want it
        r_sheet = np.nan
        r_value = 0
    else:
        fit = [a*x +b for x in x_fit]
        if plotResults:
            fig, ax = plt.subplots(1,1)
            fit_curve(ax, x_fit, fit, 0)
            plot_curve(ax, i, v, 'IV Curve', 'Current', 'Voltage', lbl, '', 0, 0)
        
        if printResults:
           print('%s: \tvdp: r_sheet: %.2e Ohm/sq, correlation: %.2e  %s' % (lbl, r_sheet, r_value, lbl_vdp))
 

    return r_sheet, r_value


def analyse_linewidth_data(path, r_sheet=np.nan, printResults=print_results, plotResults=True):
    test = 'linewidth'
    
    if path is None:
        return np.nan

    v = read_json_file(path, test, 'voltage_vsrc')
    i = read_json_file(path, test, 'current')

    lbl = assign_label(path, test)
    lbl_vdp = assign_label(path, test, vdp=True)
    t_line, a, b, x_fit, spl_dev, status = analyse_linewidth(i, v, r_sheet=r_sheet, cut_param=0.01, debug=0)

    fit = [a*x +b for x in x_fit]
    if plotResults:
        fig, ax = plt.subplots(1, 1)
        fit_curve(ax, x_fit, fit, 0)
        plot_curve(ax, i, v, 'IV Curve', 'Current [A]', 'Voltage [V]', lbl, '', 0, 0)

    if printResults:
        print('%s: \tLinewidth: %.2e um\t%s' % (lbl, t_line, lbl_vdp))
 
    return t_line
  

def analyse_cbkr_data(path, r_sheet=np.nan, printResults=print_results, plotResults=True):
    test = 'cbkr'
    
    if path is None:
        return np.nan

    v = read_json_file(path, test, 'voltage_vsrc')
    i = read_json_file(path, test, 'current')
    
    if len(v) == 0:
        return np.nan

    lbl = assign_label(path, test)
    lbl_vdp = assign_label(path, test, vdp=True)

    r_contact, a, b, x_fit, spl_dev, status = analyse_cbkr(i, v, r_sheet, cut_param=0.01, debug=0)
    fit = [a*x +b for x in x_fit]
    
    if plotResults:
        fig, ax = plt.subplots(1, 1)
        fit_curve(ax, x_fit, fit, 0)
        plot_curve(ax, i, v, 'IV Curve', 'Current [A]', 'Voltage [V]', lbl, '', 0, 0)
   
    if printResults:
       print('%s: \tcbkr: r_contact: %.2e Ohm\t%s' % (lbl, r_contact, lbl_vdp))
 

    return r_contact

   
def analyse_contact_data(path):
    test= 'contact'

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

    rho_sq, status = analyse_meander(i, v, debug=0)
    
    if print_results:
       print('%s: \tMeander: rho_sq: %.2e' % (lbl, rho_sq))
   

    return rho_sq



def analyse_breakdown_data(path, printResults=print_results, plotResults=True):
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

    v_bd, status = analyse_breakdown(v, i_elm, debug=0)
    
    if plotResults:
        fig, ax = plt.subplots(1,1)
        annotate = 'V$_{{bd}}$: {} V \n\nT$_{{avg}}$ : {} \u00B0C \nH$_{{avg}}$: {} $\%$ '.format(v_bd, round(np.mean(temp),2), round(np.mean(humidity),2))
        plot_curve(ax, v, i_elm_norm, 'IV Curve', 'Voltage [V]', 'Current [{}]'.format(i_elm_unit), lbl, annotate, x_loc, y_loc)


    if printResults:
       print('%s: \tBreakdown: v_bd: %.2e V' % (lbl, v_bd))
   
    return v_bd



def get_vdp_value(pathlist):
    """helper function to get best vdp result"""
    r_sheet = np.nan
    r_value = 0
    for f in pathlist:
        #print(f)
        rs, rv = analyse_van_der_pauw_data(f, printResults=False, plotResults=False)
        if(rv > r_value):  # we take the best value
            r_sheet = rs
            r_value = rv
    return r_sheet
    

    
    

def analyse_full_line_data(path):
    """
    This function is used to analyze various different measurements 
    and assemble one line (per fluteset) of results to be tabellized

    Parameters:
    path ... path to parent directory: subdirs for each measurment-set
    """
    dirs = glob.glob(os.path.join(path, "*"))
    dirs.sort()
    print("# serial                                 \t  vdp_poly/kOhm/sq       vdp_n/Ohm/sq     vdp_pstop/kOhm/sq   lw_n/um    lw_p2/um   lw_p4/um cbkr_poly/kOhm cbkr_n/Ohm")
    for dirc in dirs:
        for flute in ["PQCFlutesRight", "PQCFlutesLeft"]:
            label = dirc.split("/")[-1]            
            
            vdp_poly_f = get_vdp_value(find_all_files_from_path(dirc, "van_der_pauw", whitelist=[flute, "Polysilicon"], blacklist=["reverse"]))
            vdp_poly_r = get_vdp_value(find_all_files_from_path(dirc, "van_der_pauw", whitelist=[flute, "Polysilicon", "reverse"]))
            
            vdp_n_f = get_vdp_value(find_all_files_from_path(dirc, "van_der_pauw", whitelist=[flute, "n"], blacklist=["reverse"]))
            vdp_n_r = get_vdp_value(find_all_files_from_path(dirc, "van_der_pauw", whitelist=[flute, "n", "reverse"]))
            
            vdp_pstop_f = get_vdp_value(find_all_files_from_path(dirc, "van_der_pauw", whitelist=[flute, "P_stop"], blacklist=["reverse"]))
            vdp_pstop_r = get_vdp_value(find_all_files_from_path(dirc, "van_der_pauw", whitelist=[flute, "P_stop", "reverse"]))
            
            t_line_n = analyse_linewidth_data(find_all_files_from_path(dirc, "linewidth", whitelist=[flute, "n"], single=True), r_sheet=vdp_n_f, printResults=False, plotResults=False)
            t_line_pstop2 = analyse_linewidth_data(find_all_files_from_path(dirc, "linewidth", whitelist=[flute, "P_stop", "2_wire"], single=True), r_sheet=vdp_pstop_f, printResults=False, plotResults=False)
            t_line_pstop4 = analyse_linewidth_data(find_all_files_from_path(dirc, "linewidth", whitelist=[flute, "P_stop", "4_wire"], single=True), r_sheet=vdp_pstop_f, printResults=False, plotResults=False)
            
            r_contact_n = analyse_cbkr_data(find_all_files_from_path(dirc, "cbkr", whitelist=[flute, "n"], single=True), r_sheet=vdp_n_f, printResults=False, plotResults=False)
            r_contac_poly = analyse_cbkr_data(find_all_files_from_path(dirc, "cbkr", whitelist=[flute, "Polysilicon"], single=True), r_sheet=vdp_poly_f, printResults=False, plotResults=False)
            
            
            line = "{} {}  \t".format(label, flute)
            line += "{:8.2f} {:8.2f}    ".format(vdp_poly_f*1e-3, vdp_poly_r*1e-3)
            line += "{:8.2f} {:8.2f}    ".format(vdp_n_f, vdp_n_r)
            line += "{:8.2f} {:8.2f}    ".format(vdp_pstop_f*1e-3, vdp_pstop_r*1e-3)
            line += "{:8.2f} {:8.2f} {:8.2f}     ".format(t_line_n, t_line_pstop2, t_line_pstop4)
            line += "{:8.2f} {:8.2f}".format(r_contac_poly*1e-3, r_contact_n)
            
            
            
            print(line)
    
    print("")
    print("")
    print("# serial                                 \t fet       vdp_met-clov       vdp_p-cr-br/kOhm/sq  lw_cb/um  v_bd/V")
    print("# serial                                 \tv_th/V     ")
    for dirc in dirs:
        for flute in ["PQCFlutesRight", "PQCFlutesLeft"]:
            label = dirc.split("/")[-1]
            
            v_th = analyse_fet_data(find_all_files_from_path(dirc, "fet", whitelist=[flute,], single=True), printResults=False, plotResults=False)
            
            vdp_metclo_f = get_vdp_value(find_all_files_from_path(dirc, "van_der_pauw", whitelist=[flute, "metal", "clover"], blacklist=["reverse"]))
            vdp_metclo_r = get_vdp_value(find_all_files_from_path(dirc, "van_der_pauw", whitelist=[flute, "metal", "clover", "reverse"], blacklist=[]))
            
            vdp_p_cross_bridge_f = get_vdp_value(find_all_files_from_path(dirc, "van_der_pauw", whitelist=[flute, "P", "cross_bridge"], blacklist=["reverse"]))
            vdp_p_cross_bridge_r = get_vdp_value(find_all_files_from_path(dirc, "van_der_pauw", whitelist=[flute, "P", "cross_bridge", "reverse"]))
            t_line_p_cross_bridge = analyse_linewidth_data(find_all_files_from_path(dirc, "linewidth", whitelist=[flute, "P", "cross_bridge"], single=True), r_sheet=vdp_p_cross_bridge_f, printResults=False, plotResults=False)
            
            v_bd = analyse_breakdown_data(find_all_files_from_path(dirc, "breakdown", whitelist=[flute,], single=True), printResults=False, plotResults=False)
            
            line = "{} {}  \t".format(label, flute)
            line += "{:5.2f}  ".format(v_th)
            line += "{:9.2E}  {:9.2E}   ".format(vdp_metclo_f, vdp_metclo_r)
            
            line += "{:8.2f} {:8.2f}    ".format(vdp_p_cross_bridge_f*1e-3, vdp_p_cross_bridge_r*1e-3)
            line += "{:8.2f}".format(t_line_p_cross_bridge)
            
            line += "{:8.2f}".format(v_bd)
            print(line)
    
    print("")
    print("")
    print("# serial                                 \t                    mos                         fet                gcd             vdp_met-clov   vdp_p-cr-br/kOhm/sq  lw_cb/um")
    print("# serial                                 \t v_fb/V     t_ox/nm    n_ox/cm^-2 c_acc/F     v_th/V       i_surf/A    i_bulk/A")
    for dirc in dirs:
        for flute in ["PQCFlutesRight", "PQCFlutesLeft"]:
            label = dirc.split("/")[-1]
            
            v_fb1, v_fb2, t_ox, n_ox, c_acc_m = analyse_mos_data(find_all_files_from_path(dirc, "mos", whitelist=[flute,], single=True), printResults=False, plotResults=False)
            
            i_surf, i_bulk = analyse_gcd_data(find_all_files_from_path(dirc, "gcd", whitelist=[flute,], single=True), printResults=False, plotResults=False)
            
            
            line = "{} {}  \t".format(label, flute)
            line += "{:9.2E}  {:9.2E}  {:9.2E}  {:9.2E}   ".format(v_fb2, t_ox, n_ox, c_acc_m)
            line += "{:9.2E}  {:9.2E}    ".format(i_surf, i_bulk)

            print(line)
            
            

functions = {
        'iv': analyse_iv_data,
        'cv': analyse_cv_data,
        'fet': analyse_fet_data,
        'gcd': analyse_gcd_data,
        'mos': analyse_mos_data,
        'linewidth': analyse_linewidth_data,
        'van_der_pauw': analyse_van_der_pauw_data,
        'breakdown': analyse_breakdown_data,
        'cbkr': analyse_cbkr_data,
    }


def analyse_file(path, test):
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
    
    if args.test == 'full-line':
        analyse_full_line_data(args.path)
        plt.show()
        return 0
    elif args.test == 'all':
        tests = functions.keys()
    else:
        tests = [args.test]

    for test in tests:
        print(test)
        filedir = find_all_files_from_path(args.path, test)
        filedir = np.sort(filedir)
        for f in filedir:
         #  print(f)
           analyse_file(f, test)
        plt.show()

if __name__ =="__main__":
    main()
