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


def analyse_iv_data(path, plotResults=True, printResults=print_results):
    test = 'iv'

    if path is None:
        return np.nan, np.nan

    series = read_json_file(path).get('series')
    v = abs(series.get('voltage', np.array()))
    i_tot = series.get('current_elm')
    i = abs(series.get('current_hvsrc', np.array()))
    temp = series.get('temperature_chuck')
    humidity = series.get('humidity_box')

    if not v:
        return np.nan, np.nan

    x_loc = 0.3
    y_loc = 0.65

    v_max, i_max, i_800, i_600, status = analyse_iv(v, i)

    lbl = assign_label(path, test)

    if plotResults:
        annotate = 'I$_{max}$' + ': {} {} @ {}{}'.format(round(i_max, 2), "A", max(v), "V") + '\n\nT$_{avg}$' + ': {:0.2f} $^\circ C$'.format(np.mean(temp)) \
          + '\n\n H$_{avg}$:' + '{:0.2f}'.format(np.mean(humidity)) + r'$\%$'

        fig,ax = plt.subplots(1,1)
        plot_curve(ax, v, i, 'IV Curve', 'Reverse Bias Voltage [V]', 'Current [A]', lbl, annotate, x_loc, y_loc)

    if printResults:
        print('%s:  IV:\ti_600: %.3f uA\ti_800: %.3f uA' % (lbl, i_600*1e6, i_800*1e6))

    return i_600, i_800


def analyse_cv_data(path, plotResults=True, printResults=print_results):
    test = 'cv'

    if path is None:
        return np.nan, np.nan, np.nan

    series = read_json_file(path).get('series')
    v = abs(series.get('voltage_hvsrc', np.array()))
    i = series.get('current_hvsrc')
    c = series.get('capacitance')
    c2 = series.get('capacitance2')
    r = series.get('resistance')
    temp = series.get('temperature_chuck')
    humidity = series.get('humidity_box')

    if not v:
        return np.nan, np.nan, np.nan

    x_loc = 0.3
    y_loc = 0.65

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
    if plotResults:
        fig2, ax2 = plt.subplots(1,1)
        #ax2b = ax2.twinx()
        fit_curve(ax2, v_rise, a_rise * v_rise+ b_rise, color='ro')
        fit_curve(ax2, v_const, a_const * v_const+ b_const, color='kx')
        #fit_curve(ax2b, v_norm, spl_dev, color='mx')
        plot_curve(ax2, v, 1./c**2, 'Full Depletion Voltage Estimation', 'Voltage[V]', '1/C$^{2}$ [F$^{-2}$]', lbl, '', 0, 0 )

    if printResults:
    	#print(f"{lbl}: CV: v_fd: {}")
        print('%s: \tCV: v_fd: %.2e V\trho: %.2e Ohm\tconc: %.2e cm^-3' % (lbl, v_dep2, rho, conc*1e-6))

    return v_dep2, rho, conc


def analyse_mos_data(path, plotResults=True, printResults=print_results):
    test = 'mos'

    series = read_json_file(path).get('series')
    v = series.get('voltage_hvsrc')
    i = series.get('current_hvsrc')
    c = series.get('capacitance')
    c2 = series.get('capacitance2')
    r = series.get('resistance')
    temp = series.get('temperature_chuck')
    humidity = series.get('humidity_box')

    if not v:
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
        print('%s: \tMOS: v_fb2: %.2e V\tc_acc: %.2e F\tt_ox: %.3e um\tn_ox: %.2e cm^-2' % (lbl, v_fb2, c_acc_m, t_ox, n_ox))

    return v_fb1, v_fb2, t_ox, n_ox, c_acc_m


def analyse_gcd_data(path, plotResults=True, printResults=print_results):
    test = 'gcd'

    series = read_json_file(path).get('series')
    v = series.get('voltage')
    i_em = series.get('current_elm')
    i_src = series.get('current_vsrc')
    i_hvsrc = series.get('current_hvsrc')

    if not v:
        return np.nan, np.nan

    lbl = assign_label(path, test)

    i_surf, i_bulk, i_acc, i_dep, i_inv, v_acc, v_dep, v_inv, spl_dev, status = analyse_gcd(v,i_em)

    if plotResults and not math.isnan(i_surf) and not math.isnan(i_bulk):
        fig, ax = plt.subplots(1,1)
        plot_curve(ax, v, i_em, 'I-V Curve GCD', 'Voltage [V]', 'Current [{}]'.format("A"), lbl, '', 0, 0)
        fit_curve(ax, v_acc, i_acc, color='r')
        fit_curve(ax, v_dep, i_dep, color='k')
        fit_curve(ax, v_inv, i_inv, color='m')

    if printResults:
        print('%s: \tGCD: i_surf: %.2e A\t i_bulk: %.2e A' % (lbl, i_surf, i_bulk))

    return i_surf, i_bulk



def analyse_fet_data(path, plotResults=True, printResults=print_results):
    test = 'fet'

    series = read_json_file(path).get('series')
    v = series.get('voltage')
    i_em = series.get('current_elm')
    i_src = series.get('current_vsrc')
    i_hvsrc = series.get('current_hvsrc')

    if not v:
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

    series = read_json_file(path).get('series')
    v = series.get('voltage_vsrc')
    i = series.get('current')

    if not v:
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

    series = read_json_file(path).get('series')
    v = series.get('voltage_vsrc')
    i = series.get('current')

    if not v:
        return np.nan

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

    series = read_json_file(path).get('series')
    v = series.get('voltage_vsrc')
    i = series.get('current')

    if not v:
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

    if path is None:
        return np.nan

    series = read_json_file(path).get('series')
    v = series.get('voltage_vsrc')
    i = series.get('current')

    if not v:
        return np.nan

    i_norm, i_unit = normalise_parameter(i, 'A')

    lbl = assign_label(path, test)
    r_contact, a, b, x_fit, spl_dev, status = analyse_contact(i, v, cut_param=0.01, debug=0)

    fit = [a*x+b for x in x_fit]

    if print_results:
       print('%s: \tcontact: r_contact: %.2e Ohm' % (lbl, r_contact))

    return r_contact


def analyse_meander_data(path):
    test = 'meander'

    if path is None:
        return np.nan

    series = read_json_file(path).get('series')
    v = series.get('voltage_vsrc')
    i = series.get('current')

    if not v:
        return np.nan

    i_norm, i_unit = normalise_parameter(i, 'A')

    lbl = assign_label(path, test)

    rho_sq, status = analyse_meander(i, v, debug=0)

    if print_results:
       print('%s: \tMeander: rho_sq: %.2e' % (lbl, rho_sq))


    return rho_sq



def analyse_breakdown_data(path, printResults=print_results, plotResults=True):
    test = 'breakdown'

    if path is None:
        return np.nan

    series = read_json_file(path).get('series')
    v = series.get('voltage')
    i = series.get('current_hvsrc')
    i_elm = series.get('current_elm')
    temp = series.get('temperature_chuck')
    humidity = series.get('humidity_box')

    if not v:
        return np.nan

    lbl = assign_label(path, test)
    x_loc = 0.3
    y_loc = 0.5

    v_bd, status = analyse_breakdown(v, i_elm, debug=0)

    if plotResults:
        fig, ax = plt.subplots(1,1)
        annotate = 'V$_{{bd}}$: {} V \n\nT$_{{avg}}$ : {} \u00B0C \nH$_{{avg}}$: {} $\%$ '.format(v_bd, round(np.mean(temp),2), round(np.mean(humidity),2))
        plot_curve(ax, v, i_elm, 'IV Curve', 'Voltage [V]', 'Current [A]', lbl, annotate, x_loc, y_loc)


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
    flutes = ["PQCFlutesRight", "PQCFlutesLeft"]*len(dirs)
    dirs = dirs*2   # we need to double it for the two flutes
    dirs.sort()

    labels = ["n/a"]*len(dirs)
    default=0
    vdp_poly_f = [default]*len(dirs)
    vdp_poly_r = [default]*len(dirs)
    vdp_n_f = [default]*len(dirs)
    vdp_n_r = [default]*len(dirs)
    vdp_pstop_f = [default]*len(dirs)
    vdp_pstop_r = [default]*len(dirs)
    t_line_n = [default]*len(dirs)
    t_line_pstop2 = [default]*len(dirs)
    t_line_pstop4 = [default]*len(dirs)
    r_contact_n = [default]*len(dirs)
    r_contac_poly = [default]*len(dirs)

    v_th = [default]*len(dirs)
    vdp_metclo_f = [default]*len(dirs)
    vdp_metclo_r = [default]*len(dirs)
    vdp_p_cross_bridge_f = [default]*len(dirs)
    vdp_p_cross_bridge_r = [default]*len(dirs)
    t_line_p_cross_bridge = [default]*len(dirs)
    v_bd = [default]*len(dirs)

    i600 = [default]*len(dirs)
    v_fd = [default]*len(dirs)
    rho = [default]*len(dirs)
    conc = [default]*len(dirs)

    v_fb1 = [default]*len(dirs)
    v_fb2 = [default]*len(dirs)
    t_ox = [default]*len(dirs)
    n_ox = [default]*len(dirs)
    c_acc_m = [default]*len(dirs)
    i_surf = [default]*len(dirs)
    i_surf05 = [default]*len(dirs)
    i_bulk05 = [default]*len(dirs)


    print("# serial                                 \t  vdp_poly/kOhm/sq       vdp_n/Ohm/sq     vdp_pstop/kOhm/sq   lw_n/um    lw_p2/um   lw_p4/um cbkr_poly/kOhm cbkr_n/Ohm")
    for i in range(0, len(dirs)):
        labels[i] = dirs[i].split("/")[-1]

        vdp_poly_f[i] = get_vdp_value(find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "Polysilicon"], blacklist=["reverse"]))
        vdp_poly_r[i] = get_vdp_value(find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "Polysilicon", "reverse"]))

        vdp_n_f[i] = get_vdp_value(find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "n"], blacklist=["reverse"]))
        vdp_n_r[i] = get_vdp_value(find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "n", "reverse"]))

        vdp_pstop_f[i] = get_vdp_value(find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "P_stop"], blacklist=["reverse"]))
        vdp_pstop_r[i] = get_vdp_value(find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "P_stop", "reverse"]))

        t_line_n[i] = analyse_linewidth_data(find_all_files_from_path(dirs[i], "linewidth", whitelist=[flutes[i], "n"], single=True), r_sheet=vdp_n_f[i], printResults=False, plotResults=False)
        t_line_pstop2[i] = analyse_linewidth_data(find_all_files_from_path(dirs[i], "linewidth", whitelist=[flutes[i], "P_stop", "2_wire"], single=True), r_sheet=vdp_pstop_f[i], printResults=False, plotResults=False)
        t_line_pstop4[i] = analyse_linewidth_data(find_all_files_from_path(dirs[i], "linewidth", whitelist=[flutes[i], "P_stop", "4_wire"], single=True), r_sheet=vdp_pstop_f[i], printResults=False, plotResults=False)

        r_contact_n[i] = analyse_cbkr_data(find_all_files_from_path(dirs[i], "cbkr", whitelist=[flutes[i], "n"], single=True), r_sheet=vdp_n_f[i], printResults=False, plotResults=False)
        r_contac_poly[i] = analyse_cbkr_data(find_all_files_from_path(dirs[i], "cbkr", whitelist=[flutes[i], "Polysilicon"], single=True), r_sheet=vdp_poly_f[i], printResults=False, plotResults=False)


        line = "{} {}  \t".format(labels[i], flutes[i])
        line += "{:8.2f} {:8.2f}    ".format(vdp_poly_f[i]*1e-3, vdp_poly_r[i]*1e-3)
        line += "{:8.2f} {:8.2f}    ".format(vdp_n_f[i], vdp_n_r[i])
        line += "{:8.2f} {:8.2f}    ".format(vdp_pstop_f[i]*1e-3, vdp_pstop_r[i]*1e-3)
        line += "{:8.2f} {:8.2f} {:8.2f}     ".format(t_line_n[i], t_line_pstop2[i], t_line_pstop4[i])
        line += "{:8.2f} {:8.2f}".format(r_contac_poly[i]*1e-3, r_contact_n[i])

        print(line)

    print("")
    print("")
    print("# serial                                 \t fet       vdp_met-clov       vdp_p-cr-br/kOhm/sq  lw_cb/um  v_bd/V    i600/uA    V_fd/V   rho/kOhm cm   d-conc/cm^-3")
    print("# serial                                 \tv_th/V     ")
    for i in range(0, len(dirs)):
        v_th[i] = analyse_fet_data(find_all_files_from_path(dirs[i], "fet", whitelist=[flutes[i],], single=True), printResults=False, plotResults=False)

        vdp_metclo_f[i] = get_vdp_value(find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "metal", "clover"], blacklist=["reverse"]))
        vdp_metclo_r[i] = get_vdp_value(find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "metal", "clover", "reverse"], blacklist=[]))

        vdp_p_cross_bridge_f[i] = get_vdp_value(find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "P", "cross_bridge"], blacklist=["reverse"]))
        vdp_p_cross_bridge_r[i] = get_vdp_value(find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "P", "cross_bridge", "reverse"]))
        t_line_p_cross_bridge[i] = analyse_linewidth_data(find_all_files_from_path(dirs[i], "linewidth", whitelist=[flutes[i], "P", "cross_bridge"], single=True), r_sheet=vdp_p_cross_bridge_f[i], printResults=False, plotResults=False)

        v_bd[i] = analyse_breakdown_data(find_all_files_from_path(dirs[i], "breakdown", whitelist=[flutes[i],], single=True), printResults=False, plotResults=False)

        # we want this for FLute_3 and not Flute_1
        i600[i], dummy = analyse_iv_data(find_all_files_from_path(dirs[i], "iv", whitelist=[flutes[i], "3"], single=True), printResults=False, plotResults=False)
        v_fd[i], rho[i], conc[i] = analyse_cv_data(find_all_files_from_path(dirs[i], "cv", whitelist=[flutes[i], "3"], single=True), printResults=False, plotResults=False)

        line = "{} {}  \t".format(labels[i], flutes[i])
        line += "{:5.2f}  ".format(v_th[i])
        line += "{:9.2E}  {:9.2E}   ".format(vdp_metclo_f[i], vdp_metclo_r[i])

        line += "{:8.2f} {:8.2f}    ".format(vdp_p_cross_bridge_f[i]*1e-3, vdp_p_cross_bridge_r[i]*1e-3)
        line += "{:8.2f}".format(t_line_p_cross_bridge[i])

        line += "{:8.2f}     ".format(v_bd[i])
        line += "{:9.2f}  {:9.2f}  {:7.2f}  {:9.2E}   ".format(i600[i]*1e6, v_fd[i], rho[i]*1e-1, conc[i]*1e-6)
        print(line)
    print("")
    print("")
    print("# serial                                 \t                    mos                        gcd             gcd05")
    print("# serial                                 \t v_fb/V    c_acc/pF   t_ox/um n_ox/1E10cm^-2 i_surf/pA  i_surf/pA   i_bulk/pA")
    for i in range(0, len(dirs)):
        v_fb1[i], v_fb2[i], t_ox[i], n_ox[i], c_acc_m[i] = analyse_mos_data(find_all_files_from_path(dirs[i], "mos", whitelist=[flutes[i],], single=True), printResults=False, plotResults=False)
        i_surf[i], dummy = analyse_gcd_data(find_all_files_from_path(dirs[i], "gcd", whitelist=[flutes[i],], single=True), printResults=False, plotResults=False)  # only i_surf valid
        i_surf05[i], i_bulk05[i] = analyse_gcd_data(find_all_files_from_path(dirs[i], "gcd05", whitelist=[flutes[i],], single=True), printResults=False, plotResults=False)  # for i_bulk

        line = "{} {}  \t".format(labels[i], flutes[i])
        line += "{:8.2f}    {:6.2f}    {:7.3f}  {:9.2f}     ".format(v_fb2[i], c_acc_m[i]*1e12, t_ox[i], n_ox[i]*1e-10)
        line += "{:8.2f}  {:8.2f}  {:8.2f}    ".format(i_surf[i]*1e12, i_surf05[i]*1e12, i_bulk05[i]*1e12)

        print(line)


    #outfile = open("outpit.dat",'w')
    #outfile.write("#serial\tflute\t")
    #outfile.write("#VdP\trev\t")
    #
    #
    #outfile.write("#---\t---\t")
    #outfile.write("#\tflute\t")
    #for i in range(0, len(dirs)):
    #    outfile.write(lebels[i]+"\t"+flutes[i]+"\t")
    #
    #outfile.close()


functions = {
    'iv': analyse_iv_data,
    'cv': analyse_cv_data,
    'fet': analyse_fet_data,
    'gcd': analyse_gcd_data,
    'gcd05': analyse_gcd_data,
    'mos': analyse_mos_data,
    'linewidth': analyse_linewidth_data,
    'van_der_pauw': analyse_van_der_pauw_data,
    'breakdown': analyse_breakdown_data,
    'cbkr': analyse_cbkr_data,
}
"""Mapping of available analysis function names."""

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
    parser.add_argument('path', help="path to PQC analysis JSON files")
    parser.add_argument('test', help="analysis function to run")
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
