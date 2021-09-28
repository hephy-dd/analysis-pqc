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
import os
import glob
import pdb

import matplotlib.pyplot as plt
import numpy as np

from analysis_pqc import *
from pqc_rawdata import PQC_RawData
from pqc_analysis_tools import *
from itertools import repeat

__all__ = [
    'AnalysisOptions',
    'analyse_iv_data',
    'analyse_cv_data',
    'analyse_mos_data',
    'analyse_gcd_data',
    'analyse_fet_data',
    'analyse_van_der_pauw_data',
    'analyse_linewidth_data',
    'analyse_cbkr_data',
    'analyse_contact_data',
    'analyse_meander_data',
    'analyse_breakdown_data',
    'analyse_capacitor_data'
]

NOT_MEASURED = np.inf


class AnalysisOptions:

    def __init__(self, plotImgBasedir=None, label=None):
        self.plotWindow = False
        self.print = False
        self.plotImgBasedir = plotImgBasedir
        self.label = label
        self.prefixOverride = None

        self.plot = self.plotWindow or self.plotImgBasedir is not None

    # this is only temprary, we can reuse that object using once
    # used e g for VdP where one analysis function is used more than once
    def pushPrefix(self, prefixOverride):
        self.prefixOverride = prefixOverride
        return self

    def popPrefix(self, prefix):
        if self.prefixOverride is not None:
            ret = self.prefixOverride
            self.prefixOverride = None
            return ret
        else:
            return prefix

    def peekPrefix(self, defaultPrefix):
        if self.prefixOverride:
            return self.prefixOverride
        else:
            return defaultPrefix

    def savePlot(self, defaultPrefix, fig):
        if self.plotWindow:
            plt.show()
        else:
            prefix = self.popPrefix(defaultPrefix).lower()
            filename = f"{prefix}_{self.label}.png"
            fig.savefig(os.path.join(self.plotImgBasedir, filename))
            plt.close()

    def plotTitle(self, defaultPrefix):
        prefix = self.peekPrefix(defaultPrefix)
        return (f"{prefix}: {self.label}").replace('_', ' ')


def analyse_iv_data(path, options=None):
    test = 'iv'
    if path is None:
        return NOT_MEASURED, NOT_MEASURED, None

    if options is None:
        options = AnalysisOptions()
    
    series = read_json_file(path).get('series')
    timestamp = series.get('timestamp', np.array([]))
    v = abs(series.get('voltage', np.array([])))
    i_elm = -series.get('current_elm', np.array([]))
    i = abs(series.get('current_hvsrc', np.array([])))
    temp = series.get('temperature_chuck', np.array([]))
    temp_box = series.get('temperature_box', np.array([]))
    humidity = series.get('humidity_box', np.array([]))
    
    if(len(v) == 0):
        return np.nan, np.nan, None

    x_loc = 0.5
    y_loc = 0.65

    v_max, i_max, i_800, i_600, status = analyse_iv(v, i)

    lbl = assign_label(path, test)

    if options.plot:
        annotate = 'I$_{max}$' + ': {:6.2f} uA @ {:6.2f} V'.format(i_max*1e6, max(v)) + '\n'
        annotate += 'I$_{600V}$' + ': {:8.2f} uA '.format(i_600*1e6) + '\n'

        annotate +='T$_{avg}$' + ': {:0.2f} $^\circ C$'.format(np.mean(temp)) + '\n'
        annotate +='rH$_{avg}$:' + '{:0.2f}'.format(np.mean(humidity)) + r'$\%$'

        fig,ax = plt.subplots(1,1)
        plot_curve(ax, v, i*1e6, options.plotTitle("diode IV??"), 'Reverse Bias Voltage / A', 'Current / uA', 'SMU', annotate, x_loc, y_loc)
        plot_curve(ax, v, i_elm*1e6, options.plotTitle("diode IV??"), 'Reverse Bias Voltage / A', 'Current / uA', 'Electrometer', annotate, x_loc, y_loc)
        options.savePlot("diode_iv??", fig)

    if options.print:
        print('%s:  IV:\ti_600: %.3f uA\ti_800: %.3f uA' % (lbl, i_600*1e6, i_800*1e6))

    meta = read_json_file(path).get('meta')
    start_timestamp=meta.get('start_timestamp').replace('T',' ')
    rawdata=PQC_RawData(path,test,meta,series)
    #convert relative timestamp to absolute timestamp
    timestamp_abs=np.array(list(map(rel_to_abs_timestamp,repeat(start_timestamp),timestamp)))
    rawdata.add_data({'len':len(v),
                      'timestamp':timestamp,
                      'timestamp_abs':timestamp_abs,
                      'v':v,#Volt
                      'i_elm':i_elm*1e9,#A to nA
                      'i':i*1e9,#A to nA
                      'temp':temp,#degC
                      'temp_box':temp_box,#degC
                      'humidity':humidity,#percent
                      'i_600':i_600*1e12,#A to pA
                      'i_800':i_800*1e12}) #A to pA
    return i_600, i_800, rawdata


def analyse_cv_data(path, options=None):
    test = 'cv'

    if path is None:
        return NOT_MEASURED, NOT_MEASURED, NOT_MEASURED, None

    if options is None:
        options = AnalysisOptions()

    series = read_json_file(path).get('series')
    timestamp = series.get('timestamp', np.array([]))
    v = abs(series.get('voltage_hvsrc', np.array([])))
    i = series.get('current_hvsrc', np.array([]))
    c = series.get('capacitance', np.array([]))
    c2 = series.get('capacitance2', np.array([]))
    r = series.get('resistance', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    temp_box = series.get('temperature_box', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    if(len(v) == 0) or (len(c) == 0):
        return np.nan, np.nan, np.nan, None

    x_loc = 0.3
    y_loc = 0.65

    inv_c2 = 1/c**2

    lbl = assign_label(path, test)
    if "Flute_1" in path:
        area = 1.56e-6  # m^2, quarter
        return -1, -1, -1   # this does not work anyway and creates a lot of plots
    elif "Flute_3" in path:
    	area = 6.25e-6  # m^2, half (but without rounded edges)
    else:
        area = 1
        print("WARNING: clould not determine flute number - area dependent values will be wrong!")

    v_dep1, v_dep2, rho, conc, a_rise, b_rise, v_rise, a_const, b_const, v_const, spl_dev, status = analyse_cv(v, c, area=area, cut_param= 0.008, carrier='holes')


    annotate = 'V$_{{fd}}}}$: {} V\n\nT$_{{avg}}$: {} \u00B0C\nH$_{{avg}}$: {}'.format(v_dep2, round(np.mean(temp),2), round(np.mean(humidity),2)) + r'$\%$'

    #fig1, ax1 = plt.subplots(1, 1)
    #plot_curve(ax1, v_norm, c_norm, 'CV curve', 'Voltage[{}]'.format(v_unit), 'Capacitance [{}]'.format(c_unit), lbl, annotate, x_loc, y_loc)
    if options.plot:
        fig, ax2 = plt.subplots(1,1)
        #ax2b = ax2.twinx()
        fit_curve(ax2, v_rise, a_rise * v_rise+ b_rise, color='ro')
        fit_curve(ax2, v_const, a_const * v_const+ b_const, color='kx')
        #fit_curve(ax2b, v_norm, spl_dev, color='mx')
        plot_curve(ax2, v, 1./c**2, options.plotTitle("CV??"), 'Bias Voltage / V', '1/C$^{2}$ / F$^{-2}$', lbl, '', 0, 0 )
        plt.axvline(x=v_dep2, color='green', linestyle='dashed')

        resstr = f"v_fd: {v_dep2:8.2f} V\n"
        resstr += f"bulk res: {rho * 0.1:8.2f} kOhm cm\n"
        resstr += f"doping conc: {conc * 1e-18:8.2f}*1E12cm^-3"
        fig.text(0.95, 0.33, resstr, bbox=dict(facecolor='deepskyblue', alpha=0.75), horizontalalignment='right', verticalalignment='top')

        options.savePlot("diodecv??", fig)

    if options.print:
        # print(f"{lbl}: CV: v_fd: {}")
        print('%s: \tCV: v_fd: %.2e V\trho: %.2e Ohm\tconc: %.2e cm^-3' % (lbl, v_dep2, rho, conc*1e-6))
    meta = read_json_file(path).get('meta')
    start_timestamp=meta.get('start_timestamp').replace('T',' ')
    rawdata=PQC_RawData(path,test,meta,series)
    #convert relative timestamp to absolute timestamp
    timestamp_abs=np.array(list(map(rel_to_abs_timestamp,repeat(start_timestamp),timestamp)))
    '''    
    rawdata.add_data({'len':len(v),
                      'timestamp':timestamp,
                      'timestamp_abs':timestamp_abs,
                      'v':v,#Volt
                      'i_elm':i_elm*1e9,#A to nA
                      'i':i*1e9,#A to nA
                      'temp':temp,#degC
                      'temp_box':temp_box,#degC
                      'humidity':humidity,#percent
                      'i_600':i_600*1e12,#A to pA
                      'i_800':i_800*1e12}) #A to pA
    '''
    return v_dep2, rho, conc, rawdata


def analyse_mos_data(path, options=None):
    test = 'mos'

    if path is None:
        return NOT_MEASURED, NOT_MEASURED, NOT_MEASURED, NOT_MEASURED, NOT_MEASURED, None

    if options is None:
        options = AnalysisOptions()

    series = read_json_file(path).get('series')
    timestamp = series.get('timestamp', np.array([]))
    v = series.get('voltage_hvsrc', np.array([]))
    i = series.get('current_hvsrc', np.array([]))
    c = series.get('capacitance', np.array([]))
    c2 = series.get('capacitance2', np.array([]))
    r = series.get('resistance', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    temp_box = series.get('temperature_box', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    if(len(v) == 0):
        return np.nan, np.nan, np.nan, np.nan, np.nan, None

    #v_norm, v_unit = normalise_parameter(v, 'V')
    #c_norm, c_unit = normalise_parameter(c, 'F')
    try:
        v_fb1, v_fb2, c_acc, c_inv, t_ox, n_ox, a_acc, b_acc, v_acc, a_dep, b_dep, v_dep, a_inv, b_inv, v_inv,  spl_dev, status = analyse_mos(v, c)
        lbl = assign_label(path, test)
        c_acc_m = np.mean(c_acc)

    except:
        return np.nan, np.nan, np.nan, np.nan, np.nan, None
    annotate = 'V$_{{fb}}$: {} V (via intersection)\nV$_{{fb}}$: {} V (via inflection)\n\nt$_{{ox}}$: {} m\nn$_{{ox}}$: {} cm$^{{-2}}$'.format(round(v_fb2,2), round(v_fb1,2), round(t_ox, 2), round(n_ox, 2))

    x_loc = 0.15
    y_loc = 0.145

    if options.plot:
        fig, ax = plt.subplots(1,1)
        plot_curve(ax, v, c*1e12, options.plotTitle("MOS"), 'Bias Voltage / V', 'Capacitance / pF')
        #plt.ylim(0, 100)
        fit_curve(ax, v_dep, np.array([a_dep * v + b_dep for v in v_dep]) * 1e12)
        fit_curve(ax, v_acc, np.array([a_acc * v + b_acc for v in v_acc]) * 1e12)
        plt.axvline(x=v_fb2, color='black', linestyle='dashed')


        resstr = f"v_fb: {v_fb2:8.2f} V\n"
        fig.text(0.95, 0.33, resstr, bbox=dict(facecolor='deepskyblue', alpha=0.75), horizontalalignment='right', verticalalignment='top')

        options.savePlot("mos", fig)

    if options.print:
        print('%s: \tMOS: v_fb2: %.2e V\tc_acc: %.2e F\tt_ox: %.3e um\tn_ox: %.2e cm^-2' % (lbl, v_fb2, c_acc_m, t_ox, n_ox))
        
    meta = read_json_file(path).get('meta')
    start_timestamp=meta.get('start_timestamp').replace('T',' ')
    rawdata=PQC_RawData(path,test,meta,series)
    #convert relative timestamp to absolute timestamp
    timestamp_abs=np.array(list(map(rel_to_abs_timestamp,repeat(start_timestamp),timestamp)))
    rawdata.add_data({'len':len(v),
                      'timestamp':timestamp,
                      'timestamp_abs':timestamp_abs,
                      'v':v,#Volt
                      'i':i*1e9,#A to nA
                      'c':c*1e12,#F to pF
                      'r':r*1e-6,#Ohm to Mohm
                      'temp':temp,#degC
                      'temp_box':temp_box,#degC
                      'humidity':humidity,#percent
                      '_':v_fb1,
                      'v_fb2':v_fb2,#V
                      't_ox':t_ox*1e3,#um to nm
                      'n_ox':n_ox*1e-10,#1e10cm^-2
                      'c_acc_m':c_acc_m*1e12#F to pF
                      })
    
    return v_fb1, v_fb2, t_ox, n_ox, c_acc_m, rawdata


def analyse_gcd_data(path, options=None):
    test = 'gcd'

    if path is None:
        return NOT_MEASURED, NOT_MEASURED, None

    if options is None:
        options = AnalysisOptions()

    series = read_json_file(path).get('series')
    timestamp = series.get('timestamp', np.array([]))
    v = series.get('voltage', np.array([]))
    i_em = series.get('current_elm', np.array([]))
    i_src = series.get('current_vsrc', np.array([]))
    i_hvsrc = series.get('current_hvsrc', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    temp_box = series.get('temperature_box', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    if(len(v) < 3) or (len(i_em) < 3):
        return np.nan, np.nan, None


    lbl = assign_label(path, test)

    gcd_result = analyse_gcd(v,i_em, maxreldev=0.03)

    if options.plot:
        fig, ax = plt.subplots(1,1)
        plot_curve(ax, v, i_em*1e12, options.plotTitle("GCD??"), 'Gate Voltage / V', 'Leakage Current / {}'.format("pA"))
        try:
            fit_curve(ax, gcd_result.v_acc, gcd_result.i_acc*1e12, color='r')
            fit_curve(ax, gcd_result.v_dep, gcd_result.i_dep*1e12, color='k')
            fit_curve(ax, gcd_result.v_inv, gcd_result.i_inv*1e12, color='m')
        except (ValueError, ):
            pass
        resstr = "i_surf:"+" {:8.2f} pA\n".format(gcd_result.i_surf*1e12)
        resstr += "i_bulk:"+" {:8.2f} pA\n".format(gcd_result.i_bulk*1e12)
        resstr += "i_acc_relstd:"+" {:8.3f}%\n".format(gcd_result.i_acc_relstd*100)
        resstr += "i_dep_relstd:"+" {:8.3f}%\n".format(gcd_result.i_dep_relstd*100)
        resstr += "i_inv_relstd:"+" {:8.3f}%".format(gcd_result.i_inv_relstd*100)
        fig.text(0.95, 0.33, resstr, bbox=dict(facecolor='deepskyblue', alpha=0.75), horizontalalignment='right', verticalalignment='top')

        options.savePlot("fet", fig)

    if options.print:
        print('%s: \tGCD: i_surf: %.2e A\t i_bulk: %.2e A' % (lbl, gcd_result.i_surf, gcd_result.i_bulk))

    meta = read_json_file(path).get('meta')
    start_timestamp=meta.get('start_timestamp').replace('T',' ')
    rawdata=PQC_RawData(path,test,meta,series)
    #convert relative timestamp to absolute timestamp
    timestamp_abs=np.array(list(map(rel_to_abs_timestamp,repeat(start_timestamp),timestamp)))
    '''    
    rawdata.add_data({'len':len(v),
                      'timestamp':timestamp,
                      'timestamp_abs':timestamp_abs,
                      'v':v,#Volt
                      'i_elm':i_elm*1e9,#A to nA
                      'i':i*1e9,#A to nA
                      'temp':temp,#degC
                      'temp_box':temp_box,#degC
                      'humidity':humidity,#percent
                      'i_600':i_600*1e12,#A to pA
                      'i_800':i_800*1e12}) #A to pA
    '''
    return gcd_result.i_surf, gcd_result.i_bulk, rawdata


def analyse_fet_data(path, options=None):
    test = 'fet'

    if path is None:
        return NOT_MEASURED, None

    if options is None:
        options = AnalysisOptions()

    series = read_json_file(path).get('series')
    timestamp = series.get('timestamp', np.array([]))
    v = series.get('voltage', np.array([]))
    i_em = series.get('current_elm', np.array([]))
    i_src = series.get('current_vsrc', np.array([]))
    i_hvsrc = series.get('current_hvsrc', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    temp_box = series.get('temperature_box', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    if(len(v) < 3) or (len(i_em) < 3):
        return np.nan, None

    iz_em = i_em - i_em[0]

    v_th, a, b, spl_dev, status = analyse_fet(v, i_em)

    fit  = np.array([a*i +b for i in v])

    lbl = assign_label(path, test)

    if options.plot:
        fig,ax1 = plt.subplots()
        plt.title( options.plotTitle("FET"))
        lns1a = ax1.plot(v,iz_em*1e6, ls='', marker='s', ms=3, color='tab:green', label='transfer characteristics - shifted')
        lns1 = ax1.plot(v,i_em*1e6, ls='', marker='s', ms=3, label='transfer characteristics')

        ax1.set_xlabel('V$_{GS}$ [V]')
        ax1.set_ylabel('I$_{D}$ [uA]')
        ax1.set_ylabel(r'I$_\mathrm{D}$ [uA]')
        ax2 = ax1.twinx()
        lns2 = ax2.plot(v, spl_dev*1e6, ls=' ', marker='s', ms=3, color='tab:orange', label="transconductance")
        ax2.tick_params(axis='y', labelcolor='tab:orange')
        ax2.set_ylabel(r'g$_\mathrm{m}$ [S]', color='tab:orange')
        lns3 = ax1.plot(v, fit*1e6, '--r', label="tangent")
        lns4 = ax1.plot([v_th, v_th], [-3, 3], '--k', label="tangent zero crossing")
        lns = lns1 + lns1a + lns2 + lns3 + lns4
        labs = [l.get_label() for l in lns]

        plt.legend(lns, labs, loc='upper left')
        ax1.grid(linestyle='dotted')
        resstr = "V$_{th}$:"+" {:8.2f} V".format(v_th)
        fig.text(0.85, 0.85, resstr, bbox=dict(facecolor='deepskyblue', alpha=0.75), horizontalalignment='right', verticalalignment='top')

        options.savePlot("fet", fig)

    if options.print:
       print('%s: \tnFet: v_th: %.2e V' % (lbl, v_th))

    meta = read_json_file(path).get('meta')
    start_timestamp=meta.get('start_timestamp').replace('T',' ')
    rawdata=PQC_RawData(path,test,meta,series)
    timestamp_abs=np.array(list(map(rel_to_abs_timestamp,repeat(start_timestamp),timestamp)))
    rawdata.add_data({'len':len(v),
                      'timestamp':timestamp,
                      'timestamp_abs':timestamp_abs,
                      'v':v,#Volt
                      'i_elm':i_em*1e9,#A to nA
                      'temp':temp,#degC
                      'temp_box':temp_box,#degC
                      'humidity':humidity,#percent
                      'v_th':v_th#Volt
                      })
    return v_th, rawdata


def analyse_van_der_pauw_data(path, options=None, min_correlation=0.99):
    test = 'van-der-pauw'

    if path is None:
        options.popPrefix("-")
        return NOT_MEASURED, None

    if options is None:
        options = AnalysisOptions()

    series = read_json_file(path).get('series')
    timestamp = series.get('timestamp', np.array([]))
    v = series.get('voltage_vsrc', np.array([]))
    i = series.get('current', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    temp_box = series.get('temperature_box', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    if len(v) <= 3 or len(i) <= 3:
        options.popPrefix("-")
        return np.nan, None

    lbl = assign_label(path, test)
    lbl_vdp = assign_label(path, test, vdp=True)
    r_sheet, a, b, x_fit, spl_dev, status, r_value = analyse_van_der_pauw(i, v)

    if(abs(r_value) < min_correlation):  # r_value is the correlation coefficient
        r_sheet = np.nan

    fit = np.array([a*x +b for x in x_fit])
    if options.plot:

        if r_sheet > 1000.:
            resstr = "R$_{sheet}$:"+" {:8.2f} kOhm/sq\n".format(r_sheet*1e-3)
        elif r_sheet < 1.:
            resstr = "R$_{sheet}$:"+" {:8.2f} mOhm/sq\n".format(r_sheet*1e3)
        else:
            resstr = "R$_{sheet}$:"+" {:8.2f} Ohm/sq\n".format(r_sheet)
        resstr += "correlation:"+" {:5.3f}".format(r_value)


        fig, ax = plt.subplots(1,1)
        fit_curve(ax, x_fit*1e6, fit*1e3, 0)
        plot_curve(ax, i*1e6, v*1e3, options.plotTitle("VdP ???"), 'Current / uA', 'Voltage / mV', '', resstr, 0.9, 0.3)
        options.savePlot("vdp___", fig)

    if options.print:
       #print('%s: \tvdp: r_sheet: %.2e Ohm/sq, raw: %.2e Ohm, correlation: %.2e  %s' % (lbl, r_sheet, a, r_value, lbl_vdp))  # lbl_vdp
        print('%s: \tvdp: r_sheet: %.2e Ohm/sq, raw: %.2e Ohm, correlation: %.2e  %s' % (lbl, r_sheet, a, r_value, lbl_vdp))  # lbl_vdp

    meta = read_json_file(path).get('meta')
    start_timestamp=meta.get('start_timestamp').replace('T',' ')
    rawdata=PQC_RawData(path,test,meta,series)
    #convert relative timestamp to absolute timestamp
    timestamp_abs=np.array(list(map(rel_to_abs_timestamp,repeat(start_timestamp),timestamp)))
    rawdata.add_data({'len':len(v),
                      'timestamp':timestamp,
                      'timestamp_abs':timestamp_abs,
                      'v':v,#Volt
                      'i':i*1e9,#A to nA
                      'temp':temp,#degC
                      'temp_box':temp_box,#degC
                      'humidity':humidity,#percent
                      'r_sheet':r_sheet,#Ohm/Sq
                      'raw':a#Ohm
    })

    return r_sheet, rawdata


def analyse_linewidth_data(path, r_sheet=np.nan, options=None, min_correlation=0.9):
    test = 'linewidth'

    if path is None:
        options.popPrefix("-")
        return NOT_MEASURED, None

    if options is None:
        options = AnalysisOptions()

    series = read_json_file(path).get('series')
    timestamp = series.get('timestamp', np.array([]))
    v = series.get('voltage_vsrc', np.array([]))
    i = series.get('current', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    temp_box = series.get('temperature_box', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    if len(v) <= 3 or len(i) <= 3:
        options.popPrefix("-")
        return np.nan, None

    lbl = assign_label(path, test)
    lbl_vdp = assign_label(path, test, vdp=True)
    t_line, a, b, x_fit, spl_dev, r_value, status = analyse_linewidth(i, v, r_sheet=r_sheet, cut_param=-1., min_correlation=min_correlation, debug=0)

    fit = np.array([a*x +b for x in x_fit])
    if options.plot:
        if a > 1000.:
            resstr = "R:"+" {:8.2f} kOhm\n".format(a*1e-3)
        else:
            resstr = "R:"+" {:8.2f} Ohm\n".format(a)
        resstr += "lw:"+" {:5.3f} um\n".format(t_line)
        resstr += "correlation:"+" {:5.3f}".format(r_value)


        fig, ax = plt.subplots(1,1)
        fit_curve(ax, x_fit*1e6, fit*1e3, 0)
        plot_curve(ax, i*1e6, v*1e3, options.plotTitle("VdP ???"), 'Current / uA', 'Voltage / mV', '', resstr, 0.9, 0.3)
        options.savePlot("vdp___", fig)


    if options.print:
        print('%s: \tLinewidth: %.2e um\t%s' % (lbl, t_line, lbl_vdp))

    meta = read_json_file(path).get('meta')
    start_timestamp=meta.get('start_timestamp').replace('T',' ')
    rawdata=PQC_RawData(path,test,meta,series)
    #convert relative timestamp to absolute timestamp
    timestamp_abs=np.array(list(map(rel_to_abs_timestamp,repeat(start_timestamp),timestamp)))
    '''    
    rawdata.add_data({'len':len(v),
                      'timestamp':timestamp,
                      'timestamp_abs':timestamp_abs,
                      'v':v,#Volt
                      'i_elm':i_elm*1e9,#A to nA
                      'i':i*1e9,#A to nA
                      'temp':temp,#degC
                      'temp_box':temp_box,#degC
                      'humidity':humidity,#percent
                      'i_600':i_600*1e12,#A to pA
                      'i_800':i_800*1e12}) #A to pA
    '''
    
    return t_line, rawdata


def analyse_cbkr_data(path, r_sheet=np.nan, options=None, min_correlation=0.95):
    test = 'cbkr'

    if path is None:
        options.popPrefix("-")
        return NOT_MEASURED, None

    if options is None:
        options = AnalysisOptions()

    series = read_json_file(path).get('series')
    timestamp = series.get('timestamp', np.array([]))
    v = series.get('voltage_vsrc', np.array([]))
    i = series.get('current', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    temp_box = series.get('temperature_box', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    if len(v) <= 3 or len(i) <= 3:
        return np.nan, None

    lbl = assign_label(path, test)
    lbl_vdp = assign_label(path, test, vdp=True)

    r_contact, a, b, x_fit, spl_dev, r_value, status = analyse_cbkr(i, v, r_sheet, cut_param=0.01, debug=0)
    x_fit = np.array(x_fit)
    fit = np.array([a*x +b for x in x_fit])

    if options.plot:
        if r_contact > 1e3:
            resstr = "R:"+" {:8.2f} kOhm\n".format(r_contact*1e-3)
        else:
            resstr = "R:"+" {:8.2f} Ohm\n".format(r_contact)
        resstr += "correlation:"+" {:5.3f}".format(r_value)

        fig, ax = plt.subplots(1,1)
        fit_curve(ax, x_fit*1e6, fit*1e3, 0)
        #fit_curve(ax, x_fit, fit, 0)
        plot_curve(ax, i*1e6, v*1e3, options.plotTitle("CKBK ???"), 'Current / uA', 'Voltage / mV', '', resstr, 0.9, 0.3)
        options.savePlot("cbkr___", fig)

    if options.print:
       print('%s: \tcbkr: r_contact: %.2e Ohm\t%s' % (lbl, r_contact, lbl_vdp))

    meta = read_json_file(path).get('meta')
    start_timestamp=meta.get('start_timestamp').replace('T',' ')
    rawdata=PQC_RawData(path,test,meta,series)
    #convert relative timestamp to absolute timestamp
    timestamp_abs=np.array(list(map(rel_to_abs_timestamp,repeat(start_timestamp),timestamp)))
    '''    
    rawdata.add_data({'len':len(v),
                      'timestamp':timestamp,
                      'timestamp_abs':timestamp_abs,
                      'v':v,#Volt
                      'i_elm':i_elm*1e9,#A to nA
                      'i':i*1e9,#A to nA
                      'temp':temp,#degC
                      'temp_box':temp_box,#degC
                      'humidity':humidity,#percent
                      'i_600':i_600*1e12,#A to pA
                      'i_800':i_800*1e12}) #A to pA
    '''

    return r_contact, rawdata


def analyse_contact_data(path, options=None, min_correlation=0.95):
    test= 'contact'

    if path is None:
        options.popPrefix("-")
        return NOT_MEASURED, None

    if options is None:
        options = AnalysisOptions()

    series = read_json_file(path).get('series')
    timestamp = series.get('timestamp', np.array([]))
    v = series.get('voltage_vsrc', np.array([]))
    i = series.get('current', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    temp_box = series.get('temperature_box', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    if len(v) <= 3 or len(i) <= 3:
        options.popPrefix("-")
        return np.nan, None

    lbl = assign_label(path, test)
    r_contact, a, b, x_fit, spl_dev, status, r_value = analyse_contact(i, v, cut_param=0.01, debug=0)

    if abs(r_value) < min_correlation:
        r_contact = np.nan

    x_fit = np.array(x_fit)
    fit = np.array([a*x +b for x in x_fit])

    if options.plot:
        if r_contact > 1e3:
            resstr = "R:"+" {:8.2f} kOhm\n".format(r_contact*1e-3)
        else:
            resstr = "R:"+" {:8.2f} Ohm\n".format(r_contact)
        resstr += "correlation:"+" {:5.3f}".format(r_value)

        fig, ax = plt.subplots(1,1)
        fit_curve(ax, x_fit*1e6, fit*1e3, 0)
        plot_curve(ax, i*1e6, v*1e3, options.plotTitle("Contact ???"), 'Current / uA', 'Voltage / mV', '', resstr, 0.9, 0.3)
        options.savePlot("contact___", fig)

    if options.print:
       print('%s: \tcontact: r_contact: %.2e Ohm, r_value: %.2f' % (lbl, r_contact, r_value))

    meta = read_json_file(path).get('meta')
    start_timestamp=meta.get('start_timestamp').replace('T',' ')
    rawdata=PQC_RawData(path,test,meta,series)
    #convert relative timestamp to absolute timestamp
    timestamp_abs=np.array(list(map(rel_to_abs_timestamp,repeat(start_timestamp),timestamp)))
    '''    
    rawdata.add_data({'len':len(v),
                      'timestamp':timestamp,
                      'timestamp_abs':timestamp_abs,
                      'v':v,#Volt
                      'i_elm':i_elm*1e9,#A to nA
                      'i':i*1e9,#A to nA
                      'temp':temp,#degC
                      'temp_box':temp_box,#degC
                      'humidity':humidity,#percent
                      'i_600':i_600*1e12,#A to pA
                      'i_800':i_800*1e12}) #A to pA
    '''

    return r_contact, rawdata


def analyse_meander_data(path, options=None, min_correlation=0.99):
    test = 'meander'

    if path is None:
        options.popPrefix("-")
        return NOT_MEASURED, None

    if options is None:
        options = AnalysisOptions()

    series = read_json_file(path).get('series')
    timestamp = series.get('timestamp', np.array([]))
    v = series.get('voltage_vsrc', np.array([]))
    i = series.get('current', np.array([]))
    if (len(v) == 0):  # the polisilicon resistor uses a v-source and not a current source
        v = series.get('voltage', np.array([]))
        i = series.get('current_elm', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    temp_box = series.get('temperature_box', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    if len(v) <= 3 or len(i) <= 3:
        options.popPrefix("-")
        return np.nan, None

    lbl = assign_label(path, test)

    r, status, r_value = analyse_meander(i, v)

    if (r_value < min_correlation):
        r = np.nan

    #fit = np.array([r*x for x in x_fit])
    if options.plot:
        if r > 1e6:
            resstr = f"R: {r*1e-6:8.2f} MOhm\n"
        else:
            resstr = f"R: {r:8.2f} Ohm\n"
        resstr += "correlation:"+" {:5.3f}".format(r_value)


        fig, ax = plt.subplots(1,1)
        #fit_curve(ax, x_fit*1e6, fit*1e3, 0)
        plot_curve(ax, i*1e6, v*1e3, options.plotTitle("Meander ???"), 'Current / uA', 'Voltage / mV', '', resstr, 0.9, 0.3)
        options.savePlot("meander___", fig)

    if options.print:
        print(f"{lbl}: \tMeander: r: {r:.2e} r_value: {r_value:.2f}")

    meta = read_json_file(path).get('meta')
    start_timestamp=meta.get('start_timestamp').replace('T',' ')
    rawdata=PQC_RawData(path,test,meta,series)
    #convert relative timestamp to absolute timestamp
    timestamp_abs=np.array(list(map(rel_to_abs_timestamp,repeat(start_timestamp),timestamp)))
    '''    
    rawdata.add_data({'len':len(v),
                      'timestamp':timestamp,
                      'timestamp_abs':timestamp_abs,
                      'v':v,#Volt
                      'i_elm':i_elm*1e9,#A to nA
                      'i':i*1e9,#A to nA
                      'temp':temp,#degC
                      'temp_box':temp_box,#degC
                      'humidity':humidity,#percent
                      'i_600':i_600*1e12,#A to pA
                      'i_800':i_800*1e12}) #A to pA
    '''

    return r, rawdata


def analyse_breakdown_data(path, options=None):
    test = 'breakdown'

    if path is None:
        return NOT_MEASURED, None

    if options is None:
        options = AnalysisOptions()

    series = read_json_file(path).get('series')
    timestamp = series.get('timestamp', np.array([]))
    v = series.get('voltage', np.array([]))
    i = series.get('current_hvsrc', np.array([]))
    i_elm = series.get('current_elm', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    temp_box = series.get('temperature_box', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    i = np.array(i)

    lbl = assign_label(path, test)
    x_loc = 0.3
    y_loc = 0.5

    if len(v) == 0:
        return np.nan, None

    v_bd, status = analyse_breakdown(v, i, debug=0)

    if options.plot:
        fig, ax = plt.subplots(1,1)
        annotate = 'V$_{{bd}}$: {} V \n\nT$_{{avg}}$ : {} \u00B0C \nH$_{{avg}}$: {} $\%$ '.format(v_bd, round(np.mean(temp),2), round(np.mean(humidity),2))
        plot_curve(ax, v, i*1e9, 'Dielectric Breakdown', 'Voltage / V', 'Current / nA', lbl, annotate, x_loc, y_loc)
        options.savePlot("breakdown", fig)

    if options.print:
       print('%s: \tBreakdown: v_bd: %.2e V' % (lbl, v_bd))

    meta = read_json_file(path).get('meta')
    start_timestamp=meta.get('start_timestamp').replace('T',' ')
    rawdata=PQC_RawData(path,test,meta,series)
    #convert relative timestamp to absolute timestamp
    timestamp_abs=np.array(list(map(rel_to_abs_timestamp,repeat(start_timestamp),timestamp)))
    '''    
    rawdata.add_data({'len':len(v),
                      'timestamp':timestamp,
                      'timestamp_abs':timestamp_abs,
                      'v':v,#Volt
                      'i_elm':i_elm*1e9,#A to nA
                      'i':i*1e9,#A to nA
                      'temp':temp,#degC
                      'temp_box':temp_box,#degC
                      'humidity':humidity,#percent
                      'i_600':i_600*1e12,#A to pA
                      'i_800':i_800*1e12}) #A to pA
    '''

    return v_bd, rawdata


def analyse_capacitor_data(path, options=None):
    test = 'capacitor'

    if path is None:
        return NOT_MEASURED, NOT_MEASURED, NOT_MEASURED, None

    if options is None:
        options = AnalysisOptions()

    series = read_json_file(path).get('series')
    timestamp = series.get('timestamp', np.array([]))
    v = series.get('voltage_hvsrc', np.array([]))
    i = series.get('current_hvsrc', np.array([]))
    c = series.get('capacitance', np.array([]))
    r = series.get('resistance', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    temp_box = series.get('temperature_box', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    lbl = assign_label(path, test)
    x_loc = 0.3
    y_loc = 0.5

    if len(v) <= 3 or len(c) <= 3:
        return np.nan, np.nan, np.nan, None

    c_mean, c_median, d, status = analyse_capacitor(v, c, debug=0)

    if options.plot:
        pass

    if options.print:
       print('%s: \tCapacitance: %.2e F, ' % (lbl, c_median))

    meta = read_json_file(path).get('meta')
    start_timestamp=meta.get('start_timestamp').replace('T',' ')
    rawdata=PQC_RawData(path,test,meta,series)
    #convert relative timestamp to absolute timestamp
    timestamp_abs=np.array(list(map(rel_to_abs_timestamp,repeat(start_timestamp),timestamp)))    
    rawdata.add_data({'len':len(v),
                      'timestamp':timestamp,
                      'timestamp_abs':timestamp_abs,
                      'v':v,#Volt
                      'i':i*1e9,#A to nA
                      'c':c*1e12,#F to pF
                      'r':r*1e-6,#Ohm to Mohm
                      'temp':temp,#degC
                      'temp_box':temp_box,#degC
                      'humidity':humidity,#percent
                      'c_median':c_median*1e12,#F to pF
                      'd':d*1e3}) #um to nm(?)

    return c_mean, c_median, d, rawdata


def get_vdp_value(pathlist, printResults=False, plotResults=False):
    """helper function to get best vdp result"""
    r_sheet = np.nan
    r_value = 0
    rs = np.nan
    for f in pathlist:
        #print(f)
        rs = analyse_van_der_pauw_data(f, min_correlation=.7)

    return rs

def analyse_full_line_data(path):
    """
    This function is used to analyze various different measurements
    and assemble one line (per fluteset) of results to be tabellized

    Parameters:
    path ... path to parent directory: subdirs for each measurment-set
    """
    dirs = glob.glob(os.path.join(path, "*"))
    flutelist = ["PQCFlutesLeft"]
    #flutelist = ["PQCFlutesLeft", "PQCFlutesRight"]
    flutes = flutelist*len(dirs)
    dirs = dirs*len(flutelist)   # we need to double it for the two flutes
    dirs.sort()

    from pqc_resultset import PQC_resultset

    pqc_results = PQC_resultset(flutes)
    pqc_results.analyze(dirs)
    # pqc_results.prettyPrint()

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

FUNCTIONS = {
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
    'meander': analyse_meander_data,
    'capacitor': analyse_capacitor_data,
}


def analyse_file(path, test, show_plots=False):
    if test in FUNCTIONS:
        options = AnalysisOptions()
        options.print = True
        options.plot = show_plots
        options.plotWindow = show_plots
        r = FUNCTIONS.get(test)(path, options=options)
    else:
        raise ValueError(f"no such test: '{test}'")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    parser.add_argument('test')
    parser.add_argument('-p', dest='plot', action='store_true', help="plot results in window")
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.isdir(args.path):
        raise OSError(f"not a directory: {args.path}")

    if args.test == 'full-line':
        analyse_full_line_data(args.path)
        plt.show()
        return 0
    elif args.test == 'all':
        tests = FUNCTIONS.keys()
    else:
        tests = [args.test]

    for test in tests:
        print(test)
        filedir = find_all_files_from_path(args.path, None)
        filedir = np.sort(filedir)
        for f in filedir:
           analyse_file(f, test, show_plots=args.plot)
        plt.show()


if __name__ =="__main__":
    main()
