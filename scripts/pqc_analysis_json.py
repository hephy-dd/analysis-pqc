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
        if self.prefixOverride:
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
            fig.savefig(os.path.join(self.plotImgBasedir, self.popPrefix(defaultPrefix).lower()+"_"+self.label+".png"))
            plt.close()
            
    def plotTitle(self, defaultPrefix):
        return (self.peekPrefix(defaultPrefix) + ': ' + self.label).replace('_', ' ')


def analyse_iv_data(path, analysisOptions=AnalysisOptions()):
    test = 'iv'
    if path is None:
        return np.nan, np.nan

    series = read_json_file(path).get('series')
    v = abs(series.get('voltage', np.array([])))
    i_elm = -series.get('current_elm', np.array([]))
    i = abs(series.get('current_hvsrc', np.array([])))
    temp = series.get('temperature_chuck', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    if(len(v) == 0):
        return np.nan, np.nan

    x_loc = 0.5
    y_loc = 0.65

    v_max, i_max, i_800, i_600, status = analyse_iv(v, i)

    lbl = assign_label(path, test)

    if analysisOptions.plot:
        annotate = 'I$_{max}$' + ': {:6.2f} uA @ {:6.2f} V'.format(i_max*1e6, max(v)) + '\n'
        annotate += 'I$_{600V}$' + ': {:8.2f} uA '.format(i_600*1e6) + '\n'
        
        annotate +='T$_{avg}$' + ': {:0.2f} $^\circ C$'.format(np.mean(temp)) + '\n'
        annotate +='rH$_{avg}$:' + '{:0.2f}'.format(np.mean(humidity)) + r'$\%$'

        fig,ax = plt.subplots(1,1)
        plot_curve(ax, v, i*1e6, analysisOptions.plotTitle("diode IV??"), 'Reverse Bias Voltage / A', 'Current / uA', 'SMU', annotate, x_loc, y_loc)
        plot_curve(ax, v, i_elm*1e6, analysisOptions.plotTitle("diode IV??"), 'Reverse Bias Voltage / A', 'Current / uA', 'Electrometer', annotate, x_loc, y_loc)
        analysisOptions.savePlot("diode_iv??", fig)

    if analysisOptions.print:
        print('%s:  IV:\ti_600: %.3f uA\ti_800: %.3f uA' % (lbl, i_600*1e6, i_800*1e6))

    return i_600, i_800


def analyse_cv_data(path, analysisOptions=AnalysisOptions()):
    test = 'cv'

    if path is None:
        return np.nan, np.nan, np.nan

    series = read_json_file(path).get('series')
    v = abs(series.get('voltage_hvsrc', np.array([])))
    i = series.get('current_hvsrc', np.array([]))
    c = series.get('capacitance', np.array([]))
    c2 = series.get('capacitance2', np.array([]))
    r = series.get('resistance', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    if(len(v) == 0) or (len(c) == 0):
        return np.nan, np.nan, np.nan

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
    if analysisOptions.plot:
        fig, ax2 = plt.subplots(1,1)
        #ax2b = ax2.twinx()
        fit_curve(ax2, v_rise, a_rise * v_rise+ b_rise, color='ro')
        fit_curve(ax2, v_const, a_const * v_const+ b_const, color='kx')
        #fit_curve(ax2b, v_norm, spl_dev, color='mx')
        plot_curve(ax2, v, 1./c**2, analysisOptions.plotTitle("CV??"), 'Bias Voltage / V', '1/C$^{2}$ / F$^{-2}$', lbl, '', 0, 0 )
        plt.axvline(x=v_dep2, color='green', linestyle='dashed')
        
        resstr = "v_fd:"+" {:8.2f} V\n".format(v_dep2)
        resstr += "bulk res:"+" {:8.2f} kOhm cm\n".format(rho*0.1)
        resstr += "doping conc:"+" {:8.2f}*1E12cm^-3".format(conc*1e-18)
        fig.text(0.95, 0.33, resstr, bbox=dict(facecolor='deepskyblue', alpha=0.75), horizontalalignment='right', verticalalignment='top')
        
        analysisOptions.savePlot("diodecv??", fig)

    if analysisOptions.print:
    	#print(f"{lbl}: CV: v_fd: {}")
        print('%s: \tCV: v_fd: %.2e V\trho: %.2e Ohm\tconc: %.2e cm^-3' % (lbl, v_dep2, rho, conc*1e-6))

    return v_dep2, rho, conc


def analyse_mos_data(path, analysisOptions=AnalysisOptions()):
    test = 'mos'
    
    if path is None:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    series = read_json_file(path).get('series')
    v = series.get('voltage_hvsrc', np.array([]))
    i = series.get('current_hvsrc', np.array([]))
    c = series.get('capacitance', np.array([]))
    c2 = series.get('capacitance2', np.array([]))
    r = series.get('resistance', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    humidity = series.get('humidity_box', np.array([]))

    if(len(v) == 0):
        return np.nan, np.nan, np.nan, np.nan, np.nan

    #v_norm, v_unit = normalise_parameter(v, 'V')
    #c_norm, c_unit = normalise_parameter(c, 'F')
    try:
        v_fb1, v_fb2, c_acc, c_inv, t_ox, n_ox, a_acc, b_acc, v_acc, a_dep, b_dep, v_dep, a_inv, b_inv, v_inv,  spl_dev, status = analyse_mos(v, c)
        lbl = assign_label(path, test)
        c_acc_m = np.mean(c_acc)

    except:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    annotate = 'V$_{{fb}}$: {} V (via intersection)\nV$_{{fb}}$: {} V (via inflection)\n\nt$_{{ox}}$: {} m\nn$_{{ox}}$: {} cm$^{{-2}}$'.format(round(v_fb2,2), round(v_fb1,2), round(t_ox, 2), round(n_ox, 2))

    x_loc = 0.15
    y_loc = 0.145

    if analysisOptions.plot:
        fig, ax = plt.subplots(1,1)
        plot_curve(ax, v, c*1e12, analysisOptions.plotTitle("MOS"), 'Bias Voltage / V', 'Capacitance / pF')
        #plt.ylim(0, 100)
        fit_curve(ax, v_dep, np.array([a_dep*v+b_dep for v in v_dep])*1e12)
        fit_curve(ax, v_acc, np.array([a_acc*v + b_acc for v in v_acc])*1e12)
        plt.axvline(x=v_fb2, color='black', linestyle='dashed')
        
        
        resstr = "v_fb:"+" {:8.2f} V\n".format(v_fb2)
        fig.text(0.95, 0.33, resstr, bbox=dict(facecolor='deepskyblue', alpha=0.75), horizontalalignment='right', verticalalignment='top')
        
        analysisOptions.savePlot("mos", fig)
    
    if analysisOptions.print:
        print('%s: \tMOS: v_fb2: %.2e V\tc_acc: %.2e F\tt_ox: %.3e um\tn_ox: %.2e cm^-2' % (lbl, v_fb2, c_acc_m, t_ox, n_ox))

    return v_fb1, v_fb2, t_ox, n_ox, c_acc_m


def analyse_gcd_data(path, analysisOptions=AnalysisOptions()):
    test = 'gcd'
    
    if path is None:
        return np.nan, np.nan

    series = read_json_file(path).get('series')
    v = series.get('voltage', np.array([]))
    i_em = series.get('current_elm', np.array([]))
    i_src = series.get('current_vsrc', np.array([]))
    i_hvsrc = series.get('current_hvsrc', np.array([]))

    if(len(v) == 0) or (len(i_em) == 0):
        return np.nan, np.nan

    lbl = assign_label(path, test)

    gcd_result = analyse_gcd(v,i_em, maxreldev=0.03)

    if analysisOptions.plot:
        fig, ax = plt.subplots(1,1)
        plot_curve(ax, v, i_em*1e12, analysisOptions.plotTitle("GCD??"), 'Gate Voltage / V', 'Leakage Current / {}'.format("pA"))
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
        
        analysisOptions.savePlot("fet", fig)
    
    if analysisOptions.print:
        print('%s: \tGCD: i_surf: %.2e A\t i_bulk: %.2e A' % (lbl, gcd_result.i_surf, gcd_result.i_bulk))

    return gcd_result.i_surf, gcd_result.i_bulk



def analyse_fet_data(path, analysisOptions=AnalysisOptions()):
    test = 'fet'

    if path is None:
        return np.nan
        
    series = read_json_file(path).get('series')
    v = series.get('voltage', np.array([]))
    i_em = series.get('current_elm', np.array([]))
    i_src = series.get('current_vsrc', np.array([]))
    i_hvsrc = series.get('current_hvsrc', np.array([]))
    
    iz_em = i_em - i_em[0]

    if(len(v) < 3) or (len(i_em) < 3):
        return np.nan

    v_th, a, b, spl_dev, status = analyse_fet(v, i_em)

    fit  = np.array([a*i +b for i in v])

    lbl = assign_label(path, test)

    if analysisOptions.plot:
        fig,ax1 = plt.subplots()
        plt.title( analysisOptions.plotTitle("FET"))
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
        lns = lns1+lns1a+lns2+lns3+lns4
        labs = [l.get_label() for l in lns]
        
        plt.legend(lns, labs, loc='upper left')
        ax1.grid(linestyle='dotted')
        resstr = "V$_{th}$:"+" {:8.2f} V".format(v_th)
        fig.text(0.85, 0.85, resstr, bbox=dict(facecolor='deepskyblue', alpha=0.75), horizontalalignment='right', verticalalignment='top')
        
        analysisOptions.savePlot("fet", fig)

    if analysisOptions.print:
       print('%s: \tnFet: v_th: %.2e V' % (lbl, v_th))

    return v_th


def analyse_van_der_pauw_data(path, analysisOptions=AnalysisOptions(), minCorrelation=0.99):
    test = 'van-der-pauw'
    
    if path is None:
        return np.nan

    series = read_json_file(path).get('series')
    v = series.get('voltage_vsrc', np.array([]))
    i = series.get('current', np.array([]))

    if(len(v) == 0):
        return np.nan

    lbl = assign_label(path, test)
    lbl_vdp = assign_label(path, test, vdp=True)
    r_sheet, a, b, x_fit, spl_dev, status, r_value = analyse_van_der_pauw(i, v)
    
    if(abs(r_value) < minCorrelation):  # r_value is the correlation coefficient
        r_sheet = np.nan

    fit = np.array([a*x +b for x in x_fit])
    if analysisOptions.plot:
        

        if r_sheet > 1000.:
            resstr = "R$_{sheet}$:"+" {:8.2f} kOhm/sq\n".format(r_sheet*1e-3)
        elif r_sheet < 1.:
            resstr = "R$_{sheet}$:"+" {:8.2f} mOhm/sq\n".format(r_sheet*1e3)
        else:
            resstr = "R$_{sheet}$:"+" {:8.2f} Ohm/sq\n".format(r_sheet)
        resstr += "correlation:"+" {:5.3f}".format(r_value)
        
        
        fig, ax = plt.subplots(1,1)
        fit_curve(ax, x_fit*1e6, fit*1e3, 0)
        plot_curve(ax, i*1e6, v*1e3, analysisOptions.plotTitle("VdP ???"), 'Current / uA', 'Voltage / mV', '', resstr, 0.9, 0.3)
        analysisOptions.savePlot("vdp___", fig)

    if analysisOptions.print:
       #print('%s: \tvdp: r_sheet: %.2e Ohm/sq, raw: %.2e Ohm, correlation: %.2e  %s' % (lbl, r_sheet, a, r_value, lbl_vdp))  # lbl_vdp
        print('%s: \tvdp: r_sheet: %.2e Ohm/sq, raw: %.2e Ohm, correlation: %.2e  %s' % (lbl, r_sheet, a, r_value, lbl_vdp))  # lbl_vdp

    return r_sheet


def analyse_linewidth_data(path, r_sheet=np.nan, printResults=print_results, plotResults=True):
    test = 'linewidth'

    if path is None:
        return np.nan

    series = read_json_file(path).get('series')
    v = series.get('voltage_vsrc', np.array([]))
    i = series.get('current', np.array([]))
    
    if len(v) < 3:
        return np.nan

    lbl = assign_label(path, test)
    lbl_vdp = assign_label(path, test, vdp=True)
    t_line, a, b, x_fit, spl_dev, status = analyse_linewidth(i, v, r_sheet=r_sheet, cut_param=-1., debug=0)

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
    v = series.get('voltage_vsrc', np.array([]))
    i = series.get('current', np.array([]))

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


def analyse_contact_data(path, minCorrelation=0.96, printResults=print_results, plotResults=True):
    test= 'contact'

    if path is None:
        return np.nan

    series = read_json_file(path).get('series')
    v = series.get('voltage_vsrc', np.array([]))
    i = series.get('current', np.array([]))
    
    if len(v) == 0:
        return np.nan

    lbl = assign_label(path, test)
    r_contact, a, b, x_fit, spl_dev, status, r_value = analyse_contact(i, v, cut_param=0.01, debug=0)
    
    if abs(r_value) < minCorrelation:
        return np.nan

    fit = [a*x+b for x in x_fit]

    if printResults:
       print('%s: \tcontact: r_contact: %.2e Ohm, r_value: %.2f' % (lbl, r_contact, r_value))

    return r_contact



def analyse_meander_data(path, printResults=print_results, plotResults=True):
    test = 'meander'

    if path is None:
        return np.nan

    series = read_json_file(path).get('series')
    
    v = series.get('voltage_vsrc', np.array([]))
    i = series.get('current', np.array([]))
    if (len(v) == 0):  # the polisilicon resistor uses a v-source and not a current source
        v = series.get('voltage', np.array([]))
        i = series.get('current_elm', np.array([]))

    if len(v) < 3:
        return np.nan
        
    lbl = assign_label(path, test)

    r, status, r_value = analyse_meander(i, v)
    
    if (r_value < 0.95):
        return np.nan

    if printResults:
       print('%s: \tMeander: r: %.2e r_value: %.2f' % (lbl, r, r_value))


    return r



def analyse_breakdown_data(path, printResults=print_results, plotResults=True):
    test = 'breakdown'

    if path is None:
        return np.nan

    series = read_json_file(path).get('series')
    v = series.get('voltage', np.array([]))
    i = series.get('current_hvsrc', np.array([]))
    i_elm = series.get('current_elm', np.array([]))
    temp = series.get('temperature_chuck', np.array([]))
    humidity = series.get('humidity_box', np.array([]))


    lbl = assign_label(path, test)
    x_loc = 0.3
    y_loc = 0.5
    
    if len(v) < 1:
        return np.nan

    v_bd, status = analyse_breakdown(v, i_elm, debug=0)

    if plotResults:
        fig, ax = plt.subplots(1,1)
        annotate = 'V$_{{bd}}$: {} V \n\nT$_{{avg}}$ : {} \u00B0C \nH$_{{avg}}$: {} $\%$ '.format(v_bd, round(np.mean(temp),2), round(np.mean(humidity),2))
        plot_curve(ax, v, i_elm, 'IV Curve', 'Voltage [V]', 'Current [A]', lbl, annotate, x_loc, y_loc)


    if printResults:
       print('%s: \tBreakdown: v_bd: %.2e V' % (lbl, v_bd))

    return v_bd
    
    
def analyse_capacitor_data(path, printResults=print_results, plotResults=True):
    test = 'capacitor'

    if path is None:
        return np.nan, np.nan, np.nan

    series = read_json_file(path).get('series')
    v = series.get('voltage_hvsrc', np.array([]))
    c = series.get('capacitance', np.array([]))


    lbl = assign_label(path, test)
    x_loc = 0.3
    y_loc = 0.5
    
    if len(v) < 1:
        return np.nan

    c_mean, c_median, d, status = analyse_capacitor(v, c, debug=0)

    if plotResults:
        pass


    if printResults:
       print('%s: \tCapacitance: %.2e F, ' % (lbl, c_median))

    return c_mean, c_median, d



def get_vdp_value(pathlist, printResults=False, plotResults=False):
    """helper function to get best vdp result"""
    r_sheet = np.nan
    r_value = 0
    rs = np.nan
    for f in pathlist:
        #print(f)
        rs = analyse_van_der_pauw_data(f, minCorrelation=.7)

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

    pqc_results = PQC_resultset(flutes)
    pqc_results.analyze(dirs)
    pqc_results.prettyPrint()
        

    
    


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
    'meander': analyse_meander_data,
    'capacitor': analyse_capacitor_data,
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
        filedir = find_all_files_from_path(args.path, None)
        filedir = np.sort(filedir)
        for f in filedir:
         #  print(f)
           analyse_file(f, test)
        plt.show()

if __name__ =="__main__":
    main()
