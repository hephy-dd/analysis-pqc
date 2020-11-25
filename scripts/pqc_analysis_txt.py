#!/usr/bin/env python3

import os

import numpy as np
import matplotlib.pyplot as plt

from optparse import OptionParser
from analysis_pqc import *

# TODO:
# Catch incomplete files

## Single Analysis Functions
## ------------------------------------

def analyse_iv_data(fn, debug=0):
    path = '/'.join(fn.split('/')[:-1])
    fn = fn.split('/')[-1]
    if debug:
        print(" -- Analysing diode IV for file:\n\t%s" % fn)

    ## read data
    dat = np.genfromtxt(path + '/' + fn, skip_header=28)
    t, v, i_tot, i =  dat[:, 0], dat[:, 1], dat[:, 2], dat[:, 3]
    lbl = '_'.join(fn.split('_')[1:9])

    ## analyse
    v_max, i_max, i_800, i_600, status = analyse_iv(v, i)

    ## plot
    x_loc = 0.5
    y_loc = 0.2
    plt.plot(v, i, ls=' ', marker='s', ms=3, label=lbl)
    plt.annotate('I$_\mathrm{max}$: %.2e A @ %d V\n\nT$_\mathrm{avg}$: %.1f $^\circ C$\nH$_\mathrm{avg}$: %.1f ' % \
        (i_max, v_max, np.mean(dat[:, 5]), np.mean(dat[:, 6])) + r'$\%$', (x_loc,y_loc), xycoords='figure fraction', color='tab:blue', \
        bbox=dict(facecolor='white', edgecolor='tab:blue', boxstyle='round,pad=0.5'))
    plt.legend(loc='upper left')
    plt.xlabel('voltage [V]')
    plt.ylabel('current [A]')
    plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    plt.tight_layout()
    plt.savefig(path + '/fig_diode_iv_' + fn[:-4] + '.pdf')
    plt.clf()

    return v_max, i_max, i_800, i_600, fn[:-4]



def analyse_cv_data(fn, debug=0):
    path = '/'.join(fn.split('/')[:-1])
    fn = fn.split('/')[-1]
    if debug:
        print(" -- Analysing diode CV for file:\n\t%s" % fn)

    ## read data
    dat = np.genfromtxt(path + '/' + fn, skip_header=15, encoding='utf-8')
    t, v, i_tot, c, c2, r =  dat[:, 0], dat[:, 1], dat[:, 2], dat[:, 3], dat[:, 4], dat[:, 5]
    lbl = '_'.join(fn.split('_')[1:9])

    ## analyse
    v_dep1, v_dep2, rho, conc, a_rise, b_rise, v_rise, a_const, b_const, v_const, spl_dev, status = analyse_cv(v, c)

    ## Fix me first
    ## plot
    # x_loc = 0.5
    # y_loc = 0.7
    # plt.plot(v, c, ls=' ', marker='s', ms=3, label=lbl)
    # plt.annotate('V$_\mathrm{fd}$: %.1f V (via intersection)\nV$_\mathrm{fd}$: %.1f V (via inflection)\n\nT$_\mathrm{avg}$: %.1f $^\circ C$\nH$_\mathrm{avg}$: %.1f ' % \
    #     (v_dep2, v_dep1, np.mean(dat[:, 7]), np.mean(dat[:, 8])) + r'$\%$', (x_loc,y_loc), xycoords='figure fraction', color='tab:blue', \
    #     bbox=dict(facecolor='white', edgecolor='tab:blue', boxstyle='round,pad=0.5'))
    # plt.legend(loc='upper left')
    # plt.xlabel('voltage [V]')
    # plt.ylabel('capactiance [F]')
    # plt.ylim([0, +1e-11])
    # plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    # plt.tight_layout()
    # plt.savefig(path + '/fig_diode_cv' + fn[:-4] + '.pdf')
    # plt.clf()
    #
    # x_loc = 0.2
    # y_loc = 0.2
    # plt.plot(v, 1/c**2, ls=' ', marker='s', ms=3, label=lbl)
    # plt.plot(v_rise, a_rise*v_rise+b_rise, '-r')
    # plt.plot(v_const, a_const*v_const+b_const, '-r')
    # plt.plot(v, a_rise*v+b_rise, '--r')
    # plt.plot(v, a_const*v+b_const, '--r')
    # plt.annotate('V$_\mathrm{fd}$: %.1f V (via intersection)\nV$_\mathrm{fd}$: %.1f V (via inflection)\nN$_\mathrm{dop}$: %.2e cm$^{-3}$\n\nT$_\mathrm{avg}$: %.1f $^\circ C$\nH$_\mathrm{avg}$: %.1f ' % \
    #     (v_dep2, v_dep1, conc*1e-6, np.mean(dat[:, 7]), np.mean(dat[:, 8])) + r'$\%$', (x_loc,y_loc), xycoords='figure fraction', color='tab:blue', \
    #     bbox=dict(facecolor='white', edgecolor='tab:blue', boxstyle='round,pad=0.5'))
    # plt.legend(loc='upper right')
    # plt.xlabel('voltage [V]')
    # plt.ylabel('1/C^2 [1/F^2]')
    # plt.grid()
    # plt.tight_layout()
    # plt.savefig(path + '/fig_diode_1c2v_' + fn[:-4] + '.pdf')
    # plt.clf()

    return v_dep1, v_dep2, rho, conc, fn[:-4]



def analyse_mos_data(fn, debug=0):
    path = '/'.join(fn.split('/')[:-1])
    fn = fn.split('/')[-1]
    if debug:
        print(" -- Analysing MOS CV for file:\n\t%s" % fn)

    ## read data
    dat = np.genfromtxt(path + '/' + fn, skip_header=15, encoding='UTF-8')
    t, v, i_tot, c, c2, r =  dat[:, 0], dat[:, 1], dat[:, 2], dat[:, 3], dat[:, 4], dat[:, 5]
    lbl = '_'.join(fn.split('_')[1:9])

    ## analyse
    v_fb1, v_fb2, c_acc, c_inv, t_ox, n_ox, a_acc, b_acc, v_acc, a_dep, b_dep, v_dep, a_inv, b_inv, v_inv, spl_dev, status = analyse_mos(v, c)

    ## plot
    plt.plot(v, c, ls=' ', marker='s', ms=3, label=lbl)
    # plt.plot(v, spl_dev, 'g--', label='1st derivative')
    plt.plot(v_acc, a_acc*v_acc+b_acc, '-r')
    plt.plot(v_dep, a_dep*v_dep+b_dep, '-r')
    plt.plot(v_inv, a_inv*v_inv+b_inv, '-r')
    plt.plot(v, a_dep*v+b_dep, '--r')
    plt.plot(v, a_acc*v+b_acc, '--r')

    x_loc = 0.66
    y_loc = 0.145
    if v_fb1 > 4:
        x_loc = 0.11
    plt.annotate('V$_\mathrm{fb}$: %.1f V (via intersection)\nV$_\mathrm{fb}$: %.1f V (via inflection)\n\nt$_\mathrm{ox}$: %.2f um\nn$_\mathrm{ox}$: %.2e cm$^{-2}$\n\nC$_\mathrm{acc}$: %.1e F\nC$_\mathrm{inv}$: %.1e F\n\nT$_\mathrm{avg}$: %.1f $^\circ C$\nH$_\mathrm{avg}$: %.1f ' % \
        (v_fb2, v_fb1, t_ox, n_ox, np.mean(c_acc), np.mean(c_inv), np.mean(dat[:, 7]), np.mean(dat[:, 8])) + r'$\%$', (x_loc, y_loc), xycoords='figure fraction', color='tab:blue', \
        bbox=dict(facecolor='white', edgecolor='tab:blue', boxstyle='round,pad=0.5'))
    plt.legend(loc='upper left')
    plt.xlabel('voltage [V]')
    plt.ylabel('capacitance [F]')
    plt.ylim([0, 1e-10])
    plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    plt.tight_layout()
    plt.savefig(path + '/fig_mos_' + fn[:-4] + '.pdf')
    plt.clf()

    return v_fb1, v_fb2, t_ox, n_ox, fn[:-4]



def analyse_gcd_data(fn, debug=0):
    path = '/'.join(fn.split('/')[:-1])
    fn = fn.split('/')[-1]
    if debug:
        print(" -- Analysing GCD for file:\n\t%s" % fn)

    ## read data
    dat = np.genfromtxt(path + '/' + fn, skip_header=36, encoding='UTF-8')
    t, v, i_em, i_vsrc, i_hvsrc, v_bias =  dat[:, 0], dat[:, 1], dat[:, 2], dat[:, 3], dat[:, 4], dat[:, 5]
    lbl = '_'.join(fn.split('_')[1:9])

    ## analyse
    i_surf, i_bulk, i_acc, i_dep, i_inv, v_acc, v_dep, v_inv, spl_dev, status = analyse_gcd(v, i_em)

    ## plot
    plt.plot(v, i_em, ls=' ', marker='s', ms=3, label="i_em")
    plt.annotate('I$_\mathrm{surf}$: %.2e A\n\nT$_\mathrm{avg}$: %.1f $^\circ C$\nH$_\mathrm{avg}$: %.1f ' % \
        (i_surf, np.mean(dat[:, 7]), np.mean(dat[:, 8])) + r'$\%$', (0.70,0.2), xycoords='figure fraction', color='tab:blue', \
        bbox=dict(facecolor='white', edgecolor='tab:blue', boxstyle='round,pad=0.5'))
    plt.legend(loc='upper left')
    plt.xlabel('voltage [V]')
    plt.ylabel('current [A]')
    # plt.ylim([-9e-7, +1e-7])
    plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    plt.tight_layout()
    plt.savefig(path + '/fig_gcd_' + fn[:-4] + '.pdf')
    plt.clf()

    return i_surf, i_bulk, fn[:-4]



def analyse_fet_data(fn, debug=0):
    path = '/'.join(fn.split('/')[:-1])
    fn = fn.split('/')[-1]
    if debug:
        print(" -- Analysing FET for file:\n\t%s" % fn)

    ## read data
    dat = np.genfromtxt(path + '/' + fn, skip_header=34, encoding='UTF-8')
    t, v, i_em, i_vsrc, i_hvsrc, v_bias =  dat[:, 0], dat[:, 1], dat[:, 2], dat[:, 3], dat[:, 4], dat[:, 5]
    lbl = '_'.join(fn.split('_')[1:9])

    ## analyse
    v_th, a, b, spl_dev, status = analyse_fet(v, i_em)

    ## plot
    fig, ax1 = plt.subplots()
    lns1 = ax1.plot(v, i_em, ls=' ', marker='s', ms=3, label="transfer characteristics")
    ax1.set_xlabel(r'V$_\mathrm{GS}$ [V]')
    ax1.set_ylabel(r'I$_\mathrm{D}$ [A]')
    ax1.set_ylim([-1e-6, +5e-6])

    ax2 = ax1.twinx()
    lns2 = ax2.plot(v, spl_dev, ls=' ', marker='s', ms=3, color='tab:orange', label="transconductance")
    ax2.tick_params(axis='y', labelcolor='tab:orange')
    ax2.set_ylabel(r'g$_\mathrm{m}$ [S]', color='tab:orange')

    lns3 = ax1.plot(v, a*v+b, '--r', label="tangent")

    lns = lns1+lns2+lns3
    labs = [l.get_label() for l in lns]
    plt.legend(lns, labs, loc='upper left')

    plt.annotate('V$_\mathrm{th}$: %.2f + 0.05 V (via tangent)\n\nT$_\mathrm{avg}$: %.1f $^\circ C$\nH$_\mathrm{avg}$: %.1f ' % \
        (v_th, np.mean(dat[:, 7]), np.mean(dat[:, 8])) + r'$\%$', (0.13, 0.6), xycoords='figure fraction', color='tab:blue', \
        bbox=dict(facecolor='white', edgecolor='tab:blue', boxstyle='round,pad=0.5'))
    plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    plt.tight_layout()
    plt.savefig(path + '/fig_fet_' + fn[:-4] + '.pdf')
    plt.clf()

    return v_th, fn[:-4]



def analyse_van_der_pauw_data(fn, debug=0):
    path = '/'.join(fn.split('/')[:-1])
    fn = fn.split('/')[-1]
    if debug:
        print(" -- Analysing Van der Pauw for file:\n\t%s" % fn)

    ## read data
    dat = np.genfromtxt(path + '/' + fn, skip_header=13)
    t, i, v =  dat[:, 0], dat[:, 1], dat[:, 2]
    lbl = '_'.join(fn.split('_')[1:9])

    ## analyse
    r_sheet, a, b, x_fit, spl_dev, status = analyse_van_der_pauw(i, v)

    ## plot
    plt.plot(x_fit, a*x_fit+b, '-r')
    plt.plot(i, v, ls=' ', marker='s', ms=3, label=lbl)
    plt.annotate('R$_\mathrm{sheet}$: %.2e $\Omega$/sq\n\nT$_\mathrm{avg}$: %.1f $^\circ C$\nH$_\mathrm{avg}$: %.1f ' % \
        (r_sheet, np.mean(dat[:, 4]), np.mean(dat[:, 5])) + r'$\%$', (0.20, 0.6), xycoords='figure fraction', color='tab:blue', \
        bbox=dict(facecolor='white', edgecolor='tab:blue', boxstyle='round,pad=0.5'))
    plt.legend(loc='upper left')
    plt.xlabel('current [A]')
    plt.ylabel('voltage [V]')
    # plt.ylim([-9e-7, +1e-7])
    plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    plt.tight_layout()
    plt.savefig(path + '/fig_van_der_pauw_' + fn[:-4] + '.pdf')
    plt.clf()

    return r_sheet, fn[:-4]



def analyse_cross_data(fn, debug=0):
    path = '/'.join(fn.split('/')[:-1])
    fn = fn.split('/')[-1]
    if debug:
        print(" -- Analysing Cross Structure for file:\n\t%s" % fn)

    ## read data
    dat = np.genfromtxt(path + '/' + fn, skip_header=13)
    t, i, v =  dat[:, 0], dat[:, 1], dat[:, 2]
    lbl = '_'.join(fn.split('_')[1:9])

    ## analyse
    r_sheet, a, b, x_fit, spl_dev, status = analyse_cross(i, v)

    ## plot
    plt.plot(x_fit, a*x_fit+b, '-r')
    plt.plot(i, v, ls=' ', marker='s', ms=3, label=lbl)
    plt.annotate('R$_\mathrm{sheet}$: %.2e $\Omega$/sq\n\nT$_\mathrm{avg}$: %.1f $^\circ C$\nH$_\mathrm{avg}$: %.1f ' % \
        (r_sheet, np.mean(dat[:, 4]), np.mean(dat[:, 5])) + r'$\%$', (0.20, 0.6), xycoords='figure fraction', color='tab:blue', \
        bbox=dict(facecolor='white', edgecolor='tab:blue', boxstyle='round,pad=0.5'))
    plt.legend(loc='upper left')
    plt.xlabel('current [A]')
    plt.ylabel('voltage [V]')
    # plt.ylim([-9e-7, +1e-7])
    plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    plt.tight_layout()
    plt.savefig(path + '/fig_cross_' + fn[:-4] + '.pdf')
    plt.clf()

    return r_sheet, fn[:-4]



def analyse_linewidth_data(fn, r_sheet=-1, debug=0):
    path = '/'.join(fn.split('/')[:-1])
    fn = fn.split('/')[-1]
    if debug:
        print(" -- Analysing linewidth for file:\n\t%s" % fn)

    ## read data
    dat = np.genfromtxt(path + '/' + fn, skip_header=13)
    t, i, v =  dat[:, 0], dat[:, 1], dat[:, 2]
    lbl = '_'.join(fn.split('_')[1:9])

    ## analyse
    t_line, a, b, x_fit, spl_dev, status = analyse_linewidth(i, v, r_sheet, cut_param=0.01, debug=0)

    ## plot
    plt.plot(i, v, ls=' ', ms=3, marker='s', label=lbl)
    plt.plot(x_fit, a*x_fit+b, '-r')
    plt.annotate('T$_\mathrm{line}$: %.2e um\n\nT$_\mathrm{avg}$: %.1f $^\circ C$\nH$_\mathrm{avg}$: %.1f ' % \
        (t_line, np.mean(dat[:, 4]), np.mean(dat[:, 5])) + r'$\%$', (0.20, 0.6), xycoords='figure fraction', color='tab:blue', \
        bbox=dict(facecolor='white', edgecolor='tab:blue', boxstyle='round,pad=0.5'))
    plt.legend(loc='upper left')
    plt.xlabel('current [A]')
    plt.ylabel('voltage [V]')
    plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    plt.tight_layout()
    plt.savefig(path + '/fig_linewidth_' + fn[:-4] + '.pdf')
    plt.clf()

    return t_line, fn[:-4]



def analyse_cbkr_data(fn, r_sheet=-1, debug=0):
    path = '/'.join(fn.split('/')[:-1])
    fn = fn.split('/')[-1]
    if debug:
        print(" -- Analysing CBKR for file:\n\t%s" % fn)

    ## read data
    dat = np.genfromtxt(path + '/' + fn, skip_header=13)
    t, i, v =  dat[:, 0], dat[:, 1], dat[:, 2]
    lbl = '_'.join(fn.split('_')[1:9])

    ## analyse
    r_contact, a, b, x_fit, spl_dev, status = analyse_cbkr(i, v, r_sheet, cut_param=0.01, debug=0)

    ## plot
    plt.plot(i, v, ls=' ', ms=3, marker='s', label=lbl)
    plt.plot(x_fit, a*x_fit+b, '-r')
    plt.annotate('R$_\mathrm{contact}$: %.2e $\Omega$\n\nT$_\mathrm{avg}$: %.1f $^\circ C$\nH$_\mathrm{avg}$: %.1f ' % \
        (r_contact, np.mean(dat[:, 4]), np.mean(dat[:, 5])) + r'$\%$', (0.20, 0.6), xycoords='figure fraction', color='tab:blue', \
        bbox=dict(facecolor='white', edgecolor='tab:blue', boxstyle='round,pad=0.5'))
    plt.legend(loc='upper left')
    plt.xlabel('current [A]')
    plt.ylabel('voltage [V]')
    plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    plt.tight_layout()
    plt.savefig(path + '/fig_cbkr_' + fn[:-4] + '.pdf')
    plt.clf()

    return r_contact, fn[:-4]



def analyse_contact_data(fn, debug=0):
    path = '/'.join(fn.split('/')[:-1])
    fn = fn.split('/')[-1]
    if debug:
        print(" -- Analysing contact for file:\n\t%s" % fn)

    ## read data
    dat = np.genfromtxt(path + '/' + fn, skip_header=13)
    t, i, v =  dat[:, 0], dat[:, 1], dat[:, 2]
    lbl = '_'.join(fn.split('_')[1:9])

    ## analyse
    r_contact, a, b, x_fit, spl_dev, status = analyse_contact(i, v, cut_param=0.01, debug=0)

    ## plot
    plt.plot(i, v, ls=' ', ms=3, marker='s', label=lbl)
    plt.plot(x_fit, a*x_fit+b, '-r')
    plt.annotate('R$_\mathrm{contact}$: %.2e $\Omega$\n\nT$_\mathrm{avg}$: %.1f $^\circ C$\nH$_\mathrm{avg}$: %.1f ' % \
        (r_contact, np.mean(dat[:, 4]), np.mean(dat[:, 5])) + r'$\%$', (0.20, 0.6), xycoords='figure fraction', color='tab:blue', \
        bbox=dict(facecolor='white', edgecolor='tab:blue', boxstyle='round,pad=0.5'))
    plt.legend(loc='upper left')
    plt.xlabel('current [A]')
    plt.ylabel('voltage [V]')
    plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    plt.tight_layout()
    plt.savefig(path + '/fig_contact_chain_' + fn[:-4] + '.pdf')
    plt.clf()

    return r_contact, fn[:-4]



def analyse_meander_data(fn, debug=0):
    path = '/'.join(fn.split('/')[:-1])
    fn = fn.split('/')[-1]
    if debug:
        print(" -- Analysing Meander for file:\n\t%s" % fn)

    ## read data
    dat = np.genfromtxt(path + '/' + fn, skip_header=27)
    t, v, i_tot, i =  dat[0], dat[1], dat[2], dat[3]
    lbl = '_'.join(fn.split('_')[1:9])

    ## analyse
    rho_sq, status = analyse_meander(i, v, debug=0)

    return rho_sq, fn[:-4]



def analyse_breakdown_data(fn, debug=0):
    path = '/'.join(fn.split('/')[:-1])
    fn = fn.split('/')[-1]
    if debug:
        print(" -- Analysing Si02 breakdown for file:\n\t%s" % fn)

    ## read data
    dat = np.genfromtxt(path + '/' + fn, skip_header=27)
    t, v, i_tot, i =  dat[:, 0], dat[:, 1], dat[:, 2], dat[:, 3]
    lbl = '_'.join(fn.split('_')[1:9])

    ## analyse
    v_bd, status = analyse_breakdown(v, i)

    ## plot
    plt.plot(v, i, ls=' ', marker='s', ms=3, label=lbl)
    plt.annotate('V$_\mathrm{bd}$: %d V\n\nT$_\mathrm{avg}$: %.1f $^\circ C$\nH$_\mathrm{avg}$: %.1f ' % \
        (v_bd, np.mean(dat[:, 5]), np.mean(dat[:, 6])) + r'$\%$', (0.20,0.5), xycoords='figure fraction', color='tab:blue', \
        bbox=dict(facecolor='white', edgecolor='tab:blue', boxstyle='round,pad=0.5'))
    plt.legend(loc='upper left')
    plt.xlabel('voltage [V]')
    plt.ylabel('current [A]')
    plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    plt.tight_layout()
    plt.savefig(path + '/fig_breakdown_' + fn[:-4] + '.pdf')
    plt.clf()

    return v_bd, fn[:-4]



def analyse_capacitor_data(fn, debug=0):
    path = '/'.join(fn.split('/')[:-1])
    fn = fn.split('/')[-1]
    if debug:
        print(" -- Analysing capacitor for file:\n\t%s" % fn)

    ## read data
    dat = np.genfromtxt(path + '/' + fn, skip_header=15)
    t, v, i_hvsrc, c, c2, r =  dat[:, 0], dat[:, 1], dat[:, 2], dat[:, 3], dat[:, 4], dat[:, 5]
    lbl = '_'.join(fn.split('_')[1:9])

    ## analyse
    c_mean, c_median, status = analyse_capacitor(v, c)

    ## plot
    plt.plot(v, c, ls=' ', marker='s', ms=3, label=lbl)
    plt.annotate('C$_\mathrm{mean}$: %.2e F\nC$_\mathrm{median}$: %.2e F\n\nT$_\mathrm{avg}$: %.1f $^\circ C$\nH$_\mathrm{avg}$: %.1f ' % \
        (c_mean, c_median, np.mean(dat[:, 7]), np.mean(dat[:, 8])) + r'$\%$', (0.20,0.5), xycoords='figure fraction', color='tab:blue', \
        bbox=dict(facecolor='white', edgecolor='tab:blue', boxstyle='round,pad=0.5'))
    plt.legend(loc='upper left')
    plt.xlabel('voltage [V]')
    plt.ylabel('capacitance [F]')
    plt.grid(alpha=0.5, linestyle='--', linewidth=1)
    plt.tight_layout()
    plt.savefig(path + '/fig_capacitor_' + fn[:-4] + '.pdf')
    plt.clf()

    return c_mean, c_median, fn[:-4]



## Macro Functions
## ------------------------------------

def analyse_folder(path, debug=0):

    # flags
    fPrint = 0

    # init return dicts
    diode_iv_dict = {}
    diode_cv_dict = {}
    mos_dict = {}
    gcd_dict = {}
    fet_dict = {}
    van_der_pauw_dict = {}
    cross_dict = {}
    linewidth_dict = {}
    cbkr_dict = {}
    contact_dict = {}
    meander_dict = {}
    breakdown_dict = {}
    capacitor_dict = {}

    # fetch files
    file_list = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and f[-4:] == '.txt']

    r_sheet_ntop = r_sheet_pstop = r_sheet_pbulk = -1

    # loop over and look for fitting analysis
    for f in file_list:
        vals = f.split('_')

        if 1:
            if 'iv' in [v.lower() for v in f.split('_')]:
                v_max, i_max, i_800, i_600, lbl = analyse_iv_data(path + '/' + f)
                diode_iv_dict[lbl] = [v_max, i_max, i_800, i_600]
                if fPrint:
                    print('\nDiode IV: %s' % lbl)
                    print('\t%d V\t%.2e A\t%.2e A\t%.2e A' % (v_max, i_max, i_800, i_600))

        if 1:
            if 'cv' in [v.lower() for v in f.split('_')]:
                v_dep1, v_dep2, rho, conc, lbl  = analyse_cv_data(path + '/' + f)
                diode_cv_dict[lbl] = [v_dep1, v_dep2, rho, conc]
                if fPrint:
                    print('\nDiode CV: %s' % lbl)
                    print('\t%d V\t%d V\t%.2e kOhm\t%.2e cm^-3' % (v_dep1, v_dep2, rho*1e-3, conc*1e-6))

        if 1:
            if 'mos' in [v.lower() for v in f.split('_')]:
                v_fb1, v_fb2, t_ox, n_ox, lbl = analyse_mos_data(path + '/' + f)
                mos_dict[lbl] = [v_fb1, v_fb2, t_ox, n_ox]
                if fPrint:
                    print('\nMOS: %s' % lbl)
                    print('\t%.2f V\t%.2f V\t%.2f um\t%.2e cm^-3' % (v_fb1, v_fb2, t_ox, n_ox))

        if 1:
            if 'gcd' in [v.lower() for v in f.split('_')] or 'gcd05' in [v.lower() for v in f.split('_')]:
                i_surf, i_bulk, lbl = analyse_gcd_data(path + '/' + f)
                gcd_dict[lbl] = [i_surf, i_bulk]
                if fPrint:
                    print('\nGCD: %s' % lbl)
                    print('\t%.2e A\t%.2e A' % (i_surf, i_bulk))

        if 1:
            if 'fet' in [v.lower() for v in f.split('_')]:
                v_th, lbl = analyse_fet_data(path + '/' + f)
                fet_dict[lbl] = [v_th]
                if fPrint:
                    print('\nFET: %s' % lbl)
                    # FIXME: print('\t%.2f V\t%.2f V' % (v_th1, v_th2))
                    print('\t%.2f V' % v_th)

        if 1:
            if 'van-der-pauw' in [v.lower() for v in f.split('_')] or 'cross' in [v.lower() for v in f.split('_')]:
                r_sheet, lbl = analyse_van_der_pauw_data(path + '/' + f)
                van_der_pauw_dict[lbl] = [r_sheet]
                if fPrint:
                    print('\nVan-der-Pauw: %s' % lbl)
                    print('\t%.2e Ohm/sq' % (r_sheet))

        if 1:
            if 'linewidth' in [v.lower() for v in f.split('_')]:
                t_line, lbl = analyse_linewidth_data(path + '/' + f, r_sheet_ntop)
                linewidth_dict[lbl] = [t_line]
                if fPrint:
                    print('\nLinewidth: %s' % lbl)
                    print('\t%.2f um' % (t_line))

        if 1:
            if 'cbkr' in [v.lower() for v in f.split('_')]:
                r_contact, lbl = analyse_cbkr_data(path + '/' + f, r_sheet_ntop)
                cbkr_dict[lbl] = [r_contact]
                if fPrint:
                    print('\nCBKR: %s' % lbl)
                    print('\t%.2e Ohm' % (r_contact))
        if 1:
            if 'contact' in [v.lower() for v in f.split('_')]:
                 r_contact, lbl = analyse_contact_data(path + '/' + f)
                 contact_dict[lbl] = [r_contact]
                 if fPrint:
                     print('\nContact: %s' % lbl)
                     print('\t%.2e Ohm' % (r_contact))

        if 1:
            if 'meander' in [v.lower() for v in f.split('_')]:
                rho_sq, lbl = analyse_meander_data(path + '/' + f)
                meander_dict[lbl] = [rho_sq]
                if fPrint:
                    print('\nMeander: %s' % lbl)
                    print('\t%.2e specific resistance per sq' % (rho_sq))

        if 1:
            if 'breakdown' in [v.lower() for v in f.split('_')]:
                v_bd, lbl = analyse_breakdown_data(path + '/' + f)
                breakdown_dict[lbl] = [v_bd]
                if fPrint:
                    print('\nBreakdown: %s' % lbl)
                    print('\t%d V' % (v_bd))

        if 1:
            if 'capacitor' in [v.lower() for v in f.split('_')] and not 'mos' in [v.lower() for v in f.split('_')]:
                c_mean, c_median, lbl = analyse_capacitor_data(path + '/' + f)
                capacitor_dict[lbl] = [c_mean, c_median]
                if fPrint:
                    print('\nCapacitor: %s' % lbl)
                    print('\t%.2e F\t%.2e F' % (c_mean, c_median))

    return [diode_iv_dict,diode_cv_dict,mos_dict,gcd_dict,fet_dict,van_der_pauw_dict,cross_dict,linewidth_dict,cbkr_dict,contact_dict,meander_dict,breakdown_dict,capacitor_dict]



## Main Executable
## ------------------------------------

def main():
    usage = "usage: ./pqc_analyse.py -f [path_to_folder]"

    parser = OptionParser(usage=usage, version="prog 0.1")
    parser.add_option("-f", "--folder", action="store", dest="input_path", type="string", \
        help="path to input folder")

    parser.add_option("--ex", "--examples", action="store_true", dest="fExamples",  help="print examples")

    (options, args) = parser.parse_args()

    if options.fExamples:
        print("\nSome example commands for running this script\n")
        print("-- ./pqc_analyse.py")
        print("-- ./pqc_analyse.py -f /data/hgc/pqc/HPK_8in_LD_2019_TS_1001_UL_300um_5V")

    elif options.input_path:
        path = options.input_path
        analyse_folder(path)

    else:
        paths = [
            '/Users/Home/Documents/Works/pqc/pqc_analysis/data/HPK_8in_LD_2019_TS_1001_UL_300um_5V',
            '/Users/Home/Documents/Works/pqc/pqc_analysis/data/HPK_8in_LD_2019_TS_1002_UL_300um_2V'
        ]

        res_dict = {}
        for path in paths:
            key = path.split('/')[-1]
            print(key)
            res_dict[key] = analyse_folder(path)


if __name__ == "__main__":
    main()
