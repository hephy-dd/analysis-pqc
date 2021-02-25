#!/usr/bin/env python3

import pqc_analysis_json as pqc
import matplotlib.pyplot as plt
import matplotlib 
from matplotlib import gridspec
import numpy as np
import os
import re
import glob
from collections import namedtuple
import datetime
import dateutil.parser as parser
from matplotlib import colors
from itertools import islice

def params(names):
    """Function decorator returning namedtuples."""
    def params(f):
        def params(*args, **kwargs):
           return namedtuple(f.__name__, names)(*f(*args, **kwargs))
        return params
    return params



def make_chunks(data, size):
    it = iter(data)

    for i in range(0, len(data), size):
        yield [k for k in islice(it, size)]

def num2str(num, basenum=None):
    if basenum is None:
        basenum = num
    if basenum < 10:
        ret = "{:4.2f}".format(num)
    else:
        ret = "{:4.1f}".format(num)
    return ret

class PQC_value:
    def __init__(self, name='na', nicename='na', expectedValue=0., unit='', showmultiplier=1e0, stray=0.5, values=None):
        self.values = []
        self.name = name
        self.nicename = nicename
        self.unit = unit
        self.showmultiplier = showmultiplier
        self.expectedValue = expectedValue
        self.minAllowed = expectedValue * (1-stray)
        self.maxAllowed = expectedValue * (1+stray)
        self.stray = stray
        self.expectedValue = expectedValue

        if values is not None:
            self.values = values
            
    def __str__(self):
        return self.name+str(np.array(self.values)*self.showmultiplier)+self.unit
        
    def append(self, val):
        self.values.append(val)
            
    def rearrange(self, indices):
        self.values = [ self.values[indices[i]] for i in range(0, len(indices)) ]

    def getValue(self, index):
        # with multiplier to suit the unit
        return self.values[index]*self.showmultiplier
        
    def split(self, itemsperclice):
        ret = [PQC_value(self.name, self.nicename, self.expectedValue, self.unit, self.showmultiplier, self.stray, values=i) for i in make_chunks(self.values, itemsperclice)]
        return ret
        
    @params('values, nTot, nNan, nTooHigh, nTooLow, totAvg, totStd, totMed, selAvg, selStd, selMed')
    def getStats(self, minAllowed=None, maxAllowed=None):
        if minAllowed is None:
            minAllowed = self.minAllowed
        if maxAllowed is None:
            maxAllowed = self.maxAllowed
        
        nTot = len(self.values)

        selector = np.isfinite(np.array(self.values))
        
        if np.sum(selector) < 2:
            return np.array([0]), 1, 1, 0, 0, 0, 0, 0, 0, 0, 0
            
        values = np.array(self.values)[selector]*self.showmultiplier   # filter out nans
        
        if nTot < 2:
            return values, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0
        
        totMed = np.median(values)
        totAvg = np.mean(values)
        totStd = np.std(values)
        
        nNan = nTot - len(values)
        values = values[values < maxAllowed]
        nTooHigh = nTot - len(values) - nNan
        values = values[values > minAllowed]
        nTooLow = nTot - len(values) - nNan - nTooHigh
        
        selMed = np.median(values)
        selAvg = np.mean(values)
        selStd = np.std(values)
        
        return values, nTot, nNan, nTooHigh, nTooLow, totAvg, totStd, totMed, selAvg, selStd, selMed
     
    @classmethod
    def merge(new, parents, name='na', nicename='na'):
        values = np.concatenate( [t.values for t in parents])
        return new(name, nicename, parents[0].expectedValue, parents[0].unit, parents[0].showmultiplier, values=values, stray=parents[0].stray)
    

    
    # get a colorized value string for ise in latex, if the index is higher we get summary elemets
    def valueToLatex(self, index):
        if index < len(self.values):
            value = self.values[index]*self.showmultiplier
            vstr = num2str(value, self.expectedValue)
                
            if np.isnan(value):
                return "\\nanval NaN"
            elif value > self.maxAllowed:
                return "\\highval "+vstr
            elif value < self.minAllowed:
                return "\\lowval "+vstr
            return "\\okval "+vstr
        else:
            stats = self.getStats()
            sel={
                0: "\\okval "+num2str(stats.selMed, self.expectedValue),
                1: "\\okval "+num2str(stats.selAvg, self.expectedValue),
                2: "\\okval "+num2str(stats.selStd, self.expectedValue),
                3: "\\okval {}/{}".format(len(stats.values), stats.nTot),
                4: "\\okval {:2.0f}\\%".format(len(stats.values)/stats.nTot*100),
             }
        return sel.get(index-len(self.values), "\\nanval err")
        
    def summaryDesciption(self, index):
        v = ["\\hline \nMedian", "Average", "Std dev", "\\hline \nOK/Tot", "OK (rel)"]
        return v[index-len(self.values)]
    
    def headerToLatex():
        return "\\def\\nanval{\\cellcolor[HTML]{aa0000}}\n" \
               "\\def\\highval{\\cellcolor[HTML]{ff9900}}\n" \
               "\\def\\lowval{\\cellcolor[HTML]{ffff00}}\n" \
               "\\def\\okval{\\cellcolor[HTML]{ffffff}}\n\n" \
               

class PQC_resultset:
    def __init__(self, batchname, dataseries=None):
        self.batch = batchname
        self.labels = []
        self.flutes = []
        self.timestamps = []
        
        self.dataseries = {'xtimestamps': self.timestamps,
                           'xlabels':     self.labels,
                           'xflutes':     self.flutes,}

        self.dataseries['vdp_poly_f'] = PQC_value("vdp_poly", "Polysilicon VdP", 2.4, "kOhm/sq", 1e-3, stray=0.2)
        self.dataseries['vdp_poly_r'] = PQC_value("vdp_poly_rev", "Polysilicon VdP reverse", 2.4, "kOhm/sq", 1e-3, stray=0.2)

        self.dataseries['vdp_n_f'] = PQC_value("vdp_N", "N+ VdP", 35., "Ohm/sq", stray=0.2)
        self.dataseries['vdp_n_r'] = PQC_value("vdp_N_rev", "N+ VdP reverse", 35., "Ohm/sq", stray=0.2)

        self.dataseries['vdp_pstop_f'] = PQC_value("vdp_pstop", "P-stop VdP", 19., "kOhm/sq", 1e-3, stray=0.2)
        self.dataseries['vdp_pstop_r'] = PQC_value("vdp_pstop_rev", "P-stop VdP reverse", 19., "kOhm/sq", 1e-3, stray=0.2)

        self.dataseries['t_line_n'] = PQC_value("t_line_n", "Linewidth N+", 35., "um")
        self.dataseries['t_line_pstop2'] = PQC_value("t_line_pstop2", "Linewidth P-stop 2 Wire", 38., "um")
        self.dataseries['t_line_pstop4'] = PQC_value("t_line_pstop4", "Linewidth P-stop 4 Wire", 55., "um")

        self.dataseries['r_contact_n'] = PQC_value("r_contact_n", "Rcontact N+", 27., "Ohm")
        self.dataseries['r_contact_poly'] = PQC_value("r_contact_poly", "Rcontact polysilicon", 100., "kOhm", 1e-3)


        self.dataseries['v_th'] = PQC_value("fet", "FET Vth", 4., "V", stray=0.25)

        self.dataseries['vdp_metclo_f'] = PQC_value("vdp_met_clover", "Metal Cloverleaf VdP", 25., "mOhm/sq", 1e3)
        self.dataseries['vdp_metclo_r'] = PQC_value("vdp_met_clover_rev", "Metal Cloverleaf VdP reverse", 25., "mOhm/sq", 1e3)

        self.dataseries['vdp_p_cross_bridge_f'] = PQC_value("vdp_cross_bridge", "Cross Bridge VdP", 1.5, "kOhm/sq", 1e-3)
        self.dataseries['vdp_p_cross_bridge_r'] = PQC_value("vdp_cross_bridge_rev", "Cross Bridge VdP reverse", 1.5, "kOhm/sq", 1e-3)

        self.dataseries['t_line_p_cross_bridge'] = PQC_value("t_line_cb", "Linewidth cross bridge P", 35., "um")

        self.dataseries['v_bd'] = PQC_value("v_bd", "Breakdown Voltage", 215., "V")


        self.dataseries['i600'] = PQC_value("i600", "I @ 600V", 100., "uA", 1e6, stray=1.)
        self.dataseries['v_fd'] = PQC_value("v_fd", "Full depletion Voltage", 260., "V", stray=0.33)
        self.dataseries['rho'] = PQC_value("rho", "rho", 1.3, "kOhm cm", 0.1)
        self.dataseries['conc'] = PQC_value("d_conc", "Doping Concentration", 3.5, "1E12cm^-3", 1e-18)


        self.dataseries['v_fb2'] = PQC_value("v_fb", "Flatband voltage", 2.5, "V", stray=0.33)

        self.dataseries['t_ox'] = PQC_value("t_ox", "Oxide thickness", 0.67, "um", stray=0.33)
        self.dataseries['n_ox'] = PQC_value("n_ox", "Oxide concentration", 10.5, "1E10cm^-3", 1e-10)
        self.dataseries['c_acc_m'] = PQC_value("c_acc", "Accumulation capacitance", 85., "pF", 1e12, stray=0.2)

        self.dataseries['i_surf'] = PQC_value("i_surf", "Surface current", 8., "pA", -1e12, stray=1)
        self.dataseries['i_surf05'] = PQC_value("i_surf05", "Surface current 05", 11., "pA", -1e12, stray=1)
        self.dataseries['i_bulk05'] = PQC_value("i_bulk05", "Bulk current 05", 0.7, "pA", -1e12, stray=1)
        
        self.dataseries['nvdp_poly_f'] = PQC_value("nvdp_poly", "PolySi Swapped VdP", 2.4, "kOhm/sq", -1e-3, stray=0.2)
        self.dataseries['nvdp_poly_r'] = PQC_value("nvdp_poly_rev", "PolySi Swapped VdP reverse", 2.4, "kOhm/sq", -1e-3, stray=0.2)

        self.dataseries['nvdp_n_f'] = PQC_value("nvdp_N", "N+ Swapped VdP", 35., "Ohm/sq", -1., stray=0.2)
        self.dataseries['nvdp_n_r'] = PQC_value("nvdp_N_rev", "N+ Swapped VdP reverse", 35., "Ohm/sq", -1., stray=0.2)

        self.dataseries['nvdp_pstop_f'] = PQC_value("nvdp_pstop", "P-stop Swapped VdP", 19., "kOhm/sq", -1e-3, stray=0.2)
        self.dataseries['nvdp_pstop_r'] = PQC_value("nvdp_pstop_rev", "P-stop Swapped VdP rev", 19., "kOhm/sq", -1e-3, stray=0.2)
        
        self.dataseries['vdp_bulk_f'] = PQC_value("vdp_bulk", "Bulk VdP Cross", 150., "kOhm/sq", 1e-3, stray=0.8)
        self.dataseries['vdp_bulk_r'] = PQC_value("vdp_bulk_rev", "Bulk VdP Cross rev", 150., "kOhm/sq", 1e-3, stray=0.8)
        
        self.dataseries['meander_poly'] = PQC_value("meander_poly", "Polisilicon Resistor", 1.7, "MOhm", 1e-6, stray=0.5)
        self.dataseries['meander_metal'] = PQC_value("meander_metal", "Metal Meander", 260., "Ohm", 1., stray=0.5)
        
        self.dataseries['contact_poly'] = PQC_value("contact_poly", "Contact Chain PolySi", 30., "MOhm", 1e-6, stray=0.5)
        self.dataseries['contact_p'] = PQC_value("contact_p", "Contact Chain P", 85., "kOhm", 1e-3, stray=0.5)
        self.dataseries['contact_n'] = PQC_value("contact_n", "Contact Chain N", 85., "kOhm", 1e-3, stray=0.5)
        
        
        if dataseries is not None:
            self.dataseries = dataseries

    def vdp_poly_tot(self):
        return PQC_value.merge([self.dataseries['vdp_poly_f'], self.dataseries['vdp_poly_r']], "vdp_poly_tot", "PolySi VdP both")
    def vdp_n_tot(self):
        return PQC_value.merge([self.dataseries['vdp_n_f'], self.dataseries['vdp_n_f']], "vdp_N_tot", "N+ VdP both")     
    def vdp_pstop_tot(self):
        return PQC_value.merge([self.dataseries['vdp_pstop_f'], self.dataseries['vdp_pstop_f']], "vdp_pstop_tot", "P-stop VdP both")
        
    def nvdp_poly_tot(self):
        return PQC_value.merge([self.dataseries['nvdp_poly_f'], self.dataseries['nvdp_poly_r']], "nvdp_poly_tot", "PolySi Swapped VdP")
    def nvdp_n_tot(self):
        return PQC_value.merge([self.dataseries['nvdp_n_f'], self.dataseries['nvdp_n_r']], "nvdp_N_tot", "N+ Swapped VdP both")
    def nvdp_pstop_tot(self):
        return PQC_value.merge([self.dataseries['nvdp_pstop_f'], self.dataseries['nvdp_pstop_r']], "nvdp_pstop_tot", "P-stop Swapped VdP")
    
    def sortByTime(self):
        order = np.argsort(self.timestamps)
        #print(str(order))
        
        for key in self.dataseries:
            if type(self.dataseries[key]) is PQC_value:
                self.dataseries[key].rearrange(order)
            else:
                self.dataseries[key][:] = [self.dataseries[key][i] for i in order]  # we want to keep the original object so references are preserved
        
    # warning: this creates not full copies, only the dict is available then
    def split(self, itemsperclice):
        #print(str((self.dataseries["vdp_n_f"]).split(itemsperclice)))
        
        #ret = [PQC_resultset( self.value=i) for i in make_chunks(self.value, itemsperclice)]
        ret = [None] * len(self.dataseries["vdp_n_f"].split(itemsperclice))
        
        for i in range(0, len(ret)):
            ret[i] = PQC_resultset(self.batch)
        
        for key in self.dataseries:
            if type(self.dataseries[key]) is PQC_value:
                ch = self.dataseries[key].split(itemsperclice)
            else:
                ch = [i for i in make_chunks(self.dataseries[key], itemsperclice)]

            for i in range(0, len(ch)):
                ret[i].dataseries[key] = ch[i]
        
        for i in ret:
            i.timestamps = i.dataseries['xtimestamps']
            i.flutes = i.dataseries['xflutes']
            i.labels = i.dataseries['xlabels']

        return ret
    
    def analyze(self, dirs):
        print("dirs: "+str(len(dirs)))
        
        
        for i in range(0, len(dirs)):

            for currflute in ["PQCFlutesLeft", "PQCFlutesRight", "PQCFluteLeft", "PQCFluteRight", ]:   # we want both flutes if they are measured, but sometimes it's misspelled
                if len(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "cross"], blacklist=["reverse"])) < 3:
                    continue

                self.labels.append(dirs[i].split("/")[-1])
                self.flutes.append(currflute)

                x = pqc.find_all_files_from_path(dirs[i], "van_der_pauw")
                if i != []:
                    self.timestamps.append(pqc.get_timestamp(x[-1]))
                else:
                    self.timestamps.append(0)
                    
                cross = "cross"
                
                self.dataseries['vdp_poly_f'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "Polysilicon", cross], blacklist=["reverse"]), plotResults=False))
                self.dataseries['vdp_poly_r'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "Polysilicon", "reverse", cross])))

                self.dataseries['vdp_n_f'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "n", cross], blacklist=["reverse"])))
                self.dataseries['vdp_n_r'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "n", "reverse", cross])))
                
                self.dataseries['vdp_pstop_f'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "P_stop", cross], blacklist=["reverse"])))
                self.dataseries['vdp_pstop_r'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "P_stop", "reverse", cross])))
                
                self.dataseries['t_line_n'].append(pqc.analyse_linewidth_data(pqc.find_all_files_from_path(dirs[i], "linewidth", whitelist=[currflute, "n"], single=True), r_sheet=self.dataseries['vdp_n_f'].values[-1], printResults=False, plotResults=False))
                self.dataseries['t_line_pstop2'].append(pqc.analyse_linewidth_data(pqc.find_all_files_from_path(dirs[i], "linewidth", whitelist=[currflute, "P_stop", "2_wire"], single=True), r_sheet=self.dataseries['vdp_pstop_f'].values[-1], printResults=False, plotResults=False))
                self.dataseries['t_line_pstop4'].append(pqc.analyse_linewidth_data(pqc.find_all_files_from_path(dirs[i], "linewidth", whitelist=[currflute, "P_stop", "4_wire"], single=True), r_sheet=self.dataseries['vdp_pstop_f'].values[-1], printResults=False, plotResults=False))

                self.dataseries['r_contact_n'].append(pqc.analyse_cbkr_data(pqc.find_all_files_from_path(dirs[i], "cbkr", whitelist=[currflute, "n"], single=True), r_sheet=self.dataseries['vdp_n_f'].values[-1], printResults=False, plotResults=False))
                self.dataseries['r_contact_poly'].append(pqc.analyse_cbkr_data(pqc.find_all_files_from_path(dirs[i], "cbkr", whitelist=[currflute, "Polysilicon"], single=True), r_sheet=self.dataseries['vdp_poly_f'].values[-1], printResults=False, plotResults=False))

                self.dataseries['v_th'].append(pqc.analyse_fet_data(pqc.find_all_files_from_path(dirs[i], "fet", whitelist=[currflute, ], single=True), printResults=False, plotResults=False))

                self.dataseries['vdp_metclo_f'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "metal", "clover"], blacklist=["reverse"])))
                self.dataseries['vdp_metclo_r'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "metal", "clover", "reverse"], blacklist=[])))

                self.dataseries['vdp_p_cross_bridge_f'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "P", "cross_bridge"], blacklist=["reverse"])))
                self.dataseries['vdp_p_cross_bridge_r'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "P", "cross_bridge", "reverse"])))
                self.dataseries['t_line_p_cross_bridge'].append(pqc.analyse_linewidth_data(pqc.find_all_files_from_path(dirs[i], "linewidth", whitelist=[currflute, "P", "cross_bridge"], single=True), r_sheet=self.dataseries['vdp_p_cross_bridge_f'].values[-1], printResults=False, plotResults=False))

                self.dataseries['v_bd'].append(pqc.analyse_breakdown_data(pqc.find_all_files_from_path(dirs[i], "breakdown", whitelist=[currflute, ], single=True), printResults=False, plotResults=False))

                # we want this for FLute_3 and not Flute_1
                i600, dummy = pqc.analyse_iv_data(pqc.find_all_files_from_path(dirs[i], "iv", whitelist=[currflute, "3"], single=True), printResults=False, plotResults=False)
                self.dataseries['i600'].append(i600)
                
                v_fd, rho, conc = pqc.analyse_cv_data(pqc.find_all_files_from_path(dirs[i], "cv", whitelist=[currflute, "3"], single=True), printResults=False, plotResults=False)
                self.dataseries['v_fd'].append(v_fd)
                self.dataseries['rho'].append(rho)
                self.dataseries['conc'].append(conc)
                
                dummy, v_fb2, t_ox, n_ox, c_acc_m = pqc.analyse_mos_data(pqc.find_all_files_from_path(dirs[i], "mos", whitelist=[currflute, ], single=True), printResults=False, plotResults=False)
                self.dataseries['v_fb2'].append(v_fb2)
                self.dataseries['t_ox'].append(t_ox)
                self.dataseries['n_ox'].append(n_ox)
                self.dataseries['c_acc_m'].append(c_acc_m)
                
                i_surf, dummy = pqc.analyse_gcd_data(pqc.find_all_files_from_path(dirs[i], "gcd", whitelist=[currflute, ], single=True), printResults=False, plotResults=False)  # only i_surf valid
                self.dataseries['i_surf'].append(i_surf)
                
                i_surf05, i_bulk05 = pqc.analyse_gcd_data(pqc.find_all_files_from_path(dirs[i], "gcd05", whitelist=[currflute, ], single=True), printResults=False, plotResults=False)  # for i_bulk
                self.dataseries['i_surf05'].append(i_surf05)
                self.dataseries['i_bulk05'].append(i_bulk05)
                
                self.dataseries['nvdp_poly_f'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "Polysilicon", "ncross"], blacklist=["reverse"]), plotResults=False))
                self.dataseries['nvdp_poly_r'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "Polysilicon", "reverse", "ncross"])))

                self.dataseries['nvdp_n_f'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "n", "ncross"], blacklist=["reverse"])))
                self.dataseries['nvdp_n_r'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "n", "reverse", "ncross"])))

                self.dataseries['nvdp_pstop_f'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "P_stop", "ncross"], blacklist=["reverse"])))
                self.dataseries['nvdp_pstop_r'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[currflute, "P_stop", "reverse", "ncross"])))
            
                self.dataseries['vdp_bulk_f'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], None, whitelist=[currflute, "bulk", "cross"], blacklist=["reverse"])))
                self.dataseries['vdp_bulk_r'].append(pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], None, whitelist=[currflute, "bulk", "reverse", "cross"])))
                
                self.dataseries['meander_metal'].append(pqc.analyse_meander_data(pqc.find_all_files_from_path(dirs[i], "meander", whitelist=[currflute, "metal"], single=True), printResults=False, plotResults=False))
                self.dataseries['meander_poly'].append(pqc.analyse_meander_data(pqc.find_all_files_from_path(dirs[i], "meander", whitelist=[currflute, "polysilicon"], single=True), printResults=False, plotResults=False))
                
                self.dataseries['contact_poly'].append(pqc.analyse_contact_data(pqc.find_all_files_from_path(dirs[i], "contact", whitelist=[currflute, "chain", "polysilicon"], single=True), printResults=False, plotResults=False))
                self.dataseries['contact_p'].append(pqc.analyse_contact_data(pqc.find_all_files_from_path(dirs[i], "contact", whitelist=[currflute, "chain", "P"], single=True), printResults=False, plotResults=False))
                self.dataseries['contact_n'].append(pqc.analyse_contact_data(pqc.find_all_files_from_path(dirs[i], "contact", whitelist=[currflute, "chain", "N"], single=True), printResults=False, plotResults=False))
                
            
    def prettyPrint(self):
        print("# serial                                 \t  vdp_poly/kOhm/sq       vdp_n/Ohm/sq     vdp_pstop/kOhm/sq   lw_n/um    lw_p2/um   lw_p4/um cbkr_poly/kOhm cbkr_n/Ohm")
        for i in range(0, len(self.labels)):
            line = "{} {}  \t".format(self.labels[i], self.flutes[i])
            line += "{:8.2f} {:8.2f}    ".format(self.dataseries['vdp_poly_f'].getValue(i), self.dataseries['vdp_poly_r'].getValue(i))
            line += "{:8.2f} {:8.2f}    ".format(self.dataseries['vdp_n_f'].getValue(i), self.dataseries['vdp_n_r'].getValue(i))
            line += "{:8.2f} {:8.2f}    ".format(self.dataseries['vdp_pstop_f'].getValue(i), self.dataseries['vdp_pstop_r'].getValue(i))
            line += "{:8.2f} {:8.2f} {:8.2f}     ".format(self.dataseries['t_line_n'].getValue(i), self.dataseries['t_line_pstop2'].getValue(i), self.dataseries['t_line_pstop4'].getValue(i))
            line += "{:8.2f} {:8.2f}".format(self.dataseries['r_contact_poly'].getValue(i), self.dataseries['r_contact_n'].getValue(i))

            print(line)
        
        print("")
        print("")
        print("# serial                                 \t fet       vdp_met-clov       vdp_p-cr-br/kOhm/sq     lw_cb  v_bd   i600    V_fd  rho    dconc  contact chains")
        print("#                                        \tv_th/V                                                 um     V      uA       V kOhmcm 1E12cm^-3 Poly/MOhm  N/kOhm  P/kOhm")
        for i in range(0, len(self.labels)):
            line = "{} {}  \t".format(self.labels[i], self.flutes[i])
            line += "{:5.2f}  ".format(self.dataseries['v_th'].getValue(i))
            line += "{:9.2E}  {:9.2E}   ".format(self.dataseries['vdp_metclo_f'].getValue(i), self.dataseries['vdp_metclo_r'].getValue(i))

            line += "{:8.2f} {:8.2f}    ".format(self.dataseries['vdp_p_cross_bridge_f'].getValue(i), self.dataseries['vdp_p_cross_bridge_r'].getValue(i))
            line += "{:8.2f}  ".format(self.dataseries['t_line_p_cross_bridge'].getValue(i))

            line += "{:3.0f}  ".format(self.dataseries['v_bd'].getValue(i))
            line += "{:6.2f}  {:6.2f}  {:4.2f}  {:4.2f}   ".format(self.dataseries['i600'].getValue(i), self.dataseries['v_fd'].getValue(i), self.dataseries['rho'].getValue(i), self.dataseries['conc'].getValue(i))
            
            line += "{:8.2f} {:8.2f}  {:8.2f}  ".format(self.dataseries['contact_poly'].getValue(i), self.dataseries['contact_n'].getValue(i), self.dataseries['contact_p'].getValue(i))
            print(line)
        print("")
        print("")
        print("# serial                                 \t                    mos                        gcd             gcd05            vdp_bulk           poly_r   met_meand")
        print("#                                        \t v_fb/V    c_acc/pF   t_ox/um n_ox/1E10cm^-2 i_surf/pA  i_surf/pA   i_bulk/pA     kOhm/sq            MOhm      Ohm")
        for i in range(0, len(self.labels)):
            line = "{} {}  \t".format(self.labels[i], self.flutes[i])
            line += "{:8.2f}    {:6.2f}    {:7.3f}  {:9.2f}     ".format(self.dataseries['v_fb2'].getValue(i), self.dataseries['c_acc_m'].getValue(i), self.dataseries['t_ox'].getValue(i), self.dataseries['n_ox'].getValue(i))
            line += "{:8.2f}  {:8.2f}  {:8.2f}    ".format(self.dataseries['i_surf'].getValue(i), self.dataseries['i_surf05'].getValue(i), self.dataseries['i_bulk05'].getValue(i))
            line += "{:8.2f} {:8.2f}    ".format(self.dataseries['vdp_bulk_f'].getValue(i), self.dataseries['vdp_bulk_r'].getValue(i))
            line += "{:8.2f} {:8.2f}    ".format(self.dataseries['meander_poly'].getValue(i), self.dataseries['meander_metal'].getValue(i))
            print(line)
        
    def statusbar(self, pqc_value_statistics, axes, single=True, start=-0.5, stop=0.5, label=''):
        relOK = len(pqc_value_statistics.values)/pqc_value_statistics.nTot
        relNan = pqc_value_statistics.nNan/pqc_value_statistics.nTot
        relTooHigh = pqc_value_statistics.nTooHigh/pqc_value_statistics.nTot
        relTooLow = pqc_value_statistics.nTooLow/pqc_value_statistics.nTot
        
        average_delta = (start - stop) / 2
        center = stop + average_delta
        width = stop-start
        if single:
            alpha = 1.
        else:
            alpha = 0.5
        
        p1 = axes.bar(center, (relOK,), width, color="green", alpha=alpha)
        p2 = axes.bar(center, (relNan,), width, bottom=(relOK,), color="red", alpha=alpha)
        p3 = axes.bar(center, (relTooHigh,), width, bottom=(relOK+relNan,), color="orange", alpha=alpha)
        p4 = axes.bar(center, (relTooLow,), width, bottom=(relOK+relNan+relTooHigh,), color="yellow", alpha=alpha)
        
        if single:
            minthreshhold = 0.02
        else:
            minthreshhold = 0.15
            
        if single or label != '':
            axes.text(center, relOK-minthreshhold/2, 'OK', horizontalalignment='center', verticalalignment='top')
            if relNan > minthreshhold:
                axes.text(center, relOK+relNan-minthreshhold/2, 'Failed', horizontalalignment='center', verticalalignment='top')
            if relTooHigh > minthreshhold:
                axes.text(center, relOK+relNan+relTooHigh-minthreshhold/2, 'High', horizontalalignment='center', verticalalignment='top')
            if relTooLow > minthreshhold:
                axes.text(center, relOK+relNan+relTooHigh+relTooLow-minthreshhold/2, 'Low', horizontalalignment='center', verticalalignment='top')
        
        axes.text(center, 0.01, label, horizontalalignment='center', verticalalignment='bottom')
        
        plt.yticks([])
        plt.ylim([0, 1])
        if single:
            plt.xticks([])
            plt.xlim([start, stop])
            
        
    def histogram(self, pqc_values, path, stray=1.4, rangeExtension=None):
        if rangeExtension is not None:
            stats = pqc_values.getStats(minAllowed=0, maxAllowed=pqc_values.maxAllowed*rangeExtension)
        else:
            stats = pqc_values.getStats()
        
        
        if len(stats.values) == 1 and stats.nNan == 1:
            print("warning: skipping plot due to no valid results: "+pqc_values.name)
            return
        
        fig = plt.figure(figsize=(8, 6))
        
        
        gs = gridspec.GridSpec(1, 2, width_ratios=[10, 1]) 
        ax0 = plt.subplot(gs[0])
        
        if rangeExtension is not None:
            plt.title(self.batch + ": " + pqc_values.nicename + ", Ext: {:5.1E}".format(rangeExtension), fontsize=18)
            plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
        else:
            plt.title(self.batch + ": " + pqc_values.nicename + "", fontsize=18)
        
        ax1 = plt.subplot(gs[1])
                
        if rangeExtension is not None:
            #ax0.hist(pqc_values.value, bins=50, facecolor='blueviolet', alpha=1, edgecolor='black', linewidth=1)
            ax0.hist(stats.values, bins=50, range=[0, pqc_values.maxAllowed*rangeExtension], facecolor='blueviolet', alpha=1, edgecolor='black', linewidth=1)
            descNum = "Total: {}\nShown: {}\nFailed: {}\nToo high: {}\nToo low: {}".format(stats.nTot, len(stats.values),  stats.nNan, stats.nTooHigh, stats.nTooLow)
        else:
            ax0.hist(stats.values, bins=20, range=[pqc_values.minAllowed, pqc_values.maxAllowed], facecolor='blue', alpha=1, edgecolor='black', linewidth=1)
            descNum = "Total: {}\nShown: {}\nFailed: {}\nToo high: {}\nToo low: {}".format(stats.nTot, len(stats.values),  stats.nNan, stats.nTooHigh, stats.nTooLow)
        
        ax0.set_xlabel(pqc_values.unit)
        ax0.set_ylabel("occurences")
        
        #descNum = "Total number: {}\nShown: {:2.0f}%, {}\nFailed: {:2.0f}%, {}\nToo high: {:2.0f}%, {}\nToo low: {:2.0f}%, {}".format(stats.nTot, len(stats.values)/stats.nTot*1e2, len(stats.values), stats.nNan/stats.nTot*1e2, stats.nNan, stats.nTooHigh/stats.nTot*1e2, stats.nTooHigh, stats.nTooLow/stats.nTot*1e2, stats.nTooLow)
        
        fig.text(0.83, 0.85, descNum, bbox=dict(facecolor='red', alpha=0.6), horizontalalignment='right', verticalalignment='top')
        
        if abs(stats.selMed) < 9.99:
            descStat = "Total median: {2:8.2f} {3}\n".format(stats.totAvg, stats.totStd, stats.totMed, pqc_values.unit)
            descStat = descStat +"Selected avg: {0:8.2f} {3}\nSel median: {2:8.2f} {3}\nSelected Std: {1:8.2f} {3}".format(stats.selAvg, stats.selStd, stats.selMed, pqc_values.unit)
        elif abs(stats.totAvg) < 1e6:
            descStat = "Total median: {2:8.1f} {3}\n".format(stats.totAvg, stats.totStd, stats.totMed, pqc_values.unit)
            descStat = descStat +"Selected avg: {0:8.1f} {3}\nSel median: {2:8.1f} {3}\nSelected Std: {1:8.1f} {3}".format(stats.selAvg, stats.selStd, stats.selMed, pqc_values.unit)
        else:
            descStat = "Total avg: {0:9.2E} {4}\nTotal median: {1:9.2E} {4}\nSelected avg: {2:9.2E} {4}\nSel median: {3:9.2E} {4}".format(stats.totAvg, stats.totMed, stats.selAvg, stats.selMed, pqc_values.unit)
        fig.text(0.45, 0.85, descStat, bbox=dict(facecolor='yellow', alpha=0.85), horizontalalignment='right', verticalalignment='top')
        
        #for i in [(stats.totAvg, 'purple', 'solid'), (stats.totMed, 'purple', 'dashed'), (stats.selAvg, 'green', 'solid'), (stats.selMed, 'green', 'dashed')]:
        #    if (i[0] < pqc_values.maxAllowed) and (i[0] > pqc_values.minAllowed):
        #        ax0.vlines(x = i[0], ymin = 0, ymax = 3, 
        #            colors = i[1], linestyles=i[2],
        #            label = 'vline_multiple - full height') 
                    

        self.statusbar(stats, ax1)
        
        
        
        fig.tight_layout(h_pad=1.0)
        if rangeExtension is not None:
            fig.savefig(path+"/"+pqc_values.name+"_erhist.png")
        else:
            fig.savefig(path+"/"+pqc_values.name+"_hist.png")
        plt.close()

        
    def createHistograms(self, path):
        matplotlib.rcParams.update({'font.size': 14})

        histogramDir = path+"histograms_"+self.batch
        try:
            os.mkdir(histogramDir)
        except OSError:
            files = glob.glob(histogramDir+"/*")
            for f in files:
                os.remove(f)
        
        for key in self.dataseries:
            if not key.startswith('x'):
                self.histogram(self.dataseries[key], histogramDir)
        
        self.histogram(self.vdp_poly_tot(), histogramDir)
        self.histogram(self.vdp_n_tot(), histogramDir)
        self.histogram(self.vdp_pstop_tot(), histogramDir)   
        
        self.histogram(self.vdp_poly_tot(), histogramDir, rangeExtension=1.5e2)
        self.histogram(self.vdp_n_tot(), histogramDir, rangeExtension=1.5e2)
        self.histogram(self.vdp_pstop_tot(), histogramDir, rangeExtension=1.5e2)        

        self.histogram(self.nvdp_poly_tot(), histogramDir)
        self.histogram(self.nvdp_n_tot(), histogramDir)
        self.histogram(self.nvdp_pstop_tot(), histogramDir)   
        
    def shortLabel(self, i):
        fl = "x"
        try:
            lbl_list = [2,5]
            if "Left" in self.flutes[i]:
                fl = " L"
            elif "Right" in self.flutes[i]:
                fl = " R"
            else:
                fl = " err"
            
            return ' '.join([self.labels[i].split('_')[j] for j in lbl_list]) + fl
        except:
            return self.labels[i] + fl
        
    def shortBatch(self, vpx=True):
        lbl_list = [0]
        s = ' '.join([self.batch.split('_')[i] for i in lbl_list])
        if vpx:
            return s
        else:
            return re.sub(r'VPX','', s)
        
        
    # this could be done e.g. via df.to_latex() in pandas, but coloring the cells is then complicated, so done manually...
    def exportLatex(self, path):
        self.exportLatex1(path)
        self.exportLatex2(path)
        self.exportLatex3(path)
    
    
    def exportLatex1(self, path):
        f = open(path+"/histograms_"+self.batch+"/table1.tex", "w")
        f.write("% automatically created table for batch " + self.batch + "\n\n")
        
        f.write(PQC_value.headerToLatex())
        
        f.write("""\\begin{center}
        \\fontsize{4pt}{5pt}\\selectfont
    \\begin{tabular}{ |l|r|r|r|r|r|r|r|r|r|r|r| } 
        \\hline
        """+self.shortBatch()+""" & \multicolumn{2}{ c|}{PolySi VdP} & \multicolumn{2}{c|}{N+ VdP} & \multicolumn{2}{c|}{P-Stop VdP} & \multicolumn{3}{c|}{line thickness} & \multicolumn{2}{ c|}{Contact resistance} \\\\
        & \multicolumn{2}{c|}{"""+self.dataseries['vdp_poly_f'].unit+"} & \multicolumn{2}{c|}{"+self.dataseries['vdp_n_f'].unit+"} & \multicolumn{2}{c|}{"+self.dataseries['vdp_pstop_f'].unit+"} & \multicolumn{3}{c|}{"+self.dataseries['t_line_n'].unit+"} & "+self.dataseries['r_contact_poly'].unit+" & "+self.dataseries['r_contact_n'].unit+"""\\\\
         & forw & rev & forw & rev & forw & rev & N+ & p-stop2 & p-stop4 & PolySi & N+\\\\
        \\hline\n
        """)
        
        for i in range(0, len(self.labels)+5):
            if i < len(self.labels):
                line = "        \detokenize{"+self.shortLabel(i)+"}"
            else:
                line = self.dataseries['vdp_poly_f'].summaryDesciption(i)
            line = line + " & " + self.dataseries['vdp_poly_f'].valueToLatex(i)
            line = line + " & " + self.dataseries['vdp_poly_r'].valueToLatex(i)
            line = line + " & " + self.dataseries['vdp_n_f'].valueToLatex(i)
            line = line + " & " + self.dataseries['vdp_n_r'].valueToLatex(i)
            line = line + " & " + self.dataseries['vdp_pstop_f'].valueToLatex(i)
            line = line + " & " + self.dataseries['vdp_pstop_r'].valueToLatex(i)
            
            line = line + " & " + self.dataseries['t_line_n'].valueToLatex(i)
            line = line + " & " + self.dataseries['t_line_pstop2'].valueToLatex(i)
            line = line + " & " + self.dataseries['t_line_pstop4'].valueToLatex(i)
            line = line + " & " + self.dataseries['r_contact_poly'].valueToLatex(i)
            line = line + " & " + self.dataseries['r_contact_n'].valueToLatex(i)
            
            f.write(line+"\\\\\n")
        
        
        f.write("""        \\hline
    \\end{tabular}
\\end{center}""")
        
        f.close()
     
     
     
        
    def exportLatex2(self, path):
        f = open(path+"/histograms_"+self.batch+"/table2.tex", "w")
        f.write("% automatically created table for batch " + self.batch + "\n\n")
        
        f.write(PQC_value.headerToLatex())
        
        f.write("""\\begin{center}
        \\fontsize{4pt}{5pt}\\selectfont
    \\begin{tabular}{ |l|r|r|r|r|r|r|r|r|r|r|r|r| } 
        \\hline
        """+self.shortBatch()+""" & \multicolumn{2}{ c|}{MetClover VdP} & \multicolumn{3}{c|}{P cross-bridge VdP/LW} & \multicolumn{2}{ c|}{Bulk VdP cross} & Vbd & I600 & \multicolumn{3}{c|}{Diode Half CV} \\\\
        & """+" \multicolumn{2}{c|}{"+self.dataseries['vdp_metclo_f'].unit+"} & \multicolumn{2}{c|}{"+self.dataseries['vdp_p_cross_bridge_f'].unit+"} & "+self.dataseries['t_line_p_cross_bridge'].unit+ "& \multicolumn{2}{c|}{"+self.dataseries['vdp_bulk_f'].unit+"} & "+self.dataseries['v_bd'].unit+" & "+self.dataseries['i600'].unit+" & "+self.dataseries['v_fd'].unit+" & \detokenize{"+self.dataseries['rho'].unit+"} & \detokenize{"+self.dataseries['conc'].unit+"""}\\\\
          & forw & rev & forw & rev & & forw & rev & & & Vfd & rho & d. conc \\\\
        \\hline\n
        """)
        
        for i in range(0, len(self.labels)+5):
            if i < len(self.labels):
                line = "        \detokenize{"+self.shortLabel(i)+"}"
            else:
                line = self.dataseries['v_th'].summaryDesciption(i)

            line = line + " & " + self.dataseries['vdp_metclo_f'].valueToLatex(i)
            line = line + " & " + self.dataseries['vdp_metclo_r'].valueToLatex(i)
            line = line + " & " + self.dataseries['vdp_p_cross_bridge_f'].valueToLatex(i)
            line = line + " & " + self.dataseries['vdp_p_cross_bridge_r'].valueToLatex(i)
            line = line + " & " + self.dataseries['t_line_p_cross_bridge'].valueToLatex(i)
            
            line = line + " & " + self.dataseries['vdp_bulk_f'].valueToLatex(i)
            line = line + " & " + self.dataseries['vdp_bulk_r'].valueToLatex(i)
            
            line = line + " & " + self.dataseries['v_bd'].valueToLatex(i)
            line = line + " & " + self.dataseries['i600'].valueToLatex(i)
            line = line + " & " + self.dataseries['v_fd'].valueToLatex(i)
            line = line + " & " + self.dataseries['rho'].valueToLatex(i)
            line = line + " & " + self.dataseries['conc'].valueToLatex(i)
            
            f.write(line+"\\\\\n")
        
        
        f.write("""        \\hline
    \\end{tabular}
\\end{center}""")
        
        f.close()
        
        
    def exportLatex3(self, path):
        f = open(path+"/histograms_"+self.batch+"/table3.tex", "w")
        f.write("% automatically created table for batch " + self.batch + "\n\n")
        
        f.write(PQC_value.headerToLatex())
        
        f.write("""
        \\begin{center}
        \\fontsize{4pt}{5pt}\\selectfont
    \\begin{tabular}{ |l|r|r|r|r|r|r|r|r|r|r|r|r|r|} 
        \\hline
        """+self.shortBatch()+""" & FET & \multicolumn{4}{ c|}{MOS} & GCD & \multicolumn{2}{ c|}{GCD05} & \multicolumn{2}{ c|}{Meander} & \multicolumn{3}{ c|}{Contact chain}\\\\
        
        & """+self.dataseries['v_th'].unit+" & "+self.dataseries['v_fb2'].unit+" & "+self.dataseries['c_acc_m'].unit+" & "+self.dataseries['t_ox'].unit+
        " & \detokenize{"+self.dataseries['n_ox'].unit+"} & "+self.dataseries['i_surf'].unit+" & "+self.dataseries['i_surf05'].unit+" & "+
        self.dataseries['i_bulk05'].unit+" & "+self.dataseries['meander_poly'].unit+" & "+self.dataseries['meander_metal'].unit+" & "+
        self.dataseries['contact_poly'].unit+" & "+self.dataseries['contact_n'].unit+" & "+self.dataseries['contact_p'].unit+"""\\\\
        
         & Vth & Vfb &  C-acc & t-ox & n-ox & i-surf & i-surf05 & i-bulk05 & PolySi & Metal & PolySi & N+ & P \\\\
        \\hline\n
        """)
        
        for i in range(0, len(self.labels)+5):
            if i < len(self.labels):
                line = "        \detokenize{"+self.shortLabel(i)+"}"
            else:
                line = self.dataseries['v_fb2'].summaryDesciption(i)
            line = line + " & " + self.dataseries['v_th'].valueToLatex(i)
            line = line + " & " + self.dataseries['v_fb2'].valueToLatex(i)
            line = line + " & " + self.dataseries['c_acc_m'].valueToLatex(i)
            line = line + " & " + self.dataseries['t_ox'].valueToLatex(i)
            line = line + " & " + self.dataseries['n_ox'].valueToLatex(i)
            line = line + " & " + self.dataseries['i_surf'].valueToLatex(i)
            line = line + " & " + self.dataseries['i_surf05'].valueToLatex(i)
            line = line + " & " + self.dataseries['i_bulk05'].valueToLatex(i)
            
            line = line + " & " + self.dataseries['meander_poly'].valueToLatex(i)
            line = line + " & " + self.dataseries['meander_metal'].valueToLatex(i)
            
            line = line + " & " + self.dataseries['contact_poly'].valueToLatex(i)
            line = line + " & " + self.dataseries['contact_n'].valueToLatex(i)
            line = line + " & " + self.dataseries['contact_p'].valueToLatex(i)

            
            f.write(line+"\\\\\n")
        
        f.write("""        \\hline
    \\end{tabular}
\\end{center}""")
        
        f.close()
        
        
        
        
        

