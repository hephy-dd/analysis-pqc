#!/usr/bin/env python3

import pqc_analysis_json as pqc
import matplotlib.pyplot as plt
import matplotlib 
from matplotlib import gridspec
import numpy as np
import os
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

class PQC_value:
    def __init__(self, numrows, name='na', nicename='na', expectedValue=0., unit='', showmultiplier=1e0, stray=0.5, value=None):
        self.value = np.array([0.]*numrows)
        self.name = name
        self.nicename = nicename
        self.unit = unit
        self.showmultiplier = showmultiplier
        self.expectedValue = expectedValue
        self.minAllowed = expectedValue * (1-stray)
        self.maxAllowed = expectedValue * (1+stray)
        self.stray = stray

        if value is not None:
            self.value = value
            
    def rearrange(self, indices):
        self.value = self.value[indices]

    def getValue(self, index):
        # with multiplier to suit the unit
        return self.value[index]*self.showmultiplier
        
    @params('values, nTot, nNan, nTooHigh, nTooLow, totAvg, totStd, totMed, selAvg, selStd, selMed')
    def getStats(self):
        nTot = len(self.value)

        selector = np.isfinite(self.value)
        
        if np.sum(selector) < 2:
            return np.array([0]), 1, 1, 0, 0, 0, 0, 0, 0, 0, 0
            
        values = self.value[selector]*self.showmultiplier   # filter out nans
        
        if nTot < 2:
            return values, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0
        
        totMed = np.median(values)
        totAvg = np.mean(values)
        totStd = np.std(values)
        
        nNan = nTot - len(values)
        values = values[values < self.maxAllowed]
        nTooHigh = nTot - len(values) - nNan
        values = values[values > self.minAllowed]
        nTooLow = nTot - len(values) - nNan - nTooHigh
        
        selMed = np.median(values)
        selAvg = np.mean(values)
        selStd = np.std(values)
        
        return values, nTot, nNan, nTooHigh, nTooLow, totAvg, totStd, totMed, selAvg, selStd, selMed
     
    @classmethod
    def merge(new, parents, name='na', nicename='na'):
        value = np.concatenate( [t.value for t in parents])
        return new(1, name, nicename, parents[0].expectedValue, parents[0].unit, parents[0].showmultiplier, value=value, stray=parents[0].stray)

class PQC_resultset:
    def __init__(self, rows, batchname, dataseries=None):
        self.batch = batchname
        self.labels = ["na"]*rows
        self.flutes = ["na"]*rows
        self.timestamps = ["na"]*rows
        
        self.dataseries = {'xtimestamps': self.timestamps,
                           'xlabels':     self.labels,
                           'xflutes':     self.flutes,}

        self.vdp_poly_f = PQC_value(rows, "vdp_poly", "Polysilicon VdP", 2.4, "kOhm/sq", 1e-3, stray=0.2)
        self.vdp_poly_r = PQC_value(rows, "vdp_poly_rev", "Polysilicon VdP reverse", 2.4, "kOhm/sq", 1e-3, stray=0.2)
        self.dataseries['vdp_poly_f'] = self.vdp_poly_f
        self.dataseries['vdp_poly_r'] = self.vdp_poly_r
        self.vdp_n_f = PQC_value(rows, "vdp_N", "N+ VdP", 35., "Ohm/sq", stray=0.2)
        self.vdp_n_r = PQC_value(rows, "vdp_N_rev", "N+ VdP reverse", 35., "Ohm/sq", stray=0.2)
        self.dataseries['vdp_n_f'] = self.vdp_n_f
        self.dataseries['vdp_n_r'] = self.vdp_n_r
        self.vdp_pstop_f = PQC_value(rows, "vdp_pstop", "P-stop VdP", 19., "kOhm/sq", 1e-3, stray=0.2)
        self.vdp_pstop_r = PQC_value(rows, "vdp_pstop_rev", "P-stop VdP reverse", 19., "kOhm/sq", 1e-3, stray=0.2)
        self.dataseries['vdp_pstop_f'] = self.vdp_pstop_f
        self.dataseries['vdp_pstop_r'] = self.vdp_pstop_r
        
        self.t_line_n = PQC_value(rows, "t_line_n", "Linewidth N+", 35., "um")
        self.dataseries['t_line_n'] = self.t_line_n
        self.t_line_pstop2 = PQC_value(rows, "t_line_pstop2", "Linewidth P-stop 2 Wire", 38., "um")
        self.dataseries['t_line_pstop2'] = self.t_line_pstop2
        self.t_line_pstop4 = PQC_value(rows, "t_line_pstop4", "Linewidth P-stop 4 Wire", 55., "um")
        self.dataseries['t_line_pstop4'] = self.t_line_pstop4
        self.r_contact_n = PQC_value(rows, "r_contact_n", "Rcontact N+", 27., "Ohm")
        self.dataseries['r_contact_n'] = self.r_contact_n
        self.r_contac_poly = PQC_value(rows, "r_contact_poly", "Rcontact polysilicon", 100., "kOhm", 1e-3)
        self.dataseries['r_contac_poly'] = self.r_contac_poly

        self.v_th = PQC_value(rows, "fet", "FET Vth", 4., "V", stray=0.25)
        self.dataseries['v_th'] = self.v_th
        self.vdp_metclo_f = PQC_value(rows, "vdp_met_clover", "Metal Cloverleaf VdP", 25., "mOhm/sq", 1e3)
        self.dataseries['vdp_metclo_f'] = self.vdp_metclo_f
        self.vdp_metclo_r = PQC_value(rows, "vdp_met_clover_rev", "Metal Cloverleaf VdP reverse", 25., "mOhm/sq", 1e3)
        self.dataseries['vdp_metclo_r'] = self.vdp_metclo_r
        self.vdp_p_cross_bridge_f = PQC_value(rows, "vdp_cross_bridge", "Cross Bridge VdP", 1.5, "kOhm/sq", 1e-3)
        self.dataseries['vdp_p_cross_bridge_f'] = self.vdp_p_cross_bridge_f
        self.vdp_p_cross_bridge_r = PQC_value(rows, "vdp_cross_bridge_rev", "Cross Bridge VdP reverse", 1.5, "kOhm/sq", 1e-3)
        self.dataseries['vdp_p_cross_bridge_r'] = self.vdp_p_cross_bridge_r
        self.t_line_p_cross_bridge = PQC_value(rows, "t_line_cb", "Linewidth cross bridge P", 35., "um")
        self.dataseries['t_line_p_cross_bridge'] = self.t_line_p_cross_bridge
        self.v_bd = PQC_value(rows, "v_bd", "Breakdown Voltage", 215., "V")
        self.dataseries['v_bd'] = self.v_bd

        self.i600 = PQC_value(rows, "i600", "I @ 600V", 100., "uA", 1e6, stray=1.)
        self.dataseries['i600'] = self.i600
        self.v_fd = PQC_value(rows, "v_fd", "Full depletion Voltage", 260., "V", stray=0.33)
        self.dataseries['v_fd'] = self.v_fd
        self.rho = PQC_value(rows, "rho", "rho", 1.3, "kOhm cm", 0.1)
        self.dataseries['rho'] = self.rho
        self.conc = PQC_value(rows, "d_conc", "Doping Concentration", 3.5, "* 1E12 cm^-3", 1e-18)
        self.dataseries['conc'] = self.conc

        self.v_fb2 = PQC_value(rows, "v_fb", "Flatband voltage", 2.5, "V", stray=0.33)
        self.dataseries['v_fb2'] = self.v_fb2
        self.t_ox = PQC_value(rows, "t_ox", "Oxide thickness", 0.67, "um", stray=0.33)
        self.dataseries['t_ox'] = self.t_ox
        self.n_ox = PQC_value(rows, "n_ox", "Oxide concentration", 10.5, "* 1E10 cm^-3", 1e-10)
        self.dataseries['n_ox'] = self.n_ox
        self.c_acc_m = PQC_value(rows, "c_acc", "Accumulation capacitance", 85., "pF", 1e12, stray=0.2)
        self.dataseries['c_acc_m'] = self.c_acc_m
        self.i_surf = PQC_value(rows, "i_surf", "Surface current", 8., "pA", -1e12, stray=1)
        self.dataseries['i_surf'] = self.i_surf
        self.i_surf05 = PQC_value(rows, "i_surf05", "Surface current 05", 11., "pA", -1e12, stray=1)
        self.dataseries['i_surf05'] = self.i_surf05
        self.i_bulk05 = PQC_value(rows, "i_bulk05", "Bulk current 05", 0.7, "pA", -1e12, stray=1)
        self.dataseries['i_bulk05'] = self.i_bulk05
        
        self.nvdp_poly_f = PQC_value(rows, "nvdp_poly", "PolySi Swapped VdP", 2.4, "kOhm/sq", -1e-3, stray=0.2)
        self.nvdp_poly_r = PQC_value(rows, "nvdp_poly_rev", "PolySi Swapped VdP reverse", 2.4, "kOhm/sq", -1e-3, stray=0.2)
        self.dataseries['nvdp_poly_f'] = self.nvdp_poly_f
        self.dataseries['nvdp_poly_r'] = self.nvdp_poly_r
        self.nvdp_n_f = PQC_value(rows, "nvdp_N", "N+ Swapped VdP", 35., "Ohm/sq", -1., stray=0.2)
        self.nvdp_n_r = PQC_value(rows, "nvdp_N_rev", "N+ Swapped VdP reverse", 35., "Ohm/sq", -1., stray=0.2)
        self.dataseries['nvdp_n_f'] = self.nvdp_n_f
        self.dataseries['nvdp_n_r'] = self.nvdp_n_r
        self.nvdp_pstop_f = PQC_value(rows, "nvdp_pstop", "P-stop Swapped VdP", 19., "kOhm/sq", -1e-3, stray=0.2)
        self.nvdp_pstop_r = PQC_value(rows, "nvdp_pstop_rev", "P-stop Swapped VdP rev", 19., "kOhm/sq", -1e-3, stray=0.2)
        self.dataseries['nvdp_pstop_f'] = self.nvdp_pstop_f
        self.dataseries['nvdp_pstop_r'] = self.nvdp_pstop_r
        
        if dataseries is not None:
            self.dataseries = dataseries

    def vdp_poly_tot(self):
        return PQC_value.merge([self.vdp_poly_f, self.vdp_poly_r], "vdp_poly_tot", "PolySi VdP both")
    def vdp_n_tot(self):
        return PQC_value.merge([self.vdp_n_f, self.vdp_n_r], "vdp_N_tot", "N+ VdP both")     
    def vdp_pstop_tot(self):
        return PQC_value.merge([self.vdp_pstop_f, self.vdp_pstop_r], "vdp_pstop_tot", "P-stop VdP both")
        
    def nvdp_poly_tot(self):
        return PQC_value.merge([self.nvdp_poly_f, self.nvdp_poly_r], "nvdp_poly_tot", "PolySi Swapped VdP")
    def nvdp_n_tot(self):
        return PQC_value.merge([self.nvdp_n_f, self.nvdp_n_r], "nvdp_N_tot", "N+ Swapped VdP both")
    def nvdp_pstop_tot(self):
        return PQC_value.merge([self.nvdp_pstop_f, self.nvdp_pstop_r], "nvdp_pstop_tot", "P-stop Swapped VdP")
    
    def sortByTime(self):
        order = np.argsort(self.timestamps)
        print(str(order))
        
        for key in self.dataseries:
            if type(self.dataseries[key]) is PQC_value:
                self.dataseries[key].rearrange(order)
            else:
                self.dataseries[key][:] = [self.dataseries[key][i] for i in order]  # we want to keep the original object so references are preserved
        
    # warning: this creates not full copies, only the dict is available then
    def split(self, itemsperclice)
        return None
    
    def analyze(self, dirs, flutes):
        print("dirs: "+str(len(dirs))+"  "+str(len(flutes)))
        self.flutes = flutes
        for i in range(0, len(dirs)):
            self.labels[i] = dirs[i].split("/")[-1]
            
            # sometimes the wrong flute is measured or it'scalled PQCFluteLeft istead of PQCFlutesLeft
            if len(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i]])) == 0:
                if flutes[i] == "PQCFlutesLeft":
                    if len(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=["PQCFluteLeft"])) != 0:
                        flutes[i] = "PQCFluteLeft"
                    else:
                        print("Changed to right flute")
                        flutes[i] = "PQCFlutesRight"
                else:
                    if len(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=["PQCFluteRight"])) != 0:
                        flutes[i] = "PQCFluteRight"
                    else:
                        flutes[i] = "PQCFlutesLeft"
                self.flutes = flutes
            
            x = pqc.find_all_files_from_path(dirs[i], "van_der_pauw")
            if i != []:
                self.timestamps[i] = pqc.get_timestamp(x[-1])
            
            self.vdp_poly_f.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "Polysilicon", "cross"], blacklist=["reverse"]), plotResults=False)
            self.vdp_poly_r.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "Polysilicon", "reverse", "cross"]))

            self.vdp_n_f.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "n", "cross"], blacklist=["reverse"]))
            self.vdp_n_r.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "n", "reverse", "cross"]))
            
            if np.isnan(self.vdp_n_f.value[i]) and np.isnan(self.vdp_n_r.value[i]):
                print("alternate N+ name: "+dirs[i])
                self.vdp_n_f.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "n+", "cross"], blacklist=["reverse"]))
                self.vdp_n_r.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "n+", "reverse", "cross"]))
                print("now: "+str(self.vdp_n_f.value[i])+"/"+str(self.vdp_n_r.value[i]))
            
            
            self.vdp_pstop_f.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "P_stop", "cross"], blacklist=["reverse"]))
            self.vdp_pstop_r.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "P_stop", "reverse", "cross"]))

            self.t_line_n.value[i] = pqc.analyse_linewidth_data(pqc.find_all_files_from_path(dirs[i], "linewidth", whitelist=[flutes[i], "n"], single=True), r_sheet=self.vdp_n_f.value[i], printResults=False, plotResults=False)
            self.t_line_pstop2.value[i] = pqc.analyse_linewidth_data(pqc.find_all_files_from_path(dirs[i], "linewidth", whitelist=[flutes[i], "P_stop", "2_wire"], single=True), r_sheet=self.vdp_pstop_f.value[i], printResults=False, plotResults=False)
            self.t_line_pstop4.value[i] = pqc.analyse_linewidth_data(pqc.find_all_files_from_path(dirs[i], "linewidth", whitelist=[flutes[i], "P_stop", "4_wire"], single=True), r_sheet=self.vdp_pstop_f.value[i], printResults=False, plotResults=False)

            self.r_contact_n.value[i] = pqc.analyse_cbkr_data(pqc.find_all_files_from_path(dirs[i], "cbkr", whitelist=[flutes[i], "n"], single=True), r_sheet=self.vdp_n_f.value[i], printResults=False, plotResults=False)
            self.r_contac_poly.value[i] = pqc.analyse_cbkr_data(pqc.find_all_files_from_path(dirs[i], "cbkr", whitelist=[flutes[i], "Polysilicon"], single=True), r_sheet=self.vdp_poly_f.value[i], printResults=False, plotResults=False)

            self.v_th.value[i] = pqc.analyse_fet_data(pqc.find_all_files_from_path(dirs[i], "fet", whitelist=[flutes[i],], single=True), printResults=False, plotResults=False)

            self.vdp_metclo_f.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "metal", "clover"], blacklist=["reverse"]))
            self.vdp_metclo_r.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "metal", "clover", "reverse"], blacklist=[]))

            self.vdp_p_cross_bridge_f.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "P", "cross_bridge"], blacklist=["reverse"]))
            self.vdp_p_cross_bridge_r.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "P", "cross_bridge", "reverse"]))
            self.t_line_p_cross_bridge.value[i] = pqc.analyse_linewidth_data(pqc.find_all_files_from_path(dirs[i], "linewidth", whitelist=[flutes[i], "P", "cross_bridge"], single=True), r_sheet=self.vdp_p_cross_bridge_f.value[i], printResults=False, plotResults=False)

            self.v_bd.value[i] = pqc.analyse_breakdown_data(pqc.find_all_files_from_path(dirs[i], "breakdown", whitelist=[flutes[i],], single=True), printResults=False, plotResults=False)

            # we want this for FLute_3 and not Flute_1
            self.i600.value[i], dummy = pqc.analyse_iv_data(pqc.find_all_files_from_path(dirs[i], "iv", whitelist=[flutes[i], "3"], single=True), printResults=False, plotResults=False)
            self.v_fd.value[i], self.rho.value[i], self.conc.value[i] = pqc.analyse_cv_data(pqc.find_all_files_from_path(dirs[i], "cv", whitelist=[flutes[i], "3"], single=True), printResults=False, plotResults=False)
            
            dummy, self.v_fb2.value[i], self.t_ox.value[i], self.n_ox.value[i], self.c_acc_m.value[i] = pqc.analyse_mos_data(pqc.find_all_files_from_path(dirs[i], "mos", whitelist=[flutes[i],], single=True), printResults=False, plotResults=False)
            self.i_surf.value[i], dummy = pqc.analyse_gcd_data(pqc.find_all_files_from_path(dirs[i], "gcd", whitelist=[flutes[i],], single=True), printResults=False, plotResults=False)  # only i_surf valid
            self.i_surf05.value[i], self.i_bulk05.value[i] = pqc.analyse_gcd_data(pqc.find_all_files_from_path(dirs[i], "gcd05", whitelist=[flutes[i],], single=True), printResults=False, plotResults=False)  # for i_bulk
            
            self.nvdp_poly_f.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "Polysilicon", "ncross"], blacklist=["reverse"]), plotResults=False)
            self.nvdp_poly_r.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "Polysilicon", "reverse", "ncross"]))

            self.nvdp_n_f.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "n", "ncross"], blacklist=["reverse"]))
            self.nvdp_n_r.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "n", "reverse", "ncross"]))

            self.nvdp_pstop_f.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "P_stop", "ncross"], blacklist=["reverse"]))
            self.nvdp_pstop_r.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "P_stop", "reverse", "ncross"]))
            
            
            
    def prettyPrint(self):
        print("# serial                                 \t  vdp_poly/kOhm/sq       vdp_n/Ohm/sq     vdp_pstop/kOhm/sq   lw_n/um    lw_p2/um   lw_p4/um cbkr_poly/kOhm cbkr_n/Ohm")
        for i in range(0, len(self.labels)):
            line = "{} {}  \t".format(self.labels[i], self.flutes[i])
            line += "{:8.2f} {:8.2f}    ".format(self.vdp_poly_f.getValue(i), self.vdp_poly_r.getValue(i))
            line += "{:8.2f} {:8.2f}    ".format(self.vdp_n_f.getValue(i), self.vdp_n_r.getValue(i))
            line += "{:8.2f} {:8.2f}    ".format(self.vdp_pstop_f.getValue(i), self.vdp_pstop_r.getValue(i))
            line += "{:8.2f} {:8.2f} {:8.2f}     ".format(self.t_line_n.getValue(i), self.t_line_pstop2.getValue(i), self.t_line_pstop4.getValue(i))
            line += "{:8.2f} {:8.2f}".format(self.r_contac_poly.getValue(i), self.r_contact_n.getValue(i))

            print(line)
        
        print("")
        print("")
        print("# serial                                 \t fet       vdp_met-clov       vdp_p-cr-br/kOhm/sq  lw_cb/um  v_bd/V    i600/uA    V_fd/V   rho/kOhm cm   d-conc/1E12 cm^-3")
        print("#                                        \tv_th/V     ")
        for i in range(0, len(self.labels)):
            line = "{} {}  \t".format(self.labels[i], self.flutes[i])
            line += "{:5.2f}  ".format(self.v_th.getValue(i))
            line += "{:9.2E}  {:9.2E}   ".format(self.vdp_metclo_f.getValue(i), self.vdp_metclo_r.getValue(i))

            line += "{:8.2f} {:8.2f}    ".format(self.vdp_p_cross_bridge_f.getValue(i), self.vdp_p_cross_bridge_r.getValue(i))
            line += "{:8.2f}".format(self.t_line_p_cross_bridge.getValue(i))

            line += "{:8.2f}     ".format(self.v_bd.getValue(i))
            line += "{:9.2f}  {:9.2f}  {:7.2f}  {:8.2f}   ".format(self.i600.getValue(i), self.v_fd.getValue(i), self.rho.getValue(i), self.conc.getValue(i))
            print(line)
        print("")
        print("")
        print("# serial                                 \t                    mos                        gcd             gcd05")
        print("#                                        \t v_fb/V    c_acc/pF   t_ox/um n_ox/1E10cm^-2 i_surf/pA  i_surf/pA   i_bulk/pA")
        for i in range(0, len(self.labels)):
            line = "{} {}  \t".format(self.labels[i], self.flutes[i])
            line += "{:8.2f}    {:6.2f}    {:7.3f}  {:9.2f}     ".format(self.v_fb2.getValue(i), self.c_acc_m.getValue(i), self.t_ox.getValue(i), self.n_ox.getValue(i))
            line += "{:8.2f}  {:8.2f}  {:8.2f}    ".format(self.i_surf.getValue(i), self.i_surf05.getValue(i), self.i_bulk05.getValue(i))

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
            
        
    def histogram(self, pqc_values, path, stray=1.4):
        stats = pqc_values.getStats()
        if len(stats.values) == 1 and stats.nNan == 1:
            print("warning: skipping plot due to no valid results: "+pqc_values.name)
            return
        
        fig = plt.figure(figsize=(8, 6))
        
        
        gs = gridspec.GridSpec(1, 2, width_ratios=[10, 1]) 
        ax0 = plt.subplot(gs[0])
        plt.title(self.batch + ": " + pqc_values.nicename + "", fontsize=18)
        ax1 = plt.subplot(gs[1])
                
        ax0.hist(stats.values, bins=20, range=[pqc_values.minAllowed, pqc_values.maxAllowed], facecolor='blue', alpha=1, edgecolor='black', linewidth=1)
        
        
        ax0.set_xlabel(pqc_values.unit)
        ax0.set_ylabel("occurences")
        
        #descNum = "Total number: {}\nShown: {:2.0f}%, {}\nFailed: {:2.0f}%, {}\nToo high: {:2.0f}%, {}\nToo low: {:2.0f}%, {}".format(stats.nTot, len(stats.values)/stats.nTot*1e2, len(stats.values), stats.nNan/stats.nTot*1e2, stats.nNan, stats.nTooHigh/stats.nTot*1e2, stats.nTooHigh, stats.nTooLow/stats.nTot*1e2, stats.nTooLow)
        descNum = "Total: {}\nShown: {}\nFailed: {}\nToo high: {}\nToo low: {}".format(stats.nTot, len(stats.values),  stats.nNan, stats.nTooHigh, stats.nTooLow)
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
        
        for i in [(stats.totAvg, 'purple', 'solid'), (stats.totMed, 'purple', 'dashed'), (stats.selAvg, 'green', 'solid'), (stats.selMed, 'green', 'dashed')]:
            if (i[0] < pqc_values.maxAllowed) and (i[0] > pqc_values.minAllowed):
                ax0.vlines(x = i[0], ymin = 0, ymax = 3, 
                    colors = i[1], linestyles=i[2],
                    label = 'vline_multiple - full height') 
                    

        self.statusbar(stats, ax1)
        
        
        
        fig.tight_layout(h_pad=1.0)
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

        self.histogram(PQC_value.merge([self.nvdp_poly_f, self.nvdp_poly_r], "nvdp_poly_tot", "PolySi Swapped VdP"), histogramDir)
        self.histogram(PQC_value.merge([self.nvdp_n_f, self.nvdp_n_r], "nvdp_N_tot", "N+ Swapped VdP both"), histogramDir)
        self.histogram(PQC_value.merge([self.nvdp_pstop_f, self.nvdp_pstop_r], "nvdp_pstop_tot", "P-stop Swapped VdP"), histogramDir)
        
        
        
        
        
        
        
        
        
        

