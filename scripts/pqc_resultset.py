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


def params(names):
    """Function decorator returning namedtuples."""
    def params(f):
        def params(*args, **kwargs):
           return namedtuple(f.__name__, names)(*f(*args, **kwargs))
        return params
    return params

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
        if value is not None:
            self.value = value

    def getValue(self, index):
        # with multiplier to suit the unit
        return self.value[index]*self.showmultiplier
        
    @params('values, nTot, nNan, nTooHigh, nTooLow, totAvg, totStd, totMed, selAvg, selStd, selMed')
    def getStats(self):
        nTot = len(self.value)
        values = self.value[np.isfinite(self.value)]*self.showmultiplier   # filter out nans
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
        return new(1, name, nicename, parents[0].expectedValue, parents[0].unit, parents[0].showmultiplier, value=value)

class PQC_resultset:
    def __init__(self, rows, batchname):
        self.batch = batchname
        self.labels = ["na"]*rows
        self.flutes = ["na"]*rows

        self.vdp_poly_f = PQC_value(rows, "vdp_poly", "Polysilicon Van-der-Pauw", 2.5, "kOhm/sq", 1e-3)
        self.vdp_poly_r = PQC_value(rows, "vdp_poly_rev", "Polysilicon Van-der-Pauw reverse", 2.5, "kOhm/sq", 1e-3)
        self.vdp_n_f = PQC_value(rows, "vdp_N", "N+ Van-der-Pauw", 35., "Ohm/sq")
        self.vdp_n_r = PQC_value(rows, "vdp_N_rev", "N+ Van-der-Pauw reverse", 35., "Ohm/sq")
        self.vdp_pstop_f = PQC_value(rows, "vdp_pstop", "P-stop Van-der-Pauw", 20., "kOhm/sq", 1e-3)
        self.vdp_pstop_r = PQC_value(rows, "vdp_pstop_rev", "P-stop Van-der-Pauw reverse", 20., "kOhm/sq", 1e-3)
        self.t_line_n = PQC_value(rows, "t_line_n", "Linewidth N+", 35., "um")
        self.t_line_pstop2 = PQC_value(rows, "t_line_pstop2", "Linewidth P-stop 2 Wire", 38., "um")
        self.t_line_pstop4 = PQC_value(rows, "t_line_pstop4", "Linewidth P-stop 4 Wire", 55., "um")
        self.r_contact_n = PQC_value(rows, "r_contact_n", "Rcontact N+", 27., "Ohm")
        self.r_contac_poly = PQC_value(rows, "r_contact_poly", "Rcontact polysilicon", 100., "kOhm", 1e-3)

        self.v_th = PQC_value(rows, "fet", "FET Vth", 3.8, "V")
        self.vdp_metclo_f = PQC_value(rows, "vdp_met_clover", "Metal Cloverleaf VdP", 25., "mOhm/sq", 1e3)
        self.vdp_metclo_r = PQC_value(rows, "vdp_met_clover_rev", "Metal Cloverleaf VdP reverse", 25., "mOhm/sq", 1e3)
        self.vdp_p_cross_bridge_f = PQC_value(rows, "vdp_cross_bridge", "Cross Bridge VdP", 1.5, "kOhm/sq", 1e-3)
        self.vdp_p_cross_bridge_r = PQC_value(rows, "vdp_cross_bridge_rev", "Cross Bridge VdP reverse", 1.5, "kOhm/sq", 1e-3)
        self.t_line_p_cross_bridge = PQC_value(rows, "t_line_cb", "Linewidth cross bridge P", 35., "um")
        self.v_bd = PQC_value(rows, "v_bd", "Breakdown Voltage", 215., "V")

        self.i600 = PQC_value(rows, "i600", "I @ 600V", 50., "uA", 1e6)
        self.v_fd = PQC_value(rows, "v_fd", "Full depletion Voltage", 260., "V")
        self.rho = PQC_value(rows, "rho", "rho", 1.2, "kOhm cm", 0.1)
        self.conc = PQC_value(rows, "d_conc", "Doping Concentration", 3.5, "* 1E12 cm^-3", 1e-18)

        self.v_fb2 = PQC_value(rows, "v_fb", "Flatband voltage", 2.5, "V")
        self.t_ox = PQC_value(rows, "t_ox", "Oxide thickness", 0.67, "um")
        self.n_ox = PQC_value(rows, "n_ox", "Oxide concentration", 10.5, "* 1E10 cm^-3", 1e-10)
        self.c_acc_m = PQC_value(rows, "c_acc", "Accumulation capacitance", 85., "pF", 1e12)
        self.i_surf = PQC_value(rows, "i_surf", "Surface current", 8., "pA", -1e12)
        self.i_surf05 = PQC_value(rows, "i_surf05", "Surface current 05", 11., "pA", -1e12)
        self.i_bulk05 = PQC_value(rows, "i_bulk05", "Bulk current 05", 0.7, "pA", -1e12)
       
    
    
    def analyze(self, dirs, flutes):
        print("dirs: "+str(len(dirs))+"  "+str(len(flutes)))
        self.flutes = flutes
        for i in range(0, len(dirs)):
            self.labels[i] = dirs[i].split("/")[-1]
            
            # sometimes the wron flute is measured
            
            if len(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i]])) == 0:
                if flutes[i] == "PQCFlutesLeft":
                    flutes[i] = "PQCFlutesRight"
                else:
                    flutes[i] = "PQCFlutesLeft"
                self.flutes = flutes
                
            self.vdp_poly_f.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "Polysilicon"], blacklist=["reverse"]), plotResults=False)
            self.vdp_poly_r.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "Polysilicon", "reverse"]))

            self.vdp_n_f.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "n"], blacklist=["reverse"]))
            self.vdp_n_r.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "n", "reverse"]))

            self.vdp_pstop_f.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "P_stop"], blacklist=["reverse"]))
            self.vdp_pstop_r.value[i] = pqc.get_vdp_value(pqc.find_all_files_from_path(dirs[i], "van_der_pauw", whitelist=[flutes[i], "P_stop", "reverse"]))

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
        
        
        
        
    def histogram(self, pqc_values, path, stray=1.4):
    
        #font = {'family' : 'normal',
        #    'weight' : 'bold',
        #    'size'   : 22}
        #matplotlib.rc('font', **font)

        stats = pqc_values.getStats()
        if (len(stats.values) == 0):
            print("warning: skipping plot due to no valid results: "+pqc_values.name)
            return
        
        fig = plt.figure(figsize=(8, 6))
        fig.suptitle(self.batch + " - " + pqc_values.nicename+" - Histogram", fontsize=14)
        
        gs = gridspec.GridSpec(1, 2, width_ratios=[10, 1]) 
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])
                
        ax0.hist(stats.values, bins=20, range=[pqc_values.minAllowed, pqc_values.maxAllowed], facecolor='blue', alpha=0.5)

        
        ax0.set_xlabel(pqc_values.unit)
        ax0.set_ylabel("occurences")
        
        descNum = "Total number: {}\nShown: {:2.0f}%, {}\nFailed: {:2.0f}%, {}\nToo high result: {:2.0f}%, {}\nToo low result: {:2.0f}%, {}".format(stats.nTot, len(stats.values)/stats.nTot*1e2, len(stats.values), stats.nNan/stats.nTot*1e2, stats.nNan, stats.nTooHigh/stats.nTot*1e2, stats.nTooHigh, stats.nTooLow/stats.nTot*1e2, stats.nTooLow)
        fig.text(0.75, 0.85, descNum, bbox=dict(facecolor='red', alpha=0.5), horizontalalignment='right', verticalalignment='top')
        
        if abs(stats.totAvg) < 1e6:
            descStat = "Total avg: {0:8.2f} {4}\nTotal median: {1:8.2f} {4}\nSelected avg: {2:8.2f} {4}\nSelected median: {3:8.2f} {4}".format(stats.totAvg, stats.totMed, stats.selAvg, stats.selMed, pqc_values.unit)
        else:
            descStat = "Total avg: {0:9.2E} {4}\nTotal median: {1:9.2E} {4}\nSelected avg: {2:9.2E} {4}\nSelected median: {3:9.2E} {4}".format(stats.totAvg, stats.totMed, stats.selAvg, stats.selMed, pqc_values.unit)
        fig.text(0.45, 0.85, descStat, bbox=dict(facecolor='yellow', alpha=0.5), horizontalalignment='right', verticalalignment='top')
        
        for i in [(stats.totAvg, 'purple', 'solid'), (stats.totMed, 'purple', 'dashed'), (stats.selAvg, 'green', 'solid'), (stats.selMed, 'green', 'dashed')]:
            if (i[0] < pqc_values.maxAllowed) and (i[0] > pqc_values.minAllowed):
                ax0.vlines(x = i[0], ymin = 0, ymax = 3, 
                    colors = i[1], linestyles=i[2],
                    label = 'vline_multiple - full height') 
                    

        relOK = len(stats.values)/stats.nTot
        relNan = stats.nNan/stats.nTot
        relTooHigh = stats.nTooHigh/stats.nTot
        relTooLow = stats.nTooLow/stats.nTot
        
        p1 = ax1.bar(0, (relOK,), 1, color="green")
        p2 = plt.bar(0, (relNan,), 1, bottom=(relOK,), color="red")
        p3 = plt.bar(0, (relTooHigh,), 1, bottom=(relOK+relNan,), color="orange")
        p4 = plt.bar(0, (relTooLow,), 1, bottom=(relOK+relNan+relTooHigh,), color="yellow")
        
        ax1.text(0, relOK-0.01, 'OK', horizontalalignment='center', verticalalignment='top')
        if relNan > 0.02:
            ax1.text(0, relOK+relNan-0.01, 'Failed', horizontalalignment='center', verticalalignment='top')
        if relTooHigh > 0.02:
            ax1.text(0, relOK+relNan+relTooHigh-0.01, 'High', horizontalalignment='center', verticalalignment='top')
        if relTooLow > 0.02:
            ax1.text(0, relOK+relNan+relTooHigh+relTooLow-0.01, 'Low', horizontalalignment='center', verticalalignment='top')
            
        
        plt.xticks([])
        plt.yticks([])
        plt.ylim([0, 1])
        plt.xlim([-0.5, 0.5])
        
        
        
        #fig.tight_layout()
        fig.savefig(path+"/"+pqc_values.name+"_hist.png")
        plt.close()

        
    def createHistograms(self, path):
        histogramDir = path+"histograms_"+self.batch
        try:
            os.mkdir(histogramDir)
        except OSError:
            files = glob.glob(histogramDir+"/*")
            for f in files:
                os.remove(f)
        self.histogram(self.vdp_poly_f, histogramDir)
        self.histogram(self.vdp_poly_r, histogramDir)
        self.histogram(PQC_value.merge([self.vdp_poly_f, self.vdp_poly_r], "vdp_poly_tot", "Polysilicon Van-der-Pauw both"), histogramDir)
        
        self.histogram(self.vdp_n_f, histogramDir)
        self.histogram(self.vdp_n_r, histogramDir)
        self.histogram(PQC_value.merge([self.vdp_n_f, self.vdp_n_r], "vdp_N_tot", "N+ Van-der-Pauw both"), histogramDir)
        
        self.histogram(self.vdp_pstop_f, histogramDir)
        self.histogram(self.vdp_pstop_r, histogramDir)
        self.histogram(PQC_value.merge([self.vdp_pstop_f, self.vdp_pstop_r], "vdp_pstop_tot", "P-stop Van-der-Pauw both"), histogramDir)
        
        self.histogram(self.vdp_metclo_f, histogramDir)
        self.histogram(self.vdp_metclo_r, histogramDir)
        
        self.histogram(self.vdp_p_cross_bridge_f, histogramDir)
        self.histogram(self.vdp_p_cross_bridge_r, histogramDir)
        
        self.histogram(self.t_line_p_cross_bridge, histogramDir)
        
        self.histogram(self.t_line_n, histogramDir)
        self.histogram(self.t_line_pstop2, histogramDir)
        self.histogram(self.t_line_pstop4, histogramDir)
        self.histogram(self.r_contac_poly, histogramDir)
        self.histogram(self.r_contact_n, histogramDir)

        
        self.histogram(self.v_th, histogramDir)
        self.histogram(self.v_fd, histogramDir)

        self.histogram(self.i600, histogramDir)
        self.histogram(self.v_bd, histogramDir)
        self.histogram(self.rho, histogramDir)
        self.histogram(self.conc, histogramDir)
        
        self.histogram(self.v_fb2, histogramDir)
        self.histogram(self.t_ox, histogramDir)
        self.histogram(self.n_ox, histogramDir)
        self.histogram(self.c_acc_m, histogramDir)
        self.histogram(self.i_surf, histogramDir)
        self.histogram(self.i_surf05, histogramDir)
        self.histogram(self.i_bulk05, histogramDir)
        
        
        
        
        
        
        
        
        

