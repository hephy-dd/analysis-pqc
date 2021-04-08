#!/usr/bin/env python3

import argparse
import math
import os
import sys
import glob

import matplotlib.pyplot as plt
import matplotlib 
from matplotlib import gridspec
import numpy as np
from pqc_resultset import PQC_resultset
from datetime import timedelta, date
from jinja2 import Template, Environment, FileSystemLoader



def renderTemplates(pqc_resultset):
    # Create the jinja2 environment.
    # Notice the use of trim_blocks, which greatly helps control whitespace.
    template_dir = os.path.join(os.path.dirname(__file__), "templates-enabled")
    
    j2_env = Environment(loader=FileSystemLoader(template_dir),
                         trim_blocks=True)
                         
    templates = glob.glob(os.path.join(template_dir, "*"))
    for f in templates:
        filename = os.path.basename(f)
        rendered_content = j2_env.get_template(os.path.basename(filename)).render(
                batch=pqc_resultset.batch,
                dataseries=pqc_resultset.dataseries)
                
        if "stdout" in filename:
            print(rendered_content)
        else:
            with open(os.path.join(pqc_resultset.outDir, filename), "w") as fh:
                fh.write(rendered_content)
    
    
def plotTimeline(pqc_batches, path, printBatchNumbers=True):
    print(path)
    matplotlib.rcParams.update({'font.size': 14})
    
    fig = plt.figure(figsize=(8, 6))
    
    
    gs = gridspec.GridSpec(3, 1)
    
    ax0 = plt.subplot(gs[0])
    plt.title("Polysilicon VdP", fontsize=15)
    plt.yticks([])
    plt.ylim([0, 1])
    
    ax1 = plt.subplot(gs[1])   
    plt.title("N+ VdP", fontsize=15)
    plt.yticks([])
    plt.ylim([0, 1])
        
    ax2 = plt.subplot(gs[2])   
    plt.title("p-stop VdP", fontsize=15)  
       
    for res in pqc_batches:
        spoly = res.vdp_poly_tot().getStats()
        sn = res.vdp_n_tot().getStats()
        spstop = res.vdp_pstop_tot().getStats()
        lbl = res.batch
        if not printBatchNumbers:
            lbl = ''
        
        start = min(res.timestamps)
        stop = max(res.timestamps)
        cent = start + (stop-start)/2
        
        if not printBatchNumbers:
            start = cent - timedelta(days=1)
            stop = cent + timedelta(days=1)
        
        res.statusbar(spoly, ax0, single=False, start=start, stop=stop, label=lbl)
        res.statusbar(sn, ax1, single=False, start=start, stop=stop, label=lbl)
        res.statusbar(spstop, ax2, single=False, start=start, stop=stop, label=lbl)
    

    fig.autofmt_xdate()
    
    fig.tight_layout(h_pad=1.0)
    fig.savefig(path+"Timeline.png")
    plt.close()
    
    
def plotBoxplot(pqc_batches, path, values=['vdp_poly_f', 'vdp_n_f', 'vdp_pstop_f', 'vdp_poly_r', 'vdp_n_r', 'vdp_pstop_r']):
    print(path)
    matplotlib.rcParams.update({'font.size': 12})
    
    fig = plt.figure(figsize=(8, 6))
    
    gs = gridspec.GridSpec(int(len(values)/2), 2)
    labels = [ b.shortBatch(vpx=False) for b in pqc_batches ]
    
    for i in range(0, len(values)):
        ax = plt.subplot(gs[i])
        plt.title(pqc_batches[0].dataseries[values[i]].nicename, fontsize=13)
        plt.grid(axis='y', linestyle=':')
        ax.set_ylabel(pqc_batches[0].dataseries[values[i]].unit)
        
        data = [ b.dataseries[values[i]].getStats().values for b in pqc_batches ]

        ax.boxplot(data, labels=labels)
        if (i < (len(values)-2)):
            ax.set_xticklabels([])
    
    fig.tight_layout(h_pad=1.0)
    fig.savefig(path+"Boxplot.png")
    plt.close()
    
def vdpPlotBoxplot(pqc_batches, path):
    print(path)
    matplotlib.rcParams.update({'font.size': 14})
    
    fig = plt.figure(figsize=(8, 6))
    
    gs = gridspec.GridSpec(3, 1)
    labels = [ b.shortBatch() for b in pqc_batches ]
    
    ax = plt.subplot(gs[0])
    plt.grid(axis='y', linestyle=':')
    plt.title(pqc_batches[0].vdp_poly_tot().nicename, fontsize=15)
    ax.set_ylabel(pqc_batches[0].vdp_poly_tot().unit)
    data = [ b.vdp_poly_tot().getStats().values for b in pqc_batches ]
    ax.boxplot(data, labels=labels) 
    
    ax = plt.subplot(gs[1])
    plt.grid(axis='y', linestyle=':')
    plt.title(pqc_batches[0].vdp_n_tot().nicename, fontsize=15)
    ax.set_ylabel(pqc_batches[0].vdp_n_tot().unit)
    data = [ b.vdp_n_tot().getStats().values for b in pqc_batches ]
    ax.boxplot(data, labels=labels) 
    
    ax = plt.subplot(gs[2])
    plt.grid(axis='y', linestyle=':')
    plt.title(pqc_batches[0].vdp_pstop_tot().nicename, fontsize=15)
    ax.set_ylabel(pqc_batches[0].vdp_pstop_tot().unit)
    data = [ b.vdp_pstop_tot().getStats().values for b in pqc_batches ]
    ax.boxplot(data, labels=labels) 
    
    fig.tight_layout(h_pad=1.0)
    fig.savefig(path+"vdpBoxplot.png")
    plt.close()
    
    
    
def loadBatch(path, outdir=None, lazy=False, create_plots=False):
    batchname = os.path.basename(os.path.normpath(path))
    print("Batch: "+batchname)
    pqc_results = PQC_resultset(batchname)
    
    if lazy and outdir is not None:
        anatime = os.path.getmtime(pqc_results.analysisFolder(outdir))
        meastime = os.path.getmtime(path)
        if anatime > meastime:
            print("lazy mode: nothing to do")
            exit(0)
    
    pqc_results.prepareAnalysisFolder(outdir)
    pqc_results.analyze(path, outdir, create_plots)
    
    return pqc_results
    
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    parser.add_argument('-m', action='store_true', help='multibatch mode, e. g. for time analysis (experimental)')
    parser.add_argument('-o', default=None, help='override output directory location')
    parser.add_argument('-l', action='store_true', default=None, help='lazy evaluation: skip if the measurement folder is older than analysis folder')
    parser.add_argument('-H', action='store_true', default=None, help='create histograms')
    parser.add_argument('-P', action='store_true', default=None, help='create plots (for each single measurement used)')
    
    #parser.add_argument('-d', action='store_true', default=None, help='create plots with debugging infos inside (e.g. correlation coefficients)')
    args = parser.parse_args()
    
    outdir = args.o
    if outdir is None:
        outdir = args.path
    
    if not args.m:
        pqc_results = loadBatch(args.path, outdir, args.l, args.P)
        
        renderTemplates(pqc_results)
        
        if args.H:
            pqc_results.createHistograms()
        
        
    else:
        print("Multibatch mode - experimental!")
        dirs = glob.glob(os.path.join(args.path, "*"))
        dirs = [t for t in dirs if "histograms" not in t and "VPX" in t ]
        
        pqc_batches = []
        pqc_slices = []
        
        for diri in dirs:
            print("Current dir: "+str(diri))
            res = loadBatch(diri)
            res.sortByTime()
            pqc_batches.append(res)
            pqc_slices.extend(res.split(4))
        #    #res.createHistograms(args.path)
        
        print("loaded "+str(len(pqc_batches))+" batches")
        plotTimeline(pqc_batches, args.path+"histograms/batch")
        plotTimeline(pqc_slices, args.path+"histograms/fine", printBatchNumbers=False)
        
        vdpPlotBoxplot(pqc_batches, args.path+"histograms/")
        plotBoxplot(pqc_batches, args.path+"histograms/a", values=['t_line_n', 'r_contact_n', 't_line_pstop2', 't_line_pstop4', 'r_contact_poly', 'v_th'])
        plotBoxplot(pqc_batches, args.path+"histograms/b", values=['vdp_p_cross_bridge_f', 'vdp_p_cross_bridge_r', 't_line_p_cross_bridge', 'v_bd', 'i600', 'v_fd'])
        plotBoxplot(pqc_batches, args.path+"histograms/c", values=['rho', 'conc', 't_ox', 'n_ox', 'c_acc_m', 'i_surf'])
        
    plt.show()
    
        
        
if __name__ =="__main__":
    main()




