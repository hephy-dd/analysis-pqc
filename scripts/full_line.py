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

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    parser.add_argument('multibatch', nargs='?', default=None)
    return parser.parse_args()
    
    
def loadBatch(path):
    dirs = glob.glob(os.path.join(path, "*"))
    dirs = [t for t in dirs if "histograms" not in t ]
    
    dirs.sort()
    
    pqc_results = PQC_resultset(os.path.basename(os.path.dirname(path+"/")))
    pqc_results.analyze(dirs)
    
    return pqc_results
    
    
    
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
    
    
    
    
    
def main():
    args = parse_args()
    
    if args.multibatch is None:
        pqc_results = loadBatch(args.path)
        pqc_results.prettyPrint()
        pqc_results.createHistograms(args.path)
        pqc_results.exportLatex(args.path)
    else:
        print("Multibatch mode")
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
        

        for sl in pqc_slices:
            print(str(sl))
            sl.prettyPrint()
        
        print("loaded "+str(len(pqc_batches))+" batches")
        plotTimeline(pqc_batches, args.path+"histograms/batch")
        plotTimeline(pqc_slices, args.path+"histograms/fine", printBatchNumbers=False)
        
        
        

if __name__ =="__main__":
    main()




