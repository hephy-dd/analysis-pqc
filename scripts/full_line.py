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

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    parser.add_argument('multibatch', nargs='?', default=None)
    return parser.parse_args()
    
    
def loadBatch(path):
    dirs = glob.glob(os.path.join(path, "*"))
    dirs = [t for t in dirs if "histograms" not in t ]
    
    flutelist = ["PQCFlutesLeft"]
    #flutelist = ["PQCFlutesLeft", "PQCFlutesRight"]
    flutes = flutelist*len(dirs)
    dirs = dirs*len(flutelist)   # we need to double it for the two flutes
    dirs.sort()
    
    pqc_results = PQC_resultset(len(dirs), os.path.basename(os.path.dirname(path+"/")))
    pqc_results.analyze(dirs, flutes)
    
    return pqc_results
    
    
    
def plotTimeline(pqc_batches, path):
    print(path)
    matplotlib.rcParams.update({'font.size': 14})
    
    fig = plt.figure(figsize=(8, 6))
    
    
    gs = gridspec.GridSpec(3, 1)
    
    ax0 = plt.subplot(gs[0])
    plt.title("Polysilicon VdP", fontsize=15)
    for res in pqc_batches:
        stats = res.vdp_poly_tot().getStats()
        res.statusbar(stats, ax0, single=False, start=min(res.timestamps), stop=max(res.timestamps), label=res.batch)
    
    ax1 = plt.subplot(gs[1])   
    plt.title("N+ VdP", fontsize=15)     
    for res in pqc_batches:
        stats = res.vdp_n_tot().getStats()
        res.statusbar(stats, ax1, single=False, start=min(res.timestamps), stop=max(res.timestamps), label=res.batch)
        #res.boxplot(stats, ax0, vertical=True, start=min(res.timestamps), stop=max(res.timestamps), label=res.batch)
    
    ax2 = plt.subplot(gs[2])   
    plt.title("p-stop VdP", fontsize=15)     
    for res in pqc_batches:
        stats = res.vdp_pstop_tot().getStats()
        res.statusbar(stats, ax2, single=False, start=min(res.timestamps), stop=max(res.timestamps), label=res.batch)
        #res.boxplot(stats, ax0, vertical=True, start=min(res.timestamps), stop=max(res.timestamps), label=res.batch)
    

    fig.autofmt_xdate()
    
    fig.tight_layout(h_pad=1.0)
    fig.savefig(path+"timeline.png")
    plt.close()
    
    
    
    
    
def main():
    args = parse_args()
    
    if args.multibatch is None:
        pqc_results = loadBatch(args.path)
        pqc_results.prettyPrint()
        pqc_results.createHistograms(args.path)
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
            pqc_slices.append(res.split(3))
            #res.createHistograms(args.path)
        
        print("loaded "+str(len(pqc_batches))+" batches")
        plotTimeline(pqc_batches, args.path+"histograms/")
        
        
        
        

if __name__ =="__main__":
    main()




