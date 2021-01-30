#!/usr/bin/env python3

import argparse
import math
import os
import sys
import glob

import matplotlib.pyplot as plt
import numpy as np
from pqc_resultset import PQC_resultset

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    return parser.parse_args()
    
    
def main():
    args = parse_args()

        
    dirs = glob.glob(os.path.join(args.path, "*"))
    dirs = [t for t in dirs if "histograms" not in t ]
    
    flutelist = ["PQCFlutesLeft"]
    #flutelist = ["PQCFlutesLeft", "PQCFlutesRight"]
    flutes = flutelist*len(dirs)
    dirs = dirs*len(flutelist)   # we need to double it for the two flutes
    dirs.sort()

    pqc_results = PQC_resultset(len(dirs), os.path.basename(os.path.dirname(args.path)))
    pqc_results.analyze(dirs, flutes)
    pqc_results.prettyPrint()
    
    pqc_results.createHistograms(args.path)
    

if __name__ =="__main__":
    main()




