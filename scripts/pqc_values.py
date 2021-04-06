
from analysis_pqc import params
import numpy as np

def num2str(num, basenum=None):
    if basenum is None:
        basenum = num
    if basenum < 10:
        ret = "{:4.2f}".format(num)
    else:
        ret = "{:4.1f}".format(num)
    return ret
    

def make_chunks(data, size):
    it = iter(data)

    for i in range(0, len(data), size):
        yield [k for k in islice(it, size)]

"""
This Object contains one parameter obtained via the PQC measurements, but for all samples measured
the order in the array is cucial for the correct mapping, so be careful to reorder everything in the PQC_resultset.dataseries at once
"""
class PQC_Values:
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
        if index < len(self.values):
            return self.values[index]*self.showmultiplier
        else:
            stats = self.getStats()
            sel={0: stats.selMed,
                 1: stats.selAvg,
                 2: stats.selStd,
                 3: len(stats.values),
                 4: len(stats.values)/stats.nTot*100., }
            return sel.get(index-len(self.values), "error")
        
    def getValueString(self, index):
        if index < len(self.values):
            return num2str(self.values[index]*self.showmultiplier, self.expectedValue)
        else:
            stats = self.getStats()
            sel={0: num2str(stats.selMed, self.expectedValue),
                 1: num2str(stats.selAvg, self.expectedValue),
                 2: num2str(stats.selStd, self.expectedValue),
                 3: "{}/{}".format(len(stats.values), stats.nTot),
                 4: "{:2.0f}%".format(len(stats.values)/stats.nTot*100.), }
            return sel.get(index-len(self.values), "error")
            
    def getStatus(self, index):
        if index >= len(self.values):
            return 0
        value = self.values[index]*self.showmultiplier
        if np.isnan(value):
            return 4  # nan
        elif value > self.maxAllowed:
            return 3  # too high
        elif value < self.minAllowed:
            return 2  # too low
        return 1  # OK
        
        
    @staticmethod 
    def getStatLabels():
        return ["Median", "Average", "Std dev.", "OK/Tot.", "OK (rel)"]
        
    def len(self):
        return len(self.values)
        
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
        
        if (len(values)):
            selMed = np.median(values)
            selAvg = np.mean(values)
            selStd = np.std(values)
        else:
            selMed = np.nan
            selAvg = np.nan
            selStd = np.nan
        
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
               "\\def\\okval{\\cellcolor[HTML]{ffffff}}\n\n"
               
               
      
