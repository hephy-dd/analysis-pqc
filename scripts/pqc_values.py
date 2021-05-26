from itertools import islice

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
This Object contains one parameter obtained via the PQC measurements, but for
all samples measured the order in the array is cucial for the correct mapping,
so be careful to reorder everything in the PQC_resultset.dataseries at once
"""
class PQC_Values:
    def __init__(self, name='na', label='na', expected_value=0., unit='',
                 show_multiplier=1e0, stray=0.5, values=None, min_allowed=None,
                 max_allowed=None):
        self.values = []
        self.name = name
        self.label = label
        self.unit = unit
        self.show_multiplier = show_multiplier

        if values is None:
            self.values = []
        else:
            self.values = values

        if min_allowed is None:
            self.min_allowed = expected_value * (1 - stray)
        else:
            self.min_allowed = min_allowed

        if max_allowed is None:
            self.max_allowed = expected_value * (1 + stray)
        else:
            self.max_allowed = max_allowed

        self.stray = stray
        self.expected_value = expected_value

    def __len__(self):
        return len(self.values)

    def __str__(self):
        preview = format(np.array(self.values) * self.show_multiplier)
        return f"{self.name}{preview}{self.unit}"

    def append(self, val):
        self.values.append(val)

    def rearrange(self, indices):
        self.values = [ self.values[indices[i]] for i in range(0, len(indices)) ]

    def getValue(self, index):
        # with multiplier to suit the unit
        if index < len(self.values):
            return self.values[index]*self.show_multiplier
        else:
            stats = self.getStats()
            sel={0: stats.selMed,
                 1: stats.selAvg,
                 2: stats.selStd,
                 3: len(stats.values),
                 4: len(stats.values)/stats.nTot, }
            return sel.get(index-len(self.values), "error")

    def getValueString(self, index, niceText=True):
        if index < len(self.values):
            if np.isnan(self.values[index]) and niceText:
                return "failed"
            elif np.isinf(self.values[index]) and niceText:
                return "---"
            return num2str(self.values[index]*self.show_multiplier, self.expected_value)
        else:
            stats = self.getStats()
            sel={0: num2str(stats.selMed, self.expected_value),
                 1: num2str(stats.selAvg, self.expected_value),
                 2: num2str(stats.selStd, self.expected_value),
                 3: "{}/{}".format(len(stats.values), stats.nTot),
                 4: "{:3.2f}".format(len(stats.values)/stats.nTot), }
            return sel.get(index-len(self.values), "error")

    def getStatus(self, index):
        if index >= len(self.values):
            return 0
        value = self.values[index]*self.show_multiplier
        if np.isinf(value):
            return 5  # inf
        elif np.isnan(value):
            return 4  # nan
        elif value > self.max_allowed:
            return 3  # too high
        elif value < self.min_allowed:
            return 2  # too low
        return 1  # OK

    @staticmethod
    def getStatLabels():
        return ["Median", "Average", "Std dev.", "OK/Tot.", "OK (rel)"]

    def split(self, junk_size):
        kwargs = {
            'name': self.name,
            'label': self.label,
            'expected_value': self.expected_value,
            'unit': self.unit,
            'show_multiplier': self.show_multiplier,
            'stray': self.stray,
        }
        return [PQC_Values(values=i, **kwargs) for i in make_chunks(self.values, junk_size)]

    @params('values, nTot, nNan, nTooHigh, nTooLow, totAvg, totStd, totMed, selAvg, selStd, selMed')
    def getStats(self, min_allowed=None, max_allowed=None):
        if min_allowed is None:
            min_allowed = self.min_allowed
        if max_allowed is None:
            max_allowed = self.max_allowed

        nTot = len(np.array(self.values)[np.invert(np.isinf(np.array(self.values)))])

        selector = np.isfinite(np.array(self.values))

        if np.sum(selector) < 2:
            return np.array([0]), 1, 1, 0, 0, 0, 0, 0, 0, 0, 0

        values = np.array(self.values)[selector]*self.show_multiplier   # filter out nans

        if nTot < 2:
            return values, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0

        totMed = np.median(values)
        totAvg = np.mean(values)
        totStd = np.std(values)

        nNan = nTot - len(values)
        values = values[values < max_allowed]
        nTooHigh = nTot - len(values) - nNan
        values = values[values > min_allowed]
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
    def merge(new, parents, name='na', label='na'):
        values = np.concatenate( [t.values for t in parents])
        return new(name, label, parents[0].expected_value, parents[0].unit, parents[0].show_multiplier, values=values, stray=parents[0].stray)
