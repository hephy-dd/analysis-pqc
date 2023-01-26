"""Set of analysis function for PQC measurements."""

import warnings
import traceback
from collections import namedtuple

import numpy as np

from scipy.interpolate import CubicSpline
from scipy.stats import linregress
import scipy.signal

__version__ = '0.7.1'

__all__ = [
    'STATUS_NONE',
    'STATUS_PASSED',
    'STATUS_FAILED',
    'CARRIER_ELECTRONS',
    'CARRIER_HOLES',
    'analyse_iv',
    'analyse_cv',
    'analyse_mos',
    'analyse_gcd',
    'analyse_fet',
    'analyse_van_der_pauw',
    'analyse_cross',
    'analyse_linewidth',
    'analyse_cbkr',
    'analyse_contact',
    'analyse_meander',
    'analyse_breakdown',
    'analyse_capacitor'
]

## Constants
## ------------------------------------

STATUS_NONE = 'none'
STATUS_PASSED = 'passed'
STATUS_FAILED = 'failed'

CARRIER_ELECTRONS = 'electrons'
CARRIER_HOLES = 'holes'

## Helper functions
## ------------------------------------

def params(names):
    """Function decorator returning namedtuples."""
    def params(f):
        def params(*args, **kwargs):
           return namedtuple(f.__name__, names)(*f(*args, **kwargs))
        return params
    return params


@params('a, b, x_fit, spl_dev, status, r_value')
def line_regr_with_cuts(x, y, cut_param, debug=False):
    """
    Linear Regression with Cuts:
    - Normalise data set
    - Get 1st derivate of 2nd order spline fit
    - Only use data points with local slope exceeding cut_param

    Parameters:
    x ... x
    y ... y
    cut_param ... used to cut on 1st derivative of x axis

    Returns:
    i_max ... max. current
    i_800 ... current @ 800V
    i_600 ... current @ 600V
    """

    # init
    r_value = a = b = x_fit = spl_dev = -1
    status = STATUS_NONE

    # get spline fit, requires strictlty increasing array
    y_norm = y / np.max(y)
    x_norm = np.arange(len(y_norm))
    spl = CubicSpline(x_norm, y_norm)
    spl_dev = spl(x_norm, 1)

    # only use data points if local slope is above cut_param
    idx_fit = [i for i in range(len(spl_dev)) if (abs(spl_dev[i]) > cut_param)]

    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            x_fit = x[idx_fit[0]:idx_fit[-1] + 1]
            y_fit = y[idx_fit[0]:idx_fit[-1] + 1]
            a, b, r_value, p_value, std_err = scipy.stats.linregress(x_fit, y_fit)
            status = STATUS_PASSED
        except np.RankWarning:
            print("The array has too few data points. Try changing the cut_param parameter.")
            status = STATUS_FAILED
        except (ValueError, TypeError, IndexError):
            print("The array seems empty. Try changing the cut_param parameter.")
            status = STATUS_FAILED

    return a, b, x_fit, spl_dev, status, r_value


## Main Analysis Functions
## ------------------------------------

@params('v_max, i_max, i_800, i_600, i300, status')
def analyse_iv(v, i, debug=False):
    """
    Diode IV: Extract current in standard situation.

    Parameters:
    v ... voltage
    i ... current

    Returns:
    i_max ... max. current
    i_800 ... current @ 800V
    i_600 ... current @ 600V
    i_300 ... current @ 300V
    """

    ## init
    v_max = i_max = i_800 = i_600 = i_300 = np.nan
    status = STATUS_NONE

    ## init
    idx_maxv = np.argmax(np.abs(v))
    idx_maxi = np.argmax(np.abs(i))
    v_max = v[idx_maxv]
    i_max = i[idx_maxv]
    i_800 = np.array([i[k] for k in range(len(v)) if np.abs(v[k]) == 800])
    i_600 = np.array([i[k] for k in range(len(v)) if np.abs(v[k]) == 600])
    i_300 = np.array([i[k] for k in range(len(v)) if np.abs(v[k]) == 300])

    if len(i_800) != 1:
        i_800 = np.nan
    else:
        i_800 = i_800[0]
    if len(i_600) != 1:
        i_600 = np.nan
    else:
        i_600 = i_600[0]
    if len(i_300) != 1:
        i_300 = np.nan
    else:
        i_300 = i_300[0]

    status = STATUS_PASSED

    return v_max, i_max, i_800, i_600, i_300, status


@params('v_dep1, v_dep2, rho, conc, a_rise, b_rise, v_rise, a_const, b_const, v_const, spl_dev, status')
def analyse_cv(v, c, area=1.56e-6, carrier='electrons', cut_param=0.008, max_v=500, savgol_windowsize=None, min_correl=0.1, debug=False):
    """
    Diode CV: Extract depletion voltage and resistivity.

    Parameters:
    v ... voltage
    c ... capacitance
    area ... implant size in [m^2] - defaults to quarter
    carrier ... majority charge carriers ['holes', 'electrons']
    cut_param ... used to cut on 1st derivative to id voltage regions
    max_v ... for definition of fit region, only consider voltages < max_v
    savgol_windowsize ... number of points to calculate the derivative, needs to be odd
    min_correl ... minimum correlation coefficient to say that it worked

    Returns:
    v_dep1 ... full depletion voltage via inflection
    v_dep2 ... full depletion voltage via intersection
    rho ... resistivity
    conc ... bulk doping concentration
    """

    # init
    v_dep1 = v_dep2 = rho = conc = np.nan
    a_rise = b_rise = a_const = b_const = np.nan
    v_rise = []
    v_const = []
    status = STATUS_NONE

    if savgol_windowsize is None:
        # savgol_windowsize = int(len(c) / 40 + 1) * 2 + 1  # a suitable off windowsie - making 20 windows along the whole measurement
        savgol_windowsize = int(len(c) / 30 + 1) * 2 + 1  # a suitable off windowsie - making 15 windows along the whole measurement
        # the window size needs to be an odd number, therefore this strange calculation

    # invert and square
    c = 1. / c**2

    # get spline fit, requires strictlty increasing array
    y_norm = c / np.max(c)
    x_norm = np.arange(len(y_norm))
    # spl = CubicSpline(x_norm, y_norm)
    # spl_dev = spl(x_norm, 1)
    spl_dev = scipy.signal.savgol_filter(y_norm, window_length=savgol_windowsize, polyorder=1, deriv=1)

    # for definition of fit region, only consider voltages < max_v
    idv_max = max([i for i,a in enumerate(v) if abs(a) < max_v])
    spl_dev = spl_dev[:idv_max]

    idx_rise = []
    idx_const = []

    with warnings.catch_warnings():
        warnings.filterwarnings('error')

        try:
            # get regions for indexing
            idx_rise = [i for i in range(2, len(spl_dev-1)) if ((spl_dev[i]) > cut_param)]  # the first and last value seems to be off sometimes
            idx_const = [i for i in range(2, len(spl_dev-1)) if ((spl_dev[i]) < cut_param) and i > idx_rise[-1]]

            v_rise = v[idx_rise[0]:idx_rise[-1] + 1]
            v_const = v[idx_const[0]:idx_const[-1] + 1]
            c_rise = c[idx_rise[0]:idx_rise[-1] + 1]
            c_const = c[idx_const[0]:idx_const[-1] + 1]

            # line fits to each region
            a_rise, b_rise, r_value_rise, p_value_rise, std_err_rise = scipy.stats.linregress(v_rise, c_rise)
            a_const, b_const, r_value_const, p_value_const, std_err_const = scipy.stats.linregress(v_const, c_const)
            #print("c raise {:5.3f}, const {:5.3f}".format(r_value_rise,r_value_const))

            if carrier == CARRIER_HOLES:
                mu = 450 * 1e-4
            elif carrier == CARRIER_ELECTRONS:
                mu = 1350 * 1e-4
            else:
                raise ValueError('Not a valid type of majority carrier.')

            # full depletion voltage via max. 1st derivative
            v_dep1 = v[np.argmax(spl_dev)]

            # full depletion via intersection
            v_dep2 = (b_const - b_rise) / (a_rise - a_const)

            # rest
            conc = 2. / (1.6e-19 * 11.9 * 8.854e-12 * a_rise * area**2)
            rho = 1. / (mu * 1.6e-19 * conc)
            status = STATUS_PASSED

            if abs(r_value_rise) < min_correl or abs(r_value_const) < min_correl or status == STATUS_FAILED:
                status = STATUS_FAILED

            #print("v_rise: "+str(v_rise))
            #print("r rise: "+str(r_value_rise))
            #print("r const: "+str(r_value_const))
            #print("v_const: "+str(v_const))


        except np.RankWarning:
            status = STATUS_FAILED
            print("The array has too few data points. Try changing the cut_param parameter.")

        except (ValueError, TypeError, IndexError):
            status = STATUS_FAILED
            #print("The array seems empty. Try changing the cut_param parameter.")

        if status == STATUS_FAILED:
            #print("The fit didn't work as expected, returning nan")
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, STATUS_FAILED
    return v_dep1, v_dep2, rho, conc, a_rise, b_rise, v_rise, a_const, b_const, v_const, spl_dev, status


@params('v_fb1, v_fb2, c_acc, c_inv, t_ox, n_ox, a_acc, b_acc, v_acc, a_dep, b_dep, v_dep, a_inv, b_inv, v_inv,  spl_dev, status')
def analyse_mos(v, c, cut_param=0.02, debug=False, min_r_value=0.4):
    """
    Metal oxide Capacitor: Extract flatband voltage, oxide thickness and charge density.

    Parameters:
    v ... voltage (V)
    c ... capacitance (F)
    cut_param ... used to cut on 1st derivative to id voltage regions


    Returns:
    v_fb1 ... flatband voltage via inflection (V)
    v_fb2 ... flatband voltage via intersection (V)
    t_ox ... oxide thickness (um)
    n_ox ... oxide charge density (cm^-2)
    """

    ## init
    v_fb1 = v_fb2 = t_ox = n_ox = np.nan
    a_acc = b_acc = a_dep = b_dep = a_inv = b_inv = spl_dev = np.nan
    v_dep = []
    v_inv = []
    v_acc = []
    status = STATUS_NONE

    # take average of last 5 samples for accumulation and inversion capacitance
    c_acc = np.mean(c[-5:])
    c_inv = np.mean(c[:5])

    # get spline fit, requires strictlty increasing array
    y_norm = c / np.max(c)
    x_norm = np.arange(len(y_norm))
    spl = CubicSpline(x_norm, y_norm)
    spl_dev = spl(x_norm, 1)

    # get regions for indexing
    idx_acc = [i for i in range(len(spl_dev)) if (abs(spl_dev[i]) < cut_param and v[i] > v[np.argmax(spl_dev)])]
    idx_dep = [i for i in range(len(spl_dev)) if (v[i] > v[np.argmax(spl_dev)] - 0.25 and v[i] < v[np.argmax(spl_dev)] + 0.25)]
    idx_inv = [i for i in range(len(spl_dev)) if (abs(spl_dev[i]) < cut_param and v[i] < v[np.argmin(spl_dev)])]

    with warnings.catch_warnings():
        warnings.filterwarnings('error')

        try:
            v_acc = v[idx_acc[0]:idx_acc[-1] + 1]
            v_dep = v[idx_dep[0]:idx_dep[-1] + 1]
            v_inv = v[idx_inv[0]:idx_inv[-1] + 1]
            c_acc = c[idx_acc[0]:idx_acc[-1] + 1]
            c_dep = c[idx_dep[0]:idx_dep[-1] + 1]
            c_inv = c[idx_inv[0]:idx_inv[-1] + 1]

            # line fits to each region
            a_acc, b_acc, r_value_acc, p_value, std_err = scipy.stats.linregress(v_acc, c_acc)
            a_dep, b_dep, r_value_dep, p_value, std_err = scipy.stats.linregress(v_dep, c_dep)
            a_inv, b_inv, r_value_inv, p_value, std_err = scipy.stats.linregress(v_inv, c_inv)

            # print("  r: "+str(r_value_acc)+"  "+str(r_value_dep)+"  "+str(r_value_inv))

            if (np.abs(np.array([r_value_acc, r_value_dep, r_value_inv])) > min_r_value).all():
                # print("yes")
                # flatband voltage via inflection
                v_fb1 = v[np.argmax(spl_dev)]

                # flatband voltage via intersection
                v_fb2 = (b_acc - b_dep) / (a_dep - a_acc)

                # note 1: Phi_MS of -0.69V is used as standard value, this correpsonds to a p-type bulk doping of 5e12 cm^-3
                # note 2: We apply the bias voltage to the backplane while keeping the gate to ground, V_fb is therefore positive
                n_ox = np.mean(c_acc) / (1.602e-19 * (0.1290**2)) * (0.69 + v_fb2)
                t_ox = 3.9 * 8.85e-12 * (0.001290**2) / np.mean(c_acc) * 1e6
                status = STATUS_PASSED

        except np.RankWarning:
            status = STATUS_FAILED
            print("rank warning")

        except (ValueError, TypeError, IndexError):
            #traceback.print_exc()
            status = STATUS_FAILED
            #print("errrrrrrr"+str(e))

    return v_fb1, v_fb2, c_acc, c_inv, t_ox, n_ox, a_acc, b_acc, v_acc, a_dep, b_dep, v_dep, a_inv, b_inv, v_inv, spl_dev, status


@params('i_surf, i_bulk, i_acc, i_dep, i_inv, v_acc, v_dep, v_inv, i_acc_relstd, i_dep_relstd, i_inv_relstd, spl_dev, status')
def analyse_gcd(v, i, cut_param=0.01, debug=False, maxreldev=0.01):
    """
    Gate Controlled Diode: Generation currents.

    Parameters:
    v ... voltage
    i ... current
    cut_param ... used to cut on 1st derivative to id voltage regions
    maxreldev ... maximum relative (to the abs max in the three regions) standart deviation to consider measurement as good

    Returns:
    i_surf ... surface generation current
    i_bulk ... bulk generation current
    """

    # init
    i_surf = i_bulk = np.nan
    i_acc = i_dep = i_inv = spl_dev = np.nan
    status = STATUS_NONE
    v_acc = []
    v_dep = []
    v_inv = []
    # get spline fit, requires strictlty increasing array
    y_norm = np.abs(i) / np.max(np.abs(i))
    x_norm = np.arange(len(y_norm))

    spl = CubicSpline(x_norm, y_norm)
    spl_dev = spl(x_norm, 1)

    i_acc_relstd = np.nan
    i_dep_relstd = np.nan
    i_inv_relstd = np.nan


    # get regions for indexing
    try:
        vmin = v[np.argmin(i)]
        idx_acc = [i for i in range(len(spl_dev)) if (abs(spl_dev[i]) < 0.03 and v[i] < (vmin - 4.5))]
        idx_dep = [i for i in range(len(spl_dev)) if (abs(spl_dev[i]) < 0.01 and v[i] > (vmin - 2.5) and v[i] < (vmin + 2.5))]
        idx_inv = [i for i in range(len(spl_dev)) if (abs(spl_dev[i]) < 0.01 and v[i] > (vmin + 4.5))]
        v_acc = v[idx_acc[0]:idx_acc[-1] + 1]
        v_dep = v[idx_dep[0]:idx_dep[-1] + 1]
        v_inv = v[idx_inv[0]:idx_inv[-1] + 1]
        i_acc = i[idx_acc[0]:idx_acc[-1] + 1]
        i_dep = i[idx_dep[0]:idx_dep[-1] + 1]
        i_inv = i[idx_inv[0]:idx_inv[-1] + 1]

        if (len(v_acc) == 0):
            v_acc = v[1:6]
            i_acc = i[1:6]
        if (len(v_dep) == 0):
            v_dep = v[np.argmin(i):(np.argmin(i) + 5)]
            i_dep = i[np.argmin(i):(np.argmin(i) + 5)]
        if (len(v_acc) == 0):
            v_inv = v[-5:]
            i_inv = i[-5:]

        # the selection above is not stable
        # until this is fixed stay with a simpler selection
        v_acc = v[1:6]
        i_acc = i[1:6]
        v_dep = v[np.argmin(i):(np.argmin(i) + 5)]
        i_dep = i[np.argmin(i):(np.argmin(i) + 5)]
        v_inv = v[-5:]
        i_inv = i[-5:]

        i_acc_avg = np.mean(i_acc)
        i_dep_avg = np.mean(i_dep)
        i_inv_avg = np.mean(i_inv)

        i_max = np.max(np.abs([i_dep_avg, i_inv_avg, i_acc_avg]))

        i_acc_relstd = np.std(i_acc)/i_max
        i_dep_relstd = np.std(i_dep)/i_max
        i_inv_relstd = np.std(i_inv)/i_max

        # surface and bulk generation current
        i_surf = i_dep_avg - i_inv_avg
        i_bulk = i_acc_avg - i_inv_avg

        if (np.array([i_acc_relstd, i_dep_relstd, i_inv_relstd]) > maxreldev).any():
            i_surf = np.nan
            i_bulk = np.nan

        if (np.array([i_acc_avg, i_dep_avg, i_inv_avg]) > 1e-3).any():  # electrometer overange condition
            i_surf = np.nan
            i_bulk = np.nan

        status = STATUS_PASSED

    except (ValueError, TypeError, IndexError):
        status = STATUS_FAILED
        i_surf = np.nan
        i_bulk = np.nan

    return i_surf, i_bulk, i_acc, i_dep, i_inv, v_acc, v_dep, v_inv, i_acc_relstd, i_dep_relstd, i_inv_relstd, spl_dev, status


@params('v_th, a, b, spl_dev, status')
def analyse_fet(v, i, debug=False, numDev=6, thrMultDev=0.33):
    """
    Field Effect Transistor: Threshold voltage.

    Parameters:
    v ... voltage
    i ... current

    Returns:
    v_th ... threshold voltage via tangent

    sanity check: at least numDev points after the maximum derivative
     must be higher than the maximum times thrMultDev. then the measurement is considered as good
    """

    # init
    v_th = np.nan
    a = b = spl_dev = -1
    status = STATUS_NONE

    # get spline fit, requires strictlty increasing array
    i = i - i[0]  # we are havng offset problems
    y_norm = i / np.max(np.abs(i))
    x_norm = np.arange(len(y_norm))
    spl = CubicSpline(x_norm, y_norm)
    spl_dev = spl(x_norm, 1)

    # get tangent at max. of 1st derivative
    maximum = np.argmax(spl_dev)
    i_0 = i[maximum]
    v_0 = v[maximum]
    a = (i[maximum] - i[maximum - 1]) / (v[maximum] - v[maximum - 1])
    b = i_0 - a*v_0

    # threshold voltage via tangent
    if a:
        v_th = - b / a
    status = STATUS_PASSED

    if (spl_dev[maximum:(maximum+numDev)] > spl_dev[maximum]*thrMultDev).all() and maximum > numDev:
        return v_th, a, b, spl_dev, status
    return np.nan, a, b, spl_dev, status


@params('r_sheet, a, b, x_fit, spl_dev, status, r_value')
def analyse_van_der_pauw(i, v, cut_param=1e-5, debug=False):
    """
    Van der Pauw: Extract sheet resistance.

    Parameters:
    i ... current
    v ... voltage
    cut_param ... used to cut on 1st derivative to id voltage regions

    Returns:
    r_sheet ... resistance per square
    """

    a, b, x_fit, spl_dev, status, r_value = line_regr_with_cuts(i, v, cut_param, debug)
    r_sheet = np.pi / np.log(2) * a
    return r_sheet, a, b, x_fit, spl_dev, status, r_value


@params('r_sheet, a, b, x_fit, spl_dev, status')
def analyse_cross(i, v, cut_param=1e-5, debug=False):
    """
    Cross: Extract sheet resistance.

    Parameters:
    i ... current
    v ... voltage
    cut_param ... used to cut on 1st derivative to id voltage regions

    Returns:
    r_sheet ... resistance per square
    """

    a, b, x_fit, spl_dev, status, r_value = line_regr_with_cuts(i, v, cut_param, debug)
    r_sheet = np.pi / np.log(2) * a

    return r_sheet, a, b, x_fit, spl_dev, status


@params('t_line, a, b, x_fit, spl_dev, r_value, status')
def analyse_linewidth(i, v, r_sheet=np.nan, cut_param=1e-5, min_correlation=0.99, debug=False):
    """
    Linewidth: Extract linewidth.

    Parameters:
    i ... current
    v ... voltage
    r_sheet ... sheet resistance
    cut_param ... used to cut on 1st derivative to id voltage regions

    Returns:
    t_line ... linewidth in [um]
    """

    a, b, x_fit, spl_dev, status, r_value = line_regr_with_cuts(i, v, cut_param, debug)
    if abs(r_value) < min_correlation:
        return np.nan, np.nan, np.nan, x_fit, spl_dev, r_value, status


    t_line = r_sheet * 128.5 * 1. / a

    return t_line, a, b, x_fit, spl_dev, r_value, status


@params('r_contact, a, b, x_fit, spl_dev, r_value, status')
def analyse_cbkr(i, v, r_sheet=-1, cut_param=1e-5, debug=False):
    """
    Cross Bridge Kelvin Resistance Structure: Extract contact resistance.

    Parameters:
    i ... current
    v ... voltage
    r_sheet ... sheet resistance
    cut_param ... used to cut on 1st derivative to id voltage regions

    Returns:
    r_contact ... contact resistance
    """

    a, b, x_fit, spl_dev, status, r_value = line_regr_with_cuts(i, v, cut_param, debug)

    if r_sheet == -1:
        r_contact = -1
    else:
        # note: The contact isn't symmetric. It's 12.5 by 13.5 um. Solution for now is to use 13 um.
        d = 13  # contact size
        w = 33  # diffusion width
        r_contact = a - (4 * r_sheet * d**2) / (3 * w**2) * (1 + d/(2 * w - 2 * d))

    return r_contact, a, b, x_fit, spl_dev, r_value, status


@params('r_contact, a, b, x_fit, spl_dev, status, r_value')
def analyse_contact(i, v, cut_param=1e-5, debug=False):
    """
    Contact Chain: Extract metal-implant contact resistance.

    Parameters:
    i ... current
    v ... voltage
    r_sheet ... sheet resistance
    cut_param ... used to cut on 1st derivative to id voltage regions

    Returns:
    r_contact ... contact resistance
    """

    a, b, x_fit, spl_dev, status, r_value = line_regr_with_cuts(i, v, cut_param, debug)
    r_contact = a

    return r_contact, a, b, x_fit, spl_dev, status, r_value


@params('r, status, r_value')
def analyse_meander(i, v, cut_param=1e-5, debug=False):
    """
    Meander: Calculates specific resistance per square.

    Parameters:
    i ... current
    v ... voltage
    w ... strip width, use [5, 10] for [polysilicon, metal]
    nsq ... number of squares, use [476, 12853] for [polysilicon, metal]

    Returns:
    rho_sq ... specific resistance per square
    """

    status = STATUS_PASSED

    a, b, x_fit, spl_dev, status, r_value = line_regr_with_cuts(i, v, cut_param, debug)
    r = a

    return r, status, r_value


@params('v_bd, status')
def analyse_breakdown(v, i, debug=False):
    """
    Breakdown: Get oxide breakdown.

    Parameters:
    v ... voltage
    i ... current

    Returns:
    v_bd  ... breakdown voltage
    """

    status = STATUS_PASSED
    v_bd = np.nan
    if len(v) > 0:
        v_bd = v[-1]

    return v_bd, status


@params('c_mean, c_median, d, status')
def analyse_capacitor(v, c, debug=False):
    """
    Test capacitors: Get mean capacitance.

    Parameters:
    v ... voltage
    c ... capacitance

    Returns:
    c_mean ... mean capacitance
    c_median ... median capacitance
    """

    #print(str(c))

    status = STATUS_PASSED
    c_mean = np.mean(c)
    c_median = np.median(c)

    #oxide thickness
    #  d = eps_r * epx_0 * area / capacitance
    # area: (130 um)^2 = 16.9e-9 m^2
    d = 3.9 * 8.85e-12 * 16.9 * 1e-9 / c_median

    return c_mean, c_median, d, status
