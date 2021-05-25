#!/usr/bin/env python3

import argparse
import os
import glob

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec
from pqc_resultset import PQC_resultset
from datetime import timedelta
from jinja2 import Environment, FileSystemLoader


def renderTemplates(pqc_resultset, templates=None):
    # Create the jinja2 environment.
    # Notice the use of trim_blocks, which greatly helps control whitespace.
    template_dir = os.path.join(os.path.dirname(__file__), "templates")

    j2_env = Environment(loader=FileSystemLoader(template_dir), trim_blocks=True)

    filenames = set()
    for spec in templates:
        for filename in glob.glob(os.path.join(template_dir, spec)):
            # Ignore sub directories
            if os.path.isfile(filename):
                filenames.add(filename)

    for filename in filenames:
        filename = os.path.basename(filename)
        rendered_content = j2_env.get_template(filename).render(
            batch=pqc_resultset.batch,
            dataseries=pqc_resultset.dataseries
        )

        # Handle special stdout templates
        if "stdout" in filename:
            print(rendered_content)
        else:
            output_filename = os.path.join(pqc_resultset.outDir, filename)
            with open(output_filename, "w") as fh:
                print(f"rendering file {output_filename} ... ", end="")
                fh.write(rendered_content)
                print("done.")


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
        cent = start + (stop - start) / 2

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

    gs = gridspec.GridSpec(int(len(values) / 2), 2)
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



def loadBatch(path, outdir=None, lazy=False, create_plots=False, force_eval=False):
    batchname = os.path.basename(os.path.normpath(path))
    print("Batch: "+batchname)
    pqc_results = PQC_resultset(batchname)

    if lazy and outdir is not None:
        try:
            anatime = os.path.getmtime(pqc_results.analysisFolder(outdir))
            meastime = os.path.getmtime(path)

            if anatime > meastime:
                print("lazy mode: nothing to do")
                exit(0)
        except FileNotFoundError:
            print("lazy but first time")

    pqc_results.prepareAnalysisFolder(outdir)
    pqc_results.analyze(path, outdir, create_plots, force_eval)

    return pqc_results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    parser.add_argument('-m', dest='multibatch', action='store_true', help='multibatch mode, e. g. for time analysis (experimental)')
    parser.add_argument('-o', dest='outdir', metavar='DIR', default=None, help='override output directory location')
    parser.add_argument('-l', dest='lazy', action='store_true', default=None, help='lazy evaluation: skip if the measurement folder is older than analysis folder')
    parser.add_argument('-H', dest='histograms', action='store_true', default=None, help='create histograms')
    parser.add_argument('-P', dest='plots', action='store_true', default=None, help='create plots (for each single measurement used)')
    parser.add_argument('-f', dest='force', action='store_true', default=None, help='force evaluating all directories (normally, only directroies with at least one VdP measurement are evaluated to prevent blank lines if the is a wrong file or so)')
    parser.add_argument('-t', dest='templates', metavar='EXPR', action='append', default=[], help='select templates to render (eg. -t*.tex -t*.html -tall.txt)')

    #parser.add_argument('-d', action='store_true', default=None, help='create plots with debugging infos inside (e.g. correlation coefficients)')
    args = parser.parse_args()

    outdir = args.outdir or  args.path

    # Create output directory if not exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if args.multibatch:
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

        print(f"loaded {len(pqc_batches)} batches")
        plotTimeline(pqc_batches, os.path.join(args.path, 'histograms', 'batch'))
        plotTimeline(pqc_slices, os.path.join(args.path, "histograms", "fine"), printBatchNumbers=False)

        vdpPlotBoxplot(pqc_batches, os.path.join(args.path, "histograms") + '/') # TODO
        plotBoxplot(pqc_batches, os.path.join(args.path, "histograms", "a"), values=['t_line_n', 'r_contact_n', 't_line_pstop2', 't_line_pstop4', 'r_contact_poly', 'v_th'])
        plotBoxplot(pqc_batches, os.path.join(args.path, "histograms", "b"), values=['vdp_p_cross_bridge_f', 'vdp_p_cross_bridge_r', 't_line_p_cross_bridge', 'v_bd', 'i600', 'v_fd'])
        plotBoxplot(pqc_batches, os.path.join(args.path, "histograms", "c"), values=['rho', 'conc', 't_ox', 'n_ox', 'c_acc_m', 'i_surf'])
    else:
        pqc_results = loadBatch(args.path, outdir, lazy=args.lazy, create_plots=args.plots, force_eval=args.force)

        renderTemplates(pqc_results, args.templates)

        if args.histograms:
            pqc_results.createHistograms()

    plt.show()

if __name__ =="__main__":
    main()
