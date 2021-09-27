#!/usr/bin/env python3
import pdb # DEBUG

import argparse
import glob
import os
import sys
from collections.abc import Iterable
from datetime import timedelta

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec
from matplotlib.figure import Figure

from jinja2 import Environment, FileSystemLoader
import yaml

from pqc_resultset import PQC_resultset


def create_dir(dirname: str) -> None:
    """Create directory if not already exists."""
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def render_templates(pqc_resultset: PQC_resultset, templates: Iterable) -> None:
    """Render templates using a PQC resultset, templates is an iterable of glob
    statements, eg. ['*.xml', '*.txt'].
    """
    template_dir = os.path.join(os.path.dirname(__file__), "templates")

    # Create the jinja2 environment.
    # Notice the use of trim_blocks, which greatly helps control whitespace.
    j2_env = Environment(loader=FileSystemLoader(template_dir), trim_blocks=True)

    filenames = set()
    for spec in templates:
        for filename in glob.glob(os.path.join(template_dir, spec)):
            # Ignore sub directories
            if os.path.isfile(filename):
                filenames.add(filename)


    for filename in filenames:

        basename = os.path.basename(filename)
        template_id, extension = os.path.splitext(basename)

        # Skip xml templates of tests which were not carried out
        is_xml_template = extension == '.xml'
        is_valid_template = any(key for key in pqc_resultset.rawdata.keys() if key in basename)
        if is_xml_template and not is_valid_template:
            print(f"skipping XML template: {filename}")
            continue
        rendered_content = j2_env.get_template(basename).render(
            batch=pqc_resultset.batch,
            dataseries=pqc_resultset.dataseries,
            histograms=pqc_resultset.histograms,
            rawdata=pqc_resultset.rawdata[template_id]
        )

        # Handle special stdout templates
        if "stdout" in basename:
            print(rendered_content)
        else:
            create_dir(pqc_resultset.output_dir)
            if is_xml_template: output_filename = os.path.join(pqc_resultset.output_dir, pqc_resultset.rawdata[template_id].out_file_name)
            else: output_filename = os.path.join(pqc_resultset.output_dir, basename)
            
            with open(output_filename, "w") as fh:
                print(f"rendering file {output_filename} ... ", end="", flush=True)
                fh.write(rendered_content)
                print("done.")


def plot_timeline(pqc_batches: list, filename: str,
                  show_batch_numbers: bool = True) -> None:
    """Render timeline figure for batches."""
    print(f"rendering timeline {filename} ...")

    rcParams.update({'font.size': 14})

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
        spoly = res.vdp_poly_tot().get_stats()
        sn = res.vdp_n_tot().get_stats()
        spstop = res.vdp_pstop_tot().get_stats()
        lbl = res.batch
        if not show_batch_numbers:
            lbl = ''

        if not res.timestamps:
            print(f"plot_timeline: skipping empty dataseries {lbl}")
            continue

        start = min(res.timestamps)
        stop = max(res.timestamps)
        cent = start + (stop - start) / 2

        if not show_batch_numbers:
            start = cent - timedelta(days=1)
            stop = cent + timedelta(days=1)

        res.statusbar(spoly, ax0, single=False, start=start, stop=stop, label=lbl)
        res.statusbar(sn, ax1, single=False, start=start, stop=stop, label=lbl)
        res.statusbar(spstop, ax2, single=False, start=start, stop=stop, label=lbl)

    fig.autofmt_xdate()
    fig.tight_layout(h_pad=1.0)
    plot_save(fig, filename)
    plt.close()


def plot_boxplot(pqc_batches: list, filename: str, keys: list = None) -> None:
    """Render boxplot figure for selected keys of dataseries."""
    print(f"rendering boxplot {filename} ...")

    if not pqc_batches:
        return

    if keys is None:  # TODO
        keys = ['vdp_poly_f', 'vdp_n_f', 'vdp_pstop_f', 'vdp_poly_r', 'vdp_n_r', 'vdp_pstop_r']

    rcParams.update({'font.size': 12})

    fig = plt.figure(figsize=(8, 6))

    gs = gridspec.GridSpec(int(len(keys) / 2), 2)
    labels = [b.short_batch(vpx=False) for b in pqc_batches]

    for i, key in enumerate(keys):
        ax = plt.subplot(gs[i])
        plt.title(pqc_batches[0].dataseries[key].label, fontsize=13)
        plt.grid(axis='y', linestyle=':')
        ax.set_ylabel(pqc_batches[0].dataseries[key].unit)

        data = [b.dataseries[key].get_stats().values for b in pqc_batches]

        ax.boxplot(data, labels=labels)
        if i < (len(keys) - 2):
            ax.set_xticklabels([])

    fig.tight_layout(h_pad=1.0)
    plot_save(fig, filename)
    plt.close()


def plot_vdp_boxplot(pqc_batches: list, filename: str) -> None:
    """Render VdP boxplot figure for dataseries."""
    print(f"rendering boxplot {filename} ...")
    rcParams.update({'font.size': 14})

    fig = plt.figure(figsize=(8, 6))

    gs = gridspec.GridSpec(3, 1)
    labels = [b.short_batch() for b in pqc_batches]

    ax = plt.subplot(gs[0])
    plt.grid(axis='y', linestyle=':')
    plt.title(pqc_batches[0].vdp_poly_tot().label, fontsize=15)
    ax.set_ylabel(pqc_batches[0].vdp_poly_tot().unit)
    data = [b.vdp_poly_tot().get_stats().values for b in pqc_batches]
    ax.boxplot(data, labels=labels)

    ax = plt.subplot(gs[1])
    plt.grid(axis='y', linestyle=':')
    plt.title(pqc_batches[0].vdp_n_tot().label, fontsize=15)
    ax.set_ylabel(pqc_batches[0].vdp_n_tot().unit)
    data = [b.vdp_n_tot().get_stats().values for b in pqc_batches]
    ax.boxplot(data, labels=labels)

    ax = plt.subplot(gs[2])
    plt.grid(axis='y', linestyle=':')
    plt.title(pqc_batches[0].vdp_pstop_tot().label, fontsize=15)
    ax.set_ylabel(pqc_batches[0].vdp_pstop_tot().unit)
    data = [b.vdp_pstop_tot().get_stats().values for b in pqc_batches]
    ax.boxplot(data, labels=labels)

    fig.tight_layout(h_pad=1.0)
    plot_save(fig, filename)
    plt.close()


def plot_save(fig: Figure, filename: str) -> None:
    """Save plot to file, creates path if not exist."""
    path = os.path.dirname(filename)
    create_dir(path)
    fig.savefig(filename)


def load_configuration(name: str) -> dict:
    """Load dataseries configuration from YAML file in directory `config`."""
    filename = os.path.join(os.path.dirname(__file__), 'config', f'{name}.yaml')
    if not os.path.isfile(filename):
        raise ValueError(f"No such configuration: {name}")
    with open(filename) as fp:
        return yaml.safe_load(fp)


def apply_configuration(dataseries: dict, config: dict = None) -> None:
    if config is None:
        config = {}
    for key, d in config.get('pqc_values', {}).items():
        series = dataseries.get(key)
        if series is not None:
            for name, value in d.items():
                types = {
                    'name': str,
                    'label': str,
                    'expected_value': float,
                    'unit': str,
                    'show_multiplier': float,
                    'stray': float,
                    'min_allowed': float,
                    'max_allowed': float
                }
                if name in types:
                    value = types.get(name)(value)
                    setattr(series, name, value)


def load_batch(path: str, outdir: str = None, *, lazy: bool = False,
               create_plots: bool = False, create_histograms: bool = False,
               force_eval: bool = False, config: dict = None) -> PQC_resultset:
    """Create PQC resultset from batch directory and optionally creates plots
    and histograms.
    """
    has_outdir = outdir is not None
    create_plots = create_plots and has_outdir
    create_histograms = create_histograms and has_outdir
    batchname = os.path.basename(os.path.normpath(path))
    print(f"Batch: {batchname}")
    pqc_results = PQC_resultset(batchname)

    # Apply configuration
    apply_configuration(pqc_results.dataseries, config)

    if lazy and has_outdir:
        try:
            analysis_time = os.path.getmtime(pqc_results.analysis_dir(outdir))
            measure_time = os.path.getmtime(path)

            if analysis_time > measure_time:
                print("lazy mode: nothing to do")
                sys.exit(0)  # TODO!
        except FileNotFoundError:
            print("lazy but first time")

    # TODO
    # Prevent to create ouput directory if not given
    if has_outdir:
        pqc_results.prepare_analysis_dir(outdir)
    pqc_results.analyze(path, create_plots=create_plots, force_eval=force_eval)

    # Render histograms (optional)
    if create_histograms:
        print("rendering histograms... ", end="", flush=True)
        pqc_results.create_histograms()
        print("done.")

    return pqc_results


def run_multibatch(path: str, outdir: str, *, config: dict):
    print("Multibatch mode - experimental!")
    dirs = glob.glob(os.path.join(path, "*"))
    dirs = [t for t in dirs if "histograms" not in t and "VPX" in t]  # TODO!
    dirs = [t for t in dirs if os.path.isdir(t)]

    pqc_batches = []
    pqc_slices = []

    for diri in dirs:
        print(f"Current dir: {diri}")
        res = load_batch(diri, config=config)
        res.sort_by_time()
        pqc_batches.append(res)
        pqc_slices.extend(res.split(4))
    #    #res.create_histograms(args.path)

    print(f"loaded {len(pqc_batches)} batches")
    plot_timeline(pqc_batches, os.path.join(outdir, "histograms", "timeline_batch.png"))
    plot_timeline(pqc_slices, os.path.join(outdir, "histograms", "timeline_fine.png"), show_batch_numbers=False)

    plot_vdp_boxplot(pqc_batches, os.path.join(outdir, "histograms", "boxplot_vdp.png"))
    plot_boxplot(pqc_batches, os.path.join(outdir, "histograms", "boxplot_a.png"), keys=['t_line_n', 'r_contact_n', 't_line_pstop2', 't_line_pstop4', 'r_contact_poly', 'v_th'])
    plot_boxplot(pqc_batches, os.path.join(outdir, "histograms", "boxplot_b.png"), keys=['vdp_p_cross_bridge_f', 'vdp_p_cross_bridge_r', 't_line_p_cross_bridge', 'v_bd', 'i600', 'v_fd'])
    plot_boxplot(pqc_batches, os.path.join(outdir, "histograms", "boxplot_c.png"), keys=['rho', 'conc', 't_ox', 'n_ox', 'c_acc_m', 'i_surf'])


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    parser.add_argument('-m', dest='multibatch', action='store_true', help='multibatch mode, e. g. for time analysis (experimental)')
    parser.add_argument('-o', dest='outdir', metavar='DIR', help='override output directory location')
    parser.add_argument('-l', dest='lazy', action='store_true', help='lazy evaluation: skip if the measurement folder is older than analysis folder')
    parser.add_argument('-H', dest='histograms', action='store_true', help='create histograms')
    parser.add_argument('-P', dest='plots', action='store_true', help='create plots (for each single measurement used)')
    parser.add_argument('-f', dest='force', action='store_true', help='force evaluating all directories (normally, only directories with at least one VdP measurement are evaluated to prevent blank lines if the is a wrong file or so)')
    parser.add_argument('-t', dest='templates', metavar='EXPR', action='append', default=[], help='select templates to render (eg. -t*.tex -t*.html -tall.txt)')
    parser.add_argument('-c', '--config', metavar='NAME', default='default', help='select custom configuration')
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Output directory, input path if not set
    outdir = args.outdir or args.path

    # Load configuration
    config = load_configuration(args.config)

    # Create output directory
    create_dir(outdir)

    if args.multibatch:
        run_multibatch(args.path, outdir, config=config)
    else:
        pqc_results = load_batch(
            args.path,
            outdir,
            lazy=args.lazy,
            create_plots=args.plots,
            create_histograms=args.histograms,
            force_eval=args.force,
            config=config
        )
        render_templates(pqc_results, args.templates)
    #pdb.set_trace()
    plt.show()


if __name__ == "__main__":
    main()
