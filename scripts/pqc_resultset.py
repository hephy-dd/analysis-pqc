#!/usr/bin/env python3

import glob
import os
import re
import pdb

import numpy as np

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec

from analysis_pqc import params

import pqc_analysis_json as pqc
from pqc_values import PQC_Values, make_chunks
from pqc_analysis_json import AnalysisOptions

__all__ = ["Histogram", "PQC_resultset"]


class Histogram:
    def __init__(self, filename, *, title=None, description=None):
        self.filename = filename
        self.title = title or ""
        self.description = description or ""


class PQC_resultset:

    OUTPUT_PREFIX = "analysis_"

    def __init__(self, batchname, dataseries=None):
        self.batch = batchname
        self.labels = []
        self.flutes = []
        self.timestamps = []
        self.basepath = ""
        self.output_dir = None  # TODO
        self.plot_dir = None  # TODO
        self.histogram_dir = None  # TODO
        self.histograms = []
        self.rawdata = {}

        if dataseries is None:
            self.dataseries = {
                "xtimestamps": self.timestamps,
                "xlabels": self.labels,
                "xflutes": self.flutes,
            }

            # =================================================== Flute 1 ===================================================

            self.dataseries["v_th"] = PQC_Values("fet", "FET Vth", 4.0, "V", stray=0.25)

            self.dataseries["v_fb2"] = PQC_Values(
                "v_fb", "Flatband voltage", 2.5, "V", stray=0.33
            )
            self.dataseries["t_ox"] = PQC_Values(
                "t_ox", "Oxide thickness", 670.0, "nm", 1e3, stray=0.33
            )
            self.dataseries["n_ox"] = PQC_Values(
                "n_ox", "Oxide concentration", 10.5, "1E10cm^-2", 1e-10
            )
            self.dataseries["c_acc_m"] = PQC_Values(
                "c_acc", "Accumulation capacitance", 85.0, "pF", 1e12, stray=0.2
            )

            self.dataseries["cap_l"] = PQC_Values(
                "cap_l", "Capacitor", 3.0, "pF", 1e12, stray=0.5
            )
            self.dataseries["cap_r"] = PQC_Values(
                "cap_r", "Capacitor", 3.0, "pF", 1e12, stray=0.5
            )
            self.dataseries["cap_l_tox"] = PQC_Values(
                "capl_tox", "Capacitor: Oxide Thickness", 1.0, "nm", 1e9, stray=0.5
            )
            self.dataseries["cap_r_tox"] = PQC_Values(
                "capr_tox", "Capacitor: Oxide Thickness", 1.0, "nm", 1e9, stray=0.5
            )

            self.dataseries["vdp_poly_f"] = PQC_Values(
                "vdpPoly", "Polysilicon VdP", 2.2, "kOhm/sq", 1e-3, stray=0.3
            )
            self.dataseries["vdp_poly_r"] = PQC_Values(
                "vdpPoly_r", "Polysilicon VdP reverse", 2.2, "kOhm/sq", 1e-3, stray=0.3
            )

            self.dataseries["vdp_n_f"] = PQC_Values(
                "vdpN", "N+ VdP", 35.0, "Ohm/sq", stray=0.3
            )
            self.dataseries["vdp_n_r"] = PQC_Values(
                "vdpN_r", "N+ VdP reverse", 35.0, "Ohm/sq", stray=0.3
            )

            self.dataseries["vdp_pstop_f"] = PQC_Values(
                "vdpPstp", "P-stop VdP", 19.0, "kOhm/sq", 1e-3, stray=0.3
            )
            self.dataseries["vdp_pstop_r"] = PQC_Values(
                "vdpPstp_r", "P-stop VdP reverse", 19.0, "kOhm/sq", 1e-3, stray=0.3
            )

            # =================================================== Flute 2 ===================================================

            self.dataseries["i_surf"] = PQC_Values(
                "i_surf", "Surface current", 8.0, "pA", -1e12, stray=1
            )
            self.dataseries["s0"] = PQC_Values(
                "s0", "Surface generation velocity", 1.0, "cm/s", stray=1
            )

            self.dataseries["meander_poly"] = PQC_Values(
                "meand_poly", "Polisilicon Resistor", 1.7, "MOhm", 1e-6, stray=0.5
            )

            self.dataseries["t_line_n"] = PQC_Values("lw_n", "Linewidth N+", 35.0, "um")
            self.dataseries["t_line_pstop2"] = PQC_Values(
                "lw_pstp2", "Linewidth P-stop 2 Wire", 38.0, "um"
            )
            self.dataseries["t_line_pstop4"] = PQC_Values(
                "lw_pstp4", "Linewidth P-stop 4 Wire", 55.0, "um"
            )

            self.dataseries["v_bd"] = PQC_Values(
                "v_bd", "Breakdown Voltage", 215.0, "V"
            )

            # =================================================== Flute 3 ===================================================

            self.dataseries["i600"] = PQC_Values(
                "i600", "I @ 600V", 100.0, "uA", 1e6, stray=1.0
            )

            self.dataseries["i300"] = PQC_Values(
                "i300", "I @ 300V", 100.0, "nA", 1e9, stray=1.0
            )

            self.dataseries["v_fd"] = PQC_Values(
                "v_fd", "Full depletion Voltage", 260.0, "V", stray=0.33
            )
            self.dataseries["rho"] = PQC_Values(
                "rho", "rho", 3.5, "kOhm cm", 0.1, min_allowed=3.5, max_allowed=10.0
            )
            self.dataseries["conc"] = PQC_Values(
                "d_conc", "Doping Concentration", 3.5, "1E12cm^-3", 1e-18
            )

            self.dataseries["meander_metal"] = PQC_Values(
                "meand_metal", "Metal Meander", 260.0, "Ohm", 1.0, stray=0.5
            )

            self.dataseries["vdp_metclo_f"] = PQC_Values(
                "vdp_met", "Metal Cloverleaf VdP", 25.0, "mOhm/sq", 1e3
            )
            self.dataseries["vdp_metclo_r"] = PQC_Values(
                "vdp_met_r", "Metal Cloverleaf VdP reverse", 25.0, "mOhm/sq", 1e3
            )

            self.dataseries["vdp_p_cross_bridge_f"] = PQC_Values(
                "vdp_cb", "Cross Bridge VdP", 1.5, "kOhm/sq", 1e-3
            )
            self.dataseries["vdp_p_cross_bridge_r"] = PQC_Values(
                "vdp_cb_r", "Cross Bridge VdP reverse", 1.5, "kOhm/sq", 1e-3
            )
            self.dataseries["t_line_p_cross_bridge"] = PQC_Values(
                "t_line_cb", "Linewidth cross bridge P", 35.0, "um"
            )

            self.dataseries["vdp_bulk_f"] = PQC_Values(
                "vdpBulk", "Bulk VdP Cross", 66.0, "kOhm/sq", 1e-3, stray=0.8
            )
            self.dataseries["vdp_bulk_r"] = PQC_Values(
                "vdpBulk_r", "Bulk VdP Cross rev", 66.0, "kOhm/sq", 1e-3, stray=0.8
            )
            self.dataseries["vdp_bulk_rho_old"] = PQC_Values(
                "vdpBulk_rho_old",
                "Bulk VdP Resistivity (1.218)",
                3.5,
                "kOhm*cm (1.218)",
                1.0,
                min_allowed=3.5,
                max_allowed=10.0,
            )
            self.dataseries["vdp_bulk_rho"] = PQC_Values(
                "vdpBulk_rho",
                "Bulk VdP Resistivity (1.081)",
                3.5,
                "kOhm*cm (1.081)",
                1.0,
                min_allowed=3.5,
                max_allowed=10.0,
            )

            # =================================================== Flute 4 ===================================================

            self.dataseries["i_surf05"] = PQC_Values(
                "i_surf05",
                "Surface current 05",
                11.0,
                "pA",
                -1e12,
                stray=1,
                max_allowed=25.0,
            )
            self.dataseries["s0_gcd05"] = PQC_Values(
                "s0_gcd05", "Surface generation velocity GCD05", 1.0, "cm/s", stray=1
            )
            self.dataseries["i_bulk05"] = PQC_Values(
                "i_bulk05", "Bulk current 05", 0.7, "pA", -1e12, stray=1
            )

            self.dataseries["r_contact_n"] = PQC_Values(
                "r_cont_n", "Rcontact N+", 27.0, "Ohm"
            )
            self.dataseries["r_contact_poly"] = PQC_Values(
                "r_cont_poly", "Rcontact polysilicon", 100.0, "kOhm", 1e-3
            )

            self.dataseries["contact_poly"] = PQC_Values(
                "cont_poly", "Contact Chain PolySi", 20.0, "MOhm", 1e-6, stray=0.5
            )
            self.dataseries["contact_p"] = PQC_Values(
                "cont_p", "Contact Chain P", 85.0, "kOhm", 1e-3, stray=0.5
            )
            self.dataseries["contact_n"] = PQC_Values(
                "cont_n", "Contact Chain N", 85.0, "kOhm", 1e-3, stray=0.5
            )
        else:
            self.dataseries = dataseries

    def vdp_poly_tot(self):
        return PQC_Values.merge(
            [self.dataseries["vdp_poly_f"], self.dataseries["vdp_poly_r"]],
            "vdp_poly_tot",
            "PolySi VdP both",
        )

    def vdp_n_tot(self):
        return PQC_Values.merge(
            [self.dataseries["vdp_n_f"], self.dataseries["vdp_n_r"]],
            "vdp_N_tot",
            "N+ VdP both",
        )

    def vdp_pstop_tot(self):
        return PQC_Values.merge(
            [self.dataseries["vdp_pstop_f"], self.dataseries["vdp_pstop_r"]],
            "vdp_pstop_tot",
            "P-stop VdP both",
        )

    def sort_by_time(self):
        # TODO remove placeholders
        order = np.argsort([t for t in self.timestamps if t])
        # print(str(order))

        for key in self.dataseries:
            if type(self.dataseries[key]) is PQC_Values:
                self.dataseries[key].rearrange(order)
            else:
                self.dataseries[key][:] = [
                    self.dataseries[key][i] for i in order
                ]  # we want to keep the original object so references are preserved

    # warning: this creates not full copies, only the dict is available then
    def split(self, junk_size):
        # print(str((self.dataseries["vdp_n_f"]).split(junk_size)))

        # ret = [PQC_resultset( self.value=i) for i in make_chunks(self.value, junk_size)]
        ret = [None] * len(self.dataseries["vdp_n_f"].split(junk_size))

        for i in range(0, len(ret)):
            ret[i] = PQC_resultset(self.batch)

        for key in self.dataseries:
            if type(self.dataseries[key]) is PQC_Values:
                ch = self.dataseries[key].split(junk_size)
            else:
                ch = [i for i in make_chunks(self.dataseries[key], junk_size)]

            for i in range(0, len(ch)):
                ret[i].dataseries[key] = ch[i]

        for i in ret:
            i.timestamps = i.dataseries["xtimestamps"]
            i.flutes = i.dataseries["xflutes"]
            i.labels = i.dataseries["xlabels"]

        return ret

    def analyze_sample(self, path, create_plots=False, force_eval=False, config=None):
        """Analyze sample data and append results to dataseries."""

        # TODO
        # this finds out if there is an empty directory, assuming that there is at least one vdp measurement
        if (
            len(pqc.find_all_files_from_path(path, "van_der_pauw", whitelist=["cross"]))
            < 1
            and not force_eval
        ):
            print(f"   skipping incomplete measurement: {path}")
            return

        self.flutes.append("")
        label = os.path.basename(path)
        self.labels.append(label)
        self.rawdata[label] = {}

        # TODO
        x = pqc.find_all_files_from_path(path, "van_der_pauw")
        if len(x) > 0:
            self.timestamps.append(pqc.get_timestamp(x[-1]))
        else:
            self.timestamps.append(0)

        cross = "cross"
        plotImgLabel = self.dataseries["xlabels"][-1]
        if create_plots:
            plot_dir = self.plot_dir
        else:
            plot_dir = None

        options = AnalysisOptions(plot_dir, plotImgLabel)

        # =================================================== Flute 1 ===================================================

        v_th, rawdata = pqc.analyse_fet_data(
            pqc.find_most_recent_file(path, "fet", whitelist=[]),
            options=options,
            config=config,
        )
        self.dataseries["v_th"].append(v_th)
        if rawdata is not None:
            self.rawdata[label]["FET"] = rawdata

        _, v_fb2, t_ox, n_ox, c_acc_m, rawdata = pqc.analyse_mos_data(
            pqc.find_most_recent_file(path, "mos", whitelist=[]),
            options=options,
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["MOS_capacitor_HV_Source"] = rawdata
        self.dataseries["v_fb2"].append(v_fb2)
        self.dataseries["t_ox"].append(t_ox)
        self.dataseries["n_ox"].append(n_ox)
        self.dataseries["c_acc_m"].append(c_acc_m)

        c_mean, c_median, d, rawdata = pqc.analyse_capacitor_data(
            pqc.find_most_recent_file(
                path,
                "capacitor",
                whitelist=["Left", "250mV", "10kHz"],
                blacklist=["mos"],
            ),
            options=options,
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label][
                "Capacitor_test_structure_Left_10kHz_250mV_HV_Source"
            ] = rawdata
        self.dataseries["cap_l"].append(c_median)
        self.dataseries["cap_l_tox"].append(d)

        c_mean, c_median, d, rawdata = pqc.analyse_capacitor_data(
            pqc.find_most_recent_file(
                path,
                "capacitor",
                whitelist=["Right", "250mV", "10kHz"],
                blacklist=["mos"],
            ),
            options=options,
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label][
                "Capacitor_test_structure_Right_10kHz_250mV_HV_Source"
            ] = rawdata
        self.dataseries["cap_r"].append(c_median)
        self.dataseries["cap_r_tox"].append(d)

        vdp_poly_f, rawdata = pqc.analyse_van_der_pauw_data(
            pqc.find_most_recent_file(
                path,
                "van_der_pauw",
                whitelist=["Polysilicon", cross],
                blacklist=["reverse"],
            ),
            options=options.pushPrefix("VdP_poly_fwd"),
            config=config,
        )
        self.dataseries["vdp_poly_f"].append(vdp_poly_f)
        if rawdata is not None:
            self.rawdata[label]["Polysilicon_Van-der-Pauw_cross"] = rawdata

        vdp_poly_r, rawdata = pqc.analyse_van_der_pauw_data(
            pqc.find_most_recent_file(
                path, "van_der_pauw", whitelist=["Polysilicon", "reverse", cross]
            ),
            options=options.pushPrefix("VdP_poly_rev"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["Reverse_Polysilicon_Van-der-Pauw_cross"] = rawdata
        self.dataseries["vdp_poly_r"].append(vdp_poly_r)

        vdp_n_f, rawdata = pqc.analyse_van_der_pauw_data(
            pqc.find_most_recent_file(
                path, "van_der_pauw", whitelist=["n", cross], blacklist=["reverse"]
            ),
            options=options.pushPrefix("VdP_N_fwd"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["N_Van-der-Pauw_cross"] = rawdata
        self.dataseries["vdp_n_f"].append(vdp_n_f)

        vdp_n_r, rawdata = pqc.analyse_van_der_pauw_data(
            pqc.find_most_recent_file(
                path, "van_der_pauw", whitelist=["n", "reverse", cross]
            ),
            options=options.pushPrefix("VdP_N_rev"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["Reverse_N_Van-der-Pauw_cross"] = rawdata
        self.dataseries["vdp_n_r"].append(vdp_n_r)

        vdp_pstop_f, rawdata = pqc.analyse_van_der_pauw_data(
            pqc.find_most_recent_file(
                path, "van_der_pauw", whitelist=["P_stop", cross], blacklist=["reverse"]
            ),
            options=options.pushPrefix("VdP_P-Stop_fwd"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["P-stop_Van-der-Pauw_cross"] = rawdata
        self.dataseries["vdp_pstop_f"].append(vdp_pstop_f)

        vdp_pstop_r, rawdata = pqc.analyse_van_der_pauw_data(
            pqc.find_most_recent_file(
                path, "van_der_pauw", whitelist=["P_stop", "reverse", cross]
            ),
            options=options.pushPrefix("VdP_P-Stop_rev"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["Reverse_P-stop_Van-der-Pauw_cross"] = rawdata
        self.dataseries["vdp_pstop_r"].append(vdp_pstop_r)

        # =================================================== Flute 2 ===================================================

        i_surf, _, rawdata = pqc.analyse_gcd_data(
            pqc.find_most_recent_file(path, "gcd", whitelist=[]),
            options=options.pushPrefix("GCD"),
            config=config,
        )  # only i_surf valid
        self.dataseries["i_surf"].append(i_surf)
        # s0 = i_surf / q / ni[cm^-3] / Agate[cm^2]
        s0 = -i_surf / 1.602e-19 / 7.01e9 / 0.505e-2
        self.dataseries["s0"].append(s0)
        if rawdata is not None:
            rawdata.add_data({"s0": s0})
            self.rawdata[label]["GCD"] = rawdata

        t_line_n, rawdata = pqc.analyse_linewidth_data(
            pqc.find_most_recent_file(path, "linewidth", whitelist=["n"]),
            r_sheet=self.dataseries["vdp_n_f"].values[-1],
            options=options.pushPrefix("lw_n"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["N_linewidth_structure"] = rawdata
        self.dataseries["t_line_n"].append(t_line_n)

        t_line_pstop2, rawdata = pqc.analyse_linewidth_data(
            pqc.find_most_recent_file(
                path, "linewidth", whitelist=["P_stop", "2_wire"]
            ),
            r_sheet=self.dataseries["vdp_pstop_f"].values[-1],
            options=options.pushPrefix("lw_p2"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["P-stop_linewidth_structure_2-wire"] = rawdata
        self.dataseries["t_line_pstop2"].append(t_line_pstop2)

        t_line_pstop4, rawdata = pqc.analyse_linewidth_data(
            pqc.find_most_recent_file(
                path, "linewidth", whitelist=["P_stop", "4_wire"]
            ),
            r_sheet=self.dataseries["vdp_pstop_f"].values[-1],
            options=options.pushPrefix("lw_p4"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["P-stop_linewidth_structure_4-wire"] = rawdata
        self.dataseries["t_line_pstop4"].append(t_line_pstop4)

        meander_poly, rawdata = pqc.analyse_meander_data(
            pqc.find_most_recent_file(path, "meander", whitelist=["polysilicon"]),
            options=options.pushPrefix("meander_poly"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["Polysilicon_meander"] = rawdata
        self.dataseries["meander_poly"].append(meander_poly)

        v_bd, rawdata = pqc.analyse_breakdown_data(
            pqc.find_most_recent_file(path, "breakdown", whitelist=[]),
            options=options,
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["Dielectric_Breakdown_1"] = rawdata
        self.dataseries["v_bd"].append(v_bd)

        # =================================================== Flute 3 ===================================================

        # we want this for FLute_3 and not Flute_1
        i600, i300, rawdata = pqc.analyse_iv_data(
            pqc.find_most_recent_file(path, "iv", whitelist=["3"]),
            options=options.pushPrefix("IV_DiodeHalf"),
            config=config,
        )
        self.dataseries["i600"].append(i600)
        self.dataseries["i300"].append(i300)
        if rawdata is not None:
            self.rawdata[label]["Diode_IV"] = rawdata
        v_fd, rho, conc, rawdata = pqc.analyse_cv_data(
            pqc.find_most_recent_file(path, "cv", whitelist=["3"]),
            options=options.pushPrefix("CV_DiodeHalf"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["Diode_CV"] = rawdata
        self.dataseries["v_fd"].append(v_fd)
        self.dataseries["rho"].append(rho)
        self.dataseries["conc"].append(conc)

        vdp_metclo_f, rawdata = pqc.analyse_van_der_pauw_data(
            pqc.find_most_recent_file(
                path,
                "van_der_pauw",
                whitelist=["metal", "clover"],
                blacklist=["reverse"],
            ),
            options=options.pushPrefix("VdP_Metal_fwd"),
            min_correlation=0.95,
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["Metal_clover_leaf_Van-der-Pauw"] = rawdata
        self.dataseries["vdp_metclo_f"].append(vdp_metclo_f)

        vdp_metclo_r, rawdata = pqc.analyse_van_der_pauw_data(
            pqc.find_most_recent_file(
                path,
                "van_der_pauw",
                whitelist=["metal", "clover", "reverse"],
                blacklist=[],
            ),
            options=options.pushPrefix("VdP_Metal_rev"),
            min_correlation=0.95,
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["Reverse_Metal_clover_leaf_Van-der-Pauw"] = rawdata
        self.dataseries["vdp_metclo_r"].append(vdp_metclo_r)

        vdp_p_cross_bridge_f, rawdata = pqc.analyse_van_der_pauw_data(
            pqc.find_most_recent_file(
                path,
                "van_der_pauw",
                whitelist=["P", "cross_bridge"],
                blacklist=["reverse"],
            ),
            options=options.pushPrefix("VdP_P-edge_fwd"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["P_cross-bridge_Van-der-Pauw"] = rawdata
        self.dataseries["vdp_p_cross_bridge_f"].append(vdp_p_cross_bridge_f)

        vdp_p_cross_bridge_r, rawdata = pqc.analyse_van_der_pauw_data(
            pqc.find_most_recent_file(
                path, "van_der_pauw", whitelist=["P", "cross_bridge", "reverse"]
            ),
            options=options.pushPrefix("VdP_P-edge_rev"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["Reverse_P_cross-bridge_Van-der-Pauw"] = rawdata
        self.dataseries["vdp_p_cross_bridge_r"].append(vdp_p_cross_bridge_r)

        t_line_p_cross_bridge, rawdata = pqc.analyse_linewidth_data(
            pqc.find_most_recent_file(
                path, "linewidth", whitelist=["P", "cross_bridge"]
            ),
            r_sheet=self.dataseries["vdp_p_cross_bridge_f"].values[-1],
            options=options.pushPrefix("lw_P-edge"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["P_cross-bridge_linewidth"] = rawdata
        self.dataseries["t_line_p_cross_bridge"].append(t_line_p_cross_bridge)

        vdp_bulk_f, rawdata = pqc.analyse_van_der_pauw_data(
            pqc.find_most_recent_file(
                path, whitelist=["bulk", "cross"], blacklist=["reverse"]
            ),
            options=options.pushPrefix("VdP_bulk_fwd"),
            min_correlation=0.85,
            config=config,
        )
        self.dataseries["vdp_bulk_f"].append(vdp_bulk_f)
        if rawdata is not None:
            self.rawdata[label]["Bulk_cross"] = rawdata

        vdp_bulk_r, rawdata = pqc.analyse_van_der_pauw_data(
            pqc.find_most_recent_file(path, whitelist=["bulk", "reverse", "cross"]),
            options=options.pushPrefix("VdP_bulk_rev"),
            min_correlation=0.85,
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["Reverse_bulk_cross"] = rawdata
        self.dataseries["vdp_bulk_r"].append(vdp_bulk_r)

        bulk_rsheet = (
            self.dataseries["vdp_bulk_f"].values[-1]
            + self.dataseries["vdp_bulk_r"].values[-1]
        ) / 2
        bulk_rho0 = bulk_rsheet * np.log(2.0) * 2.0 * 187e-6 / (2 - np.sqrt(2.0)) * 1e-1
        self.dataseries["vdp_bulk_rho"].append(bulk_rho0 * 1.081)
        if "Bulk_cross" in self.rawdata[label]:
            self.rawdata[label]["Bulk_cross"].add_data(
                {"vdp_bulk_rho": bulk_rho0 * 1.081}
            )
        if "Reverse_bulk_cross" in self.rawdata[label]:
            self.rawdata[label]["Reverse_bulk_cross"].add_data(
                {"vdp_bulk_rho": bulk_rho0 * 1.081}
            )
        self.dataseries["vdp_bulk_rho_old"].append(bulk_rho0 * 1.218)

        meander_metal, rawdata = pqc.analyse_meander_data(
            pqc.find_most_recent_file(path, "meander", whitelist=["metal"]),
            options=options.pushPrefix("meander_metal"),
            config=config,
        )
        self.dataseries["meander_metal"].append(meander_metal)
        if rawdata is not None:
            self.rawdata[label]["Metal_meander"] = rawdata

        # =================================================== Flute 4 ===================================================

        i_surf05, i_bulk05, rawdata = pqc.analyse_gcd_data(
            pqc.find_most_recent_file(path, "gcd05", whitelist=[]),
            options=options.pushPrefix("GCD05"),
            config=config,
        )  # for i_bulk
        self.dataseries["i_surf05"].append(i_surf05)
        self.dataseries["i_bulk05"].append(i_bulk05)
        s0_gcd05 = -i_surf05 / 1.602e-19 / 7.01e9 / 0.732e-2
        self.dataseries["s0_gcd05"].append(s0_gcd05)
        if rawdata is not None:
            rawdata.add_data({"s0": s0_gcd05})
            self.rawdata[label]["GCD05"] = rawdata

        r_contact_n, rawdata = pqc.analyse_cbkr_data(
            pqc.find_most_recent_file(path, "cbkr", whitelist=["n"]),
            r_sheet=self.dataseries["vdp_n_f"].values[-1],
            options=options.pushPrefix("cbkr_n"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["N_CBKR"] = rawdata
        self.dataseries["r_contact_n"].append(r_contact_n)

        r_contact_poly, rawdata = pqc.analyse_cbkr_data(
            pqc.find_most_recent_file(path, "cbkr", whitelist=["Polysilicon"]),
            r_sheet=self.dataseries["vdp_poly_f"].values[-1],
            options=options.pushPrefix("cbkr_poly"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["Polysilicon_CBKR"] = rawdata
        self.dataseries["r_contact_poly"].append(r_contact_poly)

        contact_poly, rawdata = pqc.analyse_contact_data(
            pqc.find_most_recent_file(
                path, "contact", whitelist=["chain", "polysilicon"]
            ),
            options=options.pushPrefix("contact_chain_poly"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["Polysilicon_contact_chain"] = rawdata
        self.dataseries["contact_poly"].append(contact_poly)

        contact_p, rawdata = pqc.analyse_contact_data(
            pqc.find_most_recent_file(path, "contact", whitelist=["chain", "P"]),
            options=options.pushPrefix("contact_chain_p"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["P_contact_chain"] = rawdata
        self.dataseries["contact_p"].append(contact_p)

        contact_n, rawdata = pqc.analyse_contact_data(
            pqc.find_most_recent_file(path, "contact", whitelist=["chain", "N"]),
            options=options.pushPrefix("contact_chain_n"),
            config=config,
        )
        if rawdata is not None:
            self.rawdata[label]["N_contact_chain"] = rawdata
        self.dataseries["contact_n"].append(contact_n)

    def analyze(self, basepath, create_plots=False, force_eval=False, config=None):
        """Analyze and collect results of a batch of samples inside a directory."""

        def is_dataset(path):
            # Keep only directories
            if not os.path.isdir(path):
                return False
            basename = os.path.basename(path)
            # Skip existing analysis directories
            if basename.startswith(self.OUTPUT_PREFIX):
                return False
            return True

        self.basepath = basepath
        dirs = glob.glob(os.path.join(basepath, "*"))
        dirs = [t for t in dirs if is_dataset(t)]

        dirs.sort()

        for i, sample_path in enumerate(dirs):
            # print the currently analyzing directory
            precent = i * 100.0 / len(dirs)
            sample_name = os.path.basename(os.path.normpath(sample_path))
            print("[{:2.0f}%] Analyzing: {:}".format(precent, sample_name))
            self.analyze_sample(
                sample_path,
                create_plots=create_plots,
                force_eval=force_eval,
                config=config,
            )

        self.dataseries["xlabels"].extend(PQC_Values.get_stats_labels())
        self.dataseries["xflutes"].extend([""] * len(PQC_Values.get_stats_labels()))
        self.dataseries["xtimestamps"].extend([""] * len(PQC_Values.get_stats_labels()))

    def statusbar(
        self, pqc_value_statistics, axes, single=True, start=-0.5, stop=0.5, label=""
    ):
        relOK = len(pqc_value_statistics.values) / pqc_value_statistics.nTot
        relNan = pqc_value_statistics.nNan / pqc_value_statistics.nTot
        relTooHigh = pqc_value_statistics.nTooHigh / pqc_value_statistics.nTot
        relTooLow = pqc_value_statistics.nTooLow / pqc_value_statistics.nTot

        average_delta = (start - stop) / 2
        center = stop + average_delta
        width = stop - start
        if single:
            alpha = 1.0
        else:
            alpha = 0.5

        p1 = axes.bar(center, (relOK,), width, color="green", alpha=alpha)
        p2 = axes.bar(
            center, (relNan,), width, bottom=(relOK,), color="red", alpha=alpha
        )
        p3 = axes.bar(
            center,
            (relTooHigh,),
            width,
            bottom=(relOK + relNan,),
            color="orange",
            alpha=alpha,
        )
        p4 = axes.bar(
            center,
            (relTooLow,),
            width,
            bottom=(relOK + relNan + relTooHigh,),
            color="yellow",
            alpha=alpha,
        )

        if single:
            minthreshhold = 0.02
        else:
            minthreshhold = 0.15

        if single or label:
            axes.text(
                center,
                relOK - minthreshhold / 2,
                "OK",
                horizontalalignment="center",
                verticalalignment="top",
            )
            if relNan > minthreshhold:
                axes.text(
                    center,
                    relOK + relNan - minthreshhold / 2,
                    "Failed",
                    horizontalalignment="center",
                    verticalalignment="top",
                )
            if relTooHigh > minthreshhold:
                axes.text(
                    center,
                    relOK + relNan + relTooHigh - minthreshhold / 2,
                    "High",
                    horizontalalignment="center",
                    verticalalignment="top",
                )
            if relTooLow > minthreshhold:
                axes.text(
                    center,
                    relOK + relNan + relTooHigh + relTooLow - minthreshhold / 2,
                    "Low",
                    horizontalalignment="center",
                    verticalalignment="top",
                )

        axes.text(
            center,
            0.01,
            label,
            horizontalalignment="center",
            verticalalignment="bottom",
        )

        plt.yticks([])
        plt.ylim([0, 1])
        if single:
            plt.xticks([])
            plt.xlim([start, stop])

    def histogram(self, pqc_values, path, stray=1.4, range_extension=None):
        # Plot stats
        if range_extension is not None:
            stats = pqc_values.get_stats(
                min_allowed=0, max_allowed=pqc_values.max_allowed * range_extension
            )
        else:
            stats = pqc_values.get_stats()

        if len(stats.values) == 1 and stats.nNan == 1:
            print(f"warning: skipping plot due to no valid results: {pqc_values.name}")
            return

        fig = plt.figure(figsize=(8, 6))

        gs = gridspec.GridSpec(1, 2, width_ratios=[10, 1])
        ax0 = plt.subplot(gs[0])

        # Plot title
        if range_extension is not None:
            title = f"{self.batch}: {pqc_values.label}, Ext: {range_extension:5.1E}"
            plt.ticklabel_format(axis="x", style="sci", scilimits=(-2, 2))
        else:
            title = f"{self.batch}: {pqc_values.label}"
        plt.title(title, fontsize=18)

        ax1 = plt.subplot(gs[1])

        # Plot description
        if range_extension is not None:
            # ax0.hist(pqc_values.value, bins=50, facecolor='blueviolet', alpha=1, edgecolor='black', linewidth=1)
            ax0.hist(
                stats.values,
                bins=50,
                range=[0, pqc_values.max_allowed * range_extension],
                facecolor="blueviolet",
                alpha=1,
                edgecolor="black",
                linewidth=1,
            )
            descNum = (
                "Total: {}\nShown: {}\nFailed: {}\nToo high: {}\nToo low: {}".format(
                    stats.nTot,
                    len(stats.values),
                    stats.nNan,
                    stats.nTooHigh,
                    stats.nTooLow,
                )
            )
        else:
            ax0.hist(
                stats.values,
                bins=20,
                range=[pqc_values.min_allowed, pqc_values.max_allowed],
                facecolor="blue",
                alpha=1,
                edgecolor="black",
                linewidth=1,
            )
            descNum = (
                "Total: {}\nShown: {}\nFailed: {}\nToo high: {}\nToo low: {}".format(
                    stats.nTot,
                    len(stats.values),
                    stats.nNan,
                    stats.nTooHigh,
                    stats.nTooLow,
                )
            )

        ax0.set_xlabel(pqc_values.unit)
        ax0.set_ylabel("occurences")

        # descNum = "Total number: {}\nShown: {:2.0f}%, {}\nFailed: {:2.0f}%, {}\nToo high: {:2.0f}%, {}\nToo low: {:2.0f}%, {}".format(stats.nTot, len(stats.values)/stats.nTot*1e2, len(stats.values), stats.nNan/stats.nTot*1e2, stats.nNan, stats.nTooHigh/stats.nTot*1e2, stats.nTooHigh, stats.nTooLow/stats.nTot*1e2, stats.nTooLow)

        fig.text(
            0.83,
            0.85,
            descNum,
            bbox=dict(facecolor="red", alpha=0.6),
            horizontalalignment="right",
            verticalalignment="top",
        )

        if abs(stats.selMed) < 9.99:
            descStat = f"Total median: {stats.totMed:8.2f} {pqc_values.unit}\n"
            descStat = (
                descStat
                + "Selected avg: {0:8.2f} {3}\nSel median: {2:8.2f} {3}\nSelected Std: {1:8.2f} {3}".format(
                    stats.selAvg, stats.selStd, stats.selMed, pqc_values.unit
                )
            )
        elif abs(stats.totAvg) < 1e6:
            descStat = f"Total median: {stats.totMed:8.1f} {pqc_values.unit}\n"
            descStat = (
                descStat
                + "Selected avg: {0:8.1f} {3}\nSel median: {2:8.1f} {3}\nSelected Std: {1:8.1f} {3}".format(
                    stats.selAvg, stats.selStd, stats.selMed, pqc_values.unit
                )
            )
        else:
            descStat = "Total avg: {0:9.2E} {4}\nTotal median: {1:9.2E} {4}\nSelected avg: {2:9.2E} {4}\nSel median: {3:9.2E} {4}".format(
                stats.totAvg, stats.totMed, stats.selAvg, stats.selMed, pqc_values.unit
            )
        fig.text(
            0.45,
            0.85,
            descStat,
            bbox=dict(facecolor="yellow", alpha=0.85),
            horizontalalignment="right",
            verticalalignment="top",
        )

        # for i in [(stats.totAvg, 'purple', 'solid'), (stats.totMed, 'purple', 'dashed'), (stats.selAvg, 'green', 'solid'), (stats.selMed, 'green', 'dashed')]:
        #    if (i[0] < pqc_values.max_allowed) and (i[0] > pqc_values.min_allowed):
        #        ax0.vlines(x = i[0], ymin = 0, ymax = 3,
        #            colors = i[1], linestyles=i[2],
        #            label = 'vline_multiple - full height')

        self.statusbar(stats, ax1)

        fig.tight_layout(h_pad=1.0)
        if range_extension is not None:
            filename = os.path.join(path, f"{pqc_values.name}_erhist.png")
        else:
            filename = os.path.join(path, f"{pqc_values.name}_hist.png")
        fig.savefig(filename)
        plt.close()

        # Append histogram to resultset
        self.histograms.append(
            Histogram(
                os.path.relpath(filename, self.output_dir),  # relative path for urls
                title=title,
                description=descNum,
            )
        )

    def analysis_dir(self, base_dir=None):
        return os.path.join(base_dir or "", f"{self.OUTPUT_PREFIX}{self.batch}")

    def prepare_analysis_dir(self, base_dir=None, lazy=False):
        """either creates or empties analysis folder."""
        if base_dir is None:
            base_dir = self.basepath
        self.output_dir = self.analysis_dir(base_dir)
        self.plot_dir = os.path.join(self.output_dir, "plots")
        self.histogram_dir = os.path.join(self.output_dir, "histograms")
        try:
            os.mkdir(self.output_dir)
        except OSError:
            files = glob.glob(os.path.join(self.output_dir, "*"))
            for f in files:
                if os.path.isfile(f):
                    os.remove(f)
        try:
            os.mkdir(self.plot_dir)
        except OSError:
            files = glob.glob(os.path.join(self.plot_dir, "*"))
            for f in files:
                os.remove(f)
        try:
            os.mkdir(self.histogram_dir)
        except OSError:
            files = glob.glob(os.path.join(self.histogram_dir, "*"))
            for f in files:
                os.remove(f)

    def create_histograms(self):
        matplotlib.rcParams.update({"font.size": 14})

        histogram_dir = outdir = os.path.join(self.output_dir, "histograms")

        for key in self.dataseries:
            if not key.startswith("x"):
                self.histogram(self.dataseries[key], histogram_dir)

        self.histogram(self.vdp_poly_tot(), histogram_dir)
        self.histogram(self.vdp_n_tot(), histogram_dir)
        self.histogram(self.vdp_pstop_tot(), histogram_dir)

        self.histogram(self.vdp_poly_tot(), histogram_dir, range_extension=1.5e2)
        self.histogram(self.vdp_n_tot(), histogram_dir, range_extension=1.5e2)
        self.histogram(self.vdp_pstop_tot(), histogram_dir, range_extension=1.5e2)

    def short_label(self, i):
        fl = "x"
        try:
            lbl_list = [2, 5]
            if "Left" in self.flutes[i]:
                fl = " L"
            elif "Right" in self.flutes[i]:
                fl = " R"
            elif "RL" in self.flutes[i]:
                fl = " RL"
            elif "UL" in self.flutes[i]:
                fl = " UL"
            else:
                fl = " err"
            return " ".join([self.labels[i].split("_")[j] for j in lbl_list]) + fl
        except:
            return self.labels[i] + fl

    def short_batch(self, vpx=True):
        lbl_list = [0]
        s = " ".join([self.batch.split("_")[i] for i in lbl_list])
        if vpx:
            return s
        else:
            return re.sub(r"VPX", "", s)
