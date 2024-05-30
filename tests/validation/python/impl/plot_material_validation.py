# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# detray includes
import plotting
from .plot_material_scan import get_n_bins

# python includes
import math
import numpy as np
import os
import pandas as pd


""" Read the material scan data from file and prepare data frame """
def read_material_data(inputdir, logging):

    # Input data directory
    data_dir = os.fsencode(inputdir)

    detector_name = "default_detector"
    material_scan_file = ""

    # Find the data files by naming convention
    for file in os.listdir(data_dir):
        filename = os.fsdecode(file)

        if filename.find('navigaion_material_trace_') != -1:
            material_scan_file = inputdir + "/" + filename
            file_name = os.path.basename(material_scan_file)
            detector_name = file_name.removeprefix('material_scan_ ')
            detector_name = detector_name.removesuffix('.csv')

    detector_name = detector_name.replace('_', ' ')

    df = pd.read_csv(material_scan_file)

    return detector_name, df


""" Plot the material thickness in units of X_0 vs eta """
def compare_X0(df1, df2, detector, plotFactory, var='eta', out_format =  "pdf"):

    # Histogram bin edges
    xBinning, yBinning = get_n_bins(df1)
    lgd_ops = plotting.get_legend_options()
    lgd_ops._replace(loc = 'upper center')

    normalization = len(yBinning if var == 'eta' else xBinning) - 1 

    binning = xBinning if var == 'eta' else yBinning
    hist_data = plotFactory.hist1D(
                            x      = df[var],
                            w      = df['mat_tX0'] / normalization,
                            normalize = False,
                            label  = rf'{detector}',
                            xLabel = rf'$\{var}$',
                            yLabel = r'thickness / $X_0$',
                            bins = binning,
                            showStats = False,
                            lgd_ops = lgd_ops)

    # Move the legend ouside plo
    hist_data.lgd.set_bbox_to_anchor((0.825, 1.15))

    plotFactory.write_plot(hist_data, "t_X0",  out_format)

    hist_data = plotFactory.hist1D(
                            x      = df[var],
                            w      = df['mat_sX0'] / normalization,
                            normalize = False,
                            label  = rf'{detector}',
                            xLabel = rf'$\{var}$',
                            yLabel = r'path length / $X_0$',
                            bins = binning,
                            showStats = False,
                            lgd_ops = lgd_ops)

    # Move the legend ouside plo
    hist_data.lgd.set_bbox_to_anchor((0.825, 1.15))

    plotFactory.write_plot(hist_data, f"s_X0_{var}",  out_format)


""" Plot the material thickness in units of L_0 vs eta """
def compare_L0(df, detector, plotFactory, var='eta', out_format =  "pdf"):

    # Histogram bin edges
    xBinning, yBinning = get_n_bins(df)
    lgd_ops = plotting.get_legend_options()
    lgd_ops._replace(loc = 'upper center')

    normalization = len(yBinning if var == 'eta' else xBinning) - 1

    binning = xBinning if var == 'eta' else yBinning
    hist_data = plotFactory.hist1D(
                            x      = df[var],
                            w      = df['mat_tL0'] / normalization,
                            normalize = False,
                            label  = rf'{detector}',
                            xLabel = rf'$\{var}$',
                            yLabel = r'thickness / $\Lambda_0$',
                            bins = binning,
                            showStats = False,
                            lgd_ops = lgd_ops)

    # Move the legend ouside plo
    hist_data.lgd.set_bbox_to_anchor((0.825, 1.15))

    plotFactory.write_plot(hist_data, "t_L0",  out_format)

    hist_data = plotFactory.hist1D(
                            x      = df[var],
                            w      = df['mat_sL0'] / normalization,
                            normalize = False,
                            label  = rf'{detector}',
                            xLabel = rf'$\{var}$',
                            yLabel = r'path length / $\Lambda_0$',
                            bins = binning,
                            showStats = False,
                            lgd_ops = lgd_ops)

    # Move the legend ouside plo
    hist_data.lgd.set_bbox_to_anchor((0.825, 1.15))

    plotFactory.write_plot(hist_data, f"s_L0_{var}",  out_format)
