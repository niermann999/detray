# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# detray includes
from impl import plot_material_scan as mat_plotter
from impl import read_material_scan_data
from plotting import pyplot_factory as plt_factory
from options import (common_options, detector_io_options,
                     uniform_track_generator_options, propagation_options,
                     plotting_options)
from options import (parse_common_options, parse_detector_io_options, 
                     parse_plotting_options)

# python includes
import argparse
import os
import subprocess
import sys


def __main__():

#----------------------------------------------------------------arg parsing

    descr = "Detray Material Validation"

    # Define options
    parent_parsers = [common_options(descr),
                      detector_io_options(),
                      uniform_track_generator_options(),
                      propagation_options(),
                      plotting_options()]

    parser = argparse.ArgumentParser(description=descr, parents=parent_parsers)

    parser.add_argument("--bindir", "-bin",
                        help=("Directoy containing the validation executables"),
                        default = "./bin", type=str)
    parser.add_argument("--datadir", "-data",
                        help=("Directoy containing the data files"),
                        default = "./validation_data", type=str)
    parser.add_argument("--tolerance", "-tol",
                        help=("Tolerance for comparing the material"),
                        default = 1e-04, type=float)
    parser.add_argument("--cuda",
                        help=("Run the CUDA navigation validation."),
                        action="store_true", default=False)
    parser.add_argument("--sycl",
                        help=("Run the SYCL navigation validation."),
                        action="store_true", default=False)

    args = parser.parse_args()

    logging = parse_common_options(args, descr)
    input_dir, out_dir, out_format = parse_plotting_options(args, logging)

    # IO path for data files
    datadir = args.datadir.strip("/")

    # Check bin path
    bindir = args.bindir.strip("/")
    cpu_validation = bindir + "/detray_material_validation"
    cuda_validation = bindir + "/detray_material_validation_cuda"

    if not os.path.isdir(bindir) or not os.path.isfile(cpu_validation):
        logging.error(f"Material validation binaries were not found! ({args.bindir})")
        sys.exit(1)

#------------------------------------------------------------------------run

    # Pass on the options for the validation tools
    args_list = ["--geometry_file", args.geometry_file,
                 "--material_file", args.material_file,
                 "--tol", str(args.tolerance),
                 "--phi_steps", str(args.phi_steps),
                 "--eta_steps", str(args.eta_steps),
                 "--eta_range", str(args.eta_range[0]), str(args.eta_range[1]),
                 "--min_mask_tolerance", str(args.min_mask_tol),
                 "--max_mask_tolerance", str(args.max_mask_tol),
                 "--overstep_tolerance", str(args.overstep_tol),
                 "--path_tolerance", str(args.path_tol),
                 "--rk-tolerance", str(args.rk_error_tol),
                 "--path_limit", str(args.path_limit),
                 "--search_window", str(args.search_window[0]),
                 str(args.search_window[1])]

    if args.grid_file:
        args_list = args_list + ["--grid_file", args.grid_file]

    # Run the host validation and produce the truth data
    logging.debug("Running CPU material validation")
    #subprocess.run([cpu_validation] + args_list)

    # Run the device validation (if it has been built)
    if args.cuda and os.path.isfile(cuda_validation):
        logging.debug("Running CUDA material validation")
        #subprocess.run([cuda_validation] + args_list)

    elif args.cuda:
        logging.error("Could not find CUDA material validation executable")

    if args.sycl:
        logging.error("SYCL validation is not implemented")

#------------------------------------------------------------------------plot

    logging.info("\n\nGenerating data plots...\n")

    detector_name, df = read_material_scan_data(input_dir, logging)

    plot_factory = plt_factory(out_dir + "material_", logging)

    # The histograms are not re-weighted (if the rays are not evenly distributed
    # the material in some bins might be artificially high)!
    mat_plotter.X0_vs_eta_phi(df, detector_name, plot_factory, out_format)
    mat_plotter.L0_vs_eta_phi(df, detector_name, plot_factory, out_format)
    mat_plotter.plot_X0(df, detector_name, plot_factory, 'phi', out_format)
    mat_plotter.plot_X0(df, detector_name, plot_factory, 'eta', out_format)
    mat_plotter.plot_L0(df, detector_name, plot_factory, 'phi', out_format)
    mat_plotter.plot_L0(df, detector_name, plot_factory, 'eta', out_format)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    __main__()

#------------------------------------------------------------------------------- 
