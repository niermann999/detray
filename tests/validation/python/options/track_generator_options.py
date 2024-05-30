# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

import argparse

#-------------------------------------------------------------------------------
# Options parsing
#-------------------------------------------------------------------------------

""" Parent parser that contains generic track generator options """
def track_generator_options():

    parser = argparse.ArgumentParser(add_help=False)

    # Navigation options
    parser.add_argument("--n_tracks", "-n",
                        help=("Number of test tracks"),
                        default = 100, type=int)
    parser.add_argument("--transverse-momentum", "-p_T",
                        help=("Transverse momentum of the test particle [GeV]"),
                        default = 10, type=float)
    parser.add_argument("--eta_range", "-eta", nargs=2,
                        help=("Eta range of generated tracks"),
                        default = [-4, 4], type=float)

    return parser

""" Parent parser that contains options for the uniform track generator """
def uniform_track_generator_options():

    parser = argparse.ArgumentParser(add_help=False)

    # Navigation options
    parser.add_argument("--phi_steps",
                        help=("Number of steps in phi"),
                        default = 100, type=int)
    parser.add_argument("--eta_steps",
                        help=("Number of steps in eta"),
                        default = 100, type=int)
    parser.add_argument("--transverse-momentum", "-p_T",
                        help=("Transverse momentum of the test particle [GeV]"),
                        default = 10, type=float)
    parser.add_argument("--eta_range", "-eta", nargs=2,
                        help=("Eta range of generated tracks"),
                        default = [-4, 4], type=float)

    return parser

#-------------------------------------------------------------------------------
