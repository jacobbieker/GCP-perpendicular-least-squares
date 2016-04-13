__author__ = 'Jacob Bieker'
import os, sys, random
import numpy
from multiprocessing import Pool
from astropy.io import fits

# Fit plane or line iteratively
#   if more than one cluster, each cluster in separate STSDAS table
#
# Residuals are calculated according to resalgo
#      per:  perpendicular to the fitted relation
#            delta = (y - a*x1 - b*x2 - c)/sqrt(1+a^2+b^2)
#      y  :  in y     delta = (y - a*x1 - b*x2 - c)
#      x1 :  in x1
#      x2 :  in x2
#
# Minization algorithms:  (minalgo)
#    quartile : distance between upper (0.75) quartile point
#               and lower (0.25) quartile point
#    delta100 : SUM abs(res)
#    delta60  : SUM abs(res) of the central 60 percent of distribution
#    rms100   : SUM res^2
#    rms60    : SUM res^2 of the central 60 percent of distribution
#
# Zero points : median (use with delta and quartile)
#               mean   (use with rms)
#
# Bootstrap uncertainties on coefficients if nboot>0
# Seed for random is in each step   seed-(nboot-step)
#
# WARNING: This is a rather slow algorithm.
# FP for 200 galaxies in 10 clusters takes 5min on Sparc 10
# if restart=yes and nboot=0 (no uncertainties)
#
# Version: 21.06.95
#          16.10.95  Solaris roeskva
#          14.02.96  changed resalgo=y initial values, faster conv.
#          10.11.97  use rename instead of trename in 2.11
#          21.01.98  bootstrap changed to the proper way of doing it
#          06.03.98  fixed bug in bootstrap for 2 parameter fits
#                    improved random number generation
#          27.07.98  stop if itera>3*maxiter || niterb>3*maxiter
#                    no matter the other parameter's niter
#          30.03.04  handle INDEF automatically, works under IRAF 2.12.2 Redhat9.0
# Inger Jorgensen, Gemini Observatory
# e-mail: inger@gemini.edu

def random_number(number, seed):
    rand_nums = []
    if seed > 0:
        seed = -seed
    random.seed(a=seed)
    for i in range(number):
        rand_num = random.randint(0,1)
        rand_nums.append(rand_num)
    return rand_nums


def residuals_perpendicular():
    # TODO: Convert delta = (y - a*x1 - b*x2 - c)/sqrt(1+a^2+b^2)
    return 0


def residuals_y():
    # TODO: Convert: delta = (y - a*x1 - b*x2 - c)
    return 0


def residuals_x1():
    return 0


def residuals_x2():
    return 0


def line_solve():
    # TODO: Find the Least Squares for a line
    return 0


def plane_solve():
    # TODO: Find the least Squares for a plane
    return 0


def read_clusters(*args):
    # TODO: Read in the files containing the cluster points
    return 0


def determine_type():
    # TODO: Determine the type (i.e. Plane or Line) to solve for
    return 0


def initial_guess(largest_cluster, type_solution):
    if type_solution == "line":
        return 0
    elif type_solution == "plane":
        return 0
    # TODO: Get an inital guess from the perpendicular least squares method from the largest cluster
    return


def check_guess(cluster, type_solution, guess):
    if type_solution == "line":
        return 0
    elif type_solution == "plane":
        return 0
    # TODO: check guess with another cluster


def bootstrap_cluster(cluster):
    # TODO: bootstrap a cluster to get a different distribution to check with check_guess
    return 0


def determine_uncertainty(solutions):
    # TODO: Take a list or dict of solutions and determine the uncertainty in them
    return 0

if __name__ == "__main__":
    iterations = 0
    galaxy_name = ""
    group_name = ""
    factor_change_a = 0.05
    factor_change_b = 0.02
    restart_factor = True
    num_bootstrap = 0
    rand_seed = 1
    zeropoint_choice = "median" # or mean
    min_distance = "delta100" # or quartile, delta60, rms100, rms60
    min_residual = "per" # or y, x1, x2
    y_col = "lre_GR_sc" # prompt for it
    x1_col = "lsig_re" # prompt for it
    x2_col = "lIe_GR_sccor" # prompt for it (optional)
    list_clusters = [] # prompt "List of input STSDAS tables"

    filename = input("Enter the filename containing the cluster(s): ")
    hdulist = fits.open(filename)
    print(hdulist.info())
    print(repr(hdulist[0].header))


