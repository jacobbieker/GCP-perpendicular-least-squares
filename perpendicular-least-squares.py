__author__ = 'Jacob Bieker'
import os, sys, random
import numpy
from multiprocessing import Pool
from astropy.io import fits


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


