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

# vector sought after n = ( - 1, a, b)
# Use median zeropoint, not mean
# finding the residual: delta =(logre -a*1og(J) - b*log<I>)/(1 + a^2 + b^2)^(1/2).


def random_number(number, seed, nboot):
    rand_nums = []
    if seed <= 0 or seed <= nboot:
        seed = max(seed*seed, (seed+1)*(nboot+1))
    elif seed > 0:
        seed = -seed
    random.seed(a=seed)
    for i in range(number):
        rand_num = random.random()
        rand_nums.append(rand_num)
    return rand_nums


def min_quartile(cluster):
    return 0


def min_delta(filename, percentage):
    if percentage == 100:
        residuals = fits.getdata(filename=filename, extname="residual")
        absolute_residuals = []
        for residual in residuals:
            absolute_residual = abs(residual)
            absolute_residuals.append(absolute_residual)
        delta = numpy.mean(absolute_residuals)
        #TODO calculate sqrt( (tstat.nrows-1.)/(tstat.nrows-3.) )
        rms = numpy.std(absolute_residuals) * numpy.sqrt
        return 0
    elif percentage == 60:
        return 0


def min_rms(cluster, percentage):
    if percentage == 100:
        return 0
    elif percentage == 60:
        return 0


def zeropoint(cluster, type_solution):
    if type_solution.lower() == "median":
        # use with delta and quartile
        return 0
    elif type_solution.lower() == "mean":
        # use with rms
        return 0


def residuals_perpendicular():
    # TODO: Convert delta = (y - a*x1 - b*x2 - c)/sqrt(1+a^2+b^2)
    return 0


def residuals_y():
    # TODO: Convert: delta = (y - a*x1 - b*x2 - c)
    # The calculated residuals go into the "res" column
    return 0


def residuals_x1():
    # The calculated residuals go into the "res" column
    return 0


def residuals_x2():
    # The calculated residuals go into the "res" column
    return 0


def line_solve():
    # TODO: Find the Least Squares for a line
    return 0


def plane_solve():
    # TODO: Find the least Squares for a plane
    return 0


def read_clusters(fits_file):
    # TODO: Read in the files containing the cluster points
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


def bootstrap_cluster(cluster, nboot):
    # TODO: bootstrap a cluster to get a different distribution to check with check_guess
    return 0


def determine_uncertainty(solutions):
    # TODO: Take a list or dict of solutions and determine the uncertainty in them
    return 0

if __name__ == "__main__":
    filename = str(input("Enter the filename containing the cluster(s): ")).strip()
    tables = str(input("List of input STSDAS tables (e.g. Table1 Table2 Table3): ")).strip()
    min_choice = str(input("Distance to minimize (delta100,delta60,rms100,rms60,quartile): ")).strip() or "delta100"
    res_choice = str(input("Residual to minimize (per,y,x1,x2): ")).strip() or "per"
    y_col = str(input("Column name for y: ")).strip() or "lre_GR_sc"
    x1_col = str(input("Column name for x1: ")).strip() or "lsig_re"
    x2_col = str(input("Column name for x2 (optional): ")).strip()
    zeropoint_choice = input("Zeropoints (median, mean): ") or "median"
    galaxy_name = str(input("Column name for galaxy: ")).strip()
    group_name = str(input("Column name for group: ")).strip()
    factor_change_a = float(input("Starting factor for changes in a: ") or 0.05)
    factor_change_b = float(input("Starting factor for changes in b: ") or 0.02)
    iterations = int(input("Maximum number of iterations: ") or 0)
    restart_factor = bool(input("Restart iteration with smaller factors: ") or True)
    num_bootstrap = int(input("Number of estimates for bootstrap: ") or 0)
    rand_seed = int(input("Seed for random used in bootstrap: ") or 1)
    rand_num = int(input("Number of random numbers: ") or 1)

    # preprocess input
    hdulist = fits.open(filename)
    print(hdulist.info())
    print(hdulist[1].dump(datafile="comafit.txt", cdfile="comacolumn.txt", hfile="comaheader.txt", clobber=True))
    list_temp = tables.split(" ")
    list_clusters = [x for x in list_temp if x.strip()]
    random_numbers = random_number(number=rand_num, seed=rand_seed, nboot=num_bootstrap)
    print(random_numbers)
    # Checks for which variables and functions to call
    if not x2_col:
        # Only use two parameters
        factor_change_b = 0.0
        line_solve()
    else:
        plane_solve()

    print(hdulist[1].data.columns)
    #print(hdulist[1].data)




