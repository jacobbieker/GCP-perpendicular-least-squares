__author__ = 'Jacob Bieker'
import os, sys, random
import numpy
from astropy.table import Table
from astropy.io import fits
import copy


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
        seed = max(seed * seed, (seed + 1) * (nboot + 1))
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
        residuals = fits_table.field("residual")
        absolute_residuals = []
        for residual in residuals:
            absolute_residual = abs(residual)
            absolute_residuals.append(absolute_residual)
        delta = numpy.mean(absolute_residuals)
        # TODO calculate sqrt( (tstat.nrows-1.)/(tstat.nrows-3.) )
        rms = numpy.std(absolute_residuals) * numpy.sqrt((len(absolute_residuals) - 1) / (len(absolute_residuals) - 3))
        return 0
    elif percentage == 60:
        return 0


def min_rms(percentage):
    residuals = fits_table.field("residual")
    if percentage == 100:
        rms = numpy.std(residuals) * numpy.sqrt((len(residuals) - 1) / (len(residuals) - 3))
        return rms
    elif percentage == 60:
        lower_num = 0.2 * len(residuals) + 0.5
        higher_num = 0.8 + len(residuals) + 0.5
        residuals_60 = residuals[int(lower_num):int(higher_num)]
        rms = numpy.std(residuals_60)*numpy.sqrt((len(residuals_60) - 1.0) / (len(residuals_60) - 3.0))
        return rms


def zeropoint(cluster, type_solution, res_choice, y_col, x1_col, x2_col, a_factor, b_factor):
    """

    # derive the zero points and calculate the residuals
    # the zero points are median or mean as defined by n_zeropoint
    n_norm=1.                         # min in y
    if(n_resalgo=="per")
      n_norm=sqrt(1.+n_a**2+n_b**2)   # min perpendicular
    if(n_resalgo=="x1")
      n_norm=-1.*n_a                  # min in x1
    if(n_resalgo=="x2")
      n_norm=-1.*n_b                  # min in x2

    tcalc(tmpall,"res","0.",colfmt="f6.3")

    for(n_i=1 ; n_i<=n_clus ; n_i+=1) {
    # delta y
     n_expression = "("//n_recol//"-"//n_a//"*"//n_sigcol//"-"//n_b//"*"//n_Iecol//")"
     print(n_expression, > tmpexp)
     tcalc(tmpall,"z"//n_i,"@"//tmpexp,colfmt="f6.3")
     delete(tmpexp,verify=no)
     n_expression="z"//n_i//"*(nclus=="//n_i//")+1000.*(nclus!="//n_i//")"
     print(n_expression, >tmpexp)
     tcalc(tmpall,"z"//n_i,"@"//tmpexp,colfmt="f6.3")
     delete(tmpexp,verify=no)
     tstat(tmpall,"z"//n_i,lowlim=INDEF,highlim=100., >> "/dev/null")
     if(n_zeropoint=="median")
       n_zero=tstat.median
     else
       n_zero=tstat.mean
    # printf("Zero point for cluster  %-3d : %8.5f\n",n_i,n_zero)
    # residuals normalized
     n_expression="((z"//n_i//"-"//n_zero//")*(nclus=="//n_i//"))/"//n_norm//"+1000.*(nclus!="//n_i//")"
     print(n_expression, > tmpexp)
     tcalc(tmpall,"r"//n_i,"@"//tmpexp,colfmt="f6.3")
     delete(tmpexp,verify=no)
     n_expression="res+((z"//n_i//"-"//n_zero//")*(nclus=="//n_i//"))/"//n_norm
     print(n_expression, > tmpexp)
     tcalc(tmpall,"res","@"//tmpexp,colfmt="f6.3"
     delete(tmpexp,verify=no)
    }

    RMA coefficients: a_factor, b_factor
    :param cluster:
    :param type_solution:
    :return:
    """

    # Adds a column full of zeros to the FITS table for use in residual
    residual_column = fits.Column(name='res', format='f6.3', array=0)
    fits_table.add_column(residual_column)

    n_norm = 1.  # min in y
    if res_choice == "per":
        n_norm = numpy.sqrt(1.0 + a_factor ** 2 + b_factor ** 2)  # min perpendicular
    if res_choice == "x1":
        n_norm = -1.0 * a_factor  # min in x1
    if res_choice == "x2":
        n_norm = -1.0 * b_factor  # min in x2

    for index, galaxy in enumerate(cluster):
        zeropoint_dict = {}
        # expression "//n_recol//"-"//n_a//"*"//n_sigcol//"-"//n_b//"*"//n_Iecol//"
        # n_recol = y1, n_Iecol = x2_col, n_sigcol = x1_col
        n_recol = fits_table[y_col]
        n_Iecol = fits_table[x2_col]
        n_sigcol = fits_table[x1_col]
        expression = n_recol - a_factor * n_sigcol - b_factor * n_Iecol
        zeropoint_dict["z" + str(index)] = expression
        non_cluster_residual = fits_table[:index] + fits_table[
                                                         index + 1:]  # Potentially works, if it is a list
        zeropoint_dict["z" + str(index)] = ((zeropoint_dict["z" + str(index)]) * (fits_table[index])) + 1000.0 * non_cluster_residual
        # Ignore z values that are above 100.0
        temp_zeropoint_dict = copy.deepcopy(zeropoint_dict)
        temp_zeropoint_dict[x > 100.0] = float("NaN")
        if type_solution.lower() == "median":
            # use with delta and quartile
            n_zero = numpy.nanmedian(temp_zeropoint_dict)
        elif type_solution.lower() == "mean":
            # use with rms
            n_zero = numpy.nanmean(temp_zeropoint_dict)

            # printf("Zero point for cluster  %-3d : %8.5f\n",n_i,n_zero)
            # residuals normalized
        residuals(n_norm, index, n_zero, zeropoint_dict)


def residuals(n_norm, cluster_number, n_zero, zeropoint_dict):
    '''
     n_expression="((z"//n_i//"-"//n_zero//")*(nclus=="//n_i//"))/"//n_norm//"+1000.*(nclus!="//n_i//")"
 print(n_expression, > tmpexp)
 tcalc(tmpall,"r"//n_i,"@"//tmpexp,colfmt="f6.3")
 delete(tmpexp,verify=no)
 n_expression="res+((z"//n_i//"-"//n_zero//")*(nclus=="//n_i//"))/"//n_norm
 print(n_expression, > tmpexp)
 tcalc(tmpall,"res","@"//tmpexp,colfmt="f6.3"
 delete(tmpexp,verify=no)
    :param n_norm:
    :return:
    '''
    non_cluster_residual = fits_table[:cluster_number] + fits_table[
                                                         cluster_number + 1:]  # Potentially works, if it is a list
    residual_data = ((zeropoint_dict[cluster_number] - n_zero) * (
    fits_table[cluster_number])) / n_norm + 1000.0 * non_cluster_residual
    # fits_data[nclud] gets row, if that's what needed
    residual_number_col = fits.Column(name='r' + str(cluster_number), format='f6.3', array=residual_data)
    # TODO: tcalc(tmpall,"r"//n_i,"@"//tmpexp,colfmt="f6.3") This seems to put the result into the residual, numbered by the cluster number, so neeed to add "residuals" thing
    res_zeropoint = fits_table['res'] + ((zeropoint_dict[cluster_number] - n_zero) * (
    fits_table[cluster_number])) / n_norm
    # TODO: next tcalc puts it in the general residual line, so that is used the most
    residual_number_all = fits.Column(name='res', format='f6.3', array=residual_data)
    # new_columns = fits.ColDefs([residual_number_col, residual_number_all])
    # new_table_hdu = fits.new_table(hdulist.columns + new_columns)
    fits_table.add_column(residual_number_col)
    fits_table.add_column(residual_number_all)
    # Debug stuff

    print("\n\n\n---------- Residual Things-----------\n")
    print("Table Res\n")
    print(fits_table['res'])
    print("Table r<Cluster Number>\n")
    print(fits_table[cluster_number])
    print("\n\n\n---------- Residual Things-----------\n\n\n\n")
    status = 0
    return status


def residuals_perpendicular(y, x1, x2, ):
    n_norm = numpy.sqrt(1.0 + n_a ** 2 + n_b ** 2)
    # TODO: Convert delta = (y - a*x1 - b*x2 - c)/sqrt(1+a^2+b^2)
    return n_norm


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
    fits_data = hdulist[1]
    list_temp = tables.split(" ")
    list_clusters = [x for x in list_temp if x.strip()]
    random_numbers = random_number(number=rand_num, seed=rand_seed, nboot=num_bootstrap)
    print(random_numbers)
    fits_table = Table.read(filename, format="fits")
    # Checks for which variables and functions to call
    if not x2_col:
        # Only use two parameters
        factor_change_b = 0.0
        line_solve()
    else:
        plane_solve()

    print(hdulist[1].data.columns)
    # print(hdulist[1].data)
    # min_delta(filename="rxj1226allfit.fits", percentage=100)
