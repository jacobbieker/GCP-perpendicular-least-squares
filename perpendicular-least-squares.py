__author__ = 'Jacob Bieker'
import os, sys, random
import numpy
from astropy.table import Table, vstack
from astropy.io import fits
import copy
import scipy.odr as odr
from scipy.stats import linregress
from statsmodels.formula.api import ols


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


def min_quartile(total_galaxies):
    fits_table.sort("residual")
    low_num = total_galaxies / 4.0
    high_num = 3.0 * total_galaxies / 4.0
    fits_residual = fits_table.field("residual")
    very_low_num = fits_residual[int(low_num - 0.5)]
    tab_value = fits_residual[int(low_num + 0.5)]
    very_low_num = (very_low_num + tab_value) * 2.0
    very_high_num = fits_residual[int(high_num -0.5)]
    tab_value = fits_residual[int(high_num + 0.5)]
    very_high_num = (very_high_num + tab_value) * 2.0
    delta = very_high_num - very_low_num
    return delta


def min_delta(percentage, total_galaxies):
    residuals = fits_table.field("residual")
    absolute_residuals = []
    for residual in residuals:
        absolute_residual = abs(residual)
        absolute_residuals.append(absolute_residual)
    if percentage == 100:
        delta = numpy.mean(absolute_residuals)
        rms = numpy.std(residuals) * numpy.sqrt((len(residuals) - 1) / (len(residuals) - 3))
        return rms, delta
    elif percentage == 60:
        high_num = total_galaxies * 0.6 + 0.5
        absolute_residuals_60 = absolute_residuals[:int(high_num)]
        residuals_60 = residuals[:int(high_num)]
        delta = numpy.mean(absolute_residuals_60)
        rms = numpy.std(residuals_60)*numpy.sqrt((len(residuals_60) - 1) / (len(residuals_60) - 3))
        return rms, delta


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


def line_solve():
    # TODO: Find the Least Squares for a line
    return 0


def plane_solve():
    # TODO: Find the least Squares for a plane
    return 0


def read_clusters(filename, solve_plane, galaxy_name, group_name, y_col, x1_col, x2_col):
    '''
    Reads in the FITS file and creates one table from the different tables and counts the
    total number of galaxies, for use later
    :param filename: Name of the FITS file
    :param solve_plane: Boolean whether x2_col should be counted or not
    :param galaxy_name: Column name for the galaxy
    :param group_name: Column name for the group
    :param y_col: Column name for y
    :param x1_col: Column name for x1
    :param x2_col: Column name for x2 (optional)
    :return: The entire dataset in one table, with cluster number added to the data, and the total number
    of galaxies
    '''
    cluster_number = 0
    hdu_num = 0
    finished_table = Table()
    while (True):
        try:
            table = Table.read(filename, format="fits", hdu=hdu_num)
            cluster_number += 1
            hdu_num += 1
            if solve_plane:
                new_table = Table(table, names=(galaxy_name, group_name, y_col, x1_col, x2_col))
            else:
                new_table = Table(table, names=(galaxy_name, group_name, y_col, x1_col))
            cluster_num_column = fits.Column(cluster_number, name="cluster_number")
            new_table.add_column(cluster_num_column)
            finished_table = vstack(finished_table, new_table)
        except:
            break
    total_galaxies = len(finished_table)
    return finished_table, total_galaxies


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


def bootstrap_initial(num_bootstrap):
    flow_boot = False
    if num_bootstrap > 0:
        total_boot = num_bootstrap
        ssa = 0
        sa = 0
        ssb = 0
        sb = 0


def bootstrap_cluster():
    # Fit cluser with the most data to get an idea of where to iterate from
    '''
    for(n_i=1 ; n_i<=n_clus ; n_i+=1) {
 tstat(tmpall,"nclus",lowlim=n_i,highlim=n_i, >> "/dev/null")
    :param cluster:
    :param num_bootstrap:
    :return:
    '''
    rich_cluster = 0
    rich_members = 0
    cluster_table = fits_table.group_by('cluster_number')
    # Selects all the rows with the same cluster number and counts them to figure out which has the most
    for key, group in zip(cluster_table.groups.keys, cluster_table.groups):
        if len(group) > rich_members:
            rich_members = len(group)
            rich_cluster = key['cluster_number']

    # Fitting cluster to get rid of any NaN points
    if solve_plane:
        '''
    tselect(tmpall,tmpsel,"nclus=="//n_rich//" && "//n_recol//"<99999. && "//n_sigcol//"<99999.")
    else
        tselect(tmpall,tmpsel,"nclus=="//n_rich//" && "//n_recol//"<99999. && "//n_sigcol//"<99999. && "//n_Iecol//"<99999.")
    tinfo(tmpsel,ttout-)
    if(n_flprint)
      printf("Number of galaxies fit in this cluster:  %3d\n",tinfo.nrows)
    tfitlin(tmpsel,n_recol,n_sigcol,n_Iecol,rows="-",verbose=no, > tmpout)
        # get the RMA coefficients
        '''
    if res_choice == "y":
        # On the 6th line from tfitlin, which goes to STDIN, fields takes the 3rd and 4th whitespace
        # separated values in line 6
        # Scan scans in those values into n_a and n_b
       # head(tmpout,nlines=6) | fields("STDIN","3-4",lines="6") | \
        # scan(n_a,n_b)
    else:
        # On the 3rd to last line from tfitlin, which goes to STDIN, fields takes the 2nd and 3rd whitespace
        # separated values
        # tail(tmpout,nlines=3) | fields("STDIN","2-3",lines="1") | \ scan(n_a, n_b)
    if solve_plane:
        b_in = 0.0

    if solve_plane:
        # If two parameter fit make the face zero column
        fits_table[x2_col] = 0.0


    # TODO: bootstrap a cluster to get a different distribution to check with check_guess
    return 0


def determine_uncertainty(solutions):
    # TODO: Take a list or dict of solutions and determine the uncertainty in them
    return 0


def change_coefficients(a_factor_in, b_factor_in, a_iterations, max_iterations, b_iterations, restart_factor, flow_print, flow_a, flow_b):
    m_a = a_factor_in / 200.0
    m_b = b_factor_in / 200.0
    if (a_factor <= m_a or a_iterations > max_iterations) and ((b_factor <= m_b or b_iterations > max_iterations) or solve_plane ) or ( a_iterations > 3 * max_iterations or b_iterations > 3 * max_iterations ):
        if not restart_factor:
            flow_end = True
            flow_print = True
            a_in = a_out
            b_in = b_out
            next_res()
        else:
            restart_factor = False
            a_in = a_out
            b_in = b_out
            a_factor_in /= 100.0
            b_factor_in /= 100.0
            max_iterations /= 2.0
            if flow_print:
                '''
                if(n_flprint) {
                      print("")
                      printf("Restarting with (a,b)=(%7.4f,%7.4f)\n",n_a,n_b)
                      printf("    (n_facta,n_bfact)=(%7.4f,%7.4f)\n",n_factain,n_factbin)
                      print("")
                    }
                '''
            restart()

    # Change Coefficients
    if flow_a:
        a_in *= 1.0 + sig_a * a_factor
        a_iterations += 1
    if flow_b:
        b_in *= 1.0 * sig_b * b_factor
        b_iterations += 1
    next_res()

    return 0


def restart():
    return 0


def next_res():
    # TODO Use to restart process again?
    return 0


def determine_change_coefficients(a_iterations, b_iterations, a_factor, b_factor, delta_out, a_in, b_in, delta_in, sig_a, sig_b, very_low_in, very_high_in, minimization_algorithm, flow_a, flow_b,):
    if a_iterations == 1 and b_iterations == 1:
        a_out = a_in
        b_out = b_in
        delta_out = delta_in
        if minimization_algorithm == "quartile":
            very_low_out = very_low_in
            very_high_out = very_high_in
    else:
        if delta_in <= delta_out and a_factor > a_factor_in / 200.0 and (b_factor > b_factor_in / 200.0 or solve_plane):
            a_out = a_in
            b_out = b_in
            delta_out = delta_in
            flow_first = False
            if minimization_algorithm == "quartile":
                very_low_out = very_low_in
                very_high_out = very_high_in
            change_coefficients()
        if delta_in > delta_out and a_factor > a_factor_in / 200.0 and (b_factor > b_factor_in / 200.0 or solve_plane):
            # Change the current coefficients back to previous values
            a_in = a_out
            b_in = b_out
            if flow_a and (a_iterations == 1 or flow_first):
                sig_a = -sig_b
                flow_first = False
                change_coefficients()
            if flow_b and (b_iterations == 1 or flow_first):
                sig_b = -sig_b
                flow_first = False
                change_coefficients()
            if (a_iterations > 1 or b_iterations > 1) and not flow_first:
                if not solve_plane:
                    # Change the other coefficient
                    flow_a = not flow_a
                    flow_b = not flow_b
                flow_first = True
                if flow_a and a_iterations > 1:
                    sig_a = -sig_a
                    a_factor /= 2.0
                if flow_b and b_iterations > 1:
                    b_factor /= 2.0
                change_coefficients()


def tfitlin(table):
    # Fit line or plane by 3 or 4 different methods
    #
    # Fit:  y = a1*x+b
    #       y = a1*x+a2*x+b
    #
    # Methods: Minimizing the residuals in y
    #          Minimizing the residuals in x1
    #          Minimizing the residuals in x2 (for plane)
    #          Reduced major axis method
    #              a=(a1*a2)**0.5   b=(b1*b2)**0.5  (line)
    #              a1=(a1_1*a1_2*a1_3)**(1/3)       (plane)
    #              a2=(a2_1*a2_2*a2_3)**(1/3)       (plane)
    #              b=(b1*b2*b3)**(1/3)              (plane)
    #
    # =====> No weighting <========
    #
    # Packages used: noao, proto, testphot ttoolsx
    # Fortran program: fitplan    This is where the job is done!
    #
    # Version: 13.03.92  IRAF 2.9 sun/mira
    #          12.01.93  iraf 2.10 mira/sun
    #          16.10.95  Solaris roeskva
    #          30.03.04  Changed to use redhat binary, works under IRAF 2.12.2 Redhat 9.0
    #          19.05.08  moved to Mac OS freja, changed to use /Users/inger/bin path
    #          30.11.16  Moved to Python, Jacob Bieker
    # Inger Jorgensen, Gemini Observatory
    # e-mail: ijorgensen@gemini.edu

    def f(B, x):
        '''Linear function y = m*x + b'''
        # B is a vector of the parameters.
        # x is an array of the current x values.
        # x is in the same format as the x passed to Data or RealData.
        #
        # Return an array in the same format as y passed to Data or RealData.
        return B[0]*x + B[1]

    def f3(B, x):
        return B[0]*x + B[1]*x +B[2]

    intable = str(input("Input Table: "))
    y_col = str(input("Y column name, log r") or "logre_r")
    x1_col = str(input("X1 column name, log sigma") or "logsig_lit")
    x2_col = str(input("X2 column name, surface brightness") or "<mu_r>e")
    rows = str(input("Rows: ") or "")
    verbose = bool(input("Verbose") or False)

    if not solve_plane:
        x = table[0]
        y = table[1]
        guess = linregress(x, y)
        mod = odr.Model(f)
        dat = odr.Data(x, y)
        od = odr.ODR(dat, mod, beta0=guess[0:2])
        out = od.run()
    else:
        # For the three dimensional case, the y has to be a matrix for lstsq to work
        x = table[0]
        y = table[1]
        guess = numpy.lstsq(x, y)


if __name__ == "__main__":

    a_factor = 0
    b_factor = 0
    sig_a = 0
    sig_b = 0
    a_factor_in = 0
    b_factor_in = 0
    a_iterations = 0
    b_iterations = 0
    max_iterations = 0
    a_out = 0
    b_out = 0
    restart_factor = False

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
    fits_table, total_galaxies = read_clusters(filename)
    # Checks for which variables and functions to call
    if not x2_col:
        # Only use two parameters
        factor_change_b = 0.0
        solve_plane = False
        line_solve()
    else:
        plane_solve()

    print(hdulist[1].data.columns)
    # print(hdulist[1].data)
    # min_delta(filename="rxj1226allfit.fits", percentage=100)
